# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 19:03:19 2016
1.NPT_cutter_sphere:
    Initial atom coordinates for NPT_cutter procedure.
2.Saturation:
    Saturation by LJ atoms
3.Ads_lj:
    Single adsorption simulation case, with LJ potential
4.Isothermal:
    A series of adsorption simulations with a specific method at various pressures
@author: sleepingz
"""
import numpy as np
import matplotlib.pyplot as plt
def NPT_cutter_sphere(m,proj_name,D,T,P):
    '''
        D (nm): void diameter
        T (K): temperature, P (bar): pressure    
    '''
    m.Project(proj_name)
    nm2cm = 1e-7
    void = [D*nm2cm]*3
    void_v = 1./6*3.1416*(D*nm2cm)**3
    import update_db
    k = update_db.KeroCreate()
    import pack_mol
    N_mol_all, N_mol_net = pack_mol.Mol_Num_Estimate(m.config,\
            void, void_v,k.table)
    c = pack_mol.CubicPack(m.config, N_mol_all, N_mol_net, k.table)
    c.CubicPack()
    cell_min = pack_mol.np.min(c.cell)
    void_D = pack_mol.np.max(void)/1e-8# cm to A
    void_cell_req = int(void_D/cell_min) + 1
    nearest = pack_mol.NearestN(c.cube_length,void_cell_req)
    void_order = int((nearest + 1)/2)
    import logging
    logging.basicConfig(level = logging.DEBUG,filename = 'pack.log',\
	    filemode = 'w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)
    logging.info("Void making parameters:")
    logging.info("Void_radius: %.2f A, void_volume: %s cm^3"\
            %(void_D/2.0, format(void_v,'.2e')))
    logging.info("Required molecule number:")
    logging.info("all: %d, net: %d, adopted cube length: %d"\
            %(N_mol_all, N_mol_net, c.cube_length))
    logging.info("Suggested void order: %d"%void_order)
    logging.info("Sure to init the void?(y/n)")
    sure = raw_input()
    if sure == 'y':
        c.MakeVoid(void_order = void_order)
        c.VoidRadius()
        logging.info("Agree")
        logging.info("Final void radius:%.2f (A)"%c.void_radius)
        c.WriteData()
        datafile =  'data.%dmol.dreiding'%c.cube_length**3
        id_max = c.id_max
        import Recon
        recon = Recon.Recon_NPT_cutter(datafile,m.config,\
            void_D/2.0,void_D/2.0+1.0,T=T,P=P)
        recon.script_mod()
        recon.cutter_id_mod(id_max)
    else:
	logging.info("Deny")
    return '', 0    

def Saturation(m,proj_name,datafile,stepRun=40000):
    """
        datafile: datafile for the structure to be saturated
        stepRun: simulation # of steps for the saturation
    """
    import ads_pre
    import copy
    m.Project(proj_name)
    step = ads_pre.FinalStep('recon.lammpstrj')
    sat_dict = copy.deepcopy(ads_pre.sat_default)
    sat_dict['dump'] = 'recon.lammpstrj'
    sat_dict['data'] = datafile
    sat_dict['stepRun'] = stepRun
    sat_dict['step'] = step
    sat = ads_pre.ads_saturation(m.config,sat_dict)
    
def Ads_lj(m,proj_name,datafile,fluid,void_radius,T,P,ensemble_frame,\
    exchange_id = 8):
    '''
        Must be invoked after Saturation.
        datafile: datafile for the structure
        fluid: fluid name as 'methane'
        void_radius(A): pore void space radius
        T (K): temperature, P (MPa): pressure
        ensemble_frame: # of snapshots saved to the dump file
        exchange_id: (id of the MC exchange atoms) = (# of types in the datafile) +1
    '''
    import ads,ads_pre
    import copy,os
    m.Project(proj_name)
    convert = ads.P2mu(m.config)
    fluid_info = ads.fluid_info[fluid]
    f_mass = fluid_info['molMass']
    f_lj_e = fluid_info['atomEpsilon']
    f_lj_sig = fluid_info['atomSigma']
    mu = convert.conv(fluid,f_mass,T,P)
    
    ads_dict = copy.deepcopy(ads.ads_default)
    ads_pre.PoreRefine("ads_pre.lammpstrj",void_radius,\
        poreThickness = ads_dict['poreThickness'])
    step = ads_pre.FinalStep('refine_ads_pre.lammpstrj')
    ads_dict['fluid'] = fluid
    ads_dict['dump'] = 'refine_ads_pre.lammpstrj'
    ads_dict['data'] = datafile
    ads_dict['step'] = step
    ads_dict['atomMass'] = f_mass
    ads_dict['atomEpsilon'] = f_lj_e
    ads_dict['atomSigma'] = f_lj_sig
    ads_dict['T'] = T
    ads_dict['mu'] = mu
    ads_dict['ensembleFrame'] = ensemble_frame
    ads_dict['poreRadius'] = void_radius + ads_dict['poreThickness']
    ads_dict['dumpEvery'] = ads_dict['stepRun']/ads_dict['ensembleFrame']
    
    N_est = 4./3*3.1416*void_radius**3
    vol_mol = f_lj_sig**3
    N_max = int(N_est/vol_mol)
    PL = 5.0#!!!This is a estimated value for PL
    Ne = int(N_max*P/(PL+P))
    ads_dict['N_exchange'] = max(Ne,20)#fix_gcmc exchanges at least 20 atoms in a step 
    ads_dict['N_move'] = ads_dict['N_exchange']
    
    import datafile_mod
    datafile_mod.addAtomType(datafile,exchange_id-1,exchange_id)
    ads_dict['data'] = 'data.dreiding'#new datafile generated by addAtomType
    ads_lj = ads.ads_lj(m.config,ads_dict)
    
    method = 'lj'
    dump = 'ads_pre.lammpstrj'
    path = 'iso:%s_T_%.2f_P_%.2f'%(method,T,P)
    try:
        os.mkdir(path)
    except:
        pass
    os.system('cp -rf ff_%s refine_%s ads_%s.lmp data.dreiding %s'%\
                (dump,dump,method,path))
    os.system('rm -f ff_%s refine_%s ads_%s.lmp data.dreiding'%\
                (dump,dump,method))

def Ads_ljbiM(m,proj_name,datafile,fluid1,fluid2,x1,void_radius,T,P,\
    ensemble_frame,exchange_id = 8):
    '''
        Must be invoked after Saturation.
        datafile: datafile for the structure
        fluid1,fluid2: fluid name as 'methane'
        x1: mol ratio of fluid1
        void_radius(A): pore void space radius
        T (K): temperature, P (MPa): pressure
        ensemble_frame: # of snapshots saved to the dump file
        exchange_id: (id of the 1st type of MC exchange atoms) = (# of types in the datafile) +1
    '''
    import ads,ads_pre
    import copy,os
    m.Project(proj_name)
    convert = ads.P2mu(m.config)
    fluid_info = [ads.fluid_info[fluid1],ads.fluid_info[fluid2]]
    f_mass = [fluid_info[0]['molMass'],fluid_info[1]['molMass']]
    f_lj_e = [fluid_info[0]['atomEpsilon'],fluid_info[1]['atomEpsilon']]
    f_lj_sig = [fluid_info[0]['atomSigma'],fluid_info[1]['atomSigma']]
    fugacity = convert.biMix_f(fluid1,fluid2,T,P,x1)
    mu1 = convert.F2mu(f_mass[0],T,fugacity[0])
    mu2 = convert.F2mu(f_mass[1],T,fugacity[1])
    
    ads_dict = copy.deepcopy(ads.ads_default)
    ads_pre.PoreRefine("ads_pre.lammpstrj",void_radius,\
        poreThickness = ads_dict['poreThickness'])
    step = ads_pre.FinalStep('refine_ads_pre.lammpstrj')
    ads_dict['Name1'] = fluid1
    ads_dict['Name2'] = fluid2
    ads_dict['dump'] = 'refine_ads_pre.lammpstrj'
    ads_dict['data'] = datafile
    ads_dict['step'] = step
    ads_dict['atomMass1'],ads_dict['atomMass2'] = f_mass[0],f_mass[1]
    ads_dict['atomEpsilon1'],ads_dict['atomEpsilon2'] = f_lj_e[0],f_lj_e[1]
    ads_dict['atomSigma1'],ads_dict['atomSigma2'] = f_lj_sig[0],f_lj_sig[1]
    ads_dict['T'] = T
    ads_dict['mu1'],ads_dict['mu2'] = mu1,mu2
    ads_dict['ensembleFrame'] = ensemble_frame
    ads_dict['poreRadius'] = void_radius + ads_dict['poreThickness']
    ads_dict['dumpEvery'] = ads_dict['stepRun']/ads_dict['ensembleFrame']

    N_est = 4./3*3.1416*void_radius**3
    vol_mol = (0.5*(f_lj_sig[0]+f_lj_sig[1]))**3
    N_max = int(N_est/vol_mol)
    PL = 5.0#!!!This is a estimated value for PL
    Ne = int(N_max*P/(PL+P))
    ads_dict['N_exchange'] = max(Ne,20)
    ads_dict['N_move'] = ads_dict['N_exchange']
    ads_dict['N_swap'] = int(0.5*ads_dict['N_exchange'])
    
    import datafile_mod
    dataFile = datafile_mod.dataFile(datafile)
    dataFile.addAtomType(exchange_id-1,exchange_id+1)
    dataFile.output()
    ads_dict['data'] = 'data.dreiding'#new datafile generated by addAtomType
    ads_ljbiM = ads.ads_biM_lj(m.config,ads_dict)
    method = 'ljbiM'
    dump = 'ads_pre.lammpstrj'
    path = 'iso:%s_T_%.2f_P_%.2f_x1_%.2f'%(method,T,P,x1)
    try:
        os.mkdir(path)
    except:
        pass
    os.system('cp -rf ff_%s refine_%s ads_%s.lmp data.dreiding %s'%\
                (dump,dump,method,path))
    os.system('rm -f ff_%s refine_%s ads_%s.lmp data.dreiding'%\
                (dump,dump,method))

def Ads_molbiM(m,proj_name,datafile,fluid1,fluid2,x1,void_radius,T,P,\
    ensemble_frame,gcmcFreq,exchange_id = 8):
    import ads,ads_pre
    import copy,os
    m.Project(proj_name)
    convert = ads.P2mu(m.config)
    fluid_info = [ads.fluid_info[fluid1],ads.fluid_info[fluid2]]
    f_mass = [fluid_info[0]['molMass'],fluid_info[1]['molMass']]
    f_lj_e = [fluid_info[0]['atomEpsilon'],fluid_info[1]['atomEpsilon']]
    f_lj_sig = [fluid_info[0]['atomSigma'],fluid_info[1]['atomSigma']]
    fugacity = convert.biMix_f(fluid1,fluid2,T,P,x1)
    mu1 = convert.F2mu(f_mass[0],T,fugacity[0])
    mu2 = convert.F2mu(f_mass[1],T,fugacity[1])
    molecule_info = [ads.molecule_info[fluid1],ads.molecule_info[fluid2]]

    import datafile_mod
    dataFile = datafile_mod.dataFile(datafile)
    nMolAtomTypes = molecule_info[0]['ntypes']+molecule_info[1]['ntypes']
    molAtomID = range(exchange_id,exchange_id+nMolAtomTypes)
    molAtomID1 = molAtomID[:molecule_info[0]['ntypes']]
    molAtomID2 = molAtomID[molecule_info[0]['ntypes']:nMolAtomTypes]
    molAtomMass = molecule_info[0]['atomMasses']+molecule_info[1]['atomMasses']
    molAtomEpsilon = molecule_info[0]['atomEpsilons']+molecule_info[1]['atomEpsilons']
    molAtomSigma = molecule_info[0]['atomSigmas']+molecule_info[1]['atomSigmas']
    nMolBondTypes = len(molecule_info[0]['bondE'])+len(molecule_info[1]['bondE'])
    nMolAngleTypes = len(molecule_info[0]['angleE'])+len(molecule_info[1]['angleE'])
    molBondID = range(dataFile.nBondType+1,dataFile.nBondType+1+nMolBondTypes)
    molBondID1 = molBondID[:len(molecule_info[0]['bondE'])]
    molBondID2 = molBondID[len(molecule_info[0]['bondE']):nMolBondTypes]
    molAngleID = range(dataFile.nAngleType+1,dataFile.nAngleType+1+nMolAngleTypes)
    molAngleID1 = molAngleID[:len(molecule_info[0]['angleE'])]
    molAngleID2 = molAngleID[len(molecule_info[0]['angleE']):nMolAngleTypes]
    molBondE = molecule_info[0]['bondE']+molecule_info[1]['bondE']
    molBondLength = molecule_info[0]['bondLengh']+molecule_info[1]['bondLengh']
    molAngleE = molecule_info[0]['angleE']+molecule_info[1]['angleE']
    molAngleTheta = molecule_info[0]['angleTheta']+molecule_info[1]['angleTheta']
    
    from mol_template import mol_template
    mol_template(m,molecule_info[0]['template'],\
        bond_seq = molBondID1,angle_seq = molAngleID1)
    mol_template(m,molecule_info[1]['template'],\
        bond_seq = molBondID2,angle_seq = molAngleID2)
    addMass = '\n'.join(['mass %s %s'%(molAtomID[i],molAtomMass[i])\
        for i in range(nMolAtomTypes)])
    addPairCoef = '\n'.join(['pair_coeff %s %s %s %s'\
        %(molAtomID[i],molAtomID[i],molAtomEpsilon[i],molAtomSigma[i])\
        for i in range(nMolAtomTypes)])
    addBond = '\n'.join(['bond_coeff %s %s %s'\
        %(molBondID[i],molBondE[i],molBondLength[i]) for i in range(nMolBondTypes)])
    addAngle = '\n'.join(['angle_coeff %s %s %s'\
        %(molAngleID[i],molAngleE[i],molAngleTheta[i]) for i in range(nMolAngleTypes)])
    ads_dict = copy.deepcopy(ads.ads_default)
    ads_pre.PoreRefine("ads_pre.lammpstrj",void_radius,\
        poreThickness = ads_dict['poreThickness'])
    step = ads_pre.FinalStep('refine_ads_pre.lammpstrj')
    ads_dict['addMass'] = addMass
    ads_dict['addPairCoeff'] = addPairCoef
    ads_dict['addBond'] = addBond
    ads_dict['addAngle'] = addAngle
    ads_dict['Name1'] = fluid1
    ads_dict['Name2'] = fluid2
    ads_dict['mol1'] = molecule_info[0]['template']
    ads_dict['mol2'] = molecule_info[1]['template']
    ads_dict['dump'] = 'refine_ads_pre.lammpstrj'
    ads_dict['data'] = datafile
    ads_dict['step'] = step
    ads_dict['T'] = T
    ads_dict['mu1'],ads_dict['mu2'] = mu1,mu2
    ads_dict['molAtomTypes'] = ' '.join(map(str,molAtomID))
    ads_dict['offset1'] = molAtomID1[0]-1
    ads_dict['offset2'] = molAtomID2[0]-1
    ads_dict['ensembleFrame'] = ensemble_frame
    ads_dict['poreRadius'] = void_radius + ads_dict['poreThickness']
    ads_dict['gcmcFreq'] = gcmcFreq
    ads_dict['dumpEvery'] = gcmcFreq * 20
    ads_dict['thermoEvery'] = ads_dict['dumpEvery'] * 10 
    ads_dict['stepRun'] = ensemble_frame * ads_dict['dumpEvery']
<<<<<<< HEAD
=======
    ads_dict['coulCut'] = 2.0 * ads_dict['poreRadius']
>>>>>>> master

    N_est = 4./3*3.1416*void_radius**3
    vol_mol = (0.5*(f_lj_sig[0]+f_lj_sig[1]))**3
    N_max = int(N_est/vol_mol)
    PL = 5.0#!!!This is a estimated value for PL
    Ne = int(N_max*P/(PL+P))
    ads_dict['N_exchange'] = max(Ne,20)
    ads_dict['N_swap'] = int(0.5*ads_dict['N_exchange'])
    ads_dict['N_exchange'] = ads_dict['N_swap']
    
    dataFile.addAtomType(exchange_id-1,exchange_id-1+nMolAtomTypes)
    dataFile.addBondType(dataFile.nBondType,dataFile.nBondType+nMolBondTypes)
    dataFile.addAngleType(dataFile.nAngleType,dataFile.nAngleType+nMolAngleTypes)
    dataFile.addExtra()
    dataFile.output()
    #datafile_mod.addAtomType(datafile,exchange_id-1,exchange_id-1+nMolAtomTypes)
    ads_dict['data'] = 'data.dreiding'#new datafile generated by addAtomType
    ads_molbiM = ads.ads_biM_mol(m.config,ads_dict)
    method = 'molbiM'
    dump = 'ads_pre.lammpstrj'
    path = 'iso:%s_T_%.2f_P_%.2f_x1_%.2f'%(method,T,P,x1)
    try:
        os.mkdir(path)
    except:
        pass
    os.system('cp -rf ff_%s refine_%s ads_%s.lmp data.dreiding %s %s %s'%\
                (dump,dump,method,ads_dict['mol1'],ads_dict['mol2'],path))
    os.system('rm -f ff_%s refine_%s ads_%s.lmp data.dreiding %s %s '%\
                (dump,dump,method,ads_dict['mol1'],ads_dict['mol2']))
                
def Isotherm(m,proj_name,group_name,datafile,method,fluid,void_radius,\
    T,ensemble_frame=500,exchange_id = 8,Pmin=0.5,Pmax=30.0,num=10):
    '''
        Must be invoked after Saturation.
        datafile: datafile for the structure
        method: lj or mol
        fluid: fluid name as 'methane'
        void_radius(A): pore void space radius
        T (K): temperature
        ensemble_frame: # of snapshots saved to the dump file
        exchange_id: (id of the MC exchange atoms) = (# of types in the datafile) +1
        Pmin, Pmax and num decide the pressure values in the isotherm
    '''
    import ads_post
    import numpy as np
    P_seq = np.linspace(Pmin,Pmax,num=num,endpoint=True)
    for P in P_seq:
        if method == 'lj':
            Ads_lj(m,proj_name,datafile,fluid,void_radius,T,P,ensemble_frame,\
                exchange_id = exchange_id)
    m.Project(proj_name)
    ads_post.adsGroup('.',group_name)

def IsothermBiM(m,proj_name,group_name,datafile,method,fluid1,fluid2,\
        x1,void_radius,T,\
        ensemble_frame=500,exchange_id=8,gcmcFreq=100,Pmin=0.5,Pmax=30.0,num=10):
    import ads_post
    import numpy as np
    P_seq = np.linspace(Pmin,Pmax,num=num,endpoint=True)
    for P in P_seq:
        if method == 'ljbiM':
            Ads_ljbiM(m,proj_name,datafile,fluid1,fluid2,x1,void_radius,T,P,ensemble_frame,\
                exchange_id = exchange_id)
        if method == 'molbiM':
            Ads_molbiM(m,proj_name,datafile,fluid1,fluid2,x1,void_radius,T,P,ensemble_frame,\
                gcmcFreq,exchange_id = exchange_id)
    m.Project(proj_name)
    ads_post.adsGroup('.',group_name)

def DensityAccum(m,proj_name,accum_name,adsfile,ensemble_frame,equil_step,\
    void_radius,cutoff = 2.0,mesh = 50,exchange_id = 8,mode = 'full'):
    '''
        Must be invoked after Isotherm/Ads_lj.
        accum_name: folder name of the current DensityAccum
        adsfile: filename of the .lammpstrj file to be Accum
        ensemble_frame: # of snapshots saved to the dump file
        void_radius(A): pore void space radius
        mesh: # of mesh of the half_length, half_length = void_radius + 4.0 by default
        cutoff: space averaging cutoff radius of the coordnum
        equil_step: when step>equil_step the GCMC procedure is in equilibruim
        exchange_id: (id of the MC exchange atoms) = (# of types in the datafile) +1
        mode:
        --full: DensityAccum and Coordnum
        --partial: Coordnum only
    '''
    m.Project(proj_name)
    import ensembleDensity as eD
    
    half_length = void_radius + 4.0
    mesh = max(mesh,int(half_length/cutoff))
    
    D = eD.Density(m.config,path = accum_name,adsfile = adsfile)
    if mode == 'full':
        D.densityAccum(atomID = exchange_id, frame = ensemble_frame, equil = equil_step)
    
    D.coordNum(half_length, mesh, cutoff)
    D.finalize()
        
def IsothermAnalyze(m,proj_name,group_name,fluid,equil_step):
    """
        Must be called after Isotherm/Ads_lj AND LAMMPS runs
        group_name: name for the current isotherms
        fluid: fluid name
        equil_step: when step>equil_step the GCMC procedure is in equilibruim
    """
    m.Project(proj_name+'/%s'%group_name)
    import ads_post
    ads_post.writeAdsFile(m.config,fluid,'.',equil_step = equil_step)

def IsothermBiMAnalyze(m,proj_name,group_name,fluid1,fluid2,equil_step):
    m.Project(proj_name+'/%s'%group_name)
    import ads_post
    ads_post.writeBiMAdsFile(m.config,fluid1,fluid2,'.',equil_step = equil_step)
    
def VisualDensity(m,proj_name,case,phi=0.01):
    m.Project(proj_name)
    import ensembleDensity as eD
    D = eD.Density(m.config,path=case,mode='open')
    D.visualDensity(vis='yes',phi=phi)

def VisualDensityBiM(m,proj_name,case1,case2,phi=0.01):
    import os
    import ensembleDensity as eD
    m.Project(proj_name)
    grid_seq = []
    extent_seq = []
    for case in [case1,case2]:
        D=eD.Density(m.config,path=case,mode='open')
        grid_seq.append(D.visualDensity(vis='no',phi=phi))
        extent_seq.append((D.rxy_min,D.rxy_max,D.zmin,D.zmax))
        os.chdir('..')
    inf = 1e-3
    for grid in grid_seq:
        mask = np.ma.masked_less(grid,inf)
        grid = mask.filled(fill_value=0.0)
    vmax = max(np.max(grid_seq[0]),np.max(grid_seq[1]))
    import matplotlib.cm as cm
    plt.rc('font', family='serif')
    plt.rc(('xtick','ytick'),labelsize=15)
    plt.figure(figsize=(8,6))
    plt.subplot(121)
    plt.imshow(grid_seq[0].T,aspect='equal',origin='lower',\
        cmap=cm.Reds,alpha = 1.0,vmax=vmax,\
        extent=extent_seq[0])
    plt.title(case1)
    plt.subplot(122)
    plt.imshow(grid_seq[1].T,aspect='equal',origin='lower',\
        cmap=cm.Blues,alpha = 1.0,vmax=vmax,\
        extent=extent_seq[1])
    plt.title(case2)
    plt.show()

def VRhoSpectrum(m,proj_name,case,hist_mesh=100):
    m.Project(proj_name)
    import ensembleDensity as eD
    D = eD.Density(m.config,path=case,mode='open')
    D.readCoordfile()
    D.evalHist(mesh=hist_mesh)
    D.V_rho_Spectrum()

def Density1DRadial(m,proj_name,case,void_radius,mesh=100):
    m.Project(proj_name)
    import ensembleDensity as eD
    D = eD.Density(m.config,path=case,mode='open')
    D.density1D_radial(void_radius,mesh=mesh)
    print np.divide(D.radial_hist[1:],np.square(D.radial_r[1:]))
        
def Ads_lj93Wall(m,proj_name,void_radius,nAtoms,fluid,T,energy=1.9,sigma=3.7):
    m.Project(proj_name)
    import ads
    fluid_info = ads.fluid_info[fluid]
    f_mass = fluid_info['molMass']
    f_lj_e = fluid_info['atomEpsilon']
    f_lj_sig = fluid_info['atomSigma']
    lj93_dict = {'voidRadius':void_radius,\
            'poreRadius':void_radius + sigma,\
            'energy':energy,'nAtoms':nAtoms,'sigma':sigma,\
            'atomEpsilon':f_lj_e,'atomSigma':f_lj_sig,\
            'T':T,'atomMass':f_mass}
    lj93Wall = ads.ads_lj93Wall(m.config,lj93_dict)
    lj93Wall()

def AdsRhoFree(m,proj_name,void_radius,ads_group,fluid,T,energy=1.9,sigma=3.7,\
    hist_mesh = 100):
    m.Project(proj_name)
    import ads_post
    import ensembleDensity as eD
    ads_res_file = '%s/adsorption.data'%ads_group
    res = ads_post.readAdsFile(ads_res_file)
    nAtoms_seq = [int(r['N'][0]) for r in res]
    p_seq = [r['p'] for r in res]
    num = len(nAtoms_seq)
    new_proj_name = '%s/%s/adsRatio'%(proj_name,ads_group)
    m.Project(new_proj_name)
    f = open('rho_free.data','w')
    for i in range(num):
        Ads_lj93Wall(m,new_proj_name,void_radius,nAtoms_seq[i],\
                fluid,T,energy=energy,sigma=sigma)
        accum_name = 'accum_%.2fMPa'%p_seq[i]
        DensityAccum(m,new_proj_name,accum_name,"dump.ensembles",\
            100,200000,void_radius,exchange_id = 1)
        #curdir = new_proj
        D = eD.Density(m.config,path='VRho:%s'%(accum_name),\
            mode='open')
        D.readCoordfile()
        D.evalHist(mesh=hist_mesh)
        D.V_rho_Spectrum(plot='no')
        f.write('%8.2f,%8.4f\n'%(p_seq[i],D.RhoAd))   
    f.close()
    
def IsothermDensityAccum(m,proj_name,ads_group,method,\
        ensemble_frame,equil_step,\
        void_radius,cutoff = 2.0,mesh = 50,exchange = [8]):
    import os
    m.Project(proj_name)
    filenames = os.listdir(ads_group)
    for name in filenames:
        if name.startswith('iso:') or name.startswith('iso%'):
            new_proj_name = '%s/%s/%s'%(proj_name,ads_group,name)
            for ID in exchange:
                DensityAccum(m,new_proj_name,"AccumID%s"%ID,"ads_%s.lammpstrj"%method,\
                    ensemble_frame,equil_step,void_radius,\
                    cutoff = cutoff,mesh=mesh,exchange_id=ID,mode='full')

def AdsRatio(m,proj_name,ads_group,mode='self',hist_mesh=100,bonus=0.2):
    import os
    import ensembleDensity as eD
    m.Project(proj_name)
    filenames = os.listdir(ads_group)
    V_dict,RhoV_dict = {},{}
    rhoFree_dict = {}

    if mode == 'file':
        rhoFreefile = '%s/adsRatio/rho_free.data'%(ads_group)
        f = open(rhoFreefile)
        lines = f.readlines()
        for line in lines:
            p = line.strip().split(',')[0]
            rF = float(line.strip().split(',')[1])
            rhoFree_dict[p] = rF
        f.close()
    
    for name in filenames:
        m.Project(proj_name)
        if name.startswith('iso:') or name.startswith('iso%'):
            press = name.split('_')[-1]
            D = eD.Density(m.config,path='%s/%s/VRho:Accum'%(ads_group,name),\
                mode='open')
            D.readCoordfile()
            D.evalHist(mesh=hist_mesh)
            D.V_rho_Spectrum(plot='no')
            if mode == 'file':
                D.adsAnalysis(mode='file',rhoFreeRef=rhoFree_dict[press],bonus=bonus)
            else:
                D.adsAnalysis(mode='self',bonus=bonus)
                rhoFree_dict[press] = D.RhoFree
            RhoV_dict[press] = D.RhoV_ads
            V_dict[press] = D.V_ads
    m.Project('%s/%s'%(proj_name,ads_group))
    p_seq = sorted(V_dict.keys(),key=lambda x:float(x))
    g = open('adsRatio.data','w')
    g.write('press(MPa), RhoV, V, RhoFree\n')
    for p in p_seq:
        g.write('%8.2f,%8.4f,%8.4f,%8.4f\n'%(float(p),RhoV_dict[p],V_dict[p],rhoFree_dict[p]))
    g.close()

def AdsRatioLJ93(m,proj_name,ads_group,mode='self',hist_mesh=100,bonus=0.2):
    import os
    import ensembleDensity as eD
    m.Project(proj_name)
    filenames = os.listdir(ads_group+'/adsRatio')
    V_dict,RhoV_dict = {},{}
    rhoFree_dict = {}

    if mode == 'file':
        rhoFreefile = '%s/adsRatio/rho_free.data'%(ads_group)
        f = open(rhoFreefile)
        lines = f.readlines()
        for line in lines:
            p = line.strip().split(',')[0]
            rF = float(line.strip().split(',')[1])
            rhoFree_dict[p] = rF
        f.close()
    
    for name in filenames:
        m.Project(proj_name)
        if name.startswith('VRho'):
            press = name.split('_')[-1][:-3]
            D = eD.Density(m.config,path='%s/adsRatio/%s'%(ads_group,name),\
                mode='open')
            D.readCoordfile()
            D.evalHist(mesh=hist_mesh)
            D.V_rho_Spectrum(plot='no')
            if mode == 'file':
                D.adsAnalysis(mode='file',rhoFreeRef=rhoFree_dict[press],bonus=bonus)
            else:
                D.adsAnalysis(mode='self',bonus=bonus)
                rhoFree_dict[press] = D.RhoFree
            RhoV_dict[press] = D.RhoV_ads
            V_dict[press] = D.V_ads
    m.Project('%s/%s'%(proj_name,ads_group))
    p_seq = sorted(V_dict.keys(),key=lambda x:float(x))
    g = open('adsRatioLJ93.data','w')
    g.write('press(MPa), RhoV, V, RhoFree\n')
    for p in p_seq:
        g.write('%8.2f,%8.4f,%8.4f,%8.4f\n'%(float(p),RhoV_dict[p],V_dict[p],rhoFree_dict[p]))
    g.close()    

def VRhoSpectrumCompare(m,proj_name,ads_group,\
        case,void_radius,fluid,T,P,lj93e=2.0,hist_mesh=100):
    m.Project(proj_name)
    import ensembleDensity as eD
    import ads_post
    ads_res_file = '%s/adsorption.data'%ads_group
    res = ads_post.readAdsFile(ads_res_file)
    N_P = {}
    for r in res:
        N_P[r['p']]=r['N'][0]
    natoms = int(N_P[P])
    D = eD.Density(m.config,path=ads_group+'/'+case+'/VRho:Accum',mode='open')
    D.readCoordfile()
    D.evalHist(mesh=hist_mesh)
    D.V_rho_Spectrum(plot='no')
    new_proj_name = '%s/%s/Ads_lj93_temp'%(proj_name,ads_group)
    accum_name = 'accum_%.2fMPa_%s'%(P,lj93e)
    try:
        m.Project(new_proj_name)
        DLJ = eD.Density(m.config,path='VRho:%s'%(accum_name),\
            mode='open')
    except:
        Ads_lj93Wall(m,new_proj_name,void_radius,natoms,\
                fluid,T,energy=lj93e)
        DensityAccum(m,new_proj_name,accum_name,"dump.ensembles",\
            100,200000,void_radius,exchange_id = 1)
        DLJ = eD.Density(m.config,path='VRho:%s'%(accum_name),\
            mode='open')
    DLJ.readCoordfile()
    DLJ.evalHist(mesh=hist_mesh)
    
    import matplotlib.pyplot as plt
    plt.rc('font', family='serif')
    plt.rc('legend',numpoints=1)
    fig=plt.figure(figsize=(8,6))
    mesh_data = (2.0*DLJ.mesh)**3
    DLJ.data_hist = np.divide(DLJ.data_hist,mesh_data)
    plt.plot(np.divide(DLJ.data_x,DLJ.frame),\
        np.multiply(DLJ.data_hist,np.divide(D.data_x,D.frame)),'o-',\
        label=r"$\varepsilon_w = %s$"%lj93e)
    plt.xlabel(r'$\rho_E$',fontsize=20)
    plt.ylabel(r'$V_F$',fontsize=20)       
    plt.xlim(xmax=0.7)
    plt.plot(np.divide(D.data_x,D.frame),np.multiply(D.data_hist,np.divide(D.data_x,D.frame)),\
        'ko-',label="Kerogen Pore")
    plt.legend()
    print "LJ93 V:%.2f"%np.sum(DLJ.data_hist)
    print "Kerogen V:%.2f"%np.sum(D.data_hist)
    print "LJ93 N:%.2f"%np.sum(np.multiply(DLJ.data_hist,np.divide(DLJ.data_x,DLJ.frame)))
    print "Kerogen N:%.2f"%np.sum(np.multiply(D.data_hist,np.divide(D.data_x,D.frame)))
    fig.show()
