# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 19:03:19 2016

Initial atom coordinates for NPT_cutter procedure.

@author: sleepingz
"""
def NPT_cutter_sphere(m,proj_name,D,T,P):
    '''
        D (nm): void diameter
        T (K): temperature, P (bar): pressure    
    '''
    import sys
    sys.path.append(m.config['keros_path'])
    sys.path.append(m.config['recon_script_path'])
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
    