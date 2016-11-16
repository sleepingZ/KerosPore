# -*- coding: utf-8 -*-
import os
from dump import dump
import numpy as np
import copy

class Density():
    def __init__(self,config,adsfile = 'ads_lj.lammpstrj',\
        mode = 'new', path = 'densityTemp'):
        if mode == 'new':
            try:
                os.mkdir('densityTemp')
            except:
                pass
            os.chdir('./densityTemp')
        else:
            os.chdir(path)
            f = open('V_Rho.info','r')
            lines = f.readlines()
            self.frame = int(lines[0].strip().split(':')[1])
            self.mesh = int(lines[1].strip().split(':')[1])
            self.R_ave = float(lines[2].strip().split(':')[1])
            self.halflength = float(lines[3].strip().split(':')[1])
            
        self.mode = mode
        self.path = path
        self.config = config
        self.adsfile = '../'+adsfile
        self.engine = config['engine']+'<'
        self.dumpfile = 'ensembles.lammpstrj'
        self.coordnum = 'dump.coordnum'
    
    def finalize(self):#Must be called when a new Density class is about to construct
        f = open('V_Rho.info','w')
        f.write('frame:%d\n'%self.frame)
        f.write('mesh:%d\n'%self.mesh)
        f.write('R_ave:%.2f\n'%self.R_ave)
        f.write('halflength:%.2f\n'%self.halflength)
        f.close()        
        if self.mode == 'new':
            os.chdir('..')
            out_folder = 'VRho:%s'%(self.path)
            try:
                os.system('mv densityTemp %s'%out_folder)
            except:
                print "%s already exists, the results are in %s(NEW)"%(out_folder,out_folder)
                os.system('mv densityTemp %s'%(out_folder+'(NEW)'))
        else:
            os.chdir('..')
    
    def densityAccum(self,atomID=8,frame=100,equil=30000):
        d = dump(self.adsfile)
        times = d.time()
        times = [time for time in times if time>=equil]
        gap = times[1]-times[0]
        gap_need = (times[-1]-times[0])/frame
        if gap_need/gap >= 1:
            skip = gap_need/gap
        d.tselect.test("$t>=%d"%equil)
        d.delete()
        d.tselect.skip(skip)
        self.frame = len(d.time())
        d.aselect.all()
        d.aselect.test("$type == %d"%atomID)
        d.set("$type = 1")
        d.write('temp.ensembles.lammpstrj',0)
        f = open('temp.ensembles.lammpstrj','r')
        lines = f.readlines()
        lines_seq = [line.split() for line in lines]
        f.close()
        fp = open('ensembles.lammpstrj','w')
        fp.write('ITEM: TIMESTEP\n')
        fp.write('0\n')
        fp.write('ITEM: NUMBER OF ATOMS\n')
        fp.write('%d\n' % len(lines_seq))
        fp.write('ITEM: BOX BOUNDS ff ff ff\n')
        x = [float(item[2]) for item in lines_seq]
        y = [float(item[3]) for item in lines_seq]
        z = [float(item[4]) for item in lines_seq]
        xmax, xmin = np.max(x), np.min(x)
        ymax, ymin = np.max(y), np.min(y)
        zmax, zmin = np.max(z), np.min(z)
        fp.write('%.2f %.2f\n' % (xmin,xmax))
        fp.write('%.2f %.2f\n' % (ymin,ymax))
        fp.write('%.2f %.2f\n' % (zmin,zmax))
        fp.write('ITEM: ATOMS id type x y z \n')
        for i in range(len(lines_seq)):
            c = lines_seq[i][2:5]
            fp.write('%d %d %s %s %s\n' % (i+1,1,c[0],c[1],c[2]))
        fp.close()
        os.system('rm -f temp.ensembles.lammpstrj')
   
    def density1D_radial(self,radius,dumpfile='ensembles.lammpstrj',mesh=100):
        d = dump(dumpfile)
        t = d.time()[0]
        x,y,z = d.vecs(t,'x','y','z')
        coord = [[x[i],y[i],z[i]] for i in range(len(x))]
        r = np.linspace(0,radius,mesh,endpoint=False)
        coord_v = []
        for item in coord:
            v = [float(s) for s in item]
            rs = np.linalg.norm(v)
            coord_v.append(rs)
        hist = np.histogram(coord_v,r)
        self.radial_r = hist[1][:-1]
        self.radial_hist = hist[0]           
   
    def coordNum(self,\
            halflength,mesh,cut_off,dumpfile='ensembles.lammpstrj'):
        import lmp2py
        self.halflength = halflength
        self.mesh = mesh
        self.R_ave = cut_off
        coordnum_template = '%s/coordnum.lmpT'\
            %self.config['ads_script_path']
        d = dump(dumpfile)
        d.tselect.all()
        dump_time = d.time()
        cmd_list = lmp2py.lmp2list(coordnum_template)
        spacing = float(halflength)/float(mesh)
        for item in cmd_list:
            if item[0] == 'lattice':
                item[2] = str(spacing)
            if item[0] == 'region':
                item[3] = str(-halflength)
                item[4] = str(halflength)
                item[5] = str(-halflength)
                item[6] = str(halflength)
                item[7] = str(-halflength)
                item[8] = str(halflength)
            if item[0] == 'read_dump':
                item[2] = str(dump_time[0])
                item[1] =dumpfile
            if item[0] == 'compute':
                item[4] = str(cut_off)
        cmd_text = lmp2py.list2str(cmd_list)
        lmp2py.str2file(cmd_text,'coordnum')
        os.system(self.engine+'coordnum.lmp')
        os.system('mv log.lammps log.coordnum')
        
    
    def readCoordfile(self):
        self.hist_list = []
        d = dump(self.coordnum)
        d.tselect.all()
        hist_list = d.vecs(0,'c_DENSITY')
        hist_list = [h for h in hist_list if h!=0]
        self.hist_list = hist_list
    
    def evalHist(self,mesh=100):
        hist_list = self.hist_list
        frame = self.frame
        d_mesh = mesh
        d_max = frame
        #d_max = np.max(hist_list)
        d_range = np.linspace(0,d_max,d_mesh)
        hist = np.histogram(hist_list,d_range)
        self.data_x = hist[1][:-1]
        self.data_hist = hist[0]

    def V_rho_Spectrum(self,xlim):
        import matplotlib.pyplot as plt
        mesh_all = (2.0*self.mesh)**3
        self.data_hist = np.divide(self.data_hist,mesh_all)
        plt.rc('font', family='serif')
        plt.rc('legend',numpoints=1)
        plt.figure(figsize=(5,6))
        ax = []
        ax.append(plt.subplot(211))
        ax[0].set_position([0.2,0.6,0.7,0.38])
        ax.append(plt.subplot(212))
        ax[1].set_position([0.2,0.1,0.7,0.38])
        ax[1].plot(np.divide(self.data_x,self.frame),self.data_hist,'ko-',label="Kerogen Pore")
        plt.legend()
        plt.xlabel(r'$\rho_E$',fontsize=20)
        plt.ylabel(r'$V_F$',fontsize=20)
        plt.xlim(xmax=xlim)
        plt.show()
    
    def V_rho_SpectrumComparison(self,compData,refs,xlim=0.7):
        
#        f=open('data_regularPore_2.0','w')
#        f.write('rho_E,V_F')

        import matplotlib.pyplot as plt
        mesh_all = (2.0*self.mesh)**3
        self.data_hist = np.divide(self.data_hist,mesh_all)
        plt.rc('font', family='serif')
        plt.rc('legend',numpoints=1)
        plt.figure(figsize=(5,6))
        ax = []
        ax.append(plt.subplot(211))
        ax[0].set_position([0.2,0.6,0.7,0.38])
        for i in range(len(compData)):
            data, ref = compData[i], refs[i]
            mesh_data = (2.0*data.mesh)**3
            data.data_hist = np.divide(data.data_hist,mesh_data)
            ax[0].plot(np.divide(data.data_x,data.frame),\
                np.multiply(data.data_hist,1),'o-',\
                label=r"$\varepsilon_w = %s$"%ref)
                
#            for j in range(len(np.divide(data.data_x,data.frame))):
#                f.write('%.2f,%.6f\n'%(np.divide(data.data_x,data.frame)[j],data.data_hist[j]))

        plt.xlabel(r'$\rho_E$',fontsize=20)
        plt.ylabel(r'$V_F$',fontsize=20)       
        plt.xlim(xmax=xlim)
        plt.legend()
        ax.append(plt.subplot(212))
        ax[1].set_position([0.2,0.1,0.7,0.38])
#        ax[1].plot(np.divide(self.data_x,self.frame),self.data_hist,\
#            'ko-',label="Kerogen Pore")
        ax[1].plot(np.divide(self.data_x,self.frame),np.multiply(self.data_hist,1),\
            'ko-',label="Kerogen Pore")
        plt.legend()
        plt.xlabel(r'$\rho_E$',fontsize=20)
        plt.ylabel(r'$V_F$',fontsize=20)
        plt.xlim(xmax=xlim)
        ylim0 = ax[0].get_ylim()[1]
        ylim1 = ax[1].get_ylim()[1]
        ylim = max([ylim0,ylim1])
        ax[0].set_ylim((0,ylim))
        ax[1].set_ylim((0,ylim))
        plt.show()
#        f.close()
#        g=open('data_kerogenPore','w')
#        g.write('rho_E,V_F')
#        for j in range(len(np.divide(self.data_x,self.frame))):
#            g.write('%.2f,%.6f\n'%(np.divide(self.data_x,self.frame)[j],self.data_hist[j]))
#        g.close()  
        
    def visualDensity(self, vis = 'no', phi = 0.01):
        """
        visualDensity must be invoked after coordNum().
        The equatorial plane is the x-y plane.
        phi is the longitude angle.
        """
        d = dump(self.coordnum)
        d.aselect.all()
        time = d.time()[0]
        x, y, z, rho = d.vecs(time,"x","y","z","c_DENSITY")
        N = len(x)
        x_temp = sorted(x)
        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        zmin, zmax = np.min(z), np.max(z)
        for i in range(N):
            if x_temp[i]!=xmin:
                x_next = x_temp[i]
                break
        spacing =  x_next - xmin
        rho = np.divide(rho,self.frame*spacing**3)
        pixels = []
        pixel_x = range(int(round(xmin/spacing)),\
            int(round(xmax/spacing))+1)
        pixel_y = range(int(round(ymin/spacing)),\
            int(round(ymax/spacing))+1)
        pixel_z = range(int(round(zmin/spacing)),\
            int(round(zmax/spacing))+1)
        center = [pixel_x[len(pixel_x)/2],\
            pixel_y[len(pixel_y)/2],\
            pixel_z[len(pixel_z)/2]]
        for i in range(N):
            pixels.append([int(round(x[i]/spacing)),\
                int(round(y[i]/spacing)),\
                int(round(z[i]/spacing))])
        
        p_xmin, p_xmax = np.min(pixel_x), np.max(pixel_x)
        p_ymin, p_ymax = np.min(pixel_y), np.max(pixel_y)
        xs = center[0]-int((center[1]-p_ymin)/np.tan(phi))
        if xs < p_xmin:
            ys = center[1]-int((center[0]-p_xmin)*np.tan(phi))
            start = [p_xmin,ys]
            ye = center[1]-int((center[0]-p_xmax)*np.tan(phi))
            end = [p_xmax,ye]
        else:
            start = [xs,p_ymin]
            xe = center[0]-int((center[1]-p_ymax)/np.tan(phi))
            end = [xe,p_ymax]
        def Bresenham(start,end):
            dx = end[0]-start[0]
            dy = end[1]-start[1]
            e = -dx
            x, y = start[0],start[1]
            res = []
            for i in range(dx):
                res.append([x,y])
                x += 1
                e = e + 2*dy
                if e>=0:
                    y += 1
                    e = e - 2*dy
            return res
        plane_xy = Bresenham(start,end)
        slope = float(end[1]-start[1])/(end[0]-start[0])+0.0001
        #recover to length with units
        coord_start = np.subtract(start,center[:2])*spacing        
        def lineMap(x0,y0,coord_start,slope):
            k,ik = slope, 1.0/slope
            xs,ys = coord_start[0],coord_start[1]
            xp = (y0-ys+k*xs+ik*x0)/(k+ik)
            yp = ys + k*(xp-xs)
            return np.linalg.norm(np.subtract([xp,yp],coord_start))
        
        plane = []
        for i in range(N):
            pixel = pixels[i]
            if plane_xy.count(pixel[:2])>0:
                x0,y0 = pixel[0]*spacing,pixel[1]*spacing
                rxy = lineMap(x0,y0,coord_start,slope)
                z = pixel[2]*spacing
                plane.append([rxy,z,rho[i]])
        
        import matplotlib.pyplot as plt
        plt.rc('font', family='serif')
        plt.rc(('xtick','ytick'),labelsize=15)
        plt.figure(figsize=(8,6))
        plane = np.array(plane)
        rxy_max, z_max = tuple(np.max(plane,axis=0)[:2])
        rxy_min, z_min = tuple(np.min(plane,axis=0)[:2])
        grid_rxy,grid_z=np.mgrid[rxy_min:rxy_max:400j,z_min:z_max:400j]
        from scipy.interpolate import griddata
        grid_rho=griddata(plane[:,:2],plane[:,2],(grid_rxy,grid_z),method='cubic')
        if vis == 'no':
            return grid_rho
        else:
            plt.imshow(grid_rho.T,aspect='equal',\
                extent=(rxy_min,rxy_max,z_min,z_max),origin='lower')
            plt.show()
            

class Density_KeroPore(Density):
    def __init__(self,arg,pore_path):
        """
        pore_path: SOMEPATH/atoms_out
        tar_path: SOMEPATH/atoms_out/RDF
        The output files are all moved to tar_path when finalize is called
        """
        Density.__init__(self)
        self.para = arg
        datafile_o = 'data.9mol_du_mod'
        dumpfile_o = 'dump.atom.lammpstrj'
        sphere_center_file = 'sphere_center'
        os.system('cp -rf %s/%s %s/%s %s/%s .'%\
            (pore_path,datafile_o,pore_path,dumpfile_o,\
            pore_path,sphere_center_file))#Copy the origin 3 file to curdir
        self.method = 'KeroPore'
        self.tarpath = '%s/Density'%pore_path
        self.lmpfile = self.method+'.lmp'
        self.cmd_mod()
    
    def cmd_mod(self):
        fluidspath = '/home/sleepingz/2015Autumn/KeroPore/Adsorption'
        os.system('cp -rf %s/fluids .'%fluidspath)
        import IsoThermal as IT
        fluid = 'methane'
        fluidmass = 16.0
        ae = 0.29386
        asigma = 3.73
        from RDF import  dump_mod
        i = IT.GCMC_lj_ads(fluid=fluid,fluidmass=fluidmass,\
            T=self.para['T'],p=[self.para['p']],datafile='data.9mol_du_mod',\
            readdump='dump.atom_mod.lammpstrj',sphere_r=self.para['radius'],\
            equil=self.para['equil'],nsample=self.para['sample_num'],atom_energy=ae,atom_sigma=asigma)
        i.cmd_mod()
        self.para['mu'] = i.mu[0]#calculate mu according to p
        self.para['sphere_center'] = i.para['sphere_center']
        arg={'sphere_r':i.para['sphere_r'],'dump_step':500000}
        dump_mod(arg,'dump.atom.lammpstrj')#dump_mod will read the 'sphere_center' file
        
        for cmd in i.cmd_list:
            if cmd[:2] == ['fix','GCMC_E']:
                cmd[7] = str(i.para['fluid_id'])
                cmd[9] = str(i.T)
                cmd[10] = str(i.mu[0])
            if cmd[0] == 'dump':
                cmd[0] = '#dump'
        for cmd in i.cmd_list2:
            if cmd[:2] == ['fix','GCMC']:
                cmd[5] = '10'#Need to be improved
                cmd[7] = str(i.para['fluid_id'])
                cmd[9] = str(i.T)
                cmd[10] = str(i.mu[0])
        
        arg = self.para
        cmd_list = i.cmd_list + i.cmd_list2
        cmd_list.append(['dump','1','FLUID','custom',str(arg['freq']),'dump.ensembles','id','type','x','y','z'])
        cmd_list.append(['run',str(arg['sample_num']*arg['freq'])])
        import lmp2py
        lmp2py.str2file(lmp2py.list2str(cmd_list),self.method)

class Density_HardWall(Density):
    def __init__(self, arg, hw_path='.'):
        """
        savemode = 0: do not save the hardwall results
        pore_path: path to store the hardwall results, must tell the hw_path
        arg.keys = 'radius','sigma','E_W','T','mu','id','sample_num','equil','freq',...
        """
        Density.__init__(self)
        self.para = arg
        self.para['radius_L'] = self.para['radius'] + self.para['sigma']
        self.method = 'hardwall'
        self.lmpfile = self.method+'.lmp'
        self.lmptemplate =\
            '/home/sleepingz/2015Autumn/KeroPore/Case:hardwall/poremaker_hardwall.lmp'
        self.tarpath = hw_path
        self.cmd_mod()
            
        
    def cmd_mod(self):
        arg = self.para
        import lmp2py
        cmd_list = lmp2py.lmp2list(self.lmptemplate)
        for item in cmd_list:
            if item[0:2] == ['fix','WALL']:
                item[6] = str(arg['E_W'])
                item[7] = str(arg['sigma'])
                item[8] = str(arg['sigma']*2.5)
            if item[0:2] == ['fix','GCMC']:
                item[9] = str(arg['T'])
                item[10] = str(arg['mu'])
            if item[0:2] == ['region','SPHERE']:
                item[6] = str(arg['radius'])#between radius and radius_L to capture all adsoption structure
            if item[0:2] == ['region','SPHERE_L']:
                item[6] = str(arg['radius']+arg['sigma'])
            if item[0] == 'dump':
                item[0] = '#dump'
            if item[0] == 'run':
                item[1] = str(arg['equil'])
            if item[0] == 'create_atoms':
                item[3] = str(arg['N_iso'])
        cmd_list.append(['dump','1','all','custom',str(arg['freq']),'dump.ensembles','id','type','x','y','z'])
        cmd_list.append(['unfix','GCMC'])
        cmd_list.append(['run',str(arg['sample_num']*arg['freq'])])
        lmp2py.str2file(lmp2py.list2str(cmd_list),self.method)

class Density_StaticAtoms(Density):
    def __init__(self, arg, sa_path='.'):
        """
        savemode = 0: do not save the hardwall results
        pore_path: path to store the hardwall results, must tell the hw_path
        arg.keys = 'radius','sigma','E_W','N_iso','T','id','sample_num','equil','freq',...
        """
        Density.__init__(self)
        self.para = arg
        self.para['radius_L'] = self.para['radius'] + self.para['sigma']
        self.method = 'StaticAtoms'
        self.lmpfile = self.method+'.lmp'
        self.lmptemplate =\
            '/home/sleepingz/2015Autumn/KeroPore/Case:hardwall/poremaker_staticAtoms.lmp'
        self.tarpath = sa_path
        self.cmd_mod()
            
        
    def cmd_mod(self):
        arg = self.para
        import lmp2py
        cmd_list = lmp2py.lmp2list(self.lmptemplate)
        for item in cmd_list:
            if item[0:2] == ['region','BOX']:
                item[3] = item[5] = item[7] = str(-arg['radius']-10)
                item[4] = item[6] = item[8] = str(arg['radius']+10)
            if item[0:2] == ['region','SPHERE']:
                item[6] = str(arg['radius'])#between radius and radius_L to capture all adsoption structure
            if item[0:2] == ['region','SPHERE_L']:
                item[6] = str(arg['radius']+arg['sigma'])
            if item[0:2] == ['pair_coeff','1']:
                item[3] = str(arg['E_W'])
                item[4] = str(arg['sigma'])
            if item[0:2] == ['fix','EQUIL']:
                item[5] = item[6] = str(arg['T'])
            if item[0] == 'run':
                item[1] = str(arg['equil'])
            if item[0:2] == ['create_atoms','2']:
                item[3] = str(arg['N_iso'])
        cmd_list.append(['dump','1','all','custom',str(arg['freq']),'dump.ensembles','id','type','x','y','z'])
        cmd_list.append(['run',str(arg['sample_num']*arg['freq'])])
        lmp2py.str2file(lmp2py.list2str(cmd_list),self.method)    
#Test:
#arg={'radius':7.0,'sigma':3.7,'E_W':1.0,'T':300.0,'mu':-6.0,'id':1,'sample_num':200,\
#    'equil':10000,'freq':200,'mesh':20,'r_ave':2.0,'sphere_center':[0.0,0.0,0.0]}
#Test1: hardwall
#dh = Density_HardWall(arg)
#dh()
#dh.cmd_mod()
#dh.ensembleAccum()
#dh.density_eval()

#Test2: 3 hardwall, plot
#mu_seq = [-8.0,-7.0,-6.0]
#data = []
#for i in range(len(mu_seq)):
#    arg['mu'] = mu_seq[i]
#    dh = Density_HardWall(arg)
#    dh.cmd_mod()
#    dh()
#    dh.ensembleAccum()
#    dh.density_eval()
#    dh.read_coordfile()
#    dh.eval_hist()
#    data.append({'x':dh.data_x,'hist':dh.data_hist})
#    os.chdir('..')
#
#
#layout = 311
#
#import matplotlib.pyplot as plt
#fig = plt.figure(figsize=(8,10))
#ax = []    
#for i in range(3):
#    ax.append(plt.subplot(layout+i))
#    ax[i].plot(data[i]['x'],data[i]['hist'],'bo-')

#Test3: KeroPore prepare
#pore_path = '/home/sleepingz/2015Autumn/KeroPore/Case:dig_atoms/radius_7.0/atoms_out'
#arg['p'] = 1.0
#dk = Density_KeroPore(arg,pore_path)
#dk.cmd_mod()
#dk.ensembleAccum(mode = 1)
#dk.density_eval()
#dk.read_coordfile()
#dk.eval_hist()
#import matplotlib.pyplot as plt
#fig = plt.figure(figsize=(8,10))
#plt.plot(dk.data_x,dk.data_hist,'bo-')
