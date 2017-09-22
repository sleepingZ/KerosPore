# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:06:42 2016

@author: sleepingz
"""
import lmp2py
import os


class Recon:
    def __init__(self, data, config, method='NPT'):
        self.config = config       
        self.data = data
        self.template = 'recon_'+method
        os.system('cp -f %s/%s .' % (config['recon_script_path'], self.template+'.lmpT'))
        self.dict = {}
        self.fluid_type = 1
        
    def script_mod(self):
        cmd_list = lmp2py.lmp2list(self.template+'.lmpT')
        self.cmd_str = lmp2py.list2str(cmd_list)
        self.cmd = []
        for cmd in self.cmd_str:
            for key in self.dict.keys():
                cmd = cmd.replace('$%s$'%key,'%s'%str(self.dict[key]))        
            self.cmd.append(cmd)
        lmp2py.str2file(self.cmd,self.template)
        os.system('rm -f %s'%self.template+'.lmpT')

    def solid_atom_input(self, file_name):
        f = open(self.config['pore_data_path'] + '/%s' % file_name)
        lines = f.readlines()[1:]
        self.solid_atom_potentials = ''.join(lines)
        self.fluid_type = len(lines) + 1


class Recon_NPT(Recon):
    def __init__(self,data,config,T=300.0,P=200.0):#T:K P: bar
        Recon.__init__(self,data,config,method = 'NPT')
        self.kwds = ['data','T','P']
        self.values = [self.data,T,P]
        self.dict = dict(zip(self.kwds,self.values))
        self.script_mod()


class Recon_cutter(Recon):
    def __init__(self,config,radius,radius_L,T=300.0,data=''):
        Recon.__init__(self,data,config,method = 'cutter')
        self.kwds = ['radius','radius_L','T']
        self.values = [radius,radius_L,T]
        self.dict = dict(zip(self.kwds,self.values))
        self.script_mod()
        engine = self.config['engine']
        os.system('%s<%s'%(engine,self.template+'.lmp'))
        os.system('mv log.lammps cutter.log')


class Recon_NPT_cutter(Recon):
    def __init__(self, data, config, radius, radius_L, T=300.0, P=200.0, fluid='methane'):
        cutter = Recon_cutter(config, radius, radius_L, T=T, data=data)
        Recon.__init__(self, data, config, method='NPT_cutter')
        import ads.fluid_info as fluid_info
        fluid_data = fluid_info[fluid]
        self.solid_atom_input('recon_NPT_cutter.lmpT')
        self.kwds = ['data', 'T', 'P', 'cutter_atom_ID', 'cutter_atom_mass', 'cutter_atom_epsl', 'cutter_atom_sig']
        self.values = [self.data, T, P, self.fluid_type, fluid_data['molMass'], fluid_data['atomEpsilon'],
                       fluid_data['atomSigma']]
        self.dict = dict(zip(self.kwds, self.values))
    
    def cutter_id_mod(self,id_max):
        dump = 'cutter_atoms.lammpstrj'
        atom_num = id_max
        g = open(dump,'r')
        lines = g.readlines()
        g.close()
        g = open(dump,'w')
        atoms_lineno = 10
        index = atoms_lineno - 1
        for line in lines[atoms_lineno-1:]:
            atom_num += 1
            line_seq = line.strip().split()
            line_seq[0] = str(atom_num)
            lines[index] = ' '.join(line_seq) + ' \n'
            index += 1

        g.writelines(lines)
        g.close()
