# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:06:42 2016

@author: sleepingz
"""
import lmp2py
import os

class Recon:
    def __init__(self,data,config,method = 'NPT'):
        self.config = config       
        self.data = data
        self.template = 'recon_'+method
        os.system('cp -f %s/%s .'%(config['recon_script_path'],self.template+'.lmpT'))        
        self.dict = {}    
        
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
    def __init__(self,data,config,radius,radius_L,T=300.0,P=200.0):
        cutter = Recon_cutter(config,radius,radius_L,T=T,data=data)   
        Recon.__init__(self,data,config,method = 'NPT_cutter')
        self.kwds = ['data','T','P']
        self.values = [self.data,T,P]
        self.dict = dict(zip(self.kwds,self.values))
    
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
	

    
        
        

        
        
