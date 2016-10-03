# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 12:55:54 2016

@author: sleepingz
"""

from dump import dump
import lmp2py
import numpy as np
import os

sat_default = {'fluid':'methane','fluidmass':16.0,\
    'T':300.0,'mu':-5.0,\
    'dump':'dump.lammpstrj',\
    'data':'data.dreiding','step':500000,\
    'stepRun':20000}

def FinalStep(dumpfile):
    d = dump(dumpfile)
    time_final  = d.time()[-1]
    return time_final

def FinalFrame(dumpfile):
    d = dump(dumpfile)
    time_final = d.time()[-1]
    d.tselect.one(time_final)
    d.write('ff_'+dumpfile)
    
def PoreRefine(dumpfile, poreRadius, \
    poreThickness, global_cutoff = 15.0, cutter_id = 7):
    FinalFrame(dumpfile)
    ff_dump = 'ff_'+dumpfile
    d = dump(ff_dump)
    #x,y,z = d.vecs(time,'')
    snap = d.snaps[0]
    xprd = snap.xhi - snap.xlo
    yprd = snap.yhi - snap.ylo
    zprd = snap.zhi - snap.zlo
    snap.xlo, snap.xhi = -xprd/2.0, xprd/2.0
    snap.ylo, snap.yhi = -yprd/2.0, yprd/2.0
    snap.zlo, snap.zhi = -zprd/2.0, zprd/2.0
    atoms = snap.atoms
    x = d.names["x"]#index of the list of a single line in the dump file
    y = d.names["y"]
    z = d.names["z"]
    ix, iy, iz = [], [], []
    for atom in atoms:
        if atom[x] < snap.xlo:
            ix.append(-1)
        elif snap.xlo >= snap.xhi:
            ix.append(1)
        else:
            ix.append(0)
    for atom in atoms:
        if atom[y] < snap.ylo:
            iy.append(-1)
        elif snap.ylo >= snap.yhi:
            iy.append(1)
        else:
            iy.append(0)
    for atom in atoms:
        if atom[z] < snap.zlo:
            iz.append(-1)
        elif snap.zlo >= snap.zhi:
            iz.append(1)
        else:
            iz.append(0)
    atoms[:,x] -= np.array(ix)*xprd
    atoms[:,y] -= np.array(iy)*yprd
    atoms[:,z] -= np.array(iz)*zprd
    step_now = d.time()[0]
    d.aselect.test("$x**2+$y**2+$z**2 < %f and $type == %d"\
        %((poreRadius + poreThickness)**2,cutter_id),step_now)
    d.set("$type = %d"%(cutter_id + 1))
    d.aselect.all(step_now)
    #cut the sphere with r = r_p + global_cutoff
    d.aselect.test("$x**2+$y**2+$z**2 < %f"\
        %(poreRadius + global_cutoff)**2,step_now)
    
    d.write('refine_'+dumpfile)
    
class adsorption:
    def __init__(self,config,ads_dict=sat_default,method='saturation'):
        self.config = config
        self.dict = ads_dict
        self.template = 'ads_'+method
        os.system('cp -f %s/%s .'%(config['ads_script_path'],self.template+'.lmpT'))

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
 
class ads_saturation(adsorption):
     def __init__(self,config,ads_dict=sat_default):
         adsorption.__init__(self,config,ads_dict=ads_dict,method = 'saturation')
         self.script_mod()
         
#Test:
if __name__ == '__main__':
    PoreRefine("ads_pre.lammpstrj",30.0)