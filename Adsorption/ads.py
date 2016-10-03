# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:24:52 2016

@author: sleepingz
"""
ads_default = {'fluid':'methane','atomMass':16.0,\
    'T':300.0,'mu':-5.0,\
    'dump':'dump.lammpstrj',\
    'data':'data.dreiding',\
    'poreThickness':4.0,\
    'atomEpsilon':0.29386,'atomSigma':3.723,\
    'ensembleFrame':1,'stepRun':30000}#pressure unit: MPa

fluid_info = {'methane':{'molMass':16.0,'atomEpsilon':0.29386,'atomSigma':3.723}}
from ctypes import *

class P2mu:
    def __init__(self,config):
        self.config = config
        fluidProp = '%s/libFugacityCoef.so'%config['nist_path']
        RefProp = CDLL(fluidProp)
        self.FugacityCoef = RefProp.FugacityCoef_
        self.FugacityCoef.argtypes = [c_double,c_double,POINTER(c_char)]
        self.FugacityCoef.restype = c_double   
        
    def conv(self,fluid,molMass,T,P_real):
        import conversion as conv
        import os
        import numpy as np
        origin = os.getcwd()
        os.chdir(self.config['nist_path'])
        Lambda=conv.Lambda(T,molMass)
        fluidFile = '%s.fld' % fluid
        P_ideal = self.FugacityCoef(T,P_real*1000,fluidFile)*P_real#P(MPa)
        k_B=1.3806488e-23
        N_ideal=P_ideal*1e6/(k_B*T)
        mu=conv.k_B*T*np.log(Lambda**3*N_ideal)
        os.chdir(origin)
        return mu

from ads_pre import adsorption

class ads_lj(adsorption):
     def __init__(self,config,ads_dict):
         adsorption.__init__(self,config,ads_dict=ads_dict,method = 'lj')
         self.script_mod()

class ads_lj93Wall(adsorption):
    def __init__(self,config,lj93_dict):
        adsorption.__init__(self,config,ads_dict=lj93_dict,method = 'lj93Wall')
        self.script_mod()
    def __call__(self):
        import os
        os.system('%s<ads_lj93Wall.lmp'%(self.config['engine']))
        