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
    'ensembleFrame':1,'stepRun':10000,\
    'N_exchange':10,'N_move':10}#pressure unit: MPa

fluid_info = {'methane':{'molMass':16.0,'atomEpsilon':0.29386,'atomSigma':3.723},\
    'co2':{'molMass':44.0,'atomEpsilon':0.46879,'atomSigma':3.720}}
from ctypes import *

class P2mu:
    def __init__(self,config):
        self.config = config
        fluidProp = '%s/libFugacityCoef.so'%config['nist_path']
        RefProp = CDLL(fluidProp)
        self.FugacityCoef = RefProp.FugacityCoef_
        self.FugacityCoef.argtypes = [c_double,c_double,POINTER(c_char)]
        self.FugacityCoef.restype = c_double   
        self.biMixFugacity = RefProp.biMix_Fugacity_
        self.biMixFugacity.argtypes = [POINTER(c_char),POINTER(c_char),\
            c_double,c_double,c_double,POINTER(c_double)]
        
    def conv(self,fluid,molMass,T,P_real):
        import conversion as conv
        import os
        import numpy as np
        origin = os.getcwd()
        os.chdir(self.config['nist_path'])
        Lambda=conv.Lambda(T,molMass)
        fluidFile = '%s.fld' % fluid
        P_ideal = self.FugacityCoef(T,P_real*1000,fluidFile)*P_real#P_real(MPa)
        k_B=1.3806488e-23
        N_ideal=P_ideal*1e6/(k_B*T)
        mu=conv.k_B*T*np.log(Lambda**3*N_ideal)
        os.chdir(origin)
        return mu
    
    def F2mu(self,molMass,T,f):#f (P_ideal) in kPa
        import conversion as conv
        import numpy as np
        Lambda=conv.Lambda(T,molMass)
        k_B=1.3806488e-23
        N_ideal=f*1e3/(k_B*T)
        mu=conv.k_B*T*np.log(Lambda**3*N_ideal)
        return mu
        
    def biMix_f(self,fluid1,fluid2,T,P,x1):
        import os
        origin = os.getcwd()
        os.chdir(self.config['nist_path'])
        fluidFile1 = '%s.fld' % fluid1
        fluidFile2 = '%s.fld' % fluid2
        res = (c_double*2)()
        self.biMixFugacity(fluidFile1,fluidFile2,T,P*1000,x1,res)
        os.chdir(origin)
        return [res[0],res[1]]

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

class ads_biM_lj(adsorption):
    '''
        biM = binary Mixture
    '''
    def __init__(self,config,ads_ljbiM_dict):
        adsorption.__init__(self,config,ads_dict = ads_ljbiM_dict,method = 'ljbiM')
        self.script_mod()

class ads_biM_mol(adsorption):
    '''
        biM = binary Mixture
        mol = flexible molecule model
    '''
    def __init__(self,config,ads_molbiM_dict):
        adsorption.__init__(self,config,ads_dict = ads_molbiM_dict,method = 'molbiM')
        self.script_mod()