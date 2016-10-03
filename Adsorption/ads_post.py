# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 12:06:38 2016

@author: sleepingz
"""
import numpy as np
from scipy.optimize import leastsq

def nmolAnalysis(avefile,equil_step):
    f = open(avefile,'r')
    lines = f.readlines()[2:]
    lines_equil = []    
    for line in lines:
        seq = line.strip().split()
        seq = [int(seq[0]),float(seq[1])]
        if seq[0] >= equil_step:
            lines_equil.append(seq)
    f.close()
    nmol_seq = [line[1] for line in lines_equil]
    N_ave = np.average(nmol_seq)
    N_std = np.std(nmol_seq)
    return N_ave, N_std

def adsGroup(path,group_name):
    import os
    dirs = os.listdir(path)
    iso_dirs = [name for name in dirs if name.startswith('iso:')]
    try:
      os.system('mkdir %s'%group_name)
    except:
      pass
    for iso_dir in iso_dirs:
      try:
	os.system('mv %s %s'%(iso_dir,group_name))
      except:
	pass
    os.chdir(group_name)

def writeAdsFile(config,fluid,path,equil_step = 10000):
    import os
    dirs = os.listdir(path)
    iso_dirs = [name for name in dirs if name.startswith('iso:')]
    iso_dirs = sorted(iso_dirs,key=lambda x:float(x.split('_')[4]))
    from ads import fluid_info,P2mu
    molMass = fluid_info[fluid]['molMass']
    P2mu_ins = P2mu(config)
    f = open('adsorption.data','w')
    f.write("T(K)    p(MPa)    mu(kCal/mol)    Nave(1)    Nstd(1)\n")
    for iso_path in iso_dirs:
        name = iso_path.split('_')
        N_ave, N_std = nmolAnalysis('%s/%s/nmol.ave'%(path,iso_path),equil_step)
        T = float(name[2])
        P = float(name[4])
        mu = P2mu_ins.conv(fluid,molMass,T,P)
        f.write("%8.5f,%8.5f,%8.5f,%8.5f,%8.5f\n"%\
            (T,P,mu,N_ave,N_std))
    f.close()

def LangmuirFit(p,Nad):
    def residuals(t,q,p):
        pL,qL=t
        err=q-qL*np.array(p)/(pL+np.array(p))
        return err
    t0=[Nad[-1]/Nad[1],Nad[-1]]
    tlsq=leastsq(residuals,t0,args=(Nad,p))[0]
    #Nad_fit=tlsq[1]*np.array(p)/(tlsq[0]+np.array(p))
    pL=tlsq[0]
    qL=tlsq[1]
    return pL,qL
    
def readAdsFile(filename):
    f=open(filename)
    res = []
    f.readline()#jump the 1st line
    while True:
        line=f.readline()
        if (line==''):
            break
        line=line.split(",")
        res.append({'T':float(line[0]),'p':float(line[1]),'mu':float(line[2]),\
            'N':[float(line[3]),float(line[4])]})
    f.close()
    return res
    
class adsData():
    def __init__(self,files):
        self.keys = files.keys()
        self.keys.sort(key = lambda x : float(x))
        self.data,self.p,self.N,self.N_err={},{},{},{}
        self.pL,self.qL={},{}
        for key in self.keys:
            self.data[key] = readAdsFile(files[key])
            self.p[key]=[item['p'] for item in self.data[key]]
            self.N[key]=[item['N'][0] for item in self.data[key]]
            self.N_err[key]=[item['N'][1] for item in self.data[key]]
            self.pL[key],self.qL[key]=LangmuirFit(self.p[key],self.N[key])

def singleIsoPlot(data):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(6,5))
    plt.rc('font', family='serif')
    plt.rc('legend',numpoints=1)
    colors = ['b','g','r','c','m','k']
    ax=plt.subplot(111)
    ax.set_position([0.13,0.16,0.54,0.7])
    
    p_sample = np.linspace(0.0,20.0,100)
    i = 0
    for key in data.keys:
        ax.plot(data.p[key],data.N[key],colors[i]+'o',label = r'$\mathrm{CH}_4, r_p=%s\mathrm{\AA}$'%(float(key)+3.7))
        ax.errorbar(data.p[key],data.N[key],yerr = data.N_err[key],ls = 'None',ecolor=colors[i])
        Nad_fit=data.qL[key]*p_sample/(data.pL[key]+p_sample)
        ax.plot(p_sample,Nad_fit,colors[i]+'-')
        i += 1
    plt.xlim(0,20)
    plt.xlabel(r"$p_{\mathrm{reservior}}\mathrm{(MPa)}$",fontsize=20)
    plt.ylabel(r"$N_{\mathrm{CH_4}}}$",fontsize=20)
    plt.legend(loc='lower right',bbox_to_anchor=(1.6, 0.1),fontsize=14,numpoints=1,handletextpad=-0.5)
    plt.show()    
