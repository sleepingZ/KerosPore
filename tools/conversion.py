import numpy as np
kcal=4.186#kJ
k_B=0.0019872067
h_P=95.306976368
mvv2e=48.88821291*48.88821291
#gas_mass=16.0
#R=8314.0/gas_mass#J/(K.kg)
N_A=6.02e23
def Lambda(T,gas_mass):#(m)
    return np.sqrt(h_P*h_P/(2*np.pi*gas_mass*mvv2e*k_B*T))*1e-10
def zz(mu,T):
    return np.exp(mu/(k_B*T))
def Ndensity(rho,gas_mass):
    return rho/gas_mass*1e6*N_A#(1/m^3)
def mu_id(rho,T,gas_mass):
    return k_B*T*np.log(Lambda(T,gas_mass)**3*Ndensity(rho,gas_mass))#kcal/mole
def mu_idp(press,T,gas_mass):#press:Pa
    return k_B*T*np.log(Lambda(T,gas_mass)**3*press/1.38e-23/T)
