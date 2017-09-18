# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 19:05:42 2016

@author: sleepingz
"""

def mol_template(m,template,bond_seq=[1],angle_seq=[1]):
    path = m.config['nist_path']+'/Templates'
    f = open(path+'/'+template,'r')
    lines = f.readlines()
    f.close()
    output = []
    bond_len = len(bond_seq)
    angle_len = len(angle_seq)
    keys = ['Bond_%d'%(i+1) for i in range(bond_len)] + \
        ['Angle_%d'%(i+1) for i in range(angle_len)]
    values = bond_seq + angle_seq
    for line in lines:
        for i in range(len(keys)):
            line = line.replace('$%s$'%keys[i],'%d'%values[i])
        output.append(line)
    g = open(template,'w')
    g.writelines(output)
    g.close()
    