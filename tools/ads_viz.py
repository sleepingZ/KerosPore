# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 19:43:55 2015

@author: sleepingz
"""
import os
import sys
from dump import dump
from vmd import vmd

path = '/home/sleepingz/2015Autumn/KeroPore/Case:dig_atoms/radius_7.0/GCMC_hybrid_methane_300.0K/methane_300.000000_-6.269508'
os.chdir(path)
fp = open('../../atoms_out/sphere_center','r')
sc = str(fp.readline()).split(',')
sc = [float(item) for item in sc]
radius = float(path.split('/')[6].split('_')[1])
sigma = 3.72
radius += sigma

d = dump('dump.ads.lammpstrj')
d.aselect.all()
t = d.time()
d.tselect.one(t[-1])
d.aselect.test(" ($x -%f)**2 + ($y -%f)**2 + ($z -%f)**2 < %f"\
        %(sc[0],sc[1],sc[2],radius**2))
d.aselect.test("$type == 7 or $type==8")
d.set("$type = 10")
d.tselect.one(t[-1])
d.write('dump.ads_mod.lammpstrj')

v = vmd()
v.rep('Line')
v.new('dump.ads_mod.lammpstrj')
v('mol addrep 0')
v('mol modselect 1 0 name 10')
v('mol modselect 0 0 all not name 10')
v('mol modstyle 1 0 VDW')
v('mol modstyle 0 0 DynamicBonds')
#v.stop()