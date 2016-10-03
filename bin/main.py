# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 14:12:22 2016

@author: sleepingz
"""
import sys

class Main:
    def __init__(self):
        self.configure()
    
    def configure(self):
        self.config = {}
        f = open('config','r')
        lines = f.readlines()
        f.close()
        for line in lines:
            line_seq = line.strip().split(':')
            self.config[line_seq[0]] = line_seq[1]
	sys.path.append(self.config['keros_path'])
        sys.path.append(self.config['ads_script_path'])
        sys.path.append(self.config['recon_script_path'])
    
    def Project(self,proj_name):
        import os
        path = self.config['workspace'] + '/' + proj_name
        try:
            os.system('mkdir -p %s'%path)
        except:
            pass
        os.chdir(path)
