import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os

from sys import argv

global n_ions
############################################################
runmin=0
runmax=53


common_parent_path='/data_prace/michele/prace/cdc'
base_dir_script=common_parent_path+'/analysis/results'

paths_to_launch=[
'/cdc1200/nacl/8800_160/0.0V', 
#'/cdc1200/nacl/8800_160/1.0V'#, 
'/cdc1200/nacl/8800_80/0.0V',
'/cdc1200/nacl/8800_80/1.0V',
'/cdc800/nacl/7615_139/0.0V',
'/cdc800/nacl/7615_139/1.0V',
'/cdc800/nacl/7700_70/0.0V',
'/cdc800/nacl/7700_70/1.0V'
]

for target_path in paths_to_launch:
	dir_to_go=base_dir_script+target_path+'/2_convert_position2VMD_results'
	os.chdir(dir_to_go)
	file_1=open('traj_000.xyz','r')
	numatoms=int(file_1.readline())
	file_1.close()	
	filenames=[]
	for irun in range(runmin,runmax+1):
		filenames.append('traj_%03d.xyz'%(irun))	

	num_lin=-1
	num_fram=1
	with open('all_traj_%03d_%03d.xyz'%(runmin,runmax), 'w') as outfile:
	    for fname in filenames:
	        with open(fname) as infile:
	            for line in infile:
	            	if (num_lin%(numatoms+2)==0):
	            		line='pas\t%d\n'%(num_fram)
	            		num_fram+=1
			num_lin+=1
	                outfile.write(line)	

	outfile.close()
