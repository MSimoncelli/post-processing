import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os
from os import listdir
from os.path import isfile, join

from sys import argv

global n_ions
############################################################
if (len(argv)<=1):
        print 'USAGE: 5_plot_condictivity.py <datafile>'
        
else:
        root_name=argv[1]

from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################
conv=0.529177
max_dist=186.406 #length of the whole system
work_dir=os.getcwd()

file_list = [f for f in listdir(work_dir) if isfile(join(work_dir, f))]

data_time={}
data_msqxy={}
for files in file_list:
        if (root_name in files):
                data_time[files]=np.loadtxt(files)[:,0]
                data_msqxy[files]=np.loadtxt(files)[:,1]



name_file_save=r'5_Kubo_msd%s.pdf' %(root_name)
target_path=root_name.replace('_','/')
fig=pl.figure()
ax1=fig.add_subplot(111)

title='Mean square displacement %s' %(root_name)
title=title.replace('_', ' ')
title=title.replace('nacl', 'NaCl')
title=title.replace('kcl', 'KCl')
title=title.replace('0 0V', '0.0V')
title=title.replace('1 0V', '1.0V')
ax1.set_title(title)
#print title
ax1.grid(True)
ax1.set_xlabel(r'time $[ps]$')
ax1.set_ylabel(r'msd [\AA^2]')
list_keys=sorted(data_msqxy.keys())

for keys in list_keys:
	if keys.split('_')[-3]=='spc2':
		ion='Cl'
	else:
		if 'nacl' in root_name:
			ion='Na'
		else:
			ion='K'
	if keys.split('_')[-1][:-4]=='1':
		nam_reg=' el. 1'
	elif keys.split('_')[-1][:-4]=='2':
		nam_reg=' bulk'
	elif keys.split('_')[-1][:-4]=='3':
		nam_reg=' el. 2'
	string_name=ion+nam_reg
	ax1.plot(data_time[keys],data_msqxy[keys]*(conv*conv), label=string_name)
     

ax1.legend(loc='best')
     
fig.savefig(name_file_save, dpi=300)

