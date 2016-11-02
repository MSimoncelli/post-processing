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
data_md_x={}
data_md_y={}

for files in file_list:
        if (root_name in files):
                data_time[files]=np.loadtxt(files)[:,0]
                #take the squareroot of each displacement. Look at the equations to understand this!
                data_md_x[files]=np.sqrt(np.loadtxt(files)[:,1])
                data_md_y[files]=np.sqrt(np.loadtxt(files)[:,2])
                length=len(data_md_y[files])



name_file_save=r'5_Kubo_cond%s.pdf' %(root_name)
target_path=root_name.replace('_','/')
fig=pl.figure()
ax1=fig.add_subplot(111)


array_reg_1_x=np.zeros(length)
array_reg_2_x=np.zeros(length)
array_reg_3_x=np.zeros(length)

array_reg_1_y=np.zeros(length)
array_reg_2_y=np.zeros(length)
array_reg_3_y=np.zeros(length)

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
list_keys=sorted(data_md_x.keys())

reg_1=0
reg_2=0
reg_3=0

for keys in list_keys:
	time=data_time[keys]
	if keys.split('_')[-3]=='spc2':
		ion='Cl'
	else:
		if 'nacl' in root_name:
			ion='Na'
		else:
			ion='K'
	if keys.split('_')[-1][:-4]=='1':
		if(reg_1<1):
			array_reg_1_x=np.add(array_reg_1_x, data_md_x[keys])
			array_reg_1_y=np.add(array_reg_1_y, data_md_y[keys])
		else:
			array_reg_1_x=np.add(array_reg_1_x, -data_md_x[keys])
			array_reg_1_y=np.add(array_reg_1_y, -data_md_y[keys])
		reg_1=reg_1+10


	elif keys.split('_')[-1][:-4]=='2':
		if(reg_2<1):
			array_reg_2_x=np.add(array_reg_2_x, data_md_x[keys])
			array_reg_2_y=np.add(array_reg_2_y, data_md_y[keys])
		else:
			array_reg_2_x=np.add(array_reg_2_x, -data_md_x[keys])
			array_reg_2_y=np.add(array_reg_2_y, -data_md_y[keys])
		reg_2=reg_2+10

	elif keys.split('_')[-1][:-4]=='3':
		if(reg_3<1):
			array_reg_3_x=np.add(array_reg_3_x, data_md_x[keys])
			array_reg_3_y=np.add(array_reg_3_y, data_md_y[keys])
		else:
			array_reg_3_x=np.add(array_reg_3_x, -data_md_x[keys])
			array_reg_3_y=np.add(array_reg_3_y, -data_md_y[keys])
		reg_3=reg_3+10


print root_name
if ('8800_80' in root_name): 
        N_avg_el1=5
        N_avg_el2=5
        N_avg_bulk=70
if ('7700_70' in root_name): 
        N_avg_el1=4.5
        N_avg_el2=4.5
        N_avg_bulk=61
if ('8800_160' in root_name):
        N_avg_el1=10
        N_avg_el2=10
        N_avg_bulk=140
if ('7615_139' in root_name): 
        N_avg_el1=8.5
        N_avg_el2=8.5
        N_avg_bulk=122

ax1.plot(time,np.add(np.square(array_reg_1_x),N_avg_el1*np.square(array_reg_1_y))*(conv*conv), label='el. 1')
ax1.plot(time,np.add(np.square(array_reg_2_x),N_avg_bulk*np.square(array_reg_2_y))*(conv*conv), label='bulk')
ax1.plot(time,np.add(np.square(array_reg_3_x),N_avg_el2*np.square(array_reg_3_y))*(conv*conv), label='el. 2') 

ax1.legend(loc='best')
     
fig.savefig(name_file_save, dpi=300)

