import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os
import scipy as sp
from scipy import signal

from sys import argv
############################################################
if (len(argv)<=1):
        print 'USAGE: plor_coord_num.py runmin runmax'
        runmin = 0
        runmax = 53
else:
	runmin=int(argv[1])
	runmax=int(argv[2])

name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_dir=name_folder.split('/')
#print list_dir

name_R_charge='Rwallcharge-all_%03d_%03d.dat'%(runmin,runmax)
name_L_charge='Lwallcharge-all_%03d_%03d.dat'%(runmin,runmax)
name_kintemp='kintemp-all_%03d_%03d.dat'%(runmin,runmax)

from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################
data_R_charge=np.loadtxt(name_R_charge)
data_L_charge=np.loadtxt(name_L_charge)
data_kintemp=np.loadtxt(name_kintemp)

time_step=10**(-3)#=1 fs= 10^(-6) nanosecond. The charge is recorded every 1000 step=every ps!
tot_num_data=len(data_R_charge)


time_array=np.arange(0,int(tot_num_data))
time_array=(time_array*time_step)

fig=pl.figure()
ax1 = fig.add_subplot(111)
title= 'Charge vs time run%03d-run%03d %s %s %s %s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11] )
title=title.replace('_', '-')
print title, 'TOT_NUM_DATA=',tot_num_data
ax1.set_title(title)
ax1.grid(True)
# ax1.set_xlim([-50,50])
ax1.set_xlabel(r'time $[ns]$')
ax1.set_ylabel(r'charge [e]')

ax1.plot(time_array,data_R_charge, label='right electrode' ,color='r')	
ax1.plot(time_array,data_L_charge, label='left electrode' ,color='b')	
ax1.legend(loc='best')

#pl.show()
name_file_save='5_Charge_run%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
name_file_save=name_file_save.replace('.','_')
name_file_save+='.pdf'

fig.savefig(name_file_save, dpi=300)

##################################################################################
time_step=10**(-3)#=1 fs= 10^(-6) nanosecond. The charge is recorded every 1000 step=every ps!
tot_num_data_2=len(data_kintemp)

time_array=np.arange(0,int(tot_num_data_2))
time_array=(time_array*time_step)

filtered = sp.signal.medfilt(data_kintemp[:,1],21)

fig=pl.figure()
ax1 = fig.add_subplot(111)
title= 'Temperature vs time run%03d-run%03d %s %s %s %s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11] )
title=title.replace('_', '-')
print title, 'TOT_NUM_DATA_2=',tot_num_data_2
print 'DIFF=', tot_num_data- tot_num_data_2
print ''
ax1.set_title(title)
ax1.grid(True)
# ax1.set_xlim([-50,50])
ax1.set_xlabel(r'time $[ns]$')
ax1.set_ylabel(r'T [K]')

ax1.plot(time_array,data_kintemp[:,1],label='raw data' ,color='r')	
ax1.plot(time_array,filtered,label='median filtered (21) ' ,color='b')		
ax1.legend(loc='best')

#pl.show()
name_file_save='5_Kintemp_run%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
name_file_save=name_file_save.replace('.','_')
name_file_save+='.pdf'

fig.savefig(name_file_save, dpi=300)
