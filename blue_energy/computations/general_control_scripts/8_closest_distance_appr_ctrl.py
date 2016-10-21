#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
#import datetime
#from pylab import *
from sys import argv
############################################################
# first and last runs to include in average density
#print 'START'
#print datetime.datetime.now()
startdir = os.getcwd() #gets the current working directory

#all the files will start from this path
common_parent_path='/data_prace/michele/prace/cdc'

base_dir_script=common_parent_path+'/analysis/scripts'

'''
FULL LIST SCRIPTS.
Select the script to lauch from this list:

/1_density/density_AVG.py [runmin] [runmax]
/1_density/density_position.py [runmin] [runmax]

/2_convert_position2VMD/convert_VMD.py 

/3_RDF/control_RDF_sp.py [runmin] [runmax]

/4_coord_number/control_COOR_NUM_sp.py [runmin] [runmax]  [nome_file_cutoffs]
					0-46	0-46	  cutoff_dataset_NaCl_1.txt	
					0-46 	0-46	  cutoff_dataset_NaCl_2.txt
					0-46	0-46	  cutoff_dataset_NaCl_3.txt
/5_charge/catfiles.sh
'''
runmax=1
script_folder_1='/8_closest_distance_of_approach' #pay attention at the arguments
script_name_1='/control_closest_approach.py' #extracts only the name of the script
maxrun='%d'%(runmax)
list_arg_script_1=['0', maxrun]

'''
FULL LIST OF PATHS. 
'/cdc1200/nacl/8800_160/0.0V', 
'/cdc1200/nacl/8800_160/1.0V', 
'/cdc1200/nacl/8800_80/0.0V',
'/cdc1200/nacl/8800_80/1.0V',
'/cdc1200/kcl/8800_160/0.0V',
'/cdc1200/kcl/8800_160/1.0V',
'/cdc1200/kcl/8800_80/0.0V',
'/cdc1200/kcl/8800_80/1.0V',

'/cdc800/nacl/7615_139/0.0V',
'/cdc800/nacl/7615_139/1.0V',
'/cdc800/nacl/7700_70/0.0V',
'/cdc800/nacl/7700_70/1.0V',
'/cdc800/kcl/7615_139/0.0V',
'/cdc800/kcl/7615_139/1.0V',
'/cdc800/kcl/7700_70/0.0V',
'/cdc800/kcl/7700_70/1.0V'
'''

#put in this list all the paths in which you want to launch script[i]
paths_to_launch=[
'/cdc1200/nacl/8800_160/0.0V'#,
#'/cdc1200/nacl/8800_160/1.0V', 
#'/cdc1200/nacl/8800_80/0.0V',
#'/cdc1200/nacl/8800_80/1.0V',
#'/cdc800/nacl/7615_139/0.0V',
#'/cdc800/nacl/7615_139/1.0V',
#'/cdc800/nacl/7700_70/0.0V',
#'/cdc800/nacl/7700_70/1.0V'
]

paths_to_store=['/analysis/results'+s for s in paths_to_launch]

#copy the script in the folder where it must be executed
for target_path in paths_to_launch:
	cmd='cp -r '+base_dir_script+script_folder_1+' '+common_parent_path+target_path
	print cmd
	os.system(cmd)

proc_id=[]
pid_list=[]
indx_buf=0
for target_path in paths_to_launch:
	#go in the directory where you have just copied the script
	dir_to_go=common_parent_path+target_path+script_folder_1
	print 'cd ', dir_to_go
	os.chdir(dir_to_go)
	list_commands=['python',script_name_1[1:]]
	list_commands+=list_arg_script_1
	#list_commands.append('&')
	print list_commands
	proc_id.append(1)
	proc_id[indx_buf] = subprocess.Popen(list_commands)
	pid_list.append(proc_id[indx_buf].pid)
	indx_buf+=1

#wait all the processes launched to complete
for j in range(0, indx_buf):
	proc_id[j].wait()
proc_id=[]
indx_buf=0


#transfer data to the proper directory and clean tmp direatories
for target_path in paths_to_launch:
	cmd='mkdir -p '+common_parent_path+'/analysis/results'+target_path+'/8_closest_dist_approach_res'
	os.system(cmd)
	cmd='mv '+common_parent_path+target_path+script_folder_1+'/data_OUT '+common_parent_path+'/analysis/results'+target_path+'/8_closest_dist_approach_res/'
	print cmd
	os.system(cmd)
	
	os.chdir(common_parent_path+'/analysis/results'+target_path+'/8_closest_dist_approach_res')
	print 'cd ',common_parent_path+'/analysis/results'+target_path+'/8_closest_dist_approach_res'
	
	#cmd ='rm -r '+common_parent_path+target_path+script_folder_1
	#print cmd
	#os.system(cmd)
