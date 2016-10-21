#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
#import datetime
#from pylab import *
from sys import argv
num_max_proc=32

def get_run_max():
        data_runmax=open('current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])
        data_runmax.close()
        return dict_runmax;

dict_runmax=get_run_max()

def is_process_running(process_id):
    try:
        os.kill(process_id, 0)
        return True
    except OSError:
        return False
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



script_folder_2='/5_charge' #pay attention at the arguments
script_name_2='/control_slice.py' #extracts only the name of the script

#put in this list all the paths in which you want to launch script[i]
paths_to_launch=[
#'/cdc1200/nacl/8800_160/0.0V', 
#'/cdc1200/nacl/8800_160/1.0V',
#'/cdc1200/nacl/8800_80/0.0V',
#'/cdc1200/nacl/8800_80/1.0V',
#'/cdc800/nacl/7615_139/0.0V',
#'/cdc800/nacl/7615_139/1.0V',
#'/cdc800/nacl/7700_70/0.0V',
#'/cdc800/nacl/7700_70/1.0V',
#'/cdc1200/kcl/8800_160/0.0V',
#'/cdc1200/kcl/8800_160/1.0V',
'/cdc1200/kcl/8800_80/0.0V'
#'/cdc1200/kcl/8800_80/1.0V'
]

paths_to_store=['/analysis/results'+s for s in paths_to_launch]

#copy the script in the folder where it must be executed
for target_path in paths_to_launch:
	cmd='cp -r '+base_dir_script+script_folder_2+' '+common_parent_path+target_path
	#print cmd
	os.system(cmd)

proc_id=[]
pid_list=[]
indx_buf=0
launched=0
full_cmd=''
for target_path in paths_to_launch:
	runmax=dict_runmax[target_path]
	list_arg_script_2=['0','%d'%(runmax), '%d'%(num_max_proc)]
	print target_path, runmax
	#go in the directory where you have just copied the script
	dir_to_go=common_parent_path+target_path+script_folder_2
	print 'cd ', dir_to_go
	os.chdir(dir_to_go)
	list_commands=['python',script_name_2[1:]]
	list_commands+=list_arg_script_2
	full_cmd=''
	for string in list_commands:
		full_cmd+=string+' '	#list_commands.append('&')
	print full_cmd
	os.system(full_cmd)
	print '----------\n'
	
print 'DONE 1st part'
#transfer data to the proper directory and clean tmp direatories
for target_path in paths_to_launch:
	name_new_folder=script_folder_2+'_results'
	
	cmd='mv '+common_parent_path+target_path+script_folder_2+'/data_OUT '+common_parent_path+'/analysis/results'+target_path
	print cmd
	os.system(cmd)
	
	os.chdir(common_parent_path+'/analysis/results'+target_path)
	print 'cd ',common_parent_path+'/analysis/results'+target_path
	cmd='mkdir -p '+name_new_folder[1:]
	print cmd
	os.system(cmd)
	cmd='cp -a data_OUT/. '+name_new_folder[1:]+'/'
	print cmd
	os.system(cmd)
	#cmd='rm -r data_OUT'
	#print cmd
	#os.system(cmd)

	#cmd ='rm -r '+common_parent_path+target_path+script_folder_2
	#print cmd
	#os.system(cmd)
