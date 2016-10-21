#!/usr/bin/python
'''
this script will generate a file containing the number of available runs for each system.
I have notice that usually the lastest run available is broken, hence here I count all the runs folder-1
'''
import string, re, struct, sys, math, os, time
import numpy
import subprocess
#import datetime
#from pylab import *
from sys import argv
############################################################

startdir = os.getcwd() #gets the current working directory

#all the files will start from this path
my_common_parent_path='/data_prace/michele/prace/cdc'
common_parent_path='/data_prace/ben/prace/cdc'
base_dir_script=common_parent_path+'/analysis/scripts'


script_folder_2='/5_charge' #pay attention at the arguments
script_name_2='/catfiles.py' #extracts only the name of the script
arg_script_2=' 0 53'
#put in this list all the paths in which you want to launch script[i]
paths_to_launch=[
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
#'/cdc800/kcl/7615_139/1.0V',
'/cdc800/kcl/7700_70/0.0V'#,
#'/cdc800/kcl/7700_70/1.0V'
]

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

max_num_runs={}
#copy the script in the folder where it must be executed
for target_path in paths_to_launch:
	list_num_runs=[]
	curr_path=common_parent_path+target_path
	directories_here=get_immediate_subdirectories(curr_path)
	for nam_dir in directories_here:
		if(nam_dir.startswith( 'run' )and (len(nam_dir)<7)):
			list_num_runs.append(int(nam_dir.split('n')[1]))
	max_runs=0
	if len(list_num_runs)>0:
		max_runs=max(list_num_runs)-1	
	max_num_runs[target_path]=max_runs

	#rsync -aq --exclude 'openmpi-sessions-*' --exclude '*.tar.gz.*' /data_prace/ben/prace/cdc/cdc1200/nacl/8800_160/0.0V/ /data_prace/michele/prace/cdc/cdc1200/nacl/8800_160/0.0V/
	#/data_prace/ben/prace/cdc/cdc1200/nacl/8800_160/1.0V/run083/tmpdata/.run.out.loXCQb
	cmd='rsync -aq --exclude \'openmpi-sessions-*\' --exclude \'*.testout.rst.*\' --exclude \'*tmpdata/.run.out.*\' --exclude \'*.tar.gz.*\' '+curr_path+'/ '+my_common_parent_path+target_path+'/ '
	print cmd 
	os.system(cmd)

text=''
for data in paths_to_launch:
	text+='%s\t%d\n' %(data,max_num_runs[data])
os.chdir(startdir)
tmpf=open('current_state_simulation_b.dat', 'w')
#print 'RDF_'+int_RDF[x]+'.out'
tmpf.write(text)
tmpf.close()
