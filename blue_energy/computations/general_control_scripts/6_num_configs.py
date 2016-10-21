#!/usr/bin/python
import string, re, struct, sys, math, os, inspect
import subprocess
from sys import argv
import numpy as np
startdir = os.getcwd() #gets the current working directory

def read_numframes():
        name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
        list_dir=name_folder.split('/')
        num_run=list_dir[-1]
        n_rn=num_run.split('n')[1]
        fin=int(n_rn)
        dirrun='../run%03d'%(fin+1)
        #print dirrun
        os.chdir(dirrun)
        file_box = open( 'restart.dat', 'r' )
        line1= file_box.readline().split() #returns a list containing the lines
        file_box.close()
        numframes=int(line1[1])/1000
        dirrun='../run%03d'%(fin)
        os.chdir(dirrun)
        #print dirrun
        return numframes

def get_run_max():
        data_runmax=open('current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])
        data_runmax.close()
        return dict_runmax;

dict_runmax=get_run_max()
#all the files will start from this path
common_parent_path='/data_prace/michele/prace/cdc'
path_base_surf=['/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area_8_69_Ne',
		'/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area_9_45_Ne',
		'/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area_9_64_Ar',
		'/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area_10_02_Ar']
	
path_cdc800_el_1 ='/electrode_cdc800/surf_points_cdc800_el1.out'
path_cdc800_el_2 ='/electrode_cdc800/surf_points_cdc800_el2.out'
path_cdc1200_el_1='/electrode_cdc120/surf_points_cdc1200_el1.out'
path_cdc1200_el_2='/electrode_cdc120/surf_points_cdc1200_el2.out'

path_cdc_800_list=[path_cdc800_el_1,path_cdc800_el_2]
path_cdc_1200_list=[path_cdc1200_el_1,path_cdc1200_el_2]
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
#/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area
script_folder_2='/6_acc_surf_area' #pay attention at the arguments
script_name_2='/idens_local_surf' #extracts only the name of the script
#put in this list all the paths in which you want to launch script[i]
paths_to_launch=[
'/cdc1200/nacl/8800_160/0.0V',
'/cdc1200/nacl/8800_160/1.0V', 
'/cdc1200/nacl/8800_80/0.0V',
'/cdc1200/nacl/8800_80/1.0V',
'/cdc800/nacl/7615_139/0.0V',
'/cdc800/nacl/7615_139/1.0V',
'/cdc800/nacl/7700_70/0.0V',
'/cdc800/nacl/7700_70/1.0V',
'/cdc1200/kcl/8800_160/0.0V',
'/cdc1200/kcl/8800_160/1.0V',
'/cdc1200/kcl/8800_80/0.0V',
'/cdc1200/kcl/8800_80/1.0V'
]

runmin=0

num_max_proc=32

paths_to_store=['/analysis/results'+s for s in paths_to_launch]
sp=[]
indx_buf=0
FNULL = open(os.devnull, 'w')

#copy the script in the folder where it must be executed
text=''
for target_path in paths_to_launch:      
	runmax=dict_runmax[target_path]
	for irun in range(runmin,runmax+1):
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
		os.chdir(dir_to_go)        #changes directory
		num_frames=read_numframes()
		text+=target_path+'/run%03d'%(irun)+' n_frames=%d\n'%(num_frames)
	text+='\n'

os.chdir(startdir)
tmpf_2=open('frames_infos.dat', 'w')
#print 'RDF_'+int_RDF[x]+'.out'
tmpf_2.write(text)
tmpf_2.close()

