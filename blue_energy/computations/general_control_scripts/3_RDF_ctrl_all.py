#!/usr/bin/python
import string, re, struct, sys, math, os, inspect
import subprocess
from sys import argv
import numpy as np
startdir = os.getcwd() #gets the current working directory

#all the files will start from this path
common_parent_path='/data_prace/michele/prace/cdc'
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

def get_run_max():
        data_runmax=open('current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])-1
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()

runmin=0
num_max_proc=24



paths_to_store=['/analysis/results'+s for s in paths_to_launch]
sp=[]
indx_buf=0
FNULL = open(os.devnull, 'w')

#copy the script in the folder where it must be executed
for target_path in paths_to_launch:
	launched=0
	runmax=dict_runmax[target_path]

        for irun in range(runmin,runmax+1):
		# move to directory and uncompress files
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
		os.chdir(dir_to_go)        #changes directory
		#copy the script in the run folder
		cmd_list = ['tar', '-xzf', 'positions.tar.gz']
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(cmd_list)
		indx_buf+=1
		cmd_list_1 = ['cp', '/data_prace/michele/prace/cdc/analysis/scripts/3_RDF/calc_RDF_fort.x', dir_to_go]
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(cmd_list_1)
		indx_buf+=1
                cmd_list_1 = ['cp', '/data_prace/michele/prace/cdc/analysis/scripts/3_RDF/all_my_subroutines.mod', dir_to_go]
                sp.append(1)
                sp[indx_buf] = subprocess.Popen(cmd_list_1)
                indx_buf+=1
		cmd_list_1 = ['cp', '/data_prace/michele/prace/cdc/analysis/scripts/3_RDF/qsort_c_module.mod', dir_to_go]
                sp.append(1)
                sp[indx_buf] = subprocess.Popen(cmd_list_1)
                indx_buf+=1


		for j in range(0, indx_buf):
        		sp[j].wait()
	
	print 'untar completed'       
	pid_list=[]
	sp=[]
	indx_buf=0
	launched=0
	for irun in range(runmin,runmax+1):
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
                os.chdir(dir_to_go) 
		print dir_to_go
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(['./calc_RDF_fort.x','%03d'%(irun) ],stdout=subprocess.PIPE)
		#(out, err) = sp[indx_buf].communicate()
		pid_list.append(sp[indx_buf].pid)
		indx_buf+=1
		launched+=1
		print launched
                check=min([num_max_proc,runmax-runmin+1])
                if (launched>=check ):
                        for i in range(check,0,-1):
                                sp[indx_buf-i].wait()
                        launched=0
	for proc in sp:
		proc.wait()


	for irun in range(runmin,runmax+1):
	        #changes directory
		path_store = '/data_prace/michele/prace/cdc/analysis/results'+target_path+'/3_RDF_all'
		cmd = 'mkdir -p '+path_store
		os.system(cmd)
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
                os.chdir(dir_to_go)

		cmd= 'mv run%03dRDF* '%(irun)+path_store
		os.system(cmd)

		cmd='rm -r calc_RDF_fort.x'
		os.system(cmd)
		cmd ='rm *.out'
		#print cmd
		os.system(cmd)

		cmd ='rm *.mod'
		os.system(cmd)

