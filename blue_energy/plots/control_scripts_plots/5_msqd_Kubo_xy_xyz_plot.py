#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
from sys import argv
from os import listdir
from os.path import isfile, join
############################################################

runmin=0
def get_run_max():
        data_runmax=open('current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])-1
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()


#Irta Paths
#common_parent_path='/data_prace/michele/prace/cdc'
#base_dir_script=common_parent_path+'/analysis/plots/plot_low_level_scripts'
#base_dir_plot_storage=common_parent_path+'/analysis/plots/pdf_files'
#base_dir_results=common_parent_path+'/analysis/results/'

#latitude Paths
common_parent_path='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing'
base_dir_results=common_parent_path+'/results'
base_dir_script=common_parent_path+'/plots/plot_low_level_scripts'
base_dir_plot_storage=common_parent_path+'/plots/pdf_files'


script_name='/5_plot_conductivity.py' #extracts only the name of the script


unique_path_to_launch='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/results/general/KUBO_conductivity_xy_xyz'



#copy the script in the folder where it must be executed
cmd='cp '+base_dir_script+script_name+' '+unique_path_to_launch
os.system(cmd)

file_list = [f for f in listdir(unique_path_to_launch) if isfile(join(unique_path_to_launch, f))]

root_list=[]
for data_file in file_list:
        if (data_file!=script_name[1:]):
                root_file=data_file.split('_spc')[0]
                #print root_file
                root_file=root_file.split('C_')[1]
                if root_file not in root_list:
                        root_list.append(root_file)

os.chdir(unique_path_to_launch)
proc_id=[]
indx_buf=0
pid_list=[]
for root_name in root_list:
        list_commands=['python',script_name[1:], root_name]
        proc_id.append(1)
        proc_id[indx_buf] = subprocess.Popen(list_commands)
        pid_list.append(proc_id[indx_buf].pid)
        indx_buf+=1

#wait all the processes launched to complete
for j in range(0, indx_buf):
        proc_id[j].wait()
proc_id=[]
indx_buf=0

file_list = [f for f in listdir(unique_path_to_launch) if isfile(join(unique_path_to_launch, f))]

list_pdf_file=[]
for all_file in file_list:
        if all_file[-4:]=='.pdf':
                list_pdf_file.append(all_file)

#transfer data to the proper directory and clean tmp direatories
print list_pdf_file

for pdf_file in list_pdf_file:
        #create the path if not present!
        target_path=pdf_file.split('_Kubo_msd')[1]
        target_path=target_path.split('.pd')[0]
        target_path=target_path.replace('_','/')
        target_path=target_path.replace('1/0V','1.0V')
        target_path=target_path.replace('0/0V','0.0V')
        target_path=target_path.replace('8800/80','8800_80')
        target_path=target_path.replace('8800/160','8800_160')
        target_path=target_path.replace('7615/139','7615_139')
        target_path=target_path.replace('7700/70','7700_70')
        print target_path
        store_path=base_dir_plot_storage+'/'+target_path
        cmd='mkdir -p '+store_path
        os.system(cmd)
        here=os.getcwd()
        cmd='mv '+here+'/'+pdf_file+' '+store_path
        print store_path
        #print cmd
        os.system(cmd)