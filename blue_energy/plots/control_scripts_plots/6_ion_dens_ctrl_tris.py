#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
#import datetime
#from pylab import *
from sys import argv
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


'''
FULL LIST OF PLOTS SCRIPTS.
Select the script to lauch from this list:

/5_plot_charge.py
'''
script_name='/6_plot_ioni_dens_v2.py' #extracts only the name of the script
#'/6_plot_ioni_dens_NO_GRAD.py'
#put in this list all the paths in which you want to launch script[i]

'''
FULL LIST OF DATA FOLDERS IN EACH OF THE PATHS ABOVE
/1_density_results
/2_convert_position2VMD_results
/3_RDF_results
/4_coord_num_results
/5_charge_results	
'''
data_path='/6_acc_surf_area_5_67_Fict'
#'/6_acc_surf_area_8_69_Ne'
#'/6_acc_surf_area_5_67_FictNG'
#'/6_acc_surf_area_5_67_Fict'

end_name=data_path.split('area_')[1]
'''
FULL LIST OF PATHS ROOTS. 
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
#'/cdc1200/nacl/8800_160/0.0V',
'/cdc1200/nacl/8800_160/1.0V',
#'/cdc1200/nacl/8800_80/0.0V',
#'/cdc1200/nacl/8800_80/1.0V',
#'/cdc800/nacl/7615_139/0.0V',
#'/cdc800/nacl/7615_139/1.0V',
#'/cdc800/nacl/7700_70/0.0V',
#'/cdc800/nacl/7700_70/1.0V',
#'/cdc1200/kcl/8800_160/0.0V',
#'/cdc1200/kcl/8800_160/1.0V',
#'/cdc1200/kcl/8800_80/0.0V',
#'/cdc1200/kcl/8800_80/1.0V'
]



#copy the script in the folder where it must be executed
for target_path in paths_to_launch:
        cmd='cp '+base_dir_script+script_name+' '+base_dir_results+target_path+data_path
        #print cmd
        os.system(cmd)

proc_id=[]
pid_list=[]
indx_buf=0
for target_path in paths_to_launch:
        #go in the directory where you have just copied the script
        dir_to_go=base_dir_results+target_path+data_path
        os.chdir(dir_to_go)
        list_arg_script =[]


        list_commands=['python',script_name[1:], end_name]
        list_commands+=list_arg_script
        #list_commands.append('&')
        #print list_commands
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
        #create the path if not present!
        store_path=base_dir_plot_storage+target_path
        cmd='mkdir -p '+store_path
        os.system(cmd)

        cmd='mv '+base_dir_results+target_path+data_path+'/*.pdf '+store_path
        #print cmd
        os.system(cmd)