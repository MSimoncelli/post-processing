import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os
from collections import defaultdict


from sys import argv

global n_ions
############################################################
from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################
conv=0.529177
max_dist=186.406 #length of the whole system

def get_run_max():
        data_runmax=open('/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/plots/current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()

end_name='ERR_EnD_Name_missing'


#list of directories in which we have to launch the script
list_sub_dir_here=os.walk('.').next()[1]
list_sub_dir_here=['/'+s for s in list_sub_dir_here]

############################################
data_path={}
num_confs={}


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
#print name_check

common_parent_path='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing'


def say_me_run(time_nanosec):
        base_dir_script=common_parent_path+'/plots/plot_low_level_scripts/'
        file_frames=open(base_dir_script+'frames_infos.dat', 'r') 
        file_frames_lines= file_frames.readlines()      
        num_confs=[]    

        numframes_max=time_nanosec*1000  

        num_confs = defaultdict(list)
        for target_path in paths_to_launch:
                for irun in range(0,dict_runmax[target_path]):
                        #print target_path+'/run%03d'%(irun)
                        for line in file_frames_lines:
                                if target_path+'/run%03d'%(irun) in line:
                                        num_confs[target_path].append(int((line.split('n_frames=')[1]).split('\n')[0])) 

        num_frames={}
        run_interest={}
        for target_path in paths_to_launch:
                num_frames[target_path]=0
                for irun in range(0,len(num_confs[target_path][:])):
                        if (num_frames[target_path]<numframes_max):
                                num_frames[target_path]+=num_confs[target_path][irun]
                                run_interest[target_path]=irun  

        return run_interest

dict_runs={}

dict_runs=say_me_run(2.5)

print dict_runs

dict_runs=say_me_run(3.5)

print dict_runs