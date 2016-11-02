#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy as np
import subprocess
#import datetime
#from pylab import *
from sys import argv
import fileinput
import math
import inspect, os
############################################################
def get_run_max():
        data_runmax=open('../current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])-1
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()
#print dict_runmax

dict_runmax={'/cdc1200/nacl/8800_80/0.0V': 42, '/cdc1200/nacl/8800_160/1.0V': 44, '/cdc1200/kcl/8800_80/0.0V': 18, '/cdc1200/kcl/8800_160/1.0V': 13, '/cdc800/nacl/7615_139/0.0V': 46, '/cdc800/nacl/7700_70/1.0V': 43, '/cdc1200/kcl/8800_160/0.0V': 16, '/cdc1200/nacl/8800_160/0.0V': 44, '/cdc1200/nacl/8800_80/1.0V': 45, '/cdc1200/kcl/8800_80/1.0V': 15, '/cdc800/nacl/7615_139/1.0V': 41, '/cdc800/nacl/7700_70/0.0V': 44}

#run corresponding to 2.5 ns 
dict_runmin={'/cdc1200/nacl/8800_80/0.0V': 37, '/cdc1200/nacl/8800_160/1.0V': 38, '/cdc1200/kcl/8800_80/0.0V': 16, '/cdc1200/kcl/8800_160/1.0V': 13, '/cdc800/nacl/7615_139/0.0V': 41, '/cdc800/nacl/7700_70/1.0V': 38, '/cdc1200/kcl/8800_160/0.0V': 14, '/cdc1200/nacl/8800_160/0.0V': 37, '/cdc1200/nacl/8800_80/1.0V': 39, '/cdc1200/kcl/8800_80/1.0V': 14, '/cdc800/nacl/7615_139/1.0V': 35, '/cdc800/nacl/7700_70/0.0V': 39}




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
base_dir_plot_tables=common_parent_path+'/plots/data_table'

'''
FULL LIST OF PLOTS SCRIPTS.
Select the script to lauch from this list:
/1_plot_density_AVG.py
/5_plot_charge.py
'''
script_name='/1_plot_density_AVG.py' #extracts only the name of the script

#put in this list all the paths in which you want to launch script[i]

'''
FULL LIST OF DATA FOLDERS IN EACH OF THE PATHS ABOVE
/1_density_results
/2_convert_position2VMD_results
/3_RDF_results
/4_coord_num_results
/5_charge_results       
'''
data_path='/1_density_results'


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

#'/cdc1200/kcl/8800_160/0.0V',
#'/cdc1200/kcl/8800_160/1.0V',
#'/cdc1200/kcl/8800_80/0.0V',
#'/cdc1200/kcl/8800_80/1.0V',
#'/cdc1200/nacl/8800_160/0.0V', 
'/cdc1200/nacl/8800_160/1.0V', 
#'/cdc1200/nacl/8800_80/0.0V',
'/cdc1200/nacl/8800_80/1.0V'

#'/cdc800/nacl/7615_139/0.0V',
#'/cdc800/nacl/7615_139/1.0V',

#'/cdc800/nacl/7700_70/0.0V',

#'/cdc800/nacl/7700_70/1.0V',
]

data={}
for target_path in paths_to_launch:
        runmax=dict_runmax[target_path]
        runmin=dict_runmin[target_path]
        store_path_data=base_dir_plot_tables+target_path
        os.chdir(store_path_data)
        name_folder=os.getcwd() # script directory
        list_dir=name_folder.split('/')
        #print list_dir
        name_file_save=r'c_1_AVG_density_%03d-run%03d_%s_%s_%s_%s.txt' %(runmin,runmax, list_dir[9],list_dir[10],list_dir[11],list_dir[12].replace('.','_') )
        data[target_path]=np.loadtxt(name_file_save, usecols=(1,2,3,4,5,6))


avg_Na=np.zeros(3)
avg_Cl=np.zeros(3)
avg_O=np.zeros(3)

inc_Na=np.zeros(3)
inc_Cl=np.zeros(3)
inc_O=np.zeros(3)

def print_info(spc):
        avg_Na=np.zeros(3)
        inc_Na=np.zeros(3)
        if (spc=='Na') or (spc=='K'):
                col_num=0
        elif (spc=='Cl'):
                col_num=1
        elif (spc=='O'):
                col_num=2

        if 'cdc1200' in paths_to_launch[0]:
                text= '\hline\n & \multicolumn{3}{c}{%s, cdc1200 (KCl)} \\\ \n\hline\n' %(spc)
        else:
                text= '\hline\n & \multicolumn{3}{c}{%s, cdc800 (KCl)} \\\ \n\hline\n' %(spc)
        for target_path in paths_to_launch:
                #print target_path
                if (target_path.split('/')[-2]=='8800_160' or target_path.split('/')[-2]=='7615_139'):
                        text+='1.0M'
                else:
                        text+='0.5M'
                for index in [0,1,2]:
                        text+= ' & '
                        avg_Na[index]+=data[target_path][col_num,2*index]
                        inc_Na[index]+=data[target_path][col_num,2*index+1]**2
                        buf_print_avg=data[target_path][col_num,2*index]
                        buf_print_inc=data[target_path][col_num,2*index+1]
                        #print '%d )  %6.*f +- %6.*f'  %(index,int(('%1.0e'%(buf_print_inc)).split('e-')[-1]),buf_print_avg,int(('%1.0e'%(buf_print_inc)).split('e-')[-1]),buf_print_inc)
                        text+= r'%6.*f $\pm$ %6.*f '  %(int(('%1.0e'%(buf_print_inc)).split('e-')[-1]),buf_print_avg,int(('%1.0e'%(buf_print_inc)).split('e-')[-1]),buf_print_inc)
                text+='\\\ \n'  

        avg_Na=avg_Na/len(paths_to_launch)
        inc_Na=np.sqrt(inc_Na)/len(paths_to_launch)
        #print '---average---'
        '''
        text+='average '
        for index in [0,1,2]:
                text+=' & '
                #print '%d )  %6.*f +- %6.*f'  %(index,int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),avg_Na[index],int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),inc_Na[index])
                text+= '%6.*f $\pm$ %6.*f \t'  %(int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),avg_Na[index],int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),inc_Na[index])
        text+='\\\ \n'
        '''
        return text

text=print_info('K')
print text

text=print_info('Cl')
print text


text=print_info('O')
print text
