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

'/cdc1200/kcl/8800_160/0.0V',
#'/cdc1200/kcl/8800_160/1.0V',
'/cdc1200/kcl/8800_80/0.0V',
#'/cdc1200/kcl/8800_80/1.0V',
#'/cdc1200/nacl/8800_160/0.0V', 
#'/cdc1200/nacl/8800_160/1.0V', 
#'/cdc1200/nacl/8800_80/0.0V',
#'/cdc1200/nacl/8800_80/1.0V'

#'/cdc800/nacl/7615_139/0.0V',
        #'/cdc800/nacl/7615_139/1.0V',

#'/cdc800/nacl/7700_70/0.0V',

                #'/cdc800/nacl/7700_70/1.0V',
]

data={}

file_Oxigen='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/results/general/data_coord_num_Ion-Oxige.txt'
file_Carbon='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/results/general/data_coord_num_Ion-Carbon.txt'

data_Oxigen=np.loadtxt(name_file_save, usecols=(1,2,3,4,5,6))
data_Carbon=np.loadtxt(name_file_save, usecols=(1,2,3,4,5,6))


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
        text= '& \multicolumn{3}{c}{%s, cdc1200 (KCl)} \\\ \n' %(spc)

        for target_path in paths_to_launch:
                #print target_path
                text+=target_path.split('/')[-2]+', '+target_path.split('/')[-1]
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
        text+='average '
        for index in [0,1,2]:
                text+=' & '
                #print '%d )  %6.*f +- %6.*f'  %(index,int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),avg_Na[index],int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),inc_Na[index])
                text+= '%6.*f $\pm$ %6.*f \t'  %(int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),avg_Na[index],int(('%1.0e'%(inc_Na[index])).split('e-')[-1]),inc_Na[index])
        text+='\\\ \n'
        return text

text=print_info('Na')
print text

text=print_info('Cl')
print text


text=print_info('O')
print text
