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
dict_runmax2={'/cdc1200/nacl/8800_80/0.0V': 37, '/cdc1200/nacl/8800_160/1.0V': 38, '/cdc1200/kcl/8800_80/0.0V': 16, '/cdc1200/kcl/8800_160/1.0V': 13, '/cdc800/nacl/7615_139/0.0V': 41, '/cdc800/nacl/7700_70/1.0V': 38, '/cdc1200/kcl/8800_160/0.0V': 14, '/cdc1200/nacl/8800_160/0.0V': 37, '/cdc1200/nacl/8800_80/1.0V': 39, '/cdc1200/kcl/8800_80/1.0V': 14, '/cdc800/nacl/7615_139/1.0V': 35, '/cdc800/nacl/7700_70/0.0V': 39}
dict_runmin2={'/cdc1200/nacl/8800_80/0.0V': 32, '/cdc1200/nacl/8800_160/1.0V': 33, '/cdc1200/kcl/8800_80/0.0V': 13, '/cdc1200/kcl/8800_160/1.0V': 12, '/cdc800/nacl/7615_139/0.0V': 35, '/cdc800/nacl/7700_70/1.0V': 33, '/cdc1200/kcl/8800_160/0.0V': 12, '/cdc1200/nacl/8800_160/0.0V': 32, '/cdc1200/nacl/8800_80/1.0V': 31, '/cdc1200/kcl/8800_80/1.0V': 11, '/cdc800/nacl/7615_139/1.0V': 30, '/cdc800/nacl/7700_70/0.0V': 33}

dict_runmax3={'/cdc1200/nacl/8800_80/0.0V': 42, '/cdc1200/nacl/8800_160/1.0V': 44, '/cdc1200/kcl/8800_80/0.0V': 18, '/cdc1200/kcl/8800_160/1.0V': 13, '/cdc800/nacl/7615_139/0.0V': 46, '/cdc800/nacl/7700_70/1.0V': 43, '/cdc1200/kcl/8800_160/0.0V': 16, '/cdc1200/nacl/8800_160/0.0V': 44, '/cdc1200/nacl/8800_80/1.0V': 45, '/cdc1200/kcl/8800_80/1.0V': 15, '/cdc800/nacl/7615_139/1.0V': 41, '/cdc800/nacl/7700_70/0.0V': 44}
dict_runmin3={'/cdc1200/nacl/8800_80/0.0V': 37, '/cdc1200/nacl/8800_160/1.0V': 38, '/cdc1200/kcl/8800_80/0.0V': 16, '/cdc1200/kcl/8800_160/1.0V': 13, '/cdc800/nacl/7615_139/0.0V': 41, '/cdc800/nacl/7700_70/1.0V': 38, '/cdc1200/kcl/8800_160/0.0V': 14, '/cdc1200/nacl/8800_160/0.0V': 37, '/cdc1200/nacl/8800_80/1.0V': 39, '/cdc1200/kcl/8800_80/1.0V': 14, '/cdc800/nacl/7615_139/1.0V': 35, '/cdc800/nacl/7700_70/0.0V': 39}



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
'/cdc1200/kcl/8800_80/1.0V',
#'/cdc1200/nacl/8800_160/0.0V', 
#'/cdc1200/nacl/8800_160/1.0V', 
#'/cdc1200/nacl/8800_80/0.0V',
#'/cdc1200/nacl/8800_80/1.0V'
#'/cdc800/nacl/7615_139/0.0V',
#'/cdc800/nacl/7615_139/1.0V',

#'/cdc800/nacl/7700_70/0.0V',

#'/cdc800/nacl/7700_70/1.0V',
]

def pwu(a,b):
        c=b
        i=0
        while int(c)<1:
                c=c*10
                i+=1
        return r'%.*f$\pm%.*f$'%(i,a,i,b)

data2={}
data3={}
for target_path in paths_to_launch:
        runmax2=dict_runmax2[target_path]
        runmin2=dict_runmin2[target_path]
        runmax3=dict_runmax3[target_path]
        runmin3=dict_runmin3[target_path]
        store_path_data=base_dir_plot_tables+target_path
        os.chdir(store_path_data)
        name_folder=os.getcwd() # script directory
        list_dir=name_folder.split('/')
        #print list_dir
        name_file_save2=r'd_1_AVG_density_%03d-run%03d_%s_%s_%s_%s.txt' %(runmin2,runmax2, list_dir[9],list_dir[10],list_dir[11],list_dir[12].replace('.','_') )
        data2[target_path]=np.loadtxt(name_file_save2, usecols=(1,2,3), skiprows=1)
        name_file_save3=r'd_1_AVG_density_%03d-run%03d_%s_%s_%s_%s.txt' %(runmin3,runmax3, list_dir[9],list_dir[10],list_dir[11],list_dir[12].replace('.','_') )
        data3[target_path]=np.loadtxt(name_file_save3, usecols=(1,2,3), skiprows=1)

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



        for target_path in paths_to_launch:
                if (target_path.split('/')[-2]=='8800_160' or target_path.split('/')[-2]=='7615_139'):
                        inpt='1.0M'
                else:
                        inpt='0.5M'
                avg_Na_0=(data2[target_path][col_num,0]+data3[target_path][col_num,0])*0.5
                inc_avg_Na_0=abs(data2[target_path][col_num,0]-data3[target_path][col_num,0])*0.5
                avg_Na_1=(data2[target_path][col_num,1]+data3[target_path][col_num,1])*0.5
                inc_avg_Na_1=abs(data2[target_path][col_num,1]-data3[target_path][col_num,1])*0.5
                avg_Na_2=(data2[target_path][col_num,2]+data3[target_path][col_num,2])*0.5
                inc_avg_Na_2=abs(data2[target_path][col_num,2]-data3[target_path][col_num,2])*0.5
               
        text=inpt+' (1) & '+'%1.8f\t&\t%1.8f\t&\t%1.8f \\\ \n'%(data3[target_path][col_num,0],data3[target_path][col_num,1],data3[target_path][col_num,2])
        text+=inpt+' (2) & '+'%1.8f\t&\t%1.8f\t&\t%1.8f \\\ \n'%(data2[target_path][col_num,0],data2[target_path][col_num,1],data2[target_path][col_num,2])
        text+='avg & '+pwu(avg_Na_0,inc_avg_Na_0)+' & '+pwu(avg_Na_1,inc_avg_Na_1)+' & '+pwu(avg_Na_2,inc_avg_Na_2)+r'\\'
        return text
print 'Na or K'
text=print_info('K')
print text
print ''
print 'Cl'
text=print_info('Cl')
print text


