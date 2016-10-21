#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
#import datetime
#from pylab import *
from sys import argv
############################################################

runmin=0
runmax=53


#Irta Paths
irta_common_parent_path='/data_prace/michele/prace/cdc'
irta_base_dir_script=irta_common_parent_path+'/analysis/plots/plot_low_level_scripts'
irta_base_dir_plot_storage=irta_common_parent_path+'/analysis/plots/pdf_files'
irta_base_dir_results=irta_common_parent_path+'/analysis/results'


'''
FULL LIST OF PLOTS SCRIPTS.
Select the script to lauch from this list:

/5_plot_charge.py
'''

'''
FULL LIST OF DATA FOLDERS IN EACH OF THE PATHS ABOVE
/1_density_results
/2_convert_position2VMD_results
/3_RDF_results
/4_coord_num_results
/5_charge_results	
'''



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

paths_to_sync=[
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

folders_in_path_to_sync=[
'/1_density_results',
#'/2_convert_position2VMD_results',
'/3_RDF_results',
'/4_coord_num_results',
'/5_charge_results'
]

#create all_the_missing paths on the dell latitude
for target_path in paths_to_sync:
	current=irta_base_dir_results+target_path
	print current
	os.chdir(current)
	cmd='rm -rf 5_charge_results'
	os.system(cmd)


