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
from collections import defaultdict
from uncertainties import ufloat
############################################################
def get_run_max():
        data_runmax=open('../current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()
#print dict_runmax

def pwu(a,b):
    c=b
    i=0
    if (b<1 and b!=0):
     while int(c)<1:
             c=c*10
             i+=1
    return '%.*f$\pm$%.*f'%(i,a,i,b)

def pwu_list(lista):
    a=lista[0]
    b=lista[1]
    c=b
    i=0
    if (b<1 and b!=0):
     while int(c)<1:
             c=c*10
             i+=1
    return '%.*f$\pm$%.*f'%(i,a,i,b)

def mean_sem_disp(list_in,list_unc):
    mean_val=np.mean(list_in)
    sem_disp=0.5*abs(list_in[0]-list_in[1])
    inc_prop=0.5*(list_unc[0]**2+list_unc[0]**2)**0.5
    if sem_disp>inc_prop:
    	incert=sem_disp
    else:
    	incert=inc_prop
    return [mean_val, incert]



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
#'/cdc1200/nacl/8800_160/0.0V',
'/cdc1200/nacl/8800_160/1.0V', 
#'/cdc1200/nacl/8800_80/0.0V',
'/cdc1200/nacl/8800_80/1.0V',
#'/cdc800/nacl/7615_139/0.0V',
'/cdc800/nacl/7615_139/1.0V',
#'/cdc800/nacl/7700_70/0.0V',
'/cdc800/nacl/7700_70/1.0V',
#'/cdc1200/kcl/8800_160/0.0V',
'/cdc1200/kcl/8800_160/1.0V',
#'/cdc1200/kcl/8800_80/0.0V',
'/cdc1200/kcl/8800_80/1.0V'
]

paths_to_launch_nacl_1200=[
'/cdc1200/nacl/8800_160/1.0V', 
'/cdc1200/nacl/8800_80/1.0V'
]

paths_to_launch_nacl_800=[
'/cdc800/nacl/7615_139/1.0V',
'/cdc800/nacl/7700_70/1.0V'
]

paths_to_launch_kcl_1200=[
'/cdc1200/kcl/8800_160/1.0V',
'/cdc1200/kcl/8800_80/1.0V'
]

def compute_parameters(R_l,R_bulk,Q_1max,Q_2max):
	t_0=10**(-9)
	R_0=10**9
	C_0=1.60217662*10**(-19)
	V_0=1.0
	Q_0=1.60217662*10**(-19)	

	a_red_conv=(t_0/(R_0*C_0))
	b_red_conv=(t_0**2/(R_0**2*C_0**2))
	c_red_conv=(t_0**2*V_0/((R_0**2)*Q_0*C_0))
	d_red_conv=(t_0*V_0/(R_0*Q_0))	

	V_0=1.0 #voltage potential
	C_1=2.0*Q_1max/V_0
	C_2=2.0*Q_2max/V_0
	Q_max=Q_1max+Q_2max
	a=a_red_conv*((R_bulk+2*R_l)*C_1+(R_bulk+4*R_l)*C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	b=b_red_conv*2.0/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	c=c_red_conv*(C_1+C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	d=d_red_conv*1.0/(R_bulk+2*R_l)
	tau_1=(2.0/(a+(a**2-4.0*b)**0.5))
	tau_2=(2.0/(a-(a**2-4.0*b)**0.5))
	A_1=0.5*(1.0+((2.0*b*d-a*c)/(2*c*(a**2-4.0*b)**0.5)))
	A_2=0.5*(1.0-((2.0*b*d-a*c)/(2*c*(a**2-4.0*b)**0.5)))
	return tau_1,tau_2,A_1,A_2


def print_values(target_path):
	a=mean_sem_disp([data[target_path][0,0],data[target_path][1,0]],[data[target_path][0,1],data[target_path][1,1]])
	R_l = ufloat(a[0]*0.01, a[1]*0.01)
	R_bulk=data[target_path][0,2]*0.01
	b=mean_sem_disp([data[target_path][0,3],data[target_path][1,3]],[data[target_path][0,4],data[target_path][1,4]])
	Q_1max=ufloat(b[0], b[1])
	c=mean_sem_disp([data[target_path][0,5],data[target_path][1,5]],[data[target_path][0,6],data[target_path][1,6]])
	Q_2max=ufloat(c[0], c[1])

	tau_1,tau_2,A_1,A_2=compute_parameters(R_l,R_bulk,Q_1max,Q_2max)
	print r'$\tau_1$ [ns] & ',tau_1, ' \\\ '
	print r'$\tau_2$ [ns] & ',tau_2, ' \\\ ' 
	print '$A_1$ & ',A_1, ' \\\ ' 
	print '$A_2$ & ',A_2, ' \\\ ' 


data={}
for target_path in paths_to_launch:
	runmax=dict_runmax[target_path]
	name_file='/5_Charge_REV_sep_run000-run%03d_%s_%s_%s_%s.txt' %(runmax, target_path.split('/')[1], target_path.split('/')[2],target_path.split('/')[3], target_path.split('/')[4].replace('.','_'))
	#print name_file
	file_load=base_dir_plot_tables+target_path+name_file
	data[target_path]=np.loadtxt(file_load, usecols=(2,3,4,5,6,7,8,9,10))
text= '\hline \n' 
text +=' & '+' $R_l \;[10^7\cdot \Omega]$ &  $R_{bulk}  \;[10^7\cdot \Omega]$ & Q$_{max1}$ $[e]$ & Q$_{max2}$ $[e]$ & C$_{tot}$ '+r' $[\frac{F}{g}]$ \\'+ '\n' 
text+= '\hline \n' 
text+= '& \multicolumn{5}{c}{cdc1200 (NaCl)} \\\ \n' 
text+= '\hline \n' 
for target_path in paths_to_launch_nacl_1200:
	if '8800_160' in target_path:
		conc=1.0
	else:
		conc=0.5
	for i in [0,1]:
		text+='el %d, %1.1f M & '%(i+1, conc)+ pwu(data[target_path][i,0],data[target_path][i,1])+' & '+ '%1.6f'%data[target_path][i,2]+' & '+pwu(data[target_path][i,3],data[target_path][i,4])+' & '+pwu(data[target_path][i,5],data[target_path][i,6])+' & '+pwu(data[target_path][i,7],data[target_path][i,8])+'\\\  \n'
	text+= '\hline \n'
	text+= 'avg & '+pwu_list(mean_sem_disp([data[target_path][0,0],data[target_path][1,0]],[data[target_path][0,1],data[target_path][1,1]]))+' & '+ '%1.6f'%data[target_path][i,2]+' & '+pwu_list(mean_sem_disp([data[target_path][0,3],data[target_path][1,3]],[data[target_path][0,4],data[target_path][1,4]]))+' & '+pwu_list(mean_sem_disp([data[target_path][0,5],data[target_path][1,5]],[data[target_path][0,6],data[target_path][1,6]]))+' & '+pwu_list(mean_sem_disp([data[target_path][0,7],data[target_path][1,7]],[data[target_path][0,8],data[target_path][1,8]]))+'\\\  \n'
	text+= '\hline \n'
	text+= '\hline \n'
text+= '\hline \n' 

text+= '& \multicolumn{5}{c}{cdc800 (NaCl)} \\\ \n' 
text+= '\hline \n' 
for target_path in paths_to_launch_nacl_800:
	if '7615_139' in target_path:
		conc=1.0
	else:
		conc=0.5
	for i in [0,1]:
		text+='el %d, %1.1f M & '%(i+1, conc)+ pwu(data[target_path][i,0],data[target_path][i,1])+' & '+ '%1.6f'%data[target_path][i,2]+' & '+pwu(data[target_path][i,3],data[target_path][i,4])+' & '+pwu(data[target_path][i,5],data[target_path][i,6])+' & '+pwu(data[target_path][i,7],data[target_path][i,8])+'\\\  \n'
	text+= '\hline \n'
	text+= 'avg & '+pwu_list(mean_sem_disp([data[target_path][0,0],data[target_path][1,0]],[data[target_path][0,1],data[target_path][1,1]]))+' & '+ '%1.6f'%data[target_path][i,2]+' & '+pwu_list(mean_sem_disp([data[target_path][0,3],data[target_path][1,3]],[data[target_path][0,4],data[target_path][1,4]]))+' & '+pwu_list(mean_sem_disp([data[target_path][0,5],data[target_path][1,5]],[data[target_path][0,6],data[target_path][1,6]]))+' & '+pwu_list(mean_sem_disp([data[target_path][0,7],data[target_path][1,7]],[data[target_path][0,8],data[target_path][1,8]]))+'\\\  \n'
	text+= '\hline \n'
	text+= '\hline \n'
text+= '\hline \n' 

text+= '& \multicolumn{5}{c}{cdc1200 (KCl)} \\\ \n' 
text+= '\hline \n' 
for target_path in paths_to_launch_kcl_1200:
	if '8800_160' in target_path:
		conc=1.0
	else:
		conc=0.5
	for i in [0,1]:
		text+='el %d, %1.1f M & '%(i+1, conc)+ pwu(data[target_path][i,0],data[target_path][i,1])+' & '+ '%1.6f'%data[target_path][i,2]+' & '+pwu(data[target_path][i,3],data[target_path][i,4])+' & '+pwu(data[target_path][i,5],data[target_path][i,6])+' & '+pwu(data[target_path][i,7],data[target_path][i,8])+'\\\  \n'
	text+= '\hline \n'
	text+= 'avg & '+pwu_list(mean_sem_disp([data[target_path][0,0],data[target_path][1,0]],[data[target_path][0,1],data[target_path][1,1]]))+' & '+ '%1.6f'%data[target_path][i,2]+' & '+pwu_list(mean_sem_disp([data[target_path][0,3],data[target_path][1,3]],[data[target_path][0,4],data[target_path][1,4]]))+' & '+pwu_list(mean_sem_disp([data[target_path][0,5],data[target_path][1,5]],[data[target_path][0,6],data[target_path][1,6]]))+' & '+pwu_list(mean_sem_disp([data[target_path][0,7],data[target_path][1,7]],[data[target_path][0,8],data[target_path][1,8]]))+'\\\  \n'
	text+= '\hline \n'
	text+= '\hline \n'

text+= '\hline \n' 

print text
##########################################################################################


text= '\hline \n' 
text +=' & '+' $R_l \;[10^7\cdot \Omega]$ &  $R_{bulk}  \;[10^7\cdot \Omega]$ & Q$_{max1}$ $[e]$ & Q$_{max2}$ $[e]$ & C$_{tot}$ '+r' $[\frac{F}{g}]$ \\'+ '\n' 
text+= '\hline \n' 
text+= '& \multicolumn{5}{c}{cdc1200 (NaCl)} \\\ \n' 

for target_path in paths_to_launch_nacl_1200:
	print '\hline' 
	if '8800_160' in target_path:
		conc=1.0
	else:
		conc=0.5
	print '\multicolumn{2}{c}{cdc1200, NaCl, 1.0V, %1.2f M} \\\ ' %(conc)
	print '\hline' 
	print_values(target_path)

for target_path in paths_to_launch_nacl_800:
	if '7615_139' in target_path:
		conc=1.0
	else:
		conc=0.5
	print '\hline' 
	print '\multicolumn{2}{c}{cdc800, NaCl, 1.0V, %1.2f M} \\\ ' %(conc)
	print '\hline' 
	print_values(target_path)


for target_path in paths_to_launch_kcl_1200:
	if '8800_160' in target_path:
		conc=1.0
	else:
		conc=0.5
	print '\hline' 
	print '\multicolumn{2}{c}{cdc1200, KCl, 1.0V, %1.2f M} \\\ ' %(conc)
	print '\hline' 
	print_values(target_path)

print '\hline' 



##########################################################################################















'''



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
'''