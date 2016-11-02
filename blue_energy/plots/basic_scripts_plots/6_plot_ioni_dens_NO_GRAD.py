import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os

from sys import argv

global n_ions
############################################################
from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################
conv=0.529177
max_dist=186.406 #length of the whole system

end_name='ERR_EnD_Name_missing'

end_name=str(argv[1])
#list of directories in which we have to launch the script
list_sub_dir_here=os.walk('.').next()[1]
list_sub_dir_here=['/'+s for s in list_sub_dir_here]

############################################
data_path={}
num_confs={}
#load all the data in memory
name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_path=name_path.split('/')
name_check=name_path.split('results')[1]
name_check=name_check.split('6_acc')[0]

#print name_check

common_parent_path='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing'
base_dir_script=common_parent_path+'/plots/plot_low_level_scripts/'


file_frames=open(base_dir_script+'frames_infos.dat', 'r') 
file_frames_lines= file_frames.readlines()

num_confs=[]

for line in file_frames_lines:
	if name_check in line:
		num_confs.append(int((line.split('n_frames=')[1]).split('\n')[0]))

#print name_path
os.chdir(name_path)
file_list=os.listdir(os.getcwd())
#print file_list


data_num=[]
for name_dat_file in file_list:
	if ('surface_density_el_1_run'in name_dat_file):
		data_num.append(int(name_dat_file.split('run')[1][:3]))

#print data_num
runmin=0

runmax=max(data_num)

x_coord=np.loadtxt('surface_density_el_1_run000.out')[:,0]
x_coord*=conv


set_list_bound=[[runmin,runmax],[runmin,runmax/2],[runmax/2+1,runmax]]
for list_bound in set_list_bound: 
	data_path_el_1_Na=np.loadtxt('surface_density_el_1_run%03d.out'%(list_bound[0]))[:,2]*num_confs[list_bound[0]]
	data_path_el_1_Cl=np.loadtxt('surface_density_el_1_run%03d.out'%(list_bound[0]))[:,3]*num_confs[list_bound[0]]
	data_path_el_2_Na=np.loadtxt('surface_density_el_2_run%03d.out'%(list_bound[0]))[:,2]*num_confs[list_bound[0]]
	data_path_el_2_Cl=np.loadtxt('surface_density_el_2_run%03d.out'%(list_bound[0]))[:,3]*num_confs[list_bound[0]]
	for irun in range(list_bound[0]+1,list_bound[1]+1):
		data_path_el_1_Na+=(np.loadtxt('surface_density_el_1_run%03d.out'%(irun)))[:,2]*num_confs[irun]
		data_path_el_1_Cl+=(np.loadtxt('surface_density_el_1_run%03d.out'%(irun)))[:,3]*num_confs[irun]
		data_path_el_2_Na+=(np.loadtxt('surface_density_el_2_run%03d.out'%(irun)))[:,2]*num_confs[irun]
		data_path_el_2_Cl+=(np.loadtxt('surface_density_el_2_run%03d.out'%(irun)))[:,3]*num_confs[irun]

	#wheighted average
	somma=sum(num_confs[list_bound[0]:list_bound[1]+1])
	normalize=1.0/somma
	#print normalize
	data_path_el_1_Na*=normalize
	data_path_el_2_Na*=normalize
	data_path_el_1_Cl*=normalize
	data_path_el_2_Cl*=normalize

	#you are going from units of 1/Bohr^3 to 1/Angstrom^3
	data_path_el_1_Na/=(conv**3)
	data_path_el_2_Na/=(conv**3)
	data_path_el_1_Cl/=(conv**3)
	data_path_el_2_Cl/=(conv**3)	
	

	fig=pl.figure()
	ax1=fig.add_subplot(111)

	if list_path[9]=='nacl':
		name_ion='NaCl'
		name_label_ion_1='Na'
	else:
		name_ion='KCl'
		name_label_ion_1='K'

	name_file_save=r'6_ion_dens_%03d_%03d_%s_%s_%s_%s' %(list_bound[0],list_bound[1], list_path[8],name_ion,list_path[10],list_path[11] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+=end_name
	name_file_save+='.pdf'	
	#print name_file_save
	dim_probes=((int(end_name.split('_')[0])+0.01*int(end_name.split('_')[1]))/2.0)*conv
	info_probes=r'd(C-%s)='%(end_name.split('_')[-1])+'%1.2f'%(dim_probes)+r' $\AA$'
	title='Ionic_density_%03d_%03d_%s_%s_%s_%s_%s - no gradient' %(list_bound[0],list_bound[1], list_path[8],name_ion,list_path[10],list_path[11],info_probes)
	title=title.replace('_', ' ')
	ax1.set_title(title)
	#print title
	ax1.grid(True)
	ax1.set_xlabel(r'distance $[\AA]$')
	ax1.set_ylabel(r'Ionic density $\left[\frac{1}{\AA^3}\right]$')
	ax1.set_xlim([0, max(x_coord)])	
	ax1.set_xlim([0, 8])

	if (list_path[11]=='0.0V'):
		ax1.plot(x_coord,data_path_el_1_Na, label='el 1 '+name_label_ion_1,c='orange')
		ax1.plot(x_coord,data_path_el_1_Cl, label='el 1 Cl',c='c')	
		ax1.plot(x_coord,data_path_el_2_Na, label='el 2 '+name_label_ion_1,c='red')
		ax1.plot(x_coord,data_path_el_2_Cl, label='el 2 Cl',c='b')	
		ax1.legend(loc='best')
	else:
		ax2 = ax1.twinx()
		ax2.set_xlim([0, 8])
		ax2.set_ylabel(r'Ionic density $\left[\frac{1}{\AA^3}\right]$')
		ax1.plot(x_coord,data_path_el_1_Na, label='el 1 '+name_label_ion_1,c='orange')
		ax2.plot(x_coord,data_path_el_1_Cl, label='el 1 Cl',c='c')	
		ax2.plot(x_coord,data_path_el_2_Na, label='el 2 '+name_label_ion_1,c='red')
		ax1.plot(x_coord,data_path_el_2_Cl, label='el 2 Cl',c='b')	
		ax1.legend(loc=9)
		ax2.legend(loc=1)
		pl.gcf().subplots_adjust(right=0.87)

	#pl.plot([0, 25], [1, 1], '--', lw=1, color='black')	

	max_1=max(data_path_el_1_Na)
	max_2=max(data_path_el_1_Cl)
	max_tot=max([max_1,max_2])
	h_line=1+max_tot*0.025
	#print h_line
	#pl.show()
	fig.savefig(name_file_save, dpi=300)


x_coord=np.loadtxt('surface_probability_el_1_run000.out')[:,0]
x_coord*=conv


set_list_bound=[[runmin,runmax],[runmin,runmax/2],[runmax/2+1,runmax]]
for list_bound in set_list_bound: 
	data_path_el_1_Na=np.loadtxt('surface_probability_el_1_run%03d.out'%(list_bound[0]))[:,2]*num_confs[list_bound[0]]
	data_path_el_1_Cl=np.loadtxt('surface_probability_el_1_run%03d.out'%(list_bound[0]))[:,3]*num_confs[list_bound[0]]
	data_path_el_2_Na=np.loadtxt('surface_probability_el_2_run%03d.out'%(list_bound[0]))[:,2]*num_confs[list_bound[0]]
	data_path_el_2_Cl=np.loadtxt('surface_probability_el_2_run%03d.out'%(list_bound[0]))[:,3]*num_confs[list_bound[0]]
	for irun in range(list_bound[0]+1,list_bound[1]+1):
		data_path_el_1_Na+=(np.loadtxt('surface_probability_el_1_run%03d.out'%(irun)))[:,2]*num_confs[irun]
		data_path_el_1_Cl+=(np.loadtxt('surface_probability_el_1_run%03d.out'%(irun)))[:,3]*num_confs[irun]
		data_path_el_2_Na+=(np.loadtxt('surface_probability_el_2_run%03d.out'%(irun)))[:,2]*num_confs[irun]
		data_path_el_2_Cl+=(np.loadtxt('surface_probability_el_2_run%03d.out'%(irun)))[:,3]*num_confs[irun]

	#wheighted average
	somma=sum(num_confs[list_bound[0]:list_bound[1]+1])
	normalize=1.0/somma
	#print normalize
	data_path_el_1_Na*=(normalize)
	data_path_el_2_Na*=(normalize)
	data_path_el_1_Cl*=(normalize)
	data_path_el_2_Cl*=(normalize)	
	dx=x_coord[1]-x_coord[0]
	#I am going in the continuous limit: from the probability distribution I compute the probability density function!
	data_path_el_1_Na/=dx
	data_path_el_2_Na/=dx
	data_path_el_1_Cl/=dx
	data_path_el_2_Cl/=dx	
	

	fig=pl.figure()
	ax1=fig.add_subplot(111)
	name_file_save=r'6_ion_prob_%03d_%03d_%s_%s_%s_%s' %(list_bound[0],list_bound[1], list_path[8],name_ion,list_path[10],list_path[11] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+=end_name
	name_file_save+='.pdf'	
	#print name_file_save
	title='Ionic_probability_%03d_%03d_%s_%s_%s_%s_%s - no gradient' %(list_bound[0],list_bound[1], list_path[8],name_ion,list_path[10],list_path[11],info_probes)
	title=title.replace('_', ' ')
	ax1.set_title(title)
	#print title
	ax1.grid(True)
	ax1.set_xlabel(r'distance $[\AA]$')
	ax1.set_ylabel(r'Probability density $\left[\frac{1}{\AA}\right]$')
	ax1.set_xlim([0, max(x_coord)])
	ax1.set_xlim([0, 8])	
	area=0
	
	#for index in range(0,len(data_path_el_1_Na)):
	#	area+=dx*data_path_el_1_Cl[index]*0.01
#
	#print 'area=', area

	ax1.plot(x_coord,data_path_el_1_Na*0.01, label='el 1 '+name_label_ion_1,c='orange')
	ax1.plot(x_coord,data_path_el_1_Cl*0.01, label='el 1 Cl',c='c')	
	ax1.plot(x_coord,data_path_el_2_Na*0.01, label='el 2 '+name_label_ion_1,c='red')
	ax1.plot(x_coord,data_path_el_2_Cl*0.01, label='el 2 Cl',c='b')

	ax1.legend(loc='best')
	#pl.plot([0, 25], [1, 1], '--', lw=1, color='black')	

	max_1=max(data_path_el_1_Na)
	max_2=max(data_path_el_1_Cl)
	max_tot=max([max_1,max_2])
	h_line=1+max_tot*0.025
	#print h_line
	#pl.show()
	fig.savefig(name_file_save, dpi=300)