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
if (len(argv)<=1):
        print 'USAGE: plor_coord_num.py n_bins'
        n_bins=50   
else:
	n_bins=int(argv[1])

from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################

n_bins_histo_doc=50
step=100.0/n_bins_histo_doc
coor_x=np.arange(0,100,step)
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
common_parent_path='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing'

def build_histo_theta(array_1,coord,theta):
	#units in Angstrom directly converted in atomic units
	#cdc1200 and cdc800 have different zise, here it is enough to disc
	half_box_z=(186.406/2.0)/conv
	if (theta>0.0):
			histo_pos=int(theta/100.0*n_bins_histo_doc)
			if (coord<half_box_z):
				array_1[0,histo_pos]+=theta		
			#elif (array_1[k,0]>bulk[0]) and (array_1[k,0]<bulk[1]):
			#	histo[array_1[k,1],1]+=1		
			else: 
				array_1[1,histo_pos]+=theta

def get_run_max():
        data_runmax=open('/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/plots/current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])-1
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()


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

def say_me_time_max():
        base_dir_script=common_parent_path+'/plots/plot_low_level_scripts/'
        file_frames=open(base_dir_script+'frames_infos.dat', 'r') 
        file_frames_lines= file_frames.readlines()      
        num_confs=[]     

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
				num_frames[target_path]+=num_confs[target_path][irun] 

		for target_path in paths_to_launch:
			num_frames[target_path]=num_frames[target_path]*0.001 #divide by output frequency 

        return num_frames

interesting_time_ns=[1.5,4.0]
max_time_dict=say_me_time_max()
dict_bound_1=say_me_run(interesting_time_ns[0])
dict_bound_2=say_me_run(interesting_time_ns[1])
dict_bound_3=get_run_max()

numbers={}
data={}
x={}
y={}
indx=0
num_fram=[]
max_val_Na=[]
max_val_Cl=[]

conv=0.529177
list_colors_1=['g','#0174DF', '#08088A' ]
list_colors_2=['purple','#FF8000', '#FF0000' ]
#list of directories in which we have to launch the script
list_sub_dir_here=os.walk('.').next()[1]
list_sub_dir_here=['/'+s for s in list_sub_dir_here]

data_avg_O={}
data_avg_C={}
data_avg_Ions={}

for datadir in list_sub_dir_here:
	name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
	#print name_path
	list_path=name_path.split('/')
	target_path=name_path.split('/results')[1]
	target_path=target_path.split('/4')[0]

	ion_name_list=list_path[9].split('c')
	ion_type_1=ion_name_list[0].title()
	ion_type_2='Cl'
	ints_list=[int(s) for s in datadir.split('_') if s.isdigit()]
	runmin=0
	runmax=dict_runmax[target_path]

	n_ions=int(list_path[10].split('_')[1])*2
	if (interesting_time_ns[1]<max_time_dict[target_path]):
		list_bounds=[[0,dict_bound_1[target_path]],[dict_bound_1[target_path]+1,dict_bound_2[target_path]],[dict_bound_2[target_path]+1,dict_bound_3[target_path]] ]
		list_names_legend=['0.0-%1.1f ns'%(interesting_time_ns[0]),'%1.1f-%1.1f ns'%(interesting_time_ns[0],interesting_time_ns[1]), '%1.1f-%1.1f ns'%(interesting_time_ns[1],max_time_dict[target_path])]

	else:
		list_bounds=[[0,dict_bound_1[target_path]],[dict_bound_1[target_path]+1,dict_bound_2[target_path]]]
		list_names_legend=['0.0-%1.1f ns'%(interesting_time_ns[0]),'%1.1f-%1.1f ns'%(interesting_time_ns[0],max_time_dict[target_path])]

	array_histo_Na=np.zeros((2,n_bins_histo_doc))
	array_histo_Cl=np.zeros((2,n_bins_histo_doc))


	fig1=pl.figure()
	ax1 = fig1.add_subplot(111)
	fig2=pl.figure()
	ax2 = fig2.add_subplot(111)

	print list_path[9]
	if (list_path[9]=='nacl'):
		spc_title='NaCl'
	elif(list_path[9]=='kcl'):
		spc_title='KCl'

	if ((list_path[10]=='8800_160') or (list_path[10]=='7615_139')):
	    conC='1.0 M'
	else:
	    conC='0.5 M'

	name_file_save_Na=r'7_deg_conf_%s_%03d_%03d_%s_%s_%s_%s' %(ion_type_1,runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11] )
	name_file_save_Na= name_file_save_Na.replace('.','_')
	name_file_save_Na+='.pdf'
	title_1= r'Degree of confinement %s %s %s %s %s' %(ion_type_1, list_path[8],spc_title,conC,list_path[11])
	title_1=title_1.replace('_', ' ')
	ax1.set_title(title_1)
	ax1.grid(True)
	ax1.set_xlabel(r'Degree of confiment $[\%]$')
	ax1.set_ylabel(r'Fraction of ions')

	name_file_save_Cl=r'7_deg_conf_%s_%03d_%03d_%s_%s_%s_%s' %('Cl',runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11] )
	name_file_save_Cl= name_file_save_Cl.replace('.','_')
	name_file_save_Cl+='.pdf'
	title_2= r'Degree of confinement %s %s %s %s %s' %('Cl', list_path[8],spc_title,conC,list_path[11])
	title_2=title_2.replace('_', ' ')
	ax2.set_title(title_2)
	ax2.grid(True)
	ax2.set_xlabel(r'Degree of confiment $[\%]$')
	ax2.set_ylabel(r'Fraction of ions')


	for bound_idx in range(0,len(list_bounds)):
		runmin=list_bounds[bound_idx][0]
		runmax=list_bounds[bound_idx][1]
		#load all the data in memory
		print runmin, runmax
		for irun in range(runmin,runmax+1):
			name_file_load=r'run%03d_out_coord_num.dat' %(irun)
			data_path=r'.'+datadir+r'/'+'run%03d_out_coord_num.dat'%(irun)
			coord= np.loadtxt(data_path,usecols=(4,))
			theta= np.loadtxt(data_path,usecols=(9,))
			name= np.genfromtxt(data_path,dtype='str',usecols=(0,))
			for i in range(0,len(theta)):
				if (name[i]=='Na' or name[i]=='K'):
					build_histo_theta(array_histo_Na,coord[i],theta[i])
				elif(name[i]=='Cl'):
					build_histo_theta(array_histo_Cl,coord[i],theta[i])


		norm_el_1_Na=sum(array_histo_Na[0,:])
		norm_el_2_Na=sum(array_histo_Na[1,:])
		norm_el_1_Cl=sum(array_histo_Cl[0,:])
		norm_el_2_Cl=sum(array_histo_Cl[1,:])

		array_histo_Na[0,:]*=(1.0/(1.0*norm_el_1_Na))
		array_histo_Na[1,:]*=(1.0/(1.0*norm_el_2_Na))
		array_histo_Cl[0,:]*=(1.0/(1.0*norm_el_1_Cl))
		array_histo_Cl[1,:]*=(1.0/(1.0*norm_el_2_Cl))		

		arr_el_1_Na=array_histo_Na[0,:]
		arr_el_2_Na=array_histo_Na[1,:]
		arr_el_1_Cl=array_histo_Cl[0,:]
		arr_el_2_Cl=array_histo_Cl[1,:]

		ax1.plot(coor_x,arr_el_1_Na, label='el. - '+list_names_legend[bound_idx],color=list_colors_1[bound_idx])	
		ax1.plot(coor_x,arr_el_2_Na, label='el. + '+list_names_legend[bound_idx],color=list_colors_2[bound_idx])	
		ax2.plot(coor_x,arr_el_1_Cl, label='el. - '+list_names_legend[bound_idx],color=list_colors_1[bound_idx])	
		ax2.plot(coor_x,arr_el_2_Cl, label='el. + '+list_names_legend[bound_idx],color=list_colors_2[bound_idx])	


		max_val_Na.append(max(arr_el_1_Na))
		max_val_Na.append(max(arr_el_2_Na))
		max_val_Cl.append(max(arr_el_1_Cl))
		max_val_Cl.append(max(arr_el_2_Cl))

	max_tot_Na=max(max_val_Na)
	ax1.set_ylim(0,round(1.1*max_tot_Na,2))
	max_tot_Na=round(1.1*max_tot_Na,2)-0.02

	for j in range(1,len(arr_el_1_Na)+1):
		if((arr_el_1_Na[-j]!=0) or(arr_el_2_Na[-j]!=0)):
			max_x=len(arr_el_1_Na)-j
			break
	ax1.set_xlim(0,max_x)
	if (max_x<70):
		ax1.set_xlim(0,70)
	ax1.axvspan(0, 13, alpha=0.6, color='gray')
	ax1.axvspan(13, 39, alpha=0.2 , color='gray')
	ax1.axvspan(39, 58, alpha=0.6, color='gray')
	ax1.axvspan(58, 100, alpha=0.2 , color='gray')
	ax1.text(1, max_tot_Na, 'Edge')
	ax1.text(14, max_tot_Na, 'Plane')
	ax1.text(40, max_tot_Na, 'Hollow')
	ax1.text(59, max_tot_Na, 'Pocket')
	ax1.legend(loc=7)
	fig1.set_size_inches(6, 4.8)
	fig1.savefig(name_file_save_Na, dpi=300)

	max_tot_Cl=max(max_val_Cl[:])
	ax1.set_ylim(0,round(1.1*max_tot_Cl,2))
	max_tot_Cl=round(1.1*max_tot_Cl,2)-0.02

	for j in range(1,len(arr_el_1_Cl)+1):
		if((arr_el_1_Cl[-j]!=0) or(arr_el_2_Cl[-j]!=0)):
			max_x=len(arr_el_1_Cl)-j
			break
	ax2.set_xlim(0,max_x)
	if (max_x<70):
		ax2.set_xlim(0,70)
	ax2.axvspan(0, 13, alpha=0.6, color='gray')
	ax2.axvspan(13, 39, alpha=0.2 , color='gray')
	ax2.axvspan(39, 58, alpha=0.6, color='gray')
	ax2.axvspan(58, 100, alpha=0.2 , color='gray')
	ax2.text(1, max_tot_Cl, 'Edge')
	ax2.text(14, max_tot_Cl, 'Plane')
	ax2.text(40, max_tot_Cl, 'Hollow')
	ax2.text(59, max_tot_Cl, 'Pocket')
	ax2.legend(loc=7)
	fig2.set_size_inches(6, 4.8)
	fig2.savefig(name_file_save_Cl, dpi=300)
		

	

