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

def build_histo(array,histo_O,histo_C,histo_Na_Cl,Nbins):
	for k in range(0,len(array[:,0])):
		histo_pos=math.ceil((array[k,0]/max_dist)*(Nbins-1)) #since floor 
		if (k<n_ions/2):
			histo_O[histo_pos,0,0]+=array[k,1]
			if (array[k,1]!=0):
				histo_O[histo_pos,0,1]+=1 #count how many atoms you put!
			histo_Na_Cl[histo_pos,0]+=array[k,2] #NaCl or ClNa is the same (taking into account how you did calculations)
			if (array[k,2]!=0):
				histo_Na_Cl[histo_pos,1]+=1
			histo_C[histo_pos,0,0]+=array[k,3]
			if (array[k,3]!=0):
				histo_C[histo_pos,0,1]+=1
		else:
			histo_O[histo_pos,1,0]+=array[k,1]
			if (array[k,1]!=0):
				histo_O[histo_pos,1,1]+=1
			histo_C[histo_pos,1,0]+=array[k,3]
			if (array[k,3]!=0):
				histo_C[histo_pos,1,1]+=1

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


print max_time_dict

numbers={}
data={}
x={}
y={}
indx=0
num_fram=[]

conv=0.529177
list_colors=['#0000FF','#01DF01', '#FF8000' ]

#list of directories in which we have to launch the script
list_sub_dir_here=os.walk('.').next()[1]
list_sub_dir_here=['/'+s for s in list_sub_dir_here]

data_avg_O={}
data_avg_C={}
data_avg_Ions={}

for datadir in list_sub_dir_here:
	name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
	print name_path
	list_path=name_path.split('/')
	target_path=name_path.split('/results')[1]
	target_path=target_path.split('/4')[0]

	if 'cdc800' in target_path:
		max_dist=183.000 #length of the whole system
	else:
		max_dist=186.406 #length of the whole system


	if list_path[9]=='nacl':
		name_ion='NaCl'
		name_label_ion_1='Na'
	else:
		name_ion='KCl'
		name_label_ion_1='K'

	if ((list_path[10]=='8800_160') or (list_path[10]=='7615_139')):
		conC='1.0 M'
	else:
		conC='0.5 M'

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

	array_avg_O={}
	array_avg_C={}
	array_avg_Ions={}

	for bound_idx in range(0,len(list_bounds)):
		runmin=list_bounds[bound_idx][0]
		runmax=list_bounds[bound_idx][1]
		#load all the data in memory
		print runmin, runmax
		for irun in range(runmin,runmax+1):
			name_file_load=r'run%03d_out_coord_num.dat' %(irun)
			data_path=r'.'+datadir+r'/'+'run%03d_out_coord_num.dat'%(irun)
			data[indx] = np.loadtxt(data_path, usecols = (4,5,7,8))
			data[indx][:,0]*=conv #convert from atomic units to Angstrom
			num_lin=len(data[indx][:,0])
			num_fram.append(num_lin/n_ions)
			indx+=1


		n_cols=len(data[0][0,:])
		array_avg_O=np.zeros((n_bins,2,2)) #2=number of different ions: index 0=Na, index 1=Cl
		array_avg_C=np.zeros((n_bins,2,2))
		array_avg_Ions=np.zeros((n_bins,2))	


		for i in range(0,indx): #loop over runs
			for j in range(0,num_fram[i]):
				build_histo(data[i][j*n_ions:(j+1)*n_ions,:],array_avg_O,array_avg_C,array_avg_Ions,n_bins)		

		#perform the average
		for pos in range(0,n_bins):
			if (array_avg_Ions[pos,1]!=0):
				array_avg_Ions[pos,0]*=(1.0/array_avg_Ions[pos,1])
			for j in range(0,2):
				if (array_avg_O[pos,j,1]!=0):
					array_avg_O[pos,j,0]*=(1.0/array_avg_O[pos,j,1])
				if (array_avg_C[pos,j,1]!=0):
					array_avg_C[pos,j,0]*=(1.0/array_avg_C[pos,j,1])

		data_avg_O[bound_idx]=np.copy(array_avg_O)
		data_avg_C[bound_idx]=np.copy(array_avg_C)
		data_avg_Ions[bound_idx]=np.copy(array_avg_Ions)
			
	


	fig=pl.figure()
	ax1 = fig.add_subplot(211)
	
	run_med=indx/2-1
	name_file_save=r'4_Coord_Num_%03d_%03d_%s_%s_%s_%s_compare' %(0,dict_bound_3[target_path], list_path[8],list_path[9],list_path[10],list_path[11] )
	name_file_save= name_file_save.replace('.','_')
	
	title= r'C.N. %s %s %s %s' %(list_path[8],name_ion ,conC,list_path[11])
	title=title.replace('_', ' ')
	ax1.set_title(title)
	ax1.grid(True)
	# ax1.set_xlim([-50,50])
	
	bar_width=0.5*max_dist/n_bins
	coor_z=[]
	for x in xrange(0,n_bins):
		coor_z.append(bar_width+max_dist*x/n_bins)	

	ax2 = fig.add_subplot(212, sharex=ax1)
	ax2.grid(True)
	ax1.set_ylabel(r'coordination number')
	ax2.set_ylabel(r'coordination number')
	for bound_idx in range(0,len(list_bounds)):
		for i in range(0,len(data_avg_C[bound_idx][:,0,1])):
			if (data_avg_C[bound_idx][i,0,1]<0.1):
				data_avg_C[bound_idx][i,0,1]=float('Inf')	

		ax1.errorbar(coor_z,data_avg_C[bound_idx][:,0,0],
	     			 color=list_colors[bound_idx],
	                 yerr=np.reciprocal(np.sqrt(data_avg_C[bound_idx][:,0,1])),
	                 xerr=bar_width,
	                 label=ion_type_1+'-C '+list_names_legend[bound_idx])	
	


		for i in range(0,len(data_avg_C[bound_idx][:,1,1])):
			if (data_avg_C[bound_idx][i,1,1]<0.1):
				data_avg_C[bound_idx][i,1,1]=float('Inf')	


		ax2.set_xlabel(r'z coordinate $[\AA]$')
		ax2.errorbar(coor_z,data_avg_C[bound_idx][:,1,0],
	     			 color=list_colors[bound_idx],
	                 yerr=np.reciprocal(np.sqrt(data_avg_C[bound_idx][:,1,1])),
	                 xerr=bar_width,
	                 label='Cl-C '+list_names_legend[bound_idx])	


	ax1.legend(loc=9,prop={'size':11})
	ax2.legend(loc=9,prop={'size':11})

	max_totx1=max(data_avg_C[bound_idx][:,0,0])
	max_totx2=max(data_avg_C[bound_idx][:,1,0])	
	ax1.set_ylim(0,round(1.2*max_totx1,0))
	ax2.set_ylim(0,round(1.2*max_totx2,0))

	fig.tight_layout()
	fig.set_size_inches(6, 4.8)
	fig.savefig(name_file_save+'ion-C.pdf', dpi=150)




	figA=pl.figure()
	ax3 = figA.add_subplot(211)
	title= r'C.N. %s %s %s %s' %(list_path[8],name_ion, conC, list_path[11])
	title=title.replace('_', ' ')
	ax3.set_title(title)
	ax3.grid(True)
	ax4 = figA.add_subplot(212, sharex=ax3)
	ax3.set_ylabel(r'coordination number')
	ax4.set_ylabel(r'coordination number')
	ax4.set_xlabel(r'z coordinate $[\AA]$')
	ax4.grid(True)

	for bound_idx in range(0,len(list_bounds)):
		for i in range(0,len(data_avg_O[bound_idx][:,0,1])):
			if (data_avg_O[bound_idx][i,0,1]<0.1):
				data_avg_O[bound_idx][i,0,1]=float('Inf')	

		ax3.errorbar(coor_z,data_avg_O[bound_idx][:,0,0],
	     			 color=list_colors[bound_idx],
	                 yerr=np.reciprocal(np.sqrt(data_avg_O[bound_idx][:,0,1])),
	                 xerr=bar_width,
	                 label=ion_type_1+'-O '+list_names_legend[bound_idx])	
	


		for i in range(0,len(data_avg_O[bound_idx][:,1,1])):
			if (data_avg_O[bound_idx][i,1,1]<0.1):
				data_avg_O[bound_idx][i,1,1]=float('Inf')	


		ax4.errorbar(coor_z,data_avg_O[bound_idx][:,1,0],
	     			 color=list_colors[bound_idx],
	                 yerr=np.reciprocal(np.sqrt(data_avg_O[bound_idx][:,1,1])),
	                 xerr=bar_width,
	                 label='Cl-O '+list_names_legend[bound_idx])	

	ax3.legend(loc=8,prop={'size':11})
	ax4.legend(loc=8,prop={'size':11})
	max_totx3=max(data_avg_O[bound_idx][:,0,0])
	max_totx4=max(data_avg_O[bound_idx][:,1,0])	
	ax3.set_ylim(0,round(1.2*max_totx3,0))
	ax4.set_ylim(0,round(1.2*max_totx4,0))

	figA.tight_layout()
	figA.set_size_inches(6, 4.8)
	figA.savefig(name_file_save+'ion-O.pdf', dpi=150)
	
##################################################################################
'''
	fig=pl.figure()
	ax1 = fig.add_subplot(111)
	name_file_save=r'4_Coord_Num_%03d_%03d_%s_%s_%s_%s_dataset_%s' %(0,dict_bound_3[target_path], list_path[8],list_path[9],list_path[10],list_path[11] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+='.pdf'
	title= r'Coordination number run%03d-%03d %s %s %s %s' %(0,dict_bound_3[target_path], list_path[8],list_path[9],list_path[10],list_path[11])
	title=title.replace('_', ' ')
	ax1.set_title(title)
	ax1.grid(True)
	# ax1.set_xlim([-50,50])
	ax1.set_xlabel(r'z coordinate $[\AA]$')
	ax1.set_ylabel(r'coordination number')
	max_1=max(array_avg_O_tot[:,0,0])
	max_2=max(array_avg_O_tot[:,1,0])
	max_3=max(array_avg_C_tot[:,0,0])
	max_4=max(array_avg_C_tot[:,1,0])
	max_5=max(array_avg_Ions_tot[:,0])
	max_ax_tot=max([max_1,max_2,max_3,max_4,max_5])
	ax1.set_ylim(0, max_ax_tot*1.1)	

	coor_z=[]
	for x in xrange(0,n_bins):
		coor_z.append(max_dist*x/n_bins)	

	ax1.plot(coor_z,array_avg_O_tot[:,0,0], label=ion_type_1+'-O' ,color='r')	
	ax1.plot(coor_z,array_avg_O_tot[:,1,0], label='Cl-O',c='blue' )
	ax1.plot(coor_z,array_avg_C_tot[:,0,0], label=ion_type_1+'-C',c='black' )	
	ax1.plot(coor_z,array_avg_C_tot[:,1,0], label='Cl-C',c='green' )
	#ax1.plot(coor_z,array_avg_Ions_tot, label=ion_type_1+'-Cl',c='m' )
	art = []
	lgd = ax1.legend(loc=9, bbox_to_anchor=(1.2, 1), ncol=1)
	art.append(lgd)	

	fig.savefig(name_file_save, dpi=300,additional_artists=art, bbox_inches="tight")
'''