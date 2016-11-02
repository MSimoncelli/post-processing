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

histo_PDF_Na={}
norm_Na={}
histo_PDF_Na['bulk']=np.zeros(30)
histo_PDF_Na['edge']=np.zeros(30)
histo_PDF_Na['plane']=np.zeros(30)
histo_PDF_Na['hollow']=np.zeros(30)

histo_PDF_Cl={}
norm_Cl={}
histo_PDF_Cl['bulk']=np.zeros(30)
histo_PDF_Cl['edge']=np.zeros(30)
histo_PDF_Cl['plane']=np.zeros(30)
histo_PDF_Cl['hollow']=np.zeros(30)
coor_x=np.arange(0,30)

def build_histo_PDF(histo_PDF,coord,theta):
	if (coord>=0 and coord<29):
		if (theta<1e-2):
			histo_PDF['bulk'][coord]+=1
		elif (theta<11.0):
			histo_PDF['edge'][coord]+=1
		elif (theta<39.0):
			histo_PDF['plane'][coord]+=1
		elif(theta<58.0)		
			histo_PDF['hollow'][coord]+=1 


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

		#load all the data in memory
	for irun in range(runmin,runmax+1):
		name_file_load=r'run%03d_out_coord_num.dat' %(irun)
		data_path=r'.'+datadir+r'/'+'run%03d_out_coord_num.dat'%(irun)
		coord= np.loadtxt(data_path,usecols=(5,)) #coordination by O
		theta= np.loadtxt(data_path,usecols=(9,))
		name= np.genfromtxt(data_path,dtype='str',usecols=(0,))
		for i in range(0,len(theta)):
			if (name[i]=='Na' or name[i]=='K'):
				build_histo_PDF(histo_PDF_Na,coord[i],theta[i])
			elif(name[i]=='Cl'):
				build_histo_PDF(histo_PDF_Cl,coord[i],theta[i])

	norm_Na['bulk']=sum(histo_PDF_Na['bulk'])
	norm_Na['edge']=sum(histo_PDF_Na['edge'])
	norm_Na['plane']=sum(histo_PDF_Na['plane'])
	norm_Na['hollow']=sum(histo_PDF_Na['hollow'])

	norm_Cl['bulk']=sum(histo_PDF_Cl['bulk'])
	norm_Cl['edge']=sum(histo_PDF_Cl['edge'])
	norm_Cl['plane']=sum(histo_PDF_Cl['plane'])
	norm_Cl['hollow']=sum(histo_PDF_Cl['hollow'])

	for region in histo_PDF_Na.keys():
		histo_PDF_Na[region]=histo_PDF_Na[region]*(1.0/(1.0*norm_Na[region]))
		histo_PDF_Cl[region]=histo_PDF_Cl[region]*(1.0/(1.0*norm_Cl[region]))

	for region in ['edge', 'plane', 'hollow']:
		fig1=pl.figure()
		ax1 = fig1.add_subplot(111)
		fig2=pl.figure()
		ax2 = fig2.add_subplot(111)
		if (list_path[9]=='nacl'):
			spc_title='NaCl'
		elif(list_path[9]=='kcl'):
			spc_title='KCl'

		name_file_save_Na=r'8_solvation_prob_%s_%s_%03d_%03d_%s_%s_%s_%s' %(region,ion_type_1,runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11] )
		name_file_save_Na= name_file_save_Na.replace('.','_')
		name_file_save_Na+='.pdf'
		name_file_save_Cl=r'8_solvation_prob_%s_%s_%03d_%03d_%s_%s_%s_%s' %(region,'Cl',runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11] )
		name_file_save_Cl= name_file_save_Cl.replace('.','_')
		name_file_save_Cl+='.pdf'
		title_1= r'Coordination number by Oxigen %s run%03d-%03d %s %s %s %s' %(ion_type_1,runmin,runmax, list_path[8],spc_title,list_path[10],list_path[11])
		title_1=title_1.replace('_', ' ')
		ax1.set_title(title_1)
		ax1.grid(True)
		ax1.set_xlabel(r'Coordination number by Oxigen')
		ax1.set_ylabel(r'Probability')

		title_2= r'Coordination number by Oxigen %s run%03d-%03d %s %s %s %s' %(ion_type_2,runmin,runmax, list_path[8],spc_title,list_path[10],list_path[11])
		title_2=title_2.replace('_', ' ')
		ax1.set_title(title_2)
		ax2.grid(True)
		ax2.set_xlabel(r'Coordination number by Oxigen')
		ax2.set_ylabel(r'Probability')

		ax1.bar(coor_x,histo_PDF_Na['bulk'], label='bulk',color='b')	
		ax1.plot(coor_x,histo_PDF_Na['bulk'],color='b')	
		ax1.bar(coor_x,histo_PDF_Na[region], label=region ,color='r')	
		ax1.plot(coor_x,histo_PDF_Na[region],color='r')

		ax2.bar(coor_x,histo_PDF_Cl['bulk'], label='bulk',color='b')	
		ax2.plot(coor_x,histo_PDF_Cl['bulk'],color='b')	
		ax2.bar(coor_x,histo_PDF_Cl[region], label=region ,color='r')	
		ax2.plot(coor_x,histo_PDF_Cl[region],color='r')
			

		ax1.legend(loc='best')
		fig1.savefig(name_file_save_Na, dpi=300)
		ax1.legend(loc='best')
		fig2.savefig(name_file_save_Cl, dpi=300)


	
