import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os

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

common_parent_path='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing'
def build_histo(array_1,histo):
	#units in Angstrom directly converted in atomic units
	electrode_1	= [ 17./conv, 40./conv ]
	#bulk =[  55./conv , 131.406/conv]
	electrode_2 = [ 146.406/conv , 169.406/conv]
	for k in range(0,len(array_1[:,0])):
		if (array_1[k,1] !=0):
			if (array_1[k,0]>electrode_1[0]) and (array_1[k,0]<electrode_1[1]):
				histo[array_1[k,1],0]+=1		
			#elif (array_1[k,0]>bulk[0]) and (array_1[k,0]<bulk[1]):
				#histo[array_1[k,1],1]+=1		
			elif (array_1[k,0]>electrode_2[0]) and (array_1[k,0]<electrode_2[1]):
				histo[array_1[k,1],2]+=1


numbers={}
data={}
x={}
y={}
indx=0
num_fram=[]

conv=0.529177
max_dist=186.406 #length of the whole system






#list of directories in which we have to launch the script
list_sub_dir_here=os.walk('.').next()[1]
list_sub_dir_here=['/'+s for s in list_sub_dir_here]

for datadir in list_sub_dir_here:
	name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
	list_path=name_path.split('/')

	ion_name_list=list_path[9].split('c')
	ion_type_1=ion_name_list[0].title()
	ion_type_2='Cl'
	ints_list=[int(s) for s in datadir.split('_') if s.isdigit()]
	runmin=ints_list[0]
	runmax=ints_list[1]

	n_ions=int(list_path[10].split('_')[1])*2

	#load all the data in memory
	for irun in range(runmin,runmax+1):
		name_file_load=r'run%03d_out_coord_num.dat' %(irun)
		data_path=r'.'+datadir+r'/'+'run%03d_out_coord_num.dat'%(irun)
		##################################################################
		#FIND AND REPLACE IN A FILE: THIS IS NO LONGER NECESSARY BUT THIS IS USEFUL
		#s = open(data_path).read()
		#s = s.replace('Run ', '#Run_')
		#s = s.replace('label', '#labl')
		#f = open(data_path, 'w')
		#f.write(s)
		#f.close()
		##################################################################

		data[indx] = np.loadtxt(data_path, usecols = (4,8))
		#col 4 = z_coordinate
		#col 5 = coordination ion-O
		#col 7 = coordination Na-Cl ->here use just 1st 160 rows!!
		#col 8 = coordination ion-C	
	

		#print data[indx]
		num_lin=len(data[indx][:,0])
		num_fram.append(num_lin/n_ions)
		indx+=1


	histo_Na_C=np.zeros((n_bins,3)) 
	histo_Cl_C=np.zeros((n_bins,3))
	
	for i in range(0,indx):
		for j in range(0,num_fram[i]):
			build_histo(data[i][j*n_ions:j*n_ions+(n_ions/2),:],histo_Na_C)
			#print len(data[i][j*n_ions:j*n_ions+(n_ions/2),0])
			build_histo(data[i][j*n_ions+(n_ions/2):(j+1)*n_ions,:],histo_Cl_C)
			#print len(data[i][j*n_ions+(n_ions/2):(j+1)*n_ions,0])

	norm_el_1_Na=sum(histo_Na_C[:,0])
	#norm_bulk_Na=sum(histo_Na_C[:,1])
	norm_el_2_Na=sum(histo_Na_C[:,2])

	norm_el_1_Cl=sum(histo_Cl_C[:,0])
	#norm_bulk_Cl=sum(histo_Cl_C[:,1])
	norm_el_2_Cl=sum(histo_Cl_C[:,2])

	histo_Na_C[:,0]*=1.0/norm_el_1_Na
	#histo_Na_C[:,1]*=1.0/norm_bulk_Na
	histo_Na_C[:,2]*=1.0/norm_el_2_Na

	histo_Cl_C[:,0]*=1.0/norm_el_1_Cl
	#histo_Cl_C[:,1]*=1.0/norm_bulk_Cl
	histo_Cl_C[:,2]*=1.0/norm_el_2_Cl
	
	#start plotting data
	fig=pl.figure()
	ax1 = fig.add_subplot(111)
	name_file_save=r'4_Coord_Num_PDF_%03d_%03d_%s_%s_%s_%s_dataset_%s_Na' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+='.pdf'
	title= r'Coordination number run%03d-%03d %s %s %s %s dataset %s' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:])
	title=title.replace('_', ' ')
	ax1.set_title(title)
	ax1.grid(True)
	# ax1.set_xlim([-50,50])
	ax1.set_xlabel(r'coordination number')
	ax1.set_ylabel(r'probability')
	max_1=max(histo_Na_C[:,0])
	max_2=max(histo_Na_C[:,2 ])
	max_tot=max([max_1,max_2])
	ax1.set_ylim(0,round(1.1*max_tot,2))	
	for j in range(1,len(histo_Na_C[:,0])+1):
		if((histo_Na_C[-j,0]!=0) or(histo_Na_C[-j,2]!=0)):
			max_x=len(histo_Cl_C[:,0])-j
			break

	ax1.set_xlim(0,max_x)
	


	coor_x=np.arange(0,len(histo_Cl_C[:,0]))
	bin_size=coor_x[1]-coor_x[0]	
	ax1.plot(coor_x,histo_Na_C[:,0], label=ion_type_1+'-C el. 1',color='b')	
	ax1.bar(coor_x,histo_Na_C[:,0], color='c', alpha=0.5, align='center')	
	#ax1.plot(coor_x,histo_Na_C[:,1], label=ion_type_1+'-C bulk',color='orange')	
	ax1.plot(coor_x,histo_Na_C[:,2], label=ion_type_1+'-C el. 2',color='r')
	ax1.bar(coor_x,histo_Na_C[:,2], color='orange', alpha=0.5,align='center')
	ax1.legend(loc='best')
	fig.savefig(name_file_save, dpi=300)
###########################################################################################
	fig=pl.figure()
	ax1 = fig.add_subplot(111)
	name_file_save=r'4_Coord_Num_PDF_%03d_%03d_%s_%s_%s_%s_dataset_%s_Cl' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+='.pdf'
	title= r'Coordination number run%03d-%03d %s %s %s %s dataset %s' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:])
	title=title.replace('_', ' ')
	ax1.set_title(title)
	ax1.grid(True)
	# ax1.set_xlim([-50,50])
	ax1.set_xlabel(r'coordination number')
	ax1.set_ylabel(r'probability')
	max_1_a=max(histo_Cl_C[:,0])
	max_2_a=max(histo_Cl_C[:,2])
	max_tot_a=max([max_1_a,max_2_a])
	ax1.set_ylim(0,round(1.1*max_tot_a,2))	
	for j in range(1,len(histo_Cl_C[:,0])+1):
		if((histo_Cl_C[-j,0]!=0) or(histo_Cl_C[-j,2]!=0)):
			max_x_a=len(histo_Cl_C[:,0])-j
			break
	ax1.set_xlim(0,max_x_a)

	#for j in range(0,len(histo_Cl_C[:,2])):
	#	print j, histo_Cl_C[j,2]


	coor_x=np.arange(0,len(histo_Cl_C[:,0]))	

	ax1.plot(coor_x,histo_Cl_C[:,0], label=ion_type_2+'-C el. 1',color='b')	
	ax1.bar(coor_x,histo_Cl_C[:,0], color='c', alpha=0.5, align='center')	
	
	#ax1.plot(coor_x,histo_Cl_C[,1], label=ion_type_2+'-C bulk',color='orange')	
	ax1.plot(coor_x,histo_Cl_C[:,2], label=ion_type_2+'-C el. 2',color='r')
	ax1.bar(coor_x,histo_Cl_C[:,2], color='orange', alpha=0.5,align='center')
	ax1.legend(loc='best')
	fig.savefig(name_file_save, dpi=300)

##################################################################################

	histo_Na_C_1=np.zeros((n_bins,3)) 
	histo_Cl_C_1=np.zeros((n_bins,3))

	histo_Na_C_2=np.zeros((n_bins,3)) 
	histo_Cl_C_2=np.zeros((n_bins,3))

	for i in range(0,indx/2):
		for j in range(0,num_fram[i]):
			build_histo(data[i][j*n_ions:j*n_ions+(n_ions/2),:],histo_Na_C_1)
			#print len(data[i][j*n_ions:j*n_ions+(n_ions/2),0])
			build_histo(data[i][j*n_ions+(n_ions/2):(j+1)*n_ions,:],histo_Cl_C_1)
			#print len(data[i][j*n_ions+(n_ions/2):(j+1)*n_ions,0])

	for i in range(indx/2,indx):
		for j in range(0,num_fram[i]):
			build_histo(data[i][j*n_ions:j*n_ions+(n_ions/2),:],histo_Na_C_2)
			#print len(data[i][j*n_ions:j*n_ions+(n_ions/2),0])
			build_histo(data[i][j*n_ions+(n_ions/2):(j+1)*n_ions,:],histo_Cl_C_2)
			#print len(data[i][j*n_ions+(n_ions/2):(j+1)*n_ions,0])

	histo_Na_C_1[:,0]*=1.0/sum(histo_Na_C_1[:,0])
	histo_Na_C_1[:,2]*=1.0/sum(histo_Na_C_1[:,2])

	histo_Na_C_2[:,0]*=1.0/sum(histo_Na_C_2[:,0])
	histo_Na_C_2[:,2]*=1.0/sum(histo_Na_C_2[:,2])

	histo_Cl_C_1[:,0]*=1.0/sum(histo_Cl_C_1[:,0])
	histo_Cl_C_1[:,2]*=1.0/sum(histo_Cl_C_1[:,2])

	histo_Cl_C_2[:,0]*=1.0/sum(histo_Cl_C_2[:,0])
	histo_Cl_C_2[:,2]*=1.0/sum(histo_Cl_C_2[:,2])

	fig=pl.figure()
	ax1 = fig.add_subplot(111)
	name_file_save=r'4_Coord_Num_PDF_%03d_%03d_%s_%s_%s_%s_dataset_%s_Na_compare' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+='.pdf'
	title= r'Coordination number run%03d-%03d %s %s %s %s dataset %s' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:])
	title=title.replace('_', ' ')
	ax1.set_title(title)
	ax1.grid(True)
	# ax1.set_xlim([-50,50])
	ax1.set_xlabel(r'coordination number')
	ax1.set_ylabel(r'probability')
	max_1=max(histo_Na_C_1[:,0])
	max_2=max(histo_Na_C_1[:,2 ])
	max_3=max(histo_Na_C_2[:,0])
	max_4=max(histo_Na_C_2[:,2 ])
	max_tot=max([max_1,max_2,max_3,max_4])
	ax1.set_ylim(0,round(1.1*max_tot,2))

	for j in range(1,len(histo_Na_C[:,0])+1):
		if((histo_Na_C_1[-j,0]!=0) or(histo_Na_C_1[-j,2]!=0) or (histo_Na_C_2[-j,0]!=0) or(histo_Na_C_2[-j,2]!=0)):
			max_x=len(histo_Cl_C_1[:,0])-j
			break

	ax1.set_xlim(0,max_x)


	coor_x=np.arange(0,len(histo_Cl_C[:,0]))	

	ax1.plot(coor_x,histo_Na_C_1[:,0], label=ion_type_1+'-C el. 1 1st part',color='b')		
	ax1.plot(coor_x,histo_Na_C_1[:,2], label=ion_type_1+'-C el. 2 1st part',color='r')
	ax1.plot(coor_x,histo_Na_C_2[:,0], label=ion_type_1+'-C el. 1 2nd part',color='#7b68ee')		
	ax1.plot(coor_x,histo_Na_C_2[:,2], label=ion_type_1+'-C el. 2 2nd part',color='orange')

	ax1.legend(loc='best')
	fig.savefig(name_file_save, dpi=300)
	################################################################################
	fig=pl.figure()
	ax1 = fig.add_subplot(111)
	name_file_save=r'4_Coord_Num_PDF_%03d_%03d_%s_%s_%s_%s_dataset_%s_Cl_compare' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:] )
	name_file_save= name_file_save.replace('.','_')
	name_file_save+='.pdf'
	title= r'Coordination number run%03d-%03d %s %s %s %s dataset %s' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11],datadir[-1:])
	title=title.replace('_', ' ')
	ax1.set_title(title)
	ax1.grid(True)
	# ax1.set_xlim([-50,50])
	ax1.set_xlabel(r'coordination number')
	ax1.set_ylabel(r'probability')
	max_1=max(histo_Cl_C_1[:,0])
	max_2=max(histo_Cl_C_1[:,2 ])
	max_3=max(histo_Cl_C_2[:,0])
	max_4=max(histo_Cl_C_2[:,2 ])
	max_tot=max([max_1,max_2,max_3,max_4])
	ax1.set_ylim(0,round(1.1*max_tot,2))

	for j in range(1,len(histo_Na_C[:,0])+1):
		if((histo_Cl_C_1[-j,0]!=0) or(histo_Cl_C_1[-j,2]!=0) or (histo_Cl_C_2[-j,0]!=0) or(histo_Cl_C_2[-j,2]!=0)):
			max_x=len(histo_Cl_C_1[:,0])-j
			break

	ax1.set_xlim(0,max_x)


	coor_x=np.arange(0,len(histo_Cl_C[:,0]))	

	ax1.plot(coor_x,histo_Cl_C_1[:,0], label=ion_type_2+'-C el. 1 1st part',color='b')		
	ax1.plot(coor_x,histo_Cl_C_1[:,2], label=ion_type_2+'-C el. 2 1st part',color='r')
	ax1.plot(coor_x,histo_Cl_C_2[:,0], label=ion_type_2+'-C el. 1 2nd part',color='#3BB9FF')		
	ax1.plot(coor_x,histo_Cl_C_2[:,2], label=ion_type_2+'-C el. 2 2nd part',color='orange')

	ax1.legend(loc='best')
	fig.savefig(name_file_save, dpi=300)

