#!/usr/bin/python
import string, re, struct, sys, math, os, inspect
import subprocess
from sys import argv
import numpy as np
startdir = os.getcwd() #gets the current working directory
conv=0.529177

def read_info_specie(line_num_atom, line_name, file_positions):
	line_specie = file_positions[line_num_atom].split()
	Nspec=int(line_specie[0])
	line_name_specie=file_positions[line_name].split()
	name_spec=line_name_specie[0]
	return Nspec, name_spec

def read_numframes():
	name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
	list_dir=name_folder.split('/')
	num_run=list_dir[-1]
	n_rn=num_run.split('n')[1]
	fin=int(n_rn)
	dirrun='../run%03d'%(fin+1)
	#print dirrun
	os.chdir(dirrun)
	file_box = open( 'restart.dat', 'r' )
	line1= file_box.readline().split() #returns a list containing the lines
	file_box.close()     
	numframes=int(line1[1])/1000
        dirrun='../run%03d'%(fin)
        os.chdir(dirrun)
	#print dirrun
	return numframes

def read_sim_box_infos():
	start_dir=os.getcwd()
	dirrun = '/data_prace/michele/prace/cdc/cdc1200/nacl/8800_160/0.0V/run000'
	os.chdir(dirrun)
	file_box = open( 'restart.dat', 'r' )
	lines_r= file_box.readlines() #returns a list containing the lines
	file_box.close()     
	num_line_file=len(lines_r)
	cell=np.zeros([3])
	cell[0]=float(lines_r[num_line_file-3])
	cell[1]=float(lines_r[num_line_file-2])
	cell[2]=float(lines_r[num_line_file-1])
	os.chdir(start_dir)
	return cell

def create_input_files_for_degree_conf():
	box_size=read_sim_box_infos()
	file_rntime=open ('runtime.inpt','r')
	file_runtime_a=file_rntime.readlines() # read all the file, it contains lots of useful things	
	n_part=[]
	spc=[]
	array_index = [7,12,17,22,27,32,37,42]	

	for index in array_index:
		Nspec,name_spec=read_info_specie(index,index-1,file_runtime_a)
		n_part.append(Nspec)
		spc.append(name_spec)
		#print 'num %s=%d' %(name_spec, Nspec)	

	for k in xrange(0,len(spc)):
		if ((spc[k]=='Na')or(spc[k]=='k')):
			n_ions_1=n_part[k]	

	for k in xrange(0,len(spc)):
		if (spc[k]=='Cl'):
			n_ions_2=n_part[k]	

	name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
	if 'cdc800' in name_folder:
		z_max_electrode_1=87.75148391	

	if 'cdc1200' in name_folder:
		z_max_electrode_1=87.94962366	

	text_1='positions.out\n' #1st line: name of the positions's file 
	text_1+='wallq.out\n'
	text_1+='%d\n'%(read_numframes()) #2nd line: number of configurations in the positions' file
	text_1+= '%f %f %f\n'%(box_size[0], box_size[1], box_size[2]) #4th line length of the box in the 3 cartesian directions (bohrs)
	text_1+= '1\n'#number of type of anions
	text_1+='%d -1.0\n' %(n_ions_2) #number of Anions of type 1=n_part_[3] and charge
	text_1+= '1\n'#number of type of anions
	text_1+='%d 1.0\n' %(n_ions_1) #number of Cations of type 2=n_part[4]
	text_1+='1\n' #Number of types of solvent molecules (> NtypeS)): only H2O
	text_1+='%d\n' %(n_part[0]+n_part[1]+n_part[2]) #approx the CM of H2O with O only, then be careful to discard H1 and H2 atoms!!
	#be careful and check if you effectively the program considers also the P carbons
	text_1+='%d\n' %(n_part[5]+n_part[6]+n_part[7])  #number of carbons
	text_1+='3\n'#number of cutoffs: not interested!!
	text_1+='1 2 %f\n'%(3.27/conv) #Rcutoff O-Na
	text_1+='1 3 %f\n'%(3.95/conv) #Rcutoff O-Cl
	text_1+='2 3 %f\n'%(3.80/conv) #Rcutoff Na-Cl
	text_1+='%f\n'%(3.5/conv) #Rcutoff Na-C
	text_1+='%f\n'%(4.5/conv) #Rcutoff Cl-C
	text_1+='%f %f\n'%(0, z_max_electrode_1+0.5)#zmin and zmax (> zmin,zmax)0.0 90.0 	
	tmpf=open('input_deg_conf_el_1.dat', 'w')
	tmpf.write(text_1)
	tmpf.close()	

	text_2='positions.out\n' #1st line: name of the positions's file 
	text_2+='wallq.out\n'
	text_2+='%d\n'%(read_numframes()) #2nd line: number of configurations in the positions' file
	text_2+= '%f %f %f\n'%(box_size[0], box_size[1], box_size[2]) #4th line length of the box in the 3 cartesian directions (bohrs)
	text_2+= '1\n'#number of type of anions
	text_2+='%d -1.0\n' %(n_ions_2) #number of Anions of type 1=n_part_[3] and charge
	text_2+= '1\n'#number of type of anions
	text_2+='%d 1.0\n' %(n_ions_1) #number of Cations of type 2=n_part[4]
	text_2+='1\n' #Number of types of solvent molecules (> NtypeS)): only H2O
	text_2+='%d\n' %(n_part[0]+n_part[1]+n_part[2]) #approx the CM of H2O with O only, then be careful to discard H1 and H2 atoms!!
	#be c2reful and check if you effectively the program considers also the P carbons
	text_2+='%d\n' %(n_part[5]+n_part[6]+n_part[7]) 
	text_2+='3\n'#number of cutoffs: not interested!!
	text_2+='1 2 %f\n'%(3.27/conv) #Rcutoff O-Na
	text_2+='1 3 %f\n'%(3.95/conv) #Rcutoff O-Cl
	text_2+='2 3 %f\n'%(3.80/conv) #Rcutoff Na-Cl
	text_2+='%f\n'%(3.5/conv) #Rcutoff Na-C
	text_2+='%f\n'%(4.5/conv) #Rcutoff Cl-C
	text_2+='%f %f\n'%(box_size[2]-(z_max_electrode_1+0.5), box_size[2])#zmin and zmax (> zmin,zmax)0.0 90.0 
	tmpf=open('input_deg_conf_el_2.dat', 'w')
	tmpf.write(text_1)
	tmpf.close()

#all the files will start from this path
common_parent_path='/data_prace/michele/prace/cdc'

path_cdc800_el_1 = '/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area/electrode_cdc80/surf_points_el_1.dat'
path_cdc800_el_2 = '/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area/electrode_cdc80/surf_points_el_2.dat'
path_cdc1200_el_1='/data_prace/michele/prace/cdc//analysis/scripts/6_acc_surf_area/electrode_cdc1200/surf_points_el_1.dat'
path_cdc1200_el_2='/data_prace/michele/prace/cdc//analysis/scripts/6_acc_surf_area/electrode_cdc1200/surf_points_el_2.dat'

path_cdc_800_list=[path_cdc800_el_1,path_cdc800_el_2]
path_cdc_1200_list=[path_cdc1200_el_1,path_cdc1200_el_2]
'''
FULL LIST OF PATHS. 
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
#/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area
script_folder_2='/7_degree_of_confinement' #pay attention at the arguments
script_name_2='/degree_conf' #extracts only the name of the script
#put in this list all the paths in which you want to launch script[i]
paths_to_launch=[
'/cdc1200/nacl/8800_160/0.0V',
'/cdc1200/nacl/8800_160/1.0V', 
'/cdc1200/nacl/8800_80/0.0V',
'/cdc1200/nacl/8800_80/1.0V',
'/cdc800/nacl/7615_139/0.0V',
'/cdc800/nacl/7615_139/1.0V',
'/cdc800/nacl/7700_70/0.0V',
'/cdc800/nacl/7700_70/1.0V'
]

runmin=0
runmax=56
num_max_proc=32

paths_to_store=['/analysis/results'+s for s in paths_to_launch]
sp=[]
indx_buf=0
FNULL = open(os.devnull, 'w')

#copy the script in the folder where it must be executed
for target_path in paths_to_launch:
	launched=0
	for irun in range(runmin,runmax+1):
		# move to directory and uncompress files
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
		os.chdir(dir_to_go)        #changes directory
		#copy the script in the run folder
		cmd_list = ['tar', '-xzf', 'positions.tar.gz']
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(cmd_list)
		indx_buf+=1
		cmd_list = ['tar', '-xzf', 'wallq.tar.gz']
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(cmd_list)
		indx_buf+=1
		cmd_list_1 = ['cp', '/data_prace/michele/prace/cdc/analysis/scripts/7_degree_of_confinement/degree_conf', dir_to_go]
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(cmd_list_1)
		indx_buf+=1

		create_input_files_for_degree_conf()
		for j in range(0, indx_buf):
        		sp[j].wait()
	
	print 'untar completed'       
	pid_list=[]
	sp=[]
	indx_buf=0
	launched=0
	for irun in range(runmin,runmax+1):
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
                os.chdir(dir_to_go) 
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(['./'+script_name_2[1:]+' < input_deg_conf_el_1.dat'],stdout=subprocess.PIPE,shell=True)
		#(out, err) = sp[indx_buf].communicate()
		pid_list.append(sp[indx_buf].pid)
		indx_buf+=1
		launched+=1
		print launched
                check=min([num_max_proc,runmax-runmin+1])
                if (launched>=check ):
                        for i in range(check,0,-1):
                                sp[indx_buf-i].wait()
                        launched=0
	for proc in sp:
		proc.wait()


	for irun in range(runmin,runmax+1):
	        #changes directory
		path_store = '/data_prace/michele/prace/cdc/analysis/results'+target_path+'/7_degree_of_confinement'
		cmd = 'mkdir -p '+path_store
		os.system(cmd)
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
                os.chdir(dir_to_go)

		cmd= 'mv Chargeconfmoy2.dat '+path_store+'/Chargeconfmoy2_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Chargeconfmoy3.dat '+path_store+'/Chargeconfmoy3_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coord2.dat '+path_store+'/Charge_coord2_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coord3.dat '+path_store+'/Charge_coord3_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coordmoy2.dat '+path_store+'/Charge_coordmoy2_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coordmoy3.dat '+path_store+'/Charge_coordmoy3_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Coordcarbonmoy2.dat '+path_store+'/Coordcarbonmoy2_el1_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Coordcarbonmoy3.dat '+path_store+'/Coordcarbonmoy3_el1_run%03d.out'%(irun)
		os.system(cmd)
		

	pid_list=[]
	sp=[]
	indx_buf=0
	launched=0
	for irun in range(runmin,runmax+1):
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
                os.chdir(dir_to_go) 
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(['./'+script_name_2[1:]+' < input_deg_conf_el_2.dat'],stdout=subprocess.PIPE,shell=True)
		#(out, err) = sp[indx_buf].communicate()
		pid_list.append(sp[indx_buf].pid)
		indx_buf+=1
		launched+=1
		print launched
                check=min([num_max_proc,runmax-runmin+1])
                if (launched>=check ):
                        for i in range(check,0,-1):
                                sp[indx_buf-i].wait()
                        launched=0
	for proc in sp:
		proc.wait()


	for irun in range(runmin,runmax+1):
	        #changes directory
		path_store = '/data_prace/michele/prace/cdc/analysis/results'+target_path+'/7_degree_of_confinement'
		cmd = 'mkdir -p '+path_store
		os.system(cmd)
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
                os.chdir(dir_to_go)

		cmd= 'mv Chargeconfmoy2.dat '+path_store+'/Chargeconfmoy2_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Chargeconfmoy3.dat '+path_store+'/Chargeconfmoy3_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coord2.dat '+path_store+'/Charge_coord2_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coord3.dat '+path_store+'/Charge_coord3_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coordmoy2.dat '+path_store+'/Charge_coordmoy2_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Charge_coordmoy3.dat '+path_store+'/Charge_coordmoy3_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Coordcarbonmoy2.dat '+path_store+'/Coordcarbonmoy2_el2_run%03d.out'%(irun)
		os.system(cmd)

		cmd= 'mv Coordcarbonmoy3.dat '+path_store+'/Coordcarbonmoy3_el2_run%03d.out'%(irun)
		os.system(cmd)


		cmd='rm -r '+script_name_2[1:]
		print cmd
		os.system(cmd)
		cmd ='rm *.out'
		print cmd
		os.system(cmd)

