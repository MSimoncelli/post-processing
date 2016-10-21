#!/usr/bin/python
import string, re, struct, sys, math, os, inspect
import subprocess
from sys import argv
import numpy as np
from os import listdir
from os.path import isfile, join

startdir = os.getcwd() #gets the current working directory
flag_untar=False

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

def read_species():
	file_rntime=open ('runtime.inpt','r')
	file_runtime_a=file_rntime.readlines() # read all the file, it contains lots of useful things   
	n_part=[]
	spc=[]
	array_index = [7,12,17,22,27,32,37,42]
	for index in array_index:
		Nspec,name_spec=read_info_specie(index,index-1,file_runtime_a)
		n_part.append(Nspec)
		spc.append(name_spec)
	return n_part, spc

def get_run_max():
	data_runmax=open('current_state_simulation_b.dat', 'r')
	lines_file = data_runmax.readlines()
	dict_runmax={}
	for line in lines_file:
	        dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])-1
	data_runmax.close()
	return dict_runmax;

dict_runmax=get_run_max()
#all the files will start from this path
common_parent_path='/data_prace/michele/prace/cdc'

cmd='mkdir -p ./scripts/5_conductivity_kubo/inputs_conductivity_1'
os.system(cmd)
cmd='mkdir -p ./scripts/5_conductivity_kubo/inputs_conductivity_2'
os.system(cmd)
cmd='mkdir -p ./scripts/5_conductivity_kubo/inputs_conductivity_3'
os.system(cmd)
cmd='mkdir -p ./scripts/5_conductivity_kubo/outputs_conductivity_1'
os.system(cmd)
cmd='mkdir -p ./scripts/5_conductivity_kubo/outputs_conductivity_2'
os.system(cmd)
cmd='mkdir -p ./scripts/5_conductivity_kubo/outputs_conductivity_3'
os.system(cmd)



#/data_prace/michele/prace/cdc/analysis/scripts/6_acc_surf_area
script_folder_2='/6_acc_surf_area' #pay attention at the arguments
script_name_2='/idens_local_surf' #extracts only the name of the script
#put in this list all the paths in which you want to launch script[i]
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

runmin=0

num_max_proc=32

box_size=read_sim_box_infos()

sp=[]
indx_buf=0
FNULL = open(os.devnull, 'w')

tot_num_frames={}
for target_path in paths_to_launch:
	string_nam=target_path
	string_nam=string_nam.replace('/','_')
	string_nam=string_nam.replace('.','_')
	tot_num_frames[target_path]=0 
	dir_to_go=common_parent_path+target_path+'/run000'
	os.chdir(dir_to_go) 
	n_part_diff, spec=read_species()
	#copy the script in the folder where it must be executed
	text='%d\n'%(len(n_part_diff))
	text+='%d,.false.\n'%(n_part_diff[0]) #O
	text+='%d,.false.\n'%(n_part_diff[1]) #H1
	text+='%d,.false.\n'%(n_part_diff[2]) #H2
	text+='%d,.true.\n'%(n_part_diff[3])  #Na or K
	text+='%d,.true.\n'%(n_part_diff[4]) #Cl
	text+='%d,.false.\n'%(n_part_diff[5]) #C1
	text+='%d,.false.\n'%(n_part_diff[6]) #C2
	text+='%d,.false.\n'%(n_part_diff[7]) #P

	runmax=dict_runmax[target_path]
	text+='%d\n'%(runmax+1) #bacause it starts from run000!!!
	for irun in range(runmin,runmax+1):
		dir_to_go=common_parent_path+target_path+'/run%03d'%(irun)
		os.chdir(dir_to_go)        #changes directory
		num_frames=read_numframes()
		tot_num_frames[target_path]+=num_frames
		text+=dir_to_go+'/positions.out'+'\n%d\n'%(num_frames)
		if (flag_untar):
			cmd_list = ['tar', '-xzf', 'positions.tar.gz']
                        sp.append(1)
                        sp[indx_buf] = subprocess.Popen(cmd_list)
                        indx_buf+=1
	if (flag_untar):
		for j in range(0, indx_buf):
			sp[j].wait()

	text+='%f\n'%(box_size[0])
	text+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_1/dispfrompositions'+string_nam+'.out\n' #outfile
	text+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_1/positions'+string_nam+'.out\n' #outposfile
	text+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_1/ALL_displacement'+string_nam+'.out\n' #outposfile
	os.chdir(startdir)

	tmpf_1=open('./scripts/5_conductivity_kubo/inputs_conductivity_1/displ_f_pos'+string_nam+'.inpt', 'w')
	tmpf_1.write(text)
	tmpf_1.close()
	#now make the input file for the layer program
	for working_species in [1,2]: #specie 1=Na or K (the 1st one for which you have extracted the displacement), soecie 2 is Cl. YOU HAVE EXTRACTED THE DISPLACEMENTS OF THOSE 2 ONLY!!
                text2='%d\n'%(tot_num_frames[target_path]-1)#-1 because when you look at the displacemet you loose a frame (2 frames to compute 1 displacement!!!!)
                text2+='2\n' #2 species Na or K and Cl
                text2+='%d\n' %(working_species) #which species you are focusing on 1: Naor K, 2: Cl
                text2+='%d\n'%(n_part_diff[3])  #Na or K
                text2+='%d\n'%(n_part_diff[4]) #Cl
                text2+='3\n' #number of ziones considered!
		if 'cdc800' in target_path:
			z_max_electrode=87.75148391	
		if 'cdc1200' in target_path:
			z_max_electrode=87.94962366
                text2+='%f,%f \n'%(0.0, z_max_electrode) #boundary of zone 1
                text2+='%f,%f \n'%(z_max_electrode,max(box_size)-z_max_electrode) #boundary zone 2
                text2+='%f,%f \n'%(max(box_size)-z_max_electrode,max(box_size) ) #boundary zone 3
		text2+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_1/positions'+string_nam+'.out\n'
		text2+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_1/dispfrompositions'+string_nam+'.out\n' #outfile
		text2+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_2/sortedisp-C'+string_nam+'_spc%d.out\n'%(working_species)
		text2+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_2/keydisp-C'+string_nam+'_spc%d.out\n'%(working_species)
		tmpf_2=open('./scripts/5_conductivity_kubo/inputs_conductivity_2/layer'+string_nam+'_spc_%d.inpt'%(working_species), 'w')
		tmpf_2.write(text2)
		tmpf_2.close()


path_inputs=startdir+'/scripts/5_conductivity_kubo/inputs_conductivity_1'

file_list = [f for f in listdir(path_inputs) if isfile(join(path_inputs, f))]
#print file_list

sp=[]
indx_buf=0

for file_inpt in file_list:
	cmd_list = [startdir+'/scripts/5_conductivity_kubo/disp_from_positions.x',path_inputs+'/'+file_inpt]
	print cmd_list
	sp.append(1)
	sp[indx_buf] = subprocess.Popen(cmd_list)
	indx_buf+=1

#wait all the launched process to complete
for j in range(0, indx_buf):
	sp[j].wait()


path_inputs_2=startdir+'/scripts/5_conductivity_kubo/inputs_conductivity_2'

file_list_2 = [f for f in listdir(path_inputs_2) if isfile(join(path_inputs_2, f))]
#print file_list

sp=[]
indx_buf=0

for file_inpt_2 in file_list_2:
        cmd_list = [startdir+'/scripts/5_conductivity_kubo/layers.x',path_inputs_2+'/'+file_inpt_2]
        sp.append(1)
	print cmd_list
        sp[indx_buf] = subprocess.Popen(cmd_list)
        indx_buf+=1

for j in range(0, indx_buf):
        sp[j].wait()

dir_here=os.getcwd()
path_outputs_2=startdir+'/scripts/5_conductivity_kubo/outputs_conductivity_2'
file_list_3 = [f for f in listdir(path_outputs_2) if isfile(join(path_outputs_2, f))]

os.chdir(path_outputs_2)
num_line_useful={}
for file_out_3 in file_list_3:
        if ("keydisp-C" in file_out_3):
		file_posit_memo = np.loadtxt(file_out_3)
		num_line_useful[file_out_3]=len(file_posit_memo)


for file_out_3 in num_line_useful.keys():
	string_nam=file_out_3.split('keydisp-C_')[1]
	string_nam=string_nam.split('.out')[0]
	for layer_idx in [1,2,3]:
		text3='%d\n'%(num_line_useful[file_out_3])
		text3+='%d\n'%(layer_idx) #layer on which we do the analysis
		text3+='%d\n'%(2000) #number of steps on which the correlation is made
		text3+='1000,%f\n' %(41.341) #timestep in atomic units!
		text3+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_2/sortedisp-C_'+string_nam+'.out\n'
		text3+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_2/keydisp-C_'+string_nam+'.out\n'
		text3+=common_parent_path+'/analysis/scripts/5_conductivity_kubo'+'/outputs_conductivity_3/msd-C_'+string_nam+'_layer_%d.out\n'%(layer_idx)
		tmpf_3=open(common_parent_path+'/analysis/scripts/5_conductivity_kubo/inputs_conductivity_3/zmsd-mat_'+string_nam+'_layer_%d.inpt'%(layer_idx), 'w')
		tmpf_3.write(text3)
		tmpf_3.close()
	
os.chdir(dir_here)

path_inputs_3=startdir+'/scripts/5_conductivity_kubo/inputs_conductivity_3'

file_list_3 = [f for f in listdir(path_inputs_3) if isfile(join(path_inputs_3, f))]
#print file_list

sp=[]
indx_buf=0

for file_inpt_3 in file_list_3:
	cmd_list = [startdir+'/scripts/5_conductivity_kubo/msdlayers-mat.x',path_inputs_3+'/'+file_inpt_3]
        sp.append(1)
	print cmd_list
        sp[indx_buf] = subprocess.Popen(cmd_list)
        indx_buf+=1

for j in range(0, indx_buf):
        sp[j].wait()
