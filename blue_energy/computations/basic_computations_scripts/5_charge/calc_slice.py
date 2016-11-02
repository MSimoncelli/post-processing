import math
import numpy as np

def read_info_specie():
	file_rntime=open ('runtime.inpt','r')
	file_positions=file_rntime.readlines() # read all the file, it contains lots of useful things
	n_part=[] #those are global variables
	spc=[]
	array_line_num_atom = [7,12,17,22,27,32,37,42]
	for line_num_atom in array_line_num_atom:
		line_specie = file_positions[line_num_atom].split()
		Nspec=int(line_specie[0])
		line_name_specie=file_positions[line_num_atom-1].split()
		name_spec=line_name_specie[0]
		n_part.append(Nspec)
		spc.append(name_spec)
	return n_part, spc

n_part, spc= read_info_specie()
atoms_per_conf=sum(n_part)
base_shift_ion=sum(n_part[0:5])

to_app = np.empty((0,3), float)
list_electrode_1=[]
list_electrode_2=[]
with open('positions.out') as file_pos:
    data_line = [next(file_pos) for x in range(0, atoms_per_conf)]
    
for line_idx in range(base_shift_ion,base_shift_ion+n_part[5]):
	to_app=np.fromstring(data_line[line_idx], dtype=float, sep=' ')
	list_electrode_1.append(to_app)

electrode_1 = np.asarray(list_electrode_1)
	
for line_idx in range(base_shift_ion+n_part[5],base_shift_ion+n_part[5]+n_part[6]):
	to_app=np.fromstring(data_line[line_idx], dtype=float, sep=' ')
	list_electrode_2.append(to_app)

electrode_2 = np.asarray(list_electrode_2)

file_wallq = np.loadtxt('wallq.out')
num_line_wall_q=len(file_wallq)

atoms_C_per_conf=sum(n_part[5:8])
num_configs=num_line_wall_q/atoms_C_per_conf

length_elecrtrode=max(electrode_1[:,2])
start_electrode_2=min(electrode_2[:,2])

base_wall_q=0
text='#Q_el_1[0]\tQ_el_1[1]\tQ_el_2[0]\tQ_el_2[1]\n' 
base_shift_pos=sum(n_part[0:5])
for n_frame in range(0,num_configs):
	Q_electrode_1=np.zeros(2)
	Q_electrode_2=np.zeros(2)
	for idx_C in range(0,n_part[5]): #first el
		Q_electrode_1[int(round(electrode_1[idx_C, 2]/length_elecrtrode))]+=file_wallq[base_wall_q+idx_C]
	for idx_C in range(n_part[5],n_part[5]+n_part[6]): #second electrode
		Q_electrode_2[int(round((electrode_2[idx_C-n_part[5], 2]-start_electrode_2)/length_elecrtrode))]+=file_wallq[base_wall_q+idx_C]

	text+='%e\t%e\t%e\t%e\n' %(Q_electrode_1[0],Q_electrode_1[1],Q_electrode_2[0],Q_electrode_2[1])
	base_wall_q+=n_part[5]+n_part[6]+n_part[7]
	base_shift_pos+=atoms_per_conf

	tmpf=open('charge_slice.dat', 'w')
	tmpf.write(text)
	tmpf.close()
