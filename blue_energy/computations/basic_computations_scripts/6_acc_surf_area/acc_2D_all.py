import math
import numpy as np
np.set_printoptions(threshold='nan')
import time
import os
global num_lines_file
global Nbins
global n_dist
global histo_el_1
global histo_bulk
global histo_el_2
global n_spec_ions #number of different ions


D_obs_D_probe=9.45

bohr= 0.529177 		# in angstrom
Nbins=100
n_dist=4
n_spec_ions=2
num_regions=3 #electrode 1, bulk, electrode 2

histo=np.zeros([num_regions,n_dist,n_spec_ions,Nbins], dtype=float)
num_density=np.zeros([n_dist,n_spec_ions,num_regions], dtype=float)

#for i in xrange(0,2):
#	histo_el_1[i+1,i,0]-=1 #remove the target ion from its radial distribution function
#	histo_el_2[i+1,i,0]-=1
#	histo_bulk[i+1,i,0]-=1

def read_info_specie(line_num_atom, line_name, file_positions):
	line_specie = file_positions[line_num_atom].split()
	Nspec=int(line_specie[0])
	line_name_specie=file_positions[line_name].split()
	name_spec=line_name_specie[0]
	return Nspec, name_spec


def dist_w_PBC_opt(coor_a,coor_b,cell_size):
	halfboxxrec=2.0/cell_size[0]
	halfboxyrec=2.0/cell_size[1]

	dxcf=coor_a[0]-coor_b[0]
	dycf=coor_a[1]-coor_b[1]
	dzcf=coor_a[2]-coor_b[2]

	#minimal distance convenction
	dxcf=dxcf-cell_size[0]*int(dxcf*halfboxxrec)
	dycf=dycf-cell_size[1]*int(dycf*halfboxyrec)

	res=math.sqrt(dxcf**2+dycf**2+dzcf**2)
	return res

def read_sim_box_infos(dirrun):
	start_dir=os.getcwd()
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

dir_cdc1200='/data_prace/michele/prace/cdc/cdc1200/nacl/8800_160/0.0V/run000'
dir_cdc800='/data_prace/michele/prace/cdc/cdc800/nacl/7615_139/0.0V/run000'

lis_dir_to_go=[dir_cdc1200, dir_cdc800]

for dir_file_struct in lis_dir_to_go:
	start_dir=os.getcwd()
	os.chdir(dir_file_struct)
	cmd='tar -xzf positions.tar.gz'
	os.system(cmd)
	box_size=read_sim_box_infos(dir_file_struct)

	conv=0.529177	

	file_rntime=open ('runtime.inpt','r')
	file_runtime_a=file_rntime.readlines() # read all the file, it contains lots of useful things	

	n_part=[]
	spc=[]
	array_index = [7,12,17,22,27,32,37,42]	

	for index in array_index:
		Nspec,name_spec=read_info_specie(index,index-1,file_runtime_a)
		n_part.append(Nspec)
		spc.append(name_spec)
		print 'num %s=%d' %(name_spec, Nspec)	

	atoms_per_conf=0
	n_part_tot=0
	n_ions=0
	nspec=8	
	

	for k in xrange(0,len(spc)):
		if (((spc[k]=='Na')or(spc[k]=='Cl'))or(spc[k]=='K')):
			n_ions=n_ions+n_part[k]	

	#read all the file->slow operation	

	file_posit_memo = np.loadtxt('positions.out')
	num_line_pos=len(file_posit_memo)	

	atoms_per_conf=sum(n_part)
	num_configs=num_line_pos/atoms_per_conf	

	spc_upd=[]
	# names
	spc_upd.append('O')
	spc_upd.append(spc[3]) #Na or K
	spc_upd.append(spc[4]) #Cl 
	spc_upd.append('C')	

	n_specie_dist=[]
	#index: O, H1, H2, Na, Cl, C1, C2, P
	base_shift_ion=sum(n_part[:5])	

	electrode_1=file_posit_memo[base_shift_ion:base_shift_ion+n_part[5], :].copy()
	electrode_2=file_posit_memo[base_shift_ion+n_part[5]:base_shift_ion+n_part[5]+n_part[6], :].copy()
	cwd=os.getcwd()
	cwd=cwd.replace('/','_')
		
	name_folder='electrode_'+cwd[30:36]
	cmd= 'mkdir -p '+name_folder
	os.system(cmd)	

	nome_file_1='el_1_'+cwd[30:36]+'.dat'
	nome_file_2='el_2_'+cwd[30:36]+'.dat'	

	np.savetxt(name_folder+'/'+nome_file_1, electrode_1)
	np.savetxt(name_folder+'/'+nome_file_2, electrode_2)	

	z_max_electrode_1= max(electrode_1[:,2])# max z coordinate of the atoms of the electrode. Should be around 90.	

	text_1='%s\n'%(nome_file_1) #file comntaining the coordinates of the atoms of the electrode
	text_1+='%d\n' %(n_part[5]) #number of atoms in the electrode
	text_1+='%f\n'%(D_obs_D_probe) #third line : diameter of the obstacles in bohrs (> Dobs)
	text_1+= '%f %f\n'%(box_size[0], box_size[1]) #82.54 82.54 #length of the box in the three directions in bohrs (> Lx,Ly)
	text_1+='%f %f\n'%(0, z_max_electrode_1+0.5)#zmin and zmax (> zmin,zmax)0.0 90.0 
	text_1+='160\n' #number of bins in any direction (> nbins)
	text_1+='0.0' # diameter of the probe in bohrs (Argon atom)	

	tmpf=open(name_folder+'/input_vol_1.dat', 'w')
	tmpf.write(text_1)
	tmpf.close()	
	

	text_2='%s\n'%(nome_file_2) #file comntaining the coordinates of the atoms of the electrode
	text_2+='%d\n' %(n_part[5]) #number of atoms in the electrode
	text_2+='%f\n'%(D_obs_D_probe) #third line : diameter of the obstacles in bohrs (> Dobs)
	text_2+= '%f %f\n'%(box_size[0], box_size[1]) #82.54 82.54 #length of the box in the three directions in bohrs (> Lx,Ly)
	text_2+='%f %f\n'%(box_size[2]-(z_max_electrode_1+0.5), box_size[2])#zmin and zmax (> zmin,zmax)0.0 90.0 
	text_2+='160\n' #number of bins in any direction (> nbins)
	text_2+='0.0' # diameter of the probe in bohrs (Argon atom	

	tmpf_2=open(name_folder+'/input_vol_2.dat', 'w')
	tmpf_2.write(text_2)
	tmpf_2.close()	

	cmd='mv ./'+name_folder+' '+start_dir
	os.system(cmd)
	os.chdir(start_dir)


cmd='cp acc_vol_2D ./electrode_cdc800'
os.system(cmd)
cmd='cp acc_vol_2D ./electrode_cdc120'
os.system(cmd)

os.chdir('./electrode_cdc800')
cmd = './acc_vol_2D<input_vol_1.dat'
os.system(cmd)

cmd='mv access_surf.xyz access_surf_cdc800_el1.xyz'
os.system(cmd)
cmd='mv access_vol.xyz access_vol_cdc800_el1.xyz'
os.system(cmd)
cmd='mv surf_points.out surf_points_cdc800_el1.out'
os.system(cmd)
cmd='mv vol_points.out vol_points_cdc800_el1.out'
os.system(cmd)



cmd='./acc_vol_2D<input_vol_2.dat'
os.system(cmd)

cmd='mv access_surf.xyz access_surf_cdc800_el2.xyz'
os.system(cmd)
cmd='mv access_vol.xyz access_vol_cdc800_el2.xyz'
os.system(cmd)
cmd='mv surf_points.out surf_points_cdc800_el2.out'
os.system(cmd)
cmd='mv vol_points.out vol_points_cdc800_el2.out'
os.system(cmd)

os.chdir('../electrode_cdc120')
cmd = './acc_vol_2D<input_vol_1.dat'
os.system(cmd)

cmd='mv access_surf.xyz access_surf_cdc1200_el1.xyz'
os.system(cmd)
cmd='mv access_vol.xyz access_vol_cdc1200_el1.xyz'
os.system(cmd)
cmd='mv surf_points.out surf_points_cdc1200_el1.out'
os.system(cmd)
cmd='mv vol_points.out vol_points_cdc1200_el1.out'
os.system(cmd)

cmd='./acc_vol_2D<input_vol_2.dat'
os.system(cmd)

cmd='mv access_surf.xyz access_surf_cdc1200_el2.xyz'
os.system(cmd)
cmd='mv access_vol.xyz access_vol_cdc1200_el2.xyz'
os.system(cmd)
cmd='mv surf_points.out surf_points_cdc1200_el2.out'
os.system(cmd)
cmd='mv vol_points.out vol_points_cdc1200_el2.out'
os.system(cmd)

