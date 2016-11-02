import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os
from collections import defaultdict
from sys import argv
############################################################
if (len(argv)<=1):
        print 'USAGE: plor_coord_num.py runmin runmax'
        runmin = 0
        runmax = 66
else:
	runmin=int(argv[1])
	runmax=int(argv[2])

numbers={}

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

def get_run_max():
        data_runmax=open('/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/plots/current_state_simulation_b.dat', 'r')
        lines_file = data_runmax.readlines()
        dict_runmax={}
        for line in lines_file:
                dict_runmax[line.split('\t')[0]]=int(line.split('\t')[1])-1
        data_runmax.close()
        return dict_runmax;


dict_runmax=get_run_max()

#list of directories in which we have to launch the script
list_sub_dir_here=os.walk('.').next()[1]
list_sub_dir_here=['/'+s for s in list_sub_dir_here]
common_parent_path='/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing'
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


max_time_dict=say_me_time_max()

name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
#print name_path


name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_dir=name_folder.split('/')
#print list_dir
name_file_save=r'1_AVG_density_%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11] )
name_file_save= name_file_save.replace('.','_')
name_file_save_text=name_file_save
name_file_save+='.pdf'
name_file_load=r'density_AVG_%02d_%02d.dat' %(runmin,runmax)
if (list_dir[9]=='nacl'):
	list_dir[9]='NaCl'
elif (list_dir[9]=='kcl'):
	list_dir[9]='KCl'


if ((list_dir[10]=='8800_160') or (list_dir[10]=='7615_139')):
	conC='1.0 M'
else:
	conC='0.5 M'

max_time_dict=say_me_time_max()

name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
#print name_path
list_path=name_path.split('/')
target_path=name_path.split('/results')[1]
print target_path
target_path=target_path.split('/1_dens')[0]

list_time_max=[0,0]
max_time=max_time_dict[target_path]
if runmin==0:
	list_time_max[0]=0
	if ((dict_runmax[target_path]-3)>runmax):
		list_time_max[1]=max_time/2.0
	else:
		list_time_max[1]=max_time
else:
	list_time_max[0]=max_time/2.0
	list_time_max[1]=max_time

title= r'Average density %1.2f-%1.2f ns %s %s %s %s' %(list_time_max[0],list_time_max[1], list_dir[8],list_dir[9],conC,list_dir[11])
title=title.replace('_', ' ')

from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################

fig=pl.figure()
ax1 = fig.add_subplot(111)
data = np.loadtxt(name_file_load)
for i in xrange(0,10):
	numbers[i] = data[:,i]

ax1.set_title(title)
ax1.grid(True)
# ax1.set_xlim([-50,50])
ax1.set_xlabel(r'z coordinate $[\AA]$')
ax1.set_ylabel(r'density $\left[\frac{1}{\AA^3}\right]$')

ax1.plot(numbers[0],numbers[2],color='b', label='O ' )	
ax1.plot(numbers[0],numbers[7],color='c', label='C1')
ax1.plot(numbers[0],numbers[8],color='purple', label='C2' )

max_1=max(numbers[2])
max_2=max(numbers[7])
max_3=max(numbers[8])
max_ax_1=max([max_1,max_2,max_3])
ax1.set_ylim(0, max_ax_1*1.1)

ax2 = ax1.twinx()
ax2.set_xlim(0,200)

ax2.set_ylabel(r'density $\left[\frac{1}{\AA^3}\right]$')



for tl in ax2.get_yticklabels():
    tl.set_color('r')

for tl in ax1.get_yticklabels():
    tl.set_color('b')

for i in range(0,len(numbers[0])):
	if numbers[0][i]<17:
		index_1_1=i
	if numbers[0][i]<40:
		index_1_2=i
	if numbers[0][i]<55:
		index_2_1=i
	if numbers[0][i]<131.406:
		index_2_2=i
	if numbers[0][i]<146.406:
		index_3_1=i
	if numbers[0][i]<169.406:
		index_3_2=i

avg_O_1=sum(numbers[2][index_1_1:index_1_2])/(len(numbers[2][index_1_1:index_1_2]))
inc_O_1=(max(numbers[2][index_1_1:index_1_2])-min(numbers[2][index_1_1:index_1_2]))*0.3

avg_O_2=sum(numbers[2][index_2_1:index_2_2])/(len(numbers[2][index_2_1:index_2_2]))
inc_O_2=(max(numbers[2][index_2_1:index_2_2])-min(numbers[2][index_2_1:index_2_2]))*0.3

avg_O_3=sum(numbers[2][index_3_1:index_3_2])/(len(numbers[2][index_3_1:index_3_2]))
inc_O_3=(max(numbers[2][index_3_1:index_3_2])-min(numbers[2][index_3_1:index_3_2]))*0.3


avg_Na_1=sum(numbers[5][index_1_1:index_1_2])/(len(numbers[5][index_1_1:index_1_2]))
inc_Na_1=(max(numbers[5][index_1_1:index_1_2])-min(numbers[5][index_1_1:index_1_2]))*0.3

avg_Na_2=sum(numbers[5][index_2_1:index_2_2])/(len(numbers[5][index_2_1:index_2_2]))
inc_Na_2=(max(numbers[5][index_2_1:index_2_2])-min(numbers[5][index_2_1:index_2_2]))*0.3

avg_Na_3=sum(numbers[5][index_3_1:index_3_2])/(len(numbers[5][index_3_1:index_3_2]))
inc_Na_3=(max(numbers[5][index_3_1:index_3_2])-min(numbers[5][index_3_1:index_3_2]))*0.3

avg_Cl_1=sum(numbers[6][index_1_1:index_1_2])/(len(numbers[6][index_1_1:index_1_2]))
inc_Cl_1=(max(numbers[6][index_1_1:index_1_2])-min(numbers[6][index_1_1:index_1_2]))*0.3

avg_Cl_2=sum(numbers[6][index_2_1:index_2_2])/(len(numbers[6][index_2_1:index_2_2]))
inc_Cl_2=(max(numbers[6][index_2_1:index_2_2])-min(numbers[6][index_2_1:index_2_2]))*0.3

avg_Cl_3=sum(numbers[6][index_3_1:index_3_2])/(len(numbers[6][index_3_1:index_3_2]))
inc_Cl_3=(max(numbers[6][index_3_1:index_3_2])-min(numbers[6][index_3_1:index_3_2]))*0.3

if (list_dir[9]=='NaCl'):
	ion_name='Na'
else:
	ion_name='K'

ax2.plot(numbers[0],numbers[5],color='r', label=ion_name )
ax2.plot(numbers[0],numbers[6],color='orange', label='Cl')

pl.legend()
pl.plot([17, 17], [0, 0.01], '--', lw=1, color='b')
pl.text(18, 0.0001, 'z=17')

pl.plot([40, 40], [0, 0.01], '--', lw=1, color='b')
pl.text(41, 0.0001, 'z=40')

pl.plot([55, 55], [0, 0.01], '--', lw=1, color='b')
pl.text(56, 0.0001, 'z=55')

pl.plot([146.406, 146.406], [0, 0.01], '--', lw=1, color='b')
pl.text(142.406, 0.0001, 'z=146')

pl.plot([169.406, 169.406], [0, 0.01], '--', lw=1, color='b')
pl.text(170.406, 0.0001, 'z=169')

pl.plot([131.406, 131.406], [0, 0.01], '--', lw=1, color='b')
pl.text(112.406, 0.0001, 'z=131')

max_4=max(numbers[5])
max_5=max(numbers[6])
max_ax_2=max([max_4,max_5])

ax2.set_ylim(0, max_ax_2*1.1)

rows = [' ', ion_name, 'Cl', 'O']
columns = (r'Average $\left[1/\AA^3\right]$', 'electrode 1', 'bulk', 'electrode 2')
colors = ['w','r','orange','#2E64FE']

#cell_text=[['%1.2f'%(avg_Cl_1),'%1.2f','%1.2f'],['%d','%d','%d'],['%d','%d','%d']
cell_text=[
[ion_name,
r'%6.*f'%(int(('%1.0e'%(inc_Na_1)).split('e-')[-1]),avg_Na_1)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_Na_1)).split('e-')[-1]),inc_Na_1),
r'%6.*f'%(int(('%1.0e'%(inc_Na_2)).split('e-')[-1]),avg_Na_2)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_Na_2)).split('e-')[-1]),inc_Na_2),
r'%6.*f'%(int(('%1.0e'%(inc_Na_3)).split('e-')[-1]),avg_Na_3)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_Na_3)).split('e-')[-1]),inc_Na_3)
],
['Cl',
r'%6.*f'%(int(('%1.0e'%(inc_Cl_1)).split('e-')[-1]),avg_Cl_1)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_Cl_1)).split('e-')[-1]),inc_Cl_1),
r'%6.*f'%(int(('%1.0e'%(inc_Cl_2)).split('e-')[-1]),avg_Cl_2)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_Cl_2)).split('e-')[-1]),inc_Cl_2),
r'%6.*f'%(int(('%1.0e'%(inc_Cl_3)).split('e-')[-1]),avg_Cl_3)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_Cl_3)).split('e-')[-1]),inc_Cl_3)
],
['O',
r'%6.*f'%(int(('%1.0e'%(inc_O_1)).split('e-')[-1]),avg_O_1)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_O_1)).split('e-')[-1]),inc_O_1),
r'%6.*f'%(int(('%1.0e'%(inc_O_2)).split('e-')[-1]),avg_O_2)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_O_2)).split('e-')[-1]),inc_O_2),
r'%6.*f'%(int(('%1.0e'%(inc_O_3)).split('e-')[-1]),avg_O_3)+r'$\pm$'+'%6.*f'%(int(('%1.0e'%(inc_O_3)).split('e-')[-1]),inc_O_3)
]
]

the_table = pl.table(cellText=cell_text,
                      #colColours=colors,
                      colLabels=columns,
                      cellLoc = 'center', rowLoc = 'center',
					  loc='bottom', bbox=[0.0, -0.45, 1, 0.3])

text='#'
for col_label in columns:
	text+=col_label+'\t'
text+='\n'

for list_ion in cell_text:
	for string in list_ion:
		text+=string.replace('$\pm$', '\t') +'\t'
	text+='\n'

tmpf=open(name_file_save_text+'.txt','w')
tmpf.write(text)
tmpf.close()
table_props = the_table.properties()
table_cells = table_props['child_artists']
for cell in table_cells: cell.set_height(0.1)
for cell in table_cells: cell.set_width(0.3)
# Adjust layout to make room for the table:
#pl.subplots_adjust(left=0.2, bottom=0.2)

#fig.tight_layout()
ax1.legend(loc=2)
ax2.legend(loc=1)
#pl.gcf().subplots_adjust(right=0.87,bottom=0.3)
fig.set_size_inches(6, 4.8)
fig.savefig(name_file_save, dpi=150, bbox_inches="tight")

