import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os

from sys import argv

conv=0.529177
############################################################
if (len(argv)<=1):
        print 'USAGE: plor_coord_num.py runmin runmax'
        runmin = 0
        runmax = 67
else:
	runmin=int(argv[1])
	runmax=int(argv[2])

numbers={}

name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_dir=name_folder.split('/')
#print list_dir
name_file_save=r'd_1_AVG_density_%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11] )
name_file_save= name_file_save.replace('.','_')
name_file_save_text=name_file_save
name_file_save+='.pdf'
name_file_load=r'c_density_AVG_%02d_%02d.dat' %(runmin,runmax)
if (list_dir[9]=='nacl'):
	list_dir[9]='NaCl'
elif (list_dir[9]=='kcl'):
	list_dir[9]='KCl'

title= r'Average density from 2.5 to 3.0 ns run%03d-%03d %s %s %s %s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
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

index_1_1=0

for i in range(0,len(numbers[0])):
	#if numbers[0][i]<0.5:
	#	index_1_1=i
	if numbers[0][i]<44:
		index_1_2=i
	if numbers[0][i]<44:
		index_2_1=i
	if numbers[0][i]<142.706:
		index_2_2=i
	if numbers[0][i]<142.706:
		index_3_1=i
	if numbers[0][i]<186.406:
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
if 'cdc1200' in list_dir[8]:
	l_1=87.751484*conv
	l_2=(352.2564283784-87.751484)*conv
	l_3=352.2564283784*conv      
else:
	l_1=46.436067
	l_2=183.000000-46.436067
	l_3=183.000000   




pl.plot([l_1, l_1], [0, 0.01], '--', lw=1, color='b')
pl.text(l_1+1, 0.0001, r'z=%1.2f'%(l_1))


pl.plot([l_2, l_2], [0, 0.01], '--', lw=1, color='b')
pl.text(l_2+1, 0.0001, r'z=%1.2f'%(l_2))


pl.plot([l_3, l_3], [0, 0.01], '--', lw=1, color='b')
pl.text(l_3+1, 0.0001, r'z=%1.2f'%(l_3))

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
r'%1.8f'%(avg_Na_1),
r'%1.8f'%(avg_Na_2),
r'%1.8f'%(avg_Na_3)
],
['Cl',
r'%1.8f'%(avg_Cl_1),
r'%1.8f'%(avg_Cl_2),
r'%1.8f'%(avg_Cl_3),
],
['O',
r'%1.8f'%(avg_O_1),
r'%1.8f'%(avg_O_2),
r'%1.8f'%(avg_O_3)
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
pl.gcf().subplots_adjust(right=0.87,bottom=0.3)
#pl.show()
fig.savefig(name_file_save, dpi=300)
