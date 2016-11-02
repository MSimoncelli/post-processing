import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os
from sys import argv

run=0

global n_ions
############################################################
if (len(argv)<=1):
        print 'USAGE: plor_coord_num.py n_bins'
        n_bins=200
else:
	n_bins=int(argv[1])

from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################
def build_histo(array,histo_O,histo_C,histo_Na_Cl,Nbins):
	for k in range(0,len(array[:,0])):
		histo_pos=math.ceil((array[k,0]/max_dist)*(Nbins-1)) #since floor 
		if (k<n_ions/2):
			histo_O[histo_pos,0]+=array[k,1]
			histo_Na_Cl[histo_pos]+=array[k,2] #NaCl or ClNa is the same (taking into account how you did calculations)
			histo_C[histo_pos,0]+=array[k,3]
		else:
			histo_O[histo_pos,1]+=array[k,1]
			histo_C[histo_pos,1]+=array[k,3]


data={}
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

############################################
data_path={}
num_confs={}
#load all the data in memory
name_path=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_path=name_path.split('/')

file_list=os.listdir(os.getcwd())
data_num=[]

runmin=0

data= np.loadtxt('run%03dRDF_W_W_for.dat'%(run))


lun= len(data[:,0])

#x vector is the same for all!!      
x=data[0:lun/3,0]
lung=len(x)

#print y['C-Cl'][:,0]
r_cut_off={}
r_cut_off['H2O-H2O']=3.40
y=np.zeros((lung,3))
y[:,0]+=data[0:lun/3,1]
y[:,1]+=data[lun/3:lun*2.0/3,1]
y[:,2]+=data[lun*2.0/3:lun,1]

fig=pl.figure()
ax1=fig.add_subplot(111)
name_file_save=r'3_RDF_%03d_%s_%s_%s_%s_%s' %( run, list_path[8],list_path[9],list_path[10],list_path[11],'H2O-H2O' )
name_file_save= name_file_save.replace('.','_')
name_file_save+='.pdf'	
title='RDF_%s_%s_%s_%s_%s' %( list_path[8],list_path[9],list_path[10],list_path[11],'H2O-H2O' )
title=title.replace('_', ' ')
ax1.set_title(title)
k='H2O-H2O'
#print title
ax1.grid(True)
ax1.set_xlabel(r'radial distance $[\AA]$')
ax1.set_ylabel(r'RDF')
ax1.set_xlim([0, max(x)])
if (not(math.isnan(y[1,0]))):
	string= 'RDF '+k+' el. 1'
	#print string
	ax1.plot(x,y[:,0], label=string,c='b')
if (not(math.isnan(y[1,1]))):
	string='RDF '+k+' bulk'
	ax1.plot(x,y[:,1], label=string,c='orange')
if (not(math.isnan(y[1,2]))):
	string='RDF '+k+' el. 2'
	ax1.plot(x,y[:,2], label=string,c='r')	
ax1.legend(loc='best')
#pl.plot([0, 25], [1, 1], '--', lw=1, color='black')
l = pl.axhline(y=1,ls='--', color='black')
max_1=max(y[:,0])
max_2=max(y[:,1])
max_3=max(y[:,2])
max_tot=max([max_1,max_2,max_3])
h_line=1+max_tot*0.025
#print h_line
pl.text(0.5, h_line, '1.0')
if ((k!='C-Cl')and(k!='C-Na')):
	l = pl.axvline(x=r_cut_off[k],ls='--', color='black')
	pl.text(r_cut_off[k]*1.05, max_tot*0.8, '%1.2f'%(r_cut_off[k]))
	#pl.plot([r_cut_off[k], r_cut_off[k]], [100, 100], '--', lw=1, color='black')
#pl.show()
fig.savefig(name_file_save, dpi=300)

