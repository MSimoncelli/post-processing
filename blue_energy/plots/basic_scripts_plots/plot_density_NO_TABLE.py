import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os

from sys import argv
############################################################
if (len(argv)<=1):
        print 'USAGE: plor_coord_num.py runmin runmax'
        runmin = 0
        runmax = 53
else:
	runmin=int(argv[1])
	runmax=int(argv[2])

numbers={}

name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_dir=name_folder.split('/')
print list_dir
name_file_save=r'1_AVG_density_%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11] )
name_file_save= name_file_save.replace('.','_')
name_file_save+='.pdf'
name_file_load=r'density_AVG_%02d_%02d.dat' %(runmin,runmax)
if (list_dir[9]=='nacl'):
	list_dir[9]='NaCl'
elif (list_dir[9]=='kcl'):
	list_dir[9]='KCl'

title= r'Average density run%03d-run%03d %s %s %s %s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
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



if (list_dir[9]=='NaCl'):
	ax2.plot(numbers[0],numbers[5],color='r', label='Na' )
else:
	ax2.plot(numbers[0],numbers[5],color='r', label='K' )

ax2.plot(numbers[0],numbers[6],color='orange', label='Cl')

pl.legend()
pl.plot([17, 17], [0, 0.01], '--', lw=1, color='b')
pl.text(18, 0.0001, 'z=17')

pl.plot([40, 40], [0, 0.01], '--', lw=1, color='b')
pl.text(41, 0.0001, 'z=40')

pl.plot([55, 55], [0, 0.01], '--', lw=1, color='b')
pl.text(56, 0.0001, 'z=55')

pl.plot([146.406, 146.406], [0, 0.01], '--', lw=1, color='b')
pl.text(142.406, 0.0001, 'z=146.406')

pl.plot([169.406, 169.406], [0, 0.01], '--', lw=1, color='b')
pl.text(170.406, 0.0001, 'z=169.406')

pl.plot([131.406, 131.406], [0, 0.01], '--', lw=1, color='b')
pl.text(112.406, 0.0001, 'z=131.406')

max_4=max(numbers[5])
max_5=max(numbers[6])
max_ax_2=max([max_4,max_5])

ax2.set_ylim(0, max_ax_2*1.1)

#fig.tight_layout()
ax1.legend(loc=2)
ax2.legend(loc=1)
pl.gcf().subplots_adjust(right=0.87)
#pl.show()
fig.savefig(name_file_save, dpi=300)
