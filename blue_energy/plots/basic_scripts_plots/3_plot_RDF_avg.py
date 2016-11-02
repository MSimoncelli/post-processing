import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os

from sys import argv
from collections import defaultdict

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
for name_dat_file in file_list:
        if ('Cl_Clfin.dat'in name_dat_file):
                data_num.append(int(name_dat_file[3:6]))
runmin=0

runmax=max(data_num)


for irun in range(runmin,runmax+1):
        data_path['C-Cl']='run%03dRDF_C_Cl_fin.dat'%(irun)
        data_path['C-Na']='run%03dRDF_C_Na_fin.dat'%(irun)
        data_path['Cl-Na']='run%03dRDF_Na_Clfin.dat'%(irun)
        data_path['O-Cl']='run%03dRDF_O_Cl_fin.dat'%(irun)
        data_path['O-Na']='run%03dRDF_O_Na_fin.dat'%(irun)
        data_path['Cl-Cl']='run%03dRDF_Cl_Clfin.dat'%(irun)
        data_path['Na-Na']='run%03dRDF_Na_Nafin.dat'%(irun)
        
        ##################################################################
        #FIND AND REPLACE IN A FILE: THIS IS NO LONGER NECESSARY BUT THIS IS USEFUL
        #s = open(data_path).read()
        #s = s.replace('Run ', '#Run_')
        #s = s.replace('label', '#labl')
        #f = open(data_path, 'w')
        #f.write(s)
        #f.close()
        ##################################################################
        #load all the data
        with open(data_path['C-Cl'], 'r') as f:
                first_line = f.readline()
    
        num_confs[irun]=int(first_line.split('=')[1]) 

        for k in data_path:
                data[irun,k]=np.loadtxt(data_path[k])
                #x[k], y[k]=data[irun][k][:,0], data[irun][k][:,1]


lun= len(data[runmin,'C-Cl'][:,0])

#x vector is the same for all!!
x=data[runmin,'C-Cl'][0:lun/3,0]
lung=len(x)

#print y['C-Cl'][:,0]
r_cut_off={}

tot_weights_1={}
tot_weights_2={}
tot_weights_3={}
set_list_bound=[[runmin,runmax],[runmin,runmax/2],[runmax/2+1,runmax]]

text='#el 1,\t bulk, \t el 2\n'
for list_bound in set_list_bound: 
        cdi={}
        text+= 'runmin=%03d \trunmax=%03d\n' %(list_bound[0],list_bound[1])
        for k in data_path:
                y[k]=np.zeros((lung,3))
                tot_weights_1[k]=0
                tot_weights_2[k]=0
                tot_weights_3[k]=0

                list_num=num_confs.values()
                list_num_crop=list_num[list_bound[0]:list_bound[1]+1]

                for irun in range(list_bound[0],list_bound[1]+1):
                        if (not(math.isnan((data[irun,k][1][1])))):
                                y[k][:,0]+=(data[irun,k][0:lun/3,1]*num_confs[irun])
                                tot_weights_1[k]+=list_num[irun]

                        if (not(math.isnan((data[irun,k][lun/3+1][1])))):
                                y[k][:,1]+=(data[irun,k][lun/3:lun*2.0/3,1]*num_confs[irun])
                                tot_weights_2[k]+=list_num[irun]

                        if (not(math.isnan((data[irun,k][lun*2.0/3+1][1])))):
                                y[k][:,2]+=(data[irun,k][lun*2.0/3:lun,1]*num_confs[irun])
                                tot_weights_3[k]+=list_num[irun]

                
                #print list_num_crop
                tot_weights=sum(list_num_crop)
                if (tot_weights_1[k]):
                        y[k][:,0]*=(1.0/(tot_weights_1[k])) #average!!
                if (tot_weights_2[k]):
                        y[k][:,1]*=(1.0/(tot_weights_2[k]))
                if (tot_weights_3[k]):
                        y[k][:,2]*=(1.0/(tot_weights_3[k]))

                fig=pl.figure()
                ax1=fig.add_subplot(111)        
                w=k
                if (list_path[9]=='nacl'):
                        list_path[9]='NaCl'
                        r_cut_off['Cl-Na']=3.80
                        r_cut_off['O-Cl']=3.9
                        r_cut_off['O-Na']=3.3
                        r_cut_off['Cl-Cl']=6.40
                        r_cut_off['Na-Na']=5.10
                        r_cut_off['C-Cl']=4.20
                        r_cut_off['C-Na']=2.60
                elif (list_path[9]=='kcl'):
                        list_path[9]='KCl'
                        w=w.replace("Na", "K")
                        r_cut_off['Cl-Na']=4.20
                        r_cut_off['O-Cl']=3.90
                        r_cut_off['O-Na']=3.70
                        r_cut_off['Cl-Cl']=6.40
                        r_cut_off['Na-Na']=5.80
                        r_cut_off['C-Cl']=4.20
                        r_cut_off['C-Na']=3.10
                elif(list_path[9]=='KCl'):
                        w=w.replace("Na", "K")

                name_file_save=r'3_RDF_%03d_%03d_%s_%s_%s_%s_%s' %(list_bound[0],list_bound[1], list_path[8],list_path[9],list_path[10],list_path[11],w )
                name_file_save= name_file_save.replace('.','_')
                name_file_save+='.pdf'

                title='RDF_%03d_%03d_%s_%s_%s_%s_%s' %(list_bound[0],list_bound[1], list_path[8],list_path[9],list_path[10],list_path[11],w )
                title=title.replace('_', ' ')
                ax1.set_title(title)
                #print title
                ax1.grid(True)
                ax1.set_xlabel(r'radial distance $[\AA]$')
                ax1.set_ylabel(r'RDF')
                ax1.set_xlim([0, max(x)])
                

                cdi[k]=np.zeros(3)
                if (sum(y[k][:,0])>1e-10):
                        string= 'RDF '+w+' el. 1'
                        #print string
                        ax1.plot(x,y[k][:,0], label=string,c='b')
                        a=(np.nonzero(y[k][:,0]))
                        min_1=np.amin(a, axis=1)
                        cdi[k][0]=x[min_1[0]]
                        
                if (sum(y[k][:,1])>1e-10):
                        string='RDF '+w+' bulk'
                        ax1.plot(x,y[k][:,1], label=string,c='orange')
                        a=(np.nonzero(y[k][:,1]))
                        min_1=np.amin(a, axis=1)
                        cdi[k][1]=x[min_1[0]]

                if (sum(y[k][:,2])>1e-10):
                        string='RDF '+w+' el. 2'
                        ax1.plot(x,y[k][:,2], label=string,c='r')  
                        a=(np.nonzero(y[k][:,2]))
                        min_1=np.amin(a, axis=1)
                        cdi[k][2]=x[min_1[0]]  

                text+= '%s\t%e\t%e\t%e\n' %(k,cdi[k][0],cdi[k][1],cdi[k][2])

                ax1.legend(loc='best')
                #pl.plot([0, 25], [1, 1], '--', lw=1, color='black')
                l = pl.axhline(y=1,ls='--', color='black')
                max_1=max(y[k][:,0])
                max_2=max(y[k][:,1])
                max_3=max(y[k][:,2])
                max_tot=max([max_1,max_2,max_3])
                h_line=1+max_tot*0.025
                #print h_line
                pl.text(0.5, h_line, '1.0') 

                # major ticks every 20, minor ticks every 5                                      
                major_ticks = np.arange(0, 11, 1)                                              
                minor_ticks = np.arange(0, 10.1, 0.1)                                                              

                ax1.set_xticks(major_ticks)                                                       
                ax1.set_xticks(minor_ticks, minor=True)                                           
                
                if (cdi[k][1]>0.5):                                                                                               
                        cdi_value=sum(cdi[k][:])/3.0
                else:
                        cdi_value=(cdi[k][0]+cdi[k][2])/2.0
                # or if you want differnet settings for the grids:                               
                #ax1.grid(which='minor', alpha=0.0)                                                
                ax1.grid(which='major', alpha=0.5)  
                y_ax_lim= ax1.get_ylim()
                if ((k!='C-Cl')and(k!='C-Na')and(k!='Cl-Cl')and(k!='Na-Na')):
                        l = pl.axvline(x=r_cut_off[k],ls='--', color='black')
                        pl.text(r_cut_off[k]+0.15, y_ax_lim[1]*0.7, r'%1.2f$\pm$0.06'%(r_cut_off[k]))
                        l = pl.axvline(x=cdi_value,ls='--', color='black')
                        pl.text(cdi_value-1.3, y_ax_lim[1]*0.7, r'%1.2f$\pm$0.06'%(cdi_value))
                elif((k!='C-Cl')and(k!='C-Na')and((k=='Cl-Cl'))or(k=='Na-Na') ):
                        l = pl.axvline(x=r_cut_off[k],ls='--', color='black')
                        pl.text(r_cut_off[k]+0.15, y_ax_lim[1]*0.15, r'%1.1f$\pm$0.3'%(r_cut_off[k]))
                        l = pl.axvline(x=cdi_value,ls='--', color='black')
                        pl.text(cdi_value-1.1, y_ax_lim[1]*0.15, r'%1.1f$\pm$0.1'%(cdi_value))
                elif((k=='C-Cl')):
                        l = pl.axvline(x=r_cut_off[k],ls='--', color='black')
                        pl.text(r_cut_off[k]+0.15, y_ax_lim[1]*0.15, r'%1.1f$\pm$0.5'%(r_cut_off[k]))
                        l = pl.axvline(x=cdi_value,ls='--', color='black')
                        cdi_value=cdi[k][2]
                        pl.text(cdi_value-1.0, y_ax_lim[1]*0.15, r'%1.1f$\pm$0.1'%(cdi_value))
                elif((k=='C-Na')):
                        l = pl.axvline(x=r_cut_off[k],ls='--', color='black')
                        pl.text(r_cut_off[k]+0.15, y_ax_lim[1]*0.15, r'%1.1f$\pm$0.1'%(r_cut_off[k]))
                        cdi_value=cdi[k][0]
                        l = pl.axvline(x=cdi_value,ls='--', color='black')
                        pl.text(cdi_value-1.0, y_ax_lim[1]*0.15, r'%1.1f$\pm$0.1'%(cdi_value))
             
                fig.savefig(name_file_save, dpi=300)

name_file_save_txt=r'3_RDF_%03d_%03d_%s_%s_%s_%s.txt' %(runmin,runmax, list_path[8],list_path[9],list_path[10],list_path[11])
tmpf=open(name_file_save_txt, 'w')
tmpf.write(text)
tmpf.close()

