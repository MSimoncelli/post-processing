import matplotlib.pyplot as pl
import numpy as np
import matplotlib.cm as cm
import fileinput
import math
import inspect, os
import scipy as sp
from scipy import signal
from sys import argv
from scipy.optimize import curve_fit

##AAAAAAAAAA
#remember to set the global variables Q_max bnefore calling the function!


def pwu(a,b):
        c=b
        i=0
        if (b<1 and b!=0):
	        while int(c)<1:
	                c=c*10
	                i+=1
        return '%.*f$\pm$%.*f'%(i,a,i,b)


'''
#UNITS:
CHARGE: e=1.602*10**-19 C
TIME: nanosecond: 10**-9 s
CAPACITY: 10**-18 F
'''
############################################################

t_0=10**(-9)
R_0=10**9
C_0=1.60217662*10**(-19)
V_0=1.0
Q_0=1.60217662*10**(-19)

a_red_conv=(t_0/(R_0*C_0))
b_red_conv=(t_0**2/(R_0**2*C_0**2))
c_red_conv=(t_0**2*V_0/((R_0**2)*Q_0*C_0))
d_red_conv=(t_0*V_0/(R_0*Q_0))

print a_red_conv, b_red_conv, c_red_conv, d_red_conv


def Q_function_1(t, R_l):
	V_0=1.0 #voltage potential
	C_1=2.0*Q_1max/V_0
	C_2=2.0*Q_2max/V_0
	a=a_red_conv*((R_bulk+2*R_l)*C_1+(R_bulk+4*R_l)*C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	b=b_red_conv*2.0/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	c_1=c_red_conv*1.0/(R_l*(R_bulk+2*R_l)*C_2)
	d=d_red_conv*1.0/(R_bulk+2*R_l)
	tau_1=(2.0/(a+math.sqrt(a**2-4.0*b)))
	tau_2=(2.0/(a-math.sqrt(a**2-4.0*b)))
	tildeA_1=0.5*(1.0+((2.0*b*d-a*c_1)/(2*c_1*math.sqrt(a**2-4*b))))
	tildeA_2=0.5*(1.0-((2.0*b*d-a*c_1)/(2*c_1*math.sqrt(a**2-4*b))))
	return Q_1max*(1-tildeA_1*np.exp(- t/tau_1) -tildeA_2*np.exp(- t/tau_2))

def Q_function_2(t, R_l):
	V_0=1.0 #voltage potential
	C_1=2.0*Q_1max/V_0
	C_2=2.0*Q_2max/V_0
	a=a_red_conv*((R_bulk+2*R_l)*C_1+(R_bulk+4*R_l)*C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	b=b_red_conv*(2.0/(R_l*(R_bulk+2*R_l)*C_1*C_2))
	c_2=c_red_conv*1.0/(R_l*(R_bulk+2*R_l)*C_1)
	d=0.0
	tau_1=(2.0/(a+math.sqrt(a**2-4.0*b)))
	tau_2=(2.0/(a-math.sqrt(a**2-4.0*b)))
	hatA_1=0.5*(1.0+((2.0*b*d-a*c_2)/(2*c_2*math.sqrt(a**2-4*b))))
	hatA_2=0.5*(1.0-((2.0*b*d-a*c_2)/(2*c_2*math.sqrt(a**2-4*b))))
	#print hatA_1, hatA_2
	return Q_2max*(1-hatA_1*np.exp(- t/tau_1) -hatA_2*np.exp(- t/tau_2))

#before calling this function you need to define the variablC_1, C_2 from the previous fits!!
def Q_function_tot(t, R_l):
	V_0=1.0 #voltage potential
	C_1=2.0*Q_1max/V_0
	C_2=2.0*Q_2max/V_0
	Q_max=Q_1max+Q_2max
	a=a_red_conv*((R_bulk+2*R_l)*C_1+(R_bulk+4*R_l)*C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	b=b_red_conv*2.0/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	c=c_red_conv*(C_1+C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)
	d=d_red_conv*1.0/(R_bulk+2*R_l)
	tau_1=(2.0/(a+math.sqrt(a**2-4.0*b)))
	tau_2=(2.0/(a-math.sqrt(a**2-4.0*b)))
	A_1=0.5*(1.0+((2.0*b*d-a*c)/(2*c*math.sqrt(a**2-4*b))))
	A_2=0.5*(1.0-((2.0*b*d-a*c)/(2*c*math.sqrt(a**2-4*b))))
	return Q_max * (1-A_1*np.exp(- t/tau_1) -A_2*np.exp(- t/tau_2))

def print_result_table(columns, cell_text, name_file_save_t):
	text='#'
	for col_label in columns:
		text+=col_label+'\t'
	text+='\n'	

	for list_ion in cell_text:
		for string in list_ion:
			text+=string.replace('$\pm$', '\t') +'\t'
		text+='\n'	

	tmpf=open(name_file_save_t,'w')
	tmpf.write(text)
	tmpf.close()


if (len(argv)<=1):
        #print 'USAGE: plor_coord_num.py runmin runmax'
        runmin = 0
        runmax = 62
else:
	runmin=int(argv[1])
	runmax=int(argv[2])

name_folder=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
list_dir=name_folder.split('/')
#print list_dir

name_file_charge='tot_charge_run%03d_run%03d.dat'%(runmin,runmax)

from matplotlib import rc
rc('font',**{'family':'Texgyre'})
rc('text', usetex=True)
##############################################################
data_charge=np.loadtxt(name_file_charge)

#electrode 1, see the image in the latex file
data_slice_1=data_charge[:,0]
data_slice_2=data_charge[:,1]


data_left_charge=data_slice_1+data_slice_2

data_slice_3=data_charge[:,2]
data_slice_4=data_charge[:,3]
data_right_charge=data_slice_3+data_slice_4

time_step=10**(-3)#=1 fs= 10^(-6) nanosecond. The charge is recorded every 1000 step=every ps!
tot_num_data=len(data_left_charge)


num_data_plateau=1000

Q_2max1=-np.mean(data_slice_1[-num_data_plateau:])
inc_Q_2max1=(Q_2max1-min(-data_slice_1[-num_data_plateau:]))

Q_1max2=-np.mean(data_slice_2[-num_data_plateau:])
inc_Q_1max2=Q_1max2-min(-data_slice_2[-num_data_plateau:])

Q_1max3=np.mean(data_slice_3[-num_data_plateau:])
inc_Q_1max3=Q_1max3-min(data_slice_3[-num_data_plateau:])

Q_2max4=np.mean(data_slice_4[-num_data_plateau:])
inc_Q_2max4=Q_2max4-min(data_slice_4[-num_data_plateau:])

Q_2max1+=inc_Q_2max1
Q_1max2+=inc_Q_1max2
Q_1max3+=inc_Q_1max3
Q_2max4+=inc_Q_2max4

V_0=1.0


time_array=np.arange(0,int(tot_num_data))
time_array=(time_array*time_step)
#print list_dir[11]
if (list_dir[9]=='nacl'):
	spc_name='NaCl'
else:
	spc_name='KCl'

if ((list_dir[10]=='8800_160') or (list_dir[10]=='7615_139')):
	conC='1.0 M'
else:
	conC='0.5 M'

if spc_name=='NaCl':
    if (list_dir[10]=='8800_80'): 
        R_bulk=0.0988443862
        N_carbon=3649
    if (list_dir[10]=='7700_70'): 
        R_bulk=0.0972311627
        N_carbon=3821
    if (list_dir[10]=='8800_160'):
        R_bulk=0.0541421056
        N_carbon=3649
    if (list_dir[10]=='7615_139'): 
        R_bulk=0.0551470421  
        N_carbon=3821
if spc_name=='KCl':
    if (list_dir[10]=='8800_80'): 
        R_bulk=0.0832941274
        N_carbon=3649
    if (list_dir[10]=='8800_160'):
        R_bulk=0.0437107912
        N_carbon=3649


#pcov : 2d array
#http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
#The estimated covariance of popt. 
#The diagonals provide the variance of the parameter estimate. 
#To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov)).

fig=pl.figure()
ax1 = fig.add_subplot(111)


title= 'Charge vs time %s %s %s %s' %(list_dir[8],spc_name,conC,list_dir[11] )
title=title.replace('_', '-')
name_file_save='5_Charge_separed_run%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
name_file_save=name_file_save.replace('.','_')
name_file_save_text=name_file_save
name_file_save+='.pdf'

ax1.set_title(title)
ax1.grid(True)
# ax1.set_xlim([-50,50])
ax1.set_xlabel(r'time $[ns]$')
ax1.set_ylabel(r'charge [e]')

ax1.plot(time_array,data_slice_1, label='Q2, slice 1' ,color='#00A6FF')
ax1.plot(time_array,data_slice_2, label='Q1, slice 2' ,color='b')
ax1.plot(time_array,data_slice_3, label='Q1, slice 3' ,color='r')
ax1.plot(time_array,data_slice_4, label='Q2, slice 4' ,color='orange')


if (list_dir[11]=='1.0V'):
	#do not confuse the slice number and the index of the capacitor!!
	Q_1max=Q_1max2
	Q_2max=Q_2max1
	Rl_1, cov_Rl_1 = curve_fit(Q_function_2, time_array, -data_slice_1, p0=(0.06))#,bounds=([0.0],[1.0000]), method='trf')
	fitted_slice_1=(-1.0)*Q_function_2(time_array,Rl_1[0])
	
	Rl_2, cov_Rl_2 = curve_fit(Q_function_1, time_array,-data_slice_2, p0=(0.06))#,bounds=([0.0],[1.0000]), method='trf')
	fitted_slice_2=(-1.0)*Q_function_1(time_array,Rl_2[0])
	
	Q_1max=Q_1max3
	Q_2max=Q_2max4
	Rl_3, cov_Rl_3 = curve_fit(Q_function_1, time_array,data_slice_3, p0=(0.06))#,bounds=([0.0],[1.0000]), method='trf')
	fitted_slice_3=Q_function_1(time_array,Rl_3[0])

	Rl_4, cov_Rl_4 = curve_fit(Q_function_2, time_array, data_slice_4, p0=(0.06))#,bounds=([0.0],[1.0000]), method='trf')
	fitted_slice_4=Q_function_2(time_array,Rl_4[0])

	list_R_l=[Rl_1[0],Rl_2[0],Rl_3[0],Rl_4[0]]
	list_inc_R_l=[(cov_Rl_1[0])**0.5,(cov_Rl_2[0])**0.5,(cov_Rl_3[0])**0.5,(cov_Rl_4[0])**0.5]
	w_sum=0
	w_inc_rec=0
	mean_val=np.mean(list_R_l)
	stad_val=np.std(list_R_l)

	mean_Q_slice1=np.mean([Q_1max3,Q_1max2])
	disp_Q_slice1=0.5*abs(Q_1max3-Q_1max2)
	mean_Q_slice2=np.mean([Q_2max1,Q_2max4])
	disp_Q_slice2=0.5*abs(Q_2max1-Q_2max4)

	for i in range(0,len(list_R_l)):
		w_sum+=list_R_l[i]*(1.0/(list_inc_R_l[i]**2))
		w_inc_rec+=(1.0/(list_inc_R_l[i]**2))

	avg_tot=w_sum/w_inc_rec	
	inc_tot=(1.0/w_inc_rec)**0.5

	ax1.plot(time_array,fitted_slice_1, label='fit - Q2, slice 1' ,color='b')
	ax1.plot(time_array,fitted_slice_2, label='fit - Q1, slice 2' ,color='#2FFF00')
	ax1.plot(time_array,fitted_slice_3, label='fit - Q1, slice 3' ,color='purple')
	ax1.plot(time_array,fitted_slice_4, label='fit - Q2, slice 4' ,color='#FF00FC')

	rows = [' ', 'slice 1', 'slice 2', 'slice 3', 'slice 4']
	columns = ('Fit results', 'Qm $[e]$' , r'Rl  $\;[10^7\cdot \Omega]$')
	cell_text=[['slice 1',pwu(Q_2max1,inc_Q_2max1),pwu(Rl_1[0]*100,((cov_Rl_1[0])**0.5)*100)],
			   ['slice 2',pwu(Q_1max2,inc_Q_1max2),pwu(Rl_2[0]*100,((cov_Rl_2[0])**0.5)*100)],
			   ['slice 3',pwu(Q_1max3,inc_Q_1max3),pwu(Rl_3[0]*100,((cov_Rl_3[0])**0.5)*100)],
			   ['slice 4',pwu(Q_2max4,inc_Q_2max4),pwu(Rl_4[0]*100,((cov_Rl_4[0])**0.5)*100)],
			   ['Rl w. avg',' -- ',pwu(avg_tot*100,inc_tot*100)],
			   ['Rl avg, std','--  ',pwu(mean_val*100,stad_val*100)],
			   ['Q1m (avg s. 2,3)',pwu(mean_Q_slice1,disp_Q_slice1),' -- ' ],
			   ['Q2m (avg s. 1,4)',pwu(mean_Q_slice2,disp_Q_slice2),' -- ' ]
			   ]
	art = []
	the_table = pl.table(cellText=cell_text,
                      #colColours=colors,
                      colLabels=columns,
                      cellLoc = 'center', rowLoc = 'center',
					  loc='bottom', bbox=[0.1, -0.8, 0.8, 0.64])	
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(11)
	art.append(the_table)	
	print_result_table(columns,cell_text,name_file_save_text+'.txt')

art = []
lgd = ax1.legend(loc=9, bbox_to_anchor=(1.2, 1), ncol=1)
art.append(lgd)	
fig.set_size_inches(6, 4.8)
fig.savefig(name_file_save, dpi=150,additional_artists=art, bbox_inches="tight")

fig=pl.figure()
ax1 = fig.add_subplot(111)

title= 'Charge vs time %s %s %s %s' %(list_dir[8],spc_name,conC,list_dir[11] )
title=title.replace('_', '-')
name_file_save='5_Charge_run%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
name_file_save=name_file_save.replace('.','_')

#print title, 'TOT_NUM_DATA=',tot_num_data
ax1.set_title(title)
ax1.grid(True)
# ax1.set_xlim([-50,50])
ax1.set_xlabel(r'time $[ns]$')
ax1.set_ylabel(r'charge [e]')

ax1.plot(time_array,data_right_charge, label='right electrode ' ,color='r')
ax1.plot(time_array,data_left_charge, label='left electrode' ,color='b')	

#print 'R_bulk=', R_bulk

name_file_save_text=name_file_save
name_file_save+='.pdf'


if (list_dir[11]=='1.0V'):
	Q_max_el1=max(data_left_charge*(-1.0))
	Q_max_el2=max(data_right_charge)
	inc_Q_max_el1=Q_max_el1-min(data_left_charge[-30:]*(-1))
	inc_Q_max_el2=Q_max_el2-min((data_right_charge[-30:]))

	w1=1.0/inc_Q_max_el1**2
	w2=1.0/inc_Q_max_el2**2
	avg_Q_max=(w1*Q_max_el1+w2*Q_max_el2)/(w1+w2)
	inc_avg_Q_max=(1/(w1+w2)**0.5)

	Q_1max=Q_1max2
	Q_2max=Q_2max1
	Rl_el1, cov_Rl_el1 = curve_fit(Q_function_tot, time_array, data_left_charge*(-1))#,p0=(1),bounds=([1*low_bound],[R_l*up_bound]), method='trf')
	fitted_1=Q_function_tot(time_array,Rl_el1[0])	

	Q_1max=Q_1max3
	Q_2max=Q_2max4
	Rl_el2, cov_Rl_el2 = curve_fit(Q_function_tot, time_array, data_right_charge)#*(-1),p0=(1),bounds=([1*low_bound],[R_l*up_bound]), method='trf')
	fitted_2=Q_function_tot(time_array,Rl_el2[0])

	sigma_Rl_el1 = np.sqrt(np.diag(cov_Rl_el1))
	sigma_Rl_el2 = np.sqrt(np.diag(cov_Rl_el2))

	average=(Rl_el1[0]*(1.0/sigma_Rl_el1)**2+Rl_el2[0]*(1.0/sigma_Rl_el2)**2)/((1.0/sigma_Rl_el1)**2+(1.0/sigma_Rl_el2)**2)
	inc_avg=1.0/((1.0/sigma_Rl_el1)**2+(1.0/sigma_Rl_el2)**2)**0.5

	mean_val=np.mean([Rl_el1[0],Rl_el2[0]])
	disp_val=abs(Rl_el1[0]-Rl_el2[0])*0.5

	mean_Q_max=0.5*(Q_max_el1+Q_max_el2)
	semidispQ_max=0.5*abs(Q_max_el1-Q_max_el2)
	
	ax1.plot(time_array,fitted_1*(-1.0), label='fit - el 1' ,color='c')
	ax1.plot(time_array,fitted_2, label='fit - el 2' ,color='orange')	
	rows = [' ', 'electrode 1', 'electrode 2']
	columns = ('Fit results',  'Qm $[e]$' ,'Rl  $\;[10^7\cdot \Omega]$')
	cell_text=[['electrode 1',pwu(Q_max_el1,inc_Q_max_el1),pwu(Rl_el1[0]*100,sigma_Rl_el1[0]*100)],
			   ['electrode 2 ',pwu(Q_max_el2,inc_Q_max_el2),pwu(Rl_el2[0]*100,sigma_Rl_el2[0]*100)],
			   ['Rl w. avg', '--' ,pwu(average*100,inc_avg*100)],
			   ['Rl avg, sdisp', '--',pwu(mean_val*100,disp_val*100)],
			   ['Qm (w. avg)',pwu(avg_Q_max,inc_avg_Q_max),' -- ' ],
			   ['Qm (avg, sdisp)',pwu(mean_Q_max,semidispQ_max),' -- ' ]
			   ]
	art = []
	the_table = pl.table(cellText=cell_text,
                      #colColours=colors,
                      colLabels=columns,
                      cellLoc = 'center', rowLoc = 'center',
					  loc='bottom', bbox=[0.1, -0.56, 0.8, 0.40])	
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(11)
	art.append(the_table)	

	print_result_table(columns,cell_text,name_file_save_text+'.txt')



ax1.legend(loc='best')
fig.set_size_inches(6, 4.8)
fig.savefig(name_file_save, dpi=150,additional_artists=art, bbox_inches="tight")

