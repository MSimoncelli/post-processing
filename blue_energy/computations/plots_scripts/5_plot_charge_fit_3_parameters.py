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
import lmfit
from lmfit import minimize, Parameters, report_fit

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

def pwu_list(lista):
    a=lista[0]
    b=lista[1]
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

mass_carbon=12.0107*1.66054*10**(-24)
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


def Q_function_1(t, R_l, R_bulk, Q_1max, Q_2max):
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

def Q_function_2(t, R_l, R_bulk, Q_1max, Q_2max):
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
def Q_function_tot(t, R_l, R_bulk, Q_1max, Q_2max):
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

def weight_average(list_num, list_unc):
    w_avg=0
    err_w=0
    for indx in range(0, len(list_num)):
        w_avg+=list_num[indx]*(1.0/list_unc[indx])**2
        err_w+=(1.0/list_unc[indx])**2
    w_avg/=err_w
    err_w=1.0/(err_w)**(0.5)
    return [w_avg, err_w]

def mean_stdev(list_num, list_unc):
    mean_val=np.mean(list_R_l)
    stad_val=np.std(list_R_l)
    return [mean_val, stad_val]

def mean_variance(list_num, list_unc):
    mean_val=np.mean(list_R_l)
    tot=0
    for stddev in list_unc:
        tot+=stddev**2
    stad_val=tot**0.5
    return [mean_val, stad_val]



def Q_dataset_4(params, i, t):
    """calc gaussian from params for data set i
    using simple, hardwired naming convention"""
    R_l = params['R_l'].value
    R_bulk = params['R_bulk'].value
    Q_1max = params['Q_1max'].value
    Q_2max = params['Q_2max'].value
    if ((i==0) or (i==1)):
        res=Q_function_1(t,  R_l, R_bulk, Q_1max, Q_2max)
    elif ((i==2)or(i==3)):
        res=Q_function_2(t,  R_l, R_bulk, Q_1max, Q_2max)
    return res

def Q_dataset_2(params, i, t):
    """calc gaussian from params for data set i
    using simple, hardwired naming convention"""
    R_l = params['R_l'].value
    R_bulk = params['R_bulk'].value
    Q_1max = params['Q_1max'].value
    Q_2max = params['Q_2max'].value
    if (i==0):
        res=Q_function_1(t,  R_l, R_bulk, Q_1max, Q_2max)
    elif (i==1):
        res=Q_function_2(t,  R_l, R_bulk, Q_1max, Q_2max)
    return res


def objective_Q(params, t, data):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by Gaussian functions"""
    ndata, nx = data.shape
    resid = 0.0*data[:]
    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - Q_dataset_2(params, i, t)
    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

def fcn2min_Q1(params, t, data):
    """ model for Q1, subtract data"""
    R_l = params['R_l_1'].value
    R_bulk = params['R_bulk_1'].value
    Q_1max = params['Q_1max_1'].value
    Q_2max = params['Q_2max_1'].value
    model = Q_function_1(t,  R_l, R_bulk, Q_1max, Q_2max)
    return model - data


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
data_slice_3=data_charge[:,2]
data_slice_4=data_charge[:,3]

data_left_charge=data_slice_1+data_slice_2
data_right_charge=data_slice_3+data_slice_4
data_array_el2=[data_slice_3,data_slice_4]
data_array_el2 = np.array(data_array_el2)
data_array_el1=[(-1.0)*data_slice_2,(-1.0)*data_slice_1]
data_array_el1 = np.array(data_array_el1)

if (list_dir[9]=='nacl'):
    spc_name='NaCl'
    factor_not_plateau=1.4
else:
    spc_name='KCl'
    factor_not_plateau=2.0

if ((list_dir[10]=='8800_160') or (list_dir[10]=='7615_139')):
    conC='1.0 M'
else:
    conC='0.5 M'

time_step=10**(-3)#=1 fs= 10^(-6) nanosecond. The charge is recorded every 1000 step=every ps!
tot_num_data=len(data_slice_1)

time_array=np.arange(0,int(tot_num_data))
time_array=(time_array*time_step)


fig=pl.figure()
ax1 = fig.add_subplot(111)


title= 'Charge vs time %s %s %s %s' %(list_dir[8],spc_name,conC,list_dir[11] )
title=title.replace('_', '-')
name_file_save='5_Charge_REV_sep_run%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
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
    num_data_plateau=1000   
    Q_2max1=-np.mean(data_slice_1[-num_data_plateau:])
    inc_Q_2max1=(Q_2max1-min(-data_slice_1[-num_data_plateau:]))
    Q_1max2=-np.mean(data_slice_2[-num_data_plateau:])
    inc_Q_1max2=Q_1max2-min(-data_slice_2[-num_data_plateau:])
    Q_1max3=np.mean(data_slice_3[-num_data_plateau:])
    inc_Q_1max3=Q_1max3-min(data_slice_3[-num_data_plateau:])
    Q_2max4=np.mean(data_slice_4[-num_data_plateau:])
    inc_Q_2max4=Q_2max4-min(data_slice_4[-num_data_plateau:])       
    
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

    mass_tot_C=mass_carbon*N_carbon
    percentage_plateau=0.2
    print inc_Q_1max3
    print inc_Q_1max2

    params_1 = lmfit.Parameters()
    params_1.add( 'R_l', value=0.03910899, min=0.0,  max=0.2, vary=True)
    params_1.add( 'R_bulk', value=R_bulk, min=R_bulk*0.9,  max=R_bulk*1.1, vary=False)
    params_1.add( 'Q_1max', value=Q_1max2, min=Q_1max2-1.*inc_Q_1max2, max=Q_1max2+5.*inc_Q_1max2, vary=True) 
    params_1.add( 'Q_2max', value=Q_2max1, min=Q_2max1-1.*inc_Q_2max1, max=Q_2max1+5.*inc_Q_2max1, vary=True) 

    params_2 = lmfit.Parameters()
    params_2.add( 'R_l', value=0.05615835, min=0.0,  max=0.2, vary=True)
    params_2.add( 'R_bulk', value=R_bulk, min=R_bulk*0.9,  max=R_bulk*1.1, vary=False)
    params_2.add( 'Q_1max', value=Q_1max3, min=Q_1max3-1.*inc_Q_1max3, max=Q_1max3+5.*inc_Q_1max3) 
    params_2.add( 'Q_2max', value=Q_2max4, min=Q_2max4-1.*inc_Q_2max4, max=Q_2max4+5.*inc_Q_2max4) 

    #print 'electrode 1 method 1'
    #result_1=minimize(objective_Q, params_1, args=(time_array, data_array_el1))
    #report_fit(result_1.params)

    def residual_arr_1(params):
        return objective_Q(params, time_array, data_array_el1)

    def residual_arr_2(params):
        return objective_Q(params, time_array, data_array_el2)

    print 'electrode 1 method 2'
    mini_1 = lmfit.Minimizer(residual_arr_1, params_1)
    result_CC_1_a = mini_1.minimize(method='Nelder')
    result_1 = mini_1.minimize(method='leastsq', params=result_CC_1_a.params)
    print(lmfit.fit_report(result_1.params))

    #print 'electrode 2'
    #result_2=minimize(objective_Q, params_2, args=(time_array, data_array_el2))
    #report_fit(result_2.params)

    print 'electrode 2 method 2'
    mini_2 = lmfit.Minimizer(residual_arr_2, params_2)
    result_CC_a = mini_2.minimize(method='Nelder')
    result_2 = mini_2.minimize(method='leastsq', params=result_CC_a.params)


    print(lmfit.fit_report(result_2.params))
    #ci_2 = lmfit.conf_interval(mini_2, result_CC_2)
    #lmfit.printfuncs.report_ci(ci_2)

    R_l_1= result_1.params['R_l'].value
    R_bulk= result_1.params['R_bulk'].value
    Q_1max1= result_1.params['Q_1max'].value
    Q_2max1= result_1.params['Q_2max'].value

    inc_Rl_1=(float(str(result_1.params['R_l']).split('+/-')[1].split(',')[0]))
    #inc_R_bulk_1=(float(str(result_1.params['R_bulk']).split('+/-')[1].split(',')[0]))
    inc_Q_1max1=(float(str(result_1.params['Q_1max']).split('+/-')[1].split(',')[0]))
    inc_Q_2max1=(float(str(result_1.params['Q_2max']).split('+/-')[1].split(',')[0]))

    R_l_2= result_2.params['R_l'].value
    R_bulk= result_2.params['R_bulk'].value
    Q_1max2= result_2.params['Q_1max'].value
    Q_2max2= result_2.params['Q_2max'].value

    inc_Rl_2=(float(str(result_2.params['R_l']).split('+/-')[1].split(',')[0]))
    #inc_R_bulk_2=(float(str(result_2.params['R_bulk']).split('+/-')[1].split(',')[0]))
    inc_Q_1max2=(float(str(result_2.params['Q_1max']).split('+/-')[1].split(',')[0]))
    inc_Q_2max2=(float(str(result_2.params['Q_2max']).split('+/-')[1].split(',')[0]))


    fitted_slice_1=(-1.0)*Q_function_2(time_array,R_l_1,R_bulk,Q_1max1,Q_2max1)
    fitted_slice_2=(-1.0)*Q_function_1(time_array,R_l_1,R_bulk,Q_1max1,Q_2max1)
    fitted_slice_3=Q_function_1(time_array,R_l_2,R_bulk,Q_1max2,Q_2max2)
    fitted_slice_4=Q_function_2(time_array,R_l_2,R_bulk,Q_1max2,Q_2max2)
    
    ax1.plot(time_array,fitted_slice_1, label='fit - Q2, slice 1' ,color='b')
    ax1.plot(time_array,fitted_slice_2, label='fit - Q1, slice 2' ,color='#2FFF00')
    ax1.plot(time_array,fitted_slice_3, label='fit - Q1, slice 3' ,color='purple')
    ax1.plot(time_array,fitted_slice_4, label='fit - Q2, slice 4' ,color='#FF00FC')

    '''
    #print weight_average(list_R_bulk,list_inc_R_bulk)[0]
    rows = [' ', 'slice 1', 'slice 2']
    columns =('Fit results','$R_l \;[10^7\cdot \Omega]$','$R_{bulk}  \;[10^7\cdot \Omega]$', 'Q$_{max1}$ $[e]$','Q$_{max2}$ $[e]$' )
    cell_text=[['electrode 1',  pwu(R_l_1*100,inc_Rl_1*100),'%10.6f'%(R_bulk*100),pwu(Q_1max1,inc_Q_1max1),pwu(Q_2max1,inc_Q_2max1)],
               ['electrode 2',  pwu(R_l_2*100,inc_Rl_2*100),'%10.6f'%(R_bulk*100),pwu(Q_1max2,inc_Q_1max2),pwu(Q_2max2,inc_Q_2max2)]
               ]
    art = []
    '''

    V_0=1.0
    C_1_1=2.0*Q_1max1/V_0
    C_2_1=2.0*Q_2max1/V_0

    unc_C_1=((inc_Q_1max1**2+inc_Q_2max1**2)**0.5)*C_0/(V_0*mass_tot_C)
    C_tot_1=0.5*(C_1_1+C_2_1)*C_0/mass_tot_C

    C_1_2=2.0*Q_1max2/V_0
    C_2_2=2.0*Q_2max2/V_0

    unc_C_2=((inc_Q_1max2**2+inc_Q_2max2**2)**0.5)*C_0/(V_0*mass_tot_C)
    C_tot_2=0.5*(C_1_2+C_2_2)*C_0/mass_tot_C


    rows = [' ', 'slice 1', 'slice 2']
    columns =('Fit results','$R_l \;[10^7\cdot \Omega]$','$R_{bulk}  \;[10^7\cdot \Omega]$', 'Q$_{max1}$ $[e]$','Q$_{max2}$ $[e]$', 'C$_{tot}$ '+r' $[\frac{F}{g}]$' )
    cell_text=[['electrode 1',  pwu(R_l_1*100,inc_Rl_1*100),'%10.6f'%(R_bulk*100),pwu(Q_1max1,inc_Q_1max1),pwu(Q_2max1,inc_Q_2max1),pwu(C_tot_1,unc_C_1)],
               ['electrode 2',  pwu(R_l_2*100,inc_Rl_2*100),'%10.6f'%(R_bulk*100),pwu(Q_1max2,inc_Q_1max2),pwu(Q_2max2,inc_Q_2max2),pwu(C_tot_2,unc_C_2)]
               ]

    art = []
    the_table = pl.table(cellText=cell_text,
                      #colColours=colors,
                      colLabels=columns,
                      cellLoc = 'center', rowLoc = 'center',
                      loc='bottom', bbox=[0, -0.4, 1.35, 0.28])    
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
name_file_save='5_Charge_REV_run%03d-run%03d_%s_%s_%s_%s' %(runmin,runmax, list_dir[8],list_dir[9],list_dir[10],list_dir[11])
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
    fitted_1=Q_function_tot(time_array,R_l_1 ,R_bulk,Q_1max1,Q_2max1)
    fitted_2=Q_function_tot(time_array,R_l_2 ,R_bulk,Q_1max2,Q_2max2)
    
    ax1.plot(time_array,fitted_1*(-1.0), label='fit - el 1' ,color='c')
    ax1.plot(time_array,fitted_2, label='fit - el 2' ,color='orange')   

    V_0=1.0
    C_1_1=2.0*Q_1max1/V_0
    C_2_1=2.0*Q_2max1/V_0

    unc_C_1=((inc_Q_1max1**2+inc_Q_2max1**2)**0.5)*C_0/(V_0*mass_tot_C)
    C_tot_1=0.5*(C_1_1+C_2_1)*C_0/mass_tot_C

    C_1_2=2.0*Q_1max2/V_0
    C_2_2=2.0*Q_2max2/V_0

    unc_C_2=((inc_Q_1max2**2+inc_Q_2max2**2)**0.5)*C_0/(V_0*mass_tot_C)
    C_tot_2=0.5*(C_1_2+C_2_2)*C_0/mass_tot_C


    rows = [' ', 'slice 1', 'slice 2']
    columns =('Fit results','$R_l \;[10^7\cdot \Omega]$','$R_{bulk}  \;[10^7\cdot \Omega]$', 'Q$_{max1}$ $[e]$','Q$_{max2}$ $[e]$', 'C$_{tot}$ '+r' $[\frac{F}{g}]$' )
    cell_text=[['electrode 1',  pwu(R_l_1*100,inc_Rl_1*100),'%10.6f'%(R_bulk*100),pwu(Q_1max1,inc_Q_1max1),pwu(Q_2max1,inc_Q_2max1),pwu(C_tot_1,unc_C_1)],
               ['electrode 2',  pwu(R_l_2*100,inc_Rl_2*100),'%10.6f'%(R_bulk*100),pwu(Q_1max2,inc_Q_1max2),pwu(Q_2max2,inc_Q_2max2),pwu(C_tot_2,unc_C_2)]
               ]

    art = []
    the_table = pl.table(cellText=cell_text,
                      #colColours=colors,
                      colLabels=columns,
                      cellLoc = 'center', rowLoc = 'center',
                      loc='bottom', bbox=[-0.15, -0.4, 1.2, 0.28])   
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(11)
    art.append(the_table)   
    print_result_table(columns,cell_text,name_file_save_text+'.txt')
 
ax1.legend(loc='best')
fig.set_size_inches(6, 4.8)
fig.savefig(name_file_save, dpi=150,additional_artists=art, bbox_inches="tight")
