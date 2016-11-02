#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
#import datetime
#from pylab import *
from sys import argv
############################################################
# first and last runs to include in average density
#print 'START'
#print datetime.datetime.now()
if (len(argv)<=2):
        print 'USAGE: control_RDF.py runmin runmax'
        runmin = 0
        runmax = 46
        name_cutoff_file = 'A_cutoff_dataset.txt'
else:
	runmin=int(argv[1])
	runmax=int(argv[2])	
	name_cutoff_file=str(argv[3]) #give as argument the file containing the cutoffs

execpwd=os.getcwd()
pid_list=[]
#print "running from", execpwd
proc_id=[]
indx_buf=0

#print 'start cycle irun'
#print datetime.datetime.now()

for irun in range(runmin,runmax+1):
	# move to directory and uncompress files
	dirrun = '../run%03d' %( irun )
	os.chdir(dirrun)	#changes directory
	#copy the script in the run folder
	cmd = 'tar -xzf positions.tar.gz'
	os.system(cmd)
	cmd = 'cp ../4_coord_number/coord_num.x ./' 
	os.system(cmd)
        cmd = 'cp ../4_coord_number/'+name_cutoff_file+' ./'
        os.system(cmd)  
	proc_id.append(1)
	proc_id[indx_buf] = subprocess.Popen(['./coord_num.x',name_cutoff_file])
	pid_list.append(proc_id[indx_buf].pid)
	indx_buf+=1

#print 'OUT 1 cycle irun'
#print datetime.datetime.now()

for j in range(0, indx_buf):
	proc_id[j].wait()
	#print 'pid=', pid_list[j]
	#print datetime.datetime.now() 
	#print 'after wait ', j

#now all the calculations have been carried on.
#we have to clean the runXXX directories and store all the data in the proper folder
proc_id=[]
indx_buf=0

#create suitable folder
cmd = 'mkdir -p ../4_coord_number/data_OUT'
os.system(cmd)

for irun in range(runmin,runmax+1):
	dirrun = '../run%03d' %( irun )
	os.chdir(dirrun)	#changes directory

	#clean from the position.out
	proc_id.append(1)
	proc_id[indx_buf] = subprocess.Popen(['rm', 'positions.out'])
	indx_buf+=1

	proc_id.append(1)
	proc_id[indx_buf] = subprocess.Popen(['rm', 'coord_num.x'])
	indx_buf+=1

        proc_id.append(1)
        proc_id[indx_buf] = subprocess.Popen(['rm', name_cutoff_file])
        indx_buf+=1

	cmd='rename '+r" '"+r's/^/'+'run%03d'%(irun)+r'_/'r"' "+r'out_coord_*'
	os.system(cmd)

	start_name='run%03d*'%(irun)
	cmd = 'mv '+start_name+' ../4_coord_number/data_OUT'
	os.system(cmd)

os.chdir(execpwd)
cmd = 'mv '+name_cutoff_file+' data_OUT'
os.system(cmd)
