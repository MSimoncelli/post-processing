#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import subprocess
import shutil
#import datetime
#from pylab import *
from sys import argv

flag_untar=True
############################################################
# first and last runs to include in average density
#print 'START'
#print datetime.datetime.now()
if (len(argv)<=1):
        print 'USAGE: control_RDF.py runmin runmax'
        runmin = 0
        runmax = 20
        num_max_proc = 4
else:
	runmin=int(argv[1])
	runmax=int(argv[2])	
	num_max_proc=int(argv[3])

#build names as done in the calc_RDF_program 
spc_upd=[]
# names
spc_upd.append('O')
spc_upd.append('Na') #Na
spc_upd.append('Cl') #Cl
spc_upd.append('C')

mem=0
n_dist=4

int_RDF=[]
settings=[]
for k in range(0,n_dist):
	for i in range(0,2):
		int_RDF.append(spc_upd[k]+'-'+spc_upd[i+1])
		settings.append([spc_upd[k],spc_upd[i+1]])
		mem=mem+1


startdir = os.getcwd() #gets the current working directory
#print startdir

pid_list=[]
#print "running from", execpwd
sp=[]
indx_buf=0

#print 'start cycle irun'
#print datetime.datetime.now()
launched=0
for irun in range(runmin,runmax+1):
	# move to directory and uncompress files
	dirrun = '../run%03d' %( irun )
	os.chdir(dirrun)        #changes directory
	#copy the script in the run folder
	if flag_untar:
		cmd_list = ['tar', '-xzf', 'positions.tar.gz']
		sp.append(1)
		sp[indx_buf] = subprocess.Popen(cmd_list)
		indx_buf+=1

	cmd_list_1 = ['tar', '-xzf', 'wallq.tar.gz']
	sp.append(1)
	sp[indx_buf] = subprocess.Popen(cmd_list_1)
	indx_buf+=1
   
#print 'OUT 1 cycle irun'
#print datetime.datetime.now()

for j in range(0, indx_buf):
        sp[j].wait()
        #print 'pid=', pid_list[j]
        #print datetime.datetime.now() 
pid_list=[]
#print "running from", execpwd
sp=[]
indx_buf=0

#print 'start cycle irun'
#print datetime.datetime.now()
launched=0
for irun in range(runmin,runmax+1):
	# move to directory and uncompress files
	dirrun = '../run%03d' %( irun )
	os.chdir(dirrun)	#changes directory
	cmd = 'cp ../5_charge/calc_slice.py ./' 
	#print cmda			asdf,m  sadlknnaaasdssdasadadsfkjsdfnlkadssdaklmadslkm
	os.system(cmd)	#executes a command as it was typed in shell
	sp.append(1)
	sp[indx_buf] = subprocess.Popen(['python', 'calc_slice.py'])
	pid_list.append(sp[indx_buf].pid)
	indx_buf+=1
	launched+=1
	print launched
	if (launched>=num_max_proc):
		for i in range(num_max_proc,0,-1):
			sp[indx_buf-i].wait()
		launched=0
		
#print 'OUT 1 cycle irun'
#print datetime.datetime.now()

for j in range(0, indx_buf):
	sp[j].wait()
	#print 'pid=', pid_list[j]
	#print datetime.datetime.now() 
	#print 'after wait ', j

#now all the calculations have been carried on.
#we have to clean the runXXX directories and store all the data in the proper folder
sp=[]
indx_buf=0

#create suitable folder
cmd = 'mkdir -p '+startdir+'/data_OUT'
print cmd
os.system(cmd)

outfilename='tot_charge_run%03d_run%03d.dat'%(runmin,runmax)
here=os.getcwd()
os.chdir(startdir+'/data_OUT')
outfile=open(outfilename, 'wb')
os.chdir(here)

for irun in range(runmin,runmax+1):
	dirrun = '../run%03d' %( irun )
	os.chdir(dirrun)	#changes directory
	#clean the directory from the python script
	sp.append(1)
	sp[indx_buf] = subprocess.Popen(['rm', 'calc_slice.py'])
	indx_buf+=1

	#clean from the position.out
	sp.append(1)
	sp[indx_buf] = subprocess.Popen(['rm', 'wallq.out'])
	indx_buf+=1

	#clean from the position.out
        #sp.append(1)
        #sp[indx_buf] = subprocess.Popen(['rm', 'positions.out'])
        #indx_buf+=1
	
	readfile=open('charge_slice.dat', 'rb')
	shutil.copyfileobj(readfile, outfile)

outfile.close()
