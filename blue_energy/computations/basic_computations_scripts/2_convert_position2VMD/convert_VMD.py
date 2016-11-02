#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
#from pylab import *
from sys import argv
############################################################
# first and last runs to include in average density

if (len(argv)<=1):
        print 'USAGE: visualize_VMD.py runmin runmax'
        runmin = 0
        runmax = 53
else:
        runmin=int(argv[1])
        runmax=int(argv[2])

#print 'runmin=',runmin, ' runmax=', runmax
############################################################

conv=0.529177 #convert atomic units to Angstrom
startdir = os.getcwd() #gets the current working directory
#print "running from", execpwd
os.chdir(startdir)
info=open("info_convert2VMD.txt",'w')
info.write('#'+startdir+'\n')

"""
this script reads the file positions.out in each folder, convert it in traj.xyz and does the same operation in all the folders
it builds the file pm2xyz-info.template which is then used by the conversion program.
In the same folder from which you run this script you need the folder "conv_p2vmd" containing:
- convert2xyz.sh
- pm2xyz.x
- pm2xyz-inpt
- pm2xyz-inpt.template
(suorce pm2xyz.f90)
the script goes in each folder runXXX (range specified by runmin-runmax) decompress position.tar.gz, creates the template file, convert position file in traj.xyz file
and delete the uncompacted files.
"""

def read_info_specie(line_num_atom, line_name, file_positions):
	line_specie = file_positions[line_num_atom].split()
	Nspec=int(line_specie[0])
	line_name_specie=file_positions[line_name].split()
	name_spec=line_name_specie[0]
	return Nspec, name_spec

list_Nspec=[]
list_name_spec=[]

dirrun = '../run%03d' %( runmin )
os.chdir(dirrun)
file_runtime_IN= open( 'runtime.inpt', 'r' )
file_runtime_a=file_runtime_IN.readlines() # read all the file, it contains lots of useful things
Nspec_1,name_spec_1=read_info_specie(7,6,file_runtime_a)
#print "Nspec_1=%d, name_spec_1=%f" %(Nspec_1,name_spec_1)
list_Nspec.append(Nspec_1)
list_name_spec.append(name_spec_1)      

Nspec_2,name_spec_2=read_info_specie(12,11,file_runtime_a)
#print "Nspec_2=%d, name_spec_2=%f" %(Nspec_2,name_spec_2)
#this is Hidrogen, multiply by 2
list_Nspec.append(Nspec_2*2)
list_name_spec.append("H")      

Nspec_3,name_spec_3=read_info_specie(22,21,file_runtime_a)
#print "Nspec_3=%d, name_spec_3=%f" %(Nspec_3,name_spec_3)
list_Nspec.append(Nspec_3)
list_name_spec.append(name_spec_3)      

Nspec_4,name_spec_4=read_info_specie(27,26,file_runtime_a)
#print "Nspec_4=%d, name_spec_4=%f" %(Nspec_4,name_spec_4)
list_Nspec.append(Nspec_4)
list_name_spec.append(name_spec_4)      

Nspec_5,name_spec_5=read_info_specie(32,31,file_runtime_a)
Nspec_6,name_spec_6=read_info_specie(42,41,file_runtime_a)
#print "Nspec_5=%d, name_spec_5=%f" %(Nspec_5,name_spec_5)
#those are the carbons
list_Nspec.append(Nspec_5*2+Nspec_6)
list_name_spec.append("C")      

nskip=1
nspecies=len(list_Nspec)
os.chdir(startdir)
cmd='mkdir -p data_OUT'
os.system(cmd)

os.chdir('conv_p2vmd')
file_template= open( 'pm2xyz-inpt.template', 'w' )
file_template.write('nframes\n')
file_template.write('%d\n'%(nskip))
file_template.write('%d\n'%(nspecies))
for x in xrange(0,nspecies):
	file_template.write('%s\n'%(list_name_spec[x]))
	file_template.write('%d\n'%(list_Nspec[x]))
file_template.write('positions.out\n')
file_template.write('traj.xyz\n')
file_template.write('%.6f' %(conv))
file_template.close()


for irun in range(runmin,runmax+1):
	print 'run %d out of %d' %(irun-runmin, runmax-runmin)
	# move to directory and uncompress files
	dirrun = '../../run%03d' %( irun )
	os.system('cp pm2xyz-inpt '+dirrun)
	os.system('cp pm2xyz-inpt.template '+dirrun)
	os.system('cp pm2xyz.x '+dirrun)
	os.system('cp convert2xyz.sh '+dirrun)

	os.chdir(dirrun)	#changes directory
	cmd = 'tar -xzf positions.tar.gz'
	os.system(cmd)	#executes a command as it was typed in shell
	cmd = 'bash ./convert2xyz.sh positions.out'
	os.system(cmd)
	# open positions.out

	cmd='mv traj.xyz ../2_convert_position2VMD/data_OUT/traj_%03d.xyz'%(irun) 
	os.system(cmd)
		
	# clean directory and come back
	cmd = 'rm *.out'
	#print cmd
	os.system(cmd)
	cmd = 'rm pm2xyz*'
	os.system(cmd)
	cmd = 'rm convert2xyz.sh'
	os.system(cmd)
	
	os.chdir(startdir)	
	os.chdir('conv_p2vmd')

