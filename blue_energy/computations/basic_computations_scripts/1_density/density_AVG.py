#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
#from pylab import *
from sys import argv
############################################################
# first and last runs to include in average density

if (len(argv)<=1):
        print 'USAGE: python density_position.py runmin runmax'
        runmin = 0
        runmax = 53
else:
	runmin=int(argv[1])
	runmax=int(argv[2])	
#print 'runmin=',runmin, ' runmax=', runmax
############################################################
# density.out: this file contains the average density (as a function of the coordinate z) over 1 run, i.e. 50000 steps. 

Nbins = 500 ###
Nspec =   8 ###

# format z, total, O, H1, H2, Na, Cl, C1, C2, P: header of the file density.out
# You can see that the number of columns of that file is Nspec+2.

############################################################
# Conversions

bohr		= 0.529177 		# in angstrom
# from ion/bohr**3 to ion/A**3
conv		= 1/bohr**3

############################################################

densities = numpy.zeros([Nbins,Nspec+2], numpy.float64) #creates a matrix (of size Nbins*(Nspec+2)) of zeros, datatype float64
Ntotsteps = 0

############################################################

startdir = os.getcwd() #gets the current working directory
execpwd=os.getcwd()
#print "running from", execpwd

for irun in range(runmin,runmax+1):

	# move to directory and uncompress files

	dirrun = '../run%03d' %( irun )
	os.chdir(dirrun)	#changes directory
	cmd = 'tar -xzf otherouts.tar.gz'
	#print cmd
	os.system(cmd)	#executes a command as it was typed in shell

	# open density.out

	infile = open( 'density.out', 'r' )
	lines = infile.readlines() #returns a list containing the lines
	infile.close()

	# read number of steps for this file

	tokens = lines[2].split() #read the file density.out: the 3rd line (index 2) contains the number of steps
	Nsteps=int(tokens[3]) #the number of steps in the 3rd data in the 2nd line
	#print Nsteps, "steps"
	Ntotsteps += Nsteps

	# read densities and accumulate, remember that densities is a file which contains a coarse grained description.
	
	for i in range(Nbins):
		for j in range(Nspec+2):
			tokens = lines[i+4].split()
			densities[i,j] += float(tokens[j])*Nsteps #multiply by Nstep, it's like you were reading the same density at each step.

	# clean directory and come back

	cmd = 'rm *.out'
	#print cmd
	os.system(cmd)	
	os.chdir(execpwd)	


# normalize as appropriate

densities /= Ntotsteps #you are doing an average of the (already 1-run averaged) densities

# convert first column to angstrom

densities[:,0] *= bohr

# convert all other columns to ion/A**3
for j in range(1,Nspec+2):
	densities[:,j] *= conv

#rewrite the header and construct a file containing the average over all the runs.
text  = '#'+startdir+'\n'
text   += "# Average over %d steps\n" %( Ntotsteps )
text   += "# z (A)\t total\t O\t     H1\t    H2\t    Na\t    Cl\t    C1\t    C2\t    P\n"

for i in range(Nbins):
	for j in range(Nspec+2):
		value =	'%.6f\t' %( densities[i,j] )
		text += value
	text += "\n" 

os.chdir(startdir)
cmd='mkdir -p data_OUT'
os.system(cmd)
tmpf=open("data_OUT/c_density_AVG_%02d_%02d.dat"%(runmin,runmax),'w')
tmpf.write(text)
tmpf.close()

