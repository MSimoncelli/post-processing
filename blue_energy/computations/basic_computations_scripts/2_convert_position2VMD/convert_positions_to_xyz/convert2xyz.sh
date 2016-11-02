#!/bin/bash

if [ $# -ne 1 ]; then
	echo "Usage: file"
	exit 1;
fi

file=$1

echo "counting frames in files" ${file}

natms=$( awk 'BEGIN{n=0;start=0;}{ if ($1=="Positions") { start=1; } else if (start==1 && (($1!="Cell")&&($1!="Velocities"))) {n++;} else if (($1=="Cell")||($1=="Velocities")) { start=0; } }END{print n}' ../run000/restart.dat )

echo "found" ${natms} "atoms in ../run000/restart.dat"

nlines=$( wc ${file} | awk '{print $1}' )

echo "found" ${nlines} "lines in" ${file}

nframes=$( awk -v nl=${nlines} -v na=${natms} 'END{print nl/na}' /dev/null )

echo "number of frames in" ${file} ":" ${nframes} 

sed -e 's/nframes/'${nframes}'/g' pm2xyz-inpt.template > pm2xyz-inpt

./pm2xyz.x < pm2xyz-inpt
