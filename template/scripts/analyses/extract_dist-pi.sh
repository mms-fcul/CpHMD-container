#!/bin/bash -e

grom=/gromacs/gromacs-5.1.5/bin/gmx
ndx="../../01_boxmin/index.ndx"
tpr="../../01_boxmin/min1.tpr"


## Reduce the number of points in the xtc ##
$grom trjconv -f ../traj/traj.xtc \
      -s $tpr \
      -o tmp.xtc \
      -n $ndx \
      -skip 10 <<EOF
Solute
EOF
    
$grom mindist -f tmp.xtc \
	  -s $tpr \
	  -n $ndx \
	  -od ./dist-pi.xvg \
	  -pi <<EOF
AChE
EOF

rm -rdf tmp.xtc
    




    
