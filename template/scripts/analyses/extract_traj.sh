#!/bin/bash -e 

grom=/gromacs/gromacs-5.1.5_pH_I/bin/gmx

#concatenate using trjcat

mkdir -p Control

ndx="../../01_boxmin/index.ndx"
tpr="../../01_boxmin/min1.tpr"

for s in {001..200}
do
    if [ ! -f Control/$s.dat -a -f ../*_$s.xtc ]
    then
	$grom trjconv -f ../*_${s}.xtc \
	      -o ./traj_${s}.xtc \
	      -n $ndx \
	      -s $tpr \
	      -center -pbc mol  <<EOF
Solute
Solute
EOF

	echo "control file" > Control/$s.dat
    fi
    
done

$grom trjcat -f ./traj_*.xtc \
      -o ./traj.xtc \
      -n $ndx <<EOF
Solute
EOF


$grom trjconv -f ./traj.xtc  \
      -o ./tmp.gro \
      -n $ndx \
      -s $tpr -skip 10 <<EOF
Solute
EOF
        
rm -f ./aux* ./\#* ./*~
