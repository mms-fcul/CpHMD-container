#!/bin/bash -e 

grom=/gromacs/gromacs-5.1.5_pH_I/bin/gmx

#concatenate using trjcat

mkdir -p Control
mkdir -p fixgro
rm -rdf ./GROs/

ndx="../../01_boxmin/index.ndx"
tpr="../../01_boxmin/min1.tpr"
fix="../../AChE.fixgro"

for s in {001..200}
do
    if [ ! -f Control/$s.dat -a -f ../*_$s.xtc ]
    then
        mkdir -p GROs/	
		
	$grom trjconv -f ../*_$s.xtc \
		      -o ./GROs/frame.gro \
		      -n $ndx \
		      -s $tpr -sep  \
                      -center -pbc mol <<EOF
Solute
Solute
EOF
	
        for i in {0..1000}
        do
            if  [ -f ./GROs/frame${i}.gro ]
            then
                j=`printf "%05d\n" $i`
  
                /programs/fixbox-1.2/fixbox ./GROs/frame$i.gro ../../../../AChE.fixgro > ./GROs/corrected-frame$j.gro
            fi
        done

	cat ./GROs/corrected-*.gro > aux_cor.gro
	
	$grom trjconv -f ./aux_cor.gro  \
                      -o ./traj${s}.xtc \
                      -n $ndx <<EOF
Solute
EOF

        echo "control file" > Control/$s.dat
        rm -rdf ./GROs
    fi
done


$grom trjcat -f ./traj???.xtc \
             -o ./traj.xtc \
	     -n $ndx <<EOF
Solute
EOF
    
#$grom trjconv -f ./traj.xtc  \
#      -o ./tmp.gro \
#      -n $ndx \
#      -s $tpr -skip 10 <<EOF
#Solute
#EOF
        
rm -f ./aux* ./\#* ./*~
