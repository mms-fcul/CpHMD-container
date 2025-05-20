## Definition of variables
grom="/gromacs/gromacs-2018.6/bin/gmx"

# Variables to use the system information (.tpr), index and full trajectory of the simulation
tpr="../../init/i100.tpr"
ndx="../../boxmin/index.ndx"
traj="../traj/traj.xtc"

startres=1  # residue to start measuring the dssp partial
endres=542  #residue to end measuring dssp partial


ln -sf /home/noliveira/CpHMD/CpHMD_v3.0/template/scripts/Partial_SS .

#
# The do_dssp module provides the secondary structure of a given trajectory or structure file provided in the -f flag
# Two distinct output files are provided, an assignment of secondary structure per residue per frame (in case of a trajectory was given as input) in a matrix file (.xpm) and a count of total number of residues with each secondary structure type per frame (.xvg)
# In this case, the calculations are performed for each monomer of the protein and for the protein in its dimeric form, as 

Dir=`pwd`

pwd > ~/a.txt

ssh bio000 '
Dir=`awk '{print}' a.txt`
grom="/gromacs/gromacs-2018.6/bin/gmx"
cd $Dir

tpr="../../01_boxmin/min1.tpr"
ndx="../../01_boxmin/index.ndx"
traj="../traj/traj.xtc"


$grom do_dssp -f $traj -s $tpr -n $ndx \
      -o dssp_AChE.xpm \
      -sc ss_AChE.xvg -ver 1 <<EOF
AChE	      	      	
EOF

rm -f dd*
rm -f /home/noliveira/a.txt
'

duration=`ls ../*_???.xtc | wc -l | awk '{print $1*100}'`
#duration=1000
./Partial_SS dssp_AChE.xpm ${startres} ${endres} $duration ./anal_dssp_AChE.xvg

awk '{print $1/100, ($3/$9)*100}' ./anal_dssp_AChE.xvg > Beta_anal.xvg
awk '{print $1/100, ($2/$9)*100}' ./anal_dssp_AChE.xvg > Helix_anal.xvg


rm -f *~ *# aux* 
