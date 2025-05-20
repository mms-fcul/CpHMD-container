## Definition of variables
grom=/gromacs/gromacs-5.1.5/bin/gmx

ndx="../../01_boxmin/index.ndx"
tpr="../../01_boxmin/min1.tpr"

$grom rms -f ../traj/traj.xtc \
      -s $tpr \
      -n $ndx \
      -o rmsd_prot.xvg -fit rot+trans <<EOF
Solute
Solute
EOF


rm -f *~ *# 
