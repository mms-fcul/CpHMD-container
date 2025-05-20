#!/bin/bash -e 
if [[ $# != 1 ]] ; then
  echo 'Script needs an argument that is the starting pdb'
  exit 0
fi

grom=/gromacs/gromacs-2021.2/bin/gmx

Sys=XXX

sed -i 's/ASP/AS4/g;s/GLU/GL4/g;s/HIS/HI0/g;s/CYS/CY0/g;s/TYR/TY0/g;s/LYS/LY3/g' $1

$grom pdb2gmx -f HEWL-start.pdb -p ${Sys}.top -o ${Sys}.gro -i ${Sys}_posre.itp -ter -ignh -renum -water spc -merge all <<EOF
1
3
4
EOF


sed -i 's/Protein_chain_A/Protein/' ${Sys}.top

awk '{\
if($0=="\[ moleculetype \]"){print ";\[ begin_tit_molecule \]"}\
if($0=="; Include water topology"){print ";\[ end_tit_molecule \]"}\
print $0}' ${sys}.top > aux
mv aux ${sys}.top

##### Make new box with the correct size ####

$grom editconf -f ${Sys}.gro -o ${Sys}_box.gro -c -d 0.9 -bt dodecahedron

$grom solvate -cp ${Sys}_box.gro  -cs spc216.gro -o ${Sys}_solv.gro -p ${Sys}.top

##### Neutralize Box ####
echo "" > ${Sys}_ion.mdp 

$grom grompp -f ${Sys}_ion.mdp -c ${Sys}_solv.gro -p ${Sys}.top -o ${Sys}_ion.tpr -maxwarn 100

${grom} genion -s ${Sys}_ion.tpr -p ${Sys}.top -o ${Sys}_ion.gro -neutral <<EOF
SOL
EOF

rm -f ./\#* ./*~ ./aux*
