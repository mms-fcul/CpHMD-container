#!/bin/bash

#this script is used to calculate the average charge for all the NXXX and CXXX
#  residues in the Amber14SB FF, so that they can be used in the protstates.dic
#  for the PypKa MD program for the titrating termini (this will give us the
#  value of the charged termini)

Dir=`pwd` #/home/jsequeira/Documents/TMP/scripts_to_edit/Amber14SBpH.ff/scripts

#N-ter
#get just the lines we want to calculate the averages of
awk '/\[ N[A-Z][A-Z][A-W,Y-Z] \]/,/\[ bonds \]/ { if ($1 ~ /^(N|H1|H2|H3|CA|HA|HA1|HA2)$/) { if ($0 !~ /\[/) print } }' ${Dir}/../aminoacids.rtp > ${Dir}/tmp-nt.dat

#calculate the average for each of the
for atom in N H1 H2 H3 CA HA
do
    if [[ $atom != HA ]]; then
	echo "${atom} " `awk -v atom=${atom} '$1 == atom {sum += $3; count++} END {if (count > 0) print sum/count}' ${Dir}/tmp-nt.dat`
    else
	echo "${atom} " `awk -v atom=${atom} '$1 == atom || $1 == "HA1" || $1 == "HA2" {sum += $3; count++} END {if (count > 0) print sum/count}' ${Dir}/tmp-nt.dat`
    fi
done > nt-charges.tmp



#C-ter

#get just the lines we want to calculate the averages of
awk '/\[ C[A-Z][A-Z][A-W,Y-Z] \]/,/\[ bonds \]/ { if ($1 ~ /^(C|OC1|HC11|HC12|OC2|HC21|HC22)$/) { if ($0 !~ /\[/) print } }' ${Dir}/../aminoacids.rtp > ${Dir}/tmp-ct.dat

#calculate the average for each of the 
for atom in C OC1 HC11 HC12 OC2 HC21 HC22
do
    echo "${atom} " `awk -v atom=${atom} '$1 == atom {sum += $3; count++} END {if (count > 0) print sum/count}' ${Dir}/tmp-ct.dat`
done > ct-charges.tmp
