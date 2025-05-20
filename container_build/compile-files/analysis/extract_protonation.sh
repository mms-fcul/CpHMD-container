#!/bin/bash -e

mkdir -p protonation_analysis/

cat $1 > protonation_analysis/all.occ

i=1
cat $2 | while read line
do
    k=`printf "%02d\n" $i`
    resn=`echo $line | awk '{print $2}'`
    rnumb=`echo $line | awk '{print $1}'`
    name=`echo $line | awk '{print substr($2,0,3)}'`
    
    #### Do amines - tautomer above 2 means the residue is protonated ####
    if [[ $resn == "NTR"* ]] || [[ $resn == "LYS"* ]];
    then
	awk -v j=$i '{print ($j > 2 ? 1 : 0)}' protonation_analysis/all.occ > protonation_analysis/protonation_${rnumb}${name}.dat
    fi

    #### Do acids - tautomer above 3 means the residue is protonated  ####
    if [[ $resn == "ASP"* ]] || [[ $resn == "GLU"* ]] || [[ $resn == "CTR"* ]];
    then
	awk -v j=$i '{print ($j > 3 ? 0 : 1)}' protonation_analysis/all.occ > protonation_analysis/protonation_${rnumb}${name}.dat
    fi
    
    #### Do ARG  - tautomer above 3 means the residue is protonated ####
    if [[ $resn == "ARG"* ]] ;
    then
	awk -v j=$i '{print ($j > 3 ? 1 : 0)}' protonation_analysis/all.occ > protonation_analysis/protonation_${rnumb}${name}.dat
    fi
    
    #### Do TYR  - tautomer above 1 means the residue is protonated ####

    if [[ $resn == "TYR"* ]] ;
    then	    
	awk -v j=$i '{print ($j > 1 ? 0 : 1)}' protonation_analysis/all.occ > protonation_analysis/protonation_${rnumb}${name}.dat
    fi
    
    #### Do CYS  - tautomer above 2 means the residue is protonated ####
    if [[ $resn == "CYS"* ]] ;
    then
	awk -v j=$i '{print ($j > 2 ? 0 : 1)}' protonation_analysis/all.occ > protonation_analysis/protonation_${rnumb}${name}.dat
    fi
    
    #### Do HIS  - tautomer above 2 means the residue is protonated ####
    if [[ $resn == "HIS"* ]] ;
    then
	awk -v j=$i '{print ($j > 1 ? 1 : 0)}' protonation_analysis/all.occ > protonation_analysis/protonation_${rnumb}${name}.dat
    fi
    
    echo $resn $i
    i=$((i+1))
done

rm -rdf protonation_analysis/*~;rm -rdf protonation_analysis/*#
