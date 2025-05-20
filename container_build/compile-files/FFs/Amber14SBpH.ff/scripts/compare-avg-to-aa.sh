#!/bin/bash

#this script is used to obtain the differences from each aa block
#  to the average charge that was considered for each.

Dir=`pwd`

#for every NXXX residue other than NGLY and NPRO
for i in `awk '!/NGLY|NPRO/ && /\[ N[A-Z][A-Z][A-Z] \]/,/\[ bonds \]/ { if ($0 ~ /\[ N/) { print $2 }}' ${Dir}/../aminoacids.rtp`
do

    #first we get the information from the blocks that already have
    #  the average charge in the main chain
    awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O)$/) printf "%6s%6s%18s%6s%3s\n", $1, $2, $3, $4, $5}}' ${Dir}/../aminoacids.rtp > tmp_${i}.dat

    #and we get the information from the original blocks
    awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O)$/) printf "%9s\n", $3}}' ${Dir}/../../amber14sb_OL15.ff/aminoacids.rtp > tmp_col.dat

    #put it all together in a file
    paste tmp_${i}.dat tmp_col.dat > tmp.dat

    #and now just calculate the difference between the average
    #  and the old charges.
    awk -v amino=${i} '{diff+=$3-$6}END{print amino "  " diff}' tmp.dat >> delta_to_avg.dat

done

# CXXX residues
for i in `awk '!/CGLY|CPRO/ && /\[ C[A-Z][A-Z][A-Z] \]/,/\[ bonds \]/ { if ($0 ~ /\[ C/) { print $2 }}' ${Dir}/../aminoacids.rtp`
do

    #first we get the information from the blocks that already have
    #  the average charge in the main chain
    awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O)$/) printf "%6s%6s%18s%6s%3s\n", $1, $2, $3, $4, $5}}' ${Dir}/../aminoacids.rtp > tmp_${i}.dat

    #and we get the information from the original blocks
    awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O)$/) printf "%9s\n", $3}}' ${Dir}/../../amber14sb_OL15.ff/aminoacids.rtp > tmp_col.dat

    #put it all together in a file
    paste tmp_${i}.dat tmp_col.dat > tmp.dat

    #and now just calculate the difference between the average
    #  and the old charges.
    awk -v amino=${i} '{diff+=$3-$6}END{print amino "  " diff}' tmp.dat >> delta_to_avg.dat

done

#to finish compute the average distance to the the average charge
#  over all affected blocks
average=`awk '{sum+=$2; n++} END {print sum/n}' delta_to_avg.dat`
echo "Average diff to avg charges: ${average}" >> delta_to_avg.dat


rm -f tmp*
