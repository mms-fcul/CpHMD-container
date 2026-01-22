#!/bin/bash

#this script is used to change the aminoacids.rtp file so that the
#  charge of each block matches that of the average for the residues
#  in the format NXXX and CXXX, which are in the protstates.dic

Dir=`pwd`

#for every NXXX residue other than NGLY and NPRO
for i in `awk '!/NGLY|NPRO/ && /\[ N[A-Z][A-Z][A-Z] \]/,/\[ bonds \]/ { if ($0 ~ /\[ N/) { print $2 }}' ${Dir}/../aminoacids.rtp`
do

    #first we get the column with the average charges for the main
    #  atoms in the main chain
    awk '/NT /{printf "%2s%8s0\n", $3, $NF}' ${Dir}/../protstates.dic > tmp_avg-charge.dat

    # awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O|CB|HB.*)$/) printf "%6s%6s%18s%6s  ; %8s\n", $1, $2, $3, $4, $3; else printf "%6s%6s%18s%6s\n", $1, $2, $3, $4}}' ${Dir}/../../amber14sb_OL15.ff/aminoacids.rtp > tmp_${i}.dat

    awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O|CB|HB.*)$/) printf "%6s%6s%18s%6s  ; %8s\n", $1, $2, $3, $4, $3; else printf "%6s%6s%18s%6s\n", $1, $2, $3, $4}}' ${Dir}/../aminoacids.rtp > tmp_${i}.dat


    awk 'NR==FNR{a[$1]=$2;next} {if ($1 in a) printf "%6s%6s%18s%6s%3s%9s\n", $1, $2, a[$1], $4, $5, $6; else print}' tmp_avg-charge.dat tmp_${i}.dat > tmp-changed-charges.dat

    #now we have to compensate the new main chain charges in the side
    #  chain so that we have a correct total charge in each part of the
    #  the residues. The charge readjustment will be spredout throughout
    #  the CB and the H bond to it.
    ch_sum=`awk '{ if ($1 !~ /^(N|H1|H2|H3|CA|HA|C|O)$/) { sum +=$3 } } END {print sum}' tmp-changed-charges.dat`
    n_cb=`awk '{ if ($1 ~ /^(CB|HB.*)$/) {count++;} } END {print count;}' tmp-changed-charges.dat`
    closest_int=`echo "${ch_sum}" | awk '{printf("%d\n", $1 < 0 ? int($1-0.5) : int($1+0.5))}'`
    difference=`echo "scale=4; ${closest_int} - ${ch_sum}" | bc`

    avg_ch=`echo "scale=4; ${difference} / ${n_cb}" | bc`
    left_over=`echo "${difference} - ${avg_ch} * ${n_cb}" | bc`

    awk -v n=${avg_ch} -v left=${left_over} '{ 
    if ($1 ~ /^(CB)$/) printf("%6s%6s%18.5f%6s%3s%9s%9.5f\n", $1, $2, $3+n+left, $4, $5, $6, n+left); 
    else if ($1 ~ /^(HB.*)$/) printf("%6s%6s%18.5f%6s%3s%9s%9.5f\n", $1, $2, $3+n, $4, $5, $6, n);
    else  printf("%6s%6s%18.5f%6s%3s%9s\n", $1, $2, $3, $4, $5, $6)
}' tmp-changed-charges.dat > tmp-final-charges.dat

    #now we just have to get the changed block to where it bellongs in the
    #  aminoacids.rtp file
    #we first need to add a " [ atoms ]" line to the top of tmp-final-charges.dat
    sed -i '1s/^/ /; 1s/^ / [ atoms ]\n/' tmp-final-charges.dat
    sed -e "/\[ ${i} \]/,/\[ bonds \]/ { /\[ ${i} \]/ {p; r tmp-final-charges.dat
   }; /\[ bonds \]/ {p; d}; d }" -i ${Dir}/../aminoacids.rtp
    
done

#Note that for NGLY and NPRO the procedure to be done is different
#  than just using the average charges of the others as the atoms
#  that make up the main chain are different


#
#
### for the C-ter residues:
#
#
for i in `awk '!/CGLY|CPRO/ && /\[ C[A-Z][A-Z][A-Z] \]/,/\[ bonds \]/ { if ($0 ~ /\[ C/) { print $2 }}' ${Dir}/../aminoacids.rtp`
do

    #first we get the column with the average charges for the main
    #  atoms in the main chain
    awk '/^CT /{printf "%2s%8s0\n", $3, $NF}' ${Dir}/../protstates.dic > tmp_avg-charge.dat

    # awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H1|H2|H3|CA|HA|C|O|CB|HB.*)$/) printf "%6s%6s%18s%6s  ; %8s\n", $1, $2, $3, $4, $3; else printf "%6s%6s%18s%6s\n", $1, $2, $3, $4}}' ${Dir}/../../amber14sb_OL15.ff/aminoacids.rtp > tmp_${i}.dat

    awk '/\[ '${i}' \]/,/\[ bonds \]/{ if ($0 !~ /\[/) { if ($1 ~ /^(N|H|CA|HA|C|OC1|OC2|CB|HB.*)$/) printf "%6s%6s%18s%6s  ; %8s\n", $1, $2, $3, $4, $3; else printf "%6s%6s%18s%6s\n", $1, $2, $3, $4}}' ${Dir}/../aminoacids.rtp > tmp_${i}.dat

    
    awk 'NR==FNR{a[$1]=$2;next} {if ($1 in a) printf "%6s%6s%18s%6s%3s%9s\n", $1, $2, a[$1], $4, $5, $6; else print}' tmp_avg-charge.dat tmp_${i}.dat > tmp-changed-charges.dat

    #now we have to compensate the new main chain charges in the side
    #  chain so that we have a correct total charge in each part of the
    #  the residues. The charge readjustment will be spredout throughout
    #  the CB and the H bond to it.
    ch_sum=`awk '{ if ($1 !~ /^(N|H|CA|HA|C|OC1|OC2)$/) { sum +=$3 } } END {print sum}' tmp-changed-charges.dat`
    n_cb=`awk '{ if ($1 ~ /^(CB|HB.*)$/) {count++;} } END {print count;}' tmp-changed-charges.dat`
    closest_int=`echo "${ch_sum}" | awk '{printf("%d\n", $1 < 0 ? int($1-0.5) : int($1+0.5))}'`
    difference=`echo "scale=4; ${closest_int} - ${ch_sum}" | bc`

    avg_ch=`echo "scale=4; ${difference} / ${n_cb}" | bc`
    left_over=`echo "${difference} - ${avg_ch} * ${n_cb}" | bc`

    awk -v n=${avg_ch} -v left=${left_over} '{ 
    if ($1 ~ /^(CB)$/) printf("%6s%6s%18.5f%6s%3s%9s%9.5f\n", $1, $2, $3+n+left, $4, $5, $6, n+left); 
    else if ($1 ~ /^(HB.*)$/) printf("%6s%6s%18.5f%6s%3s%9s%9.5f\n", $1, $2, $3+n, $4, $5, $6, n);
    else  printf("%6s%6s%18.5f%6s%3s%9s\n", $1, $2, $3, $4, $5, $6)
}' tmp-changed-charges.dat > tmp-final-charges.dat

    #now we just have to get the changed block to where it bellongs in the
    #  aminoacids.rtp file
    #we first need to add a " [ atoms ]" line to the top of tmp-final-charges.dat
    sed -i '1s/^/ /; 1s/^ / [ atoms ]\n/' tmp-final-charges.dat
    #to add the comment noting the atom type change in the C
    awk -F ',' '$1 ~ / C / {print $0", atom type changed C -> CO"}; $1 !~ / C / {print $0}' tmp-final-charges.dat > tmp-final.dat
    #change in aminacids.rtp
    sed -e "/\[ ${i} \]/,/\[ bonds \]/ { /\[ ${i} \]/ {p; r tmp-final.dat
   }; /\[ bonds \]/ {p; d}; d }" -i ${Dir}/../aminoacids.rtp

    rm -f tmp*.dat
done
