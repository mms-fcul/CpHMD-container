mkdir -p tmp/

cat ../*.occ > all.occ
rm -f *.dat *.tit

## Get sites file ##

pdb=`awk '/PDBin/ {print substr($2,8,length($2)-8)}' ../CpHMD.settings`

/home/noliveira/CpHMD/CpHMD_v3.0/scripts/make_sites ../CpHMD.settings $pdb > base.sites

i=1
cat base.sites | while read line
do
    k=`printf "%02d\n" $i`
    resn=`echo $line | awk '{print $2}'`

    #### Do hetero atom ####
    if [[ $resn == "DON"* ]] || [[ $resn == "GAL"* ]]
    then	
	awk -v j=$i '{print ($j > 0 ? 1 : 0)}' all.occ > tmp/occ_$k.dat
    fi
    
    #### Do amines ####
    if [[ $resn == "NTR"* ]] || [[ $resn == "LYS"* ]];
    then
	awk -v j=$i '{print ($j > 2 ? 1 : 0)}' all.occ > tmp/occ_$k.dat
    fi

    #### Do acids ####
    if [[ $resn == "ASP"* ]] || [[ $resn == "GLU"* ]] || [[ $resn == "CTR"* ]];
    then
	awk -v j=$i '{print ($j > 3 ? 0 : 1)}' all.occ > tmp/occ_$k.dat
    fi
    
    #### Do ARG  ####
    if [[ $resn == "ARG"* ]] ;
    then
	awk -v j=$i '{print ($j > 3 ? 1 : 0)}' all.occ > tmp/occ_$k.dat
    fi
    
    #### Do TYR  ####

    if [[ $resn == "TYR"* ]] ;
    then	    
	awk -v j=$i '{print ($j > 1 ? 0 : 1)}' all.occ > tmp/occ_$k.dat
    fi
    
    #### Do CYS  ####
    if [[ $resn == "CYS"* ]] ;
    then
	awk -v j=$i '{print ($j > 2 ? 0 : 1)}' all.occ > tmp/occ_$k.dat
    fi
    
    #### Do HIS  ####
    if [[ $resn == "HIS"* ]] ;
    then
	awk -v j=$i '{print ($j > 1 ? 1 : 0)}' all.occ > tmp/occ_$k.dat
    fi
    
    echo $resn $i
    i=$((i+1))
done


rm -f *~;rm -f *#
