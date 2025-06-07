#!/bin/bash -e

### Help with usage function ###
Help()
{
    # Display Help
    echo "This function converts the regular PDB naming into CpHMD compatible naming."
    echo "Changes are made to the residue names for them to be in line with the CpHMD tautomer naming."
    echo
    echo "Syntax: singularity run --app pdb2cphmd CpHMD.sif [-h] -p <pdb file> -f <force field> -r <residue names to change for CpHMD>"
    echo "Example: singularity run --app pdb2cphmd CpHMD.sif -p start.pdb -f G54a7pH -r asp glu his "
    echo
    echo "options:"
    echo "h    Print this Help."
    echo "p    Base pdb that will be changed."
    echo "f    Choose the force-field to use, available names G54a7pH, CHARMM36pH, Amber14SBpH. "
    echo "r    Residue names to change "
    echo
}
### Reading number of variables ###
function getopts-extra () {
    declare i=1
    # if the next argument is not an option, then append it to array OPTARG
    while [[ ${OPTIND} -le $# && ${!OPTIND:0:1} != '-' ]]; do
        OPTARG[i]=${!OPTIND}
        let i++ OPTIND++
    done
}


while getopts "hp:f:r:" opt ; do
    case $opt in
        p) file=${OPTARG} ;
	   if [[ `echo $file | awk -F "." '{print $NF}'` != "pdb" ]]
	   then
	       echo "File given does not end with .pdb, please make sure to give a pdb file" >&2
	       exit 1
	   fi ;;
	f) fg=${OPTARG} ;
	   case $fg in
	       GROMOS|Gromos|gromos)\
		   ff="G54a7pH" ;;
	       CHARMM|Charmm|charmm)\
		   ff="CHARMM36pH" ;;
	       AMBER|Amber|amber)\
		   ff="Amber14SBpH" ;;
	       *|"")\
		   echo "Missing or invalid force field to prepare your pdb. Please give one of the following GROMOS, CHARMM, AMBER. " ;;
	   esac ;;
	r) getopts-extra "$@" ;
           res=( "${OPTARG[@]}" );
	   ;;
	h) Help ;
	   exit 1 ;;
	*)
            echo 'Error in command line parsing' >&2
	    Help ;
            exit 0 ;;
    esac
done
shift $(( OPTIND - 1 ))

### Copy the required FF to the working folder! ###
echo "Copying files required for pdb2gmx treatment. Force-field chosen ${ff}." 

cp -rf /FFs/${ff}.ff ./
cp -rf /FFs/residuetypes.dat ./
cp -rf /FFs/specbond.dat ./

### Get new name right ###

fname=`echo $file |  awk -F "." '{NF--;print}' | sed 's/ /./g' `

## make new CpHMD file ##

cp -rf $file ./${fname}_CpHMD.pdb

for rr in "${!res[@]}"
do
    echo "Changing residue "${res[$rr]}
    case ${res[$rr]} in
	
	asp|Asp|ASP)
	    sed -i 's/ASP/AS4/; s/ASPH/AS0 /' ./${fname}_CpHMD.pdb ;;
	glu|Glu|GLU)
	    sed -i 's/GLU/GL4/; s/GLUH/GL0 /' ./${fname}_CpHMD.pdb ;;
        tyr|Tyr|TYR)
	    sed -i 's/TYR/TY0/' ./${fname}_CpHMD.pdb ;;
	lys|Lys|LYS)
	    sed -i 's/LYS/LY3/' ./${fname}_CpHMD.pdb ;;
	cys|Cys|CYS)
	    echo
	    echo -e "\033[38;2;255;0;02m -- Warning -- \033[m"
	    echo -e "Cys is a sensitive residue to treat. If you are sure that no cysteine in your protein is oxidized then use \033[38;2;255;0;02m cys-force\033[m instead of cys to convert cysteines into CpHMD nomenclature."
	    echo ;;
	cys-force|Cys-force|CYS-force)
	    sed -i 's/CYS/CY0/' ./${fname}_CpHMD.pdb ;;
	his|His|HIS)
	    sed -i 's/HIS/HI2/; s/HSD/HI0/ ; s/HSE/HI1/ ; s/HSP/HI2/ ; s/HID/HI0/ ; s/HIE/HI1/ ; s/HIP/HI2/' ./${fname}_CpHMD.pdb ;;
	*)
	    echo 'Error: Invalid residue to change' >&2
            exit 1 ;;
    esac

    
    #echo $i
done
