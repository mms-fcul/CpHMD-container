#!/bin/bash

### Help with usage function ###
Help()
{
    # Display Help
    echo "This function adds a new molecule to the delphi databases (crg and siz). To run this script you should use:"
    echo "singularity run --app add-mol-database ||<container name>|| -f <force field ffnonbonded.itp> -i <molecule.itp file> -r <moleculetype name within itp>"
    echo
    echo "options:"
    echo "h    Print this Help."
    echo "f    force field name (GROMOS,CHARMM,AMBER) or custom location of ffnonbonded.itp."
    echo "i    itp file for the molecule to add (can be obtained by CHARMM-GUI for example"
    echo "r    Names of the moleculetypes as they appear in the .itp file. Can be more than one. (needs testing)"
    echo "e    Energy for Lennard-Jonnes radii calculation, zero, min, RT, 2RT,5RT (default 2RT) "
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


###### Defining LJ energy to use the radii calculation ######

LJenergy="2RT"


while getopts "hf:i:e:r:" opt ; do
    case $opt in
        i) file=${OPTARG} ;
	   if [[ `echo $file | awk -F "." '{print $NF}'` != "itp" ]]
	   then
	       echo "File given does not end with .itp, please make sure to give a itp file" >&2
	       exit 1
	   fi ;;
	f) ff=${OPTARG} ;
	   case $ff in
	       GROMOS|Gromos|gromos)\
		   ffloc="/FFs/G54a7pH.ff/ffnonbonded.itp" ;
		   dbloc="DataBaseT_G54a7pH" ;
		   ;;
	       CHARMM|Charmm|charmm)\
		   ffloc="/FFs/CHARMM36pH.ff/ffnonbonded.itp" ;
		   dbloc="DataBaseT_CHARMM36pH" ;
		   ;;
	       AMBER|Amber|amber)\
		   ffloc="/FFs/Amber14SBpH.ff/ffnonbonded.itp" ;
		   dbloc="DataBaseT_Amber14SBpH" ;
		   ;;
	       *)\
		   ffloc=${OPTARG}
		   if [[  `echo $ffloc | awk -F "." '{print $NF}'` != "itp" ]] 
		   then
		       echo "ffnonbonded.itp given is not an itp or is not found, please fix your path."
		       exit 1 ; 
		   fi ;;
	       "") echo "Missing or invalid force field to add the database. Please give one of the following GROMOS, CHARMM, AMBER or your custom force-field path. " ;;
	   esac ;;
	e) LJenergy=${OPTARG}
	   ;;
	r) getopts-extra "$@" ;
           res=("${OPTARG[@]}");
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

### Gather all residues in the given file ###

for rr in "${!res[@]}" ; do
    echo "Starting to build the database files for the molecule ${res[$rr]}"

    ### Creating the new .crg database to add ###
    awk -v d=${res[$rr]} '$4==d && !/^;/ {printf"%-6s%-9s%6.3f\n", $5,$4,$7}' $file >> ./charges-to-add-on-database.crg

    ### Creating the siz radii files to add on the databases ###
    ## First create the radii file from ffnonbonded.itp ##
    
    if [[ ! -z `awk '/nonbond_params/' $ffloc` ]] ; then
	echo "treating ffnonbonded as GROMOS like, identifying \[nonbond_params\] delimiter "
	awk 'NR==1,/nonbond_params/ && !/^;/' $ffloc | grep OW > params.itp	
	awk '/nonbond_params/,/\[ pair/ ' $ffloc | awk 'NF>3 && !/^;/' | grep OW | sed 's/OW//' >>params.itp 

	/analysis/GROMOS_LJ-to-radii ./params.itp ${LJenergy} > LJ_radii.dat
    elif [[ ! -z `awk '/pairtypes/' $ffloc` ]] ; then
	echo "treating ffnonbonded as CHARMM like, identifying \[pairtype\] delimiter"
	awk 'NR==1,/pairtypes/ && !/^;/' $ffloc | sed 's/\[ pairtypes \]//;/^$/d' > params.itp
	/analysis/Amber14SBpH_LJ-to-radii ./params.itp ${LJenergy} > LJ_radii.dat
	
    else
	echo "treating ffnonbonded as AMBER like"
	awk '!/^;/' $ffloc | sed '/^$/d' > params.itp
	/analysis/Amber14SBpH_LJ-to-radii ./params.itp ${LJenergy} > LJ_radii.dat
    fi	

    ### Finally reading every atom of the molecule given as .itp and writting its area for a database ###
    
    while read -r line
    do
	atomtype=`echo $line | awk '{print $2}' `
	if [ -z `awk -v i=$atomtype '$1==i  {print $1}' LJ_radii.dat ` ]
	then
	    echo "Atomtype $atomtype does not exist in ffnonbonded.itp"
	    exit 1
	fi
	atomtype=`echo $line | awk '{print $2}' `
	atom_radii=`awk -v i=$atomtype '$1==i {print $2}' LJ_radii.dat`
	atomname=`echo $line | awk '{print $5}' `
	echo $atomname ${res[$rr]} ${atom_radii} | awk '{printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' >> ./radii-to-add-on-database.siz 
    done < <(awk -v d=${res[$rr]} '$4==d && !/^;/ {print $0}' $file )
    
done

### Make the final databases to the output folder ###

#for db in siz crg
#do
#    case $db in
#	siz)     file="./radii-to-add-on-database.siz" ;;
#	siz)     file="./charges-to-add-on-database.crg" ;;
#    esac
#
#    ## check to make sure database ends with new line ##
#    if [ -z "$(tail -c 1 < /Databases/${dbloc}.${db})" ];
#    then
#	cat /Databases/${dbloc}.${db} $file  > ./Databases/${dbloc}.${db}
#    else
#	cat /Databases/${dbloc}.${db} ; echo "" ; cat $file > ./Databases/${dbloc}.${db}
#    fi
#done


rm -rf params.itp LJ_radii.dat #./radii-to-add-on-database.siz ./charges-to-add-on-database.crg
