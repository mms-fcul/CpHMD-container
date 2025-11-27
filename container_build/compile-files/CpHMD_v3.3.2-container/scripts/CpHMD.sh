#!/bin/bash  
#
prog=`basename $0` 
usage="Usage: $prog <MyProtein_001>.settings\n
       <MyProtein_001>.settings  : Constant-pH MD settings file"
#
# This allows to use extended pattern matching features
# E.g.: $PDBin == *(" ")
shopt -s extglob 
#
message ()
{
    case "$1" in
        E ) shift; echo -e "$prog: Error: $*" >&2 ; exit 1;;
        U ) shift; echo -e "$prog: Warning: $*\n$usage" >&2; exit 1;;
        W ) shift; echo -e "$prog: Warning: $*" >&2;;
        * ) message E "Wrong use of 'message' function.";;
    esac
}
#
# Parse arguments
if [ $# != 1 ]; then message U "Wrong number of arguments."; fi

#
## Defining the default locations for programs/paths inside container
## These might be overwritten if they are defined in the CpHMD.settings
## which is sourced below.

source /CpHMD/CpHMD-default.settings

# Read Simulation Parameters
source $1
#
# Name of your simulation segment (the output files will
# be generated with this name)
#export blockname=${1%.*}
export runname=${SysName}
export blockname="${SysName}_CpHrun"
#
# Check some parameters
if [[ $ffID != G54a7pH && $ffID != CHARMM36pH && $ffID != Amber14SBpH ]]; then
    message W "\n***** WARNING *****\n**** You are using a force-field which is not standard.\n**** The use of  $ffID implies that the paths for the Delphi databases, tautomer St files, and force field is given to the program."
fi

### Insert here definition of the ST default folder

if [[ -d $StDIR ]] ; then
    message W "StDIR inputed by the user at $StDIR"
else
    StDIR="/STs/St-${ffID}"
    message W "StDIR not found, using base container ST directory."
fi

# Call all functions stored in the "functions" file:
functions=$CpHDIR/scripts/functions.sh
if [[ -f $functions ]]; then
    source $functions            
else
    message E "File $functions is missing. Check CpHDIR parameter in $1."
fi

######################################################### 
# add the treatment of CpHMD.settings to mdp and fixgro #
#########################################################

awk '/#mdp# / {print}' $1 | sed 's/#mdp# //g' > ${SysName}.mdp
awk '/#fixgro# / {print}' $1 | sed 's/#fixgro# //g' > ${SysName}.fixgro

##################################################### 
# add the treatment of plumed input when plumed =1  #
#####################################################
if  [[ -n $plumed && "${plumed}" -eq "1" ]] ; then
    awk '/#plumed# / {print}' $1 | sed 's/#plumed# //g' > ${SysName}_plumed.dat

    case $plumedtype in
	grid)
	    sed -i "s/\&colvar_stride/$colvar_stride/g" ${SysName}_plumed.dat
	    sed -i "s/\&grid_name/$grid_name/g" ${SysName}_plumed.dat
	    sed -i "s/\&hills/$hills/g" ${SysName}_plumed.dat
	    sed -i "s/\&colvar_name/$colvar_name/g" ${SysName}_plumed.dat
	    ;;
	hill|static)
	    sed -i "s/\&colvar_stride/$colvar_stride/g" ${SysName}_plumed.dat
	    sed -i "s/\&hills/$hills/g" ${SysName}_plumed.dat
	    sed -i "s/\&colvar_name/$colvar_name/g" ${SysName}_plumed.dat
	    ;;
    esac
    ############################################################################
    ## Check to add the RESTART flag on plumed which is assumed to be present ##
    ############################################################################
    if [ ! -n "`awk '/RESTART/' ${SysName}_plumed.dat `" ] ;
    then
	awk -i inplace 'NR==1{print "RESTART";print;next} {print}' ${SysName}_plumed.dat
    fi	
fi

#############################################################
# add the correct cicles depending on the size of block asked
#############################################################
dt=`awk '/dt/ {print $3}' ${SysName}.mdp`
if [ -z "$dt" ] ;
then
    dt=0.002
    message W "dt was not defined in the .mdp, falling back to 0.002 default dt."
fi

EndCycle=`echo ${SimTime} $EffectiveSteps $dt | awk '{print ($1*1000)/($2*$3)}'`

message W "Running CpHMD cycles from $InitCycle until $EndCycle with 20 ps each cycle."

#################################################
## PLUMED Detection rebasing ##
#################################################
if [[ $plumed == 1 ]] ;then
    ## change /gromacs file if it is the standard one ##
    if [ $GroDIR == "/gromacs/bin/gmx" ] ; then
	export GroDIR="/gromacs-plumed/bin/gmx"
	message W "PLUMED support requested, changing the gromacs compilation to PLUMED $GroDIR"
    fi
fi

#################################################
## GPU Detection and mdrun parameters rebasing ##
#################################################

if [[ $GPU == 1 ]] ;then
    ## change GPU /gromacs file if it is the standard one ##
    if [ $GroDIR == "/gromacs/bin/gmx" ] ; then
	export GroDIR="/gromacs-gpu/bin/gmx"
	message W "GPU support requested, changing the gromacs compilation to gpu $GroDIR"
    fi
    
   OMP_NUM_THREADS=$nCPU

   ## Edit the mdrun line in the .pHmdp ##
   
   message W "$GroDIR" 
   export mdrun=`echo $GroDIR $mdrungpu `

   message W "changin mdrun to $mdrun"
   
else
    export mdrun=`echo $GroDIR $mdruncpu`
fi

## Creating a running folder to prevent clutter ##
mkdir -p ./CpHMD-run_$$

### Checking if the specific files exist and copying them to the running directory
for f in ${SysName}.mdp $GROin $TOPin $NDXin 
do
    if [ ! -f $f ]; then
        message E  "File $f is missing!!!... Program will crash"
    fi
done

cp $1  ./CpHMD-run_$$ ; cp ${SysName}.mdp ./CpHMD-run_$$ ; cp ${SysName}.fixgro ./CpHMD-run_$$ 

cp $GROin ./CpHMD-run_$$/${SysName}.gro ; GROin="./${SysName}.gro"
cp $TOPin ./CpHMD-run_$$/${SysName}.top ; TOPin="./${SysName}.top"
cp $NDXin ./CpHMD-run_$$/${SysName}.ndx ; NDXin="./${SysName}.ndx"

## Passing the Plumed settings if they exist to the folder ##
if  [[ -n $plumed && "${plumed}" -eq "1" ]] ; then
    cp ${SysName}_plumed.dat ./CpHMD-run_$$
fi

cd ./CpHMD-run_$$

#
correct_variables #SC 17/3/2010
#
# Let's keep track of the simulation place and time:
echo -e "Simulation run by $USER  @  $HOSTNAME\nInitial time: `date`" \
    > ${blockname}.info
#
# Do some housekeeping
clean_up
#
# Build ff links needed for simulations
build_forcefield
#
# Check if all files are present:
check_files
#
# Make auxiliary files
make_auxiliary_files
#
#Make sites file on the fly
make_sites $1
#
# Get .st files
"$DelphiDir"/getst ${runname}.sites "$StDIR"
#
# Make the delphi database to use
make_delphi_DB
#
#### Starts the constant-pH MD cycle ####
#
TimeStep=`awk '/dt *=/{print $3}' ${runname}.mdp`
#
if [[ $TimeStep == "" ]]; then TimeStep=0.001; fi #SC 25-11-2011
WriteStep=`awk '/nstxout-compressed *=/{print $3}' ${runname}.mdp`
WriteTime=`echo "$WriteStep*$TimeStep" | bc -l`
#
for (( Cycle=1 ; Cycle <=$EndCycle ; Cycle++ )); do
    #
    # This section is to keep track of the simulation
    sim_time=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)+$WriteTime+($StartTime*1000)" | bc -l`
    message W "debug time  $sim_time $EffectiveSteps $TimeStep $Cycle $WriteTime $StartTime "
    echo -e "\nCycle = $Cycle; time = $sim_time ps; Date: `date "+%D %T"`" \
         >> ${blockname}.info
    #
    ################### PB/MC PART #####################
    #
    # Make sure the sites file is not empty...
    sitenumball=$(($(wc -l < ${runname}-all.sites)))
    if [ $sitenumball -ne 0 ]; then
	if [ $ReduceTitration == 1 ]; then
	    if [ $((Cycle % RTInterval)) -eq 1 ]; then
		### Get rid of previous reduced sites ###
		rm -f ${runname}.sites
		### make sure the -all sites is now the new one ###
		cp -f ${runname}-all.sites ${runname}.sites
		#####################################
		# ... and call the PB/MC function,...
		echo -n "PB/MC (All) -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
		     >> ${blockname}.info
		run_PBMC red
		echo "`date "+%D %T"`" >> ${blockname}.info
		write_fractions_all_sites
		#
		# ...write fractions to files and build a new topology.
	    else
		sitenumb=$(($(wc -l < ${runname}.sites)))
		if [ $sitenumb -ne 0 ]; then
		    # ... and call the PB/MC function,...
		    echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
			 >> ${blockname}.info
		    run_PBMC
		    echo "`date "+%D %T"`" >> ${blockname}.info
		    write_fractions
		    
		fi
	    fi
	else
	    sitenumb=$(($(wc -l < ${runname}-all.sites)))
	    if [ $sitenumb -ne 0 ]; then
		# ... and call the PB/MC function,...
		echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
		     >> ${blockname}.info
		run_PBMC
		echo "`date "+%D %T"`" >> ${blockname}.info
		write_fractions_all_sites
	    fi
	fi
	
	update_topology

    else
        # ... skip the PB/MC and write fractions to files 
        message W "File ${runname}.sites is empty. PB/MC step is not performed in cycle $Cycle."
        echo "" >> "TMP_CpHMD.occ"
        echo "" >> "TMP_CpHMD.mocc"
	#migrate the charge topology back to its original name
	mv TMP_processed.top TMP_CpHMD.top
    fi
    
    #### MD PART ####
    # Call dynamics with solvent relaxation 
    if [ $RelaxSteps != 0 ]; then
        echo -n "MD relax     - Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
	run_relaxation #SC 28-11-2011
        echo "`date "+%D %T"`" >> ${blockname}.info
    else
        mv TMP_effective.gro TMP_relax.gro
    fi
    #
    # Call effective (full) dynamics
    echo -n "MD effective - Cycle = $Cycle; Date: `date "+%D %T"` - " \
        >> ${blockname}.info
    run_dynamics effective relax
    echo "`date "+%D %T"`" >> ${blockname}.info
    #
    # Call Append data function   
    data_append
done
#### Ends the constant-pH MD cycle ####
#
# Store Segment Outputs with unambigous Name
for e in gro tpr edr log xtc; do  mv -f TMP_CpHMD.$e ${blockname}.$e; done
## move the energy calculation
if [ -f Eb_calculation.dat ]; then
    mv Eb_calculation.dat ${blockname}.ene
fi
#
if [ -f TMP_CpHMD.occ ]; then
    for e in occ mocc ; do mv -f TMP_CpHMD.$e ${blockname}.$e; done
fi

if [ -f TMP_CpHMD_pullx.xvg ]; then
    for e in x f ; do mv -f TMP_CpHMD_pull$e.xvg ${blockname}_pull$e.xvg; done
fi
# Correct final timestamps of hills and colvar:
if [ -f $colvar_name ]; then
    sed -i '/^ 20.0/d' ${colvar_name}  
fi
#
# Let's keep track of the simulation end time:
echo -e "\nEnd time:     `date`" >> ${blockname}.info
#
# Clean up data function

clean_up

## Transfer production files back to original folder for run cleaning

cp -df ./${SysName}_CpHrun* ../

cp -df ./${SysName}.sites ../

if [ -f ./${SysName}-reducedtitration.sites ] ; then
    cp ./${SysName}-reducedtitration.sites ..//${SysName}_RT-sites.dat
fi
if [ -f ./${SysName}.pocc_RT ] ; then
    cp ./${SysName}.pocc_RT ./${SysName}_RT-debug.pocc_RT
fi

if [[ -n $plumed && "$plumed" == "1" ]] ; then
    if [ -f ./${colvar_name} ] && [ -f ./${hills} ] && [[ $plumedtype != "grid" ]]  ; then
	cp ./${colvar_name} ../
	cp ./${hills} ../
    fi
    #
    if [ -f ./${colvar_name} ] && [[ $plumedtype == "grid" ]]  ; then
	cp ./${colvar_name} ../
	cp  ${hills}_curr_seg ../${hills}
	cp ./${grid_name} ../
    fi
fi

if (for f in ${SysName}_CpHrun*; do diff $f ../$f; done);
then
    cd ../
    gzip  ${SysName}_CpHrun.{log,tpr}
    sleep 3
    rm -rf ./CpHMD-run_$$/*
    rm -rf ./CpHMD-run_$$
else
    message E "Error in file copy... please check local files"
fi
#
exit 0
