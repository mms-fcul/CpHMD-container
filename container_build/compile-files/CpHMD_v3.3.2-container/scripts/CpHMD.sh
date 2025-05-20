#!/bin/bash -e 
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
    message E "ffID = $ffID is not valid. Check $1."
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


## Creating a running folder to prevent clutter ##
mkdir -p ./CpHMD-run_$$
cp ./${SysName}.mdp ./CpHMD-run_$$ ; cp ./${SysName}.fixgro ./CpHMD-run_$$ ; cp $1  ./CpHMD-run_$$

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
    sitenumb=$(($(wc -l < ${runname}.sites)))
    if [ $sitenumb -ne 0 ]; then
        # ... and call the PB/MC function,...
        echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
        run_PBMC 
        echo "`date "+%D %T"`" >> ${blockname}.info	
        #
        # ...write fractions to files and build a new topology.	
	write_fractions

	update_topology
        
    #
    # Otherwise...
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
        #
        #This was passed to functions; SC 28-11-2011
        # Prepare input GRO for dynamics
	#awk -v s=$SOL1st '$1 ~ s {exit};{print $0}' TMP_effective.gro > TMP_aux.gro
        #awk -v s=$SOL1st '$1 ~ s {a=1};a'  TMP_relax.gro >> TMP_aux.gro
	
        #mv -f TMP_aux.gro TMP_relax.gro
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

#
# Let's keep track of the simulation end time:
echo -e "\nEnd time:     `date`" >> ${blockname}.info
#
# Clean up data function

clean_up

## Transfer production files back to original folder for run cleaning

cp -df ./${SysName}_CpHrun* ../

cp -df ./${SysName}.sites ../


if (for f in ${SysName}_CpHrun*; do diff $f ../$f; done);
then
    cd ../
    rmdir --ignore-fail-on-non-empty ./CpHMD-run_$$
    gzip -9 ${SysName}_CpHrun.{err,log,tpr}
else
    message E "Error in file copy... please check local files"
fi
#
exit 0
