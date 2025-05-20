#!/bin/bash -e 
#
prog=`basename $0` 
usage="Usage: $prog <MyProtein_001>.pHmdp\n
       <MyProtein_001>.pHmdp  : Constant-pH MD parameter file"
#
# This allows to use extended pattern matching features
# E.g.: $PDBin == *(" ")
shopt -s extglob #SC 14-12-2011
#
message ()
{
    case "$1" in
        E ) shift; echo -e "$prog: Error: $*" >&2; exit 1;;
        U ) shift; echo -e "$prog: Warning: $*\n$usage" >&2; exit 1;;
        W ) shift; echo -e "$prog: Warning: $*" >&2;;
        * ) message E "Wrong use of 'message' function.";;
    esac
}
#
# Parse arguments
if [ $# != 1 ]; then message U "Wrong number of arguments."; fi
#
# Read Simulation Parameters
source $1
#
# Name of your simulation segment (the output files will
# be generated with this name)
export blockname=${1%.*}
export runname=${blockname%_*}
#
# Defaults section
#
# Ionic Strength MD segment (moles/litre)
#export ionicstrMD=$ionicstr
#
# Dielectric constant of molecular interior
#export epsin=2.0
#
# Dielectric constant of the solvent 
#export epssol=80.0 
#
# Monte Carlo seed used in program petit.
# It can be used to generate different replicates.
#export seed=1234567
#
# Protein Dimension for MEAD (use "0" if you have your own .ogm & .mgm
# files or you are using the delphi module with Delphi)
#export GridSize=81
#
# Residue number OffSet to correct bug in MEAD 
# (the value should be at least twice the max # res. in the protein)
# A value of 0 (zero) bypasses this correction and allows MEAD to use
# the whole residue as the titration site. In this case, .st files
# with pKmod values from Tanford should be used.
#export offset=1000
#
# N-terminus and C-terminus conditions: CAP=capped ; CAPpro=capped 
# proline residue ; CHG=charged ; NEU=Neutral ; TAU1=tautomer1 ;
# TAU2=tautomer2 ; TAU3=tautomer3 ; TAU4=tautomer4 ; REGN=regular 
# neutral ; REGC=regular charged - when titrating (in the .sites file), 
# these settings are bypassed
# (see table in README for more information about tautomers)
#export Nterminus=CAP
#export Cterminus=CAP
#In the case of several chains write instead
#export Nterminus=(CHG CHG)
#export Cterminus=(CHG CHG)
# TODO: try insert this info into titration site list
#
# Solvent relaxation time in steps 
# E.g. 0.2 ps with a 2 fs timestep -> RelaxSteps=100
#export RelaxSteps=100
#
# Effective time in steps (real protein simulation time)
# E.g. 1 ps with a 2 fs timestep -> EffectiveSteps=500
#export EffectiveSteps=1000
#
#The name of the first solvent appearing in the input
#gro file (GROin). It separates the solvent block that will 
#be relaxed after solvent relaxation. 
#export SOL1st="SOL"
# 
# Reduced Titration Threshold (0 will use all sites)
#export RTThreshold=0.001
#
# Reduced Titration Frequency (in # of cycles)
#export RTFrequency=50
#
# The following Programs have to be downloaded and installed
# in your cluster in order to use the Constant-pH MD package
#
# Path to your CpHMD distribution #/home/machuque/gromacs/Constant_pH_MD/CpH-MD_v2.2_CHARMM
#export CpHDIR="/home/filiper/test/sequeira/FF/CpH-MD_v2.2_CHARMM"
#
# Path to your StDIR 
# (the "_Fit" files should be used together with the "offset"
# correction. The older .st work when the offset is bypassed)
#export StDIR=$CpHDIR/ST-$ffID
#
# A pdb file of the system used to know the number of chains 
# and location of the several Nter and Cter. If this variable 
# is empty, the number of chains is considered to be just one.
#export PDBin=MyStructure.pdb
#
# The Rule file for fix_topology (needed in several systems). 
# Whenever present, it will be used. If not needed, leave it empty. 
#export RULEin=""
#
# The POSRE file - To use position restraints, please provide a file;
# otherwise, leave it empty.
#export PosRe=""
# 
# Cutoff for background and pairwise interactions (nanometers) in
# membrane system (PBdim=2); in the case of a simple protein (PBdim=0)
# no cutoff is used
#export cutoff=2.5
#
# PERFIL, GRID SIZE and SCALE cannot be assigned at the same
# time. They are not independent variables.
# (for more info check the "MyProtein_XXX.pbp" file in this folder)
# export perfil=80
#export gsize=81
#export scaleP=1.0
#export scaleM=4.0
#
# Potential at Boundary (1, 2, 3 or 4):
#     1. Potential is zero (0.0) 
#     2. Dipole. The boundary potentials are approximated by the
#        Debye-Huckel potential of the equivalent dipole to the
#        molecular charge distribution.
#     3. Focusing. 
#     4. Coulombic. They are approximated by the sum of Debye-Huckel
#        potentials of all the charges.
#export bndcon=4
#
# The convergence threshold values:
# (for more info check the "MyProtein_XXX.pbp" file in this folder)
# fast: 0.01
# slow: 0.0001
#export maxc=0.01
#
# Define beginning of the Cycle:
#export InitCycle=1
#
# Define end of the Cycle:
#export EndCycle=500
# Option "-rcon" in mdrun allows the user to "decompose" smaller systems
# in several CPUs. As an example, in a small 3k atoms system the default
# option "0" estimates a value of 1.3 nm which does not allow the use of
# more than 1 CPU. The use of a value of 1.0 nm allows the use of at 
# least 4 CPUs in that simulation.
#
#export Rcon="0"
#
# Define variable for groswitch
# Currently, CpHMD only works with tautomers.
export gstaut=1
#
# Check some parameters
if [[ $ffID != G54a7pH && $ffID != CHARMM36pH && $ffID != Amber14SBpH ]]; then
    message E "ffID = $ffID is not valid. Check $1."
fi
#
# Call all functions stored in the "functions" file:
functions=$CpHDIR/scripts/functions.sh
if [[ -f $functions ]]; then
    source $functions            
else
    message E "File $functions is missing. Check CpHDIR parameter in $1."
fi
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
# Check the number of chains #SC 15-11-2011
# check_nchains #SC 15-11-2011
# Set termini states if they are not titrating
set_termini
#
# Get .st files
"$MToolsDIR"/getst ${runname}.sites "$StDIR"
#
# Make auxiliary files
make_auxiliary_files
#
#Commented by SC 30-11-2011
# Make files that remain static when not using the Reduced Titration routine
#if [ $RTThreshold == 0 ]; then
#    echo "Protein" | "$GroDIR"/trjconv -f $GROin -s TMP_CpHMD.tpr \
#        -o TMP_CpHMD_red.gro -n TMP_CpHMD.ndx -quiet
#fi
#
#### Starts the constant-pH MD cycle ####
#
TimeStep=`awk '/^dt *=/{print $3}' ${runname}.mdp`
#time step 0.001 ps is the default value in gromacs; SC 25-11-2011
if [[ $TimeStep == "" ]]; then TimeStep=0.001; fi #SC 25-11-2011
WriteStep=`awk '/^nstxout-compressed *=/{print $3}' ${runname}.mdp`
WriteTime=`echo "$WriteStep*$TimeStep" | bc -l`
#
for (( Cycle=$InitCycle ; Cycle <=$EndCycle ; Cycle++ )); do
    #
    # This section is just an index to keep track of the simulation
    sim_time=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)+$WriteTime" | bc -l`
    echo -e "\nCycle = $Cycle; time = $sim_time ps; Date: `date "+%D %T"`" \
        >> ${blockname}.info
    #
    ################### PB/MC PART #####################
    #
    # Call the Reduced Titration function when needed
    if [ $((Cycle % RTFrequency)) -eq 1 -a $RTThreshold != 0 ]; then
        echo -n "PB/MC (All) -  Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
        run_PBMC red
        echo "`date "+%D %T"`" >> ${blockname}.info
    fi
    #
    # Make sure the sites file is not empty...
    sitenumb=$(($(wc -l < ${runname}.sites)))
    if [ $sitenumb -ne 0 ]; then
        # ... and call the PB/MC function,...
        echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
        run_PBMC std
        echo "`date "+%D %T"`" >> ${blockname}.info
        #
        # ...write fractions to files and build a new topology.
        write_fractions
        build_topology
    #
    # Otherwise...
    else
        # ... skip the PB/MC and write fractions to files 
        #(when RTThreshold = 0 there are no fractions to write;
        # usefull for tests without PB/MC).
        message W "File ${runname}.sites is empty. PB/MC step is not performed in cycle $Cycle."
        if [ $RTThreshold != 0 ]; then 
            if [ $((Cycle % RTFrequency)) -eq 1 ]; then
            #Use the output from the Reduced Titration routine
                mv -f TMP_MCarlo_mod.out TMP_MCarlo_std.out
                write_fractions
                build_topology 
            else  
                write_fractions
            fi
        fi
    fi
    #
    #### MD PART ####
    # Call dynamics with solvent relaxation 
    if [ $RelaxSteps != 0 ]; then
        echo -n "MD relax     - Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
        #run_dynamics relax effective
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
# Rename original .sites file
mv -f ${runname}-all.sites ${runname}.sites
#
# Store Segment Outputs with unambigous Name
for e in gro tpr edr log xtc; do  mv -f TMP_CpHMD.$e ${blockname}.$e; done
#
if [ $RTThreshold != 0 ]; then
    for e in sites occ_red{,_mod} mocc_red ; do 
        mv -f TMP_CpHMD.$e ${blockname}.$e; done
fi
#
if [ -f TMP_CpHMD.occ ]; then
    for e in occ mocc ; do mv -f TMP_CpHMD.$e ${blockname}.$e; done
fi
#
# Let's keep track of the simulation end time:
echo -e "\nEnd time:     `date`" >> ${blockname}.info
#
# Clean up function
clean_up
#
exit 0
