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

# ------------------------------------------------------------
# Select the requested GROMACS installation.
#
# The dual-build container provides two independent executables:
#   GPU=0 -> GroDIRCPU (CPU-only GROMACS)
#   GPU=1 -> GroDIRGPU (CUDA-enabled GROMACS)
#
# Mode-specific variables take precedence over GroDIR. GroDIR remains an
# expert fallback for custom containers that provide only one executable path.
# ------------------------------------------------------------
case "${GPU:-0}" in
    0)
        selected_gromacs="${GroDIRCPU:-${GroDIR:-/gromacs/bin/gmx}}"
        selected_mode="CPU"
        ;;
    1)
        selected_gromacs="${GroDIRGPU:-${GroDIR:-/gromacs/bin/gmx-gpu}}"
        selected_mode="CUDA"
        ;;
    *)
        message E "GPU must be either 0 (CPU-only GROMACS) or 1 (CUDA GROMACS); received: ${GPU:-unset}"
        ;;
esac

# This inexpensive path check catches a missing or mistyped installation
# immediately. Compile-mode and GPU-device checks run later, after functions.sh
# has been loaded, but still before scratch setup and PB/MC.
if [[ ! -x "${selected_gromacs}" ]]; then
    message E "Selected ${selected_mode} GROMACS executable is missing or not executable: ${selected_gromacs}"
fi

export GroDIR="${selected_gromacs}"
message W "Selected ${selected_mode} GROMACS executable: ${GroDIR}"
#
# Name of your simulation segment (the output files will
# be generated with this name)
#export blockname=${1%.*}
export runname=${SysName}
export blockname="${SysName}_CpHrun"


# Explicit site/state names. The external make_sites helper still creates a
# temporary ${runname}.sites internally, but CpHMD immediately renames it.
export CPHMD_ALL_SITES_FILE="${runname}-all.sites"
export CPHMD_RT_ACTIVE_SITES_FILE="${runname}-RT-active.sites"
export CPHMD_RT_HISTORY_FILE="${runname}-RT-history.sites"
export CPHMD_RT_TEMPLATE_OCC_FILE="TMP_RT_template.occ"
export CPHMD_RT_TEMPLATE_MOCC_FILE="TMP_RT_template.mocc"
export CPHMD_RT_CURRENT_OCC_FILE="TMP_RT_current.occ"
export CPHMD_RT_CURRENT_MOCC_FILE="TMP_RT_current.mocc"

# Fixed short basename used only for disposable legacy Delphi/PETIT files.
# Keeping this independent of SysName avoids silent Fortran filename failures.
export CPHMD_PBMC_BASENAME="PBMC"

# ------------------------------------------------------------
# Application-owned output layout.
#
# When the SLURM wrapper is used, the host's node-local scratch directory is
# mounted as /out and CpHMD.sh starts in /out. These folders are therefore
# created on local scratch, not in the shared project directory. The wrapper
# later synchronizes only finalized contents back to persistent storage.
#
# These names are a fixed output contract shared by CpHMD.sh and the scheduler
# wrappers. Keeping them fixed avoids a settings/wrapper mismatch. CpHMD owns
# the scientific directories and checkpoint file. The wrapper owns
# slurm_files/, RUNNING/DONE/FAILED and resub_count.
# ------------------------------------------------------------
export CPHMD_SEGMENT_DIR="segments"
export CPHMD_RESTART_DIR="restart"
export CPHMD_LOG_DIR="log-files"
export CPHMD_RESUB_DIR="resub_maint"

for output_dir in \
    "${CPHMD_SEGMENT_DIR}" \
    "${CPHMD_RESTART_DIR}" \
    "${CPHMD_LOG_DIR}" \
    "${CPHMD_RESUB_DIR}"
do
    if [[ -z "${output_dir}" || "${output_dir}" == /* || "${output_dir}" == *".."* ]]; then
        message E "CpHMD output directories must be non-empty relative paths without '..': ${output_dir}"
    fi

    mkdir -p -- "${output_dir}" \
        || message E "Could not create CpHMD output directory ${output_dir}."
done

message W "CpHMD output layout: ${CPHMD_SEGMENT_DIR}/, ${CPHMD_RESTART_DIR}/, ${CPHMD_LOG_DIR}/ and ${CPHMD_RESUB_DIR}/."

# ------------------------------------------------------------
# Cycle-level checkpoint/restart support.
#
# Default behaviour is unchanged when CPHMD_MAX_SECONDS=0 and no
# checkpoint file is present. When enabled by a wrapper, the script
# writes compact restart files after complete CpHMD cycles and exits
# cleanly with code 10 at a cycle boundary once the walltime guard is
# reached. This is useful for long systems that must be split across
# several scheduler jobs.
# ------------------------------------------------------------
export CPHMD_CHECKPOINT_FILE="${CPHMD_CHECKPOINT_FILE:-${CPHMD_RESUB_DIR}/checkpoint.state}"
export CPHMD_WALL_START_EPOCH="${CPHMD_WALL_START_EPOCH:-$(date +%s)}"
export CPHMD_MAX_SECONDS="${CPHMD_MAX_SECONDS:-0}"

# Named application exit codes. 10 is reserved for a validated clean scheduler
# boundary; 42 denotes an MD/output contract failure and 43 a PB/MC failure.
# Neither failure code is eligible for scheduler resubmission.
readonly CPHMD_EXIT_CLEAN_RESUBMIT=10
readonly CPHMD_EXIT_OUTPUT_FAILURE=42
readonly CPHMD_EXIT_PBMC_FAILURE=43
export CPHMD_EXIT_CLEAN_RESUBMIT CPHMD_EXIT_OUTPUT_FAILURE CPHMD_EXIT_PBMC_FAILURE

export CPHMD_DONE_CYCLES=0
export CPHMD_LAST_SEGMENT_START=0
export CPHMD_LAST_SEGMENT_END=0
export CPHMD_SEGMENT_FINALIZED=0
export CPHMD_REDUCED_TITRATION="${ReduceTitration:-0}"
export CPHMD_RT_LAST_FULL_CYCLE=0
export CPHMD_RT_ACTIVE_COUNT=0

if [[ -f "${CPHMD_CHECKPOINT_FILE}" ]]; then
    # This file is generated by cphmd_write_checkpoint in functions.sh.
    # It contains simple shell assignments only.
    source "${CPHMD_CHECKPOINT_FILE}"
fi

if [[ -z "${CPHMD_DONE_CYCLES:-}" ]]; then
    export CPHMD_DONE_CYCLES=0
fi

# If restarting from a previous clean cycle-boundary checkpoint, use the
# compact restart state rather than the original input coordinates/topology.
if (( CPHMD_DONE_CYCLES > 0 )); then
    message W "Restart checkpoint found: ${CPHMD_DONE_CYCLES} completed cycles."

    # A new job may continue only from a scheduler segment that was completely
    # validated and finalized. A cycle-level scratch checkpoint is deliberately
    # insufficient because its XTC/EDR may not have been committed.
    if (( ${CPHMD_SEGMENT_FINALIZED:-0} != 1 )); then
        message E "Checkpoint ${CPHMD_CHECKPOINT_FILE} is not marked as a finalized segment. Refusing automatic continuation."
    fi

    restart_gro="${CPHMD_RESTART_DIR}/${blockname}.gro"
    restart_top="${CPHMD_RESTART_DIR}/${blockname}.top"

    if [[ -f "${restart_gro}" ]]; then
        export GROin="./${restart_gro}"
    else
        message E "Checkpoint says cycles were completed, but ${restart_gro} is missing."
    fi

    if [[ -f "${restart_top}" ]]; then
        export TOPin="./${restart_top}"
    else
        message E "Checkpoint says cycles were completed, but ${restart_top} is missing."
    fi
fi
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

############################################################
# Build optional input files embedded in the settings file. #
############################################################
# The template may always contain #fixgro# and #plumed# sections. They are
# extracted only when their corresponding workflow is active; otherwise they
# remain inert comments in the settings file and no extra working file is made.

settings_file="$1"
mdp_file="${SysName}.mdp"
fixgro_file="${SysName}.fixgro"
plumed_file="${SysName}_plumed.dat"

awk '/#mdp#[[:space:]]/ {sub(/^.*#mdp#[[:space:]]*/, ""); print}' \
    "${settings_file}" > "${mdp_file}"
if ! grep -Eq '^[[:space:]]*[^;#[:space:]]' "${mdp_file}"; then
    message E "The settings file contains no usable #mdp# block."
fi

# Validate the simple feature flags before deciding which optional blocks to use.
for flag_name in memb multiple_solutes include_itp plumed; do
    flag_value="${!flag_name:-0}"
    case "${flag_value}" in
        0|1) ;;
        *) message E "${flag_name} must be 0 or 1; received '${flag_value}'." ;;
    esac
done

rm -f -- "${fixgro_file}"
export CPHMD_FIXGRO_REQUIRED=0
if cphmd_fixgro_required; then
    export CPHMD_FIXGRO_REQUIRED=1
    awk '/#fixgro#[[:space:]]/ {sub(/^.*#fixgro#[[:space:]]*/, ""); print}' \
        "${settings_file}" > "${fixgro_file}"
    cphmd_validate_fixgro_file "${fixgro_file}"
    message W "Using the #fixgro# settings block for non-membrane PBC repair."
else
    message W "The #fixgro# settings block is inactive for this system configuration."
fi

rm -f -- "${plumed_file}"
if (( plumed == 1 )); then
    case "${plumedtype:-}" in
        grid|hill|static) ;;
        *) message E "plumedtype must be grid, hill or static when plumed=1; received '${plumedtype:-unset}'." ;;
    esac

    awk '/#plumed#[[:space:]]/ {sub(/^.*#plumed#[[:space:]]*/, ""); print}' \
        "${settings_file}" > "${plumed_file}"

    # Replace the documented template tokens. The corresponding variables are
    # required only while PLUMED is active.
    : "${colvar_stride:?Set colvar_stride when plumed=1}"
    : "${colvar_name:?Set colvar_name when plumed=1}"
    : "${hills:?Set hills when plumed=1}"

    sed -i \
        -e "s|&colvar_stride|${colvar_stride}|g" \
        -e "s|&colvar_name|${colvar_name}|g" \
        -e "s|&hills|${hills}|g" \
        "${plumed_file}"

    if [[ ${plumedtype} == "grid" ]]; then
        : "${grid_name:?Set grid_name when plumedtype=grid}"
        sed -i -e "s|&grid_name|${grid_name}|g" "${plumed_file}"
    fi

    cphmd_validate_plumed_file "${plumed_file}"

    # PLUMED restart behavior is required for continuation. For the first cycle
    # of a segment, run_dynamics writes a temporary copy without this directive.
    if ! grep -Eq '^[[:space:]]*RESTART([[:space:]]|$)' "${plumed_file}"; then
        plumed_tmp="${plumed_file}.tmp.$$"
        {
            echo "RESTART"
            cat "${plumed_file}"
        } > "${plumed_tmp}" || message E "Could not add RESTART to ${plumed_file}."
        mv -f -- "${plumed_tmp}" "${plumed_file}" \
            || message E "Could not finalize ${plumed_file}."
    fi
    message W "PLUMED is enabled with plumedtype=${plumedtype}."
else
    message W "The #plumed# settings block is inactive because plumed=0."
fi

#############################################################
# add the correct cycles depending on the size of block asked
#############################################################
dt=`awk '/dt/ {print $3}' ${SysName}.mdp`
if [ -z "$dt" ] ;
then
    dt=0.002
    message W "dt was not defined in the .mdp, falling back to 0.002 default dt."
fi

EndCycle=`echo ${SimTime} $EffectiveSteps $dt | awk '{print int(($1*1000)/($2*$3)+0.5)}'`

if (( CPHMD_DONE_CYCLES >= EndCycle )); then
    message W "CpHMD checkpoint indicates that ${CPHMD_DONE_CYCLES}/${EndCycle} cycles are already complete."
    exit 0
fi

FirstCycle=$((CPHMD_DONE_CYCLES + 1))
export CPHMD_SEGMENT_START_CYCLE="${FirstCycle}"

message W "Running CpHMD cycles from ${FirstCycle} until ${EndCycle} with 20 ps each cycle."
message W "Checkpoint restart: completed cycles = ${CPHMD_DONE_CYCLES}; starting at cycle ${FirstCycle} of ${EndCycle}."
if (( CPHMD_MAX_SECONDS > 0 )); then
    message W "CpHMD cycle-boundary walltime guard: ${CPHMD_MAX_SECONDS} seconds."
fi

#################################################
## PLUMED Detection rebasing ##
## Both CPU and CUDA GROMACS builds use the   ##
## same native PLUMED integration.             ##
#################################################
#if [[ $plumed == 1 ]] ;then
#    ## change /gromacs file if it is the standard one ##
#    if [ $GroDIR == "/gromacs/bin/gmx" ] ; then
#	export GroDIR="/gromacs-plumed/bin/gmx"
#	message W "PLUMED support requested, changing the gromacs compilation to PLUMED $GroDIR"
#    fi
#fi

#############################################################
## Build mdrun command and perform fail-fast runtime checks ##
#############################################################
# This block runs before the node-local scratch directory is created and before
# any PB/MC work. It validates:
#   * GPU/nCPU and safety-flag values;
#   * the selected GROMACS executable and compile mode;
#   * for GPU=1, actual NVIDIA device/driver visibility inside the container.
#
# With the default CPHMD_GPU_PREFLIGHT=1, a CUDA build launched on a CPU-only
# host or without Apptainer/Singularity --nv terminates here. A second,
# zero-step mdrun probe is performed after the initial TPR is generated, still
# before any PB/MC work.
cphmd_build_mdrun_command
cphmd_preflight_gromacs_runtime

## Creating a running folder to prevent clutter ##
mkdir -p ./CpHMD-run_$$

### Checking if the specific files exist and copying them to the running directory
for f in ${SysName}.mdp $GROin $TOPin $NDXin 
do
    if [ ! -f $f ]; then
        message E  "File $f is missing!!!... Program will crash"
    fi
done

cp "${settings_file}" "./CpHMD-run_$$/"
cp "${mdp_file}" "./CpHMD-run_$$/"
if (( CPHMD_FIXGRO_REQUIRED == 1 )); then
    cp "${fixgro_file}" "./CpHMD-run_$$/"
fi

cp $GROin ./CpHMD-run_$$/${SysName}.gro ; GROin="./${SysName}.gro"
cp $TOPin ./CpHMD-run_$$/${SysName}.top ; TOPin="./${SysName}.top"
cp $NDXin ./CpHMD-run_$$/${SysName}.ndx ; NDXin="./${SysName}.ndx"

## Passing the Plumed settings if they exist to the folder ##
if (( plumed == 1 )); then
    cp "${plumed_file}" "./CpHMD-run_$$/"
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
# Make auxiliary files, including the initial TPR used by the execution probe.
make_auxiliary_files || exit $?
#
# Run the actual configured mdrun command for zero integration steps. This is
# deliberately before site generation and PB/MC so invalid CPU/GPU execution
# settings fail without spending time in electrostatics/Monte Carlo work.
cphmd_preflight_mdrun_execution || exit $?
#
# Make the immutable all-site definition on the fly. The temporary legacy
# ${runname}.sites name is removed inside make_sites.
make_sites "$1"
# Restore the active RT subset/templates on checkpoint continuation, or prepare
# a clean first-run state. This occurs before any PB/MC cycle.
cphmd_initialize_or_restore_rt_state || exit $?
#
# State templates are prepared for the immutable all-site definition. Reduced
# cycles use a subset of these same site definitions.
"$DelphiDir"/getst "${CPHMD_ALL_SITES_FILE}" "$StDIR" ||     message E "getst failed for ${CPHMD_ALL_SITES_FILE}."
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
for (( Cycle=FirstCycle ; Cycle <=$EndCycle ; Cycle++ )); do
    #
    # This section is to keep track of the simulation
    sim_time=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)+$WriteTime+($StartTime*1000)" | bc -l`
    message W "debug time  $sim_time $EffectiveSteps $TimeStep $Cycle $WriteTime $StartTime "
    echo -e "\nCycle = $Cycle; time = $sim_time ps; Date: `date "+%D %T"`" \
         >> ${blockname}.info
    #
    ################### PB/MC PART #####################
    all_site_count="$(cphmd_count_nonblank "${CPHMD_ALL_SITES_FILE}")" || {
        cphmd_preserve_failure_diagnostics "${Cycle}"
        exit "${CPHMD_EXIT_OUTPUT_FAILURE}"
    }

    if (( all_site_count > 0 )); then
        if (( ReduceTitration == 1 )); then
            if cphmd_rt_full_refresh_cycle "${Cycle}"; then
                echo -n "PB/MC (All) -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
                    >> ${blockname}.info
                run_PBMC "${CPHMD_ALL_SITES_FILE}" "all-site refresh"
                rc=$?
                if (( rc != 0 )); then
                    cphmd_preserve_failure_diagnostics "${Cycle}"
                    exit "${rc}"
                fi
                cphmd_refresh_reduced_titration_state
                rc=$?
                if (( rc != 0 )); then
                    cphmd_preserve_failure_diagnostics "${Cycle}"
                    exit "${rc}"
                fi
                echo "`date "+%D %T"`" >> ${blockname}.info
            else
                active_site_count="$(cphmd_count_nonblank "${CPHMD_RT_ACTIVE_SITES_FILE}")" || {
                    cphmd_preserve_failure_diagnostics "${Cycle}"
                    exit "${CPHMD_EXIT_OUTPUT_FAILURE}"
                }

                if (( active_site_count > 0 )); then
                    echo -n "PB/MC (RT ${active_site_count}/${all_site_count}) - Cycle = $Cycle; Date: `date "+%D %T"` - " \
                        >> ${blockname}.info
                    run_PBMC "${CPHMD_RT_ACTIVE_SITES_FILE}" "reduced active-set"
                    rc=$?
                    if (( rc != 0 )); then
                        cphmd_preserve_failure_diagnostics "${Cycle}"
                        exit "${rc}"
                    fi
                    cphmd_merge_reduced_titration_cycle
                    rc=$?
                    if (( rc != 0 )); then
                        cphmd_preserve_failure_diagnostics "${Cycle}"
                        exit "${rc}"
                    fi
                    echo "`date "+%D %T"`" >> ${blockname}.info
                else
                    echo "PB/MC (RT 0/${all_site_count}; fixed-state reuse) - Cycle = $Cycle; Date: `date "+%D %T"`" \
                        >> ${blockname}.info
                    message W \
                        "Cycle ${Cycle}: active reduced-titration site set is empty; reusing validated fixed-state templates without running PB/MC."
                    cphmd_materialize_fixed_rt_cycle
                    rc=$?
                    if (( rc != 0 )); then
                        cphmd_preserve_failure_diagnostics "${Cycle}"
                        exit "${rc}"
                    fi
                fi
            fi
        else
            echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
                >> ${blockname}.info
            run_PBMC "${CPHMD_ALL_SITES_FILE}" "all-site"
            rc=$?
            if (( rc != 0 )); then
                cphmd_preserve_failure_diagnostics "${Cycle}"
                exit "${rc}"
            fi
            cphmd_materialize_nonrt_cycle
            rc=$?
            if (( rc != 0 )); then
                cphmd_preserve_failure_diagnostics "${Cycle}"
                exit "${rc}"
            fi
            echo "`date "+%D %T"`" >> ${blockname}.info
        fi

        cphmd_record_current_fractions
        rc=$?
        if (( rc != 0 )); then
            cphmd_preserve_failure_diagnostics "${Cycle}"
            exit "${rc}"
        fi

        update_topology
        rc=$?
        if (( rc != 0 )); then
            cphmd_preserve_failure_diagnostics "${Cycle}"
            exit "${rc}"
        fi
    else
        message W "All-site definition ${CPHMD_ALL_SITES_FILE} is empty. PB/MC is not performed in cycle ${Cycle}."
        echo "" >> TMP_CpHMD.occ
        echo "" >> TMP_CpHMD.mocc
        mv TMP_processed.top TMP_CpHMD.top
    fi

    #### MD PART ####
    # Each stage now propagates command failures and missing-output errors.
    # Continuing after one of these failures could reuse stale coordinates and
    # would invalidate both the structural and protonation-state trajectory.
    if [ $RelaxSteps != 0 ]; then
        echo -n "MD relax     - Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
        run_relaxation >> ${blockname}_MD.log 2>&1
        rc=$?
        if (( rc != 0 )); then
            cphmd_preserve_failure_diagnostics "${Cycle}"
            exit "${rc}"
        fi
        echo "`date "+%D %T"`" >> ${blockname}.info
    else
        mv TMP_effective.gro TMP_relax.gro
    fi

    echo -n "MD effective - Cycle = $Cycle; Date: `date "+%D %T"` - " \
        >> ${blockname}.info
    run_dynamics effective relax >> ${blockname}_MD.log 2>&1
    rc=$?
    if (( rc != 0 )); then
        cphmd_preserve_failure_diagnostics "${Cycle}"
        exit "${rc}"
    fi
    echo "`date "+%D %T"`" >> ${blockname}.info

    data_append
    rc=$?
    if (( rc != 0 )); then
        cphmd_preserve_failure_diagnostics "${Cycle}"
        exit "${rc}"
    fi

    # This performs cheap file/size checks and atomic copies only. The state is
    # marked non-finalized until segment-level binary validation succeeds.
    cphmd_write_checkpoint "${Cycle}"
    rc=$?
    if (( rc != 0 )); then
        cphmd_preserve_failure_diagnostics "${Cycle}"
        exit "${rc}"
    fi

    if (( CPHMD_MAX_SECONDS > 0 && Cycle < EndCycle )); then
        now_epoch=$(date +%s)
        elapsed=$((now_epoch - CPHMD_WALL_START_EPOCH))

        if (( elapsed >= CPHMD_MAX_SECONDS )); then
            cphmd_finalize_segment "${Cycle}"
            rc=$?
            if (( rc != 0 )); then
                cphmd_preserve_failure_diagnostics "${Cycle}"
                exit "${rc}"
            fi
            message W "Reached validated CpHMD cycle-boundary walltime guard after cycle ${Cycle}/${EndCycle}; exiting cleanly for scheduler resubmission."
            exit "${CPHMD_EXIT_CLEAN_RESUBMIT}"
        fi
    fi
done
#### Ends the constant-pH MD cycle ####
#
# Add the normal-completion timestamp before the final segment is committed.
echo -e "\nEnd time:     `date`" >> ${blockname}.info
final_log_base="../${CPHMD_LOG_DIR}/${blockname}_cycles${CPHMD_SEGMENT_START_CYCLE}-${EndCycle}"
cphmd_copy_atomic "${blockname}.info" "${final_log_base}.info"

# The last cycle was checkpointed inside the loop. Validate and finalize the
# final scheduler segment once, including optional compression. No duplicate
# *.final.* files are created; restart/ contains the authoritative endpoint.
cphmd_finalize_segment "${EndCycle}"
rc=$?
if (( rc != 0 )); then
    cphmd_preserve_failure_diagnostics "${EndCycle}"
    exit "${rc}"
fi

# Correct final timestamps of hills and colvar, preserving the historical
# PLUMED post-processing behaviour when metadynamics is active.
if [[ -n ${colvar_name:-} && -f ${colvar_name} ]]; then
    sed -i '/^ 20.0/d' ${colvar_name}
fi

#
# Clean up data function

clean_up

## Transfer production files back to original folder for run cleaning
## With the checkpoint-aware workflow, files have already been copied to the
## parent scratch directory by cphmd_write_checkpoint/finalization above.
## The scheduler wrapper can then rsync selected files from scratch to /home.

cd ../
sleep 3
rm -rf ./CpHMD-run_$$

#
exit 0
