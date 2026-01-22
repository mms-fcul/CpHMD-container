# This file contains all the functions called by CpH-MD.sh script
#It also contains the definition of some constants.
#
# Please don't delete it
#
#######################################################################
# SC 28-11-2011
# In this part some constants can be defined which may change with 
# gromacs version. When adapting to new gromacs versions, please check 
# in here.
#
# Needed to know how to delimit moleculetype section in pdb2gmx output 
# topology. 
# Careful: it may change with gromacs version
# This one works with 3.2.1 and 4.0.7 at least
export TOPSTRING="; Include water topology" 
#
#######################################################################
#
correct_variables ()
{                   
    #Turns "0.0", "0.00", "0.000", etc into "0"
    GridSize=`awk -v i=$GridSize 'BEGIN{print (i==0) ? 0 : i}'`
    offset=`awk -v i=$offset 'BEGIN{print (i==0) ? 0 : i}'` 
    #Turns "" into "0" #MM 03/10/2011
    Rcon=`awk -v i="$Rcon" 'BEGIN{print i+0}'`
}

build_forcefield ()
{
    # Call the modified forcefield:
    ln -s -f "$CpHDIR" .; ln -s -f "$CpHDIR"/*.dat .
}

check_files ()
{
    # Check for other necessary files (filenames must be respected)
    for f in ${runname}.mdp $GROin $TOPin $GroDIR
      do
      if [ ! -f $f ]; then
          echo message E  "File $f is missing!!!... Program will crash"
      fi
    done
    # TODO: create a check for the case the site list has no titrating groups
    #
    #Check topology file for delimiting strings ;SC 28-11-2011
    awk -v p=$prog '
    BEGIN{
    b_string=";[ begin_tit_molecule ]"; 
    e_string=";[ end_tit_molecule ]";
    b=0;
    end=0;
    order=1;
    warnfmt=(p": Warning: String %s appears repeated in %s file. First occurrence will be considered.\n");
    errfmt=(p": Error: String %s is missing in %s file.\n")
    }
    /^ *; *\[ +begin_tit_molecule +\]/{
    if (b) printf warnfmt, b_string, ARGV[1]; else b=1}
    /^ *; *\[ +end_tit_molecule +\]/{
    if (end) printf warnfmt, e_string, ARGV[1]; 
    else {if (!b) order=0; end=1}}
    END{
    if (!b) {printf errfmt, b_string, ARGV[1]; exit 1};
    if (!end) {printf errfmt, e_string, ARGV[1]; exit 1};
    if (!order) {printf (p": Error: String %s is before string %s in %s file.\n"), e_string, b_string, ARGV[1]; exit 1}
    }' $TOPin >&2  
    #
    # Check directories
    for d in $StDIR
      do
      if [ ! -d $d ]; then 
          message E  "Directory $d does not exist!!!... Program will crash"
      fi
    done    
    if [[ $PBdim -ne 0 && $PBdim -ne 2 ]]; then
        message E  "Number of dimensions ($PBdim) is not supported... Program will crash"
    fi

}

set_termini()
{
    # TODO: create condition for G54A7 and for CHARMM36

    
    # This is FF hack specific. 
    # Capping Conditions
    if ! egrep NT ${runname}.sites; then 
        if [[ ${#Nterminus[@]} -ne 1 ]]; then #SC 30-11-2011
            message E  "Nterminus array contains wrong number of elements (${#Nterminus[@]}). One element was expected. Please check ${blockname}.pHmdp."   
        fi
        case $Nterminus in 
                CAP )           nter=7 ;;
                CAPpro )        nter=11 ;; 
                CHG )           nter=3 ;; 
                NEU | TAU1 )    nter=0 ;;
                TAU2 )          nter=1 ;;
                TAU3 | TAU4 )    nter=2 ;;
                * ) message E "Nterminus = $Nterminus is not valid. Check .pHmdp file" ;;
            esac
        fi
        if ! egrep CT ${runname}.sites; then 
            if [[ ${#Cterminus[@]} -ne 1 ]]; then #SC 30-11-2011
                message E  "Cterminus array contains wrong number of elements (${#Cterminus[@]}). One element was expected. Please check ${blockname}.pHmdp."
            fi  
            case $Cterminus in 
                CAP | CAPpro )  cter=9 ;; 
                CHG )           cter=4 ;; 
                NEU | TAU1 )    cter=0 ;;
                TAU2 )          cter=1 ;;
                TAU3 )          cter=2 ;;
                TAU4 )          cter=3 ;;
                * ) echo message E "Cterminus = $Cterminus is not valid. Check .pHmdp file" ;;
            esac
        fi
}

make_auxiliary_files ()
{
    # Rename original files
    cp -f ${runname}.sites ${runname}-all.sites
    cp -f $TOPin TMP_CpHMD.top
    cp -f $GROin TMP_effective.gro
    #
    # Create auxiliary files with the unchanged part of the topology file; SC 28-11-2011
    awk '/^ *; *\[ +begin_tit_molecule +\]/{exit}{print}' $TOPin > TMP_CpHMD.topbegin 
    awk 'n;/^ *; *\[ +end_tit_molecule +\]/{n=1}' $TOPin > TMP_CpHMD.topend
  
    # Check if there are different ionic strength values in MD and PB
    if [[ $ionicstrMD == "" ]]; then ionicstrMD=$ionicstr; fi
    # Convert Ionic Strength from Molar to molecule/nm^3.
    ionicstrMolecule=$(awk -v i=${ionicstrMD} 'BEGIN{print i*0.6022}')
    #
    # Correct Ionic Strength and Temperature in the .mdp file according 
    # to your parameters
    sed "s/\(ionicstrength *= \).*\$/\1 $ionicstrMolecule/" \
        ${runname}.mdp > TMP_aux1.mdp
    message W "Ionic strength in ${runname}.mdp file was changed to $ionicstrMolecule molecules/nm^3 (which corresponds to $ionicstr Molar)."
    #
    tcgroups=`awk '$1=="tc-grps"{for (i=3;i<=NF;i++) if ($i !~ /^;/) {n++; if ($i ~ /;/) exit} else exit}END{print n}' ${runname}.mdp` #SC 25-11-2011
    awk -v e=$tcgroups -v t=$temp '!/^ref_t *=/{print $0}; 
        /^ref_t *=/{printf "ref_t               =  "; 
        for (i=1;i<=e;i++) printf "%-9s", t ;print ""}' \
            TMP_aux1.mdp > ${runname}.mdp #SC 25-11-2011
        message W "Temperature in ${runname}.mdp file was changed to $temp."
    #
    # Create a completely normalized index file
    #Allow the option of giving index file ; SC 28-11-2011
    if [[ $NDXin != *(" ") && -f $NDXin ]]
    then
        cp -f $NDXin TMP_CpHMD.ndx
        message W "File $NDXin will be used as index file instead of file automatically generated by make_ndx."
    else
        if [ ! -f $NDXin ]
        then
            message W "File $NDXin is missing. A file generated by make_ndx will be used instead."
        fi
        $CpHDIR/scripts/make_ndx $GROin TMP_CpHMD.ndx $GroDIR $CpHDIR/scripts/
    fi
    #
    #
    # Make relaxation .mdp gromacs file
    # Remove constraints & change NPT to NVT
    sed "s/\(nsteps *= \).*\$/\1 $RelaxSteps/
         s/\(Pcoupl *= \).*\$/\1 No/" \
        ${runname}.mdp > TMP_relax.mdp
    #
    echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
  Y Y Y\n" >> TMP_relax.mdp
    #
    # Make effective .mdp gromacs file
    sed "s/\(nsteps *= \).*\$/\1 $EffectiveSteps/" \
        ${runname}.mdp > TMP_aux2.mdp
    #
    # Modify Effective MDP file in case of using Position Restraints
    if [[ -f $PosRe && $PosRe != *(" ") ]]; then
	sed "s/\(define *= \).*\$/\1 -DPOSRES/" \
            TMP_aux2.mdp > TMP_effective.mdp
    else
        if [ ! -f $PosRe ]; then message W \
        "File $PosRe is missing. Position Restraints will not be used"; fi
        mv TMP_aux2.mdp TMP_effective.mdp
    fi
    #
    # Make initial .tpr file + TMP_processed.top
    if [[ $nchains -gt 1 ]]; then
    # use a whole gro
    "$GroDIR" grompp -f TMP_effective.mdp -po TMP_effective_out.mdp \
        -c  $GROwhole -p $TOPin -pp TMP_processed.top -n TMP_CpHMD.ndx \
        -o TMP_CpHMD.tpr -maxwarn 1000 -quiet
    else
    # use a normal gro
    "$GroDIR" grompp -f TMP_effective.mdp -po TMP_effective_out.mdp \
        -c  $GROin -p $TOPin -pp TMP_processed.top -n TMP_CpHMD.ndx \
        -o TMP_CpHMD.tpr -maxwarn 1000 -quiet
    fi
    rm -f TMP_aux*.mdp
    #
    # Build pdb of protein to be used with pdbswitch in a step 
    # of the topology updating (multi-chain) #SC 30-11-2011
    if [[ $nchains -gt 1 ]]; then
        #Here topology (-s) needs to be a pdb file so that the chain
        #information does not get lost
        echo "Protein" | "$GroDIR" trjconv -f $PDBin -s $PDBin\
            -o TMP_protein.pdb -n TMP_CpHMD.ndx -quiet
        if [ $RTThreshold == 0 ]; then 
            mv TMP_protein.pdb TMP_CpHMD_red.pdb
        fi
    else
    # Build gro of protein to be used with groswitch in a step 
    # of the topology updating (single chain) #SC 30-11-2011
        echo "Protein" | "$GroDIR" trjconv -f $GROwhole -s TMP_CpHMD.tpr \
            -o TMP_protein.gro -n TMP_CpHMD.ndx -quiet
        if [ $RTThreshold == 0 ]; then 
            cp -f $GROwhole TMP_CpHMD_red.gro # mv -> cp MM
        fi
    fi
    #
    # Make grid files .mgm and .ogm (in case you don't have your own)
    if [ $GridSize != 0 ]; then 
        echo -e "ON_GEOM_CENT 61 1.0\nON_CENT_OF_INTR 65 0.25" \
            > ${runname}.mgm
        echo -e "ON_GEOM_CENT $GridSize 1.0\nON_CENT_OF_INTR 65 0.25" \
            > ${runname}.ogm
    elif [[ $CpHModule != delphi ]]; then
        for f in ${runname}.{o,m}gm ;do
            if [ ! -f $f ];then message E \
            "File $f is missing!!!... Program will crash"; fi
        done    
    fi

    if [[ $CpHModule == delphi ]]
    then
        if [[ ! -f $GROwhole && $nchains == 1 ]]; then message E \
            "File $GROwhole is missing!!!... Program will crash"; fi
        #
        # Prepare input delphi parameters:
        if [ -f ${runname}.pbp ]; then
            mv ${runname}.pbp TMP_aux.pbp
            sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
                 s/\(export epsin*=\).*\$/\1$epsin/" \
                     TMP_aux.pbp > DELPHI.pbp
        fi
        #
        # Correct TMP_CpHMD.top with NTR and CTR
        #
        if egrep NT ${runname}.sites; then 

            awk '
            BEGIN{
                top=ARGV[1]; 
                while (getline < top)
                if ($5=="H3") ntr[++j]=$3;
                close(top)}; 
            END{
                while (getline < top) {
                write=1; 
                  for (i in ntr) {
                    if ($3==ntr[i] && $5=="N") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="H1") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="H2") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="H3") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="CA") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="C") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="O") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                  };
                  if (write) print $0}
            }' TMP_CpHMD.top > TMP_aux.top
        else
            mv TMP_CpHMD.top TMP_aux.top
        fi
        if egrep CT ${runname}.sites; then 
            awk '
            BEGIN{
                top=ARGV[1]; 
                while (getline < top)
                if ($5=="O1") ctr[++j]=$3;
                close(top)}; 
            END{
                while (getline < top) {
                write=1; 
                  for (i in ctr) {
                    if ($3==ctr[i] && $5=="C") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                    if ($3==ctr[i] && $5=="O1") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                    if ($3==ctr[i] && $5=="O2") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                    if ($3==ctr[i] && $5=="HO11") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                    if ($3==ctr[i] && $5=="HO12") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                    if ($3==ctr[i] && $5=="HO21") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                    if ($3==ctr[i] && $5=="HO22") 
                       {print substr($0,1,27),"CTR", substr($0,33); write=0};
                  };
                  if (write) print $0}
            }' TMP_aux.top > TMP_CpHMD_charge.top
        else
            mv TMP_aux.top TMP_CpHMD_charge.top
        fi
        # Generate charges database for Delphi
        "$DelphiDir"/gen_charge.awk ${runname}.sites \
            TMP_CpHMD_charge.top TMP_delphi.crg
    fi
}

run_PBMC()
{
    # Preparing files for pqr creation.
    # Do the same for dendrimer block
    if [[ $CpHModule != delphi ]]; then
    # Making protein whole removing PBC:
        echo "Protein" | "$GroDIR" trjconv -f TMP_effective.gro \
            -o TMP_aux1.gro -s TMP_CpHMD.tpr -n TMP_CpHMD.ndx -pbc mol -quiet
    #
    # The conversion can be obtained with MAKEPQR
    # NOTE: $TOPin cannot be processed top 
    # (moleculetype for the solvent cannot exist) SC 13-12-2011
	# TOPin=nasp-clys.top
	# MToolsDIR="/programs/meadTools2.0.1"
	# CpHDIR="/home/jsequeira/Projects/04_ASICs/FF/Amber14SBpH.ff"
	# ffID=Amber14SBpH
		
        # "$MToolsDIR"/makepqr W 2RT "$CpHDIR"/top/ff${ffID}nb.itp $TOPin \
        #    TMP_aux1.gro > TMP_aux.pqr

        "$MToolsDIR"/makepqr W 2RT ff${ffID}nb.itp $TOPin \
            TMP_aux1.gro > TMP_aux.pqr
    else
        if [ ${PBdim} -eq 2 ]; then
        #
        # Multi-step centering procedure:
        # 1- center with the tail of lipid 1
            echo -e "Onetail\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        # 
        # 2- center with all tails of the monolayer that has lipid 1
            echo -e "Monotail\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux1.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux2.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        # 
        # 3- center with all tails of all lipids in the system
            echo -e "Bitail\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux2.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux3.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        #
        # 4- center with all lipids in the system
            echo -e "Protein\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux3.gro -s TMP_CpHMD.tpr \
                -o TMP_${runname}.gro -n TMP_CpHMD.ndx \
                -center -pbc atom
        else
	    if [[ $nchains -gt 1 ]]; then
		#
		# In case multiple objects, use fixgro!
		echo -e "Protein\nProtein" | "$GroDIR" trjconv \
						       -f TMP_effective.gro -s TMP_CpHMD.tpr \
						       -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
						       -pbc mol -center
		#
		"$CpHDIR"/scripts/fixgro TMP_effective_aux1.gro \
			 ${runname}.fixgro > TMP_${runname}.gro
	    else
		# Making protein whole removing PBC:
		echo "Protein" | "$GroDIR" trjconv -f TMP_effective.gro \
					   -o TMP_${runname}.gro -s TMP_CpHMD.tpr \
					   -n TMP_CpHMD.ndx -pbc mol -quiet
            fi
            #
        # DEBUG - store PDB trajectory format for later analysis
        #
        #if [ $1 = "std" ]; then
        #    "$GroDIR" editconf -f TMP_${runname}.gro \
        #           -o TMP_aux_debug.pdb
        #    cat TMP_aux_debug.pdb >> ${blockname}_PQR.pdb
        #fi
        
	fi
    fi
    #

    #
    # Using the whole set of sites
    if [ $1 = "red" ]; then 
        cp -f ${runname}-all.sites ${runname}.sites
    fi 
    #
    if [[ $CpHModule != delphi ]]; then
    #
    # Offset correction
        if [ $offset = 0 ]; then
            pqr_file=TMP_aux.pqr
            all_sites_file=${runname}-all.sites
        else
        # This stmodels does not work with Cardiolipin (check stmodels_CL)
            "$MToolsDIR"/stmodels $offset TMP_aux.pqr ${runname}.sites
            mv ${runname}_stmod.sites ${runname}.sites
            if [[ $1 != red && $RTThreshold != 0 ]]; then
        # Make complete .sites and corresponding .pqr with offset
        # if reduced titration is used     
                "$MToolsDIR"/stmodels $offset TMP_aux.pqr ${runname}-all.sites
                all_sites_file=${runname}-all_stmod.sites
            fi
            pqr_file=TMP_aux_stmod.pqr
        fi
    #
    # Set sites excluded in reduced titration with right charges
        if [[ $1 != red && $RTThreshold != 0 ]]; then 
            "$MToolsDIR"/statepqr s=TMP_statepqr.in $pqr_file \
                $all_sites_file > TMP_${runname}.pqr
        else
            cp -f $pqr_file TMP_${runname}.pqr
        fi
    #
    # Set our pqr in the charged state
        "$MToolsDIR"/statepqr r=c TMP_${runname}.pqr ${runname}.sites \
            > ${runname}.pqr
    #
    # Run MeadT and Multiflex:
        "$MToolsDIR"/meadT -n $nCPU -s "$MToolsDIR" -m "$MeadDIR" \
            -b 250 -blab3 -epsin $epsin -epssol $epssol \
            -ionicstr $ionicstr  -T $temp \
            ${runname} >TMP_mead_out 2>TMP_mead_error
    else
    # Make link for DELPHI.pbp
    #   ln -s ${runname}.pbp DELPHI.pbp
    #
    # Run delphiT and Delphi: # MM 09/01/2012
	#DelphiDir="/home/jsequeira/Projects/04_ASICs/FF/DelphiTools_v2.1_AMBER/"
	
        "$DelphiDir"/delphiT $nCPU ${runname} ${PBdim} \
            >TMP_delphi.out 2>TMP_delphi.err
    fi
    #
    # Second part of the correction: We change back to the original a.a. 
    # res. numbers. In the delphi module the changes are done by delphiT.
    if [[ $offset != 0 && $CpHModule != delphi ]]; then
        mv ${runname}.pkcrg TMP_aux.pkcrg
        awk -v off=$offset '{match($0,/(^.+-)([0-9]+)$/,a);print a[1] a[2]-off*(1+($3~/^NT/)+2*($3~/^CT/))}' TMP_aux.pkcrg > ${runname}.pkcrg 
        mv ${runname}.sites TMP_aux.sites
        awk -v off=$offset '{print $1-off*(1+($2~/^NT/)+2*($2~/^CT/)), $2,$3,$4,$5,$6,$7}' TMP_aux.sites > ${runname}.sites  #MM 09/11/2010
    fi
    #
    if [[ $CpHModule == redox ]]; then
        # Mark Redox sites
        mv ${runname}.pkcrg TMP_aux2.pkcrg
        sed 's/HEMall-\([0-9]\+\)/\0    R/g' TMP_aux2.pkcrg > ${runname}.pkcrg
    fi
    #
    # Run cconvert
    "$MToolsDIR"/convert ${runname}.pkcrg ${runname}.g $temp > ${runname}.dat
    #
    # Run PETIT
    echo "$pH,$pH,1" -E "$pot,$pot,1" -T $temp -c 2 -r $seed -q 1000 100000
    "$PetitDIR"/petit -H "$pH,$pH,1" -E "$pot,$pot,1" -T $temp \
        -c 2 -r $seed -q 1000 100000 <${runname}.dat \
        >TMP_MCarlo_$1.out 2>TMP_MCarlo_$1.err 
    #   
    # Reduced Titration calculations to determine the sites to use
    # in the next cycles
    if [ $1 = "red" ]; then 
        # Clear of the .sites file before the new one is created 
        # (need to do this in order to work in case of empty .sites 
        # from reduced titration)
        rm -f ${runname}.sites ; touch ${runname}.sites
        awk -v t=$RTThreshold -v Allsites=${runname}-all.sites \
            -v Redsites=${runname}.sites '
        BEGIN{
          # Read petit output (all-site CE/MC):
          while (getline < "TMP_MCarlo_red.out")
          {
          # Read occR (from all-site CE/MC):
          if ($0 ~ /^f/)
          {
            nsites = NF - 1 ;
            for(i = 2 ; i <= NF ; i++) occR[i-1] = $i ;
          }
          # Read moccR (from all-site CE/MC):
          if ($0 ~ /^\./ && $0 !~ /tot/)
          {
            moccR[$4+1] = $5 ;
            m = 0 ;
            for (i = 6 ; i <= NF ; i++) 
              if($i > m)
              {
                m = $i ;
                # state with maximum population:
                maxstate[$4+1] = i - 6 ;
                # maximum population of that state:
                maxocc[$4+1] = $i ;
             }
            # If maximum population is above threshold, make switched=0,
            # indicating that site would be fixed (in state maxstate)
            # during next MD segments.
            switched[$4+1] = (maxocc[$4+1] > 1-t ? 0 : 1) ;
          }
          # Make input (first part of) for groswitch:
          if ($0 !~ /^f/) print $0 > "TMP_MCarlo_mod.out" ;
          }
          close("TMP_MCarlo_red.out") ;
        
          # Make input (second part of) for groswitch (among other things):
          printf("f ") > "TMP_MCarlo_mod.out" ;
          for (i = 1 ; i <= nsites ; i++)
          {
            printf ("%d ", maxstate[i]) > "TMP_MCarlo_mod.out" ;
  
            #Make input file for state_pqr with most abundant state
            printf ("%d\n", maxstate[i]) > "TMP_statepqr.in"

            printf ("%s ", switched[i] == 0 ? maxstate[i] : "-") > "TMP_template_occ" ;
            printf ("%s ", switched[i] == 0 ? moccR[i] : "-") > "TMP_template_mocc" ;
            printf ("%.6f ", moccR[i]) >> "TMP_CpHMD.mocc_red" ;
            printf ("%d ", occR[i]) >> "TMP_CpHMD.occ_red" ;
            # Modded from maxocc to maxstate
            printf ("%d ", maxstate[i]) >> "TMP_CpHMD.occ_red_mod" ;
          }
          printf("\n") > "TMP_MCarlo_mod.out" ;
          printf("\n") > "TMP_template_occ" ;
          printf("\n") > "TMP_template_mocc" ;
          printf("\n") >> "TMP_CpHMD.mocc_red" ;
          printf("\n") >> "TMP_CpHMD.occ_red" ;
          printf("\n") >> "TMP_CpHMD.occ_red_mod" ;
        
          # Make .sites:
          n = 1 ;
          while (getline < Allsites)
          {
            if (switched[n] == 1) print $0 > Redsites ;
            n++ ;
          }
          close(Allsites) ;
        }'
        #
        # Making Log for .sites
        echo "This is the .sites file at Cycle = $Cycle" >>  TMP_CpHMD.sites
        cat ${runname}.sites >> TMP_CpHMD.sites
        #
        # Writing corrected GRO file for build_topology function 
        # (1 - use tautomers) (0 - do not use tautomers)
        if [[ $nchains -gt 1 ]]; then
            "$CpHDIR"/scripts/pdbswitch $gstaut TMP_MCarlo_mod.out \
                TMP_protein.pdb > TMP_CpHMD_red.pdb #SC 11-11-2011
        elif [[ $CpHModule != delphi ]]; then
            "$CpHDIR"/scripts/groswitch $gstaut TMP_MCarlo_mod.out \
                TMP_protein.gro > TMP_CpHMD_red.gro
        else
            "$CpHDIR"/scripts/groswitch $gstaut TMP_MCarlo_mod.out \
                $GROwhole > TMP_CpHMD_red.gro
        fi
    fi
    #
    # Removing PB related auxiliary files
    rm -f ${runname}.{summ,pqr*,dat,pkcrg,g,pkint,out} \
        *.potat TMP_aux* TMP_mead_out TMP_mead_error TMP_MCarlo_red*

}

write_fractions ()
{
    # Write the occupation files
    awk -v t=$RTThreshold '
    BEGIN{
    
      # Read petit output (CE/MC):
      while (getline < "TMP_MCarlo_std.out")
      {
        # Read occ (from CE/MC):
        if ($0 ~ /^f/)
        {
          nsites = NF - 1 ;
          for(i = 2 ; i <= NF ; i++) tocc[i-1] = occ[i-1] = $i ;
        }
        # Read mocc (from CE/MC):
        if ($0 ~ /^\./ && $0 !~ /tot/) tmocc[$4+1] = mocc[$4+1] = $5 ;
      }
      close("TMP_MCarlo_std.out") ;
    
      # If t != 0 read templates and override tocc and tmocc:
      if (t != 0)
      {
        getline < "TMP_template_occ" ;
        nsites = split($0, tocc) ;
        getline < "TMP_template_mocc" ;
        split($0, tmocc) ;
        c = 0 ;
      }
    
      # Write tocc and tmocc:
      for (i = 1 ; i <= nsites ; i++)
      {
        # Substitute "-" with corresponding occ and mocc entries:
        if (t != 0 && tocc[i] == "-")
        {
          c++ ;
          tocc[i] = occ[c] ;
          tmocc[i] = mocc[c] ;
        }
        printf ("%d ", tocc[i]) >> "TMP_CpHMD.occ" ;
        printf ("%f ", tmocc[i]) >> "TMP_CpHMD.mocc" ;
      }
      printf("\n") >> "TMP_CpHMD.occ" ;
      printf("\n") >> "TMP_CpHMD.mocc" ;
    
    }'

}

build_topology ()
{
    
 # Set termini if they are titrating
    
    if [[ $nchains -gt 1 ]]; then
#################################################################
#SCAMPOS: New ter reading to account for several chains
#14-11-2011
   
        n=0
        if [[ ${NTindex[@]} != "" ]]; then
            for i in ${NTindex[@]}
            do
                NTer[$i]=`tail -n -1 TMP_CpHMD.occ | awk -v i=${NTpos[$n]} '{print $i}'`
                ((++n))
            done
        fi

        n=0
        if [[ ${CTindex[@]} != "" ]]; then
            for i in ${CTindex[@]}
            do
                CTer[$i]=`tail -n -1 TMP_CpHMD.occ | awk -v i=${CTpos[$n]} '{print $i}'`
                ((++n))
            done
        fi
#################################################################
    else
        if  egrep NT ${runname}-all.sites; then
            nter=`tail -n -1 TMP_CpHMD.occ | awk '{printf "%u\n", $1}'`
        fi
        if  egrep CT ${runname}-all.sites; then 
            cter=`tail -n -1 TMP_CpHMD.occ | awk '{printf "%u\n", $NF}'`
        fi
    fi


    # Run PDB2GMX to generate the new topology
    if [[ $nchains -gt 1 ]]; then
        # Define protonation states for the sites needed. The first flag of 
        # pdbswitch defines the use of tautomers
        "$CpHDIR"/scripts/pdbswitch $gstaut TMP_MCarlo_std.out \
                 TMP_CpHMD_red.pdb > TMP_aux1.pdb #SC 11-11-2011
        #
	# Print to the .err file the N- and C-terminus protonation states used
        (for ((i=0; i<$nchains; i++)); do 
             echo -e "${NTer[$i]}\n${CTer[$i]}"; done)
	#
        (for ((i=0; i<$nchains; i++)); do 
             echo -e "${NTer[$i]}\n${CTer[$i]}"; done) | \
                 "$GroDIR" pdb2gmx -f TMP_aux1.pdb \
			   -p TMP_aux1.top -o TMP_aux3.gro -i TMP_posre.itp \
			   -ignh -ter -ff $ffID  -water $water \
			   -quiet -merge all
    else
        # Define protonation states for the sites needed. 
        # The first flag of groswitch defines the use of tautomers
        "$CpHDIR"/scripts/groswitch $gstaut TMP_MCarlo_std.out \
                 TMP_CpHMD_red.gro > TMP_aux1.gro #SC 21-11-2011
        #
        # Removing header from GRO file (growing header can be troublesome)
        echo "GRO file from Constant pH MD" > TMP_aux2.gro 
        tail -n +2 TMP_aux1.gro >> TMP_aux2.gro 
        #
        echo -e "$nter\n$cter"  | "$GroDIR" pdb2gmx -f TMP_aux2.gro \
					    -p TMP_aux1.top -o TMP_aux3.gro -i TMP_posre.itp \
					    -ignh -ter -ff $ffID -water $water -quiet -merge all
    fi
    #################################################################################################################################################################################
    #  
    # DEBUG SC 30-11-2011  
    #   echo "Gro file from: `printf "% 3d\n" $Cycle`" >> ${blockname}_ALLGRO.gro
    #   cat TMP_aux3.gro >> ${blockname}_ALLGRO.gro
        #
        # Correct the Protein group name
        sed 's/Protein_chain_A/Protein        /g' TMP_aux1.top > TMP_CpHMD.top
    ######################################################################
    #
    #
    awk -v s="$TOPSTRING" '/moleculetype/, $0 ~ s' TMP_CpHMD.top \
        > TMP_CpHMD.toptit #SC 28-11-2011
    #
    cat TMP_CpHMD.topbegin TMP_CpHMD.toptit TMP_CpHMD.topend | \
        awk '{if($1~/Protein/){print substr($1,1,7), $2, $3}else{print $0}}' \
        > TMP_CpHMD.top  #SC 28-11-2011 MM 17-01-2012

    if [[ $RULEin != *(" ") ]]; then
        # Fix the topology with the assigned rules
        mv TMP_CpHMD.top TMP_aux2.top
        $CpHDIR/scripts/fix_topology TMP_aux2.top \
            $RULEin > TMP_CpHMD.top # SC 25-11-2011
#        nterm=`awk '$4=="NH2" && $6=="NH2" {print $3}' TMP_CpHMD_1.top`
#	if [ "$nterm" -ge 0 ]
#	then
#	    sed "s/51    53    55    61    63     1//g" TMP_CpHMD_1.top > TMP_CpHMD.top
#	fi
	
                                    #several rules need a solution; 
                                    #If RULEin is defined as an array, will it work ${RULEin[@]}? 
                                    #SC 13-12-2011
    fi

    # Fix the topology using fix_dendrimer_top. LCSF 2011-10-03
    if [[ $CpHModule == dendrimer ]]; then
        mv TMP_CpHMD.top TMP_aux2.top
        $CpHDIR/scripts/fix_dendrimer_top \
            $CpHDIR/scripts/${RULEdendr} TMP_aux2.top > TMP_CpHMD.top
    fi

    # Correct POSRE file in case of using Position Restraints
    if [[ -f $PosRe && $PosRe != *(" ") ]]; then 
        cp -f $PosRe TMP_posre.itp
    fi
    #
    # debugging the final topologies MM
    # cp TMP_CpHMD.top topology_debug_${Cycle}.top
    # Housekeeping
    rm -f TMP_aux* \#*

}

run_dynamics ()
{
    "$GroDIR" grompp -f TMP_$1.mdp -po TMP_$1_out.mdp \
        -c TMP_$2.gro -p TMP_CpHMD.top -pp TMP_processed.top \
        -n TMP_CpHMD.ndx -o TMP_$1.tpr -maxwarn 1000 -quiet

    $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
        -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr -rcon $Rcon -quiet \
        -nice 19 # SC-14-12-2011

    rm -f \#*
}

run_relaxation () #SC 28-11-2011
{
    #Solvent relaxation
    run_dynamics relax effective

    #Prepare input GRO for dynamics
    awk -v s=$SOL1st '$1 ~ s {exit};{print $0}' TMP_effective.gro > TMP_aux.gro
    awk -v s=$SOL1st '$1 ~ s {a=1};a'  TMP_relax.gro >> TMP_aux.gro
    
#    mv -f TMP_relax.gro TMP_relax_DEBUG.gro
    mv -f TMP_aux.gro TMP_relax.gro

}


data_append ()
{
    if [ $Cycle -eq $InitCycle ]; then
        InitTime=`echo $sim_time-$WriteTime | bc -l`
        initxtc=""
        initedr=""
    else
        initxtc="TMP_CpHMD.xtc"
        initedr="TMP_CpHMD.edr"
    fi
    # Append .edr files
    echo -e "$InitTime\n`echo $sim_time-$WriteTime | bc -l`" | \
        "$GroDIR" eneconv  -o TMP_CpHMD.edr -f $initedr \
                  TMP_effective.edr -settime -quiet
    #
    # Append .xtc files
    echo -e "$InitTime\n`echo $sim_time-$WriteTime | bc -l`" | \
        "$GroDIR" trjcat -f  $initxtc TMP_effective.xtc  \
                  -o TMP_CpHMD.xtc -settime
    #
    # Append and backup remaining files
    if [ -f TMP_effective.log -o -f TMP_effective0.log ]
    then
        cat TMP_effective*.log >> TMP_CpHMD.log
    fi
    cp -f TMP_effective.gro TMP_CpHMD.gro
    rm -f TMP_aux* \#* TMP_effective*.log
}

clean_up ()
{
    rm -rf ${runname}.{summ,g,pkint,pkcrg,dat,pqr,out} state*\
        *.st aminoacids.dat FF.dat specbond.dat ff* traj.trr \
        ${runname}_cpu* TMP_{statepqr.in,${runname}.pqr,effective,relax,CpHMD,aux,mead,MCarlo,template,posre,processed}*
}
