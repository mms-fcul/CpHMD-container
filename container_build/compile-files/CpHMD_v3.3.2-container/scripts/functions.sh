# This file contains all the functions called by CpH-MD.sh script
#It also contains the definition of some constants.
#
# Please don't delete it
#
correct_variables ()
{                   
    #removed offset value since CpHMD sources this value
    GridSize=`awk -v i=$GridSize 'BEGIN{print (i==0) ? 0 : i}'`
    RTThreshold=`awk -v i=$RTThreshold 'BEGIN{print (i==0) ? 0 : i}'`
}

build_forcefield ()
{
    if [ ! -d $ffDIR ]; then 
          message E  "Force field directory not defined in the .settings!!!... Program will crash"
    fi
    
    # Call the modified forcefield:
    cp -rdf "$ffDIR"/${ffID}.ff .;     
}

check_files ()
{
    # Check for other necessary files (filenames must be respected)
    for f in ${runname}.mdp $GROin $TOPin $GroDIR $NDXin
      do
      if [ ! -f $f ]; then
          message E  "File $f is missing!!!... Program will crash"
      fi
    done

    # Check topology to see if the name of the molecule is set to protein, otherwise it will crash
    mol_name=`awk '/nrexcl/{getline; print $1}' $TOPin`
    if [[ ! $mol_name =~ "Protein" ]] ; then
	message E  "Topology does not have the molecule set as Protein. Program will crash"
    fi

    #s
    # Check directories
    for d in $PetitDIR $StDIR $DelphiDir
      do
      if [ ! -d $d ]; then 
          message E  "Directory $d does not exist!!!... Program will crash"
      fi
    done    

    if [[ $PBdim -ne 0 && $PBdim -ne 2 ]]; then
        message E  "Number of dimensions ($PBdim) is not supported in ${CpHModule} module!!!... Program will crash"
    fi
    #
    # Detect cutoff-scheme to be used
    scheme=`awk '/cutoff-scheme/ {print $3}' ${runname}.mdp `
    if [[ $scheme == "verlet" ]]; then
	message W "Scheme used is verlet."
    elif [[ $scheme == "group" ]]; then
        message W "Scheme used is Group based."
	## Detect if the gromacs version supports group based
	#
	$GroDIR help  > ./gro-ver.txt 2>&1
	if [[ ! -z `awk '/error/' gro-ver.txt` ]]; then
	    message E "An error as occured with the gromacs version"
	fi	
	version=`awk '/GROMACS:/ {print $(NF)}' gro-ver.txt`
	if [[ $version != "5.1.5" ]]; then
	    message E  "Cut-off scheme ${scheme} is not supported in ${GroDIR} module!!!... Program will crash"
	fi
    else
        message W "Scheme is not group nor verlet... Performing method as group based cut off"
    fi
    
    ###########################################################################
    ## Check if gromacs version is not 2020.2 since it has freezedims bugged ##
    ###########################################################################
    $GroDIR help  > ./gro-ver.txt 2>&1 
    version=`awk '/GROMACS:/ {print $(NF)}' gro-ver.txt`
    if [[ $version == "2020.2" ]]; then
	message E  "Gromacs ${version} is not supported for CpHMD since freezedims is bugged !!... Program will crash"
    fi

    ##############################
    ### Check for Cuda Support ###
    ##############################
    cuda=`awk '/GPU support:/ {print $(NF)}' gro-ver.txt`
    if [[ $GPU == 1 && $cuda == "disabled" ]] ; then
	message E  "Gromacs ${version} was not compiled with GPU support. Please correct the gromacs compilation or remove the GPU flag"
    fi
    ##evaluate if mdp contains only one energy group ##
    if [[ `awk -F "=" '/energygrps/ {print $2}' ${runname}.mdp | awk '{print NF}'` > 1 && `awk -F "=" '/energygrps/ {print $2}' ${runname}.mdp ` != "System" && $GPU == 1 ]] ; then
	message E  "ERROR: Use of GPU and multiple energy groups is not supported, simulations would be downgraded to CPU only. Please either use only 1 energy groups or remove GPU support."
    fi
    
}

make_auxiliary_files ()
{
    # Rename original files
    cp -f $TOPin TMP_CpHMD.top
    cp -f $NDXin TMP_CpHMD.ndx
    cp -f $GROin TMP_effective.gro
    
  
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

    # Make relaxation .mdp gromacs file
    # Remove constraints & change NPT to NVT
    sed "s/\(nsteps *= \).*\$/\1 $RelaxSteps/ 
         s/\(Pcoupl *= \).*\$/\1 No/" \
        ${runname}.mdp > TMP_relax.mdp
    #
    #### Adition on 13/01/2021 to run in either group and verlet cutoff based runs ####
    ### Add Different protocol depending on the cut-off scheme ###
    # Check which cut-off is being used #
    scheme=`awk '/cutoff-scheme/ {print $3}' ${runname}.mdp `
    if [[ $scheme == "verlet" ]]
    then
        echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
        Y Y Y\n" >> TMP_relax.mdp
    elif [[ $scheme == "group" ]]
    then
        echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
        Y Y Y\n\nenergygrp_excl      =  Protein Protein" >> TMP_relax.mdp
    else
        message W "Scheme not group nor verlet... Performing method as group based cut off"
        echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
        Y Y Y\n\nenergygrp_excl      =  Protein Protein" >> TMP_relax.mdp

    fi
    #
    # Make effective .mdp gromacs file
    sed  "s/\(nsteps *= \).*\$/\1 $EffectiveSteps/" ${runname}.mdp > TMP_effective.mdp
    #
    # Modify Effective MDP file in case of using Position Restraints
    if [[ -f $PosRe && $PosRe != *(" ") ]]; then 
        sed -i "s/\(define *= \).*\$/\1 -DPOSRES/" TMP_effective.mdp
    else
        if [ ! -f $PosRe ]; then message W \
        "File $PosRe is missing. Position Restraints will not be used"; fi
    fi
    #
    ### 24/04 changed -p $TOPin to -p TMP_CpHMD.top for the grompp
    ##to run using the folder's FF
    "$GroDIR" grompp -f TMP_effective.mdp -po TMP_effective_out.mdp \
        -c  $GROin -p TMP_CpHMD.top -pp TMP_processed.top -n TMP_CpHMD.ndx \
        -o TMP_CpHMD.tpr -maxwarn 1000 -quiet
    
    rm -f TMP_aux*.mdp

    #################################################################
    ### Getting the correct pbp file for PBMC                     ###
    #################################################################
    if [[ $PBdim == 0 ]]
    then
	sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
             s/\(export epsin*=\).*\$/\1$epsin/" \
		 "$DelphiDir"/DELPHI_Prot.pbp > DELPHI.pbp
    elif [[ $PBdim == 2 ]]
    then
	 sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
             s/\(export epsin*=\).*\$/\1$epsin/" \
		 "$DelphiDir"/DELPHI_Memb.pbp > DELPHI.pbp
    fi
    
}

make_sites () 
{
    echo "Protein" | "$GroDIR" editconf -f $GROin \
				   -o TMP_protein.pdb -n TMP_CpHMD.ndx -quiet
    
    
    "$CpHDIR"/scripts/make_sites $1 TMP_protein.pdb 
    #
    # Correct TMP_CpHMD.top with NTR and CTR
    #
    #### Previously ######
    ## Would only correct it if NT was in sites, now it will try to do it
    ## Always if it has the H3 ATOM

    if [[ $ffID ==  G54a7pH ]];then 
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

    elif [[ $ffID == Amber14SBpH ]] ;then 
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
		if ($3==ntr[i] && $5=="HA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="C") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="O") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_CpHMD.top > TMP_aux.top

    elif [[ $ffID == CHARMM36pH ]] ;then 
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
		if ($3==ntr[i] && $5=="HA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_CpHMD.top > TMP_aux.top
	
    fi
	
    #if egrep CT ${runname}.sites; then
    if [[ $ffID ==  G54a7pH ]];then 
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
	
    elif [[ $ffID == Amber14SBpH ]] ;then #added AMBER case with CTR names
	awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="OC1") ctr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ctr) {
	      	if ($3==ctr[i] && $5=="N") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="H") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="CA") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HA") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="C") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="OC1") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="OC2") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC11") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC12") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC21") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC22") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_aux.top > TMP_CpHMD_charge.top
    elif [[ $ffID == CHARMM36pH ]] ;then #added AMBER case with CTR names
	awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="OT1") ctr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ctr) {
                if ($3==ctr[i] && $5=="C") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="OT1") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="O2") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT11") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT12") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT21") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT22") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_aux.top > TMP_CpHMD_charge.top
	
    fi

    ##################### Adition for offset check ############################
    ## Define Offset and verify if it is able the system is able to compute  ##
    ###########################################################################
    ## negate the pattern match of declaring an array with the chainres ##
    
    last_res=`awk '/ATOM/ {print substr($0,23,4)}' TMP_protein.pdb | tail -n 1`
    export offset=`echo $last_res | awk '{print $1 +10}'`
    ###############  Debug addition #############
    message W "Current set up offset is $offset" 

    last_tit=`awk '{print $1}' ${runname}.sites | tail -n 1 `

    #############################
    if [ $(( (offset * 2) +500 + last_tit )) -ge 10000 ] ; then
	message E "System simulated will wield a PB/MC renumbering larger than 10 000 which will crash the PBMC cycle. Your system is currently unfeasible to calculate with CpHMD."
    fi

}

make_delphi_DB()
{
    echo "atom__res_radius_" > ./DataBaseT.siz
    echo "atom__resnumbc_charge_" > ./DataBaseT.crg
    for type in crg siz ; do
	## Insert the terminals in the database
	for ter in CT NT; do
	    case $type in
		crg)   awk -v t=$ter '$2~t {printf"%-6s%-9s%6.3f\n", $1,$2,$3}' $DelphiDir/DataBaseT_${ffID}.${type} >> ./DataBaseT.${type} ;;
		siz)   awk -v t=$ter '$2~t {printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' $DelphiDir/DataBaseT_${ffID}.${type} >> ./DataBaseT.${type} ;;
	    esac
	done
	## Insert other residues
	for db_res in `awk '/ATOM/ {print substr($0,18,3)}' TMP_protein.pdb | sort | uniq` ; do	    
	    if [[ ! -z `awk -v d=$db_res '$2==d {print}' "$DelphiDir"/DataBaseT_${ffID}.${type}` ]] ; then
		## If the third character of the residue name is a number (hence a CpHMD residue most likely) ##
		if [[ `echo $db_res | awk '{if (substr($1,3,1) ~ /^[0-9]/) {print 1}}' ` == 1 ]] ;then
		    ## Check for the residue in question to not be in the database already ##
		    if [[ -z `awk  -v d=${db_res} '$2==d {print}' ./DataBaseT.${type}` ]] ; then
			case $type in
			    crg)   awk -v d=${db_res} '$2~substr(d,0,2) && substr($2,3,1)~/[0-9]/ {printf"%-6s%-9s%6.3f\n", $1,$2,$3}' $DelphiDir/DataBaseT_${ffID}.${type} >> ./DataBaseT.${type} ;;
			    siz)   awk -v d=${db_res} '$2~substr(d,0,2) && substr($2,3,1)~/[0-9]/ {printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' $DelphiDir/DataBaseT_${ffID}.${type} >> ./DataBaseT.${type} ;;
			esac
		    fi
		else
		    case $type in
			crg)   awk -v d=${db_res} '$2~d {printf"%-6s%-9s%6.3f\n", $1,$2,$3}' $DelphiDir/DataBaseT_${ffID}.${type} >> ./DataBaseT.${type} ;;
			siz)   awk -v d=${db_res} '$2~d {printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' $DelphiDir/DataBaseT_${ffID}.${type} >> ./DataBaseT.${type} ;;
		    esac
		fi
	    else
		message E  "Residue identifier $db_res is not on Delphi Database_${ffID}.${type} Program will crash"
	    fi
	done
	## Other option to avoid duplicate lines in database is using an associative array in awk##
	# awk '!a[$0]++' ./DataBaseT.${type} >> ./DataBaseT.${type} but the naming for the databases on the for cycles
	# need to be renamed
	
	## check if size of database isn't too large ##
	if [[ `cat ./DataBaseT.${type} | wc -l` -ge 1000 ]] ; then
	    message E  "DataBaseT.${type} exceeds 1000 lines. Fortran will crash!!!"
	fi
    done

    ## make CRG_FILE ##
    awk 'NF==3 {printf"%-6s%-9s%6.3f\n", $1,$2,0.000}' ./DataBaseT.crg > CRG_FILE
    
    # Generate charges database for Delphi
    "$DelphiDir"/gen_charge.awk ${runname}.sites \
		TMP_CpHMD_charge.top TMP_delphi.crg

}

run_PBMC()
{
        if [ $memb == 1 ]; then
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
	    
	    if [[ $PBdim -eq 0 ]]
	    then
		echo -e "Protein" | ${GroDIR} trjconv -f TMP_${runname}.gro -s TMP_CpHMD.tpr \
				      -n TMP_CpHMD.ndx -o TMP_aux1.gro -pbc mol -quiet

		rm -f TMP_auxcenter.dat
		for i in `awk '{print $1}' ${runname}.sites` 
		do
		    awk -v i=$i ' substr($0,1,5)+0==i {print substr($0,21,24)}'  TMP_aux1.gro >> TMP_auxcenter.dat
		done
		
		protx=`awk '{x+=$1;n++}END{print x/n}' TMP_auxcenter.dat`
		proty=`awk '{x+=$2;n++}END{print x/n}' TMP_auxcenter.dat`
		
		halfsizex=`tail -n 1 TMP_aux1.gro |awk '{print $1/2}'`
		halfsizey=`tail -n 1 TMP_aux1.gro |awk '{print $2/2}'`
		
		XCoor=`echo $protx $halfsizex | awk '{print $2-$1}'`
		YCoor=`echo $proty $halfsizey | awk '{print $2-$1}'`

		${GroDIR} editconf -f TMP_aux1.gro -o TMP_aux2.gro \
		  -translate ${XCoor} ${YCoor} 0 -quiet
    
		echo -e "Protein" | ${GroDIR} trjconv -f TMP_aux2.gro -s TMP_CpHMD.tpr \
				      -n TMP_CpHMD.ndx -o TMP_${runname}.gro -pbc atom -quiet

	    fi
        else
	    if [[ $multiple_solutes -eq 1 || $include_itp -eq 1 ]] ; then 
		#
		# In case multiple objects, use fixgro!
		# A fixgro variable was created on pHmdp for rare cases
		# with multiple objects but a single chain
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

	    ##DEBUG - store PDB trajectory format for later analysis
            #
            #"$GroDIR" editconf -f TMP_${runname}.gro \
            #       -o TMP_aux_debug.pdb
            #cat TMP_aux_debug.pdb >> ${blockname}_PQR.pdb
        
	fi
    #
    # Run delphiT and Delphi: # MM 09/01/2012
        "$DelphiDir"/delphiT $nCPU ${runname} ${PBdim} \
		    >TMP_delphi.out 2>TMP_delphi.err

    ### Check for fortran STOP within the output ###
	if grep -q "FORTRAN STOP" TMP_delphi.err
	then
	    message E "Fortran Stop on delphi step, something went wrong! on PB/MC cycle."
	fi
    # Run cconvert
    "$DelphiDir"/convert ${runname}.pkcrg ${runname}.g $temp > ${runname}.dat
    #
    # Run PETIT
    echo "$pH,$pH,1" -E "$pot,$pot,1" -T $temp -c 2 -r $seed -q 1000 100000
    "$PetitDIR"/petit -H "$pH,$pH,1" -E "$pot,$pot,1" -T $temp \
        -c 2 -r $seed -q 1000 100000 <${runname}.dat \
        >TMP_MCarlo.out 2>TMP_MCarlo.err 
    #
    #### ADD reduce titration code here! ####

    #########################################
    # Removing PB related auxiliary files

    rm -f ${runname}.{summ,pqr*,dat,pkcrg,g,pkint,out} \
        *.potat TMP_aux* TMP_mead_out TMP_mead_error 

}

write_fractions ()
{
    # Write the occupation files
    awk '
    BEGIN{
    
      # Read petit output (CE/MC):
      while (getline < "TMP_MCarlo.out")
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
      close("TMP_MCarlo.out") ;
    
      # Write tocc and tmocc:
      for (i = 1 ; i <= nsites ; i++)
      {
        printf ("%d ", tocc[i]) >> "TMP_CpHMD.occ" ;
        printf ("%f ", tmocc[i]) >> "TMP_CpHMD.mocc" ;
      }
      printf("\n") >> "TMP_CpHMD.occ" ;
      printf("\n") >> "TMP_CpHMD.mocc" ;
    
    }'

}

update_topology ()
{
    mv TMP_CpHMD.top TMP_CpHMD-pre.top
    
    #"$CpHDIR"/scripts/groswitch 1 TMP_MCarlo.out \
    #             TMP_effective-pre.gro > TMP_effective.gro
    
    ### First we create a state_file similar to the output of MEAD ###
    
    # Getting res number and res name #
    awk '{print $1, substr($2,1,3)}' ${runname}.sites > aux_sites.sites

    # Getting the prot states from the montecarlo output #

    awk '$1=="f" {for (i=2; i<=NF ; i++) print $i}' TMP_MCarlo.out > tmp_states.out


    paste tmp_states.out aux_sites.sites > states_file.dat


    #### Run update topology script ####

    "$CpHDIR"/scripts/update-top states_file.dat ${ffID}.ff/protstates.dic TMP_CpHMD-pre.top > TMP_CpHMD.top
    
    #rm -rf tmp_states.out aux_sites.sites
}


run_dynamics ()
{
    "$GroDIR" grompp -f TMP_$1.mdp -po TMP_$1_out.mdp \
        -c TMP_$2.gro -p TMP_CpHMD.top -pp TMP_processed.top \
        -n TMP_CpHMD.ndx -o TMP_$1.tpr -maxwarn 1000 -quiet

    $mdrun -s TMP_$1.tpr -x TMP_$1.xtc -c TMP_$1.gro \
        -e TMP_$1.edr -g TMP_$1.log -o TMP_$1.trr \
        -nice 19 # SC-14-12-2011

    rm -f \#*
}

run_relaxation () #SC 28-11-2011
{
    #Solvent relaxation
    run_dynamics relax effective

    #Prepare input GRO for dynamics
    # Solvent marker is hardcoded in this version
    SOL1st="SOL"
    awk -v s=$SOL1st '$1 ~ s {exit};{print $0}' TMP_effective.gro > TMP_aux.gro
    awk -v s=$SOL1st '$1 ~ s {a=1};a'  TMP_relax.gro >> TMP_aux.gro
    
#    mv -f TMP_relax.gro TMP_relax_DEBUG.gro
    mv -f TMP_aux.gro TMP_relax.gro

}

data_append ()
{
    if [ $Cycle -eq 1 ]; then
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
    if [ -f pullx.xvg ]; then
	for e in x f ; do
	    awk -v t=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)" | bc -l` \
		'!/^\#/ && !/^\@/ && $1!=0.0000{print $1+t,$2}' pull$e.xvg \
		>> TMP_CpHMD_pull$e.xvg
	done
    fi
    rm -f TMP_aux* \#* TMP_effective*.log pull?.xvg
}

clean_up ()
{    
    rm -rf ${runname}.{summ,g,pkint,pkcrg,dat,pqr,out} state*\
        *.st aminoacids.dat FF.dat specbond.dat ff* traj.trr \
        ${runname}_cpu* TMP_{statepqr.in,${runname}.pqr,effective,relax,CpHMD,aux,mead,MCarlo,template,posre,processed,pqr}*

   
}
