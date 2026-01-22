#!/usr/bin/awk -f

##############################################################
#                                                            #
#         Vitor H. Teixeira, Lisboa, 2013.03.26              #
#                                                            #
#     This program generates a file with a variable number   #
#  of columns and a number of lines equal to the number of   #
#  atoms. The first column has the reference charge of the   #
#  "site" and the following columns correspond to the atomic #
#  charges in each tautomer. The number of columns may vary  #
#  along the file.                                           #
#                                                            #
##############################################################

BEGIN{
  cmd = "gen_charge.awk" ;
  usage = "Usage: "cmd" <sites> <top> <out> " ;
  if (ARGC != 4) error("Wrong number of arguments.\n" usage) ;

  filecheck(sites = ARGV[1]) ;
  filecheck(top = ARGV[2]) ;
  out = ARGV[3] ;


###  Reading sites file 
  while (getline < sites) {

    n++ ;
    Taut[n] = NF - 1 ;
    site_i[n] = $1 ;
    for (i = 1; i <= Taut[n]; i++) {

      site[n,i] = $(i+1) ;

##########  Remove below  ##########
#
#      if (site[n,i] ~ /^(NT|CT)/)
#        tempS[n,i] = substr(site[n,i],3,2) ;
#      else 
#
##########  Remove above  ##########
      tempS[n,i] = substr(site[n,i],1,2) ;


    }
  }
  close(sites) ;




###  Initializing number of atoms
  natoms = 0 ;

###  Reading TOP file (Only atoms block)
  while (getline < top) {

    if ($1 ~ /^(;|#)/) continue ;
    if (NF == 0) {

      if (moleculetype) moleculetype = 0 ;
      if (atom) protein = atom = 0 ;
    }
    if (moleculetype) {

      if ($1 ~ /^Protein/) protein = 1 ;
    }
    if (atom && protein) {

      natoms++ ;
      res_i[natoms] = $3 ;
      tempR[natoms] = substr($4,1,2) ;
      atom_n[natoms] = $5 ;
      charge[natoms] = $7 ;
    }
    if ($0 ~ /\[ moleculetype \]/) moleculetype = 1 ;
    if ($0 ~ /\[ atoms \]/) atom = 1 ;
  }
  close (top_file) ;


###  Reading st files
  for (s = 1; s <= n; s++) {

## Read/Check only the first file
    filecheck(st_file = site[s,1] ".st") ;
    getline < st_file ;     # Do not care about pKmod, yet

###  Get number of atoms to read charges later
    l = 0 ;
    Tot_p = Tot_d = 0 ;
    while (getline < st_file) {

      l++ ;
      Tot_p += $3 ;
      Tot_d += $4 ;
    }
###  Check if site is anionic or cationic
    if (Tot_p > 0.99) type = "C" ;   # site is cationic
    else type = "A" ;  # site is anionic
    close(st_file) ;

## Read all tauts (check all taut files)
    for (t = 1 ; t <= Taut[s]; t++) {

      filecheck(st_file = site[s,t] ".st") ;
      getline < st_file ; 

###  Read charges from the non ref state (Tautomer case)
### IF site is cationic read column 4 else read column 3
      for (i = 1; i <= l; i++) {

        if (type == "C") {

	  getline < st_file ;
	  s_resn_n[s,t,i] = $1 ;
	  s_atom_n[s,t,i] = $2 ;
	  s_chg_Ref[s,t,i] = $3 ;
	  s_chg_noRef[s,t,i] = $4 ;
	} else if(type == "A") {

	  getline < st_file ;
	  s_resn_n[s,t,i] = $1 ;
	  s_atom_n[s,t,i] = $2 ;
	  s_chg_noRef[s,t,i] = $3 ;
	  s_chg_Ref[s,t,i] = $4 ;
	}
      }
      close(st_file) ;
    }

###  Assigning number of atoms of each site
    s_atoms[s] = l ;
  }


###  Writing charges for all states 
  for (pt_at = 1; pt_at <= natoms; pt_at++) {

    found[pt_at] = 0 ;
    for (s = 1; s <= n; s++) {

      if ((res_i[pt_at] == site_i[s]) &&	\
	  (tempR[pt_at] == tempS[s,1])) {

	for (mdl = 1; mdl <= s_atoms[s]; mdl++) {

	  if (atom_n[pt_at] == s_atom_n[s,1,mdl]) {

	    found[pt_at] = 1 ;
	    for (t = 1; t <= Taut[s]; t++) {

	      if (t == 1) process_atom(out,0) ;
	      else  process_atom(out,1) ;
	      if (t == Taut[s]) printf("\n") > out ;
	    }
	  }
	}   # [model atoms - mdl]
      }
    }

    if (found[pt_at] == 0) process_atom(out,-1) ;
  }  # [protein atoms - pt_at]

}



#######################################################
#                      FUNCTIONS                      #
#######################################################


function error(msg) {

  print cmd ": Error: " msg | "cat 1>&2" ;
  close ("cat 1>&2") ;
  exit 1 ;
}


function warning(msg) {

  print cmd ": Warning: " msg | "cat 1>&2" ;
  close ("cat 1>&2") ;
}


###  Function to check existence of files
function filecheck(file) {

  if (system("test -f "file))
    error("File "file" does not exist.") ;
  if (system("test -r "file))
    error("File "file" exists but is not readable.") ;
}


###  Writes atom coordinates to pdb
function process_atom(file,flag) {

  if (flag == 0)  # Exists in sites and first is ref

    printf("%9.3f %8.3f",s_chg_Ref[s,t,mdl],
	   s_chg_noRef[s,t,mdl]) >file ;

  else if (flag == 1)  # Exists in sites and is non Ref

    printf("%9.3f",s_chg_noRef[s,t,mdl]) >file ;
  else if (flag == -1)  # Does not exist in sites (equal to pqr)

    printf("%9.3f\n",charge[pt_at]) >file ;
}



