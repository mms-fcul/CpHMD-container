#!/usr/bin/awk -f

##############################################################
#                                                            #
#         Vitor H. Teixeira, Lisboa, 2013.03.26              #
#                                                            #
#     Program to generate files to run delphi. It only       #
#  supports calculations with tautomers, one focus and has   #
#  periodic boundary conditions capabilities.                #
#     It is assumed (for computational simplicity) that PQR  #
#  has ITS SITE ATOMS in the charged reference state. This   #
#  was not fully tested with proteins.                       #
#     It does not matter if residue name is ARG or AR[0-4],  #
#  the program uses the first 2 characters to compare with   #
#  the residue in sites file and then it generates all the   #
#  necessary files with the needed names.                    #
#                                                            #
#     Small explanation of delphi parameters can be found    #
#  in the DELPHI.pbp file.                                   #
#                                                            #
#   Periodic boundary conditions: When doing calculations    #
# with PBC; non-linear iteractions, relaction factor for     #
# linear and non-linear are only used in the coarse grid.    #
# Focus and model are not calculated with PBC nor with       #
# non-linear PB, so their relative parameters are not        #
# considered. Number of linear PB iteractions is fixed in    #
# the focus/model runs.                                      #
#                                                            #
##############################################################

BEGIN{
  cmd = "gen_files.awk" ;
  usage = "Usage: "cmd" <sites> <pqr> <crg_file> <DELPHI.pbp> <vectors> <centers> <Dim>" ;
  if (ARGC != 8) error("Wrong number of arguments.\n" usage) ;

  filecheck(sites = ARGV[1]) ;
  filecheck(pqr_file = ARGV[2]) ;
  filecheck(crg = ARGV[3]) ;
  filecheck(pbp = ARGV[4]) ;
  vec = ARGV[5] ;
  filecheck(centers = ARGV[6]) ;
  dimension = ARGV[7] ;

  split(vec,aux,",") ;
  if (dimension == "0") {
    CentX = CentY = CentZ = 0.0 ;
  } else if (dimension == "2") {
    CentX = CentY = aux[1] / 2 ;
  } else error("Wrong definition of dimension.\n" usage) ;


  for (i == 1; i <= 12; i++) rd_pbp[i] = 0 ;

### Reading pbp file (DELPHI.pbp)

  while (getline < pbp) {

    if ($1 !~ "#") {

      if (substr($2,1,4) == "ioni") {
	split($2,rd,"=") ;
	ion = rd[2] ;
	rd_pbp[1] = 1 ;
      }
      if (substr($2,1,4) == "epsi") {
        split($2,rd,"=") ;
	epsin = rd[2] ;
	rd_pbp[2] = 1 ;
      }
      if (substr($2,1,6) == "perfil") {
	if (substr($2,7,2) == "Pb") {
	  split($2,rd,"=") ;
	  pfP = rd[2] ;
	} else if (substr($2,7,2) == "Pf") {
	  split($2,rd,"=") ;
	  pfPf = rd[2] ;
	} else if (substr($2,7,1) == "M") {
	  split($2,rd,"=") ;
	  pfM = rd[2] ;
	}
	rd_pbp[3] = 1 ;
      }
      if (substr($2,1,5) == "gsize") {
	if (substr($2,6,1) == "=") {
	  split($2,rd,"=") ;
	  gsP = gsPf = gsM = rd[2] ;
	} else if (substr($2,6,1) == "P") {
          split($2,rd,"=") ;
	  gsP = rd[2] ;
	} else if (substr($2,6,1) == "M") {
	  split($2,rd,"=") ;
	  gsPf = gsM = rd[2] ;
	}
	rd_pbp[4] = 1 ;
      }
      if (substr($2,1,4) == "scal") {
	if (substr($2,6,1) == "P") {
          split($2,rd,"=") ;
	  scP = rd[2] ;
	} else if (substr($2,6,1) == "M") {
	  split($2,rd,"=") ;
	  scPf = scM = rd[2] ;
	}
	rd_pbp[5] = 1 ;
      }
      if (substr($2,1,4) == "bndc") {
        split($2,rd,"=") ;
	bd = rd[2] ;
	rd_pbp[6] = 1 ;
      }
      if (substr($2,1,4) == "maxc") {
        split($2,rd,"=") ;
	maxc = rd[2] ;
	rd_pbp[7] = 1 ;
      }
      if (substr($2,1,4) == "rmsc") {
        split($2,rd,"=") ;
	rmsc = rd[2] ;
	rd_pbp[8] = 1 ;
      }
      if (substr($2,1,4) == "lini") {
        split($2,rd,"=") ;
	linit = rd[2] ;
	rd_pbp[9] = 1 ;
      }
      if (substr($2,1,4) == "noni") {
	split($2,rd,"=") ;
	nonit = rd[2] ;
	rd_pbp[10] = 1 ;
      }
      if (substr($2,1,4) == "relf") {
	split($2,rd,"=") ;
	relaxL = rd[2] ;
	rd_pbp[11] = 1 ;
      }
      if (substr($2,1,4) == "relp") {
	split($2,rd,"=") ;
	relaxN = rd[2] ;
	rd_pbp[12] = 1 ;
      }
      if (substr($2,1,4) == "fcrg") fcrg = 1 ;
      if (substr($2,1,3) == "pbx") pbx = 1 ;
      if (substr($2,1,3) == "pby") pby = 1 ;
      if (substr($2,1,3) == "pbz") pbz = 1 ;
    }
  }
  close(pbp) ;

###  Checking grid parameters readed
  if (rd_pbp[3] == 1 && rd_pbp[4] == 1 && rd_pbp[5] == 1)
    error("\nYou cannot use perfil, gsize and scale at same time\n\n") ;


###  Reading sites file (Tautomer case)
  while (getline < sites) {

    n++ ;
    Taut[n] = NF - 1 ;
    site_i[n] = $1 ;
    for (i = 1; i <= Taut[n]; i++) {

      site[n,i] = $(i+1) ;
      tempS[n,i] = substr(site[n,i],1,2) ;
    }
  }
  close(sites) ;



###  Reading zeroed crg_file
  while (getline < crg) {

    k++ ;
    if (NF == 3) {
      crg_a[k] = $1 ;
      crg_r[k] = $2 ;
      crg_c[k] = $3 ;
    } else if (NF == 2) {
      crg_a[k] = $1 ;
      crg_c[k] = $2 ;
    }
  }
  close (crg) ;

###  Reading PQR file
  while (getline < pqr_file) {

    m++ ;
    atom_n[m] = $3 ;
    res_n[m] = substr($4,1,3) ;
    tempR[m] = substr($4,1,2) ;
    res_i[m] = $5 ;
    x[m] = $6 ;
    y[m] = $7 ;
    z[m] = $8 ;
    ChgRef[m] = $9 ;     # Assuming current site in PQR is in Reference state
    Chg_noRef[m] = $9 ;
    Rad[m] = $10 ;
  }
  close(pqr_file) ;

###  Reading each site center (calculated previously)
  if (dimension == 2) {
    while (getline < centers) {

      cent++ ;
      c_site[cent] = $1 ;
      c_z[cent] = $4 * 10 ;
    }
    close (centers) ;
  } else if (dimension == 0) {
    while (getline < centers) {

      cent++ ;
      c_site[cent] = $1 ;
      c_x[cent] = $2 * 10 ;
      c_y[cent] = $3 * 10 ;
      c_z[cent] = $4 * 10 ;
    }
    close (centers) ;
  }

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

###  Read charges from the non ref state 
### IF site is cationic read column 4 else read column 3
      for (i = 1; i <= l; i++) {

        if (type == "C") {

	  getline < st_file ;
	  s_resn_n[s,t,i] = $1 ;
	  s_atom_n[s,t,i] = $2 ;
	  s_chg_noRef[s,t,i] = $4 ;
	} else if(type == "A") {

	  getline < st_file ;
	  s_resn_n[s,t,i] = $1 ;
	  s_atom_n[s,t,i] = $2 ;
	  s_chg_noRef[s,t,i] = $3 ;
	}
      }
      close(st_file) ;
    }

###  Assigning number of atoms of each site
    s_atoms[s] = l ;
  }



###  Processing PQR file to create PDB and PRM files for all sites
### The number of PRM and PDB files will be the same and correspond to
### n+1, where n is the number of taut of each site. A File with 
### charges was readed and appended with the site protonation state.
  for (s = 1; s <= n; s++) {

###  Initializing variables for model coordinate center
    for (t = 1; t <= Taut[s]; t++) {

      a = length(site[s,t]) ;
      out = sprintf("%d-%s%s%d",site_i[s],substr(site[s,t],1,1),
		    tolower(substr(site[s,t],2,a-5)),
		    substr(site[s,t],a,1)) ;


###   Write parameter file. This will generate 2*(n+1) prm files, 
###  corresponding to the big grid and focus grid parameters. One prm
###  file for each taut +1 will also be generated for the models.
###  Arguments for write_parm - Filename (other parameters are global)

      if (t == Taut[s]) {

	write_parmM(out"N.prm") ;
	write_parmM(out"C.prm") ;
	write_parmPf("P_"out"Nf.prm") ;
	write_parmPb("P_"out"Nb.prm") ;
	write_parmPf("P_"out"Cf.prm") ;
	write_parmPb("P_"out"Cb.prm") ;
	write_crg_header("P_"out"N.crg") ;
	write_crg_header("P_"out"C.crg") ;
      } else {

	write_parmM(out"N.prm") ;
	write_parmPb("P_"out"Nb.prm") ;
	write_parmPf("P_"out"Nf.prm") ;
	write_crg_header("P_"out"N.crg") ;
      }

      for (p = 1; p <= m; p++) {

	pp[p] = ppc[p] = 0 ;

###  For model compounds

	if ((res_i[p] == site_i[s]) && (tempR[p] == tempS[s,t])) {

	  for (mdl = 1; mdl <= s_atoms[s]; mdl++) {

	    if (atom_n[p] == s_atom_n[s,t,mdl]) {

	      Chg_noRef[s,t,mdl] = s_chg_noRef[s,t,mdl] ;
	      if (t == Taut[s]) 
		process_atom(out"C.pdb",t,0) ;
	      process_atom(out"N.pdb",t-1,0) ;
	      break ;
	    }
	  }  # [model atoms - mdl]
	}

###  Creating Protein PDB file with only site s not in reference
### All other sites are writen in their usual name 

##  If atom is from current site write its particular taut state

	if ((tempS[s,t] == tempR[p]) && (res_i[p] == site_i[s])) {

	  for (mdl = 1; mdl <= s_atoms[s]; mdl++) {

	    if (atom_n[p] == s_atom_n[s,t,mdl]) {

	      pp[p] = 1 ;
	      process_atom("P_"out"N.pdb",t-1,0) ;
	      printf("%-4s  %-2s%1d%4d %7.3f\n",atom_n[p],
		     substr(site[s,t],1,2),t-1,
		     res_i[p],Chg_noRef[s,t,mdl]) >"P_"out"N.crg" ;
	      if (t == Taut[s]) {

		ppc[p] = 1 ;
		process_atom("P_"out"C.pdb",t,0) ;
		printf("%-4s  %-2s%1d%4d %7.3f\n",atom_n[p],
		       substr(site[s,t],1,2),t,
		       res_i[p],ChgRef[p]) >"P_"out"C.crg" ;
	      }
	    }
	  }  # [model atoms - mdl]
	}

###  If atom does not belong to site s nor is from a titrable terminal
### write usual sites name (for pdb writting)
	if (pp[p] == 0)
	  process_atom("P_"out"N.pdb",t,-1) ;
	if (ppc[p] == 0 && t == Taut[s])
	  process_atom("P_"out"C.pdb",t,-1) ;

      }  # [protein atoms - p]


###  If atom does not belong to site s write 0.00 charge
      for (c = 1; c <= k; c++)
	printf("%-4s  %-3s %11.3f\n",crg_a[c],
	       crg_r[c],crg_c[c]) >"P_"out"N.crg" ;

      if (t == Taut[s]) {
	for (c = 1; c <= k; c++) {

	  printf("%-4s  %-3s %11.3f\n",crg_a[c],
		 crg_r[c],crg_c[c]) >"P_"out"C.crg" ;
	}
      }

    }  # [Tautomers - t]

  }  # [sites - s]
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


# Write parameter for the big file (Before focusing) 
# Filename as argument 
# ionic_strength, dielectric_constant, points, spacing, boundary (Global)
function write_parmPb(file) {

  if (rd_pbp[3] == 1 && rd_pbp[5] == 1)
    printf("perfil=%f\nscale=%.10f\nindi=%f\n",pfP,scP,epsin) > file ;
  if (rd_pbp[4] == 1 && rd_pbp[5] == 1)
    printf("gsize=%d\nscale=%.10f\nindi=%f\n",gsP,scP,epsin) > file ;
  if (rd_pbp[3] == 1 && rd_pbp[4] == 1)
    printf("perfil=%f\ngsize=%d\nindi=%f\n",pfP,gsP,epsin) > file ;
  printf("exdi=80.0\nprbrad=1.4\nbndcon=%d\nmaxc=%f\n",bd,maxc) > file ;
  if (rd_pbp[8] == 1) printf("rmsc=%f\n",rmsc) >file ;
  printf("linit=%d\nenergy(s,c)\nsalt=%f\n",linit,ion) > file ;
  if (rd_pbp[10] == 1) printf("nonit=%d\n",nonit) > file ;
  if (rd_pbp[11] == 1) printf("relfac=%f\n",relaxL) > file ;
  if (rd_pbp[12] == 1) printf("relpar=%f\n",relaxN) > file ;
  if (fcrg == 1) printf("fcrg=true\n") > file ;
  if (pbx == 1) printf("pbx=true\n") > file ;
  if (pby == 1) printf("pby=true\n") > file ;
  if (pbz == 1) printf("pbz=true\n") > file ;
  printf("in(siz,file=\"DataBaseT.siz\")\n") > file ;

  comp = length(file) ;
  fl1 = substr(file,1,comp-5) ;

  printf("in(crg,file=\"%s.crg\")\n",fl1) > file ;
  printf("in(pdb,file=\"%s.pdb\")\n",fl1) > file ;
  printf("out(phi,file=\"%s.phi\")\n",fl1) > file ;
  printf("in(frc,file=\"%s.pdb\")\n",fl1) > file ;
  printf("site(a,q,p)\nout(frc,file=\"%sb.frc\")\n",fl1) > file ;
  if (dimension == "2") {
    for (a = 1; a <= cent; a++)
      if (c_site[a] == site_i[s]) 
        printf("ACent(%f,%f,%f)\n",CentX,CentY,c_z[a]) > file ;
  } else if (dimension == "0") {
    for (a = 1; a <= cent; a++)
      if (c_site[a] == site_i[s]) 
	printf("ACent(%f,%f,%f)\n",c_x[a],c_y[a],c_z[a]) > file ;
  }
}


# Write parameter for the focus file (Reads phi from big run)
# Filename as argument 
# ionic_strength, dielectric_constant, points, spacing (Global)
# Boundary is 3 (focus) and number of linear interactions is 500
function write_parmPf(file) {

  if (rd_pbp[3] == 1 && rd_pbp[5] == 1)
    printf("perfil=%f\nscale=%.10f\nindi=%f\n",pfPf,scPf,epsin) > file ;
  if (rd_pbp[4] == 1 && rd_pbp[5] == 1)
    printf("gsize=%d\nscale=%.10f\nindi=%f\n",gsPf,scPf,epsin) > file ;
  if (rd_pbp[3] == 1 && rd_pbp[4] == 1)
    printf("perfil=%f\ngsize=%d\nindi=%f\n",pfPf,gsPf,epsin) > file ;
  printf("exdi=80.0\nprbrad=1.4\nbndcon=%d\nmaxc=%f\n",3,maxc) > file ;
  if (rd_pbp[8] == 1) printf("rmsc=%f\n",rmsc) >file ;
  printf("linit=%d\nenergy(s,c)\nsalt=%f\n",500,ion) > file ;
#  if (rd_pbp[10] == 1) printf("nonit=%d\n",nonit) > file ;
#  if (rd_pbp[11] == 1) printf("relfac=%f\n",relaxL) > file ;
#  if (rd_pbp[12] == 1) printf("relpar=%f\n",relaxN) > file ;
#  if (fcrg == 1) printf("fcrg=true\n") > file ;
#  if (pbx == 1) printf("pbx=true\n") > file ;
#  if (pby == 1) printf("pby=true\n") > file ;
#  if (pbz == 1) printf("pbz=true\n") > file ;
  printf("in(siz,file=\"DataBaseT.siz\")\n") > file ;

  comp = length(file) ;
  fl2 = substr(file,1,comp-5) ;

  printf("in(crg,file=\"%s.crg\")\n",fl2) > file ;
  printf("in(pdb,file=\"%s.pdb\")\n",fl2) > file ;
  printf("in(phi,file=\"%s.phi\")\n",fl2) > file ;
  printf("in(frc,file=\"%s.pdb\")\n",fl2) > file ;
  printf("site(a,q,p)\nout(frc,file=\"%sf.frc\")\n",fl2) > file ;
  if (dimension == "2") {
    for (a = 1; a <= cent; a++)
      if (c_site[a] == site_i[s]) 
        printf("ACent(%f,%f,%f)\n",CentX,CentY,c_z[a]) > file ;
  } else if (dimension == "0")
    for (a = 1; a <= cent; a++)
      if (c_site[a] == site_i[s]) 
	printf("ACent(%f,%f,%f)\n",c_x[a],c_y[a],c_z[a]) > file ;
}


###  Writes parameters for the model compound 
##  Center is written in another function
# number of linear interactions is 500
function write_parmM(file,z) {

  if (rd_pbp[3] == 1 && rd_pbp[5] == 1)
    printf("perfil=%f\nscale=%.10f\nindi=%f\n",pfM,scM,epsin) > file ;
  if (rd_pbp[4] == 1 && rd_pbp[5] == 1)
    printf("gsize=%d\nscale=%.10f\nindi=%f\n",gsM,scM,epsin) > file ;
  if (rd_pbp[3] == 1 && rd_pbp[4] == 1)
    printf("perfil=%f\ngsize=%d\nindi=%f\n",pfM,gsM,epsin) > file ;
  printf("exdi=80.0\nprbrad=1.4\nbndcon=%d\nmaxc=%f\n",bd,maxc) > file ;
  if (rd_pbp[8] == 1) printf("rmsc=%f\n",rmsc) >file ;
  printf("linit=%d\nenergy(s,c)\nsalt=%f\n",500,ion) > file ;
#  if (rd_pbp[10] == 1) printf("nonit=%d\n",nonit) > file ;
#  if (fcrg == 1) printf("fcrg=true\n") > file ;
#  if (pbx == 1) printf("pbx=true\n") > file ;
#  if (pby == 1) printf("pby=true\n") > file ;
#  if (pbz == 1) printf("pbz=true\n") > file ;
  printf("in(crg,file=\"DataBaseT.crg\")\n") > file ;
  printf("in(siz,file=\"DataBaseT.siz\")\n") > file ;

  cmp = length(file) ;
  fl5 = substr(file,1,cmp-4) ;

  printf("in(pdb,file=\"%s.pdb\")\n",fl5) >file ;
  if (dimension == "2") {
    for (a = 1; a <= cent; a++)
      if (c_site[a] == site_i[s]) 
        printf("ACenter(%f,%f,%f)\n",CentX,CentY,c_z[a]) > file ;
  } else if (dimension == "0")
    for (a = 1; a <= cent; a++)
      if (c_site[a] == site_i[s]) 
	printf("ACent(%f,%f,%f)\n",c_x[a],c_y[a],c_z[a]) > file ;
}


###  Writes the header of the crg data file
function write_crg_header(file) {

  printf("!crg file created by gen_files.awk\n") >file ;
  printf("atom__resnumbc_charge_\n") >file ;
}


###  Writes atom coordinates to pdb
function process_atom(file,tautnumb,resflag) {

  if (resflag < 0)   # equal PQR

    printf("ATOM %6d %-4s %-4s %4d %11.3f%8.3f%8.3f\n",
	     p,atom_n[p],res_n[p],res_i[p],
	     x[p],y[p],z[p]) >file ;

  else if (resflag == 0)

    printf("ATOM %6d %-4s %-2s%d  %4d %11.3f%8.3f%8.3f\n",
	   p,atom_n[p],substr(site[s,t],1,2),tautnumb,res_i[p],
	   x[p],y[p],z[p]) >file ;

  else if (resflag > 0)
    printf("ATOM %6d %-4s %-2s%d  %4d %11.3f%8.3f%8.3f\n",
	   p,atom_n[p],substr(site[ss,tempt],1,2),tautnumb,res_i[p],
	   x[p],y[p],z[p]) >file ;
}


