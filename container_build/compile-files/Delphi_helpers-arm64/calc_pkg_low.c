
/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "my_mem.h"


/* Constants */
#define MAXLINESIZE         256
#define MAXSITES           1024
#define MAXPSEUDOSITES        8
#define MAXSITENAME          32
#define MAXSITEATOMS        128
#define MAXATOMNAME           5

/*  #define MAXATOMBUF       131072       For the huge line (buffer) of atoms per site  */
/*  #define MAXBACKATOM       32768       For the indexes of those atoms  */


/* Memory Functions */
/* In file my_mem.c */



/* Template Functions */
void usage                   (void) ;
void error                   (char *text) ;
int lines_in_file            (FILE *ifile) ;
int big_lines_in_file        (FILE *ifile) ;
void read_sites              (FILE *ifile, int n_sites, int *n_taut,
			      int *site_i, char ***site_name, 
			      char ***solv_site) ;
void read_pqr                (FILE *ifile, int n_atoms, 
			      float **charge) ;
void read_solv               (FILE *ifile, int n_sites, int *site_i, 
			      int n_solv, float **solv_s, 
			      float **solv_p, char ***write, 
			      char ***solv_site) ;
void read_sts                (int n_sites, char ***site_name, int *n_taut,
			      char **site_type, float **site_prot, 
			      char ***st_atom_n, int *st_atoms) ;
void read_pkmod              (int n_sites, int *n_taut, 
			      char ***site_name, float **pkmod) ;
void read_pair_inter         (FILE *ifile, int n_pair_inter, int *write_si,
			      int *write_sj, int *site_si, int *site_sj,
			      int *site_ti, int *site_tj) ;
void read_back_inter         (FILE *ifile, int n_back_inter,
			      int *site_bk, int *line_bk,
			      int **atom_bk) ;
void read_charges            (FILE *ifile, int n_atoms, 
			      float **charge_states) ;
void read_frc                (int n_sites, int *n_taut, 
			      char ***site_name, char ***solv_site,
			      int *site_i, int frc_lines, 
			      char ****frc_An, int ***frc_Ri, 
			      float ***frc_Pt) ;


/*  Global variables  */
char focus[MAXSITENAME] ;



/* Main Program */

int main(int argc, char **argv) {

  FILE *isites, *ipqr, *isolv, *ipairlist ;
  FILE *ibacklist, *icharges, *oout, *itemp ;

  char sites[MAXSITENAME], pqr[MAXSITENAME], solv[MAXSITENAME] ;
  char pairlist[MAXLINESIZE], backlist[MAXLINESIZE], charges[MAXLINESIZE] ;
  char out[MAXSITENAME], ftemp[MAXSITENAME] ;
  char ***site_name, ***solv_site, **site_type, ***st_atom_n ;
  char ***write_pk, ****frc_An ;

  int s, t, l, n_sites, n_solv, frc_lines, n_atoms, pt_at, st_at ;
  int n_pair_inter, n_back_inter, temp_sj, sj, temp_si, si, pair ;
  int *n_taut, *site_i, *st_atoms, ***frc_Ri ;
  int *write_si, *write_sj, *site_si, *site_sj, *site_ti ;
  int *site_tj, *site_bk, *line_bk ;
  int **atom_bk ;

  float **chg, **pkmod, **solv_s, **solv_p, **back, **site_prot ;
  float ***frc_Pt, Temperature, bg ;
  float site_solv, prot_solv, pkint, sum_temp, interaction ;
  float conv_pk , conv_g, kBoltz_au = 5.98435e-6 ; /* e^2/(Angstrom*K) */


/* Read Arguments */
  if (argc != 10) {

    error("Wrong number of arguments\n\n") ;
    usage() ; exit(1) ;
  }

  strcpy(sites,argv[1]) ;
  strcpy(pqr,argv[2]) ;
  strcpy(solv,argv[3]) ;
  strcpy(pairlist,argv[4]) ;
  strcpy(backlist,argv[5]) ;
  strcpy(charges,argv[6]) ;
  strcpy(focus,argv[7]) ;
  strcpy(ftemp,argv[8]) ;
  Temperature = atof(ftemp) ;
  strcpy(out,argv[9]) ;

/*  Conversions from delphi units (KT) to the MEAD out
 * normal units: pk (for pKint) and e^2/A (for interactions)
 * */

  conv_pk = 2.30259 ;  /* from KT to pK units */
  conv_g = kBoltz_au * Temperature ;




/*  Testing focus option  */
  if ((strcmp(focus,"yes") != 0) && 
      (strcmp(focus,"no") != 0)) {

    error("Please choose <yes> or <no> for focus\n\n") ;
    usage() ; exit(1) ;
  }

/*  Test opening files */
  if ((isites = fopen(sites, "r")) == NULL) {
    error("Can not open sites file\n") ;
    exit(1) ;
  }
  if ((ipqr = fopen(pqr, "r")) == NULL) {
    error("Can not open pqr file\n") ;
    exit(1) ;
  }
  if ((isolv = fopen(solv, "r")) == NULL) {
    error("Can not open solvation data file\n") ;
    exit(1) ;
  }
  if ((ipairlist = fopen(pairlist, "r")) == NULL) {
    error("Can not open pairwise interaction list file\n") ;
    exit(1) ;
  }
  if ((ibacklist = fopen(backlist, "r")) == NULL) {
    error("Can not open background interaction list file\n") ;
    exit(1) ;
  }
  if ((icharges = fopen(charges, "r")) == NULL) {
    error("Can not open charges file\n") ;
    exit(1) ;
  }


/***************************************************
 *
 *                  READING DATA
 *
 ***************************************************/


/**********************************************
 *          Reading sites file
 **********************************************/

/*  Count the number of lines in sites file */
  n_sites = lines_in_file(isites) ;
  fclose(isites) ;

/*  Allocate memory for all sites. Since a variable  
 * number of tautomers may exist, allocate for the 
 * maximum number in each line (each site) */ 
  n_taut = calloc(n_sites,sizeof(int)) ;
  site_i = calloc(n_sites,sizeof(int)) ;
  site_name = alloc_3D_char(n_sites,MAXPSEUDOSITES,MAXSITENAME) ;
  solv_site = alloc_3D_char(n_sites,MAXPSEUDOSITES,3) ;

/*  Count the number of words (tauts) in each line
 * and read them. Save also information of the site
 * name in a temp variable to compare with when
 * reading solvation energies.  */
  isites = fopen(sites, "r") ;
  read_sites(isites,n_sites,n_taut,site_i,site_name,solv_site) ;
  fclose(isites) ;



/**********************************************
 *          Reading pqr file
 **********************************************/

/*  Count the number of lines (atoms) in pqr file */
  n_atoms = lines_in_file(ipqr) ;
  fclose(ipqr) ;

/*  Allocate memory for the charge field in pqr file */
  chg = alloc_2D_float(n_atoms,MAXPSEUDOSITES) ;
  /*  chg_ref = calloc(n_atoms,sizeof(float)) ;*/

/*  Read pqr file. Only interested in the charge column */
  ipqr = fopen(pqr, "r") ;
  read_pqr(ipqr,n_atoms,chg) ;
  fclose(ipqr) ;



/**********************************************
 *          Reading solvation file
 **********************************************/

/*  Count lines in solvation energy file */
  n_solv = lines_in_file(isolv) ;
  fclose(isolv) ;

/*  The amount of memory to be allocated corresponds to
 * the double of the total number of states (site + Protein).
 * Also save the varible with the name to write in pk output. */
  solv_s = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;
  solv_p = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;
  write_pk = alloc_3D_char(n_sites,MAXPSEUDOSITES,16) ;

/*  Read file  */
  isolv = fopen(solv, "r") ;
  read_solv(isolv,n_sites,site_i,n_solv,solv_s,solv_p,write_pk,solv_site) ;
  fclose(isolv) ;



/**********************************************
 *    Reading st files according to sites info
 **********************************************/

/*  Allocate memory for information on st files, namelly 
 * site type, atom name, number of atoms. All info is taken 
 * from the first taut. */
  site_type = alloc_2D_char(n_sites,2) ;
  st_atom_n = alloc_3D_char(n_sites,MAXSITEATOMS,MAXATOMNAME) ;
  st_atoms = calloc(n_sites,sizeof(int)) ;
  site_prot = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;

/*  Read st files and their pKmod's  */
  read_sts(n_sites,site_name,n_taut,site_type,site_prot,st_atom_n,st_atoms) ;



/**********************************************
 *    Reading st files according to sites info
 *          EXTRACTING PKMOD ONLY
 **********************************************/

/*  Alloc memory for pkmod  */
  pkmod = alloc_2D_float(n_sites,*n_taut) ;

/*  Read pkmod from model compounds */
  read_pkmod(n_sites,n_taut,site_name,pkmod) ;



/**********************************************
 *     Reading pairwise interaction list
 **********************************************/

/*  Count the number of lines in interaction file */
  n_pair_inter = lines_in_file(ipairlist) ;
  fclose(ipairlist) ;

/*  Allocate memory for interaction list file */
  write_si = calloc(n_pair_inter,sizeof(int)) ;
  write_sj = calloc(n_pair_inter,sizeof(int)) ;
  site_si = calloc(n_pair_inter,sizeof(int)) ;
  site_ti = calloc(n_pair_inter,sizeof(int)) ;
  site_sj = calloc(n_pair_inter,sizeof(int)) ;
  site_tj = calloc(n_pair_inter,sizeof(int)) ;

/*  Read interaction pair list file. Reading everything 
 * except 1st, 4th and 7th columns */
  ipairlist = fopen(pairlist, "r") ;
  read_pair_inter(ipairlist,n_pair_inter,write_si,write_sj,site_si,site_sj,site_ti,site_tj) ;
  fclose(ipairlist) ;

/**********************************************
 *  Possible transformation of these data into
 * a data structure.
 * This would simplify the calling function */



/**********************************************
 *     Reading background interaction list
 **********************************************/


/*  Count the number of lines in back interaction file */
  n_back_inter = lines_in_file(ibacklist) ;
  fclose(ibacklist) ;

/*  Allocate memory for interaction list file */
  site_bk = calloc(n_back_inter,sizeof(int)) ;
  line_bk = calloc(n_back_inter,sizeof(int)) ;
  atom_bk = alloc_2D_int(n_back_inter,n_atoms) ;

  ibacklist = fopen(backlist, "r") ;
  read_back_inter(ibacklist,n_back_inter,site_bk,line_bk,atom_bk) ;




  /********  DEBUG *************


  {
    int i, j ;

    for (i = 0; i < n_back_inter; i++) {

      printf("%3d ",*(site_bk + i)) ;

      for (j = 0; j < *(line_bk + i); j++) {

	printf(" %d",atom_bk[i][j]) ;
      }
      printf("\n") ;
    }
  }

  
  ********  DEBUG *************/







/**********************************************
 *          Reading charges file
 **********************************************/

/*  Count the number of lines (atoms) in charges file.
 * Already counted in the pqr. */

/*  Allocate memory for all charges in file */
/*  chg = alloc_2D_float(n_atoms,MAXPSEUDOSITES) ;*/
/* Already allocated in the pqr reading */

/*  Read charges file ... everything  */
  icharges = fopen(charges, "r") ;
  read_charges(icharges,n_atoms,chg) ;
  fclose(icharges) ;



/**********************************************
 *    Reading frc files according to sites info
 **********************************************/

/*  Assuming the above read pqr file has all the needed atoms,
 * that number will be used to read frc files ignoring all the
 * atoms added by the slice program. */

  frc_lines = n_atoms ;

/*  Allocate memory for the Atom name, Residue number, 
 * Charge and Potential. */
  frc_An = alloc_4D_char(n_sites,MAXPSEUDOSITES,frc_lines,5) ;
  frc_Ri = alloc_3D_int(n_sites,MAXPSEUDOSITES,frc_lines) ;
  frc_Pt = alloc_3D_float(n_sites,MAXPSEUDOSITES,frc_lines) ;

  read_frc(n_sites,n_taut,site_name,solv_site,site_i,
	   frc_lines,frc_An,frc_Ri,frc_Pt) ;

/**********************************************
 *  Possible transformation of these data into
 * a data structure.
 * This would simplify the calling function */



/***************************************************
 *
 *        Calculating background contribution
 *
 ***************************************************/
/* Allocating memmory for background interactions and 
 * for the switch that will retain if atom site was found */
  back = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;



  for (s = 0; s < n_sites; s++) {
    for (t = 0; t <= *(n_taut + s); t++) {

      back[s][t] = 0.0 ;
    }
    for (l = 0; l < *(line_bk + s); l++) {

      for (t = 0; t <= *(n_taut + s); t++) {

	back[s][t] += chg[atom_bk[s][l]][0] * frc_Pt[s][t][atom_bk[s][l]] ;
      }
    }
  }


/*  Setting name for contributions file  */
  sprintf(ftemp,"Contributions") ;

  itemp = fopen(ftemp, "w") ;



/***************************************************
 *
 *      Combining solvation and background energies
 *    and calculate pkint values 
 *
 ***************************************************/

  sprintf(ftemp,"%s.pkint",out) ;
  oout = fopen(ftemp, "w") ;

  for (s = 0; s < n_sites; s++) {
    for (t = 0; t < *(n_taut + s); t++) {

      site_solv = solv_s[s][0] - solv_s[s][t+1] ;
      prot_solv = solv_p[s][0] - solv_p[s][t+1] ;

/*  Converting the previously calculated background energies
 * to pk units. This may actually be temperature dependent */
      bg = (back[s][0] - back[s][t+1]) / conv_pk ;

/*  Printing the different contributions to the pKint.
 */
      fprintf(itemp,"%6s-%-3d   %12.6f   %12.6f\n",write_pk[s][t],
	      site_i[s],(prot_solv - site_solv),bg) ;

      if (strcmp(site_type[s],"A") == 0)
	pkint = pkmod[s][t] + (prot_solv - site_solv + bg) ;

      else 
	pkint = pkmod[s][t] - (prot_solv - site_solv + bg) ;

      fprintf(oout,"%12.6e %s %c%c%c%d-%d %.0f\n",pkint,
	      site_type[s],write_pk[s][t][0],
	      toupper(write_pk[s][t][1]),
	      toupper(write_pk[s][t][2]),
	      t+1,site_i[s],site_prot[s][t]) ;
    }
  }


/***************************************************
 *
 *        Calculating pairwise interactions based
 *     on the list previously read ... and write them
 *
 ***************************************************/

  sprintf(ftemp,"%s.g",out) ;
  oout = fopen(ftemp, "w") ;

/*  For each pair in the list, go through all protein atoms 
 * and then through all site atoms until finding site */
  for (pair = 0; pair < n_pair_inter; pair++) {    /* for each pair */

    sum_temp = si = sj = 0.0 ;
    for (temp_sj = 0; temp_sj < n_sites; temp_sj++) {
      if (*(site_i + temp_sj)  == *(site_sj + pair)) {

	sj = temp_sj ;
	break ;
      }
    }
    for (temp_si = 0; temp_si < n_sites; temp_si++) {
      if (*(site_i + temp_si)  == *(site_si + pair)) {

	si = temp_si ;
	break ;
      }
    }
    for (pt_at = 0; pt_at < frc_lines ; pt_at++) {  /* protein atom */
      for (st_at = 0; st_at < st_atoms[sj] ; st_at++) {  /* site atoms */

	if ((site_i[sj] == frc_Ri[si][*(site_ti + pair)][pt_at]) && 
	    (strcmp(st_atom_n[sj][st_at],frc_An[si][*(site_ti + pair)][pt_at]) == 0)) {

	  sum_temp += ((frc_Pt[si][0][pt_at] -	\
			frc_Pt[si][*(site_ti + pair)][pt_at]) * \
		       (chg[pt_at][0] - chg[pt_at][*(site_tj + pair)])) ;
	  break ;
	}
      }  /*  site atoms (st_at)  */
    }  /*  protein atoms (pt_at)  */

/*  The interaction should be actually sum_temp / 2, but since
 * we are using only one site, we had to multiplying per 2 which
 * actually simplifies with the division */
    interaction = fabs(sum_temp * conv_g) ;
    fprintf(oout,"%d %d %15.6e\n",*(write_si + pair),*(write_sj + pair),
	    interaction) ;
  }


  return 0 ;
}



/*********************************************/
/*********************************************/
/*********************************************/
/***************FUNCTIONS*********************/
/*********************************************/
/*********************************************/
/*********************************************/


/* Show usage */
void usage(void) {

  fprintf(stderr, "Usage: calc_pka <sites> <pqr> <solv> <plist> <blist>"
	  " <charges> <flag> <temperature> <out>\n") ;
  fprintf(stderr, " * sites contain the groups whose pKint is to be "
	  "calculated.\n * solv is a file containing the solvation energy "
	  "of each \n    site in the \"system\" and as a model compound.\n"
	  " * pqr is a standard pqr file.\n * solv is a table of the "
	  "solvation energy of the site in the\n    system and as a model "
	  "in both protonated and deprotonated\n    states.\n") ;
  fprintf(stderr," * plist is a discriminated list of all pairwise "
	  "interactions\n    to be calculated.\n * blist is the "
	  "discriminated list of atoms to include in the\n"
	  "    background calculation, using the same cutoff as in "
	  "the\n    pairwise interactions.\n * charges is the "
	  "database of atomic charges set to zero.\n") ;
  fprintf(stderr," * flag is for the presence or not of "
	  "focus. use yes if you\n    are using focus, no otherwise.\n"
	  " * temperature is used for the units conversion from KT to \n"
	  "    atomic units\n * out is the name without extension that "
	  "pkint and g file\n    will acquire.\n\n") ;
  return ;
}  /*  usage  */


/*********************************************/


void error(char *text) {

  fprintf(stderr, "\ncalc_pkg: ERROR: %s",text) ;
  return ;
}  /*  error  */



/*********************************************/



/*  Reads the number of lines in a file.
 *  http://stackoverflow.com/questions/4278845/count-the-lines-of-a-file-in-c
 *    This is a way to read a file and get the number of lines 
 * without any interest for the contents */
int lines_in_file(FILE *ifile) {

  int lines ;

  lines = 0 ;
  while (EOF != (fscanf(ifile,"%*[^\n]"), fscanf(ifile,"%*c"))) 
    lines++ ;

  return lines ;
}  /*  lines_in_file  */



/*********************************************/



/*  Count the number of words (tauts) in each line
 * and read them. */
void read_sites(FILE *ifile, int n_sites, int *taut, int *site_i, 
		char ***site_name, char ***solv_site) {

  char buf[MAXLINESIZE] ;
  int wordcount, i, j, k;

  for (j = 0; j < n_sites; j++) {

    if (fgets(buf, sizeof(buf), ifile) != 0) {

      wordcount = i = 0 ;
      while(buf[i]) {

	if (isgraph(buf[i]))
	  if (!isgraph(buf[i-1]))
	    wordcount++ ;
	i++ ;
      }
      *(taut + j) = wordcount - 1 ;
    }
  }

  rewind(ifile) ;
  for (j = 0; j < n_sites; j++) {

    if (fscanf(ifile,"%d",(site_i + j)) == 1) {

      for (k = 0; k < *(taut + j); k++) {

	if (fscanf(ifile,"%s",site_name[j][k]) == 1) {

	  strncpy(solv_site[j][k],site_name[j][k],1) ;
	  solv_site[j][k][1] = tolower(site_name[j][k][1]) ;
	  solv_site[j][k][2] = '\0' ;

	} 
      }
    }
  }
  return ;
}  /*  read_sites  */



/*********************************************/



/*  Read charges in the pqr file  */
void read_pqr(FILE *ifile, int n_atoms, float **chg) {

  int i ;

  for (i = 0; i < n_atoms; i++) {

    if (fscanf(ifile,"%*s%*s%*s%*s%*s%*s%*s%*s%f%*s",&chg[i][0]) == 1)
      {} ;  /*  NULL STATEMENT  */
  }
  return ;
}  /*  read_pqr  */



/*********************************************/



/* Reads the solvation energy data that was greped from 
 * the delphi output run. ch1 will have the residue name 
 * and number as well as tautomer to be split. pch, temp, 
 * aux_s, aux_i, aux_t are the temporary variables containing 
 * the return of strtok, for the tautomer number temp, 
 * state (N|C) residue number and tautomer number.
 */
void read_solv(FILE *f, int n_sites, int *site_i, int n_solv,
	       float **solv_s, float **solv_p, char ***write,
	       char ***solv_site) {

  char ch1[MAXSITENAME] ;
  char *pch, *temp, *aux_s ;
  int s, dim, i, aux_i, aux_t ;
  float aux;

/*  Reading solvation file one line at a time  */
  for (i = 0; i < n_solv; i++) {

    if (fscanf(f,"%s%f",ch1,&aux) == 2) {

      if (strncmp(ch1,"P_",2) == 0) {

	pch = strtok (ch1,"_-") ;
	pch = strtok (NULL,"_-") ;
	aux_i = atoi (pch) ;
	pch = strtok (NULL,"_-") ;
	temp = &pch[3] ;
	aux_t = atoi (temp) ;
	dim = strlen (pch) ;
	if (dim == 5) aux_s = &pch[4] ;
	else aux_s = &pch[5] ;

	for (s = 0; s < n_sites; s++) {

	  if ((strncmp(pch,solv_site[s][0],2) == 0) &&
	      (site_i[s] == aux_i)) {

	    if (strcmp(aux_s,"C") == 0) {
	      
	      solv_p[s][0] = aux ;
	    } else {

	      solv_p[s][aux_t] = aux ;
	      strncpy(write[s][(aux_t - 1)],pch,4) ;
	    }
	    break ;
	  }
	}

      } else {
	
	pch = strtok (ch1,"-") ;
	aux_i = atoi (pch) ;
	pch = strtok (NULL,"-") ;
	temp = &pch[3] ;
	aux_t = atoi (temp) ;
	dim = strlen (pch) ;
	if (dim == 5) aux_s = &pch[4] ;
	else aux_s = &pch[5] ;

	for (s = 0; s < n_sites; s++) {

	  if ((strncmp(pch,solv_site[s][0],2) == 0) &&
	      (site_i[s] == aux_i)) {

	    if (strcmp(aux_s,"C") == 0) {
	      
	      solv_s[s][0] = aux ;
	    } else {

	      solv_s[s][aux_t] = aux ;
	      strncpy(write[s][(aux_t - 1)],pch,4) ;
	    }
	    break ;
	  }
	}
      }
    }
  }

}  /*  read_solv  */



/*********************************************/



/* Reads st files, assigning a site type according to total
 * charge in the protonated state, and the number of protons
 * according to the difference between the protonated and 
 * deprotonated columns. Only reads atom name. The site 
 * type and number of protons are read for every tautomer.  */
void read_sts(int n_sites, char ***site_name, int *n_taut,
	      char **site_type, float **site_prot,
	      char ***st_atom_n, int *st_atoms) {

  FILE *f ;
  char filename[MAXSITENAME], buf[MAXLINESIZE] ;
  int s, t, l, st_line ;
  float tot_p, tot_d, tempP, tempD ;

  tot_p = 0.5 ;

  for (s = 0; s < n_sites; s++) {

/*  Opening first tautomer to know number of atoms and save 
 * that number which is the same in every tautomer of the 
 * same site. */
    strcpy(filename,site_name[s][0]) ;
    strcat(filename,".st") ;
    if ((f = fopen(filename,"r")) == NULL) { 

      error("Can not open st file\n") ;
      exit(1) ;
    }
/*  Get ride of first line .... for now  */
    if (fgets(buf, sizeof(buf), f) != 0)
      {} ;  /*  NULL STATEMENT  */
    l = 0 ;
/*  Count number of atoms  */
    while (fgets(buf, sizeof(buf), f) != NULL) l++ ;
/*  Saving number of atoms for this site  */
    *(st_atoms + s) = l ;
    fclose(f) ;


/*  Going through all tautomers of the present site  */
    for (t = 0; t < *(n_taut + s); t++) {

/*  Defining st file name  */
      strcpy(filename,site_name[s][t]) ;
      strcat(filename,".st") ;
      if ((f = fopen(filename,"r")) == NULL) { 

	error("Can not open st file\n") ;
	exit(1) ;
      }

      tot_p = tot_d = 0.0 ;

/*  Get ride of first line .... for now  */
      if (fgets(buf, sizeof(buf), f) != 0)
	{} ;  /*  NULL STATEMENT  */

/*  Read atom name and find total protonated and 
 * deprotonated charge  */
      for (st_line = 0; st_line < l; st_line++) {

	if (fscanf(f,"%*s%s%f%f",st_atom_n[s][st_line],&tempP,&tempD) == 3) {

	  tot_p += tempP ;
	  tot_d += tempD ;
	}
      }

/*  Determining the difference in total charge to know 
 * how many protons are changing. */
      site_prot[s][t] = (tot_p - tot_d) ;

      fclose(f) ;
    }

/*  Check the type of site (Anionic or Cationic). Site type is 
 * defined by the total charge in the protonated column. 
 * If tot_p > 0.99 site is cationic. This is done after the last
 * tautomer to prevent it to make the same check t times. */
    if (tot_p > 0.99) strcpy(site_type[s],"C") ;
    else if (tot_p < 0.01) strcpy(site_type[s],"A") ;
    else {
      error("Can not assign site type (cationic/anionic)\n") ;
      printf("Check charges on file %s\n",filename) ;
      exit(2) ;
    }
  }
}  /*  read_sts  */



/*********************************************/



/*  Reading pkmod from all tautomers since for the same
 * site, there are cases in which the tautomers have different
 * pkmods (example asp). */
void read_pkmod(int numb_sites, int *numb_taut, char ***site_name,
		float **pkmod) {

  FILE *file ;
  char filename[MAXSITENAME] ;
  int s, t ;


  for (s = 0; s < numb_sites; s++) {
    for (t = 0; t < *(numb_taut + s); t++) {

/*  Defining st file name  */
      strcpy(filename,site_name[s][t]) ;
      strcat(filename,".st") ;
      if ((file = fopen(filename,"r")) == NULL) { 

	error("Can not open st file\n") ;
	exit(1) ;
      }
/*  Read pkmod value from every st/taut file  */

      if (fscanf(file,"%f",&pkmod[s][t]) == 1) 
	{} ;  /*  NULL STATEMENT  */
      fclose(file) ;
    }
  }
}  /*  read_pkmod  */



/*********************************************/



/* Reads specific data from the list of pairwise interactions.
 * Reads col 2 and 3 as pseudo site for writing in .g
 * Reads col 5,6 and 7,8 as site and taut of site 1 and site
 *  2, respectivelly */
void read_pair_inter(FILE *ifile, int n_pair_inter, int *write_si,
		int *write_sj, int *site_si, int *site_sj,
		int *site_ti, int *site_tj) {

  int i ;

  for (i = 0; i < n_pair_inter; i++) {

    if (fscanf(ifile,"%d%d%*s%d%d%*s%d%d",&*(write_si + i),
	       &*(write_sj + i),&*(site_si + i),&*(site_ti + i),
	       &*(site_sj + i),&*(site_tj + i)) == 6)
      { } ;  /*  NULL STATEMENT  */
  }
  return ;
}  /*  read_pair_inter  */



/*********************************************/



/* Reads a list of atoms per pseudosite where the
 * first column is the site index, the second 
 * column is the tautomer index and the rest is 
 * a variable number of atom indexes 
 */
void read_back_inter(FILE *ifile, int n_back_inter, 
		     int *site_bk, int *line_bk, 
		     int **atom_bk ) {

  int wordcount, i, j, eval_ch, keep_ch ;


  for (i = 0; i < n_back_inter; i++) {

    wordcount = eval_ch = 0 ;
    keep_ch = getc(ifile) ;
    while (eval_ch != '\n') {

      eval_ch = getc(ifile) ;
      if (isgraph(eval_ch))
	if (!isgraph(keep_ch))
	  wordcount++ ;
      keep_ch = eval_ch ;
    }
    *(line_bk + i) = wordcount - 1 ;
  }

  rewind(ifile) ;
  for (i = 0; i < n_back_inter; i++) {

    if (fscanf(ifile,"%d",&*(site_bk + i)) == 1) {

      for (j = 0; j < *(line_bk + i); j++) {

	if (fscanf(ifile,"%d",&atom_bk[i][j]) == 1) 
	  {}  /*  NULL STATEMENT  */
      }
    }
  }
  return ;
}  /*  read_back_inter  */



/*********************************************/



/*  Read charges in the pqr file  */
void read_charges(FILE *ifile, int n_atoms, float **chg) {

  int i ;

  for (i = 0; i < n_atoms; i++) {
    /*    for ( j = 1; j < 7; j++) {*/

    if (fscanf(ifile,"%*f%f%f%f%f%f%f",&chg[i][1],&chg[i][2],
	       &chg[i][3],&chg[i][4],&chg[i][5],&chg[i][6]) == 1)
      {} ;  /*  NULL STATEMENT  */
      /*    }*/
  }
  return ;
}  /*  read_charges  */



/*********************************************/



/* Reads frc files ... all of them, only ignoring the first 
 * 12 lines and the last line. It considers different ways of
 * reading the potentials, since as an example, the charge and 
 * potential columns can appear together. 
 * It can read frc files from a focus run, case where it reads 
 * first the focus and then it reads the big box and fills in 
 * values where the focus as potential 0.0 */
void read_frc(int n_sites, int *taut, char ***site_name, 
	      char ***solv_site, int *site_i, int frc_lines, 
	      char ****frc_An, int ***frc_Ri, 
	      float ***frc_Pt) {

  FILE *f ;
  char filename[MAXSITENAME], buf[MAXLINESIZE] ;
  char temp_qPt[MAXSITENAME], temp[MAXSITENAME] ;
  int s, t, l, len ;
  float temp_PT, epsilon = 0.0001 ;


/*  Going through all states  */
  for (s = 0; s < n_sites; s++) {
    for (t = 0; t <= *(taut + s); t++) {

      if (t == *(taut + s)) {
	break ;
      }
/*  Defining frc file name
 * If focus exists, start with it, else read normal box
 */
      if (strcmp(focus,"yes") == 0) {
	if (strcmp(solv_site[s][t],"Nt") == 0) {

	  sprintf(filename,"P_%d-Ntr%dNf.frc",site_i[s],(t+1)) ;
	} else if (strcmp(solv_site[s][t],"Ct") == 0) {

	  sprintf(filename,"P_%d-Ctr%dNf.frc",site_i[s],(t+1)) ;
	} else {

	  sprintf(filename,"P_%d-%s%c%dNf.frc",site_i[s],solv_site[s][t],
		  tolower(site_name[s][t][2]),(t+1)) ;
	}
      } else if (strcmp(focus,"no") == 0) {
	if (strcmp(solv_site[s][t],"Nt") == 0) {

	  sprintf(filename,"P_%d-Ntr%dN.frc",site_i[s],(t+1)) ;
	} else if (strcmp(solv_site[s][t],"Ct") == 0) {

	  sprintf(filename,"P_%d-Ctr%dN.frc",site_i[s],(t+1)) ;
	} else {

	  sprintf(filename,"P_%d-%s%c%dN.frc",site_i[s],solv_site[s][t],
		  tolower(site_name[s][t][2]),(t+1)) ;
	}
      }

      if ((f = fopen(filename,"r")) == NULL) { 

	error("Can not open frc file\n") ;
	exit(1) ;
      }

/*  Do not care about the first 12 lines  */
      for (l = 0; l < 12; l++) {
	if (fgets(buf, sizeof(buf), f) != 0)
	  {} ;  /*  NULL STATEMENT  */
      }


      for (l = 0; l < frc_lines; l++) {

	if (fscanf(f,"%s%*s%d%s",frc_An[s][t+1][l], 
		   &frc_Ri[s][t+1][l],temp_qPt) == 3) {

/*  See if this column is actually 2 or not by the length of its string.
 * If it is larger than 10, then its two columns and I need to split it, 
 * else its only one and I need to read another argument.
 */
	  len = strlen(temp_qPt) ;
	  if (len > 10) {

	    if (strncmp(temp_qPt,"-",1) == 0) {

	      strncpy(temp,temp_qPt,7) ;
	      temp[7] = '\0' ;
/*	      frc_q[s][t+1][l] = atof(temp) ; */
	      strncpy(temp,temp_qPt+7,(len - 7)) ;
	      temp[(len-7)] = '\0' ;
	      frc_Pt[s][t+1][l] = atof(temp) ;
	      
	    } else {

	      strncpy(temp,temp_qPt,6) ;
	      temp[6] = '\0' ;
/*	      frc_q[s][t+1][l] = atof(temp) ;  */
	      strncpy(temp,temp_qPt+6,len) ;
	      frc_Pt[s][t+1][l] = atof(temp) ;
	    }
	  } else {

/*	    frc_q[s][t+1][l] = atof(temp_qPt) ; */
	    if (fscanf(f,"%f",&frc_Pt[s][t+1][l]) == 1) 
	      { } ;  /*  NULL STATEMENT  */

	  }
	}
      }
      fclose (f) ;


/*  Only do this block if we are using focus. Defining frc 
 * file name (Filling zeros from focus with big box values) */
      if (strcmp(focus,"yes") == 0) {
	if (strcmp(solv_site[s][t],"Nt") == 0) {

	  sprintf(filename,"P_%d-Ntr%dNb.frc",site_i[s],(t+1)) ;
	} else if (strcmp(solv_site[s][t],"Ct") == 0) {

	  sprintf(filename,"P_%d-Ctr%dNb.frc",site_i[s],(t+1)) ;
	} else {

	  sprintf(filename,"P_%d-%s%c%dNb.frc",site_i[s],solv_site[s][t],
		  tolower(site_name[s][t][2]),(t+1)) ;
	}

	if ((f = fopen(filename,"r")) == NULL) { 

	  error("Can not open frc file\n") ;
	  exit(1) ;
	}

/*  Do not care about the first 12 lines  */
	for (l = 0; l < 12; l++) {
	  if (fgets(buf, sizeof(buf), f) != 0)
	    {} ;  /*  NULL STATEMENT  */
	}


	for (l = 0; l < frc_lines; l++) {

/*  Only read new Potential value if previously read from focus is zero  */
	  if (fscanf(f,"%*s%*s%*d%s",temp_qPt) == 1) {

/*  See if this column is actually 2 or not by the length of its string.
 * If it is larger than 10, then its two columns and I need to split it, 
 * else its only one and I need to read another argument.
 */
	    len = strlen(temp_qPt) ;
	    if (len > 10) {

	      if (strncmp(temp_qPt,"-",1) == 0) {

		strncpy(temp,temp_qPt,7) ;
		temp[7] = '\0' ;
/*		frc_q[s][t+1][l] = atof(temp) ; */
		strncpy(temp,temp_qPt+7,(len - 7)) ;
		temp[(len-7)] = '\0' ;
		temp_PT = atof(temp) ;
		if (fabs(frc_Pt[s][t+1][l] - 0.0) < epsilon) 
		  frc_Pt[s][t+1][l] = temp_PT ;
	      
	      } else {

		strncpy(temp,temp_qPt,6) ;
		temp[6] = '\0' ;
/*	        frc_q[s][t+1][l] = atof(temp) ; */
		strncpy(temp,temp_qPt+6,len) ;
		temp_PT = atof(temp) ;
		if (fabs(frc_Pt[s][t+1][l] - 0.0) < epsilon) 
		  frc_Pt[s][t+1][l] = temp_PT ;

	      }
	    } else {

/*	      frc_q[s][t+1][l] = atof(temp_qPt) ; */
	      if (fscanf(f,"%f",&temp_PT) == 1) {

		if (fabs(frc_Pt[s][t+1][l] - 0.0) < epsilon) 
		  frc_Pt[s][t+1][l] = temp_PT ;

	      }
	    }
	  }
	}
	fclose (f) ;
      }

    }  /*  Ending tautomer cycle (t)  */

/*  Like above. Defining frc file name for the reference state 
 * If focus exists start with it, else 
 */
    if (strcmp(focus,"yes") == 0) {
      if (strcmp(solv_site[s][0],"Nt") == 0) {

	sprintf(filename,"P_%d-Ntr%dCf.frc",site_i[s],*(taut + s)) ;
      } else if (strcmp(solv_site[s][0],"Ct") == 0) {

	sprintf(filename,"P_%d-Ctr%dCf.frc",site_i[s],*(taut + s)) ;
      } else {

	sprintf(filename,"P_%d-%s%c%dCf.frc",site_i[s],solv_site[s][0],
		tolower(site_name[s][0][2]),*(taut + s)) ;
      }
    } else if (strcmp(focus,"no") == 0) {
      if (strcmp(solv_site[s][0],"Nt") == 0) {

	sprintf(filename,"P_%d-Ntr%dC.frc",site_i[s],*(taut + s)) ;
      } else if (strcmp(solv_site[s][0],"Ct") == 0) {

	sprintf(filename,"P_%d-Ctr%dC.frc",site_i[s],*(taut + s)) ;
      } else {

	sprintf(filename,"P_%d-%s%c%dC.frc",site_i[s],solv_site[s][0],
		tolower(site_name[s][0][2]),*(taut + s)) ;
      }
    }

    if ((f = fopen(filename,"r")) == NULL) { 

      error("Can not open frc file\n") ;
      exit(1) ;
    }

/*  Do not care about the first 12 lines  */
    for (l = 0; l < 12; l++) {
      if (fgets(buf, sizeof(buf), f) != 0) 
	{} ;  /*  NULL STATEMENT  */
    }


    for (l = 0; l < frc_lines; l++) {

      if (fscanf(f,"%s%*s%d%s",frc_An[s][0][l], 
		 &frc_Ri[s][0][l],temp_qPt) == 3) {

/*  See if this column is actually 2 or not by the length of its string.
 * If it is larger than 10, then its two columns and I need to split it, 
 * else its only one and I need to read another argument.
 */
	len = strlen(temp_qPt) ;
	if (len > 10) {

	  if (strncmp(temp_qPt,"-",1) == 0) {

	    strncpy(temp,temp_qPt,7) ;
	    temp[7] = '\0' ;
/*	    frc_q[s][0][l] = atof(temp) ;  */
	    strncpy(temp,temp_qPt+7,(len - 7)) ;
	    temp[(len-7)] = '\0' ;
	    frc_Pt[s][0][l] = atof(temp) ;
	      
	  } else {

	    strncpy(temp,temp_qPt,6) ;
	    temp[6] = '\0' ;
/*	    frc_q[s][0][l] = atof(temp) ;  */
	    strncpy(temp,temp_qPt+6,len) ;
	    frc_Pt[s][0][l] = atof(temp) ;
	  }
	} else {

/*	  frc_q[s][0][l] = atof(temp_qPt) ;  */
	  if (fscanf(f,"%f",&frc_Pt[s][0][l]) == 1) 
	    { } ;  /*  NULL STATEMENT  */

	}
      }
    }
    fclose (f) ;



/*  Only do this block if we are using focus. Defining frc 
 * file name for the reference state 
 * (Filling zeros from focus with big box values) */
    if (strcmp(focus,"yes") == 0) {
      if (strcmp(solv_site[s][0],"Nt") == 0) {

	sprintf(filename,"P_%d-Ntr%dCb.frc",site_i[s],*(taut + s)) ;
      } else if (strcmp(solv_site[s][0],"Ct") == 0) {

	sprintf(filename,"P_%d-Ctr%dCb.frc",site_i[s],*(taut + s)) ;
      } else {

	sprintf(filename,"P_%d-%s%c%dCb.frc",site_i[s],solv_site[s][0],
		tolower(site_name[s][0][2]),*(taut + s)) ;
      }


      if ((f = fopen(filename,"r")) == NULL) { 

	error("Can not open frc file\n") ;
	exit(1) ;
      }

/*  Do not care about the first 12 lines  */
      for (l = 0; l < 12; l++) {
	if (fgets(buf, sizeof(buf), f) != 0)
	  {} ;  /*  NULL STATEMENT  */
      }


      for (l = 0; l < frc_lines; l++) {

/*  Only read new Potential value if previously read from focus is zero  */
	if (fscanf(f,"%*s%*s%*d%s",temp_qPt) == 1) {

/*  See if this column is actually 2 or not by the length of its string.
 * If it is larger than 10, then its two columns and I need to split it, 
 * else its only one and I need to read another argument.
 */
	  len = strlen(temp_qPt) ;
	  if (len > 10) {

	    if (strncmp(temp_qPt,"-",1) == 0) {

	      strncpy(temp,temp_qPt,7) ;
	      temp[7] = '\0' ;
/*	      frc_q[s][0][l] = atof(temp) ; */
	      strncpy(temp,temp_qPt+7,(len - 7)) ;
	      temp[(len-7)] = '\0' ;
	      temp_PT = atof(temp) ;
	      if (fabs(frc_Pt[s][0][l] - 0.0) < epsilon) 
		frc_Pt[s][0][l] = temp_PT ;
	      
	    } else {

	      strncpy(temp,temp_qPt,6) ;
	      temp[6] = '\0' ;
/*	      frc_q[s][0][l] = atof(temp) ; */
	      strncpy(temp,temp_qPt+6,len) ;
	      temp_PT = atof(temp) ;
	      if (fabs(frc_Pt[s][0][l] - 0.0) < epsilon) 
		frc_Pt[s][0][l] = temp_PT ;

	    }
	  } else {

/*	    frc_q[s][0][l] = atof(temp_qPt) ; */
	    if (fscanf(f,"%f",&temp_PT) == 1) {

	      if (fabs(frc_Pt[s][0][l] - 0.0) < epsilon) 
		frc_Pt[s][0][l] = temp_PT ;

	    }
	  }
	}
      }
      fclose (f) ;
    }
      
  }
}  /*  read_frc  */



