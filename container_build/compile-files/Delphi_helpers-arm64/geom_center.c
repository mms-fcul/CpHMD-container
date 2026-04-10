
/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "my_mem.h"

/* Constants */
#define MAXLINESIZE      256
#define MAXSITES        4096
#define MAXPSEUDOSITES    32
#define MAXSITENAME       32
#define MAXSITEATOMS     128
#define MAXATOMNAME        5


/* Memory Functions */
/* In file my_mem.c */



/* Template Functions */
void usage                   (void) ;
void error                   (char *text) ;
int lines_in_file            (FILE *ifile) ;
void read_sites              (FILE *ifile, int n_sites, int *n_taut,
			      int *site_i, char ***site_name, 
			      char ***solv_site) ;
void read_sts                (int n_sites, char ***site_name, 
                              char **site_type, char ***st_atom_n, 
                              int *st_atoms) ;
void read_stsb               (int n_sites, int *n_taut, char ***site_name, 
                              char ****st_atom_nb, int *st_atoms, 
			      int *st_atomsb) ;
void read_gro                (FILE *ifile, int n_atoms, int *gro_resi,
			      char **gro_resn, char **gro_atomn, 
			      float *gro_x, float *gro_y, float *gro_z) ;
float distPBC                 (float a,float b,float c,float u,
			       float v,float w) ;


/*  Global variables  */
float cutoff, vector, half_box ;
int dimension ;


/* Main Program */

int main(int argc, char **argv) {

  FILE *isites, *igro, *ocent, *opairs, *oback ;

  char sites[MAXLINESIZE], gro[MAXLINESIZE], temp[MAXLINESIZE] ;
  char buf[MAXLINESIZE], file_pairs[MAXLINESIZE] ;
  char file_center[MAXLINESIZE], file_back[MAXLINESIZE] ;
  char ***site_name, ***solv_site, **site_type, ***st_atom_n ;
  char **gro_resn, **gro_atomn ;
  char ****st_atom_nb ;

  int s, si, sj, t, ti, tj, n_sites, n_atoms, pt_at, st_at, cnt ;
  int *n_taut, *site_i, *st_atoms, *st_atomsb, *gro_resi, *found ;
  int **corr ;

  float **t_x, **t_y, ** t_z ;
  float *gro_x, *gro_y, *gro_z ;
  float max_x, max_y, max_z, min_x, min_y, min_z, distance ;
  float tempx, tempy, tempz, centX, centY, centZ ;



/* Read Arguments */
  if (argc != 6) {

    error("Wrong number of arguments\n\n") ;
    usage() ; exit(1) ;
  }

  strcpy(sites,argv[1]) ;
  strcpy(gro,argv[2]) ;
  strcpy(temp,argv[3]) ;
  cutoff = atof(temp) ;
  strcpy(temp,argv[4]) ;
  vector = atof(temp) ;
  strcpy(temp,argv[5]) ;
  dimension = atoi(temp) ;


/*  Testing parameters and initialize vectors  */
  if (vector > 0.0) 
    half_box = vector / 2 ;
  if ((cutoff >= half_box) && (dimension > 1)) 
    error("cut-off larger than half box vector\n") ;



/*  Test opening files */
  if ((isites = fopen(sites, "r")) == NULL) {
    error("Can not open sites file\n") ;
    exit(1) ;
  }
  if ((igro = fopen(gro, "r")) == NULL) {
    error("Can not open gro file\n") ;
    exit(1) ;
  }

/*  Setting name for geometric center file  */
  sprintf(file_center,"cent") ;

/*  Setting name for pairwise interaction list file  */
  sprintf(file_pairs,"List_Inter") ;

/*  Setting name for site-atom distance file  */
  sprintf(file_back,"Back_list") ;



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
 *    Reading st files according to sites info
 **********************************************/

/*  Allocate memory for information on st files, namelly 
 * site type, atom name, number of atoms. All info is taken 
 * from the first taut. */
  site_type = alloc_2D_char(n_sites,2) ;
  st_atom_n = alloc_3D_char(n_sites,MAXSITEATOMS,MAXATOMNAME) ;
  st_atoms = calloc(n_sites,sizeof(int)) ;

/*  Read st files and their pKmod's  */
  read_sts(n_sites,site_name,site_type,st_atom_n,st_atoms) ;



/**********************************************
 *    Reading st files but only hydrogens
 **********************************************/

/*  Allocate memory for information on st files, namelly 
 * site type, atom name, number of atoms. All info is taken 
 * from the first taut. */
  st_atom_nb = alloc_4D_char(n_sites,MAXPSEUDOSITES,MAXSITEATOMS,MAXATOMNAME) ;
  st_atomsb = calloc(n_sites,sizeof(int)) ;

/*  Read st files and their pKmod's  */
  read_stsb(n_sites,n_taut,site_name,st_atom_nb,st_atoms,st_atomsb) ;



/**********************************************
 *          Reading gro file
 **********************************************/

/*  Reads the second line in gro file which corresponds 
 * to the number of atoms */
  igro = fopen(gro, "r") ;
  fgets(buf, sizeof(buf), igro) ;
  if (fscanf(igro,"%d",&n_atoms) == 1) 
    {} ;  /*  NULL STATEMENT  */
  fclose(igro) ;

/*  Allocate memory for the data to be read:
 * residue name, number, atom name, x, y and z */
  gro_resi = alloc_1D_int(n_atoms) ;
  gro_resn = alloc_2D_char(n_atoms,MAXATOMNAME) ;
  gro_atomn = alloc_2D_char(n_atoms,MAXATOMNAME) ;
  gro_x = alloc_1D_float(n_atoms) ;
  gro_y = alloc_1D_float(n_atoms) ;
  gro_z = alloc_1D_float(n_atoms) ;

/*  Read gro file. */
  igro = fopen(gro, "r") ;
  read_gro(igro,n_atoms,gro_resi,gro_resn,gro_atomn,gro_x,gro_y,gro_z) ;
  fclose(igro) ;



/***************************************************
 *         Part I. geometric center
 *
 *   Calculating the geometric center for each site
 *  Do not care about tautomers (same for all)
 *
 ***************************************************/

/*  Opening file for geometric center calculation  */
  ocent = fopen(file_center, "w") ;

  for (s = 0; s < n_sites; s++) { /* sites */

    max_x = max_y = max_z = -1000.0 ;
    min_x = min_y = min_z = 1000.0 ;  
    for (pt_at = 0; pt_at < n_atoms; pt_at++) {  /* Protein atoms */

      for (st_at = 0; st_at < *(st_atoms + s); st_at++) {  /* Site atoms */

	if ((site_i[s] == gro_resi[pt_at]) && 
	    (strcmp(st_atom_n[s][st_at],gro_atomn[pt_at]) == 0)) {

	  if (gro_x[pt_at] < min_x) min_x = gro_x[pt_at] ;
	  if (gro_x[pt_at] > max_x) max_x = gro_x[pt_at] ;
	  if (gro_y[pt_at] < min_y) min_y = gro_y[pt_at] ;
	  if (gro_y[pt_at] > max_y) max_y = gro_y[pt_at] ;
	  if (gro_z[pt_at] < min_z) min_z = gro_z[pt_at] ;
	  if (gro_z[pt_at] > max_z) max_z = gro_z[pt_at] ;

	  if (st_at == (*(st_atoms + s) - 1))

/*  Center is determined by finding the max and min in the 
 * three axes and then average */
	    fprintf(ocent,"%4d %12.5f %12.5f %12.5f\n",site_i[s],
		    (max_x + min_x) / 2 ,
		    (max_y + min_y) / 2 ,
		    (max_z + min_z) / 2 ) ;
	  
	}
      }  /*  end site atoms  */
    }  /*  end protein atoms  */
  }  /*  end sites  */


/***************************************************
 *         Part II. Search pairwise groups
 * 
 * Getting coordinates for distance calculation:
 *  + Anionic -> use coordinates from charged hydrogen
 *  + Cationic -> use average coordinate of all hydrogens
 ****************************************************/


  t_x = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;
  t_y = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;
  t_z = alloc_2D_float(n_sites,MAXPSEUDOSITES) ;
  corr = alloc_2D_int(n_sites,MAXPSEUDOSITES) ;

  for (s = 0; s < n_sites; s++) { /* sites */

    tempx = tempy = tempz = 0.0 ;
    for (pt_at = 0; pt_at < n_atoms; pt_at++) {  /* Protein atoms */

      if (site_i[s] == gro_resi[pt_at]) {

	for (st_at = 0; st_at < *(st_atomsb + s); st_at++) {  /* Site atoms */
	  if (strcmp(st_atom_nb[s][0][st_at],gro_atomn[pt_at]) == 0) {

	    tempx += gro_x[pt_at] ;
	    tempy += gro_y[pt_at] ;
	    tempz += gro_z[pt_at] ;

	    if (st_at == (*(st_atomsb + s) - 1)) {

	      centX = tempx / *(st_atomsb + s) ;
	      centY = tempy / *(st_atomsb + s) ;
	      centZ = tempz / *(st_atomsb + s) ;

	      for (t = 0; t < *(n_taut + s); t++) {

		t_x[s][t] = centX ;
		t_y[s][t] = centY ;
		t_z[s][t] = centZ ;
	      }
	    }
	  }
	}  /*  st_at  */
      }
    }  /*  pt_at  */
  }  /*  sites  */

/*  Correcpondence from site with tauts to total 
 * number of pseudo sites */
  cnt = 0 ;
  for (s = 0; s < n_sites; s++) {  /*  sites  */
    for (t = 0; t < *(n_taut + s); t++) {  /*  tautomers  */

      cnt++ ;
      corr[s][t] = cnt ;
    }
  }

/***************************************************
 * Measure all distances between all hydrogens from
 * all sites. The treatment is different for cationic
 * and anionic.
 ***************************************************/
/*  Opening file for pairwise interaction list  */
  opairs = fopen(file_pairs, "w") ;

  for (si = 0 ; si < n_sites; si++) {
    for (sj = 0; sj < n_sites; sj++) {
      for (ti = 0; ti < *(n_taut + si); ti++) {  
	for (tj = 0; tj < *(n_taut + sj); tj++) {

	  if (si == sj) {
	    if (ti == tj)

	      fprintf(opairs,"%5d %5d %6.1f \t%5d %-3d - %5d %-3d\n",
		      corr[si][ti],corr[sj][tj],0.0,site_i[si],(ti+1),
		      site_i[sj],(tj+1)) ;
	    else

	      fprintf(opairs,"%5d %5d %6.1f \t%5d %-3d - %5d %-3d\n",
		      corr[si][ti],corr[sj][tj],1.0,site_i[si],(ti+1),
		      site_i[sj],(tj+1)) ;
	  } else {

/*  The distance calculated below is a changed distance which
 * calculates the smallest distance between two objects considering
 * periodic boundary conditions */

	    if (cutoff < 0.0) {

	      fprintf(opairs,"%5d %5d %6s \t%5d %-3d - %5d %-3d\n",
		      corr[si][ti],corr[sj][tj],"x",site_i[si],(ti+1),
		      site_i[sj],(tj+1)) ;
	    } else {	    

	      distance = distPBC(t_x[si][ti],t_y[si][ti],t_z[si][ti],	\
				 t_x[sj][tj],t_y[sj][tj],t_z[sj][tj]) ;

	      if (distance <= cutoff) 
		fprintf(opairs,"%5d %5d %6s \t%5d %-3d - %5d %-3d\n",
			corr[si][ti],corr[sj][tj],"x",site_i[si],(ti+1),
			site_i[sj],(tj+1)) ;
	      else
		fprintf(opairs,"%5d %5d %6.1f \t%5d %-3d - %5d %-3d\n",
			corr[si][ti],corr[sj][tj],0.0,site_i[si],(ti+1),
			site_i[sj],(tj+1)) ;
	    }
	  }
	}  /*  tautomer cycle (tj)  */
      }  /*  tautomer cycle (ti)  */
    }  /*  site cycle (sj)  */
  }  /*  site cycle (si)  */


/***************************************************
 *   Part III. Gen distance matrix for background
 * 
 *  Going through each site and using the above 
 * tautomer center, find all protein atoms inside
 * a given cutoff (same as used for pairwise
 * interactions). Atoms from the current site are
 * excluded 
 ***************************************************/

/*  Opening file for background interaction list  */
  oback = fopen(file_back, "w") ;

/* Allocating memmory */
  found = calloc(n_atoms,sizeof(int)) ;

  for (s = 0; s < n_sites; s++) {  /*  sites  */

    fprintf(oback,"%3d",s) ;
    for (pt_at = 0; pt_at < n_atoms; pt_at++) {  /* Protein atoms */

      found[pt_at] = 0 ;
      for (st_at = 0; st_at < *(st_atoms + s); st_at++) {  /* Site atoms */

/*  If atom belongs to own site, is excluded 
 * (in this case, found means it is not accounted)
 */
	if ((site_i[s] == gro_resi[pt_at]) &&
	    (strcmp(st_atom_n[s][st_at],gro_atomn[pt_at]) == 0)) {

	  found[pt_at] = 1 ;
	  break ;
	}
      }  /* Site atoms */

/*  If atom does not belong to own site, it is 
 * accounted if it is inside cutoff.
 */
      if (found[pt_at] == 0) {

	if (cutoff < 0.0) {

	  fprintf(oback," %d",pt_at) ;
	} else {
	  
	  distance = distPBC(t_x[s][0],t_y[s][0],t_z[s][0],		\
			     gro_x[pt_at],gro_y[pt_at],gro_z[pt_at]) ;

	  if (distance <= cutoff) 
	    fprintf(oback," %d",pt_at) ;
	}
      }

    }  /*  Protein atoms  */
    fprintf(oback,"\n") ;
  }  /*  sites  */


  return 0 ;
}



/*********************************************/
/*********************************************/
/*********************************************/
/************** FUNCTIONS ********************/
/*********************************************/
/*********************************************/
/*********************************************/



/* Show usage */
void usage(void) {

  fprintf(stderr, "Usage: geom_center <sites> <gro> <cutoff> <vector> <dim>\n") ;
  fprintf(stderr, " * sites is a simple MEAD sites file\n"
	  " * gro is a typical gro file to find the geometric center\n"
	  "    of each site.\n * cutoff is the radii inside which background"
	  " and pairwise\n    interactions are going to be calculated (-1 for") ;
  fprintf(stderr, " all atoms).\n * vector is the box size in the X axis (-1 to "
	  "be ignored).\n * dim is the number of dimensions in calculation (2 or "
	  "3\n    for bilayers; 0 for proteins or peptides alone).\n\n") ;

  return ;
}  /*  usage  */



/*********************************************/



void error(char *text) {

  fprintf(stderr, "\ngeom_center: ERROR: %s",text) ;
  return ;
}  /*  error  */



/*********************************************/



/*  Reads the number of lines in a file */
int lines_in_file(FILE *ifile) {

  char buf[MAXLINESIZE] ;
  int line_num = 0 ;

  while (fgets(buf, sizeof(buf), ifile) != NULL) {

    line_num++ ;
  }
  return line_num ;
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



/* Reads st files, assigning a site type according to total
 * charge in the protonated state. Only reads atom name.
 * The site type and atom name are only read for the first
 * tautomer.  */
void read_sts(int n_sites, char ***site_name, 
	      char **site_type, char ***st_atom_n, 
	      int *st_atoms) {

  FILE *f ;
  char filename[MAXSITENAME], buf[MAXLINESIZE] ;
  int s, l, st_line ;
  float tot_p, temp ;

  for (s = 0; s < n_sites; s++) {

/*  Defining st file name  */
    strcpy(filename,site_name[s][0]) ;
    strcat(filename,".st") ;
    if ((f = fopen(filename,"r")) == NULL) { 

      error("Can not open st file\n") ;
      exit(1) ;
    }

/*  Get ride of first line .... for now  */
    if (fgets(buf, sizeof(buf), f) != 0)
      {} ;  /*  NULL STATEMENT  */

    l = 0 ; tot_p = 0.0 ;

/*  Count number of atoms  */
    while (fgets(buf, sizeof(buf), f) != NULL) l++ ;

    rewind(f) ;
/*  Get ride of first line .... for now  */
    if (fgets(buf, sizeof(buf), f) != 0)
      {} ;  /*  NULL STATEMENT  */

/*  Read atom name and find total protonated charge  */
    for (st_line = 0; st_line < l; st_line++) {

      if (fscanf(f,"%*s%s%f%*f",st_atom_n[s][st_line],&temp) == 2) 
	tot_p += temp ;
    }

/*  Check type of site (Anionic or Cationic)  */    
    if (tot_p > 0.99) strcpy(site_type[s],"C") ;
    else strcpy(site_type[s],"A") ;

/*  Saving number of atoms for this site  */
    *(st_atoms + s) = l ;

    fclose(f) ;

  }
}  /*  read_sts  */



/* Reads st files, according  a site type according to total
 * charge in the protonated state. Only reads atom name.
 * The site type and atom name are only read for the first
 * tautomer.  */
void read_stsb(int n_sites, int *taut, char ***site_name, 
	       char ****st_atom_nb, int *st_atoms,
	       int *st_atomsb) {

  FILE *f ;
  char filename[MAXSITENAME], buf[MAXLINESIZE], name[MAXLINESIZE] ;
  int s, t, st_line, line_h = 0 ;

  for (s = 0; s < n_sites; s++) {

    for (t = 0; t < *(taut + s) ; t++) {

/*  Defining st file name  */
      strcpy(filename,site_name[s][t]) ;
      strcat(filename,".st") ;
      if ((f = fopen(filename,"r")) == NULL) { 

	error("Can not open st file\n") ;
	exit(1) ;
      }

/*  Get ride of first line  */
      if (fgets(buf, sizeof(buf), f) != 0)
	{} ;  /*  NULL STATEMENT  */

      line_h = 0 ;

/*  Read hydrogen name if site anionic, else read all hydrogens  */
      for (st_line = 0; st_line < *(st_atoms + s); st_line++) {

	if (fscanf(f,"%*s%s%*f%*f",name) == 1) 
	  {} ;  /*  NULL STATEMENT  */

	if (strncmp(name,"H",1) == 0) {

	  strcpy(st_atom_nb[s][t][line_h],name) ;
	  line_h++ ;
	}
      }
      fclose(f) ;

    }  /*  tautomer cycle  */

/*  Saving number of atoms for this site  */
    *(st_atomsb + s) = line_h ;


  }  /*  sites cycle  */
}  /*  read_stsb  */



/*********************************************/



/*  Read charges in the pqr file  */
void read_gro(FILE *ifile, int n_atoms, int *gro_resi,
	      char **gro_resn, char **gro_atomn, 
	      float *gro_x, float *gro_y, float *gro_z) {

  char buf[MAXLINESIZE], temp[MAXSITENAME] ;
  int i, j, k ;


  fgets(buf, sizeof(buf), ifile) ;
  fgets(buf, sizeof(buf), ifile) ;

/*  Read line by line each field and 
 *    dump velocities if they exist
 */
  for (i = 0; i < n_atoms; i++) {

    if (fgets(buf, sizeof(buf), ifile) != 0) {

      strncpy(temp, buf, 5) ; temp[5] = '\0' ;
      gro_resi[i] = atoi(temp) ;
      strncpy(gro_resn[i], buf + 5, 3) ;
      gro_resn[i][3] = '\0' ;
      strncpy(temp, buf + 11, 4) ; temp[4] = '\0' ; 

      for (j = 0, k = 0; j < (int)strlen(temp); j++, k++) {
	if (temp[j] != ' ') 
            gro_atomn[i][k] = temp[j];                     
	else 
	  k-- ;                                     
      }
      gro_atomn[i][k] = '\0' ;

      strncpy(temp, buf + 20, 8) ; temp[8] = '\0' ;
      gro_x[i] = atof(temp) ;
      strncpy(temp, buf + 28, 8) ; temp[8] = '\0' ;
      gro_y[i] = atof(temp) ;
      strncpy(temp, buf + 36, 8) ; temp[8] = '\0' ;
      gro_z[i] = atof(temp) ;
    }
  }
  return ;
}  /*  read_gro  */



/*********************************************/



/* Calculates the distance between two objects using
 * periodic boundary conditions */
float distPBC(float a,float b, float c,float u,
	      float v,float w) {

  float tp_x, tp_y, tp_z, dx, dy, dz ;


  tp_x = fabs(a - u) ; tp_y = fabs(b - v) ;
  if (tp_x > half_box) dx = vector - tp_x ;
  else dx = tp_x ;
  if (tp_y > half_box) dy = vector - tp_y ;
  else dy = tp_y ;
  if (dimension == 2) dz = c - w ;
  else {
    tp_z = fabs(c - w) ;
    if (tp_z > half_box) dz = vector - tp_z ;
    else dz = tp_z ;
  }

  return (sqrt((dx)*(dx)+(dy)*(dy)+(dz)*(dz))) ;
}  /*  distPBC  */



