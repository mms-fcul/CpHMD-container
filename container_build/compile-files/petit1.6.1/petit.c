/*
This file is part of PETIT, version 1.6.1.

Copyright (c) 2001-2012, Antonio M. Baptista, Instituto de Tecnologia
Quimica e Biologica, Universidade Nova de Lisboa, Portugal.

PETIT is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 2 of the License, or (at your
option) any later version.

PETIT is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with PETIT.  If not, see <http://www.gnu.org/licenses/>.

For further details and info check the README file.

You can get PETIT at <http://www.itqb.unl.pt/simulation>.
*/


#define VERSION "1.6.1"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>

#define STRSIZE 1000
#define LINESIZE 10000
#define MAXDELTA 50
#define MAXNPKHALFS 5
#define LN10 2.302585092994

void parse_arguments(int argc, char **argv) ;
void usage(void) ;
void read_data(void) ;
void make_initializations(void) ;
void select_pairs(void) ;
void initialize_EpH_point(float E, float pH) ;
void initialize_E_line(float E) ;
void mc_step(void) ;
int  accept_move(float dU) ;
void compute_statistics(int t) ;
void write_EpH_point(float E, float pH) ;
void write_E_line(float E) ;
float invexp(float x) ;
void clean_memory(void) ;
double sqrtp(double x) ;
void strsplit(char *s, const char *delim) ;
void error(char errtype, char *format, ...) ;


int nsites, npairs, pHsteps, Esteps, taumax, nset, nmicro,
    compute_set, compute_pair_corr, compute_errors,
    compute_energetics, nsitesP, nsitesR, nsubstr ;
int  *nstates, *occ, **stateocc, *state, *pair1, *pair2, *titration,
     *npkhalfs, **cf, **cc, **buf, *setsite, *binP, *binR ;
long mcsteps, eqsteps, seed, avgP, avgR, avgP2, avgR2, avgPR ;
long  *avg, **subavg, **avgpair, *avgmicro ;
float pHmin, pHmax, dpH, Emin, Emax, dE, T, couple_min, min_corr, energy0 ;
float **pkhalf, *pmean, ****gg, **g, **u ;
double avgU, avgU2, *avgUn ;
char cmd[STRSIZE], **site_name, *site_type, sset[STRSIZE], **substr ;

float kBoltz_au  = 5.98435e-6 ;  /*  e^2/(Angstrom*K)  */
float kBoltz_meV = 0.0861734 ;   /*  meV/K  */


int main(int argc, char **argv)
{
  int i, j, t ;
  float pH, E ;

  parse_arguments(argc, argv) ;
  read_data() ;
  select_pairs() ;

  for (i = 0 ; i < Esteps ; i++)
  {
    E = Emin + i * dE ;
    initialize_E_line(E) ;
    for (j = 0 ; j < pHsteps ; j++)
    {
      pH = pHmin + j * dpH ;
      initialize_EpH_point(E, pH) ;
      for (t = 0 ; t < eqsteps ; t++) mc_step() ;
      for (t = 0 ; t < mcsteps ; t++)
      {
        mc_step() ;
        compute_statistics(t) ;
      }
      write_EpH_point(E, pH) ;
    }
    write_E_line(E) ;
  }

  printf("f") ;
  for (i = 0 ; i < nsites ; i++) printf(" %1d", state[i]) ;
  printf("\n#\n") ;

  /*  clean_memory() ;  */
  return 0 ;
}


void parse_arguments(int argc, char **argv)
{
  int c ;

  strcpy(cmd, argv[0]) ;

  T = 300.0 ;
  couple_min = 2.0 ;
  pHmin =  -5.0 ;
  pHmax =  20.0 ;
  dpH   =   0.5 ;
  /* E is actually F*E, given in meV */
  Emin =  0.0 ;
  Emax =  0.0 ;
  dE   =  1.0 ;
  seed = 1234567 ;
  eqsteps = 1000 ;
  taumax = 20 ;
  nset = 0 ;
  compute_set = 0 ;
  compute_pair_corr = 0 ;
  compute_errors = 0 ;

  while ((c = getopt(argc, argv, "H:E:T:c:q:r:p:s:t:evd")) != -1)
  {
    switch(c)
    {
    case 'H':
      strsplit(optarg, ",") ;
      if (nsubstr != 3) error('U', "Wrong argument for option -H.\n") ;
      pHmin = atof(substr[0]) ;
      pHmax = atof(substr[1]) ;
      dpH = atof(substr[2]) ;
      break ;
    case 'E':
      strsplit(optarg, ",") ;
      if (nsubstr != 3) error('U', "Wrong argument for option -E.\n") ;
      Emin = atof(substr[0]) ;
      Emax = atof(substr[1]) ;
      dE = atof(substr[2]) ;
      break ;
    case 'T':
      T = atof(optarg) ;
      break ;
    case 'c':
      couple_min = atof(optarg) ;
      break ;
    case 'q':
      eqsteps = atol(optarg) ;
      break ;
    case 'r':
      seed = atol(optarg) ;
      break ;
    case 'p':
      compute_pair_corr = 1 ;
      min_corr = atof(optarg) ;
      break ;
    case 's':
      compute_set = 1 ;
      strcpy(sset, optarg) ;
      strsplit(sset, ",") ;
      if ((nset = nsubstr) < 2)
        error('U', "Option -s requires a minimum of 2 sites.\n") ;
      break ;
    case 't':
      compute_errors = 1 ;
      taumax = atoi(optarg) ;
      break ;
    case 'e':
      compute_energetics = 1 ;
      break ;
    case 'v':
      fprintf(stderr, "Version: %s\n", VERSION) ;
      exit(0) ;
      break ;
    case 'd':
      fprintf(stderr, "Defaults:\n") ;
      fprintf(stderr, "  -H %f,%f,%f   -E %f,%f,%f\n",
              pHmin, pHmax, dpH, Emin, Emax, dE) ;
      fprintf(stderr, "  -T %f   -c %f   -q %ld   -r %ld\n",
              T, couple_min, eqsteps, seed) ;
      exit(0) ;
      break ;
    default:
      usage() ;
      exit(1) ;
    }
  }
  if (argc - optind != 1) error('U', "Wrong number of arguments\n") ;
  mcsteps = atol(argv[optind]) ;
}


void usage(void)
{
  fprintf(stderr, "Usage: %s [options] MCsteps  < input_file  > ...\n", cmd) ;
  fprintf(stderr, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
    "Options:\n",
    "  -H pHmin,pHmax,dpH\t : solution pH range and increment\n",
    "  -E Emin,Emax,dE   \t : solution potential range and increment (meV)\n",
    "  -T temperature    \t : temperature (Kelvin) \n",
    "  -c couple_min     \t : couple threshold for double moves (pK units) \n",
    "                    \t   (full calculation if zero)\n",
    "  -q eqsteps        \t : number of equilibration steps\n",
    "  -r seed           \t : seed for random numbers\n",
    "  -p min_corr       \t : cutoff for printing pair correlations\n",
    "  -e                \t : compute energetics\n",
    "  -s site1,site2,...\t : set of sites for microstate statistics\n",
    "                    \t   (whose calculation is switched on)\n",
    "  -t taumax         \t : maximum correlation time\n",
    "                    \t   (switches on calculation of errors)\n",
    "  -d                \t : shows defaults (if alone)\n",
    "  -v                \t : shows version\n") ;
}


void read_data(void)
{
  int i, j, si, sj ;
  float auxfloat ;
  char stit[STRSIZE] ;

  /* Read number of sites */
  scanf("%d\n", &nsites) ;

  make_initializations() ;
  g = calloc(nsites, sizeof(float *)) ;

  /* Read individual info */
  nsitesP = nsitesR = 0 ;
  for (i = 0 ; i < nsites ; i++)
  {
    scanf("%s %d %c %s\n", site_name[i], &nstates[i], &site_type[i], stit) ;

    if (site_type[i] != 'P' && site_type[i] != 'R')
      error('E', "Site %s has wrong ligand type.\n", site_name[i]) ;
    if (site_type[i] == 'P') nsitesP++ ;
    else nsitesR++ ;

    switch(stit[0])
    {
    case '*':  /* titrable through single and double moves */
      titration[i] = 3 ;
      break ;
    case 'd':  /* titrable through double moves only */
      titration[i] = 2 ;
      break ;
    case 's':  /* titrable through single moves only */
      titration[i] = 1 ;
      break ;
    case 'n':  /* non-titrable but multistate */
      titration[i] = 4 ;
      break ;
    default:
      j = atoi(stit) ;
      if (j >= 0 && j < nstates[i])   /* fixed state */
      {
	titration[i] = 0 ;
	state[i] = j ;
	break ;
      }
      else error('E', "Pseudo-site %s has wrong titration mode %c.\n",
		 site_name[i], stit) ;
    }
    if (titration[i] != 0) state[i] = 0 ;

    g[i] = calloc(nstates[i], sizeof(float)) ;
    stateocc[i] = calloc(nstates[i], sizeof(int)) ;

    for (si = 0 ; si < nstates[i] ; si++)
    {
      scanf("%d %f\n", &stateocc[i][si], &auxfloat) ;
      g[i][si] = auxfloat / (kBoltz_au * T) ;
    }
  }

  if (nsitesR == 0 && Esteps > 1)
    error('E', "Range of E does not make sense without redox sites.\n") ;

  u = calloc(nsites, sizeof(float *)) ;
  for (i = 0 ; i < nsites ; i++) u[i] = calloc(nstates[i], sizeof(float)) ;
  subavg = calloc(nsites, sizeof(long *)) ;
  for (i = 0 ; i < nsites ; i++) subavg[i] = calloc(nstates[i], sizeof(long)) ;
  gg = calloc(nsites, sizeof(float ***)) ;
  for (i = 0 ; i < nsites ; i++)
  {
    gg[i] = calloc(nstates[i], sizeof(float **)) ;
    for (si = 0 ; si < nstates[i] ; si++)
    {
      gg[i][si] = calloc(nsites, sizeof(float *)) ;
      for (j = 0 ; j < nsites ; j++)
        gg[i][si][j] = calloc(nstates[j], sizeof(float)) ;
    }
  }

  /* Read pairwise info */
  for (i = 0 ; i < nsites - 1 ; i++)
  for (j = i + 1 ; j < nsites ; j++)
  {
    for (si = 0 ; si < nstates[i] ; si++)
    for (sj = 0 ; sj < nstates[j] ; sj++)
    {
      scanf("%*d %*d %*d %*d %f\n", &auxfloat) ;
      gg[i][si][j][sj] = gg[j][sj][i][si] = auxfloat / (kBoltz_au * T) ;
    }
  }
}


void make_initializations(void)
{
  int i ;

  pHsteps = rint(1 + (pHmax - pHmin) / dpH) ;
  Esteps  = rint(1 + (Emax - Emin) / dE) ;
  srand48(seed) ;

  site_name = calloc(nsites, sizeof(char *)) ;
  for (i = 0 ; i < nsites ; i++)
    site_name[i] = calloc(STRSIZE, sizeof(char)) ;
  site_type = calloc(nsites, sizeof(char)) ;
  titration = calloc(nsites, sizeof(int)) ;
  state = calloc(nsites, sizeof(int)) ;
  stateocc = calloc(nsites, sizeof(int *)) ;

  nstates = calloc(nsites, sizeof(int)) ;
  occ = calloc(nsites, sizeof(int)) ;
  pair1 = calloc(nsites * (nsites - 1) / 2, sizeof(int)) ;
  pair2 = calloc(nsites * (nsites - 1) / 2, sizeof(int)) ;
  avg = calloc(nsites, sizeof(long)) ;
  pmean = calloc(nsites, sizeof(float)) ;
  binP = calloc(nsites + 1, sizeof(int)) ;
  binR = calloc(nsites + 1, sizeof(int)) ;
  pkhalf = calloc(nsites, sizeof(float *)) ;
  for (i = 0 ; i < nsites ; i++)
    pkhalf[i] = calloc(MAXNPKHALFS, sizeof(float)) ;   
  npkhalfs = calloc(nsites, sizeof(int *)) ;
  if (compute_pair_corr)
  {
    avgpair = calloc(nsites, sizeof(long *)) ;
    for (i = 0 ; i < nsites ; i++) avgpair[i] = calloc(nsites, sizeof(long)) ;
  }
  if (compute_energetics) avgUn = calloc(nsites, sizeof(double)) ;
  if (compute_set)
  {
    setsite = calloc(nset, sizeof(int)) ;
    strsplit(sset, ",") ;
    for (i = 0 ; i < nset ; i++)
    {
      setsite[i] = atoi(substr[i]) ;
      if (setsite[i] < 0 || setsite[i] >= nsites)
        error('E', "Site number given to -s is out of range.\n") ;
    }
    nmicro = 1 << nset ;
    avgmicro = calloc(nmicro, sizeof(long)) ;
  }
  if (compute_errors)
  {
    cf =  calloc(nsites, sizeof(int *)) ;
    for (i = 0 ; i < nsites ; i++) cf[i] = calloc(taumax, sizeof(int)) ;
    cc = calloc(nsites, sizeof(int *)) ;
    for (i = 0 ; i < nsites ; i++) cc[i] = calloc(taumax, sizeof(int)) ;
    buf = calloc(nsites, sizeof(int *)) ;
    for (i = 0 ; i < nsites ; i++) buf[i] = calloc(taumax, sizeof(int)) ;
  }

  /* Write header */
  printf("# PETIT: Proton and Electron TITration \n") ;
  printf("# version: %s \n#\n", VERSION) ;
  printf("# Total number of sites = %d\n", nsites) ;
  printf("# pH (min,max,delta): %f,%f,%f\n", pHmin, pHmax, dpH) ;
  printf("# E  (min,max,delta): %f,%f,%f\n", Emin, Emax, dE) ;
  printf("# Temperature = %f K\n", T) ;
  printf("# Production MC steps per (E,pH) point = %ld\n", mcsteps) ;
  printf("# Equilibration MC steps per (E,pH) point = %ld\n", eqsteps) ;
  printf("# Couple treshold for double moves = %f pK units\n",
         couple_min) ;
  if (compute_pair_corr)
    printf("# Site-site correlations r computed and written if |r| >= %f.\n",
           min_corr) ;
  else printf("# Computation of site-site correlations switched off.\n") ;
  if (compute_set)
  {
    printf("# Set of sites selected for microstate statistics: %d",
           setsite[0]) ;
    for (i = 1 ; i < nset ; i++) printf(",%d", setsite[i]) ;
    printf("\n") ;
  }
  else printf("# Computation of microstate statistics switched off.\n") ;
  if (compute_errors)
    printf("# Maximum correlation time used for error calculation = %d\n",
           taumax) ;
  else printf("# Error calculation switched off.\n") ;
  if (compute_energetics)
    printf("# Energetics calculations switched on.\n") ;
  else printf("# Energetics calculations switched off.\n") ;
  printf("# Seed for random number generator = %ld\n", seed) ;
  printf("#\n") ;
}


void select_pairs(void)
{
  int i, j, si, sj, coup ;
  float ggmax ;

  printf("### Coupled site pairs, with max|gg| >= %f pK units:\n",
         couple_min) ;
  npairs = 0 ;
  for (i = 0 ; i < nsites - 1 ; i++)
  for (j = i + 1 ; j < nsites ; j++)
  {
    coup = 0 ;
    ggmax = 0 ;
    for (si = 0 ; si < nstates[i] ; si++)
    for (sj = 0 ; sj < nstates[j] ; sj++)
    {
      if (titration[i] >= 2 && titration[j] >= 2)
        if (fabs(gg[i][si][j][sj]) >= LN10 * couple_min)
        {
          coup = 1 ;
          if (fabs(gg[i][si][j][sj]) > ggmax) ggmax = gg[i][si][j][sj] ;
        }
    }

    if (coup)
    {
      pair1[npairs] = i ;
      pair2[npairs] = j ;
      printf("###   %3d %3d  %10.6f\n", i, j, ggmax / LN10) ;
      npairs++ ;
    }
  }
  printf("### Total number of coupled site pairs = %d\n#\n", npairs) ;
}


void initialize_EpH_point(float E, float pH)
{
  int i, si ;
  float pLi ;

  for (i = 0 ; i < nsites ; i++)
  {
    if (site_type[i] == 'P') pLi = pH ;
    /* - 2.3 kT pLi <-> E (where E is really F*E) */
    else pLi = - E / (LN10 * kBoltz_meV * T) ;
    for (si = 0 ; si < nstates[i] ; si++)
      u[i][si] = LN10 * pLi * stateocc[i][si] + g[i][si] ;
  }
}


void initialize_E_line(float E)
{
  int i ;

  for (i = 0 ; i < nsites ; i++) npkhalfs[i] = 0 ;
}


void mc_step(void)
{
  int i, j, k, si, sinew, sj, sjnew, sk, p ;
  float dU, **gi, **ginew, **gj, **gjnew ;

  for (i = 0 ; i < nsites ; i++)
  {
    if (titration[i] == 0) continue ;
    si = state[i] ;
    sinew = (int) (drand48() * nstates[i]) ;
    gi = gg[i][si] ;
    ginew = gg[i][sinew] ;
    dU = u[i][sinew] - u[i][si] ;
    for (j = 0 ; j < nsites ; j++)
      if (j != i)
      {
        sj = state[j] ;
        dU += ginew[j][sj] - gi[j][sj] ;
      }
    if (accept_move(dU)) state[i] = sinew ;
  }

  for (p = 0 ; p < npairs ; p++)
  {
    i = pair1[p] ;
    j = pair2[p] ;
    if (titration[i] == 0 || titration[j] == 0) continue ;
    si = state[i] ;
    sj = state[j] ;
    sinew = (int) (drand48() * nstates[i]) ;
    sjnew = (int) (drand48() * nstates[j]) ;
    gi = gg[i][si] ;
    ginew = gg[i][sinew] ;
    gj = gg[j][sj] ;
    gjnew = gg[j][sjnew] ;
    dU = u[i][sinew] - u[i][si] + u[j][sjnew] - u[j][sj]
         + ginew[j][sjnew] - gi[j][sj] ;
    for (k = 0 ; k < nsites ; k++)
      if (k != i && k != j)
      {
        sk = state[k] ;
        dU += ginew[k][sk] - gi[k][sk] + gjnew[k][sk] - gj[k][sk] ;
      }
    if (accept_move(dU))
    {
      state[i] = sinew ;
      state[j] = sjnew ;
    }
  }

}


int accept_move(float dU)
{
  if (dU > MAXDELTA) return 0 ;
  if (dU <= 0) return 1 ;
  if (invexp(dU) > drand48()) return 1 ;
  else return 0 ;
}


/* This function is a piecewise rational polynomial approximation of
   exp(-x) in the intervals [0,6], [6,13] and [13,70]. Thus, it may
   give significant (even drastic) errors outside the range [0,70]. */
float invexp(float x)
{
  float x2, x3, x4 ;

  if (x > 13)
  {
    return 8.194236147130614e-10 - 1.3290994520804703e-11 * x ;
  }
  else
  {
    x2 = x * x ;
    x3 = x2 * x ;
    if (x > 6)
    {
      return (-0.0013245823657199278 + 0.00027464252539452071 * x -
              0.000019314947607346905 * x2 + 4.598224667374957e-7 * x3) /
             (1.0 - 0.5165170691890946 * x + 0.09211442135429947 * x2 -
              0.006143102546214945 * x3) ;
    }
    else
    {
      x4 = x2 * x2 ;
      return (0.9999965470613797 - 0.3960827416191208 * x +
              0.06303500815508939 * x2 - 0.00476617578304489 * x3 +
              0.00014392025197088043 * x4)/
             (1.0 + 0.6038220689877429 * x + 0.16732494517488303 * x2 +
              0.026354026827091058 * x3 + 0.00289071552898347 * x4) ;
    }
  }
}


void compute_statistics(int t)
{
  int i, nP, nR, s, occi ;
  int j, tt, tau0, tau ;
  float energy ;

  if (t == 0)
  {
    avgP = avgR = avgP2 = avgR2 = avgPR = 0 ;
    if (compute_energetics) avgU = avgU2 = 0.0 ;
    for (i = 0 ; i < nsites ; i++)
    {
      avg[i] = 0 ;
      for (s = 0 ; s < nstates[i] ; s++) subavg[i][s] = 0 ;
      if (compute_energetics) avgUn[i] = 0.0 ;
      if (compute_errors)
        for (tau = 0 ; tau < taumax ; tau++)
          cf[i][tau] = cc[i][tau] = buf[i][tau] = 0 ;
    }
    if (compute_pair_corr)
      for (i = 0 ; i < nsites - 1 ; i++)
      for (j = i + 1 ; j < nsites ; j++) avgpair[i][j] = 0 ;
    for (i = 0 ; i < nsitesP + 1 ; i++) binP[i] = 0 ;
    for (i = 0 ; i < nsitesR + 1 ; i++) binR[i] = 0 ;
    if (compute_set)
      for (j = 0 ; j < nmicro ; j++) avgmicro[j] = 0 ;
  }

  nR = nP = 0 ;
  energy = 0.0 ;
  for (i = 0 ; i < nsites ; i++)
  {
    avg[i] += occ[i] = occi = stateocc[i][state[i]] ;
    subavg[i][state[i]]++ ;
    if (site_type[i] == 'P') nP += occi ;
    else nR += occi ;
  }
  /* The auto-correlation function is computed using the on-the-run
     method described by Allen & Tildesley (last paragraph of
     sec. 6.3.1) */
  if (compute_errors)
  {
    for (i = 0 ; i < nsites ; i++)
    {
      tt = t % taumax ;
      buf[i][tt] = occi = occ[i] ;
      for (tau0 = 0 ; tau0 < taumax ; tau0++)
      {
        tau = t - tau0 ;
        if (tau < 0) continue ;
        tau = tau % taumax ;
        cf[i][tau] += buf[i][tau0] * occi ;
        cc[i][tau] += buf[i][tau0] + occi ;
      }
    }
  }
  if (compute_pair_corr)
  {
    for (i = 0 ; i < nsites - 1 ; i++)
    for (j = i + 1 ; j < nsites ; j++) avgpair[i][j] += occ[i] * occ[j] ;
  }
  if (compute_energetics)
  {
    for (i = 0 ; i < nsites ; i++)
    {
      energy += u[i][state[i]] ;
      for (j = i + 1 ; j < nsites ; j++)
        energy += gg[i][state[i]][j][state[j]] ;
    }
  }
  avgP += nP ;
  avgR += nR ;
  avgP2 += nP * nP ;
  avgR2 += nR * nR ;
  avgPR += nP * nR ;
  binP[nP]++ ;
  binR[nR]++ ;

  if (compute_energetics)
  {
    /* offset energy to minimize roundoff errors */
    if (t == 0) energy0 = energy ;
    energy -= energy0 ;
    avgU += energy ;
    avgU2 += energy * energy ;
    for (i = 0 ; i < nsites ; i++) avgUn[i] += occ[i] * energy ;
  }

  if (compute_set)
  {
    j = 0 ;
    for (i = 0 ; i < nset ; i++)
      if (occ[setsite[i]]) j += 1 << i ;
    avgmicro[j]++ ;
  }

}


void write_EpH_point(float E, float pH)
{
  int i, j, s, si, tmp_stateocc ;
  float p, meanP, meanR, stdevP, stdevR, meanU,
        stdevU, dGi, dHi, TdSi, corr, pLi, meani, meanj ;
  int tau, tcor ;
  float corf, corf0 = 0, err ;

  printf("#### E=%f   pH=%f\n", E, pH) ;

  for (i = 0 ; i < nsites ; i++)
  {
    meani = avg[i] / (float) mcsteps ;
    printf(". %8.3f  %7.3f  %4d  %14.8e  ", E, pH, i, meani) ;
    for (s = 0 ; s < nstates[i] ; s++)
      printf(" %14.8e", subavg[i][s] / (float) mcsteps) ;
    printf("\n") ;
    if (site_type[i] == 'P' && pH != pHmin)
    {
      p = pmean[i] ;

/*  Added by VT 2013.04.17 to account for sites that may loose/gain
 * more than 1 proton/electron. As is, it works for a site that captures
 * 2 protons. */
      tmp_stateocc = 1 ;
      for (si = 0 ; si < nstates[i] ; si++) {

	if (stateocc[i][si] > tmp_stateocc) 
	  tmp_stateocc = stateocc[i][si] ;
      }
      switch (tmp_stateocc) {

      case 1 :
	if (((p > 0.5) && (meani <= 0.5)) || ((p < 0.5) && (meani >= 0.5))) {

	  if (npkhalfs[i] == MAXNPKHALFS)
	    error('W', "Number of pKhalf values exceeds %d for site %s.\n",
		  MAXNPKHALFS, site_name[i]) ;
	  else
	    pkhalf[i][npkhalfs[i]++] = pH - dpH * (meani - 0.5) / (meani - p) ;
	}
	break ;
      case 2 :
	if (((p > 0.5) && (meani <= 0.5)) || ((p < 0.5) && (meani >= 0.5))) {

	  if (npkhalfs[i] == MAXNPKHALFS)
	    error('W', "Number of pKhalf values exceeds %d for site %s.\n",
		  MAXNPKHALFS, site_name[i]) ;
	  else
	    pkhalf[i][npkhalfs[i]++] = pH - dpH * (meani - 0.5) / (meani - p) ;
	}
	if (((p > 1.5) && (meani <= 1.5)) || ((p < 1.5) && (meani >= 1.5))) {

	  if (npkhalfs[i] == MAXNPKHALFS)
	    error('W', "Number of pKhalf values exceeds %d for site %s.\n",
		  MAXNPKHALFS, site_name[i]) ;
	  else
	    pkhalf[i][npkhalfs[i]++] = pH - dpH * (meani - 1.5) / (meani - p) ;
	}
	break ;
      default :
	error('E', "Site %s may have unsupported occupancy %d.\n",
	      site_name[i],tmp_stateocc) ;
      }
/* End of change added by VT 2013.04.17 */

    }
    pmean[i] = meani ;
  }
  meanP = avgP / (float) mcsteps ;
  meanR = avgR / (float) mcsteps ;
  stdevP = sqrtp(avgP2 / (float) mcsteps - meanP * meanP) ;
  stdevR = sqrtp(avgR2 / (float) mcsteps - meanR * meanR) ;
  printf(". %8.3f  %7.3f  totP  %14.8e %14.8e\n", E, pH, meanP, stdevP) ;
  printf(". %8.3f  %7.3f  totR  %14.8e %14.8e\n", E, pH, meanR, stdevR) ;

  printf("P %8.3f  %7.3f ", E, pH) ;
  for (i = 0 ; i < nsitesP + 1 ; i++) printf(" %d", binP[i]) ;
  printf("\n") ;
  printf("R %8.3f  %7.3f ", E, pH) ;
  for (i = 0 ; i < nsitesR + 1 ; i++) printf(" %d", binR[i]) ;
  printf("\n") ;

  if (compute_errors)
  {
    for (i = 0 ; i < nsites ; i++)
    {
      meani = avg[i] / (float) mcsteps ;
      tcor = taumax ;
      for (tau = 0 ; tau < taumax ; tau++)
      {
        corf = (cf[i][tau] - meani * cc[i][tau]) / (float) (mcsteps - tau) +
               meani * meani ;
        if (tau == 0) corf0 = corf ;
        if (corf0 == 0) corf = 1 ;
        else corf = corf / corf0 ;
        /* tcor is estimated as the time at which corf < 0.1, as done
           by Beroza et al. */
        if ((tcor == taumax) && (fabs(corf) < 0.10)) tcor = tau ;
      }
      /* Eq. 6 of Beroza et al. See also Eq. 6.13 of Allen & Tildesley */
      err = sqrt(corf0 * tcor / (float) mcsteps) ;
      printf("t %8.3f  %7.3f  %4d  %14.8e (%d)\n", E, pH, i, err, tcor) ;
    }
  }

  if (compute_energetics)
  {
    meanU = avgU / (float) mcsteps ;
    for (i = 0 ; i < nsites ; i++)
    {
      meani = avg[i] / (float) mcsteps ;
      /* kT units */
      if (meani != 0 && meani != 1)
      {
        if (site_type[i] == 'P') pLi = pH ;
        else pLi = - E / (LN10 * kBoltz_meV * T) ;
        dGi = - log(meani / (1.0 - meani)) - LN10 * pLi ;
        TdSi = (avgUn[i] / (float) mcsteps - meani * meanU) / 
               (meani - meani * meani) + log(meani / (1.0 - meani)) ;
        dHi = dGi + TdSi ;
        printf("e %8.3f  %7.3f  %4d  %15.8e %15.8e %15.8e\n",
               E, pH, i, dGi, dHi, TdSi) ;
      }
    }
    stdevU = sqrtp(avgU2 / (float) mcsteps - meanU * meanU) ;
    /* energy offset has to be added to meanU but not to stdevU */
    printf("e %8.3f  %7.3f   tot  %15.8e %14.8e\n",
           E, pH, meanU + energy0, stdevU) ;
  }

  if (compute_pair_corr)
  {
    for (i = 0 ; i < nsites - 1 ; i++)
    {
      meani = avg[i] / (float) mcsteps ;
      if (meani == 0 || meani == 1) continue ;
      for (j = i + 1 ; j < nsites ; j++)
      {
        meanj = avg[j] / (float) mcsteps ;
        if (meanj == 0 || meanj == 1) continue ;
        corr = (avgpair[i][j] / (float) mcsteps - meani * meanj) /
               sqrtp((meani - meani * meani) * (meanj - meanj * meanj)) ;
        if (fabs(corr) >= min_corr)
          printf(": %8.3f  %7.3f  %4d %4d  %15.8e\n", E, pH, i, j, corr) ;
      }
    }
  }
  if (stdevP == 0 || stdevR == 0) corr = 0 ;
  else corr = (avgPR / (float) mcsteps - meanP * meanR) / (stdevP * stdevR) ;
  printf(": %8.3f  %7.3f  totP totR  %15.8e\n", E, pH, corr) ;

  if (compute_set)
  {
    for (j = 0 ; j < nmicro ; j++)
    {
      printf("m %8.3f  %7.3f  ", E, pH) ;
      for (i = 0 ; i < nset ; i++) printf("%1d", 1 & (j >> i)) ;
      printf("  %14.8e\n", avgmicro[j] / (float) mcsteps) ;
    }
  }

}


void write_E_line(float E)
{
  int i, k ;

  printf("#\n# pKhalfs for E = %f :\n", E) ;
  printf("#   site       number     E      pKhalf(s)\n") ;
  printf("#-----------------------------------------------\n") ;
  for (i = 0 ; i < nsites ; i++)
  {
    if (titration[i] == 4) continue ;
    printf("> %-13s %4d %8.3f  ", site_name[i], i, E) ;
    if (titration[i] == 0) printf("NonTitrable\n") ;
    else if (site_type[i] == 'R') printf("    ----\n") ;
    else if (npkhalfs[i] == 0) printf("NotInRange\n") ;
    else
    {
      printf("%9.4f", pkhalf[i][0]) ;
      for (k = 1 ; k < npkhalfs[i] ; k++) printf(", %.4f", pkhalf[i][k]) ;
      printf("\n") ;
    }
  }
  printf("#\n") ;
}


/* This treats negative numbers as zero. It is used because of rounding
   errors in numbers that should be non-negative by definition. */
double sqrtp(double x)
{
  return (x>0 ? sqrt(x) : 0) ;
}

/*
void clean_memory(void)
{
  int i, j, s ;

  free(site_type) ;
  for (i = 0 ; i < nsites ; i++) free(site_name[i]) ;
  free(site_name) ;
  free(titration) ;
  free(state) ;
  free(occ) ;
  free(pair1) ;
  free(pair2) ;
  free(avg) ;
  free(pmean) ;
  free(binP) ;
  free(binR) ;
  for (i = 0 ; i < nsites ; i++) free(pkhalf[i]) ;
  free(pkhalf) ;
  free(npkhalfs) ;
  for (i = 0 ; i < nsites ; i++)
  {
    for (s = 0 ; s < nstates[i] ; s++)
    {
      for (j = 0 ; j < nsites ; j++) free(gg[i][s][j]) ;
      free(gg[i][s]) ;
    }
    free(gg[i]) ;
  }
  free(gg) ;
  free(nstates) ;
  for (i = 0 ; i < nsites ; i++) free(g[i]) ;
  free(g) ;
  for (i = 0 ; i < nsites ; i++) free(stateocc[i]) ;
  free(stateocc) ;
  for (i = 0 ; i < nsites ; i++) free(u[i]) ;
  free(u) ;
  for (i = 0 ; i < nsites ; i++) free(subavg[i]) ;
  free(subavg) ;
  if (compute_pair_corr)
  {
    for (i = 0 ; i < nsites ; i++) free(avgpair[i]) ;
    free(avgpair) ;
  }
  if (compute_energetics) free(avgUn) ;
  if (compute_set)
  {
    free(setsite) ;
    free(avgmicro) ;
  }
  if (compute_errors)
  {
    for (i = 0 ; i < nsites ; i++) free(cf[i]) ;
    free(cf) ;
    for (i = 0 ; i < nsites ; i++) free(cc[i]) ;
    free(cc) ;
    for (i = 0 ; i < nsites ; i++) free(buf[i]) ;
    free(buf) ;
  }
  strsplit(NULL, ",") ;
  }
*/



/* This splits string s, using any character in delim, assigning the
   resulting parts to substr and the number of parts to nsubstr. Both
   substr and nsubstr must be shared by this function and the calling
   one. When called with s == NULL it frees the memory of substr. It
   does not work very well... */
void strsplit(char *s, const char *delim)
{
  int i ;
  char *c, aux[LINESIZE] ;

  if (substr != NULL)
  {
    for (i = 0 ; i < nsubstr ; i++) free(substr[i]) ;
    free(substr) ;
    if (s == NULL) return ;
  }
  strcpy(aux, s) ;
  c = strtok(aux, delim) ;
  for (nsubstr = 0 ; c != NULL ; nsubstr++) c = strtok(NULL, delim) ;
  if (nsubstr > 0)
  {
    substr = calloc(nsubstr, sizeof(char *)) ;
    for (i = 0 ; i < nsubstr ; i++) substr[i] = calloc(STRSIZE, sizeof(char)) ;
    strcpy(aux, s) ;
    strcpy(substr[0], strtok(aux, delim)) ;
    for (i = 1 ; i < nsubstr ; i++) strcpy(substr[i], strtok(NULL, delim)) ;
  }
}


void error(char errtype, char *format, ...)
{
  va_list args;

  va_start(args, format);
  if (errtype != 'W' && errtype != 'E' && errtype != 'U')
    error('E', "Wrong use of error function.\n") ;
  if (errtype == 'W') fprintf(stderr, "%s: WARNING: ", cmd);
  else fprintf(stderr, "%s: ERROR: ", cmd);
  vfprintf(stderr, format, args);
  va_end(args);
  if (errtype == 'U') usage() ;
  if (errtype != 'W') exit(1) ;
  return ;
}

