#!/usr/bin/awk -f

#########################################################################
# Averaging_g: Makes an average of the calculated pairwise interactions.
#
#   The program reads a file with calculated pairwise interactions and 
# makes the average of them. This is important because the pairwise
# interactions are calculated in both sites as a-b and b-a and are not 
# equal.
#########################################################################


BEGIN{
  cmd = "Filer.awk" ;
  usage = "Usage: "cmd" <Concatenated_g> <out_g>\n" ;
  if (ARGC != 3) error("Wrong number of arguments.\n" usage) ;

  filecheck(input = ARGV[1]) ;
  output = ARGV[2] ;



####  Reading calculated pairwise interactions
  while (getline < input) {

    n++ ;
    a[n] = $1 ;
    b[n] = $2 ;
    inter[$1,$2] = $3 ;
  }
  close (input) ;

####  Symetrizing interactions
  for (i = 1; i <= n; i++) {
      
    junk = (inter[a[i],b[i]] + inter[b[i],a[i]]) / 2 ;
    printf("%d %d %15.6e\n",a[i],b[i],junk) >output ;
  }
}


function warning(msg) {

  print cmd ": Warning: " msg | "cat 1>&2" ;
  close ("cat 1>&2") ;
}

function error(msg) {

  print cmd ": Error: " msg | "cat 1>&2" ;
  close ("cat 1>&2") ;
  exit 1 ;
}

function filecheck(file) {

  if (system("test -f "file))
    error("File "file" does not exist.") ;
  if (system("test -r "file))
    error("File "file" exists but is not readable.") ;
}

