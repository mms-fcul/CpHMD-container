#!/usr/bin/awk -f

#############################################################################
#                                                                           #
# This reads a file with coordinates from a gromacs run and adds a slice    #
# of atoms which corresponds to a percentage of the box vector. It may be   #
# in X and Y directions or all 3 (X, Y and Z). The coordinate file is pdb   #
# like. It can be a pqr file with charge and radius from which only the     #
# radius is read (charge is placed 0.000 for new added atoms), or a pdb     #
# file without the occupation and b-factor columns as the input generated   #
# for delphi. The vector is the last line from a gro gromacs file (in A).   #
# Slice is a percentage of the vector which is going to be built by         #
# periodicity in X and Y, or in X, Y and Z.                                 #
#                                                                           #
#                                                                           #
# Vitor H Teixeira, Univ. Lisboa, 2013.03.26                                #
#                                                                           #
#############################################################################

BEGIN{
  cmd = "slice.awk" ;
  usage = "Usage: "cmd" <coords> <vector> <slice> <out> <2/3>\n" \
      " * coords is a pqr file\n" \
      " * vector is the vector from the gromacs gro file\n" \
      " * slice is the percent of vector to add in x,-x,y and -y\n" \
      " * 2/3 For PBC in 2D or 3D\n." ;

  if (ARGC != 6) 
    error("Wrong number of arguments.\n" usage) ;

  filecheck(coords = ARGV[1]) ;
  vector = ARGV[2] ;
  slice = ARGV[3] ;
  out = ARGV[4] ;
  dimension = ARGV[5] ;

  if ((dimension != 2) && (dimension != 3))
    error("Bad option for dimension\n" usage) ;

####  Setting higher and lower limit according to slice
  split (vector,vec,",") ;
  lowerX = slice * vec[1] ;
  lowerY = slice * vec[2] ;
  lowerZ = slice * vec[3] ;
  higherX = (1 - slice) * vec[1] ;
  higherY = (1 - slice) * vec[2] ;
  higherZ = (1 - slice) * vec[3] ;

####  Numner of atoms
  line = last_atom = 0 ;


###  Reading coords file and processing
  while (getline < coords) {

    if ($1 == "ATOM") {

      line++ ;
      residue = $5 + 1 ;
    }
  }
  close (coords) ;

### Assign total number of lines to atoms
  atoms = line ;

  while (getline < coords) {

    if($1 == "ATOM") {
      if (dimension == 3) {
        if ($6 <= lowerX) {
	  if ($7 <= lowerY) {
	    if ($8 <= lowerZ) {

	      add_atom_3("plus","plus","plus") ;
	      add_atom_2("x","y","plus","plus") ;
	      add_atom_2("x","z","plus","plus") ;
	      add_atom_2("y","z","plus","plus") ;
	      add_atom_1("x","plus") ;
	      add_atom_1("y","plus") ;
	      add_atom_1("z","plus") ;
	    } else if ($8 >= higherZ) {

	      add_atom_3("plus","plus","minus") ;
	      add_atom_2("x","y","plus","plus") ;
	      add_atom_2("x","z","plus","minus") ;
	      add_atom_2("y","z","plus","minus") ;
	      add_atom_1("x","plus") ;
	      add_atom_1("y","plus") ;
	      add_atom_1("z","minus") ;
	    }
	  } else if ($7 >= higherY) {
	    if ($8 <= lowerZ) {

	      add_atom_3("plus","minus","plus") ;
	      add_atom_2("x","y","plus","minus") ;
	      add_atom_2("x","z","plus","plus") ;
	      add_atom_2("y","z","minus","plus") ;
	      add_atom_1("x","plus") ;
	      add_atom_1("y","minus") ;
	      add_atom_1("z","plus") ;
	    } else if ($8 >= higherZ) {

	      add_atom_3("plus","minus","minus") ;
	      add_atom_2("x","y","plus","minus") ;
	      add_atom_2("x","z","plus","minus") ;
	      add_atom_2("y","z","minus","minus") ;
	      add_atom_1("x","plus") ;
	      add_atom_1("y","minus") ;
	      add_atom_1("z","minus") ;
	    }
	  }
	} else if ($6 >= higherX) {
	  if ($7 <= lowerY) {
	    if ($8 <= lowerZ) {

	      add_atom_3("minus","plus","plus") ;
	      add_atom_2("x","y","minus","plus") ;
	      add_atom_2("x","z","minus","plus") ;
	      add_atom_2("y","z","plus","plus") ;
	      add_atom_1("x","minus") ;
	      add_atom_1("y","plus") ;
	      add_atom_1("z","plus") ;
	    } else if ($8 >= higherZ) {

	      add_atom_3("minus","plus","minus") ;
	      add_atom_2("x","y","minus","plus") ;
	      add_atom_2("x","z","minus","minus") ;
	      add_atom_2("y","z","plus","minus") ;
	      add_atom_1("x","minus") ;
	      add_atom_1("y","plus") ;
	      add_atom_1("z","minus") ;
	    }
	  } else if ($7 >= higherY) {
	    if ($8 <= lowerZ) {

	      add_atom_3("minus","minus","plus") ;
	      add_atom_2("x","y","minus","minus") ;
	      add_atom_2("x","z","minus","plus") ;
	      add_atom_2("y","z","minus","plus") ;
	      add_atom_1("x","minus") ;
	      add_atom_1("y","minus") ;
	      add_atom_1("z","plus") ;
	    } else if ($8 >= higherZ) {

	      add_atom_3("minus","minus","minus") ;
	      add_atom_2("x","y","minus","minus") ;
	      add_atom_2("x","z","minus","minus") ;
	      add_atom_2("y","z","minus","minus") ;
	      add_atom_1("x","minus") ;
	      add_atom_1("y","minus") ;
	      add_atom_1("z","minus") ;
	    }
	  }
	}
	if (($6 <= lowerX && $7 <= lowerY) && ($8 > lowerZ && $8 < higherZ)) {

	  add_atom_2("x","y","plus","plus") ;
	  add_atom_1("x","plus") ;
	  add_atom_1("y","plus") ;
	} else if (($6 <= lowerX && $7 >= higherY) && ($8 > lowerZ && $8 < higherZ)) {

	  add_atom_2("x","y","plus","minus") ;
	  add_atom_1("x","plus") ;
	  add_atom_1("y","minus") ;
	} else if (($6 >= higherX && $7 <= lowerY) && ($8 > lowerZ && $8 < higherZ)) {

	  add_atom_2("x","y","minus","plus") ;
	  add_atom_1("x","minus") ;
	  add_atom_1("y","plus") ;
	} else if (($6 >= higherX && $7 >= higherY) && ($8 > lowerZ && $8 < higherZ)) {

	  add_atom_2("x","y","minus","minus") ;
	  add_atom_1("x","minus") ;
	  add_atom_1("y","minus") ;
	} else if (($6 <= lowerX && $8 <= lowerZ) && ($7 > lowerY && $7 < higherY)) {

	  add_atom_2("x","z","plus","plus") ;
	  add_atom_1("x","plus") ;
	  add_atom_1("z","plus") ;
	} else if (($6 <= lowerX && $8 >= higherZ) && ($7 > lowerY && $7 < higherY)) {

	  add_atom_2("x","z","plus","minus") ;
	  add_atom_1("x","plus") ;
	  add_atom_1("z","minus") ;
	} else if (($6 >= higherX && $8 <= lowerZ) && ($7 > lowerY && $7 < higherY)) {

	  add_atom_2("x","z","minus","plus") ;
	  add_atom_1("x","minus") ;
	  add_atom_1("z","plus") ;
	} else if (($6 >= higherX && $8 >= higherZ) && ($7 > lowerY && $7 < higherY)) {

	  add_atom_2("x","z","minus","minus") ;
	  add_atom_1("x","minus") ;
	  add_atom_1("z","minus") ;
	} else if (($7 <= lowerY && $8 <= lowerZ) && ($6 > lowerX && $6 < higherX)) {

	  add_atom_2("y","z","plus","plus") ;
	  add_atom_1("y","plus") ;
	  add_atom_1("z","plus") ;
	} else if (($7 <= lowerY && $8 >= higherZ) && ($6 > lowerX && $6 < higherX)) {

	  add_atom_2("y","z","plus","minus") ;
	  add_atom_1("y","plus") ;
	  add_atom_1("z","minus") ;
	} else if (($7 >= higherY && $8 <= lowerZ) && ($6 > lowerX && $6 < higherX)) {

	  add_atom_2("y","z","minus","plus") ;
	  add_atom_1("y","minus") ;
	  add_atom_1("z","plus") ;
	} else if (($7 >= higherY && $8 >= higherZ) && ($6 > lowerX && $6 < higherX)) {

	  add_atom_2("y","z","minus","minus") ;
	  add_atom_1("y","minus") ;
	  add_atom_1("z","minus") ;
	} else if ($6 <= lowerX) {
	  if (($7 > lowerY && $7 < higherY) && ($8 > lowerZ && $8 < higherZ)) 

	    add_atom_1("x","plus") ;
	} else if ($6 >= higherX) {
	  if (($7 > lowerY && $7 < higherY) && ($8 > lowerZ && $8 < higherZ)) 

	    add_atom_1("x","minus") ;
	} else if ($7 <= lowerY) {
	  if (($6 > lowerX && $6 < higherX) && ($8 > lowerZ && $8 < higherZ)) 

	    add_atom_1("y","plus") ;
	} else if ($7 >= higherY) {
	  if (($6 > lowerX && $6 < higherX) && ($8 > lowerZ && $8 < higherZ)) 

	    add_atom_1("y","minus") ;
	} else if ($8 <= lowerZ) {
	  if (($6 > lowerX && $6 < higherX) && ($7 > lowerY && $7 < higherY)) 

	    add_atom_1("z","plus") ;
	} else if ($8 >= higherZ) {
	  if (($6 > lowerX && $6 < higherX) && ($7 > lowerY && $7 < higherY)) 

	    add_atom_1("z","minus") ;
	}
      } else if (dimension == 2) {

	if ($6 <= lowerX) {
	  if ($7 <= lowerY) {

	    add_atom_2("x","y","plus","plus") ;
	    add_atom_1("x","plus") ;
	    add_atom_1("y","plus") ;
	  } else if ($7 >= higherY) {

	    add_atom_2("x","y","plus","minus") ;
	    add_atom_1("x","plus") ;
	    add_atom_1("y","minus") ;
	  } else {

	    add_atom_1("x","plus") ; 
	  }
	} else if ($6 >= higherX) {
	  if ($7 <= lowerY) {

	    add_atom_2("x","y","minus","plus") ;
	    add_atom_1("x","minus") ;
	    add_atom_1("y","plus") ;
	  } else if ($7 >= higherY) {

	    add_atom_2("x","y","minus","minus") ;
	    add_atom_1("x","minus") ;
	    add_atom_1("y","minus") ;
	  } else {

	    add_atom_1("x","minus") ;
	  }
	} else if ($7 <= lowerY) {

	  add_atom_1("y","plus") ; 
	} else if ($7 >= higherY) {

	  add_atom_1("y","minus") ;
	}
      }
      print $0 >out ;
###  Save number of atoms
      last_atom++ ;
    }
  }
  close(coords) ;

####  Writing new coords
  for (i = (last_atom + 1); i <= atoms; i++)

    printf("ATOM %6d %-4s %-4s %4d %11.3f%8.3f%8.3f %6.3f %6.3f\n",
	   i,atom_n[i],res_n[i],res_i[i],x[i],y[i],z[i],
	   0.000,radii[i]) >out ;
}


function warning(msg) {

  print cmd ": Warning: " msg | "cat 1>&2" ;
}


function error(msg) {

  print cmd ": Error: " msg | "cat 1>&2" ;
  exit 1 ;
}


function filecheck(file) {

  if (system("test -f "file))
    error("File "file" does not exist.") ;
  if (system("test -r "file))
    error("File "file" exists but is not readable.") ;
}


function add_atom_3(a,b,c) {

  atoms++ ;
  if ((atoms % 100) == 0) residue++ ;
  atom_n[atoms] = $3 ;
  res_i[atoms] = residue ;
  res_n[atoms] = $4 ;
  z[atoms] = $8 ;
  radii[atoms] = $10 ;

  if (a == "plus" && b == "plus" && c == "plus") {

    x[atoms] = $6 + vec[1] ;
    y[atoms] = $7 + vec[2] ;
    z[atoms] = $8 + vec[3] ;
  } else if (a == "plus" && b == "plus" && c == "minus") {

    x[atoms] = $6 + vec[1] ;
    y[atoms] = $7 + vec[2] ;
    z[atoms] = $8 - vec[3] ;
  } else if (a == "plus" && b == "minus" && c == "plus") {

    x[atoms] = $6 + vec[1] ;
    y[atoms] = $7 - vec[2] ;
    z[atoms] = $8 + vec[3] ;
  } else if (a == "plus" && b == "minus" && c == "minus") {

    x[atoms] = $6 + vec[1] ;
    y[atoms] = $7 - vec[2] ;
    z[atoms] = $8 - vec[3] ;
  } else if (a == "minus" && b == "plus" && c == "plus") {

    x[atoms] = $6 - vec[1] ;
    y[atoms] = $7 + vec[2] ;
    z[atoms] = $8 + vec[3] ;
  } else if (a == "minus" && b == "plus" && c == "minus") {

    x[atoms] = $6 - vec[1] ;
    y[atoms] = $7 + vec[2] ;
    z[atoms] = $8 - vec[3] ;
  } else if (a == "minus" && b == "minus" && c == "plus") {

    x[atoms] = $6 - vec[1] ;
    y[atoms] = $7 - vec[2] ;
    z[atoms] = $8 + vec[3] ;
  } else if (a == "minus" && b == "minus" && c == "minus") {

    x[atoms] = $6 - vec[1] ;
    y[atoms] = $7 - vec[2] ;
    z[atoms] = $8 - vec[3] ;
  }
}


function add_atom_2(a,b,c,d) {

  atoms++ ;
  if ((atoms % 100) == 0) residue++ ;
  atom_n[atoms] = $3 ;
  res_i[atoms] = residue ;
  res_n[atoms] = $4 ;
  z[atoms] = $8 ;
  radii[atoms] = $10 ;

  if (a == "x" && b == "y") {
    if (c == "plus" && d == "plus") {

      x[atoms] = $6 + vec[1] ;
      y[atoms] = $7 + vec[2] ;
      z[atoms] = $8 ;
    } else if (c == "plus" && d == "minus") {

      x[atoms] = $6 + vec[1] ;
      y[atoms] = $7 - vec[2] ;
      z[atoms] = $8 ;
    } else if (c == "minus" && d == "plus") {

      x[atoms] = $6 - vec[1] ;
      y[atoms] = $7 + vec[2] ;
      z[atoms] = $8 ;
    } else if (c == "minus" && d == "minus") {

      x[atoms] = $6 - vec[1] ;
      y[atoms] = $7 - vec[2] ;
      z[atoms] = $8 ;
    }
  } else if (a == "x" && b == "z") {
    if (c == "plus" && d == "plus") {

      x[atoms] = $6 + vec[1] ;
      y[atoms] = $7 ;
      z[atoms] = $8 + vec[3] ;
    } else if (c == "plus" && d == "minus") {

      x[atoms] = $6 + vec[1] ;
      y[atoms] = $7 ;
      z[atoms] = $8 - vec[3] ;
    } else if (c == "minus" && d == "plus") {

      x[atoms] = $6 - vec[1] ;
      y[atoms] = $7 ;
      z[atoms] = $8 + vec[3] ;
    } else if (c == "minus" && d == "minus") {

      x[atoms] = $6 - vec[1] ;
      y[atoms] = $7 ;
      z[atoms] = $8 - vec[3] ;
    }
  } else if (a == "y" && b == "z") {
    if (c == "plus" && d == "plus") {

      x[atoms] = $6 ;
      y[atoms] = $7 + vec[2] ;
      z[atoms] = $8 + vec[3] ;
    } else if (c == "plus" && d == "minus") {

      x[atoms] = $6 ;
      y[atoms] = $7 + vec[2] ;
      z[atoms] = $8 - vec[3] ;
    } else if (c == "minus" && d == "plus") {

      x[atoms] = $6 ;
      y[atoms] = $7 - vec[2] ;
      z[atoms] = $8 + vec[3] ;
    } else if (c == "minus" && d == "minus") {

      x[atoms] = $6 ;
      y[atoms] = $7 - vec[2] ;
      z[atoms] = $8 - vec[3] ;
    }
  }
}


function add_atom_1(a,sign) {

  atoms++ ;
  if ((atoms % 100) == 0) residue++ ;
  atom_n[atoms] = $3 ;
  res_i[atoms] = residue ;
  res_n[atoms] = $4 ;
  z[atoms] = $8 ;
  radii[atoms] = $10 ;

  if (a == "x") {
    if (sign == "plus") {

      x[atoms] = $6 + vec[1] ;
      y[atoms] = $7 ;
      z[atoms] = $8 ;
    } else if (sign == "minus") {

      x[atoms] = $6 - vec[1] ;
      y[atoms] = $7 ;
      z[atoms] = $8 ;
    }
  } else if (a == "y") {
    if (sign == "plus") {

      x[atoms] = $6 ;
      y[atoms] = $7 + vec[2] ;
      z[atoms] = $8 ;
    } else if (sign == "minus") {

      x[atoms] = $6 ;
      y[atoms] = $7 - vec[2] ;
      z[atoms] = $8 ;
    }
  } else if (a == "z") {
    if (sign == "plus") {

      x[atoms] = $6 ;
      y[atoms] = $7 ;
      z[atoms] = $8 + vec[3] ;
    } else if (sign == "minus") {

      x[atoms] = $6 ;
      y[atoms] = $7 ;
      z[atoms] = $8 - vec[3] ;
    }
  }
}

