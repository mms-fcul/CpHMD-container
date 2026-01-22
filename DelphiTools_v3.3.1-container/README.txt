#################################################
###### New wrapper for delphi calculations ######
#################################################
Date 17/09/2024

- Reduced the database for CHARMM to test out the database limits

Date 31/05/2024

- Added the previour refered check, when pkint is too large then it outputs an error for wrong calc_pkg_low calculation.

--------------------------------------
Date: 29/05/24

It was found that for large systems the delphi calculation started to break
leading for pkints that were way too large to be ok.
This problem arose on PGP calculations where acidic residues at pH 7.0 and
fully solvated, were yelding a completely protonated state. Super odd and it
turns out something went wrong on the calc_pkg_low.

- Miguel recompiled calc_pkg_low which was then replaced in version Delphi_v3.3

- additionally was suggested to add a error line when values above 100 were
outputed in the pkint (clearly indicating an unhealthy simulation) - This will
be worked out in the meantime

---------------------------------------------

Date: --/--/24

- Change the N and C termini offset from 2 and 3* to +100 and +200 
