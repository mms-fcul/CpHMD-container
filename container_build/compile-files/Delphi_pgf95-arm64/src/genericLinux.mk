SRC1=qlog.f95 pointers.f95 oper_coord.f95 misc.f95
#
SRC2=misc2.f95
#
SRC3=crgarrmod.f95  epsmakmod.f95  qqintmod.f95  setbcmod.f95   setrcmod.f95 \
encalcmod.f95  extrmmod.f95   rdhmod.f95    setcrgmod.f95  wrtsitmod.f95
#
#----------------------------------------------------
OBJ1=$(SRC1:.f95=.o)
OBJ2=$(SRC2:.f95=.o)
OBJ3=$(SRC3:.f95=.o)
#----------------------------------------------------
default: delphi95
#----------------------------------------------------
clean:
#	echo $(FC) $(FLAGS)
	rm -f *.o *.mod
#----------------------------------------------------

#----------------------------------------------------
delphi95:$(OBJ1) $(OBJ2) $(OBJ3) $(SRC1) $(SRC2) $(SRC3)
#	$(FC) $(FLAGS) -o $@ $(VPATH)/qdiff4v.F95 $(OBJ1) $(OBJ2) $(OBJ3) 
	$(FC) $(FLAGS) -o delphi95 $(VPATH)/qdiff4v.F95 $(OBJ1) $(OBJ2) $(OBJ3) 
#---------------------------------------------------
#compile:
#       $(FC) $(FLAGS) -c  $(SRC1) $(SRC2) $(SRC3)

%.o:%.f95
#	$(FC) $(FLAGS) -c $(VPATH) /$< 
	$(FC) $(FLAGS) -c $(VPATH)/$*.f95 
#
