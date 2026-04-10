VPATH=./src
#LINUX INTEL
  FC=pgf95
  CC=pgcc

  FLAGS= -O3 -I$(VPATH)/ -DLINUX -DARM64 -Mcray=pointer -Mvect -Mcache_align -Munroll -Mnoframe 
  LFLAGS = -O3 -L/usr/lib64/ 

#  FC=gfortran
#  FLAGS=-g -fcray-pointer -ffixed-line-length-none -O3 -DLINUX -DIFC -DDP -fdefault-real-8


include $(VPATH)/genericLinux.mk


