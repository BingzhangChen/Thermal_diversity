.SUFFIXES: .o .i .f90 .f
.f90.o:
	$(CFT) -c $(FFLAGS) $*.f90 -o $*.o
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o
CFT  = ifort
#CFT  = gfortran
FFLAGS=-O3 -xHost -r8 -i4
#FFLAGS=-O0 -r8 -i4 -g -check all -fpe-all=0 -warn -traceback -debug extended -nogen-interface -warn interfaces
#FFLAGS=-O3 -fdefault-real-8 -fdefault-double-8 -march=native -fimplicit-none -Wall  -Wline-truncation  -fwhole-file  -std=f2008
#FFLAGS=-O0 -g -fdefault-real-8 -fdefault-double-8 -fimplicit-none -Wall -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace
MAIN = Topt
#Source file
SRCS = slab.f90 siteTopt.f90

OBJS =$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
#Executable file
$(MAIN): $(OBJS)
	$(CFT) $(FFLAGS) $(LFLAGS) -o $(MAIN) $(OBJS)

clean:
	rm -rf *.o *.mod
