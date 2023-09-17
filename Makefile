#
#  Makefile to create IntTest
#
COMPILER = gfortran #ifort
CFLAGS   = -mcmodel=medium -ffree-line-length-512 -fdec #-double-size 64 -integer-size 64 -shared-intel
#CFLAGS   =  -check all
#CFLAGS   =  -no-ipo -r8 -i8
#CFLAGS    =  -mcmodel=medium -r8 -i8 -double-size 64 -integer-size 64 -traceback 
LIBS     = 

#  Objects to create
#
OBJECTS = nrtype.o \
	  shared.o \
	  funcs.o \
          FindEigenvalue.o

#	  DCUHREm.o \
#	  quadpack.o \
#	  quadpack2.o \

#	  nrutil.o \
#	  myIntegrate.o \

#  Pattern rule(s)  
#
%.o : %.f90
	$(COMPILER) $(CFLAGS) -c $<

#  Ultimate target
#
Dirac: $(OBJECTS)
	$(COMPILER) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS)

#  Dependencies
#
#myKinds.o: 

#newtonRaph.o:

#rungeKutta.o: myKinds.o

#odeInt.o: rungeKutta.o

#SolveSE.o: odeInt.o newtonRaph.o 


# Clean up rule
#
clean: 
	rm -f *.o *.mod *~ Dirac
