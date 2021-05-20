
OFILES =   ranfun.o vector.o kernel.o

HEAD = example.h  twalk.h
OPTI =  -O -D_FILE_OFFSET_BITS=64 # to use 64 bit file access for very big files
#-g debugger, -O to optimize

#***** If you want to save only the ACCEPTED iterations include
#DEFS = -DSAVEACCONLY #-DFLUSHEVERYIT #and this to flush points at every iteration


### THE GSL LIBRARY IS NO LONGER NEEDED, ALL UTILITIES ARE COPIED FROM THE GSL SOURCE FILES IN ranfun.cpp
LIBS = -lm


## If you want to still use the gsl library change ranfun_gsl.c to ranfun.cpp and use:
##LIBS = -lm -lgsl -lgslcblas 

### The gsl library is included in most modern linux distributions,
### for most Unixes, like Mac and (even!) for windows it may be downloaded from:
### http://www.gnu.org/software/gsl/



COMPILER = g++ 

CFLAGS = $(DEFS) $(OPTI)  -Wno-write-strings


example : $(HEAD) example.o $(OFILES) 
	$(COMPILER)  $(CFLAGS) -o example example.o $(OFILES) $(LIBS)


%.o: %.cpp $(HEAD)
	$(COMPILER) $(CFLAGS) -c -o $@ $< 

