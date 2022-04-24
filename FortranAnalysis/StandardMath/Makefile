mFlags=../makeflags  
include $(mFlags)

all: $(StandardMathObj) install


BasicConstants.o:
	$(F77) $(FLAGS) -c -I$(LocModLoc) BasicConstants.f

Interpolation.o:
	$(F77) $(FLAGS) -c -I$(LocModLoc) Interpolation.f

SortArray.o:
	$(F77) $(FLAGS) -c -I$(LocModLoc) SortArray.f

SimpleSmoothing.o: BasicConstants.o
	$(F77) $(FLAGS) -c -I$(LocModLoc) SimpleSmoothing.f


TabulatedIntegral.o: Interpolation.o
	$(F77) $(FLAGS) -c -I$(LocModLoc) TabulatedIntegral.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
