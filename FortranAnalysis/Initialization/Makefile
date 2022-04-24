mFlags=../makeflags  
include $(mFlags)

all: $(IniObj) install

AsymmetryInitializations.o:
	$(MakeFObj) AsymmetryInitializations.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
