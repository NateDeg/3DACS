mFlags=../makeflags  
include $(mFlags)

all: $(InputObj) install

InputUnitConversions.o :
	$(MakeFObj) InputUnitConversions.f

DataCubeInput.o : InputUnitConversions.o
	$(MakeFObj) DataCubeInput.f

AsymmetryIn.o : DataCubeInput.o
	$(MakeFObj) AsymmetryIn.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
