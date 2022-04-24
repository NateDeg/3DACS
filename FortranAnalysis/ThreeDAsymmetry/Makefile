mFlags=../makeflags  
include $(mFlags)

all: $(ThreeDAsymObj) install


AsymTypes.o :
	$(MakeFObj) AsymTypes.f

ThreeDAsym_Core.o : AsymTypes.o
	$(MakeFObj) ThreeDAsym_Core.f

ThreeDAsym_Minimum.o : ThreeDAsym_Core.o
	$(MakeFObj) ThreeDAsym_Minimum.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
