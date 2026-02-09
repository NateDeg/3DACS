mFlags=../makeflags  
include $(mFlags)

all: $(LopsideObj) install


ProfileLopsidedNessCalc.o :
	$(MakeFObj) ProfileLopsidedNessCalc.f



clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
