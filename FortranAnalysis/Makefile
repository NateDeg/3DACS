mFlags=makeflags 
include $(mFlags)


all: LocalLib Programs

Full:
	cd ../third_party/cfitsio; make clean; ./configure; make
	make clean
	make


LocalLib: 
	cd StandardMath; make;
	cd ObjectDefinitions; make;
	cd Globals; make;
	cd SNEstimate; make;
	cd Inputs; make all;
	cd Initialization; make all;
	cd Minimizers; make;
	cd ThreeDAsymmetry; make; 
	cd Outputs; make;
	cd Mains; make;


Programs:
	cd $(LocBinLoc); make; make install



clean:
	cd StandardMath; make clean;
	cd ObjectDefinitions; make clean;
	cd Globals; make clean;
	cd SNEstimate; make clean;
	cd Inputs; make clean;
	cd Minimizers; make clean;
	cd Initialization; make clean;
	cd ThreeDAsymmetry; make clean;
	cd Outputs; make clean;
	cd Mains; make clean;


	cd $(LocBinLoc); rm *.o;
	cd $(LocModLoc); rm *.mod
	cd ../Programs; rm $(AllProgs)

	

