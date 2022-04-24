mFlags=../makeflags 
include $(mFlags)

VPATH :=$(BinLoc)

ProgramList= $(AllProgs)

all: $(ProgramList)


Asymmetry3D : Asymmetry3D.o
	$(F77) $(FLAGS) -o Asymmetry3D  Asymmetry3D.o $(AllObj) $(FitsLibLoc) 

Min_Asymmetry3D : Min_Asymmetry3D.o
	$(F77) $(FLAGS) -o Min_Asymmetry3D  Min_Asymmetry3D.o $(AllObj) $(FitsLibLoc)

install:
	mv $(AllProgs) $(PARENT_DIR)/Programs


clean:
	rm -f $(ProgramList);
	rm -f *.o;


