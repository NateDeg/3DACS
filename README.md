3D Asymmetries in data CubeS
aka 3DACS

version: 1.0

June 5, 2023

----
Description:

This package contains 2 codes for calculating the 3D asymmetry.  The first is:
Asymmetry3D
which calculates the 1D, 2D, and 3D asymmetry of a cube (plus mask) about a specific point.

The second code is:
Min_Asymmetry3D
which attempts to find the point that minimizes the asymmetry for the 1D, 2D, and 3D calculations.    

----
authors: N. Deg

--
Paper Reference:
Deg et al., 2023, MNRAS


-----

DEPENDENCIES:
gfortran, cfitsio*
*For convenience, versions of the cfitsio libraries are contained in:
/third_party

-----
Installation

There are 2 sets of installation options.  The simplest is:
1) In terminal run:
    cd FortranAnalysis
    make clean
    make Full
This will configure and install the third party version of cfitsio in the third_party directory.  Note that this will not overwrite your existing cfitsio installation.  If you wish to adjust the main code body afterwards, simply run 
    make clean
    make
In the FortranAnalysis directory as it is not necessary to rebuild cfitsio directory.

The second installation option is:
1) Adjust the line 10 in FortranAnalysis/makeflags to point to your cfitsio installation.
2) In terminal:
    cd FortranAnalysis
    make clean
    make
    
    
2) In terminal run:
    cd src/
    make clean
    make 

Once installed, the program executables should be found in /Programs
    
----
Running

Both codes require a single text input file.  A sample input file can be found in:
Inputs/

The input file specifies both the cube and the mask being used.  Note that the calculation requires a mask as many cubes can be low S/N.  It is possible to use the cube itself as a mask (for noiseless cubes as an example) by adjusting the mask switch in the input file.

Either of the codes can be run in terminal using:
$(PathToPrograms)/ProgramName $(PathToInputFile)/InputFile.txt

For example, in the main folder, the core code can be run:

    ./Programs/Asymmetry3D Inputs/AsymmetryInputs.in 

The codes can be run to either calculate the asymmetries using the traditional absolute value method, or using a modified square difference method.  This can be switched in the input file.

---
Inputs

There are a variety of required inputs in the input file.  The file format structure is


       Name of the Data cube
filename
    Mask file name
filename
    Mask file type (1 == normal mask, 2 == construct mask from file using some fraction of the total flux [specify in this line])
1
    The name of the output text file containing the asymmetry information
filename
    Centre Location (Only used in the Asymmetry3D calculation, but this line is always required.)
x   y   z
    Background Switch ( 1 == Estimate the asymmetry due to the background noise, 2 == no background)
2
    Asymmetry Method Switch (0==Absolute sign method, 1== Squared Difference Method)
1
    That random see to be used for any calculations (if >0, the time will be used to get a seed)
-1


Do note that all the lines are required in this specific order.



