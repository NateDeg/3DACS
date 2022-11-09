3D Asymmetries in data CubeS

version: 1.
----
Description:

This package contains 2 codes for calculating the 3D asymmetry.  The first is:
Asymmetry3D
which calculates the 1D, 2D, and 3D asymmetry of a cube (plus mask) about a specific point.

The second code is:
Min_Asymmetry3D
which attempts to find the point that minimizes the asymmetry for the 1D, 2D, and 3D calculations.

    **THE MINIMIZATION CODE IS STILL UNDER DEVELOPMENT AND IS NOT CURRENTLY FUNCTIONAL
    



----
authors: N. Deg

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





