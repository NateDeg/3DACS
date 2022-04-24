ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This is a program for getting the 3D Asymmetry for some
c       cube.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program ThreeDAsymmetry
c       Modules neeeded at this level

      use AsymmetryGlobals
      use AsymmetryInputMod
      use AsymmetryIniMod
      use ThreeDAsym_MinimumMod
      use SNRCalcMod
      use AsymmetryOutputsMod
      implicit none

      real SNPeak,SNA,SN_Int
      integer i,j,k
      Type(DataCubeAsymmetry) MapAsym,ProfileAsym
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "3D Asymmetry minimum Analysis"

c       Get the main inputs
      call Asym3DInput()
c           Initialize everything needed
      call Asym3DIni()
      print*, "Done input and Initialization"
c

c       Since we have the noise and beam and have masked the cube, we can
c           get some S/N measures
      call SNREstimate(ObservedDC,DataCubeMask,Beam
     &                    ,SNPeak,SNA,SN_Int)


c           Do the full 3D asymmetry calculation
      WorkingDC=> ObservedDC
      call Calc_MinThreeDAsym()


      WorkingDC=> ObservedMap
      call Calc_MinThreeDAsym()

      WorkingDC=> ObservedProfile
      call Calc_MinThreeDAsym()

      print*, " "
      print*, "Done Asymmetry calculation"
      print*, "Cell rms and integrated S/N measurement: "
     &          ,ObservedDC%DH%Uncertainty, SN_Int
      print*, "Asymmetry calculated about point: "
     &          , ObservedDC%DA%RotationPoint
      print*, "Signal Asymmetry =", ObservedDC%DA%Signal_Asym
      print*, "Background Asymmetry =", ObservedDC%DA%Back_Asym
      print*, "Total Asymmetry =", ObservedDC%DA%Asym
      print*, " "
      print*, "2D Asymmetry calculated about point:"
     &              ,ObservedMap%DA%RotationPoint(1:2)
      print*, "2D Signal Asymmetry", ObservedMap%DA%Asym

      print*, " "
      print*, "1D Asymmetry calculated about point:"
     &              ,ObservedProfile%DA%RotationPoint(3)
      print*, "1D Signal Asymmetry", ObservedProfile%DA%Asym


      call OutputAsymmetry(SNPeak,SNA,SN_Int)

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

