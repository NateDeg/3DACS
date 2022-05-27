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
      use AsymTypesMod

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

c       SET THE SWITCH FOR ABSOLUTE OR SQUARED DIFFERENCE ASYMMETRY
      if (ObservedDC%DA%AsymMethodSwitch .eq. 0) then
        AsymPtPoint => AbsoluteAsymPtSum
        AsymCombPoint => AbsoluteAsymTermComb
        MeasAsymPoint => AbsoluteMeasAsymCalc
        BackAsymPoint => AbsoluteBackAsymCalc
      elseif (ObservedDC%DA%AsymMethodSwitch .eq. 1) then
        AsymPtPoint => SquaredAsymPtSum
        AsymCombPoint => SquaredAsymTermComb
        MeasAsymPoint => SquaredMeasAsymCalc
        BackAsymPoint => SquaredBackAsymCalc
      endif


c       Since we have the noise and beam and have masked the cube, we can
c           get some S/N measures
      call SNREstimate(ObservedDC,DataCubeMask,Beam
     &                    ,SNPeak,SNA,SN_Int)

      print*, "Symmetric Flag",ObservedDC%DH%SymmetricMaskSwitch
c       The minimzer can't do symmetric masking every step with the 1D and 2D, so the best idea is to use the asymmetric mask to make the moment maps then symmetrize when doing the 3D minimization.  This is built into the minimzer, so we don't need to do any work on that here.

      call QuickMom0Construction()
      call QuickProfileConstruction()



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
      print*, "Total Numerator =", ObservedDC%DA%TotAbsDiff
      print*, "Total Denominator =", ObservedDC%DA%TotFlux


      print*, " "
      print*, "2D Asymmetry calculated about point:"
     &              ,ObservedMap%DA%RotationPoint(1:2)
      print*, "2D Signal Asymmetry", ObservedMap%DA%Asym
      print*, "Total Numerator =", ObservedMap%DA%TotAbsDiff
      print*, "Total Denominator =", ObservedMap%DA%TotFlux

      print*, " "
      print*, "1D Asymmetry calculated about point:"
     &              ,ObservedProfile%DA%RotationPoint(3)
      print*, "1D Signal Asymmetry", ObservedProfile%DA%Signal_Asym
      print*, "Total Numerator =", ObservedProfile%DA%TotAbsDiff
      print*, "Total Denominator =", ObservedProfile%DA%TotFlux


      call OutputAsymmetry(SNPeak,SNA,SN_Int)

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

