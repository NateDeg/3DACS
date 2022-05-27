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
      use ThreeDAsymCoreMod
      use SNRCalcMod
      use AsymmetryOutputsMod
      use AsymTypesMod
      implicit none

      real SNPeak,SNA,SN_Int
      integer i,j,k
      Type(DataCubeAsymmetry) MapAsym,ProfileAsym
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "3D Asymmetry Analysis"

c       Get the main inputs
      call Asym3DInput()
c       Initialize everything needed
      call Asym3DIni()

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
c           Symmeterize the mask about the rotation point
      if(ObservedDC%DH%SymmetricMaskSwitch.eq. 1) then
        call MakeSymmetricMask(ObservedDC%DA%RotationPoint
     &          ,DataCubeMask,SymmetricMask)
c           Remask the data
        call MaskCube(SymmetricMask,ObservedDC)

        call QuickMom0Construction()

        call QuickProfileConstruction()
      else
        call QuickMom0Construction()

        call QuickProfileConstruction()

      endif

c      call MakeSymmetricMask(ObservedDC%DA%RotationPoint
c     &          ,DataCubeMask,SymmetricMask)
cc           Remask the data
c      call MaskCube(SymmetricMask)

c           Do the full 3D asymmetry calculation
      call ThreeDAsym_CenterSet(ObservedDC)

c           Do the moment map signal asymmetry
      call GetSignalAsym(ObservedMap,MapAsym)
      ObservedMap%DA%Signal_Asym=MapAsym%Asym
      ObservedMap%DA%Asym=MapAsym%Asym
      ObservedMap%DA%TotAbsDiff=MapAsym%TotAbsDiff
      ObservedMap%DA%TotFlux=MapAsym%TotFlux


c           Also do the moment map signal asymmetry
      call GetSignalAsym(ObservedProfile,ProfileAsym)
      ObservedProfile%DA%Signal_Asym=ProfileAsym%Asym
      ObservedProfile%DA%TotAbsDiff =ProfileAsym%TotAbsDiff
      ObservedProfile%DA%TotFlux=ProfileAsym%TotFlux


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





