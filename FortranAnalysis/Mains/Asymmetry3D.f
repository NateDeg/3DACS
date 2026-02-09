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
      use ProfileLopsidednessMod
      implicit none

      real SNPeak,SNA,SN_Int,SN_Med
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
     &                    ,SNPeak,SNA,SN_Int,SN_Med)

      print*, "Symmetric Flag",ObservedDC%DH%SymmetricMaskSwitch
c      ObservedDC%DH%SymmetricMaskSwitch=0
c           Symmeterize the mask about the rotation point
      if(ObservedDC%DH%SymmetricMaskSwitch.eq. 1) then
        call MakeSymmetricMask(ObservedDC%DA%RotationPoint
     &          ,DataCubeMask,SymmetricMask)
c           Remask the data
c        SymmetricMask=DataCubeMask
        call MaskCube(SymmetricMask,ObservedDC)

        call QuickMom0Construction()

        call QuickProfileConstruction()
      else
        call QuickMom0Construction()

        call QuickProfileConstruction()

      endif

      if (ObservedDC%DA%AsymMethodSwitch .eq. 0) then
        ObservedProfile%DH%Uncertainty
     &          =ObservedProfile%DH%ExpectedUncertainty
        ObservedMap%DH%Uncertainty
     &          =ObservedMap%DH%ExpectedUncertainty
      elseif (ObservedDC%DA%AsymMethodSwitch .eq. 1) then
        ObservedProfile%DH%Uncertainty
     &          =ObservedProfile%DH%ExpectedUncertaintySquared
        ObservedMap%DH%Uncertainty
     &          =ObservedMap%DH%ExpectedUncertaintySquared
      endif

c      call MakeSymmetricMask(ObservedDC%DA%RotationPoint
c     &          ,DataCubeMask,SymmetricMask)
cc           Remask the data
c      call MaskCube(SymmetricMask)

c           Do the full 3D asymmetry calculation
      call ThreeDAsym_CenterSet(ObservedDC)

c           Do the moment map signal asymmetry
      call ThreeDAsym_CenterSet(ObservedMap)

c           Also do the moment map signal asymmetry
      call ThreeDAsym_CenterSet(ObservedProfile)

c           Add in a calculation for the profile lopsidedness
      call CalcProfileLopsidedness(ObservedProfile)


      print*, " "
      print*, "Done Asymmetry calculation"
      print*, "Cell rms and integrated S/N measurement: "
     &          ,ObservedDC%DH%Uncertainty, SN_Int,SN_Med
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
      print*, "2D RMS Estimate: "
     &          ,ObservedMap%DH%Uncertainty
      print*, "2D Total Asymmetry", ObservedMap%DA%Asym
      print*, "2D Signal Asymmetry =", ObservedMap%DA%Signal_Asym
      print*, "2D Background Asymmetry =", ObservedMap%DA%Back_Asym
      print*, "Total Numerator =", ObservedMap%DA%TotAbsDiff
      print*, "Total Denominator =", ObservedMap%DA%TotFlux

      print*, " "
      print*, "1D Asymmetry calculated about point:"
     &              ,ObservedProfile%DA%RotationPoint(3)
      print*, "1D RMS Estimate: "
     &          ,ObservedProfile%DH%Uncertainty
      print*, "1D Total Asymmetry", ObservedProfile%DA%Asym
      print*, "1D Signal Asymmetry =", ObservedProfile%DA%Signal_Asym
      print*, "1D Background Asymmetry =", ObservedProfile%DA%Back_Asym
      print*, "Total Numerator =", ObservedProfile%DA%TotAbsDiff
      print*, "Total Denominator =", ObservedProfile%DA%TotFlux

      print*, " "
      print*, "1D Lopsided calculated about point:"
     &              ,ObservedProfile%DA%LopChan
      print*, "1D lopsided RMS Estimate: "
     &          ,ObservedProfile%DA%LopRMS
      print*, "1D Total lopsidedness"
     &      , ObservedProfile%DA%LopA
      print*, "1D Signal lopsidedness ="
     &      , ObservedProfile%DA%LopC
      print*, "1D Background lopsidedness ="
     &      , ObservedProfile%DA%LopB
      print*, "Total lopsidedness Numerator ="
     &      , ObservedProfile%DA%LopP
      print*, "Total lopsidedness Denominator ="
     &      , ObservedProfile%DA%LopQ


      call OutputAsymmetry(SNPeak,SNA,SN_Int,SN_Med)

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





