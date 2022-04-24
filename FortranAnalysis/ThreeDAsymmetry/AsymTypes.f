cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for calculating
c       the 3D asymmetry
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module AsymTypesMod

      use AsymmetryGlobals

      implicit none

      contains

cccc
      subroutine AbsoluteAsymPtSum(Asym,F1,F2)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      real, INTENT(IN) :: F1,F2
      Asym%TotAbsDiff=Asym%TotAbsDiff+abs(F1-F2)
      Asym%TotFlux=Asym%TotFlux+F1+F2
      return
      end subroutine
ccccc

cccc
      subroutine AbsoluteAsymTermComb(Asym)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      Asym%Asym=Asym%TotAbsDiff/Asym%TotFlux
      return
      end subroutine
ccccc

cccc
      subroutine AbsoluteMeasAsymCalc(Asym)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      Asym%Asym=Asym%Signal_Asym
     &                      -Asym%Back_Asym
      return
      end subroutine
ccccc


cccc
      subroutine SquaredAsymPtSum(Asym,F1,F2)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      real, INTENT(IN) :: F1,F2
      Asym%TotAbsDiff=Asym%TotAbsDiff+(F1-F2)**2.
      Asym%TotFlux=Asym%TotFlux+(F1+F2)**2.
      return
      end subroutine
ccccc

cccc
      subroutine SquaredAsymTermComb(Asym)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      Asym%TotAbsDiff=sqrt(Asym%TotAbsDiff)
      Asym%TotFlux=sqrt(Asym%TotFlux)
      Asym%Asym=Asym%TotAbsDiff/Asym%TotFlux
      return
      end subroutine
ccccc

cccc
      subroutine SquaredMeasAsymCalc(Asym)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      Asym%Asym=sqrt(Asym%Signal_Asym**2.
     &                      -Asym%Back_Asym**2.)
      return
      end subroutine
ccccc


      end module
