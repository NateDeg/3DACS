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
     &                      -Asym%Back_Asym/Asym%TotFlux
      return
      end subroutine
ccccc

cccc
      subroutine AbsoluteBackAsymCalc(BackSum,nPairs,RMS)
c      use BasicRanNumGen
      implicit none
      integer, INTENT(IN) :: nPairs
      real, INTENT(IN) :: RMS
      real, INTENT(INOUT) :: BackSum
c      integer i
c      real T1,T2

      BackSum=real(nPairs)*2.*RMS/sqrt(3.1415)
c      BackSum=0.
c      do i=1, nPairs
c        T1=gasdev(idum)*RMS
c        T2=gasdev(idum)*RMS
c        BackSum=BackSum+abs(T1-T2)
c      enddo
c      print*, "AbsoluteBackground summation",BackSum,RMS
c     &          , 2.*nPairs*RMS
c     &          , 2.*nPairs*RMS/sqrt(3.1415921)
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
c      Asym%TotAbsDiff=sqrt(Asym%TotAbsDiff)
c      Asym%TotFlux=sqrt(Asym%TotFlux)
      Asym%Asym=sqrt(Asym%TotAbsDiff/Asym%TotFlux)
      return
      end subroutine
ccccc

cccc
      subroutine SquaredMeasAsymCalc(Asym)
      implicit none
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
      real numer, denom
      Asym%Asym=sqrt(Asym%Signal_Asym**2.
     &                      -Asym%Back_Asym**2.)

      numer=Asym%TotAbsDiff-Asym%Back_Asym
c      denom=Asym%TotFlux-Asym%Back_Asym
      denom=Asym%TotFlux-Asym%Back_Asym

      Asym%Asym=sqrt(numer/denom)
      print*, "Squared Diff Check",Asym%TotAbsDiff
     &          ,Asym%TotFlux,Asym%Back_Asym
     &          , Asym%Asym


      return
      end subroutine
ccccc

ccccc
      subroutine SquaredBackAsymCalc(BackSum,nPairs,RMS)
      use BasicRanNumGen
      implicit none
      integer, INTENT(IN) :: nPairs
      real, INTENT(IN) :: RMS
      real, INTENT(INOUT) :: BackSum
      integer i
      real T1,T2

      BackSum=2.*real(nPairs)*RMS**2.
      print*, "BackSum Calculation", nPairs, RMS,nPairs/2
c      BackSum=0.
c      do i=1, nPairs
c        T1=gasdev(idum)*RMS
c        T2=gasdev(idum)*RMS
c        BackSum=BackSum+(T1-T2)**2.
c      enddo
c      print*, "Squared Background summation",BackSum,RMS
c     &          , 2.*nPairs*RMS**2.
c     &          , 4.*nPairs*RMS**2./3.1415921

      return
      end subroutine
cccccc


      end module
