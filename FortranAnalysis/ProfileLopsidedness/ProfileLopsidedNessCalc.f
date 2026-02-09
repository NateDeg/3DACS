cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for calculating
c       the 3D asymmetry
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ProfileLopsidednessMod

      use AsymmetryGlobals
      use DataCubeMod
      use InterpolateMod
      use AsymTypesMod
      use SimpleSmoothMod
      use IntegrationMod

      implicit none
      contains

cccc

cccccc
c
      subroutine CalcProfileLopsidedness(DC)
      implicit none
      Type(DataCube),INTENT(INOUT) ::  DC

      real VSys
      integer nBins(2)
      real IntFlux(2)

      integer nTot

ccc
c      print*, "Calc Profile Lopsidedness"

c       Start by getting VSys
      call EstimateVSys(VSys,DC)
c       Now get the integraded flux
      call CalcIntFlux(nBins,IntFlux,VSys,DC)
c       Now we have all the pieces we need to get the
c       lopsidedness using the normal 3D asymmetry equations
      nTot=sum(nBins)
      DC%DA%LopRMS=sqrt(real(nTot)/2.)
     &      *DC%DH%ExpectedUncertaintySquared

      DC%DA%LopP=(IntFlux(1)-IntFlux(2))**2.
      DC%DA%LopQ=(IntFlux(1)+IntFlux(2))**2.
      DC%DA%LopB=2.*DC%DA%LopRMS**2.
      DC%DA%LopC=sqrt(DC%DA%LopP/DC%DA%LopQ)
      DC%DA%LopA=sqrt((DC%DA%LopP-DC%DA%LopB)
     &          /(DC%DA%LopQ-DC%DA%LopB))
      if (DC%DA%LopP-DC%DA%LopB .lt. 0.) then
        DC%DA%LopA=-1.
      endif
      DC%DA%LopVSys=VSys
      DC%DA%LopChan=(VSys-DC%DH%RefVal(2))
     &          /(DC%DH%ChannelSize)
     &          +DC%DH%RefLocation(2)


      return
      end subroutine
ccccc

cccc
c
      subroutine CalcIntFlux(nBins,IntFlux,VSys,DC)
      implicit none
      integer, INTENT(INOUT) :: nBins(2)
      real, INTENT(INOUT) :: IntFlux(2)
      real, INTENT(IN) :: VSys
      Type(DataCube), INTENT(IN) :: DC
      real IntLims(2), CurrVel,CurrFlux
      integer i
ccc


c       Start be doing the lower velocity integral
c       To do this, we need to get the lower limit
c        as well as the number of unmasked channels
      nBins=0
      do i=0,DC%DH%nChannels-1
        CurrFlux=DC%Flux(0,0,i)
        CurrVel=DC%Channels(i)
        if (CurrVel .lt. VSys .and.
     &          CurrFlux .ne. 0.) then
            nBins(1)=nBins(1)+1
            if (nBins(1) .eq. 1) then
                IntLims(1)=CurrVel
            endif
        elseif(CurrVel .gt. VSys) then
            goto 100
        endif
      enddo
100   continue
      IntLims(2)=VSys
c      print*, "Bin check1", nBins,IntLims


      call TabulatedTrapzIntegral(IntFlux(1)
     &          ,DC%DH%nChannels
     &      ,DC%Channels(0:DC%DH%nChannels-1)
     &      ,DC%Flux(0,0,0:DC%DH%nChannels-1),IntLims)
c      print*, "Low Flux Check", IntFlux(1)


c       Now do the upper velocity integral
c           Again, we need to get the upper limit
c           number of unmasked channels
c      do i=0,DC%DH%nChannels-1
c        print*, "Flux Check", i,DC%Flux(0,0,i)
c      enddo
      do i=DC%DH%nChannels-1,0,-1
        CurrFlux=DC%Flux(0,0,i)
        CurrVel=DC%Channels(i)
        if (CurrVel .gt. VSys .and.
     &          CurrFlux .ne. 0.) then
c            print*, "High Loop",i
c     &      , CurrVel,CurrFlux,nBins(2)
            nBins(2)=nBins(2)+1
            if (nBins(2) .eq. 1) then
c                print*, "Upper check",i
                IntLims(2)=CurrVel
            endif
        elseif(CurrVel .lt. VSys) then
            goto 200
        endif
      enddo
200   continue
      IntLims(1)=VSys
c      print*, "Bin check2", nBins, IntLims

      call TabulatedTrapzIntegral(IntFlux(2)
     &          ,DC%DH%nChannels
     &      ,DC%Channels(0:DC%DH%nChannels-1)
     &      ,DC%Flux(0,0,0:DC%DH%nChannels-1),IntLims)
c      print*, "High Flux Check", IntFlux(2)
c      print*, "Int Check", IntFlux, sum(DC%Flux)
c      print*, sum(IntFlux/(DC%Channels(1)-DC%Channels(0)))
c       Renormalize the integrated flux by the total flux
c       This is simply to make the calculation of the RMS later easier
      IntFlux=IntFlux/(sum(IntFlux))*sum(DC%Flux)
c      print*, IntFlux

      return
      end subroutine


ccccc
c
      subroutine EstimateVSys(VSys,DC)
      implicit none
      real, INTENT(INOUT) :: VSys
      Type(DataCube),INTENT(IN) :: DC

      real VSysSmooth(3)
      integer nSmooth, i
      real,ALLOCATABLE :: SmoothedFlux(:)
ccc
      ALLOCATE(SmoothedFlux(0:DC%DH%nChannels-1))


      do i=1,3
        nSmooth=2*i+1
        call Simple1DSmooth(DC%DH%nChannels
     &      ,nSmooth
     &      ,DC%Flux(0,0,0:DC%DH%nChannels-1)
     &      ,SmoothedFlux(0:DC%DH%nChannels-1))

        call EstimateVSysFromCumSum(
     &       DC%DH%nChannels
     &      ,DC%Channels(0:DC%DH%nChannels-1)
     &      ,SmoothedFlux(0:DC%DH%nChannels-1)
     &      ,VSysSmooth(i))
      enddo
c      print*, VSysSmooth
c      print*, (VSysSmooth-DC%DH%RefVal(2))
c     &          /(DC%DH%ChannelSize)
c     &          +DC%DH%RefLocation(2)
      VSys=sum(VSysSmooth)/3.
c      print*, VSys

      DEALLOCATE(SmoothedFlux)
      return
      end subroutine
ccccc

cccc
c
      subroutine EstimateVSysFromCumSum(nChannels
     &          ,Vels,Flux,VSys)
      implicit none
      real, INTENT(INOUT) :: VSys
      integer, INTENT(IN) :: nChannels
      real, INTENT(IN) :: Vels(0:nChannels-1)
     &              ,Flux(0:nChannels-1)
      real TotFlux
      real CumNormFlux(0:nChannels-1)
      real VLow,VHigh
      integer i
      real TargFrac
cccc

      TargFrac=0.1
      TotFlux=sum(Flux)
      CumNormFlux(0)=Flux(0)/TotFlux
      do i=1, nChannels-1
        CumNormFlux(i)=CumNormFlux(i-1)+Flux(i)/TotFlux
c        print*, i, CumNormFlux(i)
        
        if (CumNormFlux(i) .ge. TargFrac
     &      .and. CumNormFlux(i-1) .le. TargFrac) then
            VLow=Vels(i)
        endif
        if (CumNormFlux(i) .ge. 1.-TargFrac
     &      .and. CumNormFlux(i-1) .le. 1.-TargFrac) then
            VHigh=Vels(i)
        endif
      enddo
c      print*, CumNormFlux
c      print*, VLow,VHigh
      VSys=(VLow+VHigh)/2.
c      print*, "VSys estimate", VSys

      return
      end subroutine

cccc
      

      end module
