cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines to calculate the moments 
c     of some a passed array of n-body particles per 
c     pixel element
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ProfileMod
      implicit none

      Type ProfileHeader
        integer nChannels
        real ChannelSize
        real RefChannel, RefChannelValue
        real Start

        integer nSmooth
        real CubeRMS

      end Type

      Type ProfileCharacteristics
        real VSys
        real PeakFluxFrac
        real VEdges(2),FluxEdges(2)
        integer EdgeChannelLocs(2),VSysChannelLoc
        integer PeakChannels(2)
        real PeakFlux(2),PeakVel(2)
        integer nPeaks
        real FluxWeightedVel
        logical BadProfileFlag
        real RMS, SN_Peak
      end Type

      Type ProfileAsymmetry
        logical BackgroundAsymFlag
        real AsymMeas,SignalAsym,BackAsym
        real Asym_Sys,Signal_Sys,BackAsym_Sys
        real DiffSum,FluxSum,BackSum
        real V_Min
        
        real Lopsidedness,LowVelIntegral,HighVelIntegral
        real DelV

      end Type

      Type Profile
        Type(ProfileHeader) PH
        Type(ProfileCharacteristics) PC
        Type(ProfileAsymmetry) PA
        real,dimension(:),ALLOCATABLE :: Channels
        real,dimension(:),ALLOCATABLE :: Flux
        real,dimension(:),ALLOCATABLE :: SmoothedFlux
        real,dimension(:),ALLOCATABLE :: Slope
        real,dimension(:),ALLOCATABLE :: NoiseFlux
        integer,dimension(:),ALLOCATABLE :: nPixelsInChan
      end Type

      contains
ccccccc
      subroutine SimpleAllocateProfile(P)
c
      implicit none
      Type(Profile),INTENT(INOUT) :: P
      integer i, j,mPix
      real Delta

c       Initialize profile arrays
      ALLOCATE(P%Channels(0:P%PH%nChannels-1))
      ALLOCATE(P%Flux(0:P%PH%nChannels-1))
      ALLOCATE(P%SmoothedFlux(0:P%PH%nChannels-1))
      ALLOCATE(P%Slope(0:P%PH%nChannels-1))
      ALLOCATE(P%NoiseFlux(0:P%PH%nChannels-1))
      ALLOCATE(P%nPixelsInChan(0:P%PH%nChannels-1))
c       Set the initial flux to zero
      P%Flux=0.


      return
      end subroutine
cccccccc

cccccc
      subroutine DeAllocateProfile(P)
      implicit none
      Type(Profile) P

      DEALLOCATE(P%Channels)
      DEALLOCATE(P%Flux)
      DEALLOCATE(P%SmoothedFlux)
      DEALLOCATE(P%Slope)
      DEALLOCATE(P%NoiseFlux)
      DEALLOCATE(P%nPixelsInChan)

      return
      end subroutine
cccccccc

      end module
