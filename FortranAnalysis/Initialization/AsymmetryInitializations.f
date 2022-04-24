cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for processing
c       and initializing necessary quantities for the asymmetry
c       calculations
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module AsymmetryIniMod

      use AsymmetryGlobals
      use EstimateCubeNoiseMod
      use DataCubeMod

      implicit none

      contains



ccccc
c     This routine initializes all quantities needed for the cube analysis
      subroutine Asym3DIni()
      implicit none
      integer timearray(3)

      print*, "Process and Initializes objects"
c       First set the background asymmetry flag switch
      if(BackgroundSwitch .eq. 1) then
        ObservedDC%DA%BackgroundAsymFlag=.True.
      else
        ObservedDC%DA%BackgroundAsymFlag=.False.
      endif
c       Next set the asymmetry method switch
      ObservedDC%DA%AsymMethodSwitch=AsymMethodSwitch
c       Now set the symmetric mask switch
      ObservedDC%DH%SymmetricMaskSwitch=1
c       Next estimate the noise in the cube
      call EstimateNoise(ObservedDC
     &                  ,ObservedDC%DH%Uncertainty)
c       Now mask the data
      call MaskData()

c       Finally set the random seed
      if(idum .ge. 0) then
        call itime(timeArray)
        idum=abs(timeArray(1)*idum)
        idum=-int(idum*timeArray(2))-timeArray(3)
      endif

      return
      end subroutine
ccccccccc

cc
c       This subroutine masks the data used for the asymmetry calculation
      subroutine MaskData()
      implicit none
c
c       Before masking the data, make sure the mask data cube is zeros and ones
c       Depending on the mask option, the cube mask may need to be adjusted
      if(MaskFileTypeSwitch .eq. 2) then
c           For this switch, the mask is a noiseless cube and set so that
c           99% of the total flux is contained in the mask
        call AdjustMaskByFlux()
      endif
c       Once the mask is adjusted, mask the data cube
      call MaskCube(DataCubeMask,ObservedDC)

      return
      end subroutine
ccccccc




cccc
c       This routine takes the masked cube and sets the cells that sum
c           up to 99% of the flux to 1 and the rest to 0
      subroutine AdjustMaskByFlux()
      use SNRCalcMod
      use sort
      implicit none

      integer i, j, k,ii,jj


      real TotFlux,fSum,fTerm
      integer TotCells
      real,ALLOCATABLE :: FlatFlux(:)
      integer,ALLOCATABLE :: FlatIndx(:),ThreeDIndx(:,:)

      real FTotFrac

c       First get the total flux and total number of cells
      TotFlux=sum(DataCubeMask%Flux)
      TotCells=DataCubeMask%DH%nPixels(0)*DataCubeMask%DH%nPixels(1)
     &          *DataCubeMask%DH%nChannels
c      print*, TotFlux, TotCells
c       Now allocate the flattened arrays
      ALLOCATE(FlatFlux(TotCells))
      ALLOCATE(FlatIndx(TotCells))
      ALLOCATE(ThreeDIndx(TotCells,3))

c       Loop through the cells and fill the flattened flux array
      call MakeIndexedFlatFluxArr(DataCubeMask,FlatFlux
     &          ,TotCells,ThreeDIndx)       !/SNEstimate/SNR_Estimate.f
c       Sort the flattened array
      call indexx(TotCells,FlatFlux,FlatIndx)   !/StandardMath/SortArray.f


      FTotFrac=MaskFluxFrac
c       Loop through all cells by brightness
      fSum=0.
      do ii=TotCells,1,-1
c           Assign the indices
        jj=FlatIndx(ii)
        i=ThreeDIndx(jj,1)
        j=ThreeDIndx(jj,2)
        k=ThreeDIndx(jj,3)
c           Sum up the flux
        fSum=fSum+FlatFlux(jj)
        if(fSum .lt. FTotFrac*TotFlux) then
            fTerm=1.
        else
            fTerm=0.
        endif
c           Reassign the flux
        DataCubeMask%Flux(i,j,k)=fTerm
      enddo
c       Clear the memory
      DEALLOCATE(FlatFlux)
      DEALLOCATE(FlatIndx)
      DEALLOCATE(ThreeDIndx)

      return
      end subroutine
ccccccc


cccc
c           This routine masks a cube
c
      subroutine MaskCube(MaskUsed,WorkingDC)
      implicit none
      Type(DataCube),INTENT(IN) :: MaskUsed
      Type(DataCube),INTENT(INOUT) :: WorkingDC
      integer i,j,k
c           Loop through the cube
      do i=0,WorkingDC%DH%nPixels(0)-1
        do j=0,WorkingDC%DH%nPixels(1)-1
            do k=0,WorkingDC%DH%nChannels-1
                WorkingDC%MaskedFlux(i,j,k)=WorkingDC%Flux(i,j,k)
     &                  *MaskUsed%Flux(i,j,k)
            enddo
        enddo
      enddo

      return
      end subroutine
cccccc



cccccc
c       This routine calculates the Mom0 map from the observed datacube object
c           The Mom0 map is stored in the ObservedMap object
      subroutine QuickMom0Construction()
      implicit none
      integer i,j,k
c       Copy the observed data cube header to the map
      ObservedMap%DH=ObservedDC%DH
c       The map is only 1 channel wide
      ObservedMap%DH%nChannels=1
c       Also copy over the asymmetry object to keep the rotation point
      ObservedMap%DA=ObservedDC%DA
c       The channel rotation point is 0
      ObservedMap%DA%RotationPoint(3)=0
c       The mask can't be adjusted for a moment 0 map
      ObservedMap%DH%SymmetricMaskSwitch=0
c       Allocate the map
      call AllocateDataCube(ObservedMap)

c           Loop through all pixels
      do i=0,ObservedDC%DH%nPixels(0)-1
        do j=0,ObservedDC%DH%nPixels(1)-1
c               Sum up the cube flux along the channel axis
            ObservedMap%Flux(i,j,0)=sum(ObservedDC%MaskedFlux(i,j,:))
        enddo
      enddo
      ObservedMap%MaskedFlux=ObservedMap%Flux
c      print*, "Map Flux Check", sum(ObservedMap%Flux(:,:,0))
c      print*, "Cube Flux Check",sum(ObservedDC%Flux)

      return
      end subroutine
cccccccc

cccccc
c       This routine calculates the profile from the observed datacube object
c           The profile is stored in the ObservedProfile object
      subroutine QuickProfileConstruction()
      implicit none
      integer i

c       Copy the observed data cube header to the profile
      ObservedProfile%DH=ObservedDC%DH
c       The profile is only 1 pixel in size
      ObservedProfile%DH%nPixels(0:1)=1
c       Also copy over the asymmetry object to keep the rotation point
      ObservedProfile%DA=ObservedDC%DA
c       The pixel rotation point is 0
      ObservedProfile%DA%RotationPoint(1:2)=0
c       The mask can't be adjusted for a profile
      ObservedProfile%DH%SymmetricMaskSwitch=0
c           Allocate the profile
      call AllocateDataCube(ObservedProfile)
c           Loop over all channels
      do i=0,ObservedDC%DH%nChannels-1
c           Sum the flux in all pixels for a given channel
        ObservedProfile%Flux(0,0,i)=sum(ObservedDC%MaskedFlux(:,:,i))
      enddo
      ObservedProfile%MaskedFlux=ObservedProfile%Flux
c      print*,"Profile Flux Check", sum(ObservedProfile%Flux(0,0,:))
      return
      end subroutine
cccccccc



      end module
