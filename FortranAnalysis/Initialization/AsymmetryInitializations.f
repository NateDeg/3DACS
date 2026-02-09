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
      integer nPixel, nCells,TCells
      real CellCount

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
c      print*, "Moment Calculation Checks", ObservedMap%DH%Uncertainty

c       Allocate the map
      call AllocateDataCube(ObservedMap)

c           Loop through all pixels
      CellCount=0.
      nPixel=0
      TCells=0
      do i=0,ObservedDC%DH%nPixels(0)-1
        do j=0,ObservedDC%DH%nPixels(1)-1
c               Sum up the cube flux along the channel axis
            ObservedMap%Flux(i,j,0)=sum(ObservedDC%MaskedFlux(i,j,:))
            nCells=0
            do k=0,ObservedDC%DH%nChannels-1
                if (ObservedDC%MaskedFlux(i,j,k) .ne. 0.0) then
                    nCells=nCells+1
                endif
            enddo
            CellCount=CellCount+sqrt(real(nCells))
            TCells=TCells+nCells
            if(ObservedMap%Flux(i,j,0) .ne. 0.0) then
                nPixel=nPixel+1
            endif
        enddo
      enddo
      ObservedMap%MaskedFlux=ObservedMap%Flux
c      print*, "Total number of cells", TCells,nPixel
      ObservedMap%DH%ExpectedUncertainty
     &          =ObservedMap%DH%Uncertainty
     &          *CellCount/real(nPixel)
      ObservedMap%DH%ExpectedUncertaintySquared
     &          =ObservedMap%DH%Uncertainty
     &          *sqrt(TCells/real(nPixel))

c      ObservedMap%DH%Uncertainty=ObservedMap%DH%Uncertainty
c     &          *CellCount/real(nPixel)
c      print*, "New average uncertainty",ObservedMap%DH%Uncertainty

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
      integer i,j,k
      integer nChan, nCells,TCells
      real CellCount
      real,ALLOCATABLE :: SigmaChannel(:)
      real ConvolutionFactor

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


      ALLOCATE(SigmaChannel(0:ObservedDC%DH%nChannels-1))
      ConvolutionFactor=4.*Pi
     &          *Beam%BeamSigmaVector(0)
     &          *Beam%BeamSigmaVector(1)
c      print*, "Beam Area Check"
c     &      ,Beam%BeamSigmaVector(0)
c     &      ,Beam%BeamSigmaVector(1)
c     &      ,ConvolutionFactor

c           Loop over all channels
      CellCount=0.
      TCells=0
      nChan=0
      do i=0,ObservedDC%DH%nChannels-1
c           Sum the flux in all pixels for a given channel
        ObservedProfile%Flux(0,0,i)=sum(ObservedDC%MaskedFlux(:,:,i))
        
        nCells=0
        do j=0,ObservedDC%DH%nPixels(0)-1
            do k=0,ObservedDC%DH%nPixels(1)-1
                if (ObservedDC%MaskedFlux(j,k,i) .ne. 0.0) then
                    nCells=nCells+1
                endif
            enddo
        enddo
c        print*, nCells
        SigmaChannel(i)=(ObservedDC%DH%Uncertainty**2.
     &              *real(nCells)
     &              *ConvolutionFactor)**0.5
c        print*, i, SigmaChannel(i)
        CellCount=CellCount+sqrt(real(nCells))
        TCells=TCells+nCells
        if(ObservedProfile%Flux(0,0,i) .ne. 0.0) then
            nChan=nChan+1
        endif
      enddo
      ObservedProfile%MaskedFlux=ObservedProfile%Flux
c      print*, "Orignal uncertainty",ObservedProfile%DH%Uncertainty
c      print*, CellCount,TCells,nChan
      ObservedProfile%DH%ExpectedUncertainty
     &          =ObservedMap%DH%Uncertainty
     &          *CellCount/real(nChan)
      ObservedProfile%DH%ExpectedUncertaintySquared
     &          =ObservedMap%DH%Uncertainty
     &          *sqrt(TCells/real(nChan))


c       Get the profile uncertainty
      ObservedProfile%DH%Uncertainty=ObservedProfile%DH%Uncertainty
     &          *CellCount/real(nChan)
c       An extra factor must be added to the noise
c           to account for the beam correlation

c      print*, "Initial Estimate"
c     &      ,ObservedProfile%DH%Uncertainty
c     &      ,ObservedProfile%DH%ExpectedUncertainty
c     &  ,ObservedProfile%DH%ExpectedUncertaintySquared

      ObservedProfile%DH%Uncertainty=
     &         ObservedProfile%DH%Uncertainty
     &          *sqrt(ConvolutionFactor)

      ObservedProfile%DH%ExpectedUncertainty=
     &         ObservedProfile%DH%ExpectedUncertainty
     &          *sqrt(ConvolutionFactor)

      ObservedProfile%DH%ExpectedUncertaintySquared=
     &   ObservedProfile%DH%ExpectedUncertaintySquared
     &          *sqrt(ConvolutionFactor)

c      print*, "Convolved Estimate"
c     &      ,ObservedProfile%DH%Uncertainty
c     &      ,ObservedProfile%DH%ExpectedUncertainty
c     &  ,ObservedProfile%DH%ExpectedUncertaintySquared


c      ObservedProfile%DH%Uncertainty
c     &    =(ObservedProfile%DH%Uncertainty**2.
c     &      *(4.*Pi*Beam%BeamSigmaVector(0)
c     &          *Beam%BeamSigmaVector(1)))**0.5
c       A similar factor needs to be added to the 'expected uncertainty'
c      print*, "Ori 1D"
c     &  ,ObservedProfile%DH%ExpectedUncertainty
c     &  ,4.*Pi*Beam%BeamSigmaVector(0)
c     &          *Beam%BeamSigmaVector(1)
c     &  ,Beam%BeamSigmaVector(0)
    
c      ObservedProfile%DH%ExpectedUncertainty
c     &    =((ObservedProfile%DH
c     &          %ExpectedUncertainty)
c     &      *(4.*Pi*Beam%BeamSigmaVector(0)
c     &          *Beam%BeamSigmaVector(1)))
c       And to the expected uncertainty squared

c      ObservedProfile%DH%ExpectedUncertaintySquared
c     &    =ObservedProfile%DH%ExpectedUncertainty**2.


c      print*, "Avg Noise",
c     &          sqrt(sum(SigmaChannel**2.)/real(nChan))

c      print*, "Simple Check", ObservedDC%DH%Uncertainty
c     &          ,ObservedDC%DH%Uncertainty
c     &          *sqrt(CellCount*ConvolutionFactor)
c     *          /sqrt(real(nChan))
      ObservedProfile%DH%ExpectedUncertaintySquared=
     &           ObservedDC%DH%Uncertainty
     &          *sqrt(TCells*ConvolutionFactor)
     &          /sqrt(real(nChan))
c      ObservedProfile%DH%ExpectedUncertaintySquared=
c     &   sqrt(sum(SigmaChannel**2.))/real(nChan)

c      open(20,file="TestProfile.txt")
c      do i=0,ObservedDC%DH%nChannels-1
c        write(20,*) i
c     &          , ObservedProfile%MaskedFlux(0,0,i)
c     &          , SigmaChannel(i)
c     & ,ObservedProfile%DH%ExpectedUncertaintySquared
c     & ,ObservedProfile%DH%ExpectedUncertainty
c      enddo
c      close(20)

      DEALLOCATE(SigmaChannel)

c      print*,Beam%BeamSigmaVector
c      print*,"Profile Flux Check", sum(ObservedProfile%Flux(0,0,:))
      return
      end subroutine
cccccccc



      end module
