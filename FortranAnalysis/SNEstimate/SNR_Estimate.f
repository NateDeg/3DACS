
      module SNRCalcMod
      use DataCubeMod
      use BeamMod
      use sort


      contains

cccc
      subroutine SNREstimate(Cube,MaskCube,Beam
     &              ,SNPeak,SN_Avg,SNR99)
c       This routine estimates the S/N ratio of a cube using a variety of methods
c           It assumes that the cube has already been masked
      implicit none
      Type(DataCube), INTENT(IN) :: Cube,MaskCube
      real,INTENT(INOUT) :: SNPeak,SN_Avg,SNR99
      Type(Beam2D), INTENT(IN) :: Beam

      real TotFlux, FluxSum
      integer TotCells

      real,ALLOCATABLE :: FlatFlux(:)
      integer,ALLOCATABLE :: FlatIndx(:)
      integer i, j,count


      print*, "Estimate S/N"


c       Get the peak S/N
      SNPeak=maxval(Cube%MaskedFlux)/Cube%DH%Uncertainty
c       Get the average S/N
      TotFlux=sum(Cube%MaskedFlux)
      TotCells=int(Sum(MaskCube%Flux))
      SN_Avg=TotFlux/TotCells

c       Allocate the flattened flux array
c      ALLOCATE(FlatFlux(TotCells))
c      ALLOCATE(FlatIndx(TotCells))
c       Loop through the cells and fill the flattened flux array
c      call MakeMaskedFlatFluxArr(Cube,MaskCube,FlatFlux,TotCells)
c       Sort the flattened array
c      call indexx(TotCells,FlatFlux,FlatIndx)
c       Go through the flattened array until we reach 99% of the total flux
c      fluxSum=0.
c      do i=TotCells,1,-1
c        j=FlatIndx(i)
c        FluxSum=fluxSum+FlatFlux(j)
c        if(fluxSum .gt. 0.99*TotFlux) then
c            count=i
c            goto 100
c        endif
c      enddo
c100   continue
c       Now get the S/N_99 using the equation from Appendix A of Westmeier et al. 2021

      SNR99=TotFlux/((sqrt(TotCells*Beam%BeamAreaPixels))
     &                  *Cube%DH%Uncertainty)
c      print*, "SN99 Test", SNR99, FluxSum,count,Beam%BeamAreaPixels
c     &              ,Cube%DH%Uncertainty
c     &              ,sqrt(count*Beam%BeamAreaPixels)
c     &              ,sqrt(count*Beam%BeamAreaPixels)*Cube%DH%Uncertainty

c       Remember to deallocate the flattened flux
c      DEALLOCATE(FlatFlux)
c      DEALLOCATE(FlatIndx)


      return
      end subroutine
ccccc


ccccc
      subroutine MakeMaskedFlatFluxArr(Cube,MaskCube,FlatFlux,TotCells)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube,MaskCube
      integer, INTENT(IN) :: TotCells
      real,INTENT(INOUT) :: FlatFlux(TotCells)
      integer i,j,k,count

c       Initialize the counter
      count=0
      do i=0,Cube%DH%nPixels(0)-1
        do j=0,Cube%DH%nPixels(1)-1
            do k=0,Cube%DH%nChannels-1
                if(MaskCube%Flux(i,j,k) .gt. 0.) then
                    count=count+1
                    FlatFlux(count)=Cube%Flux(i,j,k)
                endif
            enddo
        enddo
      enddo
      return
      end subroutine
cccccccccc

ccccc
      subroutine MakeIndexedFlatFluxArr(Cube,FlatFlux,TotCells,Indx)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube
      integer, INTENT(IN) :: TotCells
      integer,INTENT(INOUT) :: Indx(TotCells,3)
      real,INTENT(INOUT) :: FlatFlux(TotCells)
      integer i,j,k,count

c       Initialize the counter
      count=0
      do i=0,Cube%DH%nPixels(0)-1
        do j=0,Cube%DH%nPixels(1)-1
            do k=0,Cube%DH%nChannels-1
                count=count+1
                FlatFlux(count)=Cube%Flux(i,j,k)
                Indx(count,1)=i
                Indx(count,2)=j
                Indx(count,3)=k
            enddo
        enddo
      enddo
      return
      end subroutine
cccccccccc


      end module

