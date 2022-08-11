cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for calculating
c       the 3D asymmetry
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ThreeDAsymCoreMod

      use AsymmetryGlobals
      use DataCubeMod
      use InterpolateMod
      use AsymTypesMod

      implicit none


      PROCEDURE(AsymPtSumInterface)
     &              ,POINTER :: AsymPtPoint =>null()
      PROCEDURE(AsymTermCombineInterface)
     &              ,POINTER :: AsymCombPoint =>null()
      PROCEDURE(MeasAsymInterface)
     &              ,POINTER :: MeasAsymPoint =>null()
      PROCEDURE(BackgroundInterface)
     &              ,POINTER :: BackAsymPoint =>null()


      ABSTRACT INTERFACE
        subroutine AsymPtSumInterface(
     &              Asym,F1,F2)
            import :: DataCubeAsymmetry
            IMPLICIT NONE
            Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
            real, INTENT(IN) :: F1,F2
        END subroutine AsymPtSumInterface
      END INTERFACE

      ABSTRACT INTERFACE
        subroutine AsymTermCombineInterface(Asym)
            import :: DataCubeAsymmetry
            IMPLICIT NONE
            Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
        END subroutine AsymTermCombineInterface
      END INTERFACE

      ABSTRACT INTERFACE
        subroutine MeasAsymInterface(Asym)
            import :: DataCubeAsymmetry
            IMPLICIT NONE
            Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym
        END subroutine MeasAsymInterface
      END INTERFACE

      ABSTRACT INTERFACE
        subroutine BackgroundInterface(BackSum,nPairs,RMS)
            IMPLICIT NONE
            integer, INTENT(IN) :: nPairs
            real, INTENT(IN) :: RMS
            real, INTENT(INOUT) :: BackSum
        END subroutine BackgroundInterface
      END INTERFACE

      contains

cccc


cccccc
c
      subroutine ThreeDAsym_CenterSet(DC)
      implicit none

      integer Extent(3)
      Type(DataCubeAsymmetry) CubeAsym
      real BackSum
      Type(DataCube),INTENT(INOUT) ::  DC


      print*, "Core 3D Asymmetry"
c
c      Get the signal asymmetry value
      call GetSignalAsym(DC,CubeAsym)
c       Set the cube signal asymmetry value
      DC%DA%Signal_Asym=CubeAsym%Asym
      DC%DA%TotAbsDiff=CubeAsym%TotAbsDiff
      DC%DA%TotFlux=CubeAsym%TotFlux
c
      if(DC%DA%BackgroundAsymFlag) then
        call BackAsymPoint(BackSum,CubeAsym%nPairs
     &          ,DC%DH%Uncertainty)
        DC%DA%Back_Asym=BackSum
      else
        DC%DA%Back_Asym=0.
      endif
c       Combine the signal and background asymmetry to get the total asymmetry
      call MeasAsymPoint(DC%DA)
c      DC%DA%Asym=DC%DA%Signal_Asym
c     &                      -DC%DA%Back_Asym

      return
      end subroutine

cccccccc

cccc
      subroutine GetSignalAsym(DC,CubeAsym)
      implicit none
      integer Extent(3)
      Type(DataCubeAsymmetry),INTENT(INOUT) :: CubeAsym
      Type(DataCube), INTENT(INOUT) :: DC

c       Determine the extent using the cube's 'rotation point'
      call DetermineExtent(Extent,DC%DA%RotationPoint
     &               ,DC)
c       Run the main asymmetry loop
      call Main3DAsymLoop(DC%DA%RotationPoint
     &                  ,Extent,DC,CubeAsym)

      return
      end subroutine
ccccc




ccccc
      subroutine DetermineExtent(Extent,RotPt,DC)
      implicit none
      real, INTENT(IN) ::  RotPt(3)
      Type(DataCube), INTENT(IN) :: DC
      integer,INTENT(INOUT) :: Extent(3)

      integer i, CentIndx(3)
      do i=1,3
        CentIndx(i)=int(RotPt(i))
      enddo
c      print*, CentIndx

      Extent(1)=2*max(CentIndx(1),DC%DH%nPixels(0)-CentIndx(1)-1)+1
      Extent(2)=2*max(CentIndx(2),DC%DH%nPixels(1)-CentIndx(2)-1)+1
      Extent(3)=2*max(CentIndx(3),DC%DH%nChannels-CentIndx(3)-1)+1

c      print*, "RotPot", RotPt
c      print*, "Cent IIndx", CentIndx
c      print*, "Size Check", Extent
c      print*, "Cube Size", DC%DH%nPixels,DC%DH%nChannels

      return
      end subroutine
ccccccccc



ccccccccc
      subroutine Main3DAsymLoop(RotPt,Extent,DC,Asym)
      implicit none
      real,INTENT(IN) :: RotPt(3)
      integer, INTENT(INOUT) :: Extent(3)
      Type(DataCube), INTENT(IN) :: DC
      Type(DataCubeAsymmetry), INTENT(INOUT) :: Asym

      integer i,j,k
      integer CentIndx(3)
      integer NegIndx(3),CurrIndx(3),Size(3)
      integer nPairs
      real Ai,Bi

      real PosPt(3),NegPt(3)

c      print*, "Main  3D Aysm Loop",shape(DC%Flux)

      Size(1)=DC%DH%nPixels(0)
      Size(2)=DC%DH%nPixels(1)
      Size(3)=DC%DH%nChannels



      Asym%A=0.
      Asym%B=0.
      Asym%TotFlux=0.
      Asym%TotAbsDiff=0.
      nPairs=0

      do i=1,3
        CentIndx(i)=int(RotPt(i))
      enddo

c      print*, "main asym loop", CentIndx, Extent

c      i=37
c      j=37
c      k=0
c      print*, "Cent Pt Check", RotPt, CentIndx
c      print*, "Extent Check", Extent
c      print*, "Flux", DC%Flux(CentIndx(1),CentIndx(2),CentIndx(3))
c      Extent(1:3)=Extent(1:3)-50

c      print*, "loop low lims", CentIndx-Extent
c      print*, "loop high lims", CentIndx+Extent
c           Loop through all pixels
c      do i=0, Extent(1)/2   !One dimension can be cut in half to get all pairs
c      do i=60,60
      do i=CentIndx(1)-Extent(1)/2,CentIndx(1)+Extent(1)/2

c        print*, i,Asym%TotFlux,nPairs
c           Use i,j,k to identify the current index
        CurrIndx(1)=i
c        print*, i,Extent(1)-1
c        do j=0,Extent(2)-1
c        do j=47,47
        do j=CentIndx(2)-Extent(2)/2,CentIndx(2)+Extent(2)/2
            CurrIndx(2)=j
c            do k=0,Extent(3)-1
c            do k=49,49
            do k=CentIndx(3)-Extent(3)/2,CentIndx(3)+Extent(3)/2
                CurrIndx(3)=k
c          Get precise position for the pair of points
                call GetPtPairs(CentIndx,CurrIndx,RotPt,PosPt,NegPt)
c          Get the flux for each pixel by index
c                print*, "Pt Coords", PosPt,NegPt
                call GetFluxByPos(PosPt,Size,DC,Ai)
                call GetFluxByPos(NegPt,Size,DC,Bi)
c                print*, "Point Flux", PosPt,NegPt,Ai, Bi
c           Get the absolute difference between pixels and the sum of the fluxes
c                Asym%TotAbsDiff=Asym%TotAbsDiff+abs(Ai-Bi)
c                Asym%TotFlux=Asym%TotFlux+Ai+Bi       !The normal definition
c               THIS USING OUR NEW SQUARED DIFFERENCE METHOD
c                Asym%TotAbsDiff=Asym%TotAbsDiff+(Ai-Bi)**2.
c                Asym%TotFlux=Asym%TotFlux+(Ai+Bi)**2.       !The normal definition
c                call AsymPtPoint(Asym,Ai,Bi)
c                   Check that both points have flux in case of masking
c                if(Ai .gt. 0. .and. Bi .gt. 0.) then
                if(Ai .ne. 0. .or. Bi .ne. 0.) then
c                    if(Ai .lt. 0. .or. Bi .lt. 0.) then
c                        print*, "Point below zero flux", i,j,K
c     &                  ,PosPt,Ai,NegPt,Bi
c                    endif
c                    if(Ai .eq. 0. .or. Bi .eq. 0.) then
c                        print*, "hmmm", Ai,PosPt,Bi,NegPt
c                        print*, i,j,k
c            call GetFluxByInterpolation_Spec(PosPt,DC,Ai,Size)
c            call GetFluxByInterpolation_Spec(NegPt,DC,Bi,Size)
c                    endif

c                   Keep track of the number of unique pairs
                    call AsymPtPoint(Asym,Ai,Bi)
                    nPairs=nPairs +1
c                    print*, i,j,k,PosPt,NegPt,Ai,Bi
                endif

            enddo
        enddo
      enddo

c           THIS IS FOR TRYING SQUARES
c      Asym%TotAbsDiff=sqrt(Asym%TotAbsDiff)
c      Asym%TotFlux=sqrt(Asym%TotFlux)

c      Asym%Asym=Asym%TotAbsDiff/Asym%TotFlux
      call AsymCombPoint(Asym)
      Asym%nPairs=nPairs
c       Try this
c       print*, "Counted nPairs",nPairs
c      nPairs=Extent(1)*Extent(2)*Extent(3)/2
c      Asym%nPairs=nPairs
c      print*, "Asym Test",Asym%Asym
c     &          ,Asym%TotAbsDiff,Asym%TotFlux,nPairs
c     &          , "Squared Method"

      return
      end subroutine
ccccc


cccc
      subroutine GetPtPairs(CentIndx,CurrIndx,RotPt,PosPt,NegPt)
      implicit none
      integer,INTENT(IN) :: CentIndx(3),CurrIndx(3)
      real,INTENT(IN) :: RotPt(3)
      real,INTENT(INOUT) :: PosPt(3),NegPt(3)
      integer i,delta
    

      do i=1,3
        delta=CentIndx(i)-CurrIndx(i)
        PosPt(i)=RotPt(i)-delta
        NegPt(i)=RotPt(i)+delta
c        print*, "Delta check", i, delta, 2*CentIndx(i)-CurrIndx(i)
c     &              ,CurrIndx(i),PosPt(i), NegPt(i)
      enddo

      return
      end subroutine
cccccc


cc
c
      subroutine GetFluxByPos(Pt,Size,DC,Flux)
      implicit none
      Type(DataCube),INTENT(IN) :: DC
      integer, INTENT(IN) :: Size(3)
      real,INTENT(IN) :: Pt(3)
      real,INTENT(INOUT) :: Flux
      logical BoundCheck
c           First check if the indx is in the cube
      call GetBoundCheck(Pt,Size,BoundCheck)
      if(BoundCheck) then
c        Flux=DC%Flux(Indx(1),Indx(2),Indx(3))
        call GetFluxByInterpolation(Pt,DC,Flux,Size)
      else
        Flux=0.
      endif
      return
      end subroutine
cccccc



cccc
c
      subroutine GetBoundCheck(Pt,Size,BoundCheck)
      implicit none
      integer,INTENT(IN) :: Size(3)
      real,INTENT(IN) ::Pt(3)
      logical,INTENT(INOUT) :: BoundCheck
      integer i
      real,parameter :: TINY=1.e-20

      BoundCheck=.True.
      do i=1,3
c        print*, "Boundcheck", Pt(i), Size(i)
        if(Pt(i) .lt. 0. .or. Pt(i) .ge. Size(i)+TINY) then
            BoundCheck=.False.
c            print*, "False check", i, Pt(i),Size(i)
        endif
      enddo

      return
      end subroutine
ccccccc


ccccc
c
      subroutine GetFluxByInterpolation(Pt,DC,Flux,Size)
      implicit none
      integer,INTENT(IN) :: Size(3)
      real,INTENT(IN) :: Pt(3)
      Type(DataCube),INTENT(IN) :: DC
      real,INTENT(INOUT)::Flux

      integer i,j,k,ii,jj,kk,ll,CurrIndx(3)
      real CornerPts(8,4),InterpolatePt(4)
      logical BoundCheck

c      print*,"DC shape", shape(DC%Flux)
c       Get the indices of the lower corner of the box
      do i=1,3
        CurrIndx(i)=int(Pt(i))
      enddo

c       Set the interpolated point position
      InterpolatePt(1:3)=Pt(1:3)
c       Place the cube values that surround the point into an array
c               Note that the order has to be lower vel to upper vel
      do i=1,2
        ii=CurrIndx(3)+(i-1)
        do j=1,2
            jj=CurrIndx(2)+(j-1)
            do k=1,2
                kk=CurrIndx(1)+(k-1)
c                   Get the count for the corners
                ll=k+(j-1)*2+(i-1)*2*2
c               Set the corner point values
                CornerPts(ll,1)=real(kk)
                CornerPts(ll,2)=real(jj)
                CornerPts(ll,3)=real(ii)
c                print*, kk,jj,ii,ll
                call GetBoundCheck(CornerPts(ll,1:3),Size,BoundCheck)
c                print*, "Corner Pts", ll,CornerPts(ll,1:3)
c     &                  ,BoundCheck,Size
                if(BoundCheck) then
c                    CornerPts(ll,4)=DC%MaskedFlux(kk,jj,ii)
                    CornerPts(ll,4)=DC%Flux(kk,jj,ii)
c                    print*,"Corner checck", DC%Flux(kk,jj,ii)
                else
                    CornerPts(ll,4)=0.
                endif
c               If any corner points flux is zero, set the total flux to zero...
c                if (CornerPts(ll,4) .eq. 0.) then
c                    Flux=0.
c                    return
c                endif
c                print*, "Corner Pts", ll,CornerPts(ll,1:4)
            enddo
        enddo
      enddo
c           Use trilinear interpolation to get the flux at the specified point
      call TriLinearInterpolation(InterpolatePt,CornerPts) !/src/StandardMath/Interpolation.f
c      print*, "Interpolated flux", InterpolatePt
      Flux=InterpolatePt(4)
      kk=CurrIndx(1)
      jj=CurrIndx(2)
      ii=CurrIndx(3)
c      print*, CurrIndx,Pt,Flux
      if(DC%MaskedFlux(kk,jj,ii) .eq. 0.) then
        Flux=0.
      endif

      return
      end subroutine
ccccccc


ccccc
      subroutine EstimateNoiseSignal(BackSum,nPairs,RMS)
      use BasicRanNumGen
      implicit none
      integer, INTENT(IN) :: nPairs
      real, INTENT(IN) :: RMS
c      integer,INTENT(INOUT) :: idum
      real, INTENT(INOUT) :: BackSum
      integer i
      real T1,T2

c      BackSum=real(nPairs)*2.*RMS/sqrt(3.1415)
      BackSum=0.
      do i=1, nPairs
        T1=gasdev(idum)*RMS
        T2=gasdev(idum)*RMS
        BackSum=BackSum+(T1-T2)**2.
      enddo
      print*, "Background summation",BackSum,RMS
     &          , 2.*nPairs*RMS**2.
     &          , 4.*nPairs*RMS**2./3.1415921
      return
      end subroutine
cccccc


cccc
      subroutine MakeSymmetricMask(RotPt,MaskedDC,SymMaskedDC)
      implicit none
      real,INTENT(IN) :: RotPt(3)
      Type(DataCube),INTENT(IN) :: MaskedDC
      Type(DataCube),INTENT(INOUT) :: SymMaskedDC

      integer i,j,k
      integer CentIndx(3),Size(3)
      integer CurrIndx(3)
      integer PosIndx(3),NegIndx(3)
      integer Extent(3)

      real TempPt(3)
      logical bCheck1,bCheck2
      real F1,F2

      print*, "Making a mask symmetric about the rotation point"
     &          , RotPt
      print*, "Initial number of cells", sum(MaskedDC%Flux)
c       First set the Symmeterized mask flux to match the mask
      SymMaskedDC%Flux=MaskedDC%Flux
c           First determine the extent that we're going to need
      call  DetermineExtent(Extent,RotPt,MaskedDC)
c      print*, "Extent",Extent
c           Set the cube size for bound checking
      Size(1)=MaskedDC%DH%nPixels(0)
      Size(2)=MaskedDC%DH%nPixels(1)
      Size(3)=MaskedDC%DH%nChannels

      do i=1,3
        CentIndx(i)=int(RotPt(i))
      enddo

c           Loop through all pixels
      do i=0, Extent(1)-1
c           Use i,j,k to identify the current index
        CurrIndx(1)=i
c        print*, i,Extent(1)-1
        do j=0,Extent(2)-1
            CurrIndx(2)=j
            do k=0,Extent(3)-1
                CurrIndx(3)=k
                call GetPtIndices(CentIndx,CurrIndx,RotPt
     &                  ,PosIndx,NegIndx)
                TempPt=real(PosIndx)
                call GetBoundCheck(TempPt,Size,bCheck1)
                TempPt=real(NegIndx)
                call GetBoundCheck(TempPt,Size,bCheck2)
c               Check that both points are inbounds
                if(bCheck1 .and. bCheck2) then
c                    print*, "points are in bound",PosIndx,NegIndx
c                       Get the flux for each point
                    F1=SymMaskedDC%Flux(PosIndx(1)
     &                      ,PosIndx(2),PosIndx(3))
                    F2=SymMaskedDC%Flux(NegIndx(1)
     &                       ,NegIndx(2),NegIndx(3))
c                       Check if the fluxes are different
                    if(F1 .ne. F2) then
c                        print*, "Unequal flux",PosIndx,F1
c     &                              ,NegIndx,F2
c                       Set the two cells to 1 in the symmetric mask
                        SymMaskedDC%Flux(PosIndx(1)
     &                          ,PosIndx(2),PosIndx(3))=1.
                        SymMaskedDC%Flux(NegIndx(1)
     &                          ,NegIndx(2),NegIndx(3))=1.
                    endif
                endif
c           Check if the points are inbound
c                print*, i,j,k,PosIndx,NegIndx,RotPt
            enddo
        enddo
      enddo
      print*, "Final numbler of cells", sum(SymMaskedDC%Flux)

      return
      end subroutine
cccccc

cccccc
c           This subroutine gets indices matching points
      subroutine GetPtIndices(CentIndx,CurrIndx,RotPt
     &          ,PosIndx,NegIndx)
      implicit none
      integer,INTENT(IN) :: CentIndx(3),CurrIndx(3)
      real,INTENT(IN) :: RotPt(3)
      integer,INTENT(INOUT) :: PosIndx(3),NegIndx(3)
      integer i,delta


      do i=1,3
        delta=CentIndx(i)-CurrIndx(i)
        PosIndx(i)=int(RotPt(i)-delta)
        NegIndx(i)=int(RotPt(i)+delta)
      enddo

      return
      end subroutine
cccccc



ccccc
c
      subroutine GetFluxByInterpolation_Spec(Pt,DC,Flux,Size)
      implicit none
      integer,INTENT(IN) :: Size(3)
      real,INTENT(IN) :: Pt(3)
      Type(DataCube),INTENT(IN) :: DC
      real,INTENT(INOUT)::Flux

      integer i,j,k,ii,jj,kk,ll,CurrIndx(3)
      real CornerPts(8,4),InterpolatePt(4)
      logical BoundCheck

c      print*,"DC shape", shape(DC%Flux)
c       Get the indices of the lower corner of the box
      do i=1,3
        CurrIndx(i)=int(Pt(i))
      enddo

c       Set the interpolated point position
      InterpolatePt(1:3)=Pt(1:3)
c       Place the cube values that surround the point into an array
c               Note that the order has to be lower vel to upper vel
      do i=1,2
        ii=CurrIndx(3)+(i-1)
        do j=1,2
            jj=CurrIndx(2)+(j-1)
            do k=1,2
                kk=CurrIndx(1)+(k-1)
c                   Get the count for the corners
                ll=k+(j-1)*2+(i-1)*2*2
c               Set the corner point values
                CornerPts(ll,1)=real(kk)
                CornerPts(ll,2)=real(jj)
                CornerPts(ll,3)=real(ii)
c                print*, kk,jj,ii,ll
                call GetBoundCheck(CornerPts(ll,1:3),Size,BoundCheck)
c                print*, "Corner Pts", ll,CornerPts(ll,1:3)
c     &                  ,BoundCheck,Size
                if(BoundCheck) then
                    CornerPts(ll,4)=DC%MaskedFlux(kk,jj,ii)
c                    CornerPts(ll,4)=DC%Flux(kk,jj,ii)
c                    print*,"Corner checck", DC%Flux(kk,jj,ii)
                else
                    CornerPts(ll,4)=0.
                endif
c               If any corner points flux is zero, set the total flux to zero...
c                if (CornerPts(ll,4) .eq. 0.) then
c                    Flux=0.
c                    return
c                endif
                print*, "Corner Pts", ll,CornerPts(ll,1:4)
            enddo
        enddo
      enddo
c           Use trilinear interpolation to get the flux at the specified point
      call TriLinearInterpolation(InterpolatePt,CornerPts) !/src/StandardMath/Interpolation.f
      print*, "Interpolated flux", InterpolatePt
      Flux=InterpolatePt(4)

      return
      end subroutine
ccccccc

      end module
