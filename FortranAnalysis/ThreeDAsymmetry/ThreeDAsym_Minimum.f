cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for finding the
c       point that minimizes the 3D asymmetry.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ThreeDAsym_MinimumMod

      use AsymmetryGlobals
      use DataCubeMod
      use ThreeDAsymCoreMod
      use DownhillSimplexMod

      implicit none

      contains


cccccc
c
c       This is the full asymmetry minimization routine
c
c       This routine assumes that you have set the 'WorkingDC' pointer
c           to the correct cube for all calculations.
      subroutine Calc_MinThreeDAsym()
      implicit none

      real,ALLOCATABLE:: AsymTests(:), Params(:,:)

      real rtol
      integer iter, nParams, mParams
      logical NoConvergenceFlag
      integer i

      Type(DataCubeAsymmetry) CubeAsym



      print*, "Calculating minimum Asymmetry"
c       First figure out how many free parameters we'll be using
      call DetermineNParams(nParams,WorkingDC)
      print*, "number of free parameters",nParams
      mParams=nParams+1
c       Now allocate the parameter vector
      ALLOCATE(Params(nParams+1,nParams))
      ALLOCATE(AsymTests(nParams+1))

c           First initialize the parameter vectors needed
c               for the downhill simplex function
      call InitializeParams(nParams,Params,AsymTests,1)

      rtol=5.e-2
      iter=0

c       Now we can call the downhill simplex routine
      call amoeba(Params,AsymTests,mParams,nParams
     &          ,rtol
     &          ,CubeAsymmetryFromParamVector
     &          ,iter,NoConvergenceFlag)
      print*, "After first pass"
      do i=1, mParams
        print*, "Params Set", i, Params(i,:),AsymTests(i)
      enddo

c           Reinitialize the parameter vectors around the best point
      call InitializeParams(nParams,Params,AsymTests,2)
c       Re-run ameoba with the new starting points
      call amoeba(Params,AsymTests,mParams,nParams
     &          ,rtol
     &          ,CubeAsymmetryFromParamVector
     &          ,iter,NoConvergenceFlag)
      do i=1, mParams
        print*, "Final Params Set", i, Params(i,:),AsymTests(i)
      enddo

c       Re-calculate everything with the best parameter
c           Set the rotation point to the best point
      call SetRotationPointFromParam(nParams,Params(1,:))
      call ThreeDAsym_CenterSet(WorkingDC)


      return
      end subroutine

cccccccc

ccccc
c
      subroutine DetermineNParams(nParams,DC)
      implicit none
      integer,INTENT(INOUT) :: nParams
      Type(DataCube),INTENT(IN) :: DC
    
      nParams=0
      if(DC%DH%nPixels(0) .gt. 1) then
        nParams=nParams+1
      endif
      if(DC%DH%nPixels(1) .gt. 1) then
        nParams=nParams+1
      endif
      if(DC%DH%nChannels .gt. 1) then
        nParams=nParams+1
      endif
      return
      end subroutine
cccccc

cccccc
      subroutine InitializeParams(nParams,Params,Asyms,PassCount)
      use BasicRanNumGen
      implicit none
      integer, INTENT(IN) :: nParams,PassCount
      real,INTENT(INOUT) :: Params(nParams+1,nParams)
      real,INTENT(INOUT) :: Asyms(nParams+1)
      integer i, j,count
      integer Size(3)

      real AvgFlux, AvgFluxFactor,PtFlux
      real Lims(2,nParams)
      real range(nParams),Start(nParams)
      logical BoundCheck


c       Start by getting the average flux as a limit
      if(nParams .eq. 1) then
        AvgFluxFactor=1.25
      else
        AvgFluxFactor=2.
      endif
      Size(1:2)=WorkingDC%DH%nPixels(0:1)
      Size(3)=WorkingDC%DH%nChannels

      AvgFlux=sum(abs(WorkingDC%Flux))/abs(Product(Size(1:3)))

      count=0
c       Set the parameter limits
      call SetLims(nParams,Lims)
      call SetRangeAndStartPoints(nParams,Lims,range,Start
     &                  ,Params(1,:),PassCount)


c       Now set the initial parameter parmeter vectors
      do i=PassCount,nParams+1
c        i=1
100     call SetParamVector(nParams,i,Params(i,:), count,Lims
     &              ,range,Start)
c        print*, "Ini Params val", i, count, Params(i,:)
        call SetRotationPointFromParam(nParams,Params(i,:))
c        print*, "Rotation Point", WorkingDC%DA%RotationPoint
c           Get the flux at this point
        call GetFluxByPos(WorkingDC%DA%RotationPoint,Size
     &          ,WorkingDC,PtFlux)
c        print*, "Point flux", i,Params(i,:)
c     &              ,PtFlux,count,nParams

c       Do a check to see if the count gets too high
        if(count .ge. 5000) then
            print*, "Unable to get a point with a high"
     &          ,"enough flux over the average"
            print*, "Avg flux and factor", AvgFlux,AvgFluxFactor
            stop
        endif
c           Make sure there is enough flux at this point
c               Note that this must happen after the count check
        if(PtFlux .lt. AvgFluxFactor*AvgFlux) goto 100
c           Once aa point is accepted, calculate the asymmetry at that point
        call CubeAsymmetryFromParamVector(nParams
     &              ,Params(i,:),Asyms(i))
        print*, "Initial parameter", Params(i,:),Asyms(i)

      enddo

      return
      end subroutine
ccccccc

cccccc
c       This routine sets the limits based on the number of parameters
c
      subroutine SetLims(nParams,Lims)
      implicit none
      integer,INTENT(IN):: nParams
      real,INTENT(INOUT):: Lims(2,nParams)
c       First set the parameter limits
      if(nParams .eq. 1) then
        Lims(1,1)=0.
        Lims(2,1)=real(WorkingDC%DH%nChannels)-1.
      elseif(nParams .eq. 2) then
        Lims(1,1:2)=0.
        Lims(2,1:2)=real(WorkingDC%DH%nPixels(0:1))-1.
      elseif(nParams .eq. 3) then
        Lims(1,1:3)=0.
        Lims(2,1:2)=real(WorkingDC%DH%nPixels(0:1))-1.
        Lims(2,3)=real(WorkingDC%DH%nChannels)-1.
      endif
      return
      end subroutine
cccccc

cccccc
c       This routine sets the range and start points
c           for selecting random parameters.  The precise values depend
c           on whether it is the 1st pass or second pass
c
      subroutine SetRangeAndStartPoints(nParams,Lims
     &                  ,range,Start,Param,PassCount)
      implicit none
      integer,INTENT(IN):: nParams,PassCount
      real,INTENT(IN):: Lims(2,nParams),Param(nParams)
      real,INTENT(INOUT) :: range(nParams),Start(nParams)
      integer i
      real Width
c       First set the parameter limits

      Width=5.
      if(PassCount .eq. 1) then
c           For the first pass use the full range of the parameter dimension
        do i=1,nParams
            range(i)=Lims(2,i)-Lims(1,i)
            Start(i)=Lims(1,i)
        enddo
      elseif(PassCount .eq. 2) then
        do i=1, nParams
            range(i)=2.*Width
            Start(i)=Param(i)-Width
        enddo
      endif
c      print*, 'Range', range
c      print*, "Start point", Start
        
      return
      end subroutine
cccccc

cccccccc
c           This subroutine sets a specific parameter vector
c
      subroutine SetParamVector(nParams,step,Param,count
     &          ,Lims,range,Start)
      implicit none
      integer,INTENT(IN):: nParams,step
      integer,INTENT(INOUT) :: count
      real,INTENT(INOUT) :: Param(nParams)
      real,INTENT(IN) :: Lims(2,nParams)
      real,INTENT(IN) ::range(nParams),Start(nParams)




      integer i
      real RandVal
c       If this is the first iteration of a point, use the cube
c           center

      if(count .eq. 0 .and. step .eq. 1) then
        do i=1,nParams
            Param(i)=Lims(2,i)/2.  !Note that Lims(2,i)==upper limit of the parameter
        enddo
      else
c           Otherwise, set the point to be somewhere in the cube
        do i=1,nParams
            call RANDOM_NUMBER(RandVal)!RANDOM_NUMBER is a gfortran built-in random number
c            Param(i)=(Lims(2,i)-Lims(1,i))*RandVal+Lims(1,i)
            Param(i)=range(i)*RandVal+Start(i)
        enddo
      endif
      count=count+1


      return
      end subroutine
cccccc


ccccc
c       This routine sets the cubes rotation point based on
c           a parameter vector
c
      subroutine SetRotationPointFromParam(nParams,Param)
      implicit none
      integer,INTENT(IN) :: nParams
      real,INTENT(IN) :: Param(nParams)

c       Initialize all rotation points to zero
      WorkingDC%DA%RotationPoint(1:3)=0.
      if(nParams .eq. 1) then
        WorkingDC%DA%RotationPoint(3)=Param(1)
      elseif(nParams .eq. 2) then
        WorkingDC%DA%RotationPoint(1:2)=Param(1:2)
      elseif(nParams .eq. 3) then
        WorkingDC%DA%RotationPoint(1:3)=Param(1:3)
      endif

      return
      end subroutine
ccccc

cccccccc
c
      subroutine CubeAsymmetryFromParamVector(nParams,ParamVector,Asym)
      use AsymmetryIniMod
c           This routine calculates 3D asymmetry from a given parameter
c               vector
      implicit none
      integer,INTENT(IN) ::  nParams
      real,INTENT(IN) :: ParamVector(nParams)
      real,INTENT(INOUT) :: Asym

      Type(DataCubeAsymmetry) CubeAsym
      logical BoundCheck
      integer CubeSize(3)

c      print*, "Initial cube asym from param"
c      print*,"Paramvector",ParamVector
      CubeSize(1:2)=WorkingDC%DH%nPixels(0:1)
      CubeSize(3)=WorkingDC%DH%nChannels
c       Now set the rotation point based on the parameter vector
      call SetRotationPointFromParam(nParams,ParamVector)
c      print*, "Getting Asym from param vector"
c      print*, "Rotation point is ", WorkingDC%DA%RotationPoint
c     &          ,nParams
c       Check that the rotation point
      call GetBoundCheck(WorkingDC%DA%RotationPoint,CubeSize,BoundCheck)
c      print*, "Bound Check", CubeSize,BoundCheck
      if(BoundCheck) then
c       With the rotation point set, we can get the signal asymmetry
c           If this is 3D, we can make the mask symmetric
      if(WorkingDC%DH%SymmetricMaskSwitch.eq. 1) then
        call MakeSymmetricMask(WorkingDC%DA%RotationPoint
     &          ,DataCubeMask,SymmetricMask)
c           Remask the data
        call MaskCube(SymmetricMask,WorkingDC)
      endif

        call GetSignalAsym(WorkingDC,CubeAsym)
        Asym=CubeAsym%Asym
c        Asym=3.
      else
c           If the point is outside the cube, set the asymmetry to be too large
        Asym=2.
      endif

      return
      end subroutine
ccccccccc





      end module
