cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading in the
c       controlling inputs for the 3D asymmetry code.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module AsymmetryInputMod

      use AsymmetryGlobals
      use DataCubeInputMod

      implicit none

      contains



ccccc
c     This routine converts angles to radians
      subroutine Asym3DInput()
      implicit none
      character(200) infile
      character(1) junk
      character(10) ErrMsg
      integer ReadUnit

      print*, "Get Asymmetry inputs"
c       First get the runtime input file
      call getarg(1,infile)
c       Check that the file exists
      ErrMsg="Input"
      call FileExistsCheck(infile,ErrMsg)

      ReadUnit=10
c       If the file exists, open it and get the inputs
      open(ReadUnit, file=trim(infile),status='old')
c       Read in the name of the data cube to analyze
      read(ReadUnit,*) junk
      ErrMsg="Data cube"
      call FileNameInput(DataCubeName,ReadUnit,ErrMsg)
c       Read in the name of the data cube mask to use
      read(ReadUnit,*) junk
      ErrMsg="Mask cube"
      call FileNameInput(MaskName,ReadUnit,ErrMsg)
c       Now get the type of mask file we'll be using
      read(ReadUnit,*) junk
      read(ReadUnit,*) MaskFileTypeSwitch
      if(MaskFileTypeSwitch .eq. 2) then
        Backspace(ReadUnit)
        read(ReadUnit,*) MaskFileTypeSwitch, MaskFluxFrac
        print*, "Mask Flux fraction", MaskFluxFrac
        if(MaskFluxFrac .le. 0. .or. MaskFluxFrac .ge. 1.) then
            print*, "The fraction of the flux to be used "
     &              ,"in the mask must be between 0 and 1"
            stop
        endif
      endif
c       Get the name of the output file
      read(ReadUnit,*) junk
      read(ReadUnit,'(a)')AsymmetryFileName

c       Read in the center indices to use for the Asymmetry3D code
      read(ReadUnit,*) junk
      read(ReadUnit,*) ObservedDC%DA%RotationPoint(1:3)
c       And now for the background switch
      read(ReadUnit,*) junk
      read(ReadUnit,*) BackgroundSwitch
c       And now for the Method switch
      read(ReadUnit,*) junk
      read(ReadUnit,*) AsymMethodSwitch
      if(AsymMethodSwitch .lt. 0
     &          .or. AsymMethodSwitch .ge. 2) then
        print*, "Select Valid Asymmetry Method"
      endif
            
c       Finally read in the random seed
      read(ReadUnit,*) junk
      read(ReadUnit,*) idum
c       Close the input file
      close(ReadUnit)
      print*, "Done reading main input file"
      
c       Once everything from the main input is finished, read in the data cube
c           and the mask from their respective files
c           Read in the mask
      call ReadFullDataCube(DataCubeMask
     &                      ,Beam,MaskName)
      print*, "Done reading in cube"
c       Initialize the symmeterized mask
      SymmetricMask%DH=DataCubeMask%DH
      call AllocateDataCube(SymmetricMask)
c           And the actual cube
      call ReadFullDataCube(ObservedDC
     &                      ,Beam,DataCubeName)
c       In order to get a S/N measure, we will need the beam to be allocated
c           to get the beam area measurement
      call Allocate_Beam2D(Beam
     &            ,ObservedDC%DH%nPixels(0:1))

      return
      end subroutine
ccccccccc


ccccc
c       This routine reads in a filename from an input file and checks if it exists
      subroutine FileNameInput(filenameVar,unit,ErrMsg)
      implicit none
      integer, INTENT(IN) :: unit
      character(*), INTENT(INOUT) :: filenameVar
      character(*), INTENT(IN) :: ErrMsg

      read(unit,'(a)') filenameVar
      call FileExistsCheck(trim(filenameVar),ErrMsg)

      return
      end subroutine
ccccc

ccccc
c           This routine checks if a specified file exists
c               and stops the code if it's missing
      subroutine FileExistsCheck(filename,errMsg)
      implicit none
      character(*),INTENT(IN):: filename
      character(*),INTENT(IN):: errMsg
      logical file_exists


      INQUIRE(FILE=filename, EXIST=file_exists)
      if(file_exists .eqv. .False.) then
        print*, trim(errMsg)," file supplied named: ", trim(filename)
        print*, "does not exist.  The code will stop here."
        stop
      endif
      return
      end subroutine
ccccccc
      end module
