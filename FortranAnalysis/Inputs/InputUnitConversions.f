cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading in a
c       datacube.  It also contains routines for getting
c       a text version of a data cube object header.
c       Since the beam parameters are usually contained in the
c       datacube headers, this file also gets those in combined
c       input files
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module InputUnitConversionsMod

      use DataCubeMod
      use BeamMod
      use CommonConsts

      implicit none

      contains



ccccc
c     This routine converts angles to radians
      subroutine GeneralAngularConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=degrees to radians
        A=A*Pi/180.
      endif

      return
      end subroutine
ccccccccc

cccccccc
c       This routine converts distances (angular) to arcseconds
      subroutine GeneralDistanceConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=degrees to arcseconds
        A=A*3600.
      endif
      return
      end subroutine
cccccccc

cccccc
c       This routine does velocity conversions
      subroutine GeneralVelocityConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

c      print*, "Vel Conversion", A, Switch
      if(Switch .eq. 0) then        !0=m/s to km/s
        A=A/1000.
      endif
      return
      end subroutine
ccccccc
ccccccc
c     This routine does brightness conversion
      subroutine GeneralBrightnessConversion(Switch,A,ChannelSize
     &              ,BeamArea,PixelSize)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A
      real, INTENT(IN) :: ChannelSize,BeamArea,PixelSize

c       The surface brightness needs to be converted to Jy pixel^2
      if(Switch .eq. 0) then        !0=Jy km/s arcsec^-2 to Jy arcsec^-2
        A=A/abs(ChannelSize)
        A=A*abs(PixelSize)**2.           !Then Jy arcsec^2 to Jy pixel
      elseif(Switch .eq. 1) then    !1=mJy/beam to Jy/Beam
        A=A/1000.
      elseif(Switch .eq. 2) then    !2=Jy arcsec^-2 to Jy/Beam
        A=A*BeamArea
      endif
      return
      end subroutine
cccccc




ccccc
c
c       This routine converts the tilted ring
c           radial values (radius and widths) to pixel values
c       It assumes that the radius is in arcseconds
      subroutine Radii_PixelConversion(R
     &                  ,PixelSize)
      implicit none
      real, INTENT(IN) :: PixelSize
      real, INTENT(INOUT) :: R
c
      R=R/abs(PixelSize)
      return
      end subroutine
ccccccc



ccccc
c
c       This routine converts the tilted ring
c           center coordinates to a pixel value
c
c       For this routine to work, StartVal must be the
c       pixel location (0) and have units of arcseconds
      subroutine PixelPositionConversion(Switch,A
     &                  ,PixelSize,StartVal)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(IN) :: StartVal,PixelSize
      real, INTENT(INOUT) :: A
c
      if(Switch .eq. 2) then        !2==pixels already
        A=A
      elseif(Switch .eq. 1) then    !1==arcseconds to pixels
        A=(A-StartVal)/PixelSize
      elseif(Switch .eq. 0) then    !0==degrees to pixels
        A=A*3600.
        A=(A-StartVal)/PixelSize
      endif

      return
      end subroutine
ccccccc



cccccccc
c
      subroutine DataCubeUnitConversions(DC,Beam)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      Type(Beam2D), INTENT(INOUT) :: Beam
      real z
      real TempFreq,V1,V2
      integer i


c      print*, "Converting Units"

c       Do the velocity conversion first --- goal is to have the units in km/s
c      print*, "Data cube vel units ",trim(DC%DH%Units(2))
c      print*, "Data cube vel type ",trim(DC%DH%AxisType(2))
      if(trim(DC%DH%AxisType(2)) .eq. 'VELO-LSR'
     &       .or. trim(DC%DH%AxisType(2))
     &        .eq. 'VELOHEL'
     &    .or. trim(DC%DH%AxisType(2)) .eq. 'VOPT') then
c        print*, "DC Vel Type trigger"
        if(trim(DC%DH%Units(2)) .eq. 'm/s'
     &       .or. trim(DC%DH%Units(2)) .eq. 'm s-1') then
c            print*, "DC Vel Unit trigger"
            DC%DH%ChannelSize=DC%DH%ChannelSize/1000.
            DC%DH%ChannelCent=DC%DH%ChannelCent/1000.
            DC%DH%RefVal(2)=DC%DH%RefVal(2)/1000.
        endif
      endif


c       Do the angule conversion next -- initially want units of arcseconds
      if(trim(DC%DH%Units(0)) .eq. 'DEGREE') then
        Beam%BeamMajorAxis=Beam%BeamMajorAxis*3600.
        Beam%BeamMinorAxis=Beam%BeamMinorAxis*3600.
        DC%DH%PixelSize=DC%DH%PixelSize*3600.
        Beam%PixelSize=DC%DH%PixelSize
        DC%DH%RefVal(0:1)=DC%DH%RefVal(0:1)*3600.
      endif
c       Finally, put the beam units into pixels as the TR portion of the
c           code will work in pixel space
      Beam%BeamMajorAxis=Beam%BeamMajorAxis
     &                  /abs(DC%DH%PixelSize(0))
      Beam%BeamMinorAxis=Beam%BeamMinorAxis
     &                  /abs(DC%DH%PixelSize(0))



      return
      end subroutine
cccccccc



ccccccc
c     This routine does brightness conversion for a data cube to
c       get final units of Jy/pixel
      subroutine DCBrightnessConversion(DC,Beam)
      implicit none
      Type(Beam2D),INTENT(IN) :: Beam
      Type(DataCube),INTENT(INOUT) :: DC

      real BeamArea

c       The cell brighness needs to be in Jy/pixel
c      print*, "DataCube Units ", DC%DH%FUnit
      if(DC%DH%FUnit .eq. 'Jy/beam') then
c       The major and minor beam axis should already be converted to pixels
c           so get the area in pixels and divide it out from the flux
        BeamArea=2.*Pi
     &          *Beam%BeamMajorAxis*Beam%BeamMinorAxis
        DC%Flux=DC%Flux/BeamArea
      endif


      return
      end subroutine
ccccccc

cccccc
c
      subroutine RedshiftCalc(z,RestFreq,Freq)
      implicit none
      real,INTENT(IN) :: RestFreq,Freq
      real,INTENT(OUT) :: z
      real numer,denom

c       Using the optical definition
      numer=RestFreq-Freq
      denom=Freq
      z=numer/denom
c      print*, numer,denom,z

      return
      end subroutine
cccccccc




      end module
