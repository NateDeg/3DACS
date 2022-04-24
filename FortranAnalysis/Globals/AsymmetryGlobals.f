

      module AsymmetryGlobals
c               Use the various object definition routines
      use DataCubeMod
      use BeamMod
      use ProfileMod


      implicit none

      character(500) DataCubeName, MaskName
      character(500) AsymmetryFileName

      integer MaskFileTypeSwitch, BackgroundSwitch
      integer CenterSwitch
      integer AsymMethodSwitch
      real MaskFluxFrac

      Type(DataCube),target:: ObservedDC, DataCubeMask
      Type(DataCube),target:: SymmetricMask
      Type(DataCube),target:: ObservedMap, ObservedProfile
      Type(Beam2D) Beam

      Type(DataCube),pointer :: WorkingDC => null()


      integer idum


      integer nSmooth
      real PeakFrac

      end module AsymmetryGlobals
