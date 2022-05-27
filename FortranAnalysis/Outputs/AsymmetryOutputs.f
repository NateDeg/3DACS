cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines for outputing the results
c       from asymmetry measurements.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module AsymmetryOutputsMod
      use AsymmetryGlobals

      contains
ccccc
c
      subroutine OutputAsymmetry(SNPeak,SNAvg,SNInt)
      implicit none
      real,INTENT(IN) :: SNPeak,SNAvg,SNInt
      
      character(100) LineStr
      character(50) VarName
      character(20) ValStr
      integer WriteUnit, i

      print*, "Outputing asymmetry measurements to file"

      WriteUnit=10
      open(WriteUnit,file=trim(AsymmetryFileName),status='replace')
      write(WriteUnit,'(a)') "Asymmetry for cube:"
      write(WriteUnit,'(a)') trim(DataCubeName)

c           Write out the cube measurements
      write(WriteUnit,'(a)') ""
      VarName="Cube_RMS"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedDC%DH%Uncertainty,LineStr)
    
      VarName="SN"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,SNInt,LineStr)

      VarName="SN_Peak"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,SNPeak,LineStr)

      VarName="SN_Avg"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,SNAvg,LineStr)

c           Write out the 3D asymmetry measures
      write(WriteUnit,'(a)') ""
      VarName="Rot_Point_3D"
      LineStr=trim(VarName)//"=["
      do i=1,3
        write(ValStr,'(F8.3)') ObservedDC%DA%RotationPoint(i)
        LineStr=trim(LineStr)//trim(ValStr)
        if(i .lt. 3) then
            LineStr=trim(LineStr)//","
        else
            LineStr=trim(LineStr)//"]"
        endif
      enddo
      write(WriteUnit,'(a)') trim(LineStr)

      VarName="C3D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedDC%DA%Signal_Asym,LineStr)

      VarName="P3D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedDC%DA%TotAbsDiff,LineStr)

      VarName="Q3D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedDC%DA%TotFlux,LineStr)

      VarName="B3D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedDC%DA%Back_Asym,LineStr)

      VarName="A3D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedDC%DA%Asym,LineStr)



c           Write out the 2D asymmetry measures
      write(WriteUnit,'(a)') ""
      VarName="Rot_Point_2D"
      LineStr=trim(VarName)//"=["
      do i=1,2
        write(ValStr,'(F8.3)') ObservedMap%DA%RotationPoint(i)
        LineStr=trim(LineStr)//trim(ValStr)
        if(i .lt. 2) then
            LineStr=trim(LineStr)//","
        else
            LineStr=trim(LineStr)//"]"
        endif
      enddo
      write(WriteUnit,'(a)') trim(LineStr)

      VarName="Map_RMS"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedMap%DH%Uncertainty,LineStr)

      VarName="A2D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedMap%DA%Asym,LineStr)

      VarName="C2D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedMap%DA%Signal_Asym,LineStr)

      VarName="B2D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedMap%DA%Back_Asym,LineStr)

      VarName="P2D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedMap%DA%TotAbsDiff,LineStr)

      VarName="Q2D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedMap%DA%TotFlux,LineStr)



c           Write out the 1D asymmetry measures
      write(WriteUnit,'(a)') ""
      VarName="Rot_Point_1D"
      LineStr=trim(VarName)//"=["
      write(ValStr,'(F8.3)') ObservedProfile%DA%RotationPoint(3)
      LineStr=trim(LineStr)//trim(ValStr)
      LineStr=trim(LineStr)//"]"
      write(WriteUnit,'(a)') trim(LineStr)

      VarName="Profile_RMS"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedProfile%DH%Uncertainty,LineStr)

      VarName="A1D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedProfile%DA%Asym,LineStr)

      VarName="C1D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedProfile%DA%Signal_Asym,LineStr)

      VarName="B1D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedProfile%DA%Back_Asym,LineStr)

      VarName="P1D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedProfile%DA%TotAbsDiff,LineStr)

      VarName="Q1D"
      call WriteRealOutputStr(WriteUnit,VarName
     &          ,ObservedProfile%DA%TotFlux,LineStr)

      close(10)

      return
      end subroutine
ccccccc

cccc
c       This subroutine makes an output string from a real value
      subroutine WriteRealOutputStr(Unit,VarName,Val,Line)
      implicit none
      integer,INTENT(IN) :: Unit
      real,intent(IN) :: Val
      character(50),INTENT(IN):: VarName
      character(100),INTENT(INOUT):: Line
      character(20) ValStr

      write(ValStr,'(F15.7)') Val
      Line=trim(VarName)//"="//trim(ValStr)
      write(Unit,'(a)') trim(Line)
    
      return
      end subroutine
ccccc

c

      end module AsymmetryOutputsMod
