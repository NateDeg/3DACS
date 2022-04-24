
      module BasicDerivativeMod
      contains

ccccc
      subroutine CalculateSimpleArrayDerivative(nPoints
     &              ,dX,Y,YPrime)
      implicit none
      integer, INTENT(IN) :: nPoints
      real, INTENT(IN) :: dX, Y(0:nPoints-1)
      real, INTENT(INOUT) :: YPrime(0:nPoints-1)

      integer i

c      print*, "get simple derivs", nPoints
c       Get the majority of the slopes
      do i=1, nPoints-2
        YPrime(i)=(Y(i+1)-Y(i-1))/(2.*dX)
      enddo

c       Get the first and last slope points
      i=0
      YPrime(i)=(Y(i+1)-Y(i))/(dX)
      i=nPoints-1
      YPrime(i)=(Y(i)-Y(i-1))/(dX)

      return
      end subroutine
cccccc


      end module

