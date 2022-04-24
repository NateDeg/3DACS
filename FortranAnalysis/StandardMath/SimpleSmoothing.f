
      module SimpleSmoothMod
      use CommonConsts
      contains



ccccc
      subroutine Simple1DSmooth(nPoints,nSmooth,Y,YSmooth)
      implicit none
      integer, INTENT(IN) :: nPoints,nSmooth
      real,INTENT(IN) :: Y(0:nPoints-1)
      real,INTENT(INOUT) :: YSmooth(0:nPoints-1)

      real Kernel(nSmooth),Val
      integer i,j,k

c       First calculate a Gaussian kernel based on the number of smoothing points
      call Calculate1DGaussianKernel(nSmooth,Kernel)
c       Now initialize the smoothed values to zero
      YSmooth(0:nPoints-1)=0.
c       The next bit is to loop through the full profile
      do i=0,nPoints-1
c       For each point, calculate it's contribution to all smoothed points
        do j=1,nSmooth
c           First get the correct index
            k=i-nSmooth/2+j-1
            if(k .ge. 0 .and. k .le. nPoints-1) then
                Val=Y(i)
            else
                Val=0.
            endif
            if (k .ge. 0 .and. k .le. nPoints-1) then
                YSmooth(k)=Val*Kernel(j)+YSmooth(k)
            endif
        enddo
      enddo

      return
      end subroutine

cccc


ccccc
      subroutine Calculate1DGaussianKernel(n,K)
c       This subroutine calculates a simple 1D gaussian kernel
      implicit none
      integer, INTENT(IN) :: n
      real,INTENT(INOUT) :: K(n)

      real Mean,Width,Sigma
      integer i

c       Get the mean, the 95% limits (based on the size/2), and sigma
      Mean=(n+1)/2
      Width=n/2.
      Sigma=Width/2.

c       Calculate the Guassian kernel
      do i=1,n
        K(i)=exp(-(real(i)-Mean)**2./(2.*Sigma**2.))
     &          /(2.*Pi*Sigma**2.)**0.5
      enddo
c       Normalize the kernel to 1
      K=K/sum(K)
      return
      end subroutine
ccccccc



      end module

