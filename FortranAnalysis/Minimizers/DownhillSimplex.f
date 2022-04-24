cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains modified versions of the amoeba and
c       amotry routines from Press et al.  They implement the
c       Nelder and Mead method of minimizing a multi-dimensional function
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DownhillSimplexMod



      implicit none

      contains

ccccccc
      subroutine amoeba(p,y,mp,np,ftol,funk,iter
     &          ,NoConvergenceFlag)
      implicit none
      integer, INTENT(IN) :: mp,np
      integer,INTENT(INOUT) ::  iter
      real,INTENT(IN) :: ftol
      logical,intent(INOUT):: NoConvergenceFlag
      real,INTENT(INOUT) :: p(mp,np),y(mp)

      interface
        subroutine funk(nump, paramvector,val)
            integer,intent(IN) :: numP
            real,INTENT(IN) :: paramvector(nump)
            real,INTENT(INOUT) :: val
        end subroutine
      end interface

      integer NMAX,ITMAX

      real rtol,TINY
      parameter(ITMAX=5000,TINY=1.e-10)

      integer i,ihi,ilo,inhi,j,m,n
      real sum,swap,ysave,ytry
      real  psum(np)

c
c      ALLOCATE(psum(ndim))

      NoConvergenceFlag=.False.
c      print*, "Starting downhill simplex",PID,y
      iter=0
 1    do n=1,np
        sum=0.
        do m=1,np+1
            sum=sum+p(m,n)
        enddo
        psum(n)=sum
      enddo

c      print*, "amoeba fnk call", psum
c      call funk(np,psum,ytry)

2     ilo=1
      if(y(1) .gt. y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do i=1,np+1
        if(y(i) .le. y(ilo)) ilo=i
        if(y(i) .gt. y(ihi)) then
            inhi=ihi
            ihi=i
        elseif(y(i) .gt. y(inhi)) then
            if(i .ne. ihi) inhi=i
        endif
      enddo

      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
      print*, "Current tolerance", iter,rtol,y(ihi),y(ilo)
      if(rtol .lt. ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do n=1,np
            swap=p(1,n)
            p(1,n)=p(ilo,n)
            p(ilo,n)=swap
        enddo
        return
      endif
c
      if(iter .ge. ITMAX) then
        NoConvergenceFlag=.True.
        return
      endif

c      print*, "amoeba", psum
      iter=iter+2
      call amotry(ytry, p,y,psum,mp,np,funk,ihi,-1.0)
      if(ytry .le. y(ilo)) then
        call amotry(ytry, p,y,psum,mp,np,funk,ihi,2.0)
      elseif(ytry .ge. y(inhi)) then
        ysave=y(ihi)
        call amotry(ytry, p,y,psum,mp,np,funk,ihi,0.5)
        if(ytry .ge. ysave) then
            do i=1,np+1
                if(i .ne. ilo) then
                    do j=1,np
                        psum(j)=0.5*(p(i,j)+p(ilo,j))
                        p(i,j)=psum(j)
                    enddo
c                    print*, "amoeba fnk call", psum
                    call funk(np,psum,ytry)
                endif
            enddo
            iter=iter+np
            goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2

      return
      end subroutine
ccccccccccc


      subroutine amotry(ynew,p,y,psum,mp,np,funk,ihi,fac)
      implicit none
      integer,INTENT(IN) :: mp,np
      integer,INTENT(IN) :: ihi
      real,INTENT(IN) :: fac
      real, INTENT(INOUT) :: ynew,p(mp,np),y(mp),psum(np)
      interface
        subroutine funk(nump, paramvector,val)
            integer,intent(IN) :: numP
            real,INTENT(IN) :: paramvector(nump)
            real,INTENT(INOUT) :: val
        end subroutine
      end interface

      integer j
      real fac1,fac2,ytry,ptry(np)
c
c      print*, "in amotry",shape(p),mp,np,ihi
c
      fac1=(1.-fac)/real(np)
      fac2=fac1-fac
      do j=1,np
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
      enddo
c      print*, p(ihi,:)
c      print*, ptry
c      print*, shape(ptry)

      call funk(np,ptry,ytry)
      if(ytry .lt. y(ihi)) then
        y(ihi)=ytry
        do j=1,np
            psum(j)=psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j)=ptry(j)
        enddo
      endif
      ynew=ytry


      return
      end subroutine

      end module
