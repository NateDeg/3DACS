
      module IntegrationMod
      use InterpolateMod
      contains

cccc
      subroutine TabulatedTrapzIntegral(Integral,n,X,Y,XEdge)
c       This routine calculates the integral of some profile with
c           regular X intervals using a trapezoid rul
      implicit none
      integer,INTENT(IN) :: n
      real,INTENT(IN) :: X(n),Y(n),XEdge(2)
      real,INTENT(INOUT) :: Integral

      integer i,j, IDLow,IDHigh
      real dX
      integer EdgeIDs(2),nBody,BodyLims(2)
      real BodyTerm,EdgeTerm
      real P1(2),P2(2),YEdge

c       print*, "Tab Int", XEdge,n,X(n)
c       First figure out dX
      dX=X(2)-X(1)
c       Figure out the bin ids of the two edges
      do i=1,2
        EdgeIDs(i)=(XEdge(i)-X(1))/dX+1
c           Check to make sure the edges are inside the profile
        if(EdgeIDs(i) .le. 0 .or. EdgeIDs(i) .gt. n) then
            print*, "Edge is not included in the integral"
            Integral=0.
            return
        endif
c        print*, i, EdgeIDs(i), X(EdgeIDs(i)),X(EdgeIDs(i)+1)
c     &          ,XEdge(i)
      enddo
c       Do the body integration
      nBody=EdgeIDs(2)-EdgeIDs(1)-1
      BodyLims(1)=EdgeIDs(1)+1
      BodyLims(2)=EdgeIDs(2)
c      print*, "Tab Check", EdgeIDs
c     &      ,BodyLims,n,nBody
c      print*, "number of points in body integration", nBody
      call TrapzBodyIntegration(BodyTerm,nBody+1
     &          ,X(BodyLims(1):BodyLims(2))
     &          ,Y(BodyLims(1):BodyLims(2)))

c      print*, "Body integration", BodyTerm
c       Now do the 2 edge integrations
    
      if (EdgeIDs(1) .lt. 0) then
        EdgeIDs(1)=0
      endif
      if (EdgeIDs(2) .eq. n) then
        EdgeIDs(2)=n-1
      endif
    
      call TrapzEdgeIntegration(EdgeTerm,nBody+3
     &                          ,X(EdgeIDs(1):EdgeIDs(2)+1)
     &                          ,Y(EdgeIDs(1):EdgeIDs(2)+1)
     &                          ,XEdge)

c      print*, "Edge Integration", EdgeTerm

      Integral=BodyTerm+EdgeTerm
c      Integral=BodyTerm
c      print*, "Table Integration", BodyTerm, EdgeTerm
      return
      end subroutine
ccccc


ccccc
      subroutine TrapzBodyIntegration(BodyInt,n,X,Y)
c           This routine does the trapezoid integration for the body (full bins)
c           of the profile
      implicit none
      integer,INTENT(IN):: n
      real,INTENT(IN) :: X(n),Y(n)
      real,INTENT(INOUT) :: BodyInt
      integer i
      real Area

c      print*, "Body Edge Check", X(1),X(n)
      BodyInt=0.
      do i=1,n-1
        call TrapRule(X(i:i+1),Y(i:i+1),Area)
        BodyInt=BodyInt+Area
c        print*, "Body Int", i, n, X(i:i+1)
      enddo
c      print*, "Final Body Int", BodyInt
      return
      end subroutine
ccccccc

cccccc
      subroutine TrapzEdgeIntegration(EdgeInt,n,X,Y,XEdges)
      implicit none
      integer,INTENT(IN):: n
      real,INTENT(IN) :: X(n),Y(n),XEdges(2)
      real,INTENT(INOUT) :: EdgeInt
      integer i,j
      real P1(2), P2(2), YEdge
      real XPair(2), YPair(2), PairInt

      EdgeInt=0.
c      print*, "Edge Check", X(1), X(n)
      do i=1, 2
c           Set up a pair of points for interpolation
        j=i
        if(i .eq. 2) then
            j=n-1
        endif
        P1(1)=X(j)
        P1(2)=Y(j)
        P2(1)=X(j+1)
        P2(2)=Y(j+1)
c           Get the interpolated points
        call SimpleInterpolateY(P1,P2,XEdges(i),YEdge)
c           Set up 2 element arrays for trapezoid integration
        if( i.eq. 1) then
            XPair(1)=XEdges(i)
            YPair(1)=YEdge
            XPair(2)=P2(1)
            YPair(2)=P2(2)
        elseif(i .eq. 2) then
            XPair(2)=XEdges(i)
            YPair(2)=YEdge
            XPair(1)=P1(1)
            YPair(1)=P1(2)
        endif
c       Run the trapezoid rule on this smaller section
c        print*, "Edge Pairs", XPair, YPair
c     &          , XPair(2)-XPair(1),X(2)-X(1)
c     &          , (XPair(2)-XPair(1))/(X(2)-X(1))
        call TrapRule(XPair,YPair, PairInt)
        EdgeInt=EdgeInt+PairInt
      enddo

      return
      end subroutine
ccccccc

cccc
      subroutine TrapRule(X,Y,Area)
c           The basic trapezoid rule
      implicit none
      real,INTENT(IN) :: X(2), Y(2)
      real,INTENT(INOUT) :: Area

      real dX,dY,YMin
      dX=X(2)-X(1)
      dY=abs(Y(2)-Y(1))
      YMin=minval(Y)
      Area=dX*(YMin+0.5*dY)
      return
      end subroutine
ccccccc


      end module

