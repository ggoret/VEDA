c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program rota
      implicit none
      integer i
      double precision dtor,mx,pi,r,r1,r2,r3,rtod,twopi,u1,u2,u3,xhi1,
     & xhi2,xhi3
      dimension mx(9),r(9),r1(9),r2(9),r3(9),u1(3),u2(3),u3(3)
      common/angkte/ pi,twopi,dtor,rtod
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
       do i=1,3
      u1(i)=0.
      u2(i)=0.
      u3(i)=0.
       enddo
      u1(1)=1.
      u2(2)=1.
      u3(3)=1.
      call redmx(r)
      call prtmx(r)
      read(5,*) xhi1,xhi2,xhi3
      call rmxp(xhi1,u1,r1)
      call rmxp(xhi2,u2,r2)
      call rmxp(xhi3,u3,r3)
      write(6,*) ' '
      call pro2mx(r,r3,r1)
      call pro2mx(mx,r2,r)
      call prtmx(mx)
      write(6,*) ' '
      call pro2mx(r,r1,r3)
      call pro2mx(mx,r2,r)
      call prtmx(mx)
      stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine redmx(r)
c this subroutine reads a matrix
      implicit none
      integer i,j
      double precision r
      dimension r(3,3)
      read(5,*) ((r(i,j),j=1,3),i=1,3)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prtmx(r)
c this subroutine prints a matrix
      implicit none
      integer i,j
      double precision r
      dimension r(3,3)
      write(6,'(3(/3f10.5))') ((r(i,j),j=1,3),i=1,3)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pro2mx(ro,r2,r1)
c this subroutine calculates the product of two matrices: ro=r2*r1
      implicit none
      integer i,j,k
      double precision add,r1,r2,ro
      dimension r1(3,3),r2(3,3),ro(3,3)
      do 30 i=1,3
      do 20 j=1,3
      add=0.d0
      do 10 k=1,3
      add=add+r2(i,k)*r1(k,j)
10    continue
      ro(i,j)=add
20    continue
30    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmx2p(rotx,xhi,vn)
      implicit none
      integer i
      double precision ang,cang,dtor,fuzz,norm,pi,rotx,rtod,sang,scale,
     & twopi,vecx,vn,xhi
      dimension rotx(3,3),vecx(3),vn(3)
      common/angkte/ pi,twopi,dtor,rtod
      fuzz=0.01
      cang=min(max((rotx(1,1)+rotx(2,2)+rotx(3,3)-1)/2,-1.d0),1.d0)
      ang=acos(cang)
      sang=sin(ang)
         if(sang.gt.fuzz) then
      scale=2*sang
      vecx(1)=(rotx(3,2)-rotx(2,3))/scale
      vecx(2)=(rotx(1,3)-rotx(3,1))/scale
      vecx(3)=(rotx(2,1)-rotx(1,2))/scale
         else
            if(cang.gt.0.) then
      vecx(1)=0.d0
      vecx(2)=0.d0
      vecx(3)=1.d0
            else
       do i=1,3
      vecx(i)=sqrt(max((rotx(i,i)+1)/2,0.d0))
       enddo
               if(vecx(3).gt.fuzz) then
      vecx(1)=sign(vecx(1),rotx(1,3))
      vecx(2)=sign(vecx(2),rotx(2,3))
               else if(vecx(2).gt.fuzz) then
      vecx(1)=sign(vecx(1),rotx(1,2))
               endif
            endif
         endif
      xhi=rtod*ang
       do i=1,3
      vn(i)=vecx(i)
       enddo
      norm=vecx(1)**2+vecx(2)**2+vecx(3)**2
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmxp(xhi,vn,rotx)
      implicit none
      integer eijk,i,j,k
      double precision ang,cang,dtor,norm,pi,rotx,rtod,sang,twopi,vn,xhi
      dimension eijk(3,3,3),rotx(3,3),vn(3)
      common/angkte/ pi,twopi,dtor,rtod
      data eijk/
     & 0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
      norm=sqrt(vn(1)**2+vn(2)**2+vn(3)**2)
       do i=1,3
      vn(i)=vn(i)/norm
       enddo
      ang=dtor*xhi
      cang=cos(ang)
      sang=sin(ang)
      do 30 j=1,3
      do 20 i=1,3
      rotx(i,j)=(1-cang)*vn(i)*vn(j)
      if(i.eq.j) rotx(i,j)=rotx(i,j)+cang
      do 10 k=1,3
         if(eijk(i,k,j).gt.0) then
      rotx(i,j)=rotx(i,j)+sang*vn(k)
         else if(eijk(i,k,j).lt.0) then
      rotx(i,j)=rotx(i,j)-sang*vn(k)
         endif
10    continue
20    continue
30    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trans(r)
c this subroutine transposes a matrix
      implicit none
      integer i,j
      double precision r,r1
      dimension r(3,3),r1(3,3)
       do i=1,3
        do j=1,3
      r1(i,j)=r(j,i)
        enddo
       enddo
       do i=1,3
        do j=1,3
      r(i,j)=r1(i,j)
        enddo
       enddo
      return
      end
