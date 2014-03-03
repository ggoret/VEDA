c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program eul2mtx
      implicit none
      integer i,length
      double precision alf,bet,dtor,fuzz,gam,pi,r,rtod,twopi,u
      character card*80
      dimension r(9),u(3)
      common/angkte/ pi,twopi,dtor,rtod
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      fuzz=0.1
      i=0
10    read(5,'(a)',end=20) card
      call compact(card)
      if(length(card).eq.0) go to 20
      read(card,*) alf,bet,gam,u
      call rmxe(alf,bet,gam,r)
      i=i+1
      write(6,'(9f10.5,3f10.3,a,i3)') r,u,'  #',i
      go to 10
20    stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmxe(alpha,beta,gamma,rotx)
c this subroutine writes the rotation matrix, given the eulerian parameters
      implicit none
      double precision alf,alpha,bet,beta,cosa,cosb,cosg,dtor,gam,gamma,
     & pi,rotx,rtod,sina,sinb,sing,twopi
      dimension rotx(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      alf=dtor*alpha
      bet=dtor*beta
      gam=dtor*gamma
      cosa=cos(alf)
      cosb=cos(bet)
      cosg=cos(gam)
      sina=sin(alf)
      sinb=sin(bet)
      sing=sin(gam)
      rotx(1,1)= cosa*cosb*cosg-sina*sing
      rotx(1,2)=-cosa*cosb*sing-sina*cosg
      rotx(1,3)= cosa*sinb
      rotx(2,1)= sina*cosb*cosg+cosa*sing
      rotx(2,2)=-sina*cosb*sing+cosa*cosg
      rotx(2,3)= sina*sinb
      rotx(3,1)=     -sinb*cosg
      rotx(3,2)=      sinb*sing
      rotx(3,3)=      cosb
c     write(6,'(a,e15.5)') ' trace ',(1+cosb)*(1+cos(alf+gam))-1
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function length(card)
      implicit none
      integer i,length
      character card*80
       do i=80,1,-1
      if(card(i:i).ne.' ') go to 10
       enddo
10    length=i
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine compact(card)
      implicit none
      integer i,j,length
      character card*80
      logical blanco,white
      external length
      blanco=card(1:1).eq.' '
      j=0
      do 10 i=1,length(card)
      white=card(i:i).eq.' '
      if(white.and.blanco) go to 10
         if(white) then
      blanco=.true.
         else
      blanco=.false.
         endif
      j=j+1
      card(j:j)=card(i:i)
10    continue
       do i=j+1,length(card)
      card(i:i)=' '
       enddo
      return
      end
