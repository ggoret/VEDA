c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 DOS VERSIONES DEL MAIN
      program rt
      implicit none
      integer i,iecr,ilec,ind,is,length,lun,ms,ns
      double precision al,be,dtor,ga,mx,pi,rms,rmx,rtod,trs,trx,twopi,tx
      character card*80
      external length
      parameter(ms=1000)
      dimension mx(9),rms(9,ms),rmx(9),trs(3,ms),trx(3),tx(3)
      common/angkte/ pi,twopi,dtor,rtod
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      ilec=5
      iecr=6
      lun=1
      open(unit=lun,file='SYM',form='formatted',status='unknown')
      ns=0
10    read(lun,'(a)',end=20) card
      if(length(card).eq.0) go to 20
      ns=ns+1
      read(card,*) al,be,ga,(trs(i,ns),i=1,3)
      call rmxe(al,be,ga,rms(1,ns))
      go to 10
20    continue
c     write(iecr,'(/a)') ' RT of the fixed molecule:'
      read(ilec,*) al,be,ga,trx
      call rmxe(al,be,ga,rmx(1))
c     write(iecr,'(/a)') ' symmetry identifiers:'
30    read(ilec,'(a)',end=40) card
      if(length(card).eq.0) go to 40
      ind=index(card,'#')
      card(ind:ind)=' '
      read(card,*) is
      call pro2mx(mx(1),rms(1,is),rmx(1))
      call rmx2e(mx(1),al,be,ga)
      call promv(tx(1),rms(1,is),trx(1))
       do i=1,3
      tx(i)=tx(i)+trs(i,is)
       enddo
      write(iecr,'(a,i5)') ' SymMate #',is
      write(iecr,'(a,3f15.5)') ' Rotation :   ',al,be,ga
      write(iecr,'(a,3f15.3)') ' Translation :',tx
      go to 30
40    stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program rt
      implicit none
      integer i,iecr,ilec,ind,is,length,lun,ms,sym
      double precision al,be,dtor,ga,mx,pi,rms,rmx,rtod,trs,trx,twopi,tx
      character card*80
      external length
      parameter(ms=1000)
      dimension mx(9),rms(9,ms),rmx(9),sym(ms),trs(3,ms),trx(3),tx(3)
      common/angkte/ pi,twopi,dtor,rtod
      common/ioprg/ ilec,iecr
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      ilec=5
      iecr=6
      lun=1
      open(unit=lun,file='sym',form='formatted',status='unknown')
       do is=1,ms
      sym(is)=0
        do i=1,9
      rms(i,is)=0.
        enddo
        do i=1,3
      trs(i,is)=0.
        enddo
       enddo
10    read(lun,'(a)',end=20) card
      if(length(card).eq.0) go to 20
      ind=index(card,'#')
         if(ind.le.0) then
      write(iecr,'(/a/)') ' >> ERROR << missing sym-ID'
      stop
         else
      read(card(ind+1:80),*) is
            if(is.gt.ms) then
      write(iecr,'(/a/)') ' >> ERROR << increase ms'
      stop
            endif
         endif
      sym(is)=1
      read(card,*) al,be,ga,(trs(i,is),i=1,3)
      call rmxe(al,be,ga,rms(1,is))
      go to 10
20    continue
c     write(iecr,'(/a)') ' RT of the fixed molecule:'
      read(ilec,*) al,be,ga,trx
      call rmxe(al,be,ga,rmx(1))
c     write(iecr,'(/a)') ' symmetry identifiers:'
30    read(ilec,'(a)',end=40) card
      if(length(card).eq.0) go to 40
      ind=index(card,'#')
      card(ind:ind)=' '
      read(card,*) is
         if(sym(is).ne.1) then
      write(iecr,'(/a/)') ' sym-ID not in list'
      go to 30
         endif
      call pro2mx(mx(1),rms(1,is),rmx(1))
      call rmx2e(mx(1),al,be,ga)
      call promv(tx(1),rms(1,is),trx(1))
       do i=1,3
      tx(i)=tx(i)+trs(i,is)
       enddo
      write(iecr,'(a,i5)') ' SymMate #',is
      write(iecr,'(a,3f15.5)') ' Rotation :   ',al,be,ga
      write(iecr,'(a,3f15.3)') ' Translation :',tx
      go to 30
40    stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmxe(alpha,beta,gamma,rotx)
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
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmx2e(rotx,alpha,beta,gamma)
      implicit none
      double precision alpha,ang,beta,cang,dtor,fuzz,gamma,pi,rotx,rtod,
     & sang,twopi
      dimension rotx(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      fuzz=0.01
c     cang=min(max(rotx(3,3),-1.d0),1.d0)
c     ang=acos(cang)
      cang=sqrt(min(max(2.d0-(rotx(1,3)**2+rotx(2,3)**2+rotx(3,1)**2+
     . rotx(3,2)**2+rotx(3,3)**2),0.d0),1.d0))
      ang=acos(sign(cang,rotx(3,3)))
      sang=sin(ang)
         if(sang.gt.fuzz) then
      alpha=atan2(rotx(2,3),rotx(1,3))
      gamma=atan2(rotx(3,2),-rotx(3,1))
         else
      alpha=atan2(-rotx(1,2),rotx(1,1)*rotx(3,3))
      gamma=0.
         endif
      alpha=mod(twopi+alpha,twopi)*rtod
      beta=ang*rtod
      gamma=mod(twopi+gamma,twopi)*rtod
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pro2mx(ro,r2,r1)
      implicit none
      integer i,j,k
      double precision add,r1,r2,ro
      dimension r1(3,3),r2(3,3),ro(3,3)
       do i=1,3
        do j=1,3
      add=0.d0
         do k=1,3
      add=add+r2(i,k)*r1(k,j)
         enddo
      ro(i,j)=add
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine promv(uo,r2,u1)
      implicit none
      integer i,j
      double precision add,r2,u1,uo
      dimension r2(3,3),u1(3),uo(3)
       do i=1,3
      add=0.d0
        do j=1,3
      add=add+r2(i,j)*u1(j)
        enddo
      uo(i)=add
       enddo
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
