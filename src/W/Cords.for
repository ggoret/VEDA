c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program cords
      implicit none
      integer i,i1,i2,iecr,ilec,iout,j,kprt,length,na
      double precision al,be,com,dtor,ga,pi,rmx,rtod,trx,twopi,xi,xo
      character card*80,cbid1*1,forml*80,inv*1,lun*80,mas*1
      external length
      dimension com(3),rmx(3,3),trx(3),xi(3),xo(3)
      common/angkte/ pi,twopi,dtor,rtod
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      data i1,i2/98,99/
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      ilec=5
      iecr=6
      kprt=1
       do i=1,80
      forml(i:i)=' '
       enddo
      forml(1:11)='(30x,3f8.3)'
10    write(iecr,'(/a)') ' enter the input-coord filename:'
      read(ilec,'(a)',end=70) lun
      if(length(lun).eq.0) go to 80
      open(unit=i1,file=lun,form='formatted',status='old')
      na=0
       do i=1,3
      com(i)=0.
       enddo
20    read(i1,'(a)',end=30) card
      if(card(1:4).ne.'ATOM') go to 20
      na=na+1
      read(card,fmt=forml) xi
       do i=1,3
      com(i)=com(i)+xi(i)
       enddo
      go to 20
30     do i=1,3
      com(i)=com(i)/na
       enddo
      write(iecr,'(a)') ' enter the output-coord filename:'
      read(ilec,'(a)') lun
      open(unit=i2,file=lun,form='formatted',status='new')
      write(iecr,'(a)') ' rotation in the center-of-mass [N/y]:'
      read(ilec,'(a)') mas
      if(length(mas).eq.0) mas='n'
      write(iecr,'(a)') ' enter rotation and translation:'
      read(ilec,*) al,be,ga,trx
      write(iecr,'(a)') ' inversion [N/y]:'
      read(ilec,'(a)') inv
      if(length(inv).eq.0) inv='n'
      call rmxe(al,be,ga,rmx)
         if(inv.eq.'y') then
       do i=1,3
        do j=1,3
      rmx(i,j)=-rmx(i,j)
        enddo
       enddo
         endif
      write(iecr,'(/a,3(/3f15.5))') ' matrix M',((rmx(i,j),j=1,3),i=1,3)
      write(iecr,'(/a/3f15.3)') ' translation T',(trx(i),i=1,3)
         if(mas.eq.'y') then
       do i=1,3
      trx(i)=trx(i)+com(i)
       enddo
         endif
      rewind(unit=i1)
40    read(i1,'(a)') card
         if(card(1:4).ne.'ATOM') then
      write(i2,'(a)') card
      go to 40
         endif
      go to 60
50    read(i1,'(a)',end=70) card
         if(card(1:4).ne.'ATOM') then
      if(card(1:3).eq.'TER') write(i2,'(a)') 'TER'
      go to 50
         endif
60    read(card,fmt=forml) xi
         if(mas.eq.'y') then
       do i=1,3
      xi(i)=xi(i)-com(i)
       enddo
         endif
       do i=1,3
      xo(i)=0.
        do j=1,3
      xo(i)=xo(i)+rmx(i,j)*xi(j)
        enddo
      xo(i)=xo(i)+trx(i)
       enddo
      write(card(31:54),'(3f8.3)') xo
      write(i2,'(a)') card
      go to 50
70    close(unit=i1)
      close(unit=i2)
      go to 10
80    stop
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
      subroutine rmxe(alpha,beta,gamma,rotx)
      implicit none
      integer iecr,ilec,iout,kprt
      double precision alf,alpha,bet,beta,cosa,cosb,cosg,dtor,gam,gamma,
     & pi,rotx,rtod,sina,sinb,sing,twopi
      character cbid1*1
      dimension rotx(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      if(kprt.ne.0) write(iecr,'(a,3f10.3)')
     . ' euler matrix; alpha,beta,gamma',alpha,beta,gamma
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
