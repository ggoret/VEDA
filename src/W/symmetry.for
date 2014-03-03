c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program symmetry
      implicit none
      double precision dtor,pi,rtod,twopi
      character card*80
      external length
      common/angkte/ pi,twopi,dtor,rtod
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=pi*2.d0
      dtor=pi/180.d0
      rtod=180.d0/pi
      write(6,'(/a)') ' enter symmetry: ico hel cn dn p1 other'
      read(5,'(a)') card
      call compact(card)
         if(card(1:3).eq.'ico') then
      call icosim
         else if(card(1:3).eq.'hel') then
      write(6,'(/a)') ' enter parameterization: ctu ele'
      read(5,'(a)') card
            if(card(1:3).eq.'ctu') then
      call helix
            else
      call tubes
            endif
         else if(card(1:2).eq.'cn') then
      call cn
         else if(card(1:2).eq.'dn') then
      call dn
         else if(card(1:5).eq.'other') then
      call other
         else
      write(6,'(a)') ' symmetry set to P1'
      call p1
         endif
      stop
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
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine icosim
      implicit none
      integer S,i,ii,ine,inp,is,length,n,ns,o1,o2
      double precision alpha,ang,ang2,ang3,angles,beta,dtor,gamma,omega,
     & pi,pio2,rn,rtod,rv,sim,theta,twopi,un,uv,xhi,zero
      character card*80,cart*80,formo*80
      logical unit
      external length
      dimension angles(3,60),ii(60,60),rn(3,3),rv(9),sim(3,3,60),un(3),
     & uv(3)
      common/angkte/ pi,twopi,dtor,rtod
      data o1,o2/1,2/
c========================
      zero=0.d0
      pio2=pi/2.d0
      ang2=acos(1.d0/sqrt(5.d0))/2.d0
      ang3=acos(sin(pi*3.d0/5.d0)*(1.d0+sqrt(5.d0))/sqrt(15.d0))
c========================
      call euler(angles,ns)
      write(6,'(6(/a))') ' ENTER option:',
     . ' 5Z2Y.1|2 5Z2X.1|2 5Y2Z.1|2 5Y2X.1|2 5X2Z.1|2 5X2Y.1|2',
     . ' 3Z2Y.1|2 3Z2X.1|2 3Y2Z.1|2 3Y2X.1|2 3X2Z.1|2 3X2Y.1|2',
     . ' 2Z2Y.1|2',
     . ' (ROTA BTV = 3Z2X.1 ; SFV TBE = 5Y2Z.1 ; IBDV = 2Z2Y.1)',
     . ' or rotation: [p ang cos/x/y/z/] or [e alpha/beta/gamma]'
      read(5,'(a)') card
         if(length(card).eq.0) then
      unit=.true.
      stop
         else
      unit=.false.
            if(card(1:6).eq.'5Z2Y.1') then
      unit=.true.
            else if(card(1:6).eq.'5Z2Y.2') then
      call rmxe(pi,pi,zero,rv)
            else if(card(1:6).eq.'5Z2X.1') then
      call rmxe(pio2,zero,zero,rv)
            else if(card(1:6).eq.'5Z2X.2') then
      call rmxe(-pio2,zero,zero,rv)
            else if(card(1:6).eq.'5Y2Z.1') then
      call rmxe(pio2,pio2,pio2,rv)
            else if(card(1:6).eq.'5Y2Z.2') then
      call rmxe(pio2,pio2,-pio2,rv)
            else if(card(1:6).eq.'5Y2X.1') then
      call rmxe(pio2,pio2,pi,rv)
            else if(card(1:6).eq.'5Y2X.2') then
      call rmxe(pio2,pio2,zero,rv)
            else if(card(1:6).eq.'5X2Z.1') then
      call rmxe(zero,pio2,pio2,rv)
            else if(card(1:6).eq.'5X2Z.2') then
      call rmxe(zero,pio2,-pio2,rv)
            else if(card(1:6).eq.'5X2Y.1') then
      call rmxe(zero,pio2,pi,rv)
            else if(card(1:6).eq.'5X2Y.2') then
      call rmxe(zero,pio2,zero,rv)
            else if(card(1:6).eq.'3Z2Y.1') then
      call rmxe(pi,ang3,pi,rv)
            else if(card(1:6).eq.'3Z2Y.2') then
      call rmxe(zero,ang3,pi,rv)
            else if(card(1:6).eq.'3Z2X.1') then
      call rmxe(pio2,ang3,pi,rv)
            else if(card(1:6).eq.'3Z2X.2') then
      call rmxe(-pio2,ang3,pi,rv)
            else if(card(1:6).eq.'3Y2Z.1') then
      call rmxe(pio2-ang3,pio2,pio2,rv)
            else if(card(1:6).eq.'3Y2Z.2') then
      call rmxe(3*pio2+ang3,pio2,-pio2,rv)
            else if(card(1:6).eq.'3Y2X.1') then
      call rmxe(pio2,pio2+ang3,pi,rv)
            else if(card(1:6).eq.'3Y2X.2') then
      call rmxe(pio2,pio2-ang3,zero,rv)
            else if(card(1:6).eq.'3X2Z.1') then
      call rmxe(pi-ang3,pio2,pio2,rv)
            else if(card(1:6).eq.'3X2Z.2') then
      call rmxe(ang3,pio2,-pio2,rv)
            else if(card(1:6).eq.'3X2Y.1') then
      call rmxe(pi,pio2-ang3,zero,rv)
            else if(card(1:6).eq.'3X2Y.2') then
      call rmxe(zero,pio2-ang3,zero,rv)
            else if(card(1:6).eq.'2Z2Y.1') then
      call rmxe(zero,ang2,zero,rv)
            else if(card(1:6).eq.'2Z2Y.2') then
      call rmxe(pi,pio2-ang2,pi,rv)
            else
      ine=index(card,'e')
      inp=index(card,'p')
               if(ine.ge.1) then
      card(ine:ine)=' '
      read(card,*) alpha,beta,gamma
      call rmxe(alpha*dtor,beta*dtor,gamma*dtor,rv)
               else if(inp.ge.1) then
      card(inp:inp)=' '
      read(card,*) ang,uv
      call rmxp(ang*dtor,uv,rv)
               else
      write(6,'(a)') ' invalid option'
      stop
               endif
            endif
         endif
         if(unit) then
      alpha=0.
      beta=0.
      gamma=0.
         else
      call rmx2e(rv,alpha,beta,gamma)
      alpha=alpha*rtod
      beta=beta*rtod
      gamma=gamma*rtod
         endif
      write(6,'(2a,3f12.6)') card(1:6),'; rotation:',alpha,beta,gamma
c========================
      open(unit=o1,file='sym',form='formatted',status='unknown')
      open(unit=o2,file='gs.sym',form='formatted',status='unknown')
cjn   open(unit=3,file='SYM',form='formatted',status='unknown')
       do n=1,80
      card(n:n)=' '
      cart(n:n)=' '
      formo(n:n)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      S=0
      do 10 is=1,ns
      alpha=angles(1,is)
      beta=angles(2,is)
      gamma=angles(3,is)
      call rmxe(alpha,beta,gamma,rn)
      call rmx2p(rn,xhi,un)
         if(.not.unit) then
      call promv(uv,rv,un)
       do n=1,3
      un(n)=uv(n)
       enddo
      call rmxp(xhi,un,rn)
      call rmx2e(rn,alpha,beta,gamma)
         endif
      alpha=alpha*rtod
      beta=beta*rtod
      gamma=gamma*rtod
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') alpha,beta,gamma,0.,0.,0.
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
cjn O usa la transpuesta
       do n=1,3
      write(o2,'(3f12.6)') (rn(i,n),i=1,3)
       enddo
      write(o2,'(3f12.6)') 0.,0.,0.
      call v2a(un,theta,omega)
      xhi=xhi*rtod
      theta=theta*rtod
      omega=omega*rtod
         if(xhi.lt.1.d-4) then
      theta=0.d0
      omega=0.d0
      un(1)=0.d0
      un(2)=0.d0
      un(3)=0.d0
         endif
         if(is.lt.10) then
      formo="('#',i1,'  ',3f10.5,' ; ',f8.3,3f10.5,' #',i1)"
         else
      formo="('#',i2,' ',3f10.5,' ; ',f8.3,3f10.5,' #',i2)"
         endif
cjn   write(6,fmt=formo) is,xhi,theta,omega,xhi,un,is
       do n=1,3
        do i=1,3
      sim(n,i,is)=rn(n,i)
        enddo
       enddo
10    continue
c========================
      call simtes(ns,sim,ii)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dn
      implicit none
      integer N,S,i,j,length,o1,o2
      double precision alpha,beta,cosa,dang,dtor,gamma,pi,rotx,rtod,
     & sina,twopi,u
      character card*80,cart*80
      external length
      dimension rotx(3,3),u(3)
      common/angkte/ pi,twopi,dtor,rtod
      data o1,o2/1,2/
c========================
      write(6,'(a)') ' Enter order of n-fold symmetry:'
      read(5,*) N
c========================
      open(unit=1,file='sym',form='formatted',status='unknown')
      open(unit=2,file='gs.sym',form='formatted',status='unknown')
       do i=1,80
      card(i:i)=' '
      cart(i:i)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      S=0
      do 10 j=0,6
      dang=j*twopi/N
      cosa=cos(dang)
      sina=sin(dang)
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') j*360.d0/N,0.,0.,0.,0.,0.
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else if(S.lt.100) then
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         else
      write(card(66:68),'(i3)') S
      write(cart(12:14),'(i3)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') cosa,sina,0.
      write(o2,'(3f12.6)') -sina,cosa,0.
      write(o2,'(3f12.6)') 0.,0.,1.
      write(o2,'(3f12.6)') 0.,0.,0.
10    continue
      do 20 j=0,6
      dang=j*twopi/N
      u(1)=cos(dang)
      u(2)=sin(dang)
      u(3)=0.d0
      call rmxp(pi,u,rotx)
      call rmx2e(rotx,alpha,beta,gamma)
      alpha=alpha*rtod
      beta=beta*rtod
      gamma=gamma*rtod
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') alpha,beta,gamma,0.,0.,0.
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else if(S.lt.100) then
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         else
      write(card(66:68),'(i3)') S
      write(cart(12:14),'(i3)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') (rotx(1,i),i=1,3)
      write(o2,'(3f12.6)') (rotx(2,i),i=1,3)
      write(o2,'(3f12.6)') (rotx(3,i),i=1,3)
      write(o2,'(3f12.6)') 0.,0.,0.
20    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tubes
      implicit none
      integer N,S,i,j,length,o1,o2,signo,start
      double precision cosa,dang,dtor,pi,rtod,sina,twopi,w1,w2,x
      character card*80,cart*80
      dimension x(2)
      common/angkte/ pi,twopi,dtor,rtod
      data o1,o2/1,2/
c========================
      write(6,'(/a,a,5(/a))') ' enter elementary helix shifts',
     . ' (degrees , Angstroms), number of starts + N:',
     . ' tubuline (dimero) :     128.583 3.005 1'
      read(5,*) x,start,N
      write(6,'(/a)') ' base vector:'
      write(6,'(a,2f15.5,a,i5)') ' x1,x2:',x(1),x(2),' ; starts:',start
      x(1)=x(1)*dtor/twopi
c========================
      open(unit=o1,file='sym',form='formatted',status='unknown')
      open(unit=o2,file='gs.sym',form='formatted',status='unknown')
       do i=1,80
      card(i:i)=' '
      cart(i:i)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      write(card(1:60),'(6f10.3)') 0.,0.,0.,0.,0.,0.
      S=1
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') 1.,0.,0.
      write(o2,'(3f12.6)') 0.,1.,0.
      write(o2,'(3f12.6)') 0.,0.,1.
      write(o2,'(3f12.6)') 0.,0.,0.
       do j=0,start-1
        do signo=1,-1,-2
         do i=1,N
      w1=signo*i*x(1)+(twopi/start)*j
      w2=signo*i*x(2)
      w1=mod(1.d0+mod(w1,1.d0),1.d0)
      dang=twopi*w1
      cosa=cos(dang)
      sina=sin(dang)
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') w1*360,0.,0.,0.,0.,w2
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else if(S.lt.100) then
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         else if(S.lt.1000) then
      write(card(66:68),'(i3)') S
      write(cart(12:14),'(i3)') S
         else
      write(card(66:69),'(i4)') S
      write(cart(12:15),'(i4)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') cosa,sina,0.
      write(o2,'(3f12.6)') -sina,cosa,0.
      write(o2,'(3f12.6)') 0.,0.,1.
      write(o2,'(3f12.6)') 0.,0.,w2
         enddo
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine cn
      implicit none
      integer N,S,i,j,length,o1,o2
      double precision cosa,dang,dtor,pi,rtod,sina,twopi
      character card*80,cart*80
      external length
      common/angkte/ pi,twopi,dtor,rtod
      data o1,o2/1,2/
c========================
      write(6,'(a)') ' Enter order of n-fold symmetry:'
      read(5,*) N
c========================
      open(unit=o1,file='sym',form='formatted',status='unknown')
      open(unit=o2,file='gs.sym',form='formatted',status='unknown')
       do i=1,80
      card(i:i)=' '
      cart(i:i)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      S=0
      do 10 j=0,6
      dang=j*twopi/N
      cosa=cos(dang)
      sina=sin(dang)
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') j*360.d0/N,0.,0.,0.,0.,0.
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else if(S.lt.100) then
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         else
      write(card(66:68),'(i3)') S
      write(cart(12:14),'(i3)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') cosa,sina,0.
      write(o2,'(3f12.6)') -sina,cosa,0.
      write(o2,'(3f12.6)') 0.,0.,1.
      write(o2,'(3f12.6)') 0.,0.,0.
10    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine helix
      implicit none
      integer A,B,C,D,M,N,S,T,U,i,j,k,length,o1,o2
      double precision cosa,dang,det,dtor,norme,pi,rtod,scale,sina,
     & twopi,u1,u2,w1,w2,x
      character card*80,cart*80
      dimension u1(2),u2(2),x(2)
      common/angkte/ pi,twopi,dtor,rtod
      data o1,o2/1,2/
c========================
1     write(6,'(6(/a))') ' enter the helix parameters (c,T,U) + N:',
     . ' petit tubes :  820.0    -9   103',
     . ' gros tubes  : 1220.0   -35   293',
     . ' ron tubes   : 1440.0   -25    54',
     . ' ken tubes   :  770.0   -13    28',
     . ' tubuline (dimero) :      840.0   -31    76'
      read(5,*) scale,T,U,N
c========================
      u1(1)=T/real(U)
      u1(2)=1/real(U)
      i=abs(U/T)
      u2(1)=i*u1(1)-sign(1,T)
      u2(2)=i*u1(2)
      det=u1(1)*u2(2)-u1(2)*u2(1)
         if(det.eq.0.0) then
      write(6,'(/20x,a)') ' >>> incorrect base vectors <<<'
      go to 1
         else if(det.lt.0.0) then
      write(6,'(20x,a)') ' >>> vectors permuted <<<'
      w1=u1(1)
      w2=u1(2)
      u1(1)=u2(1)
      u1(2)=u2(2)
      u2(1)=w1
      u2(2)=w2
         endif
      write(6,'(/a)') ' base vectors:'
      write(6,'(a,f15.5,a,f15.5,a/a,f15.5,a,f15.5,a)')
     . ' u1 , u2 :',u1(1)*360,' degrees , ',u1(2)*scale, ' Angstroms',
     . ' v1 , v2 :',u2(1)*360,' degrees , ',u2(2)*scale, ' Angstroms'
c========================
      det=u1(1)*u2(2)-u1(2)*u2(1)
      A=nint(u2(2)/det)
      B=nint(-u1(2)/det)
      write(6,'(/a/a/a,f10.5/a,f10.5/a,i5,a,i5)')
     . ' COEFFS that produce an ANGULAR repeat:',
     . ' >>> ZERO CHECK <<<',
     . ' A * u1 + B * v1 - 1 =',A*u1(1)+B*u2(1)-1,
     . ' A * u2 + B * v2     =',A*u1(2)+B*u2(2),
     . ' A =',A,' , B =',B
         if(B.eq.0) then
      write(6,'(/20x,a)') ' >>> null B <<<'
      go to 1
         endif
c------
      i=A
      j=B
2        if(mod(i,j).ne.0) then
      k=mod(i,j)
      i=j
      j=k
      go to 2
         endif
      M=j
c------
       do i=1,A,sign(1,A)
         if(mod(M+B*i,A).eq.0) then
      C=i
      D=(M+B*C)/A
         endif
       enddo
      k=nint(1.e9)
      norme=1.e10
       do i=-50,50
      w1=u1(1)*(C+A*i)+u2(1)*(D+B*i)
         if(abs(w1).le.norme) then
      k=i
      norme=abs(w1)
         endif
       enddo
      C=C+A*k
      D=D+B*k
         if(C*u1(2)+D*u2(2).lt.0.0) then
      C=-C
      D=-D
      M=-M
         endif
c------
      write(6,'(a,i5)') ' GREATEST common divisor: M =',M
      x(1)=C*u1(1)+D*u2(1)
      x(2)=C*u1(2)+D*u2(2)
      write(6,'(/a/a,f15.5,a/a,f15.5,a/a,i5,a,i5)')
     . ' ELEMENTARY helix:',
     . ' x1 = C * u1 + D * v1 =',x(1)*360,' degrees',
     . ' x2 = C * u2 + D * v2 =',x(2)*scale,' Angstroms',
     . ' C =',C,' , D =',D
      write(6,'(a/a,i5)') ' >>> ZERO CHECK <<<',
     . ' A * D - B * C - M =',A*D-B*C-M
      write(6,'(/a,2f15.5)') ' ELEMENTARY helix:',x(1)*360,x(2)*scale
c========================
      open(unit=o1,file='sym',form='formatted',status='unknown')
      open(unit=o2,file='gs.sym',form='formatted',status='unknown')
       do i=1,80
      card(i:i)=' '
      cart(i:i)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      S=0
      do 10 i=0,2*U
      w1=i*x(1)
      w2=i*x(2)
      w1=mod(1.d0+mod(w1,1.d0),1.d0)
      w2=mod(w2,1.d0)
      if(mod(i,U).eq.0.and.i.ne.0) go to 10
      if(i.gt.U) w2=w2-1.d0
      dang=twopi*w1
      cosa=cos(dang)
      sina=sin(dang)
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') w1*360,0.,0.,0.,0.,w2*scale
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else if(S.lt.100) then
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         else if(S.lt.1000) then
      write(card(66:68),'(i3)') S
      write(cart(12:14),'(i3)') S
         else
      write(card(66:69),'(i4)') S
      write(cart(12:15),'(i4)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') cosa,sina,0.
      write(o2,'(3f12.6)') -sina,cosa,0.
      write(o2,'(3f12.6)') 0.,0.,1.
      write(o2,'(3f12.6)') 0.,0.,w2*scale
10    continue
c========================
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine p1
      implicit none
      integer i,length,o1,o2
      character card*80,cart*80
      external length
      data o1,o2/1,2/
c========================
      open(unit=o1,file='sym',form='formatted',status='unknown')
      open(unit=o2,file='gs.sym',form='formatted',status='unknown')
       do i=1,80
      card(i:i)=' '
      cart(i:i)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      write(card(1:60),'(6f10.3)') 0.,0.,0.,0.,0.,0.
      write(card(66:66),'(i1)') 1
      write(cart(12:12),'(i1)') 1
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
      write(o2,'(3f12.6)') 1.,0.,0.
      write(o2,'(3f12.6)') 0.,1.,0.
      write(o2,'(3f12.6)') 0.,0.,1.
      write(o2,'(3f12.6)') 0.,0.,0.
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
      subroutine rmxp(ang,un,rotx)
      implicit none
      integer eijk,i,j,k
      double precision ang,cang,rotx,sang,un
      dimension eijk(3,3,3),rotx(3,3),un(3)
      data eijk/
     & 0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
      cang=cos(ang)
      sang=sin(ang)
       do j=1,3
        do i=1,3
      rotx(i,j)=(1.d0-cang)*un(i)*un(j)
      if(i.eq.j) rotx(i,j)=rotx(i,j)+cang
         do k=1,3
         if(eijk(i,k,j).eq.1) then
      rotx(i,j)=rotx(i,j)+sang*un(k)
         else if(eijk(i,k,j).eq.-1) then
      rotx(i,j)=rotx(i,j)-sang*un(k)
         endif
         enddo
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmx2e(rotx,alpha,beta,gamma)
      implicit none
      double precision alpha,ang,beta,cang,dtor,fuzz,gamma,pi,rotx,rtod,
     & sang,twopi
      dimension rotx(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      fuzz=0.001
      cang=min(max(rotx(3,3),-1.d0),1.d0)
      ang=acos(cang)
      sang=sin(ang)
         if(sang.gt.fuzz) then
      alpha=atan2(rotx(2,3),rotx(1,3))
      gamma=atan2(rotx(3,2),-rotx(3,1))
         else
      alpha=atan2(-rotx(1,2),rotx(1,1)*rotx(3,3))
      gamma=0.
         endif
      alpha=mod(twopi+alpha,twopi)
      beta=ang
      gamma=mod(twopi+gamma,twopi)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cjn eje 5 || Z , 2 || Y , 3 NX
      subroutine euler(angles,ns)
      implicit none
      integer alpha,beta,gamma,i,ns
      double precision angles,pio5,val
      dimension alpha(60),angles(3,60),beta(60),gamma(60),val(4)
      data alpha/
     & 0,0,0,1,5,5,9,9,3,7,1,3,7,0,0,2,4,6,8,0,2,8,0,4,6,1,3,7,9,3,
     & 5,5,7,9,1,0,4,6,0,2,6,4,8,8,2,3,7,1,9,5,2,8,4,6,0,0,0,0,0,0/
      data beta/
     & 1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,
     & 2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4/
      data gamma/
     & 0,2,8,0,4,6,0,2,6,4,8,8,2,4,6,1,3,7,9,3,5,5,7,9,1,2,4,6,8,0,
     & 2,8,0,4,6,1,5,5,9,9,3,7,1,3,7,2,8,4,6,0,3,7,1,9,5,0,2,8,4,6/
      ns=60
      val(1)=0.d0
      val(2)=acos(1.d0/sqrt(5.d0))
      val(3)=atan2(1.d0,1.d0)*4.d0-acos(1.d0/sqrt(5.d0))
      val(4)=atan2(1.d0,1.d0)*4.d0
      pio5=atan2(1.d0,1.d0)*4.d0/5.d0
       do i=1,60
      angles(1,i)=alpha(i)*pio5
      angles(2,i)=val(beta(i))
      angles(3,i)=gamma(i)*pio5
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmxe(alpha,beta,gamma,rotx)
      implicit none
      double precision alpha,beta,cosa,cosb,cosg,gamma,rotx,sina,sinb,
     & sing
      dimension rotx(3,3)
      cosa=cos(alpha)
      sina=sin(alpha)
      cosb=cos(beta)
      sinb=sin(beta)
      cosg=cos(gamma)
      sing=sin(gamma)
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
      subroutine rmx2p(rotx,ang,un)
      implicit none
      integer i
      double precision ang,cang,fuzz,rotx,sang,scale,un
      dimension rotx(3,3),un(3)
      fuzz=0.001
      cang=min(max((rotx(1,1)+rotx(2,2)+rotx(3,3)-1.d0)/2.d0,-1.d0),
     . 1.d0)
      ang=acos(cang)
      sang=sin(ang)
         if(sang.gt.fuzz) then
      scale=2*sang
      un(1)=(rotx(3,2)-rotx(2,3))/scale
      un(2)=(rotx(1,3)-rotx(3,1))/scale
      un(3)=(rotx(2,1)-rotx(1,2))/scale
         else
            if(cang.gt.0.) then
      un(1)=0.d0
      un(2)=0.d0
      un(3)=1.d0
            else
       do i=1,3
      un(i)=sqrt(max((rotx(i,i)+1)/2,0.d0))
       enddo
               if(un(3).gt.fuzz) then
      un(1)=sign(un(1),rotx(1,3))
      un(2)=sign(un(2),rotx(2,3))
               else if(un(2).gt.fuzz) then
      un(1)=sign(un(1),rotx(1,2))
               endif
            endif
         endif
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine v2a(un,theta,phi)
      implicit none
      double precision cthe,fuzz,phi,sthe,theta,twopi,un,unorm
      dimension un(3)
      twopi=atan2(1.d0,1.d0)*8.d0
      fuzz=0.0001
      unorm=sqrt(un(1)**2+un(2)**2+un(3)**2)
      cthe=min(max(un(3)/unorm,-1.d0),1.d0)
      sthe=sqrt(1-cthe**2)
         if(sthe.gt.fuzz) then
      phi=atan2(un(2),un(1))
         else
      phi=0.
         endif
      theta=acos(cthe)
      phi=mod(phi+twopi,twopi)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine promv(uo,r1,u1)
      implicit none
      integer i,j
      double precision r1,sum,u1,uo
      dimension r1(3,3),u1(3),uo(3)
       do i=1,3
      sum=0.d0
        do j=1,3
      sum=sum+r1(i,j)*u1(j)
        enddo
      uo(i)=sum
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine simtes(ns,sim,mult)
      implicit none
      integer i,j,k,l,mult,ns
      double precision ang,dist,error,rn,sim
      dimension mult(ns,ns),rn(9),sim(9,ns)
      error=0.d0
       do i=1,ns
        do j=1,ns
      call pro2mx(rn(1),sim(1,j),sim(1,i))
      dist=1.d0
      l=0
         do k=1,ns
      call dista(rn(1),sim(1,k),ang)
         if(ang.le.dist) then
      dist=ang
      l=k
         endif
         enddo
      mult(i,j)=l
      if(dist.ge.error) error=dist
        enddo
       enddo
      write(6,'(/a,f8.5)') ' Symmetry violation:',error
       do i=1,ns
        do j=1,ns
         do k=1,ns
         if(j.ne.k.and.mult(i,j).eq.mult(i,k)) then
      write(6,'(a,i3,a,2i3)') ' elementos iguales en linea',i,':',j,k
         endif
         enddo
        enddo
       enddo
       do i=1,ns
        do j=1,ns
         do k=1,ns
         if(j.ne.k.and.mult(j,i).eq.mult(k,i)) then
      write(6,'(a,i3,a,2i3)') ' elementos iguales en columna',i,':',j,k
         endif
         enddo
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pro2mx(ro,r2,r1)
      implicit none
      integer i,j,k
      double precision r1,r2,ro,sum
      dimension r1(3,3),r2(3,3),ro(3,3)
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+r2(i,k)*r1(k,j)
         enddo
      ro(i,j)=sum
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dista(r0,r1,sum)
      implicit none
      integer i
      double precision r0,r1,sum
      dimension r0(9),r1(9)
      sum=0.d0
       do i=1,9
      sum=sum+(r0(i)-r1(i))**2
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine other
      implicit none
      integer S,i,length,n,o1,o2
      double precision alpha,beta,dtor,gamma,pi,rn,rtod,twopi
      character card*80,cart*80,lect*80
      external length
      dimension rn(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      data o1,o2/1,2/
      open(unit=o1,file='sym',form='formatted',status='unknown')
      open(unit=o2,file='gs.sym',form='formatted',status='unknown')
       do n=1,80
      card(n:n)=' '
      cart(n:n)=' '
       enddo
      card(65:65)='#'
      cart(1:40)='.LSQ_RT_SYM       R     12  (3F12.6)    '
      S=0
c========================
      write(6,'(a)') ' Enter euler angles:'
10    read(5,'(a)',end=20) lect
      if(length(lect).eq.0) go to 20
      read(lect,*) alpha,beta,gamma
      card(66:70)='     '
      cart(12:16)='     '
      write(card(1:60),'(6f10.3)') alpha,beta,gamma,0.,0.,0.
      alpha=alpha*dtor
      beta=beta*dtor
      gamma=gamma*dtor
      call rmxe(alpha,beta,gamma,rn)
      S=S+1
         if(S.lt.10) then
      write(card(66:66),'(i1)') S
      write(cart(12:12),'(i1)') S
         else
      write(card(66:67),'(i2)') S
      write(cart(12:13),'(i2)') S
         endif
      write(o1,'(a)') card(1:length(card))
      write(o2,'(a)') cart(1:length(cart))
cjn O usa la transpuesta
       do n=1,3
      write(o2,'(3f12.6)') (rn(i,n),i=1,3)
       enddo
      write(o2,'(3f12.6)') 0.,0.,0.
      go to 10
20    return
      end
