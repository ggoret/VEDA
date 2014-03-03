c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program refer
      implicit none
      integer i,i1,iecr,ierr,ilec,ind,j,nat,o1
      double precision al,be,cm,det,dtor,ga,pi,rot,rota,rr,rtod,t,t1,t2,
     & t3,tin,twopi,val,ver,xi,xo
      character atnam*4,card*80,forml*20,lun*500
      dimension cm(3),rot(3,3),rota(3,3),rr(6),tin(3,3),val(3),ver(3),
     & xi(3),xo(3)
      common/angkte/ pi,twopi,dtor,rtod
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      ilec=5
      iecr=6
      i1=1
      o1=2
      forml=' (12x,a4,14x,3f8.3) '
      write(iecr,'(a)') ' enter the input and the output pdb filenames:'
      read(ilec,'(a)') lun
      call compakt(lun,500)
      ind=index(lun,' ')-1
      open(unit=i1,file=lun(1:ind),form='formatted',status='old')
      rewind(unit=i1)
       do i=1,ind
      lun(i:i)=' '
       enddo
      call compakt(lun,500)
      ind=index(lun,' ')-1
      open(unit=o1,file=lun(1:ind),form='formatted',status='unknown')
      rewind(unit=o1)
       do i=1,3
      cm(i)=0.
        do j=1,3
      tin(i,j)=0.
        enddo
       enddo
      nat=0
10    read(i1,'(a)',end=20) card
      if(card(1:4).ne.'ATOM') go to 10
      read(card,fmt=forml) atnam,xi
      if(atnam(2:2).eq.'H'.or.atnam(2:2).eq.'D'.or.atnam(2:2).eq.'E'.or.
     . atnam(2:2).eq.'X'.or.atnam(1:1).eq.'H') go to 10
         if(nat.eq.0) then
       do i=1,3
      xo(i)=xi(i)
       enddo
         endif
      nat=nat+1
       do i=1,3
      cm(i)=cm(i)+xi(i)
        do j=1,3
      tin(i,j)=tin(i,j)+xi(i)*xi(j)
        enddo
       enddo
      go to 10
20    continue
       do i=1,3
      cm(i)=cm(i)/nat
       enddo
       do i=1,3
        do j=1,3
      tin(i,j)=tin(i,j)/nat-cm(i)*cm(j)
        enddo
       enddo
      write(iecr,'(a,3f15.5)') ' center-of-mass:',cm
      write(iecr,'(a,3(/3f10.2))') ' inertia tensor:',
     . ((tin(i,j),j=1,3),i=1,3)
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call rs(3,3,tin,val(1),rota,rr(1),rr(4),1,ierr)
      t=val(1)
      val(1)=val(3)
      val(3)=t
      call det3(rota,det)
      t=-sign(1.d0,det)
       do i=1,3
      ver(i)=0.
        do j=1,3
      ver(i)=ver(i)+(xo(j)-cm(j))*rota(j,i)
        enddo
       enddo
      t1=sign(1.d0,ver(1))
      t2=sign(1.d0,ver(2))
      t3=sign(1.d0,ver(3))
      t1=t*t2*t3
       do i=1,3
      rot(i,1)=rota(i,3)*t3
      rot(i,2)=rota(i,2)*t2
      rot(i,3)=rota(i,1)*t1
       enddo
      call rmx2e(rot,al,be,ga)
      call rmxe(180-ga,be,180-al,rota)
      write(iecr,'(a,i10,a,3f10.3)')
     . ' center of mass of the',nat,' input atoms =',cm
      write(iecr,'(a,3e15.5)') ' inertia moments of input coords. =',val
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(o1,'(a,6f10.3)') 'ROTTRA',al,be,ga,cm
      rewind(unit=i1)
30    read(i1,'(a)',end=40) card
         if(card(1:4).ne.'ATOM') then
      if(card(1:6).ne.'ROTTRA') write(o1,'(a)') card
      go to 30
         endif
      read(card,fmt=forml) atnam,xi
       do i=1,3
      xo(i)=0.
        do j=1,3
      xo(i)=xo(i)+rota(i,j)*(xi(j)-cm(j))
        enddo
       enddo
      write(card(31:54),'(3f8.3)') xo
      write(o1,'(a)') card
      go to 30
40    close(unit=i1)
      close(unit=o1)
      stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine compakt(card,nfield)
      implicit none
      integer i,j,lungo,nfield
      character card*1
      logical blanco,white
      external lungo
      dimension card(nfield)
      blanco=card(1).eq.' '
      j=0
       do i=1,lungo(card,nfield)
      white=card(i).eq.' '
         if((.not.white).or.(.not.blanco)) then
            if(white) then
      blanco=.true.
            else
      blanco=.false.
            endif
      j=j+1
      card(j)=card(i)
         endif
       enddo
       do i=j+1,lungo(card,nfield)
      card(i)=' '
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rs(nm,n,a,w,z,fv1,fv2,matz,ierr)
      implicit none
      integer ierr,matz,n,nm
      double precision a,fv1,fv2,w,z
      dimension a(nm,n),fv1(n),fv2(n),w(n),z(nm,n)
      if(n.le.nm) go to 10
      ierr=10*n
      go to 50
10    if(matz.ne.0) go to 20
      call tred1(nm,n,a,w,fv1,fv2)
      call tqlrat(n,w,fv2,ierr)
      go to 50
20    call tred2(nm,n,a,w,fv1,z)
      call tql2(nm,n,w,fv1,z,ierr)
50    return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine det3(r,d)
      implicit none
      double precision a,b,d,r
      dimension r(3,3)
      a=r(1,1)*r(2,2)*r(3,3)+r(1,2)*r(2,3)*r(3,1)+r(3,2)*r(2,1)*r(1,3)
      b=r(1,3)*r(2,2)*r(3,1)+r(2,1)*r(1,2)*r(3,3)+r(2,3)*r(3,2)*r(1,1)
      d=a-b
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
      write(6,'(a,3f10.3)') ' euler; alpha,beta,gamma',alpha,beta,gamma
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmxe(alpha,beta,gamma,rotx)
      implicit none
      double precision alf,alpha,bet,beta,cosa,cosb,cosg,dtor,gam,gamma,
     & pi,rotx,rtod,sina,sinb,sing,twopi
      dimension rotx(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      write(6,'(a,3f10.3)') ' euler matrix; alpha,beta,gamma',
     . alpha,beta,gamma
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
      subroutine tred2(nm,n,a,d,e,z)
      implicit none
      integer i,ii,j,jp1,k,l,n,nm
      double precision a,d,e,f,g,h,hh,scale,z
      dimension a(nm,n),d(n),e(n),z(nm,n)
      do 100 i=1,n
      do 80 j=i,n
      z(j,i)=a(j,i)
80    continue
      d(i)=a(n,i)
100   continue
      if(n.eq.1) go to 510
      do 300 ii=2,n
      i=n+2-ii
      l=i-1
      h=0.0d0
      scale=0.0d0
      if(l.lt.2) go to 130
      do 120 k=1,l
      scale=scale+abs(d(k))
120   continue
      if(scale.ne.0.0d0) go to 140
130   e(i)=d(l)
      do 135 j=1,l
      d(j)=z(l,j)
      z(i,j)=0.0d0
      z(j,i)=0.0d0
135   continue
      go to 290
140   do 150 k=1,l
      d(k)=d(k)/scale
      h=h+d(k)*d(k)
150   continue
      f=d(l)
      g=-sign(sqrt(h),f)
      e(i)=scale*g
      h=h-f*g
      d(l)=f-g
      do 170 j=1,l
      e(j)=0.0d0
170   continue
      do 240 j=1,l
      f=d(j)
      z(j,i)=f
      g=e(j)+z(j,j)*f
      jp1=j+1
      if(l.lt.jp1) go to 220
      do 200 k=jp1,l
      g=g+z(k,j)*d(k)
      e(k)=e(k)+z(k,j)*f
200   continue
220   e(j)=g
240   continue
      f=0.0d0
      do 245 j=1,l
      e(j)=e(j)/h
      f=f+e(j)*d(j)
245   continue
      hh=f/(h+h)
      do 250 j=1,l
      e(j)=e(j)-hh*d(j)
250   continue
      do 280 j=1,l
      f=d(j)
      g=e(j)
      do 260 k=j,l
      z(k,j)=z(k,j)-f*e(k)-g*d(k)
260   continue
      d(j)=z(l,j)
      z(i,j)=0.0d0
280   continue
290   d(i)=h
300   continue
      do 500 i=2,n
      l=i-1
      z(n,l)=z(l,l)
      z(l,l)=1.0d0
      h=d(i)
      if(h.eq.0.0d0) go to 380
      do 330 k=1,l
      d(k)=z(k,i)/h
330   continue
      do 370 j=1,l
      g=0.0d0
      do 340 k=1,l
      g=g+z(k,i)*z(k,j)
340   continue
      do 360 k=1,l
      z(k,j)=z(k,j)-g*d(k)
360   continue
370   continue
380   do 400 k=1,l
      z(k,i)=0.0d0
400   continue
500   continue
510   do 520 i=1,n
      d(i)=z(n,i)
      z(n,i)=0.0d0
520   continue
      z(n,n)=1.0d0
      e(1)=0.0d0
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tred1(nm,n,a,d,e,e2)
      implicit none
      integer i,ii,j,jp1,k,l,n,nm
      double precision a,d,e,e2,f,g,h,scale
      dimension a(nm,n),d(n),e(n),e2(n)
      do 100 i=1,n
      d(i)=a(n,i)
      a(n,i)=a(i,i)
100   continue
      do 300 ii=1,n
      i=n+1-ii
      l=i-1
      h=0.0d0
      scale=0.0d0
      if(l.lt.1) go to 130
      do 120 k=1,l
      scale=scale+abs(d(k))
120   continue
      if(scale.ne.0.0d0) go to 140
      do 125 j=1,l
      d(j)=a(l,j)
      a(l,j)=a(i,j)
      a(i,j)=0.0d0
125   continue
130   e(i)=0.0d0
      e2(i)=0.0d0
      go to 300
140   do 150 k=1,l
      d(k)=d(k)/scale
      h=h+d(k)*d(k)
150   continue
      e2(i)=scale*scale*h
      f=d(l)
      g=-sign(sqrt(h),f)
      e(i)=scale*g
      h=h-f*g
      d(l)=f-g
      if(l.eq.1) go to 285
      do 170 j=1,l
      e(j)=0.0d0
170   continue
      do 240 j=1,l
      f=d(j)
      g=e(j)+a(j,j)*f
      jp1=j+1
      if(l.lt.jp1) go to 220
      do 200 k=jp1,l
      g=g+a(k,j)*d(k)
      e(k)=e(k)+a(k,j)*f
200   continue
220   e(j)=g
240   continue
      f=0.0d0
      do 245 j=1,l
      e(j)=e(j)/h
      f=f+e(j)*d(j)
245   continue
      h=f/(h+h)
      do 250 j=1,l
      e(j)=e(j)-h*d(j)
250   continue
      do 280 j=1,l
      f=d(j)
      g=e(j)
      do 260 k=j,l
      a(k,j)=a(k,j)-f*e(k)-g*d(k)
260   continue
280   continue
285   do 290 j=1,l
      f=d(j)
      d(j)=a(l,j)
      a(l,j)=a(i,j)
      a(i,j)=f*scale
290   continue
300   continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function lungo(card,nfield)
      implicit none
      integer i,lungo,nfield
      character card*1
      dimension card(nfield)
       do i=nfield,1,-1
      if(card(i).ne.' ') go to 10
       enddo
10    lungo=i
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tql2(nm,n,d,e,z,ierr)
      implicit none
      integer i,ierr,ii,j,k,l,l1,l2,m,mml,n,nm
      double precision c,c2,c3,d,dl1,e,el1,f,g,h,p,pythag,r,s,s2,tst1,
     & tst2,z
      external pythag
      dimension d(n),e(n),z(nm,n)
cjn uninitialized values
      data c3,s2/0.,0./
      ierr=0
      if(n.eq.1) go to 1001
      do 100 i=2,n
      e(i-1)=e(i)
100   continue
      f=0.0d0
      tst1=0.0d0
      e(n)=0.0d0
      do 240 l=1,n
      j=0
      h=abs(d(l))+abs(e(l))
      if(tst1.lt.h) tst1=h
      do 110 m=l,n
      tst2=tst1+abs(e(m))
      if(tst2.eq.tst1) go to 120
110   continue
120   if(m.eq.l) go to 220
130   if(j.eq.30) go to 1000
      j=j+1
      l1=l+1
      l2=l1+1
      g=d(l)
      p=(d(l1)-g)/(2.0d0*e(l))
      r=pythag(p,1.0d0)
      d(l)=e(l)/(p+sign(r,p))
      d(l1)=e(l)*(p+sign(r,p))
      dl1=d(l1)
      h=g-d(l)
      if(l2.gt.n) go to 145
      do 140 i=l2,n
      d(i)=d(i)-h
140   continue
145   f=f+h
      p=d(m)
      c=1.0d0
      c2=c
      el1=e(l1)
      s=0.0d0
      mml=m-l
      do 200 ii=1,mml
      c3=c2
      c2=c
      s2=s
      i=m-ii
      g=c*e(i)
      h=c*p
      r=pythag(p,e(i))
      e(i+1)=s*r
      s=e(i)/r
      c=p/r
      p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 180 k=1,n
      h=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*h
      z(k,i)=c*z(k,i)-s*h
180   continue
200   continue
      p=-s*s2*c3*el1*e(l)/dl1
      e(l)=s*p
      d(l)=c*p
      tst2=tst1+abs(e(l))
      if(tst2.gt.tst1) go to 130
220   d(l)=d(l)+f
240   continue
      do 300 ii=2,n
      i=ii-1
      k=i
      p=d(i)
      do 260 j=ii,n
      if(d(j).ge.p) go to 260
      k=j
      p=d(j)
260   continue
      if(k.eq.i) go to 300
      d(k)=d(i)
      d(i)=p
      do 280 j=1,n
      p=z(j,i)
      z(j,i)=z(j,k)
      z(j,k)=p
280   continue
300   continue
      go to 1001
1000  ierr=l
1001  return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tqlrat(n,d,e2,ierr)
      implicit none
      integer i,ierr,ii,j,l,l1,m,mml,n
      double precision b,c,d,e2,epslon,f,g,h,p,pythag,r,s,t
      external epslon,pythag
      dimension d(n),e2(n)
cjn uninitialized values
      data b,c/0.,0./
      ierr=0
      if(n.eq.1) go to 1001
      do 100 i=2,n
      e2(i-1)=e2(i)
100   continue
      f=0.0d0
      t=0.0d0
      e2(n)=0.0d0
      do 290 l=1,n
      j=0
      h=abs(d(l))+sqrt(e2(l))
      if(t.gt.h) go to 105
      t=h
      b=epslon(t)
      c=b*b
105   do 110 m=l,n
      if(e2(m).le.c) go to 120
110   continue
120   if(m.eq.l) go to 210
130   if(j.eq.30) go to 1000
      j=j+1
      l1=l+1
      s=sqrt(e2(l))
      g=d(l)
      p=(d(l1)-g)/(2.0d0*s)
      r=pythag(p,1.0d0)
      d(l)=s/(p+sign(r,p))
      h=g-d(l)
      do 140 i=l1,n
      d(i)=d(i)-h
140   continue
      f=f+h
      g=d(m)
      if(g.eq.0.0d0) g=b
      h=g
      s=0.0d0
      mml=m-l
      do 200 ii=1,mml
      i=m-ii
      p=g*h
      r=p+e2(i)
      e2(i+1)=s*r
      s=e2(i)/r
      d(i+1)=h+s*(h+d(i))
      g=d(i)-e2(i)/g
      if(g.eq.0.0d0) g=b
      h=g*p/r
200   continue
      e2(l)=s*g
      d(l)=h
      if(h.eq.0.0d0) go to 210
      if(abs(e2(l)).le.abs(c/h)) go to 210
      e2(l)=h*e2(l)
      if(e2(l).ne.0.0d0) go to 130
210   p=d(l)+f
      if(l.eq.1) go to 250
      do 230 ii=2,l
      i=l+2-ii
      if(p.ge.d(i-1)) go to 270
      d(i)=d(i-1)
230   continue
250   i=1
270   d(i)=p
290   continue
      go to 1001
1000  ierr=l
1001  return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function pythag(a,b)
      implicit none
      double precision a,b,p,pythag,r,s,t,u
      p=max(abs(a),abs(b))
      if(p.eq.0.0d0) go to 20
      r=(min(abs(a),abs(b))/p)**2
10    continue
      t=4.0d0+r
      if(t.eq.4.0d0) go to 20
      s=r/t
      u=1.0d0+2.0d0*s
      p=u*p
      r=(s/u)**2*r
      go to 10
20    pythag=p
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function epslon(x)
      implicit none
      double precision a,b,c,eps,epslon,x
      a=4.0d0/3.0d0
10    b=a-1.0d0
      c=b+b+b
      eps=abs(c-1.0d0)
      if(eps.eq.0.0d0) go to 10
      epslon=eps*abs(x)
      return
      end
