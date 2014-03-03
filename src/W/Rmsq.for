c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program rmsq
      implicit none
      integer i1,i2,iecr,ilec,iout,kprt,kprtt,length,ma,n,na,na1,nprt
      double precision dtor,lxy,pi,qxy,rr,rtod,twopi
      character aa*4,cbid1*1,lun*80
      external length
      parameter(ma=1000000)
      dimension aa(2*ma),lxy(6),qxy(6,6),rr(6*ma)
      common/angkte/ pi,twopi,dtor,rtod
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      common/oprt/ kprtt(10)
      data i1,i2/98,99/
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      ilec=5
      iecr=6
      iout=9
      nprt=10
       do n=1,nprt
      kprtt(n)=0
       enddo
c     write(iecr,'(a)') ' enter printing flags [N  1 0 ... ]:'
c     read(ilec,'(a)') lun
c        if(length(lun).ne.0) then
c     read(lun,*) nprt,(kprtt(n),n=1,nprt)
c        endif
10    write(iecr,'(/a)') ' enter the target-coord [y] filename:'
      read(ilec,'(a)',end=20) lun
      if(lun(1:3).eq.'end'.or.lun(1:3).eq.'   ') go to 20
      open(unit=i1,file=lun,form='formatted',status='old')
      write(iecr,'(/a,a)') ' target file:',lun(1:40)
      write(iecr,'(a)') ' enter the search-coord [x] filename:'
      read(ilec,'(a)') lun
      open(unit=i2,file=lun,form='formatted',status='old')
      write(iecr,'(a,a)') ' search file:',lun(1:40)
      call cordin(i2,ma,aa(1),rr(1),na)
      close(unit=i2)
      call cordin(i1,ma,aa(1+na),rr(1+3*na),na1)
      close(unit=i1)
         if(na1.ne.na) then
      write(iecr,2010) na1,na
      go to 10
         endif
      call quanta(na,rr(1),rr(1+3*na),qxy,lxy)
      call lsqfit(lxy,qxy)
      go to 10
20    stop
2010  format(/' >> warning << different nb. of atoms:',2i10)
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
      subroutine cordin(i1,ma,axyz,cxyz,na)
      implicit none
      integer i,i1,iecr,ilec,iout,kprt,ma,na
      double precision biso,cxyz,xi
      character atnam*4,axyz*4,cbid1*1
      logical first,last
      dimension axyz(ma),cxyz(3,ma),xi(3)
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      first=.true.
      last=.false.
      na=0
10    call lecatc(i1,atnam,xi,biso,first,last)
      if(last) go to 20
      na=na+1
      if(na.gt.ma) go to 901
      axyz(na)=atnam
       do i=1,3
      cxyz(i,na)=xi(i)
       enddo
      go to 10
20    return
901   write(iecr,3010)
      stop
3010  format('stop >> cordin << increase ma')
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine quanta(na,cxyz1,cxyz2,qxy,lxy)
      implicit none
      integer i,iecr,ilec,iout,j,kprt,kprtt,n,na
      double precision cxyz1,cxyz2,lxy,qd,qr,qt,qxy,sl,sq,xy
      character cbid1*1
      dimension cxyz1(3,na),cxyz2(3,na),lxy(6),qxy(6,6),sl(6),sq(6,6),
     & xy(6)
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      common/oprt/ kprtt(10)
      kprt=kprtt(2)
       do i=1,6
        do j=1,6
      sq(i,j)=0.d0
        enddo
      sl(i)=0.d0
       enddo
       do n=1,na
        do i=1,3
      xy(i)=cxyz1(i,n)
      xy(3+i)=cxyz2(i,n)
        enddo
        do i=1,6
         do j=1,6
      sq(i,j)=sq(i,j)+xy(i)*xy(j)
         enddo
      sl(i)=sl(i)+xy(i)
        enddo
       enddo
       do i=1,6
        do j=1,6
      qxy(i,j)=(sq(i,j)/na)-(sl(i)/na)*(sl(j)/na)
        enddo
      lxy(i)=sl(i)/na
       enddo
         if(kprt.ne.0) then
      write(iecr,'(/a,2(/3f15.5))') ' <x> and <y>',(lxy(i),i=1,6)
      write(iecr,'(/a,3(/3f15.5))') ' <x*x>',((qxy(i,j),j=1,3),i=1,3)
      write(iecr,'(/a,3(/3f15.5))') ' <y*y>',((qxy(i,j),j=4,6),i=4,6)
      write(iecr,'(/a,3(/3f15.5))') ' <x*y>',((qxy(i,j),j=4,6),i=1,3)
         endif
      qr=0.d0
      qt=0.d0
        do i=1,3
      qr=qr+qxy(i,i)+qxy(3+i,3+i)-qxy(i,3+i)-qxy(3+i,i)
      qt=qt+(lxy(i)-lxy(3+i))**2
        enddo
      qd=qr+qt
      if(qd.gt.0.d0) qd=sqrt(qd)
      if(qr.gt.0.d0) qr=sqrt(qr)
      if(qt.gt.0.d0) qt=sqrt(qt)
      write(iecr,'(/a,f10.5)') ' rms shift = sqrt(<[y-x]^2>):',qd
      write(iecr,'(a,f10.5)') ' rms shift centered coordinates:',qr
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lsqfit(lxy,qxy)
      implicit none
      integer i,iecr,ierr,ilec,iout,j,k,kprt,kprtt
      double precision al,ang,be,det,ga,lxy,mat,qmx,qxy,rms,rmx,rr,som,
     & sum,t,trx,val
      character cbid1*1
      dimension lxy(6),mat(3,3),qmx(3,3),qxy(6,6),rmx(3,3),rr(6),trx(3),
     & val(3)
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      common/oprt/ kprtt(10)
      kprt=kprtt(1)
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+qxy(i,3+k)*qxy(3+k,j)
         enddo
      qmx(i,j)=sum
        enddo
       enddo
      call rs(3,3,qmx,val,mat,rr(1),rr(4),1,ierr)
      call det3(mat,det)
      t=sign(1.d0,det)
      sum=0.d0
       do i=1,3
      mat(i,1)=mat(i,1)*t
      sum=sum+qxy(i,i)+qxy(3+i,3+i)-2*sqrt(val(i))
       enddo
      if(sum.gt.0.d0) sum=sqrt(sum)
      rms=sum
      call det3(mat,det)
         if(kprt.ne.0) then
      write(iecr,'(/a,e15.5)') ' determinant',det
      write(iecr,'(/a,3e15.5)') ' eigenvalues',(val(i),i=1,3)
         endif
       do i=1,3
      val(i)=sqrt(1/val(i))
       enddo
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+mat(i,k)*val(k)*mat(j,k)
         enddo
      qmx(i,j)=sum
        enddo
       enddo
       do i=1,3
      som=0.d0
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+qxy(3+i,k)*qmx(k,j)
         enddo
      som=som+sum*lxy(j)
      rmx(i,j)=sum
        enddo
      trx(i)=lxy(3+i)-som
       enddo
      call det3(rmx,det)
      write(iecr,'(/a,f15.5)')
     . ' intrinsic rms = sqrt(<[y-(M.x+T)]^2>):',rms
      write(iecr,'(a,3(/3f15.5))') ' matrix M',((rmx(i,j),j=1,3),i=1,3)
      write(iecr,'(a/3f15.3)') ' translation T',(trx(i),i=1,3)
         if(det.lt.0.) then
       do i=1,3
        do j=1,3
      rmx(i,j)=-rmx(i,j)
        enddo
       enddo
         endif
      kprt=0
      call lsqmx(rmx,mat)
      call rmx2e(mat,al,be,ga)
      call rmx2p(mat,ang,val)
      write(iecr,'(a/f15.5)') ' determinant(M)',det
         if(det.lt.0.) then
      write(iecr,'(a/f15.3,3f15.5)') ' rotation R = (-1) M',ang,val
         else
      write(iecr,'(a/f15.3,3f15.5)') ' rotation R = M',ang,val
         endif
      write(iecr,'(/a,6f9.3)') ' transformation (R,T):',al,be,ga,trx
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lecatc(i1,carat,xi,biso,first,last)
      implicit none
      integer i,i1,j,nort
      double precision a,alpha,b,beta,biso,c,dtor,gamma,pi,rbid1,rf,ro,
     & rtod,twopi,xi,xo
      character atnam*4,carat*4,card*80,forml*80,tcar*80,tipo*4
      logical first,last
      dimension xi(3),xo(3)
      common/angkte/ pi,twopi,dtor,rtod
      common/cell/ a,b,c,alpha,beta,gamma,rbid1(16)
      common/ortm/ ro(3,3),rf(3,3),nort
      save forml,tipo
      carat='    '
         if(first) then
      first=.false.
       do i=1,80
      forml(i:i)=' '
       enddo
      forml(1:26)='(12X,A4,14X,3F8.3,6X,F6.2)'
      rewind(unit=i1)
      read(i1,fmt='(a)') card
            if(card(1:4).eq.'CELL') then
      tipo='frac'
      tcar(1:76)=card(5:80)
c lecturainterna
      read(tcar,*) a,b,c,alpha,beta,gamma
      call celda
      nort=1
      call ortho
            else
      tipo='orto'
      go to 30
            endif
         endif
20    read(i1,fmt='(a)',end=70) card
      if(card(1:8).eq.'END-DATA') go to 70
30    if(card(1:4).ne.'ATOM') go to 20
      read(card,fmt=forml) atnam,xi,biso
c     if(atnam.ne.' CA '
c    . .and.atnam.ne.' C  '.and.atnam.ne.' N  '.and.atnam.ne.' O  '
c    . ) go to 20
      if(atnam(2:2).eq.'H'.or.atnam(2:2).eq.'D'.or.atnam(2:2).eq.'E'.or.
     . atnam(2:2).eq.'X'.or.atnam(1:1).eq.'H') go to 20
      carat(1:2)=atnam(1:2)
         if(tipo.ne.'orto') then
      do 50 i=1,3
      xo(i)=0.
      do 40 j=1,3
      xo(i)=xo(i)+ro(i,j)*xi(j)
40    continue
50    continue
      do 60 i=1,3
      xi(i)=xo(i)
60    continue
      biso=biso*8*pi**2
         endif
      return
70    last=.true.
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
      subroutine lsqmx(r1,r2)
      implicit none
      integer i,ierr,j,k
      double precision mat,r1,r2,r3,rr,sum,val
      dimension mat(3,3),r1(3,3),r2(3,3),r3(3,3),rr(6),val(3)
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+r1(i,k)*r1(j,k)
         enddo
      r3(i,j)=sum
        enddo
       enddo
      call rs(3,3,r3,val,mat,rr(1),rr(4),1,ierr)
       do i=1,3
      val(i)=sqrt(1/val(i))
       enddo
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+mat(i,k)*val(k)*mat(j,k)
         enddo
      r3(i,j)=sum
        enddo
       enddo
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+r3(i,k)*r1(k,j)
         enddo
      r2(i,j)=sum
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmx2e(rotx,alpha,beta,gamma)
      implicit none
      integer iecr,ilec,iout,kprt
      double precision alpha,ang,beta,cang,dtor,fuzz,gamma,pi,rotx,rtod,
     & sang,twopi
      character cbid1*1
      dimension rotx(3,3)
      common/angkte/ pi,twopi,dtor,rtod
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
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
      if(kprt.ne.0) write(iecr,'(a,3f10.3)')
     . ' euler; alpha,beta,gamma',alpha,beta,gamma
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmx2p(rotx,xhi,vn)
      implicit none
      integer i,iecr,ilec,iout,kprt
      double precision ang,cang,dtor,fuzz,norm,pi,rotx,rtod,sang,scale,
     & twopi,vecx,vn,xhi
      character cbid1*1
      dimension rotx(3,3),vecx(3),vn(3)
      common/angkte/ pi,twopi,dtor,rtod
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      fuzz=0.001
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
      if(kprt.ne.0) write(iecr,'(a,f10.3,3f10.6,f10.6)')
     . ' polar; xhi,versor,norm',xhi,vn,norm
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
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine celda
      implicit none
      integer iecr,ilec,iout,kprt
      double precision a,a11,a21,a22,a31,a32,a33,alpha,ast,b,beta,bst,c,
     & ca,cast,cb,cbst,cg,cgst,cst,dtor,gamma,pi,r,rtod,sa,sast,sb,sbst,
     & sg,sgst,twopi,vol
      character cbid1*1
      common/angkte/ pi,twopi,dtor,rtod
      common/cell/ a,b,c,alpha,beta,gamma,ca,cb,cg,sa,sb,sg,ast,bst,cst,
     & cast,cbst,cgst,sast,sbst,sgst,vol
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      common/star/ a11,a21,a22,a31,a32,a33
      r=alpha*dtor
      ca=cos(r)
      sa=sin(r)
      r=beta*dtor
      cb=cos(r)
      sb=sin(r)
      r=gamma*dtor
      cg=cos(r)
      sg=sin(r)
      vol=a*b*c*sqrt(1.+2.*ca*cb*cg-ca*ca-cb*cb-cg*cg)
      cast=(cb*cg-ca)/(sb*sg)
      cbst=(cg*ca-cb)/(sg*sa)
      cgst=(ca*cb-cg)/(sa*sb)
      sast=sqrt(1.-cast*cast)
      sbst=sqrt(1.-cbst*cbst)
      sgst=sqrt(1.-cgst*cgst)
      ast=1/(a*sb*sgst)
      bst=1/(b*sg*sast)
      cst=1/(c*sa*sbst)
      a11=ast*sbst*sg
      a21=-ast*sbst*cg
      a22=bst*sast
      a31=ast*cbst
      a32=bst*cast
      a33=cst
      if(kprt.ne.0) write(iecr,2010) a,b,c,alpha,beta,gamma,vol
      return
2010  format(/' cell parameters =',6f10.3/' volume =',e15.5)
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ortho
      implicit none
      integer i,iecr,ij,ilec,iout,ip,iq,is,isym,it,itr,j,kprt,mss,n0,
     & neq,nort,nsym,nts
      double precision a,alpha,arg,ast,b,beta,bst,c,c0,c1,c2,c3,ca,cast,
     & cb,cbst,cg,cgst,cst,det,fuzz,gamma,gd,gi,p,p1,p2,p3,proy,rbid1,
     & rf,ro,rs,rx,sa,sast,sb,sbst,sg,sgst,unor,v1,v2,v3,vol
      character cbid1*1
      dimension c0(3,48),c1(3),c2(3),c3(3),gd(3,3),gi(3,3),n0(48),rs(9),
     & rx(3,3),v1(3),v2(3),v3(3)
      common/cell/ a,b,c,alpha,beta,gamma,ca,cb,cg,sa,sb,sg,ast,bst,cst,
     & cast,cbst,cgst,sast,sbst,sgst,vol
      common/ioprg/ ilec,iecr,kprt,iout,cbid1(81)
      common/ortm/ ro(3,3),rf(3,3),nort
      common/simt/ rbid1(153),mss(9,48),neq,nts
cjn uninitialized values
      data iq,proy/0,0./
      fuzz=0.001
      if(nort.eq.0.and.neq.eq.1) nort=1
      if(nort.ne.0) go to 190
      gd(1,1)=a*a
      gd(1,2)=a*b*cg
      gd(1,3)=a*c*cb
      gd(2,2)=b*b
      gd(2,3)=b*c*ca
      gd(3,3)=c*c
      gd(2,1)=gd(1,2)
      gd(3,1)=gd(1,3)
      gd(3,2)=gd(2,3)
      gi(1,1)=ast*ast
      gi(1,2)=ast*bst*cgst
      gi(1,3)=ast*cst*cbst
      gi(2,2)=bst*bst
      gi(2,3)=bst*cst*cast
      gi(3,3)=cst*cst
      gi(2,1)=gi(1,2)
      gi(3,1)=gi(1,3)
      gi(3,2)=gi(2,3)
      nsym=1
      do 60 is=2,neq
      n0(is)=0
      do 10 ij=1,9
      rs(ij)=mss(ij,is)
10    continue
      call det3(rs,det)
      itr=mss(1,is)+mss(5,is)+mss(9,is)
      if(nint(det).ne.1.or.itr.gt.2) go to 60
      isym=itr+3+abs(itr)/2
         if(isym.ge.nsym) then
      nsym=isym
      call pro2mx(rx,gd,rs)
            if(itr.gt.-1) then
      unor=sqrt(3.-itr*(itr-2.))*vol
      c3(1)=(rx(3,2)-rx(2,3))/unor
      c3(2)=(rx(1,3)-rx(3,1))/unor
      c3(3)=(rx(2,1)-rx(1,2))/unor
            else
      do 20 i=1,3
      arg=(rx(i,i)+gd(i,i))/2
      if(arg.le.0.) arg=0.
      v3(i)=sqrt(arg)
20    continue
               if(v3(3).ne.0.) then
      v3(1)=sign(v3(1),rx(1,3)+gd(1,3))
      v3(2)=sign(v3(2),rx(2,3)+gd(2,3))
               else if(v3(2).ne.0.) then
      v3(1)=sign(v3(1),rx(1,2)+gd(1,2))
               endif
      do 40 i=1,3
      c3(i)=0.
      do 30 j=1,3
      c3(i)=c3(i)+gi(j,i)*v3(j)
30    continue
40    continue
            endif
      n0(is)=isym
      do 50 i=1,3
      c0(i,is)=c3(i)
50    continue
      p3=abs(c3(3))/cst
      p2=abs(c3(2))/bst
      p1=abs(c3(1))/ast
      proy=max(p3,p2,p1)
         endif
60    continue
         if(nsym.eq.1) then
      nort=1
      go to 190
         endif
      ip=0
      do 70 is=2,neq
      if(n0(is).ne.nsym) go to 70
      p3=abs(c0(3,is))/cst
      p2=abs(c0(2,is))/bst
      p1=abs(c0(1,is))/ast
      p=max(p3,p2,p1)
      if(proy-p.gt.fuzz) go to 70
         if(p-p3.le.fuzz) then
      it=3
         else if(p-p2.le.fuzz) then
      it=2
         else if(p-p1.le.fuzz) then
      it=1
         else
      go to 901
         endif
      if(ip.eq.3.and.it.ne.3) go to 70
      proy=p
      ip=it
      iq=is
70    continue
      do 80 i=1,3
      c3(i)=c0(i,iq)
      if(c0(ip,iq).lt.0.) c3(i)=-c3(i)
80    continue
      do 100 i=1,3
      v3(i)=0.
      do 90 j=1,3
      v3(i)=v3(i)+gd(i,j)*c3(j)
90    continue
100   continue
      do 110 i=1,3
      c1(i)=0.
      c2(i)=0.
      v1(i)=0.
110   continue
      p1=abs(v3(1))/a
      p2=abs(v3(2))/b
      p3=abs(v3(3))/c
      p=min(p1,p2,p3)
         if(p1-p.le.fuzz) then
      c1(1)=1/a
      proy=v3(1)/a
         else if(p2-p.le.fuzz) then
      c1(2)=1/b
      proy=v3(2)/b
         else if(p3-p.le.fuzz) then
      c1(3)=1/c
      proy=v3(3)/c
         else
      go to 901
         endif
      do 120 i=1,3
      c1(i)=c1(i)-proy*c3(i)
120   continue
      unor=0.
      do 140 i=1,3
      do 130 j=1,3
      v1(i)=v1(i)+gd(i,j)*c1(j)
130   continue
      unor=unor+c1(i)*v1(i)
140   continue
      unor=sqrt(unor)
      do 150 i=1,3
      c1(i)=c1(i)/unor
      v1(i)=v1(i)/unor
150   continue
      v2(1)=(c3(2)*c1(3)-c3(3)*c1(2))*vol
      v2(2)=(c3(3)*c1(1)-c3(1)*c1(3))*vol
      v2(3)=(c3(1)*c1(2)-c3(2)*c1(1))*vol
      do 170 i=1,3
      do 160 j=1,3
      c2(i)=c2(i)+gi(j,i)*v2(j)
160   continue
170   continue
      do 180 i=1,3
      ro(1,i)=v1(i)
      ro(2,i)=v2(i)
      ro(3,i)=v3(i)
      rf(i,1)=c1(i)
      rf(i,2)=c2(i)
      rf(i,3)=c3(i)
180   continue
      go to 220
190   do 210 j=1,3
      do 200 i=1,3
      ro(i,j)=0.
      rf(i,j)=0.
200   continue
210   continue
         if(nort.eq.1) then
      ro(1,1)=a
      ro(1,2)=b*cg
      ro(1,3)=c*cb
      ro(2,2)=b*sg
      ro(2,3)=-c*sb*cast
      ro(3,3)=c*sa*sbst
      rf(1,1)=ast*sbst*sg
      rf(1,2)=-ast*sbst*cg
      rf(1,3)=ast*cbst
      rf(2,2)=bst*sast
      rf(2,3)=bst*cast
      rf(3,3)=cst
         else if(nort.eq.2) then
      ro(1,2)=b
      ro(1,3)=c*ca
      ro(1,1)=a*cg
      ro(2,3)=c*sa
      ro(2,1)=-a*sg*cbst
      ro(3,1)=a*sb*sgst
      rf(2,1)=bst*sgst*sa
      rf(2,2)=-bst*sgst*ca
      rf(2,3)=bst*cgst
      rf(3,2)=cst*sbst
      rf(3,3)=cst*cbst
      rf(1,3)=ast
         else if(nort.eq.3) then
      ro(1,3)=c
      ro(1,1)=a*cb
      ro(1,2)=b*ca
      ro(2,1)=a*sb
      ro(2,2)=-b*sa*cgst
      ro(3,2)=b*sg*sast
      rf(3,1)=cst*sast*sb
      rf(3,2)=-cst*sast*cb
      rf(3,3)=cst*cast
      rf(1,2)=ast*sgst
      rf(1,3)=ast*cgst
      rf(2,3)=bst
         endif
220      if(kprt.ne.0) then
      write(iecr,2010) nort
      write(iecr,2020) ((ro(i,j),j=1,3),i=1,3),((rf(i,j),j=1,3),i=1,3)
         endif
      return
901   write(iecr,3010)
      write(iout,3010)
      stop
2010  format(/' transformation of coordinates (direct or reciprocal)'/
     . ' nort=0, orthog axes have z along highest symmetry axis'/
     . ' nort=1, orthog axes have x along a, z along cstar'/
     . ' nort=2, orthog axes have x along b, z along astar'/
     . ' nort=3, orthog axes have x along c, z along bstar'/
     . ' here nort is',i5)
2020  format(' orthogonalizing matrix for xyz'/3(/3f10.3)/
     . /' orthogonalizing matrix for hkl'/3(/3f10.6))
3010  format('stop >> ortho << increase fuzz')
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
      subroutine pro2mx(ro,r2,r1)
      implicit none
      integer i,j,k
      double precision add,d1,d2,r1,r2,ro
      dimension d1(3,3),d2(3,3),r1(3,3),r2(3,3),ro(3,3)
       do i=1,3
        do j=1,3
      d1(i,j)=r1(i,j)
      d2(i,j)=r2(i,j)
        enddo
       enddo
      do 30 i=1,3
      do 20 j=1,3
      add=0.d0
      do 10 k=1,3
      add=add+d2(i,k)*d1(k,j)
10    continue
      ro(i,j)=add
20    continue
30    continue
      return
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
