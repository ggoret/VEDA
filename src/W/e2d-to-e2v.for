c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program e2d
      implicit none
      integer i1,io1,length,mr,mx,my,mz,n,nx,ny,nz,o1,px,py,pz,xl,xlow,
     & xu,xup,yl,ylow,yu,yup,zl,zlow,zu,zup
      real a,alpha,b,beta,c,dx,dy,dz,gamma,magno,rr,scale,ss
      character card*80,lun*80
      external length
      parameter(mr=1000000)
      dimension rr(mr),ss(mr)
      data i1,o1,io1/1,2,3/
      write(6,'(/a)') ' enter the input ezd-map filename'
      read(5,'(a)') lun
      open(unit=i1,file=lun,form='formatted',status='old')
      rewind(unit=i1)
10    read(i1,'(a)') card
      if(card(1:4).ne.'CELL') go to 10
      read(card(5:80),*) a,b,c,alpha,beta,gamma
      read(i1,'(a)') card
      read(card(7:80),*) xlow,ylow,zlow
      read(i1,'(a)') card
      read(card(7:80),*) mx,my,mz
      read(i1,'(a)') card
      read(card(5:80),*) nx,ny,nz
      read(i1,'(a)') card
      read(card(6:80),*) scale
      read(i1,*)
      dx=a/nx
      dy=b/ny
      dz=c/nz
      xup=xlow+mx-1
      yup=ylow+my-1
      zup=zlow+mz-1
      write(6,'(/a/6f10.2/3i5)') ' input parameters:',
     . a,b,c,alpha,beta,gamma,nx,ny,nz
      write(6,'(a,6i5/a,3i5)') ' input indices [x],[y],[z]:',
     . xlow,xup,ylow,yup,zlow,zup,' extent:',mx,my,mz
         if(mx*my.gt.mr) then
      write(6,'(/a,i10)') ' STOP : e2d : set mr =',mx*my
      stop
         endif
      px=mx
      py=my
      pz=mz
      call control(px)
      call control(py)
      call control(pz)
      write(6,'(/a)') ' enter the magnification'
      read(5,*) magno
      dx=dx*magno
      dy=dy*magno
      dz=dz*magno
         if(px.ne.mx.or.py.ne.my.or.pz.ne.mz.or.magno.ne.1.0) then
      write(6,'(/a,3i5)') ' extent and|or magnification changed'
      xl=xlow+(mx-px)/2
      xu=xup-(mx-px+1)/2
      yl=ylow+(my-py)/2
      yu=yup-(my-py+1)/2
      zl=zlow+(mz-pz)/2
      zu=zup-(mz-pz+1)/2
      mx=xu-xl+1
      my=yu-yl+1
      mz=zu-zl+1
      nx=mx
      ny=my
      nz=mz
      a=mx*dx
      b=my*dy
      c=mz*dz
       do n=1,80
      card(n:n)=' '
       enddo
      write(card,'(6f12.3)') a,b,c,alpha,beta,gamma
      read(card,*) a,b,c,alpha,beta,gamma
      write(6,'(/a)') ' output ezd-map filename: d/emap.ezd'
       do n=1,80
      lun(n:n)=' '
       enddo
      lun(1:10)='d/emap.ezd'
      open(unit=io1,file=lun,form='formatted',status='unknown')
      rewind(unit=io1)
      write(io1,'(a)') 'EZD_MAP suitable for Ten Eyck FFT'
      write(card,'(a)') '! created by e2d'
      call compact(card)
      write(io1,'(a)') card(1:length(card))
      write(card,'(a,6f10.3)') 'CELL',a,b,c,alpha,beta,gamma
      call compact(card)
      write(io1,'(a)') card(1:length(card))
      write(card,'(a,3i5)') 'ORIGIN',xl,yl,zl
      call compact(card)
      write(io1,'(a)') card(1:length(card))
      write(card,'(a,3i5)') 'EXTENT',mx,my,mz
      call compact(card)
      write(io1,'(a)') card(1:length(card))
      write(card,'(a,3i5)') 'GRID',nx,ny,nz
      call compact(card)
      write(io1,'(a)') card(1:length(card))
      write(card,'(a,f15.3)') 'SCALE',scale
      call compact(card)
      write(io1,'(a)') card(1:length(card))
      write(io1,'(a)') 'MAP'
      call shrink(i1,io1,xlow,xup,ylow,yup,zlow,zup,xl,xu,yl,yu,zl,zu,
     . rr(1),ss(1))
      write(io1,'(a)') 'END'
      i1=io1
      rewind(unit=i1)
       do n=1,8
      read(i1,'(a)') card
       enddo
         else
      xl=xlow
      xu=xup
      yl=ylow
      yu=yup
      zl=zlow
      zu=zup
      nx=mx
      ny=my
      nz=mz
      a=mx*dx
      b=my*dy
      c=mz*dz
       do n=1,80
      card(n:n)=' '
       enddo
      write(card,'(6f12.3)') a,b,c,alpha,beta,gamma
      read(card,*) a,b,c,alpha,beta,gamma
         endif
      write(6,'(/a/6f10.2/3i5)') ' output parameters:',
     . a,b,c,alpha,beta,gamma,nx,ny,nz
      write(6,'(a,6i5/a,3i5)') ' output indices [x],[y],[z]:',
     . xl,xu,yl,yu,zl,zu,' extent:',mx,my,mz
      write(6,'(a)') ' enter the output map filename'
      read(5,'(a)') lun
      open(unit=o1,file=lun,form='unformatted',status='unknown')
      rewind(unit=o1)
      write(o1) mx,my,mz,a,b,c,alpha,beta,gamma
      call leesc(i1,o1,mx,my,mz,rr(1),scale)
      write(6,'(6f10.3,a)') a,b,c,alpha,beta,gamma,'     CELL'
      stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine control(n)
      implicit none
      integer i,maxs,mcd,mp1,n,ndiv,p1
      parameter(mp1=266,maxs=1000)
      dimension p1(mp1)
      data p1/
     &  46,  58,  62,  74,  82,  86,  92,  94, 106, 116, 118, 122,
     & 124, 134, 138, 142, 146, 148, 158, 164, 166, 172, 174, 178,
     & 184, 186, 188, 194, 202, 206, 212, 214, 218, 222, 226, 230,
     & 232, 236, 244, 246, 248, 254, 258, 262, 268, 274, 276, 278,
     & 282, 284, 290, 292, 296, 298, 302, 310, 314, 316, 318, 322,
     & 326, 328, 332, 334, 344, 346, 348, 354, 356, 358, 362, 366,
     & 368, 370, 372, 376, 382, 386, 388, 394, 398, 402, 404, 406,
     & 410, 412, 414, 422, 424, 426, 428, 430, 434, 436, 438, 444,
     & 446, 452, 454, 458, 460, 464, 466, 470, 472, 474, 478, 482,
     & 488, 492, 496, 498, 502, 506, 508, 514, 516, 518, 522, 524,
     & 526, 530, 534, 536, 538, 542, 548, 552, 554, 556, 558, 562,
     & 564, 566, 568, 574, 580, 582, 584, 586, 590, 592, 596, 598,
     & 602, 604, 606, 610, 614, 618, 620, 622, 626, 628, 632, 634,
     & 636, 638, 642, 644, 652, 654, 656, 658, 662, 664, 666, 668,
     & 670, 674, 678, 682, 688, 690, 692, 694, 696, 698, 706, 708,
     & 710, 712, 716, 718, 724, 730, 732, 734, 736, 738, 740, 742,
     & 744, 746, 752, 754, 758, 762, 764, 766, 772, 774, 776, 778,
     & 782, 786, 788, 790, 794, 796, 804, 806, 808, 812, 814, 820,
     & 822, 824, 826, 828, 830, 834, 844, 846, 848, 852, 854, 856,
     & 860, 868, 870, 872, 874, 876, 888, 890, 892, 894, 902, 904,
     & 906, 908, 916, 920, 928, 930, 932, 938, 940, 942, 944, 946,
     & 948, 954, 956, 962, 964, 966, 970, 976, 978, 984, 986, 992,
     & 994, 996/
      mcd=2
      n=(n/mcd)*mcd
      ndiv=0
10       if(n.gt.maxs) then
      n=n/mcd
      ndiv=ndiv+1
      go to 10
         endif
      n=(n/mcd)*mcd
20    continue
       do i=1,mp1
         if(n.eq.p1(i)) then
      n=n-mcd
      go to 20
         endif
       enddo
      if(ndiv.ne.0) n=n*mcd**ndiv
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine shrink(i1,o1,xlow,xup,ylow,yup,zlow,zup,xl,xu,yl,yu,zl,
     & zu,rr,ss)
      implicit none
      integer i1,ix,iy,iz,mxy,nxy,o1,xl,xlow,xu,xup,yl,ylow,yu,yup,zl,
     & zlow,zu,zup
      real rr,ss
      logical last
      dimension rr(xlow:xup,ylow:yup),ss(xl:xu,yl:yu)
      mxy=(xup-xlow+1)*(yup-ylow+1)
      nxy=(xu-xl+1)*(yu-yl+1)
      last=.false.
       do iz=zlow,zup
      if(iz.eq.zup) last=.true.
      call lisec(i1,mxy,rr,last)
         if(iz.ge.zl.and.iz.le.zu) then
      if(iz.eq.zu) last=.true.
        do iy=yl,yu
         do ix=xl,xu
      ss(ix,iy)=rr(ix,iy)
         enddo
        enddo
      call essec(o1,nxy,ss,last)
         endif
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine essec(o1,mxy,rr,last)
      implicit none
      integer i,il,iu,ixy,j,len,length,lm,mxy,o1
      real line,rr
      character card*80
      logical first,last
      external length
      parameter(len=7)
      dimension line(len),rr(mxy)
      data first/.true./
      save il,iu,line,lm
         if(first) then
      first=.false.
      lm=len
      ixy=1
         else
            if(lm.ne.len) then
        do i=1,len-lm
      line(i+lm)=rr(i)
        enddo
      write(card,'(7f10.1)') (line(i),i=1,len)
      call compact(card)
      write(o1,'(a)') card(1:length(card))
            endif
      ixy=1+len-lm
         endif
      il=1-lm
       do j=ixy,mxy,len
      il=il+len
      iu=min(il+len-1,mxy)
      lm=iu-il+1
        do i=1,lm
      line(i)=rr(i+il-1)
        enddo
         if(lm.eq.len.or.last) then
      write(card,'(7f10.1)') (line(i),i=1,lm)
      call compact(card)
      write(o1,'(a)') card(1:length(card))
         endif
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lisec(i1,mxy,rr,last)
      implicit none
      integer i,i1,il,iu,ixy,j,len,lm,ml,mxy
      real line,rr
      logical first,last
      parameter(len=7)
      dimension line(len),rr(mxy)
      data first/.true./
      save il,iu,line,lm,ml
         if(first) then
      first=.false.
      ml=len
      lm=len
      ixy=1
         else
            if(lm.ne.len) then
        do i=1,len-lm
      rr(i)=line(i+lm)
        enddo
            endif
      ixy=1+len-lm
         endif
      il=1-lm
       do j=ixy,mxy,len
      il=il+len
      iu=min(il+len-1,mxy)
      lm=iu-il+1
      if(last) ml=lm
      read(i1,*) (line(i),i=1,ml)
        do i=il,iu
      rr(i)=line(i-il+1)
        enddo
       enddo
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
      subroutine leesc(i1,o1,mx,my,mz,rr,scale)
      implicit none
      integer i,i1,il,iu,ixyz,ms,mx,mxy,mxyz,my,mz,ns,o1,xy
      real map,rr,scale,sigma,ss,sum1,sum2
      dimension rr(mx*my),ss(7)
      sum1=0.
      sum2=0.
      mxy=mx*my
      mxyz=mx*my*mz
      xy=0
      il=-6
       do ixyz=1,mxyz,7
      il=il+7
      iu=min(il+6,mxyz)
      ns=iu-il+1
      read(i1,*) (ss(i),i=1,ns)
      ms=mxy-xy
         if(ms.le.ns) then
        do i=1,ms
      map=ss(i)/scale
      sum1=sum1+map
      sum2=sum2+map**2
      rr(xy+i)=map
        enddo
      write(o1) rr
            if(ms.lt.ns) then
        do i=1,ns-ms
      map=ss(ms+i)/scale
      sum1=sum1+map
      sum2=sum2+map**2
      rr(i)=map
        enddo
            endif
      xy=ns-ms
         else
        do i=1,ns
      map=ss(i)/scale
      sum1=sum1+map
      sum2=sum2+map**2
      rr(xy+i)=map
        enddo
      xy=xy+ns
         endif
       enddo
      sum1=sum1/mxyz
      sum2=sum2/mxyz
      sigma=sqrt(sum2-sum1**2)
      write(6,'(a,2f15.5/a,f12.5)') ' average , sigma :',sum1,sigma,
     . ' set ezd-map SCALE to = ',sigma*scale
      return
      end
