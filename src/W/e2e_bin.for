c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program e2e_bin
      implicit none
      integer i1,iecr,ilec,kprt,length,mr,mx,my,mz,nx,ny,nz,o1,step,xl,
     & xlow,xu,xup,yl,ylow,yu,yup,zl,zlow,zu,zup
      real a,alpha,b,beta,c,gamma,rr,scale,ss
      character card*80,lun*80
      logical inv
      external length
      parameter(mr=100 000 000)
      dimension rr(mr),ss(mr)
      common/ioprg/ ilec,iecr,kprt
      ilec=5
      iecr=6
      kprt=1
      i1=98
      o1=99
      write(iecr,'(/a)') ' ENTER the input ezd-map filename'
      read(ilec,'(a)') lun
      if(length(lun).eq.0) go to 20
      open(unit=i1,file=lun,form='formatted',status='old')
      read(i1,*)
10    read(i1,'(a)') card
      if(card(1:1).eq.'!') go to 10
      read(card(5:80),*) a,b,c,alpha,beta,gamma
      read(i1,'(a)') card
      read(card(7:80),*) xl,yl,zl
      read(i1,'(a)') card
      read(card(7:80),*) mx,my,mz
      read(i1,'(a)') card
      read(card(5:80),*) nx,ny,nz
      read(i1,'(a)') card
      read(card(6:80),*) scale
      read(i1,*)
      write(iecr,'(/a/6f10.2/3i5)') ' input parameters:',
     . a,b,c,alpha,beta,gamma,nx,ny,nz
      xu=xl+mx-1
      yu=yl+my-1
      zu=zl+mz-1
      write(iecr,'(/a,6i5/a,3i5)') ' input indices [x],[y],[z]:',
     . xl,xu,yl,yu,zl,zu,' extent:',mx,my,mz
         if(mx*my*mz.gt.mr) then
      write(iecr,'(/a,i10)') 'stop >> e2e_bin << set mr =',mx*my*mz
      stop
         endif
      call lee(i1,mx,my,mz,rr(1))
      write(iecr,'(a)') ' ENTER the step'
      read(ilec,*) step
      xlow=xl/step
      if(xlow*step.lt.xl) xlow=xlow+1
      ylow=yl/step
      if(ylow*step.lt.yl) ylow=ylow+1
      zlow=zl/step
      if(zlow*step.lt.zl) zlow=zlow+1
      xup=xu/step
      if(xup*step.gt.xu) xup=xup-1
      yup=yu/step
      if(yup*step.gt.yu) yup=yup-1
      zup=zu/step
      if(zup*step.gt.zu) zup=zup-1
      mx=xup-xlow+1
      my=yup-ylow+1
      mz=zup-zlow+1
      inv=.true.
      call bin(xl,xu,yl,yu,zl,zu,xlow,xup,ylow,yup,zlow,zup,rr(1),ss(1),
     & step)
c ----------------------------------------------------------------------
      inv=.false.
         if(inv) then
      xu=-xlow
      yu=-ylow
      zu=-zlow
      xlow=-xup
      ylow=-yup
      zlow=-zup
      xup=xu
      yup=yu
      zup=zu
         endif
c ----------------------------------------------------------------------
      a=mx*(a/nx)*step
      b=my*(b/ny)*step
      c=mz*(c/nz)*step
      nx=mx
      ny=my
      nz=mz
      xl=xlow
      yl=ylow
      zl=zlow
      write(iecr,'(/a)') ' ENTER the output ezd-map filename'
      read(ilec,'(a)') lun
      open(unit=o1,file=lun,form='formatted',status='unknown')
      write(iecr,'(/a/6f10.2/3i5)') ' output parameters:',
     . a,b,c,alpha,beta,gamma,nx,ny,nz
      write(iecr,'(/a,6i5/a,3i5)') ' output indices [x],[y],[z]:',
     . xlow,xup,ylow,yup,zlow,zup,' extent:',mx,my,mz
      write(o1,'(a)') 'EZD_MAP'
      write(card,'(a)') '! created by e2e_bin'
      call compact(card)
      write(o1,'(a)') card(1:length(card))
      write(card,'(a4,6f10.3)') 'CELL',a,b,c,alpha,beta,gamma
      call compact(card)
      write(o1,'(a)') card(1:length(card))
      write(card,'(a6,3i5)') 'ORIGIN',xl,yl,zl
      call compact(card)
      write(o1,'(a)') card(1:length(card))
      write(card,'(a6,3i5)') 'EXTENT',mx,my,mz
      call compact(card)
      write(o1,'(a)') card(1:length(card))
      write(card,'(a4,3i5)') 'GRID',nx,ny,nz
      call compact(card)
      write(o1,'(a)') card(1:length(card))
      write(card,'(a5,f15.3)') 'SCALE',scale
      call compact(card)
      write(o1,'(a)') card(1:length(card))
      write(o1,'(a)') 'MAP'
      call escribe(o1,mx,my,mz,ss(1))
      write(o1,'(a)') 'END'
20    stop
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lee(i1,mx,my,mz,rr)
      implicit none
      integer i,i1,iecr,il,ilec,iu,ixyz,kprt,mx,mxyz,my,mz
      real map,mmax,mmin,rr,sigma,sum,sum2
      dimension rr(mx*my*mz)
      common/ioprg/ ilec,iecr,kprt
      common/stats/ sigma
      write(iecr,'(/a)') ' >>>>>>>>>> READING MAP <<<<<<<<<<'
      mxyz=mx*my*mz
      il=-6
       do ixyz=1,mxyz,7
      il=il+7
      iu=min(il+6,mxyz)
      read(i1,*) (rr(i),i=il,iu)
       enddo
      mmin=1.e10
      mmax=-1.e10
      sum=0.
      sum2=0.
       do ixyz=1,mxyz
      map=rr(ixyz)
      if(map.lt.mmin) mmin=map
      if(map.gt.mmax) mmax=map
      sum=sum+map
      sum2=sum2+map**2
       enddo
      sum=sum/mxyz
      sum2=sum2/mxyz
      sigma=sqrt(sum2-sum**2)
      write(iecr,'(/a,3e15.5)') ' >> min, max, sigma <<',mmin,mmax,sigma
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine bin(xl,xu,yl,yu,zl,zu,xlow,xup,ylow,yup,zlow,zup,rr,ss,
     & step)
      implicit none
      integer iecr,ilec,ix,iy,iz,jx,jy,jz,kprt,step,sx,sy,sz,xl,xlow,xu,
     & xup,yl,ylow,yu,yup,zl,zlow,zu,zup
      real rr,ss
      dimension rr(xl:xu,yl:yu,zl:zu),ss(xlow:xup,ylow:yup,zlow:zup)
      common/ioprg/ ilec,iecr,kprt
      sx=xlow*step
      sy=ylow*step
      sz=zlow*step
      write(iecr,'(/a)') ' >>>>>>>>>> SKIPPING <<<<<<<<<< '
cjn  . sx-xl,sy-yl,sz-zl
      jz=zlow-1
       do iz=sz,zu,step
      jz=jz+1
      jy=ylow-1
        do iy=sy,yu,step
      jy=jy+1
      jx=xlow-1
         do ix=sx,xu,step
      jx=jx+1
      ss(jx,jy,jz)=rr(ix,iy,iz)
         enddo
        enddo
       enddo
         if(jx.ne.xup.or.jy.ne.yup.or.jz.ne.zup) then
      write(iecr,'(/a/3i5,5x,3i5)') ' j[xyz] .ne. [xyz]up :',
     . jx,jy,jz,xup,yup,zup
         endif
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine escribe(o1,mx,my,mz,rr)
      implicit none
      integer i,il,iu,ixyz,length,mx,mxyz,my,mz,o1
      real rr
      character card*80
      external length
      dimension rr(mx*my*mz)
      mxyz=mx*my*mz
      il=-6
       do ixyz=1,mxyz,7
      il=il+7
      iu=min(il+6,mxyz)
      write(card,'(7f10.1)') (rr(i),i=il,iu)
      call compact(card)
      write(o1,'(a)') card(1:length(card))
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
