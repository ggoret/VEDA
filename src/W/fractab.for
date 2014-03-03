c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program fractab
      implicit none
      integer mc
      parameter(mc=10000000)
      integer ilec,iecr,i1,nort,kmax,lmax,mt,i,hmax,o1
      real a,b,c,alpha,beta,gamma,sqhmax,fract
      complex cc(mc)
      character lun1*80,lun2*80
      ilec=5
      iecr=6
      i1=98
      o1=99
      write(iecr,'(a)') ' enter the input tab filename :'
      read(ilec,'(a)') lun1
      open(unit=i1,file=lun1,form='unformatted',status='old')
      read(i1) a,b,c,alpha,beta,gamma,nort,hmax,kmax,lmax,sqhmax
      mt=(hmax+2)*(2*kmax+1)*(2*lmax+1)
      if(mt.gt.mc) go to 901
      read(i1) (cc(i),i=1,mt)
      write(iecr,'(a)') ' enter the fraction [tabout=tabin/fraction]:'
      read(ilec,*) fract
       do i=1,mt
      cc(i)=cc(i)/fract
       enddo
      write(iecr,'(a)') ' enter the output tab filename :'
      read(ilec,'(a)') lun2
         if(lun2.ne.lun1) then
      open(unit=o1,file=lun2,form='unformatted',status='unknown')
         else
      o1=i1
         endif
      rewind(o1)
      write(o1) a,b,c,alpha,beta,gamma,nort,hmax,kmax,lmax,sqhmax
      write(o1) (cc(i),i=1,mt)
      stop
901   write(iecr,3010) mt
      stop
3010  format('stop >> fractab << set mc=',i10)
      end
