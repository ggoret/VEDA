c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program basic
      implicit none
      integer i,i1,j,length,nat
      double precision bx,cm,diag,dist,dm,tin,xyz
      character atnam*4,card*80,forml*80,lun*80
      dimension bx(2,3),cm(3),tin(3,3),xyz(3)
      external length
      i1=1
       do i=1,80
      forml(i:i)=' '
       enddo
      forml(1:26)='(12X,A4,14X,3F8.3,6X,F6.2)'
10    write(6,'(/a)') ' enter the input-coord filename:'
      read(5,'(a)') lun
      if(lun(1:1).eq.' ') go to 100
      open(unit=i1,file=lun,form='formatted',status='old')
      rewind(unit=i1)
      nat=0
       do i=1,3
      cm(i)=0.
        do j=1,3
      tin(i,j)=0.
        enddo
       enddo
20    read(i1,'(a)',end=30) card
      if(card(1:4).ne.'ATOM') go to 20
      read(card,fmt=forml) atnam,xyz
      nat=nat+1
       do i=1,3
      cm(i)=cm(i)+xyz(i)
        do j=1,3
      tin(i,j)=tin(i,j)+xyz(i)*xyz(j)
        enddo
       enddo
      go to 20
30    continue
       do i=1,3
      cm(i)=cm(i)/nat
       enddo
      diag=0.
       do i=1,3
        do j=1,3
      tin(i,j)=tin(i,j)/nat-cm(i)*cm(j)
        enddo
      diag=diag+tin(i,i)
       enddo
       do i=1,3
      tin(i,i)=diag-tin(i,i)
       enddo
      write(6,'(a,a)') ' coordinates filename: ',lun(1:length(lun))
      write(6,'(a,3F15.5)') ' center-of-mass:',cm
      write(6,'(a,3(/3F10.2))') ' (true) inertia tensor:',
     . ((tin(i,j),j=1,3),i=1,3)
      rewind(unit=i1)
       do i=1,3
      bx(1,i)=0.
      bx(2,i)=0.
       enddo
      dm=0.
40    read(i1,'(a)',end=50) card
      if(card(1:4).ne.'ATOM') go to 40
      read(card,fmt=forml) atnam,xyz
      dist=0.
       do i=1,3
      xyz(i)=xyz(i)-cm(i)
      if(xyz(i).lt.bx(1,i)) bx(1,i)=xyz(i)
      if(xyz(i).gt.bx(2,i)) bx(2,i)=xyz(i)
      dist=dist+xyz(i)**2
       enddo
      if(dist.gt.dm) dm=dist
      go to 40
50    continue
      write(6,'(a,3(/2f10.2))') ' minimal bx:',(bx(1,i),bx(2,i),i=1,3)
      write(6,'(a,f10.2)') ' maximal distance from c-o-m:',sqrt(dm)
      close(unit=i1)
      go to 10
100   stop
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
