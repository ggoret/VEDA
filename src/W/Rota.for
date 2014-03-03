c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program rota
      implicit none
      integer i,length
      double precision alf,ang,bet,det,dtor,em,error,fuzz,gam,ome,pi,r1,
     & r2,r3,rtod,s,s2,seq,the,twopi,u1,u2,u3
      character card*80
      logical vec1,vec2
      dimension r1(9),r2(9),r3(9),u1(3),u2(3),u3(3)
      common/angkte/ pi,twopi,dtor,rtod
      pi=atan2(1.d0,1.d0)*4.d0
      twopi=atan2(1.d0,1.d0)*8.d0
      dtor=atan2(1.d0,1.d0)/45.d0
      rtod=45.d0/atan2(1.d0,1.d0)
      fuzz=0.1
10    write(6,'(/a)') ' fala > '
      read(5,'(a)') card
      call compact(card)
      if(length(card).eq.0) stop
         if(card(1:1).eq.'h'.or.card(1:1).eq.'H') then
      write(6,2000)
         else if(card(1:1).eq.'e'.or.card(1:1).eq.'E') then
c escribe la matriz de rotacion dada en parametrizacion de euler
      read(card(3:80),*) alf,bet,gam
            if(card(2:2).eq.'-') then
      ang=alf
      alf=180.-gam
      gam=180.-ang
            endif
      call rmxe(alf,bet,gam,r1)
      call rmx2e(r1,alf,bet,gam)
      write(6,'(/a,3f8.3)') ' euler:',alf,bet,gam
      call rmx2p(r1,ang,u1)
      call v2a(u1,the,ome)
      write(6,'(a,f8.3,3f10.5,5x,a,3f8.3)') ' polar:',ang,u1,
     . ';',ang,the,ome
      call prtmx(r1)
         else if(card(1:1).eq.'p'.or.card(1:1).eq.'P') then
c escribe la matriz de rotacion dada en parametrizacion polar
      read(card(3:80),*) ang,u1
      call rmxp(ang,u1,r1)
      call rmx2e(r1,alf,bet,gam)
      write(6,'(/a,3f8.3)') ' euler:',alf,bet,gam
      call rmx2p(r1,ang,u1)
      call v2a(u1,the,ome)
      write(6,'(a,f8.3,3f10.5,5x,a,3f8.3)') ' polar:',ang,u1,
     . ';',ang,the,ome
      call prtmx(r1)
         else if(card(1:1).eq.'r'.or.card(1:1).eq.'R') then
c escribe las formas parametricas de una matriz de rotacion
      write(6,'(a)') ' lee la matriz por lineas [igual que el eco]:'
      call redmx(r2)
            if(card(3:3).eq.'t'.or.card(3:3).eq.'T') then
      write(6,'(/a)') ' eco (transpuesta):'
      call trans(r2)
            else
      write(6,'(/a)') ' eco:'
            endif
      call prtmx(r2)
      call lsqmx(r2,r1)
      call det3(r1,det)
      write(6,'(/a)') ' matriz de cuadrados minimos:'
      call prtmx(r1)
      seq=0.
      em=0.
       do i=1,9
      error=abs(r2(i)-r1(i))
      seq=seq+error**2
      if(error.gt.em) em=error
       enddo
      write(6,'(/a,2f10.5)') ' r-m-s y error max.',sqrt(seq/9.),em
            if(det.le.0.) then
      write(6,'(/a)') ' determinante negativo'
            else
      call rmx2e(r1,alf,bet,gam)
      write(6,'(/a,3f8.3)') ' euler:',alf,bet,gam
      call rmx2p(r1,ang,u1)
      call v2a(u1,the,ome)
      write(6,'(a,f8.3,3f10.5,5x,a,3f8.3)') ' polar:',ang,u1,
     . ';',ang,the,ome
            endif
         else if(card(1:1).eq.'v'.or.card(1:1).eq.'V') then
c escribe un versor dadas sus proyecciones o sus coordenadas polares
            if(card(2:2).eq.'c'.or.card(2:2).eq.'C'.or.
     . card(2:2).eq.' ') then
      read(card(3:80),*) u1
      call v2a(u1,the,ome)
      call a2v(the,ome,u1)
            else if(card(2:2).eq.'a'.or.card(2:2).eq.'A') then
      read(card(3:80),*) the,ome
      call a2v(the,ome,u1)
            else
      go to 901
            endif
      write(6,'(/a,3f10.5,5x,a,2f8.3)') ' versor:',
     . u1,';',the,ome
         else if(card(1:1).eq.'x'.or.card(1:1).eq.'X') then
c multiplica (matriz o vector) * (matriz o vector)
      vec1=.false.
      vec2=.false.
      write(6,'(2a)') ' primer elemento (vector o matriz ;',
     . ' notacion polaca inversa):'
      read(5,'(a)') card
      call compact(card)
            if(card(1:1).eq.'v'.or.card(1:1).eq.'V') then
      vec1=.true.
               if(card(2:2).eq.'c'.or.card(2:2).eq.'C'.or.
     . card(2:2).eq.' ') then
      read(card(3:80),*) u1
               else if(card(2:2).eq.'a'.or.card(2:2).eq.'A') then
      read(card(3:80),*) the,ome
      call a2v(the,ome,u1)
               else
      go to 901
               endif
            else if(card(1:1).eq.'e'.or.card(1:1).eq.'E') then
      read(card(3:80),*) alf,bet,gam
               if(card(2:2).eq.'-') then
      ang=alf
      alf=180.-gam
      gam=180.-ang
               endif
      call rmxe(alf,bet,gam,r1)
            else if(card(1:1).eq.'p'.or.card(1:1).eq.'P') then
      read(card(3:80),*) ang,u1
      call rmxp(ang,u1,r1)
            else
      go to 902
            endif
      write(6,'(a)') ' segundo elemento (vector o matriz):'
      read(5,'(a)') card
      call compact(card)
            if(card(1:1).eq.'v'.or.card(1:1).eq.'V') then
      if(.not.vec1) go to 903
      vec2=.true.
               if(card(2:2).eq.'c'.or.card(2:2).eq.'C'.or.
     . card(2:2).eq.' ') then
      read(card(3:80),*) u2
               else if(card(2:2).eq.'a'.or.card(2:2).eq.'A') then
      read(card(3:80),*) the,ome
      call a2v(the,ome,u2)
               else
      go to 901
               endif
            else if(card(1:1).eq.'e'.or.card(1:1).eq.'E') then
      read(card(3:80),*) alf,bet,gam
               if(card(2:2).eq.'-') then
      ang=alf
      alf=180.-gam
      gam=180.-ang
               endif
      call rmxe(alf,bet,gam,r2)
            else if(card(1:1).eq.'p'.or.card(1:1).eq.'P') then
      read(card(3:80),*) ang,u2
      call rmxp(ang,u2,r2)
            else
      go to 902
            endif
            if(vec1.and.vec2) then
      write(6,'(a)') ' producto escalar [.] o vectorial [x]:'
      read(5,'(a)') card
      call compact(card)
               if(card(1:1).eq.'.') then
      call dotvv(s,u2,u1)
      write(6,'(/a,f10.5)') ' producto escalar:',s
               else if(card(1:1).eq.'x') then
      call vecvv(u3,u2,u1)
      write(6,'(/a,3f10.5)') ' vector producto:',u3
               else
      go to 904
               endif
            else if(vec1) then
      call promv(u3,r2,u1)
      call v2a(u3,the,ome)
      write(6,'(/a,3f10.5,5x,a,2f8.3)') ' vector girado:',
     . u3,';',the,ome
            else
      call pro2mx(r3,r2,r1)
      write(6,'(/a)') ' rotacion producto'
      call rmx2e(r3,alf,bet,gam)
      write(6,'(a,3f8.3)') ' euler:',alf,bet,gam
      call rmx2p(r3,ang,u3)
      call v2a(u3,the,ome)
      write(6,'(a,f8.3,3f10.5,5x,a,3f8.3)') ' polar:',ang,u3,
     . ';',ang,the,ome
      call prtmx(r3)
            endif
         else if(card(1:1).eq.'i'.or.card(1:1).eq.'I') then
c interpreta R*x+T como R*(x-o)+o+u con o.u=0 y R*u=u
      write(6,'(a)') ' lee la matriz de rotacion:'
      read(5,'(a)') card
      call compact(card)
            if(card(1:1).eq.'e'.or.card(1:1).eq.'E') then
      read(card(3:80),*) alf,bet,gam
               if(card(2:2).eq.'-') then
      ang=alf
      alf=180.-gam
      gam=180.-ang
               endif
      call rmxe(alf,bet,gam,r1)
      call rmx2p(r1,ang,u1)
            else if(card(1:1).eq.'p'.or.card(1:1).eq.'P') then
      read(card(3:80),*) ang,u1
      call rmxp(ang,u1,r1)
            else
      go to 905
            endif
      write(6,'(a)') ' lee la translacion:'
      read(5,'(a)') card
      call compact(card)
            if(card(1:2).eq.'vc'.or.card(1:2).eq.'VC'.or.
     . card(1:2).eq.'v ') then
      read(card(3:80),*) u2
      call v2a(u2,the,ome)
      call a2v(the,ome,u3)
      call dotvv(s,u1,u3)
      s=1.d0-s*s
               if(abs(s).le.1.d-6) then
       do i=1,3
      u3(i)=0.d0
      u1(i)=u2(i)
       enddo
      go to 20
               endif
            else
      go to 906
            endif
      write(6,'(/a/a,f8.3/a,3f10.5)') ' rotacion',
     . ' angulo de rotacion:',ang,' eje de rotacion:',u1
            if(ang.lt.fuzz) then
      call v2a(u2,the,ome)
      call a2v(the,ome,u1)
      ang=180.
            endif
      call vecvv(u3,u2,u1)
      s=0.
       do i=1,3
      s=s+u2(i)*u1(i)
       enddo
       do i=1,3
      u1(i)=s*u1(i)
       enddo
      s2=sin(ang*dtor)/(1-cos(ang*dtor))
       do i=1,3
      u3(i)=0.5*(u2(i)-u1(i)+s2*u3(i))
       enddo
      write(6,'(/a/a,3f10.5)') ' verificacion',' translacion (input):',
     . u2
      call promv(u2,r1,u3)
       do i=1,3
      u3(i)=-u2(i)+u3(i)+u1(i)
       enddo
      write(6,'(a,3f10.5)') ' translacion (output):',u3
20    write(6,'(/a,3f10.5)') ' centro de rotacion:',u3
      write(6,'(a,3f10.5)') ' translacion       :',u1
         else if(card(1:1).eq.'s'.or.card(1:1).eq.'S') then
c exit *************************************************************************
      stop
         else
      write(6,2000)
         endif
      go to 10
901   write(6,'(a)') ' debe ser: v , vc, va'
      stop
902   write(6,'(a)') ' debe ser: v , vc, va, e, p'
      stop
903   write(6,'(a)') ' ahora debe ser: e, p'
      stop
904   write(6,'(a)') ' ahora debe ser: ., x'
      stop
905   write(6,'(a)') ' debe ser: e, p'
      stop
906   write(6,'(a)') ' debe ser: v , vc'
      stop
2000  format(
     . '--------------------------------------------------------------'/
     . '    POSIBLES OPCIONES'//
     . 'e  {alf} {beta} {gama}'/
     . 'e- inversa de {alf} {beta} {gama}'/
     . 'p  {ang} {cosx} {cosy} {cosz}'/
     . 'r  {matriz}'/
     . 'v  {vector}'/
     . 'vc {cosx} {cosy} {cosz}'/
     . 'va {colatitud} {longitud}'/
     . 'x  (calcula productos binarios varios; notacion polonesa).'/
     . 'i  (interpreta rotacion-translacion como rotacion alrededor de'/
     . '    un centro seguida de translacion paralela al eje).'/
     . 's  (exit).'/
     . 'h  (kelp).'/
     . 'blank  (exit).'/
     . '--------------------------------------------------------------')
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rmx2e(rotx,alpha,beta,gamma)
c this subroutine determines the eulerian representation of a rotation matrix
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
      subroutine rmx2p(rotx,xhi,vn)
c this subroutine determines the polar representation of a rotation matrix
      implicit none
      integer i
      double precision ang,cang,dtor,fuzz,norm,pi,rotx,rtod,sang,scale,
     & twopi,vecx,vn,xhi
      dimension rotx(3,3),vecx(3),vn(3)
      common/angkte/ pi,twopi,dtor,rtod
      fuzz=0.01
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
      return
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
      subroutine rmxp(xhi,vn,rotx)
c this subroutine writes the rotation matrix, given its polar parameters
c the vector defining the axis may be not normalized
      implicit none
      integer eijk,i,j,k
      double precision ang,cang,dtor,norm,pi,rotx,rtod,sang,twopi,vn,xhi
      dimension eijk(3,3,3),rotx(3,3),vn(3)
      common/angkte/ pi,twopi,dtor,rtod
      data eijk/
     & 0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
      norm=sqrt(vn(1)**2+vn(2)**2+vn(3)**2)
       do i=1,3
      vn(i)=vn(i)/norm
       enddo
      ang=dtor*xhi
      cang=cos(ang)
      sang=sin(ang)
      do 30 j=1,3
      do 20 i=1,3
      rotx(i,j)=(1-cang)*vn(i)*vn(j)
      if(i.eq.j) rotx(i,j)=rotx(i,j)+cang
      do 10 k=1,3
         if(eijk(i,k,j).gt.0) then
      rotx(i,j)=rotx(i,j)+sang*vn(k)
         else if(eijk(i,k,j).lt.0) then
      rotx(i,j)=rotx(i,j)-sang*vn(k)
         endif
10    continue
20    continue
30    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pro2mx(ro,r2,r1)
c this subroutine calculates the product of two matrices: ro=r2*r1
      implicit none
      integer i,j,k
      double precision add,r1,r2,ro
      dimension r1(3,3),r2(3,3),ro(3,3)
      do 30 i=1,3
      do 20 j=1,3
      add=0.d0
      do 10 k=1,3
      add=add+r2(i,k)*r1(k,j)
10    continue
      ro(i,j)=add
20    continue
30    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine promv(uo,r2,u1)
c this subroutine calculates the product of a matrix times a vector: uo=r2*u1
      implicit none
      integer i,j
      double precision add,r2,u1,uo
      dimension r2(3,3),u1(3),uo(3)
      do 20 i=1,3
      add=0.d0
      do 10 j=1,3
      add=add+r2(i,j)*u1(j)
10    continue
      uo(i)=add
20    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine dotvv(s,u2,u1)
c this subroutine calculates the dot product of two vectors
      implicit none
      integer i
      double precision s,u1,u2
      dimension u1(3),u2(3)
      s=0.
      do 10 i=1,3
      s=s+u2(i)*u1(i)
10    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine vecvv(u3,u2,u1)
c this subroutine calculates the skew product of two vectors: u1xu2=u3
      implicit none
      integer eijk,i,j,k
      double precision u1,u2,u3
      dimension eijk(3,3,3),u1(3),u2(3),u3(3)
      data eijk/
     & 0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
      do 30 i=1,3
      u3(i)=0.
      do 20 j=1,3
      do 10 k=1,3
         if(eijk(i,j,k).gt.0) then
      u3(i)=u3(i)+u1(j)*u2(k)
         else if(eijk(i,j,k).lt.0) then
      u3(i)=u3(i)-u1(j)*u2(k)
         endif
10    continue
20    continue
30    continue
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine a2v(thed,omed,u1)
c this subroutine determines the unit vector, given the colatitude and longitud
      implicit none
      double precision come,cthe,dtor,ome,omed,pi,rtod,some,sthe,the,
     & thed,twopi,u1
      dimension u1(3)
      common/angkte/ pi,twopi,dtor,rtod
      ome=omed/rtod
      the=thed/rtod
      sthe=sin(the)
      cthe=cos(the)
      some=sin(ome)
      come=cos(ome)
      u1(1)=sthe*come
      u1(2)=sthe*some
      u1(3)=cthe
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine v2a(u1,thed,omed)
c this subroutine determines the colatitude and longitud, given the vector
      implicit none
      double precision cthe,dtor,fuzz,ome,omed,pi,rtod,sthe,thed,twopi,
     & u1,unorm
      dimension u1(3)
      common/angkte/ pi,twopi,dtor,rtod
      fuzz=0.001
      unorm=sqrt(u1(1)**2+u1(2)**2+u1(3)**2)
      cthe=min(max(u1(3)/unorm,-1.),1.)
      sthe=sqrt(1-cthe**2)
         if(sthe.gt.fuzz) then
      ome=atan2(u1(2),u1(1))
         else
      ome=0.
         endif
      thed=acos(cthe)*rtod
      omed=mod(twopi+ome,twopi)*rtod
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine redmx(r)
c this subroutine reads a matrix
      implicit none
      integer i,j
      double precision r
      dimension r(3,3)
      read(5,*) ((r(i,j),j=1,3),i=1,3)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trans(r)
c this subroutine transposes a matrix
      implicit none
      integer i,j
      double precision r,r1
      dimension r(3,3),r1(3,3)
       do i=1,3
        do j=1,3
      r1(i,j)=r(j,i)
        enddo
       enddo
       do i=1,3
        do j=1,3
      r(i,j)=r1(i,j)
        enddo
       enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine det3(r,d)
c this subroutine computes the determinant of a 3x3 matrix
      implicit none
      double precision a,b,d,r
      dimension r(3,3)
      a=r(1,1)*r(2,2)*r(3,3)+r(1,2)*r(2,3)*r(3,1)+r(3,2)*r(2,1)*r(1,3)
      b=r(1,3)*r(2,2)*r(3,1)+r(2,1)*r(1,2)*r(3,3)+r(2,3)*r(3,2)*r(1,1)
      d=a-b
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine prtmx(r)
c this subroutine prints a matrix
      implicit none
      integer i,j
      double precision r
      dimension r(3,3)
      write(6,'(3(/3f10.5))') ((r(i,j),j=1,3),i=1,3)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lsqmx(r1,r2)
c this subroutine determines the rotation matrix closest to a given matrix
      implicit none
      integer i,j,k
      double precision mat,r1,r2,r3,rr,sum,val
      dimension mat(3,3),r1(3,3),r2(3,3),r3(3,3),rr(3),val(3)
       do i=1,3
        do j=1,3
      sum=0.d0
         do k=1,3
      sum=sum+r1(i,k)*r1(j,k)
         enddo
      r3(i,j)=sum
        enddo
       enddo
      call rs(3,3,r3,val,mat,rr)
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
      subroutine rs(nm,n,a,w,z,fv1)
      implicit none
      integer n,nm
      double precision a,fv1,w,z
      dimension a(nm,n),fv1(n),w(n),z(nm,n)
      call tred2(nm,n,a,w,fv1,z)
      call tql2(nm,n,w,fv1,z)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tql2(nm,n,d,e,z)
      implicit none
      integer i,ii,j,k,l,l1,l2,m,mml,n,nm
      double precision c,c2,c3,d,dl1,e,el1,f,g,h,p,pythag,r,s,s2,tst1,
     & tst2,z
      external pythag
      dimension d(n),e(n),z(nm,n)
cjn uninitialized values
      data c3,s2/0.,0./
      if(n.eq.1) go to 400
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
130   if(j.eq.30) stop ' > tql2 < '
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
400   return
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
