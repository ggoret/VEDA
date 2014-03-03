# to recover individual items, csh this file
#=======================================================================
echo diagstd.f
cat >diagstd.f <<"ENDOF diagstd.f"
CoC====================================================================
CoC   Diagstd:
CoC
CoC   Diagonalization of a matrix, real, symmetrical.
CoC
CoC   Diagonalization routine: TQLI (EISPACK).
CoC  -simple, public domain, but slow (ALL eigenvectors are computed).
CoC
CoC====================================================================
CoC
CoC   INPUT matrix filename (expected from, e.g., PDBMAT program):
CoC   ************************************************************
CoC   CERFACS -Formatted : matrix.sdijf
CoC   CERFBIN -Binary    : matrix.sdijb
CoC
CoC   Input matrix format: i, j, non-zero-ij-element.
CoC
CoC   It can start with a title, recognized by:
CoC   !,# in first column, or 'program-name>' as first word.
CoC
CoC   OUTPUT:
CoC   *******
CoC   Eigenvector filename: matrix.eigenfacs (in CERFACS format)
CoC
CoC.....................................................................
      program diagstd
      implicit none
      integer natmax, ndim, nvecout
CoC
CoC   MEMORY LIMITS:
CoC   **************
CoC   NATMAX: Maximum number of atoms allowed.
CoC
      parameter( natmax=5000,
     .           ndim=3*natmax,
     .           nvecout=26 )
CoC.....................................................................
c     YHS-Sep-2002: Premiere version, from Diagijr v1.13 (Bordeaux).
c.......................................................................
      logical qcrois, qexist, qinterr
      integer evord(ndim), i, ii, j, jj, k, lmot, natom, nbig, nmots,
     .        nmotsmx, nord, nredond, ntit, ntrace, nvec, nunit, rdunit,
     .        unmess, unmodes
      double precision amat(ndim,ndim), ev(ndim), evsort(ndim),
     .       matrd, trace, work(ndim), som, x, y, z, dist, dmax
      parameter(nmotsmx=80)
      character cformat*20, cstatus*20, eige*4,
     .       lignrd*50, lign80*80, matrice*20, mots(nmotsmx)*80,
     .       nomfich*20,
     .       program*9, progrer*12, progrwn*12, version*32
c.......................................................................
      version=' Version 1.10, March 2008.'
      program=' Diagstd>'
      progrer='%Diagstd-Er:'
      progrwn='%Diagstd-Wn:'
c     Sortie standard:
      unmess=6
      write(unmess,'(2A)') program,version
CoC   Eigenvector are given by increasing eigenvalues (LOWE).
c     Or by decreasing values (HIGH).
      eige='LOWE'
      nvec=3*natmax
      nunit=10
      rdunit=nunit
      nunit=nunit+1
c     Detection de la matrice d'entree:
c     --------------------------------
      cformat='FORMATTED'
      matrice='CERFACS'
      cstatus='OLD'
      nomfich='matrix.sdijf'
      inquire(file=nomfich,exist=qexist)
      if (qexist) goto 50
      nomfich='matrice.sdijf'
      inquire(file=nomfich,exist=qexist)
      if (qexist) goto 50
      nomfich='pdbmat.sdijf'
      inquire(file=nomfich,exist=qexist)
      if (qexist) goto 50
      cformat='UNFORMATTED'
      matrice='CERFBIN'
      nomfich='matrix.sdijb'
      inquire(file=nomfich,exist=qexist)
      if (qexist) goto 50
      nomfich='matrice.sdijb'
      inquire(file=nomfich,exist=qexist)
      if (qexist) goto 50
      nomfich='pdbmat.sdijb'
      inquire(file=nomfich,exist=qexist)
      if (qexist) goto 50
      write(unmess,'(/2A)') progrer,' Matrix not found.'
      write(unmess,'(A)')
     .' Expected filenames are: ',
     .' matrix.sdijf  (CERFACS) -formatted, free format.',
     .' matrice.sdijf ',
     .' pdbmat.sdijf  ',
     .' matrix.sdijb  (CERFBIN) -unformatted, free format.',
     .' matrice.sdijb ',
     .' pdbmat.sdijb  '
      stop '*Required*'
  50  continue
      call openam1(nomfich,cformat,cstatus,rdunit,.false.,
     .     qinterr,qexist)
      if (qinterr) stop '*Matrix file could not be opened*'
      write(unmess,'(3A)') program,
     .    ' Matrix to be read from file: ',nomfich
c     ============================================
c     Lecture matrice d'entree (CERFACS, CERFBIN):
c     ============================================
      write(unmess,'(/2A)') program,
     .    ' Matrix to be read is in CERFACS Format.'
c     1) La matrice a un titre ?
c     --------------------------
      ntit=0
  55  continue
      if (matrice.eq.'CERFACS') then
          read(rdunit,'(A)',end=60) lign80
      else
          read(rdunit,end=60) lign80
      endif
      if (lign80(1:1).eq.'!'.or.lign80(1:1).eq.'#') then
          ntit=ntit+1
          write(unmess,'(2A)') lign80(1:50),' ...'
          goto 55
      else
          lignrd=lign80(1:50)
          call string_split(lign80,80,' ',mots,nmotsmx,nmots)
          call stringcl(mots(1),lmot)
          if (mots(1)(lmot:lmot).eq.'>') then
              ntit=ntit+1
              write(unmess,'(2A)') lignrd,' ...'
              goto 55
          endif
      endif
  60  continue
      rewind(rdunit)
      write(unmess,'(2A,I6,A)') program,' It has ',ntit,' title lignes.'
c     2) Ordre de la matrice, nombre de lignes:
c     -----------------------------------------
      if (ntit.gt.0) then
          do i=1,ntit
             if (matrice.eq.'CERFACS') then
                 read(rdunit,*)
             else
                 read(rdunit)
             endif
          enddo
      endif
      k=0
      nord=0
  90  continue
      if (matrice.eq.'CERFACS') then
      read(rdunit,*,end=100) i,j
      else
      read(rdunit,end=100) i,j
      endif
      k=k+1
      if (i.le.0.or.j.le.0) then
          write(unmess,'(/2A,I9,2(A,I6))')
     .    progrer,' in ligne: ',k,' I= ',i,' J= ',j
          stop
      endif
      if (i.gt.nord) nord=i
      if (j.gt.nord) nord=j
      goto 90
 100  continue
      write(unmess,'(/2A,I9)')
     .     program,' Matrix dimension  (Nord)  =',nord
      write(unmess,'(2A,I9)')
     .     program,' Number of non-zero elements',k
      natom=nord/3
      if (natom.gt.natmax.or.nord.gt.ndim) then
          write(unmess,'(2A)')
     .    progrer,' Matrix can not be read.'
          if (natom.gt.natmax) write(unmess,'(2(A,I9))')
     .   ' Natom= ',natom,' > natmax= ',natmax
          if (nord.gt.ndim) write(unmess,'(2(A,I9))')
     .   ' Nord=  ',nord,' > Ndim=  ',ndim
          stop
      endif
c     3) Lecture de la matrice:
c     -------------------------
      rewind(rdunit)
      if (ntit.gt.0) then
          do i=1,ntit
          if (matrice.eq.'CERFACS') then
              read(rdunit,*)
          else
              read(rdunit)
          endif
          enddo
      endif
      nredond=0
      ntrace=0
      trace=0.d0
      nbig=0
      do i=1,nord
        do j=1,nord
         amat(i,j)=0.d0
        enddo
      enddo
      do jj=1,k
         if (matrice.eq.'CERFACS') then
         read(rdunit,*,err=95) i,j,matrd
         else
         read(rdunit,err=95) i,j,matrd
         endif
         if (dabs(matrd).gt.0.d0) then
             amat(i,j)=matrd
             amat(j,i)=matrd
             if (i.eq.j) then
                trace=trace+matrd
                ntrace=ntrace+1
             endif
             if (matrd.gt.1E+10) then
                 nbig=nbig+1
                 if (nbig.lt.10) then
                     write(unmess,'(2A,2I12,A,G12.3)')
     .               progrwn,' Element: ',i,j,' = ',matrd
                 else
                     if (nbig.eq.10) write(unmess,*) '...'
                 endif
             endif
         else
             nredond=nredond+1
         endif
      enddo
      goto 105
  95  continue
      write(unmess,'(2A,I6)')
     .     progrer,' while reading ligne ',k
      write(unmess,'(2I6,F16.8)') ' i, j, matrd= ',i,j,matrd
      stop
 105  continue
      write(unmess,'(2A,I9)') program,
     .    ' Nb of elements found twice:',nredond
      if (nredond.gt.0)
     .write(unmess,'(2A/)') progrwn,' Ok ?'
      write(unmess,'(2A,I9)') program,
     .    ' Nb of elements    > 1E+10 :',nbig
      if (nbig.gt.0)
     .write(unmess,'(2A/)') progrwn,' Ok ?'
      write(unmess,'(2A,F31.7)') program,
     .    ' Matrix trace:',trace
      if (nord-ntrace.gt.0)
     .write(unmess,'(2A,I7)') progrwn,
     .    ' Nb of zero elements there:',nord-ntrace
c     Diagonalisation:
c     ----------------
      nomfich='matrix.eigenfacs'
      cformat='FORMATTED'
      cstatus='ove'
      unmodes=nunit
      nunit=nunit+1
      call openam1(nomfich,cformat,cstatus,
     .     unmodes,.true.,
     .     qinterr,qexist)
      write(unmess,'(/2A)') program,' Diagonalization.'
      if (nvec.gt.nord) nvec=nord
      write(unmess,'(A,I6,A)') program,
     .      nvec,' eigenvectors are about to be computed. '
c      nvecout=nvec
c     Initialisations:
      do i=1,ndim
         ev(i)=0.d0
      enddo
c     Eigenvalues/Matrix Diagonalization
CoC
CoC   TRED2 and TQLI (based on the original EISPACK library)
CoC   perform a diagonalization of a real symmetric matrix based
CoC   on the QL algorithm.
      CALL TRED2(amat,nord,ndim,ev,work)
      CALL TQLI(ev,work,nord,ndim,amat)
      trace=0.d0
      do i=1,nvec
         trace=trace+ev(i)
      enddo
      write(unmess,'(/2A,F24.7)') program,
     .     ' Sum of eigenvalues =',trace
c     Trier par ordre croissant ou decroissant:
      qcrois=.true.
      if (eige.eq.'HIGH') qcrois=.false.
      call trieri(ev,nvec,ndim,evsort,evord,qcrois)
      write(unmess,'(/2A/(5F15.7))') program,
     .    ' Eigenvalues: ',(ev(evord(i)),i=1,nvecout)
      WRITE(unmess,'(/2A/(5F15.7))') program,
     .    ' Frequencies (cm-1, '//
     .     'if the matrix is a hessien in CHARMM units):',
     .    (sqrt(dabs(ev(evord(i))))*108.591365,i=1,nvecout)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Normalisation des modes normaux:
c     -----------------------------------------------
      if (nord.ne.3*natom) stop
      do j=1,nvecout
         i=evord(j)
         som=0.
         do k=1,nord
         som=som+amat(k,i)**2
         enddo
         som=sqrt(som/natom)
         do k=1,nord
         amat(k,i)=amat(k,i)/som
         enddo
      enddo
c     Calcul du deplacement maximal:
c     -----------------------------------------------
      do j=7,nvecout
         i=evord(j)
         dmax=0.
         do k=1,nord
            ii=3*(k-1)
            x=amat(ii+1,i)
            y=amat(ii+2,i)
            z=amat(ii+3,i)
            dist=x**2+y**2+z**2
            if (dist.ge.dmax) dmax = dist
         enddo
      write(6,'(a,i5,f10.5)') ' Maximal displacement', j,sqrt(dmax)
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Ecriture des modes normaux au format 'CERFACS':
c     -----------------------------------------------
      do j=1,nvecout
         i=evord(j)
         write(unmodes,'(A,I5,7X,A,1PG12.4)') ' VECTOR',j,'VALUE',ev(i)
         write(unmodes,'(1X,35A)') ('-',k=1,35)
         write(unmodes,'(3(1PG12.4))') (amat(k,i),k=1,nord)
      enddo
      write(unmess,'(/2A)')
     .      program,' Normal end.'
      stop
      end
"ENDOF diagstd.f"
#=======================================================================
echo getchi.f
cat >getchi.f <<"ENDOF getchi.f"
c-------------------------------------------
      subroutine getchi(message,numlu,qok)
c
c     NUMLU obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-oct-96
c
      implicit none
cI/O:
      double precision numlu
      logical qok
      character*(*) message
cLocal:
      integer ntry, ntrymx
      double precision iread,istep
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getchi> ',message
      read(5,*,end=200,err=100) istep,iread
      numlu=istep*iread
c
      write(6,*) 'Getchi> ',numlu
c
      qok=.true.
      return
 200  continue
      return
      end
"ENDOF getchi.f"
#=======================================================================
echo getnam.f
cat >getnam.f <<"ENDOF getnam.f"
c---------------------------------------------------
      subroutine getnam(message,nomlu,lnomlu,qok)
c
c     NOMLU obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-oct-96
c
      implicit none
cI/O:
      integer lnomlu
      logical qok
      character*(*) message, nomlu
cLocal:
      integer ntry, ntrymx
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getnam> ',message
      read(5,'(A)',end=200,err=100) nomlu
c
      call stringcl(nomlu,lnomlu)
      write(6,'(A,A)') ' Getnam> ',nomlu(1:lnomlu)
c
      qok=.true.
      return
 200  continue
      return
      end
"ENDOF getnam.f"
#=======================================================================
echo getnum.f
cat >getnum.f <<"ENDOF getnum.f"
      subroutine getnum(message,numlu,nummin,nummax,qok)
c
c     NUMLU obtenu en reponse au MESSAGE.
c     NUMLU doit etre inferieur a nummax et superieur a nummin.
c
c     NTRYMX essais en cas de probleme.
c     qok=.false. => Probleme a la lecture.
c
c     YHS-oct-1996: version 1.0
c     YHS-nov-2000: version 3.0
c
      implicit none
cI/O:
      integer numlu, nummax, nummin
      logical qok
      character*(*) message
cLocal:
      integer ntry, ntrymx, iread
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getnum> ',message
      read(5,*,end=200,err=100) iread
      numlu=iread
c
      write(6,*) 'Getnum> ',numlu
      if (nummin.le.nummax) then
      if (numlu.gt.nummax) then
          write(6,'(A,I6,A)')
     .  '%Getnum-Err: number larger than ',nummax,
     .  ' This is not allowed. Sorry.'
          numlu=nummax
          return
      else if (numlu.lt.nummin) then
          write(6,'(A,I6,A)')
     .  '%Getnum-Err: number smaller than ',nummin,
     .  ' This is not allowed. Sorry.'
          numlu=nummin
          return
      endif
      endif
c
      qok=.true.
      return
 200  continue
      return
      end
"ENDOF getnum.f"
#=======================================================================
echo getrep.f
cat >getrep.f <<"ENDOF getrep.f"
c---------------------------------------------------
      subroutine getrep(message,qinfo,qok)
c
c     qinfo obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-jan-00
c     YHS-oct-00
c
      implicit none
cI/O:
      logical qok, qinfo
      character*(*) message
cLocal:
      integer ntry, ntrymx
      character*1 cread
cBegin:
      ntrymx=2
c
      qinfo=.false.
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) then
          goto 200
          return
      endif
c
      write(6,'(A,A)') ' Getrep> ',message
      read(5,'(A)',end=200,err=100) cread
c
      if (cread.eq.'T'.or.cread.eq.'t'.or.
     .    cread.eq.'Y'.or.cread.eq.'y'.or.
     .    cread.eq.'O'.or.cread.eq.'o') then
          qinfo=.true.
      else
      if (cread.ne.'F'.and.cread.ne.'f'.and.
     .    cread.ne.'N'.and.cread.ne.'n')
     .    write(6,'(3A)') '%Getrep-W> Unexpected answer:',cread,
     . '. Assumed answer is: NO.'
      endif
c
      write(6,*) 'Getrep> ',qinfo
c
      qok=.true.
      return
 200  continue
      write(6,'(A)') '%Getrep-W> No answer.'//
     .    ' Assumed answer is: NO.'
      return
      end
"ENDOF getrep.f"
#=======================================================================
echo mintomaj.f
cat >mintomaj.f <<"ENDOF mintomaj.f"
c-----------------------------------------------------------------------
      subroutine mintomaj(chaine)
c     Les caracteres minuscules sont mis en MAJUSCULES.
c     Les autres ne sont pas touches.
c     YHS-Oct-98: Premiere version (Toulouse).
c     YHS-Sep-03: Dernieres modifications (Lyon).
      character*(*) chaine
c Local:
      integer icar, ilettre, taille
      character*26  carmaj, carmin
      carmin='qwertyuiopasdfghjklzxcvbnm'
      carmaj='QWERTYUIOPASDFGHJKLZXCVBNM'
      taille=len(chaine)
      if (taille.le.0) return
      do icar=1,taille
         ilettre=index(carmin,chaine(icar:icar))
         if (ilettre.gt.0) then
             chaine(icar:icar)=carmaj(ilettre:ilettre)
         endif
      enddo
      return
      end
"ENDOF mintomaj.f"
#=======================================================================
echo openam.f
cat >openam.f <<"ENDOF openam.f"
c
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT.
c
c     input:
c        namfil: nom du fichier a ouvrir.
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output:
c        qexist: flag / existence du fichier
c        qinterr: Pas de nom pour le fichier cherche.
c
c     YHS-oct-1993: Premiere version.
c     YHS-jan-2000: Derniere modification.
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*10 cformat, cstatus
      character*64 namfil
c Local
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then
         write(*,*) 'Openam> Program is stopping on user request.'
         stop
      endif
c
c     Checks if filename is consistent with the opening:
c
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          write(6,'(A/A)') '%Openam-Err> File not found. Filename: ',
     .    namfil
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(*,'(/A)')
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c
      if (qverbos) then
         write(*,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(*,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status='UNKNOWN',unit=unit)
c
      return
      end
"ENDOF openam.f"
#=======================================================================
echo openam1.f
cat >openam1.f <<"ENDOF openam1.f"
c----------------------------------------------------------------
      SUBROUTINE openam1(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT,
c     a priori suite a une interrogation...
c
c     input:
c        namfil: nom du fichier a ouvrir.
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output:
c        qexist: flag / existence du fichier
c        qinterr: Pas de nom pour le fichier cherche.
c
c     YHS-oct-93
c     YHS-jan-95
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*(*) namfil, cformat, cstatus
c Local
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then
         write(*,*) 'Openam> Program is stopping on user request.'
         stop
      endif
c
c     Checks if filename is consistent with the opening:
c
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) write(6,'(A)') '%Openam-Err> File not found.'
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(*,'(/A)')
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c
      if (qverbos) then
         write(*,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(*,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)
c
      return
      end
"ENDOF openam1.f"
#=======================================================================
echo openam2.f
cat >openam2.f <<"ENDOF openam2.f"
c-----------------------------------------------------------------------
      subroutine openam2(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT,
c     a priori suite a une interrogation...
c
c     input:
c        namfil: nom du fichier a ouvrir.
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output:
c        qexist: flag / existence du fichier
c        qinterr: Pas de nom pour le fichier cherche.
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*(*) namfil, cformat, cstatus
c Local
      integer lnom
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then
         write(6,'(2A)') 'Openam> Program is stopping on user request.'
         stop
      endif
c     Checks if filename is consistent with the opening:
      call stringcl2(namfil,lnom)
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) then
              write(6,'(/2A)') '%Openam-Err> File: ',namfil(1:lnom)
              write(6,'(A)')
     .      ' Expected in the current directory, but not found.'
          endif
          return
      endif
      if (qexist.and.cstatus.eq.'NEW') then
         write(6,'(/2A)')
     .      '%Openam-Err> This file exists:',namfil(1:lnom)
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)
      if (qverbos) then
         write(6,'(/2A)') ' Openam> File opened: ',namfil(1:lnom)
      endif
      return
      end
"ENDOF openam2.f"
#=======================================================================
echo pdbmat.f
cat >pdbmat.f <<"ENDOF pdbmat.f"
c=====================================================================
c   Pdbmat:
c
c   Computes the mass-weighted second derivatives energy matrix,
c   using Tirion's model, that is, an elastic network model (ENM).
c   In such models, close particles (atoms) are linked by springs.
c
c   In Tirions's model, springs are set between atoms less than
c   CUToff Angstroms away from each others.
c
c   In Hinsen's version, springs are weighted by exp(-dij/rkh)**2.
c   Herein, both kind of typical distances (cutoff & rkh) can be mixed.
c
c   Goal: Computing the low-frequency (collective) normal modes
c   of vibration of the system.
c
c   To do so, the matrix produced by pdbmat has to be diagonalized
c   with a program like DIAGSTD or, for large systems, DIAGRTB,
c   BLZPACK, etc.
c
c=====================================================================
c
c   INPUT:
c   ******
c  -A parameter file named pdbmat.dat.
c
c   Note that each run of pdbmat produces a pdbmat.dat_run file,
c   where parameter values are shown (and shortly commented).
c   pdbmat.dat_run can be modified and used as a pdbmat.dat file,
c   for further runs.
c   So, the simplest is to compile and run pdbmat...
c   and have a look at the pdbmat.dat_run file so produced.
c   If you want to see there the syntax of other available commands,
c   just raise the PRINting level.
c
c  -Among the parameters:
c
c  *The name of a file with the coordinates of the system,
c   in FREE or PDB (protein data bank) format.
c
c   Free format: x, y, z, mass.
c   PDB  format: Only lignes with the ATOM keyword are considered.
c        Masses can be given in the Bfactors column.
c   Note that masses can be read, but they are not required.
c   In the later case, they are all set to 1.0 (as often done).
c
c  *A way to identify pairs of neighbors:
c   Either a CUTOFF value (standard Tirion's model)
c   Or a file with a list of pairs of neighbors.
c   Format of each ligne of this later file: atom-number atom-number
c
c   Alternatively, Hinsen's version can be used.
c
c   For one-atom-per-residue protein models, typical values are:
c   Tirion's distance cutoff = 10-12 Angstroms.
c   Hinsen's typical range   = 3 Angstroms.
c
c   OUTPUT:
c   *******
c
c  -Default output matrix filename:
c   Formatted file: pdbmat.sdijf
c   Binary    file: pdbmat.sdijb
c
c   This is a matrix in format: i, j, non-zero-i-j-matrix-element
c
c  -Output coordinate filename (in free format):
c   pdbmat.xyzm
c
c   This is a coordinate file with, for each atom:
c   x, y, z, mass, block-number
c   For a pdb file, the block-number is the amino-acid residue number
c  (it is of use only in other programs like DIAGRTB).
c
c   More specialized ones:
c   ----------------------
c
c  -VMD-command filename:
c   pdbmat.vmd
c
c   This is a command file for vizualising the elastic network
c   with VMD, the Visual Molecular Dynamics program (v1.8):
c   vmd -e pdbmat.vmd
c   But if you just want to vizualise a standard ENM (defined with
c   a cutoff distance), there are more efficient (faster) ways.
c
c  -Molscript-command-file:
c   pdbmat.molscript
c
c   This is a command file for vizualising the elastic network
c   with Molscript (v2.1).
c   molscript < pdbmat.molscript > pdbmat.ps
c   But if you just want to vizualise a standard ENM (defined with
c   a cutoff distance) with Molscript, there are much faster ways.
c
c.....................................................................
      program pdbmat
      implicit none
      integer natmax, nresmx, nvoismx
c   This is a fortran 77 program, so it has predefined:
c
c   **************
c   MEMORY LIMITS:
c   **************
c
c   NATMAX  :  Maximum number of atoms (particles).
c   NRESMX  :  Maximum number of residues (groups of particles).
c   NVOISMX :  Maximum number of pairs of neighbors.
      parameter( NATMAX=50000 )
      parameter( NRESMX=50000 )
      parameter( NVOISMX=1000000 )
c   Increase them if needed. To (re)compile pdbmat, type:
c   make pdbmat
c   or:
c   g77 -o pdbmat pdbmat.f
c   or use your favorite fortran compiler (instead of g77).
c.....................................................................
c
c   ABOUT Tirion's model:
c   *********************
c
c   Principe du modele (Tirion, 1996):
c
c   Tous les atomes a moins de "cutoff" les uns des autres
c   sont supposes lies par des ressorts, qui ont tous
c   la meme raideur.
c   Simplification supplementaire par rapport au modele initial:
c   les atomes sont supposes avoir tous la meme taille
c  (le cutoff est le meme pour toutes les paires d'atomes).
c   On peut de plus poser qu'ils ont tous la meme masse.
c   Sinon, celles-ci sont lues dans la colonne des
c   facteurs B du fichier pdb.
c
c   Principaux resultats:
c
c   Les modes de vibration de basse frequence obtenus
c   a partir d'un tel modele sont tres voisins de ceux
c   obtenus avec un modele beaucoup plus detaille, tels
c   ceux utilises lors des etudes de Dynamique Moleculaire.
c
c   Dans le cas ou le mouvement fonctionnel d'une proteine est un
c   mouvement d'ensemble (collectif), on constate qu'il peut tres
c   souvent etre decrit comme une combinaison lineaire de quelques
c   uns de ces modes (de un a trois).
c
c   Principaux avantages:
c
c   Pas besoin de prendre en compte tous les atomes.
c   Pas besoin de minimisation d'energie prealablement
c   au calcul des modes de vibration (E=0 par construction).
c
c.....................................................................
c
c   MAIN REFERENCES:
c   ****************
c
c   1) M.M. Tirion (1996):
c  "Large amplitude elastic motions in proteins from
c   a single-parameter, atomic analysis",
c   Phys. Rev. letters vol.77(9), p1905-1908.
c
c   2) K. Hinsen (1998):
c  "Analysis of domain motions by approximate normal mode calculations"
c   Proteins vol.33, p417-429.
c
c   3) F. Tama, Y.H. Sanejouand (2001):
c  "Conformational change of proteins arising
c   from normal modes calculations"
c   Protein Engineering vol.14, p1-6.
c
c   4) A.R. Atilgan, S.R. Durell, R.L. Jernigan, M.C. Demirel,
c   O. Keskin, I. Bahar (2001):
c  "Anisotropy of Fluctuation Dynamics of Proteins with an Elastic
c   Network Model"
c   Biophys Journal vol.80, p.505-515.
c
c   5) S. Nicolay, Y.H. Sanejouand (2006):
c  "Functional modes of proteins are among the most robust"
c   Phys. Rev. letters vol.96, p078104.
c
c   In case of problem, feel free to contact:
c   Yves-Henri.Sanejouand@ens-lyon.fr
c  (bug reports may help you, but also others)
c
c.....................................................................
c     YHS-Nov-1996: Version 1.00.
c     Versions released at http://ecole-modelisation.free.fr/modes.html
c     YHS-Mar-2001: Version 3.31.
c     YHS-Feb-2004: Version 3.50.
c     YHS-Feb-2008: Version 3.73.
c     Version used by the ELNEMO Web site (http://www.elnemo.org):
c     YHS-Feb-2004: Version 3.46.
c
c   New features since the previously released version (v3.50):
c  -The list of pairs of neighbors can be read in a file.
c  -The elastic network can be visualized, using VMD or Molscript.
c  -The output format has been extended.
c  -Hinsen's weight added.
c  -Bond-dihedral energy term added (1-4 bond, dihedral-like).
c  -cumentation included in the code.
c
c.....................................................................
      integer nmotsmax, ntopmax
      parameter(ntopmax=10*natmax,nmotsmax=100)
      integer fatres(nresmx+1),
     .        i, idmax, idres(nresmx), ii, imax, imin, ires,
     .        iresat(natmax), iseed, ivois(nvoismx),
     .        j, jangle(ntopmax), jat, jbond(ntopmax), jdihe(ntopmax),
     .        jj, jsomang(ntopmax), jvois(nvoismx),
     .        k, kcom, kk, klist,
     .        ll, lnom, lnomlst, lnommls, lnommtx, lnompdb, lnomvmd,
     .        namax, namin, nangle(ntopmax), nangles, natom,
     .        nbig, nbmax, nbmin, nbond(ntopmax), nbonds, ndat,
     .        ndihe(ntopmax), ndihs, nl, nmax, nmin, nmots, nntr,
     .        nnzero, nres, nunit, nunknown, nvois, nvoisat(natmax),
     .        prtlev, uninp, unlst, unmol, unout, unpdb, unrsd, unvmd
      double precision cutbnd, cutoff, ddf, der2(3,3*natmax),
     .        dist, dist2, dmax, dmin, dmoy, drms,
     .        elemnt, elmax, fvois(natmax),
     .        kangle, kbond, kdihe, kfce, kij, knonb, kvois(natmax),
     .        levelshft, massat(natmax), nmoy, nrms,
     .        random, rave, rbig, rdev, rinput, rkh, rmax, rmin, rsmall,
     .        rx, ry, rz,
     .        trace, unknown, xat(natmax), yat(natmax), zat(natmax)
      logical qbinary, qerror, qexist, qfread, qinter, qlist,
     .        qmasse, qmtx, qok, qpdb, qvois(natmax)
      character atonam(natmax)*4, cformat*32, csep*1, cstatus*32,
     .        lign80*80, motinp*80, mots(nmotsmax)*132,
     .        nomfich*64, nomlst*64, nommls*64, nommtx*64, nompdb*64,
     .        nomvmd*64,
     .        program*8, progrer*11, progrwn*11,
     .        residus_standards*132, residus_stshort*21,
     .        resnam(natmax)*4, segid(natmax)*4, ssunam(natmax)*1,
     .        ssusel*1, typbond*80, typmas*80, typout*80, version*32
      parameter(rbig=1e10,rsmall=1e-10,unknown=9999.d9)
c.......................................................................
      version=' Version 3.73, February 2008.'
      unrsd=0
      nangles=0
      ndihs=0
c.......................................................................
      idmax=21
      residus_standards='   ILE PHE TRP LEU CYS VAL MET TYR ALA HIS '//
     .                     'GLY THR SER PRO ARG GLN ASN ASP GLU LYS '
      residus_stshort='IFWLCVMYAHGTSPRQNDEKX'
      program=' Pdbmat>'
      progrer='%Pdbmat-Er>'
      progrwn='%Pdbmat-Wn>'
      write(6,'(2A)') program,
     .' Computes the Hessian matrix, using an Elastic Network Model.'
      write(6,'(2A)') program,version
c     ===============================================
c     Ouverture et Lecture du fichier d'instructions:
c     ===============================================
c     Default values:
      nompdb='pdbmat.ent'
      cutoff=10.d0
c     Neighbor list:
      nomlst='NONE'
      qlist=.false.
c     Typical distance for Hinsen's weigth:
c    (negative value means 1/rkh=0)
      rkh=-1.d0
      qfread=.false.
      nommls='NONE'
      nomvmd='NONE'
      typmas='CONS'
      qmasse=.false.
      knonb=1.0d0
c     Topological terms (like in ref.6):
      typbond='NONE'
c     Distance-cutoff value for defining covalent (chemical) bonds:
      cutbnd=4.0d0
      kangle=0.0d0
      kbond=1000.0d0
      kdihe=0.0d0
c     Others:
      typout='  FREE'
      qmtx=.false.
      qbinary=.false.
      prtlev=0
      levelshft=1e-8
c     Used only with a levelshift, to add some noise (not useful ?).
      iseed=27041961
      nunit=10
      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat'
      cformat="FORMATTED"
      cstatus="old"
      call openam2(nomfich,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)
      if (qinter.or..not.qexist) then
          write(6,'(/2A/(2A))') progrwn,
     .  ' No pdbmat.dat command file found.',
     .    progrwn,' Defaults assumed for all options. ',
     .    progrwn,' See the pdbmat.dat_run file if you need an example.'
          goto 110
      else
          write(6,'(/2A)') program,
     .  ' Options to be read in pdbmat.dat file.'
      endif
 50   continue
      read(uninp,'(A)',end=100) lign80
      kcom=index(lign80,'!')
      k=index(lign80,'=')
      motinp=' '
      if (k.gt.0.and.(kcom.le.0.or.kcom.gt.k)) then
          motinp=lign80(1:k)
      else
          if (k.le.0.and.kcom.gt.1) then
          write(6,'(/2A/A)') progrwn,
     .  ' No separator (=) in command ligne:',
     .    lign80
          write(6,'(2A)') progrwn,' This ligne is skipped.'
          endif
          goto 50
      endif
      call mintomaj(motinp)
      kcom=index(lign80(k+1:80),'!')
      if (kcom.gt.0) lign80(k+kcom:80)=' '
      klist=index(lign80(k+1:80),'?')
      if (index(motinp,' FILENAME').gt.0.or.
     .    index(motinp,' SCRIPT').gt.0) then
          if (index(motinp,'MATRI').gt.0) then
              nommtx=lign80(k+1:80)
              qmtx=.true.
          elseif (index(motinp,'LIST ').gt.0.or.
     .            index(motinp,'NEIGH').gt.0) then
              nomlst=lign80(k+1:80)
              qlist=.true.
          elseif (index(motinp,'MOLS').gt.0) then
              nommls=lign80(k+1:80)
          elseif (index(motinp,'VMD').gt.0) then
              nomvmd=lign80(k+1:80)
          else
              nompdb=lign80(k+1:80)
          endif
      else if (index(motinp,' DEFINITION').gt.0) then
          typbond=lign80(k+1:80)
          call mintomaj(typbond)
          call stringcl2(typbond,lnom)
          if (typbond(1:3).eq.'ALL') then
              typbond=' ALL'
          else if (typbond(1:3).eq.'NON') then
              typbond='NONE'
          else if (typbond(1:3).eq.'CON') then
              typbond='CONSECUTIF'
          else
              write(6,'(/3A)') progrwn,' Bond definition :',
     .        typbond(1:4)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: NONe, ALL, CONsecutive.'
              write(6,'(A)') ' Default assumed.'
              typbond='NONE'
          endif
      else if (index(motinp,'MASS').gt.0) then
          typmas=lign80(k+1:80)
          call mintomaj(typmas)
          call stringcl2(typmas,lnom)
          if (typmas(1:3).eq.'PDB'.or.typmas(1:3).eq.'COO') then
              qmasse=.true.
              typmas='COOR'
          else if (typmas(1:3).ne.'CON') then
              write(6,'(/3A)') progrwn,' Origin of mass values :',
     .        typmas(1:3)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: CONstant, COOr, PDB.'
              write(6,'(A)') ' Default assumed.'
              qmasse=.false.
              typmas='CONS'
          endif
      else if (index(motinp,'FORMAT').gt.0) then
          typout=lign80(k+1:80)
          call mintomaj(typout)
          call stringcl2(typout,lnom)
          if (typout(1:1).eq.'B'.or.typout(1:1).eq.'U') then
              qbinary=.true.
              typout='BINARY'
          else if (typout(1:1).ne.'F') then
              write(6,'(/3A)') progrwn,' Kind of matrix format :',
     .        typout(1:1)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: Free, Binary, Formatted, Unformatted.'
              write(6,'(A)') ' Default assumed.'
              qbinary=.false.
              typout='  FREE'
          else
              qbinary=.false.
              typout='  FREE'
          endif
      else
          qok=.false.
          read(lign80(k+1:80),*,end=90,err=90) rinput
          if (index(motinp,'SHIFT ').gt.0) then
               qok=.true.
               levelshft=rinput
          else if (index(motinp,'CUTOF').gt.0.or.
     .             index(motinp,'DISTANCE').gt.0) then
               qok=.true.
               cutoff=rinput
          else if (index(motinp,'HINSEN').gt.0) then
               qok=.true.
               rkh=rinput
          else if (index(motinp,'INTERAC').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   knonb=rinput
               endif
               if (index(motinp,' CUTOF').gt.0.or.
     .             index(motinp,' DIST').gt.0) then
                   qok=.true.
                   cutoff=rinput
               endif
          else if (index(motinp,'BOND').gt.0.and.
     .        (index(motinp,' FORCE ').gt.0.or.
     .         index(motinp,' CONST').gt.0)) then
               qok=.true.
               kbond=rinput
          else if (index(motinp,' LENGTH').gt.0) then
               qok=.true.
               cutbnd=rinput
          else if (index(motinp,'PRINT').gt.0) then
               qok=.true.
               prtlev=int(rinput)
          else if (index(motinp,'ANGLE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kangle=rinput
               endif
          else if (index(motinp,'DIHE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kdihe=rinput
               endif
          endif
  90      continue
          if (.not.qok) then
               write(6,'(/2A/A)') progrwn,
     .       ' No known or incomplete set of keywords in ligne:',
     .         motinp
               write(6,'(2A)') progrwn,
     .       ' This command ligne is skipped.'
          endif
      endif
      goto 50
 100  continue
      close(uninp)
 110  continue
      call stringcl2(nompdb,lnompdb)
      call stringcl2(nomlst,lnomlst)
      call stringcl2(nommls,lnommls)
      call stringcl2(nomvmd,lnomvmd)
      if (nomlst.eq.'none'.or.nomlst.eq.'NONE') qlist=.false.
      if (nommls.eq.'none') nommls='NONE'
      if (nomvmd.eq.'none') nomvmd='NONE'
      if (.not.qmtx) then
      if (qbinary) then
        nommtx="pdbmat.sdijb"
      else
        nommtx="pdbmat.sdijf"
      endif
      endif
      call stringcl2(nommtx,lnommtx)
c     Resume des commandes:
c     ---------------------
      write(6,'(/3A)') program,' Coordinate filename     = ',
     .      nompdb(1:lnompdb)
      if (qlist) then
      write(6,'(3A)') program,' Neighbor-list filename  = ',
     .      nomlst(1:lnomlst)
      else
      write(6,'(/2A,F10.2)') program,
     .        ' Distance cutoff         = ',cutoff
      endif
      if (rkh.gt.0.d0)
     .write(6,'(8X,A,F10.2)') " Hinsen's typical range  = ",rkh
      write(6,'(A,F10.2)')
     .'         Force constant          = ',knonb
      if (typbond.ne.'NONE') then
      write(6,'(A,6X,A)')
     .'         Kind of bond definition = ',typbond(1:4)
      write(6,'(A,F10.2)')
     .'         Maximum bond length     = ',cutbnd,
     .'         Bond force constant     = ',kbond,
     .'         Angle force constant    = ',kangle,
     .'         Dihedral force constant = ',kdihe
      endif
      write(6,'(A,6X,A)')
     .'         Origin of mass values   = ',typmas(1:4)
      write(6,'(3A)') program,' Matrix filename         = ',
     .      nommtx(1:lnommtx)
      if (prtlev.gt.0) then
      write(6,'(2A,1PG10.1)') program,
     .        ' Levelshift              = ',levelshft
      write(6,'(A,3X,I7)')
     .'         PRINTing level          = ',prtlev
      endif
      if (nommls.ne.'NONE')
     .write(6,'(3A)') program,' Molscript filename      = ',
     .      nommls(1:lnommls)
      if (nomvmd.ne.'NONE')
     .write(6,'(3A)') program,' VMD script filename     = ',
     .      nomvmd(1:lnomvmd)
c     Sauvegarde du fichier de commandes complet:
c     -------------------------------------------
      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat_run'
      cformat="FORMATTED"
      cstatus="ove"
      call openam2(nomfich,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)
c     Pour etre plus clair:
      if (typbond.eq.'NONE') then
          cutbnd=0.d0
          kbond=0.d0
          kangle=0.d0
          kdihe=0.d0
      endif
      if (qlist) cutoff=-1
      write(uninp,'(2A)')
     .'! This file can be modified and used as a command file',
     .' (named pdbmat.dat) for pdbmat.'
      write(uninp,'(2A)') ' Coordinate FILENAME        = ',
     .      nompdb(1:lnompdb)
      write(uninp,'(2A)') ' MATRIx FILENAME            = ',
     .      nommtx(1:lnommtx)
      write(uninp,'(A,F10.3,A)') ' INTERACtion DISTance CUTOF = ',
     .      cutoff,' ! For defining the list of interacting atoms.'
      write(uninp,'(A,F10.3,A)') ' INTERACtion FORCE CONStant = ',
     .      knonb,' ! For specifying frequency units.'
      if (prtlev.gt.0.or.rkh.gt.0.d0)
     .write(uninp,'(A,F10.3,A)') " HINSEN's typical range     = ",
     .      rkh,' ! Force constant weighting (if negative: none).'
      if (prtlev.gt.0.or.nomlst(1:lnomlst).ne.'NONE')
     .write(uninp,'(3A)') ' NEIGHbor-list FILENAME     = ',
     .      nomlst(1:lnomlst),
     .  ' ! For defining this list yourself.'
      write(uninp,'(A,6X,2A)') ' Origin of MASS values      = ',
     .      typmas(1:4),' ! CONstant, or from COOrdinate file.'
      write(uninp,'(A,8X,I2,A)') ' Output PRINTing level      = ',
     .      prtlev,' ! =1: more detailled. =2: debug level.'
c     Not often used:
      if (prtlev.gt.0.or.nommls(1:lnommls).ne.'NONE')
     .write(uninp,'(3A)') ' MOLScript command FILEname = ',
     .      nommls(1:lnommls),
     .  ' ! To draw the network with Molscript.'
      if (prtlev.gt.0.or.nomvmd(1:lnomvmd).ne.'NONE')
     .write(uninp,'(3A)') ' VMD command FILEname       = ',
     .      nomvmd(1:lnomvmd),
     .  ' ! vmd -e this-file (to visualize the network with VMD).'
c     Rarely used:
      if (prtlev.gt.0.or.typbond.ne.'NONE') then
      write(uninp,'(A,6X,2A)') ' Bond DEFINITION            = ',
     .      typbond(1:4),' ! NONe, ALL, or between CONsecutive atoms.'
      write(uninp,'(A,F10.3)') ' Maximum bond LENGTH        = ',cutbnd
      write(uninp,'(A,F10.3)') ' BOND FORCE CONStant        = ',kbond
      write(uninp,'(A,F10.3)') ' ANGLE FORCE CONStant       = ',kangle
      write(uninp,'(A,F10.3)') ' DIHEdral FORCE CONStant    = ',kdihe
      write(uninp,'(A,1PG10.1,A)') ' LevelSHIFT                 = ',
     .      levelshft,
     .  ' ! Non-zero value often required (numerical reasons).'
      write(uninp,'(A,4X,2A)') ' Matrix FORMAT              = ',
     .      typout(1:6),' ! Free, or Binary, matrix saved.'
      endif
      close(uninp)
c     Tests:
      if ((.not.qlist.and.cutoff.lt.0.d0.and.rkh.lt.0.d0).or.
     .    knonb.lt.0.d0.or.
     .   (typbond(1:4).ne.'NONE'.and.(cutbnd.lt.0.d0.or.kbond.lt.0.d0
     .   .or.kangle.lt.0.d0.or.kdihe.lt.0.d0))) then
          write(6,'(/2A)') progrer,
     .  ' Distances and force constants can not have negative values !'
          stop '*Commands are not consistent*'
      endif
c     On recherche l'information/sous-unite:
      call string_split(nompdb,lnompdb,":",
     .                  mots,nmotsmax,nmots)
      call stringcl2(mots(1),lnom)
      if (nmots.gt.1) then
          call stringcl2(mots(2),lnom)
          ssusel=mots(nmots)(1:1)
          write(6,'(3A)') program,' Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,' The end of filename, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)(1:64)
      call stringcl2(nompdb,lnompdb)
c
c     Lecture du fichier de coordonnees:
c     ==================================
      if (prtlev.gt.0)
     .  write(6,'(/(4A))') program,
     .' Coordinate file ',nompdb(1:lnompdb),' to be opened.'
      unpdb=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam2(nompdb,cformat,cstatus,unpdb,.true.,
     .     qinter,qexist)
      if (qinter) stop '*No readable coordinate file found*'
c     Format pdb ?
      nl=0
      qpdb=.false.
 120  continue
      read(unpdb,'(A)',end=130) lign80
      if (lign80(1:5).eq.'ATOM '.or.lign80(1:6).eq.'HETATM ') then
          qpdb=.true.
          goto 130
      else
          nl=nl+1
      endif
      goto 120
 130  continue
      rewind(unpdb)
      do i=1,natmax
         xat(i)=unknown
         yat(i)=unknown
         zat(i)=unknown
         massat(i)=unknown
         iresat(i)=i
      enddo
      if (qpdb) then
          write(6,'(/2A)') program,
     .  ' Coordinate file in PDB format.'
          call rdatompdbp(unpdb,ssusel,xat,yat,zat,massat,
     .         atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .         fatres,nresmx,nres,qerror,prtlev)
          if (natom.eq.nres) then
              write(6,'(/2A)') program,' Study of a standard ENM model.'
          else
              write(6,'(/2A)') progrwn,
     .      ' Study of a several-atom-per-residue model '//
     .      '(this is not that standard).'
          endif
      else
          if (nl.eq.0) then
              write(6,'(/2A)') progrer,' Empty coordinate file.'
              stop
          endif
          write(6,'(/2A)') program,
     .  ' Coordinate file in Free format.'
          call readxyz(unpdb,xat,yat,zat,massat,iresat,natmax,natom,
     .         ndat,qerror,prtlev)
          if (qmasse.and.ndat.lt.4) then
              write(6,'(/2A)') progrer,
     .      ' Masses were not all found, as expected.'
              qmasse=.false.
          endif
      endif
c     Tests:
c     ======
      if (qerror) stop
      if (natom.le.1) then
          write(6,'(2A)') progrer,
     .  ' Not enough atoms found in file. Nothing done.'
          stop
      endif
      write(6,'(/2A)') program,' Coordinate statistics: '
      call vecstat(xat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))')
     .' <x>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax
      call vecstat(yat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))')
     .' <y>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax
      call vecstat(zat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))')
     .' <z>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax
      if (qmasse) then
          write(6,'(/2A)') program,' Mass statistics: '
          call vecstat(massat,natom,rmin,rmax,rave,rdev)
          write(6,'(4(A,F12.6))')
     .  ' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax
          if (rmin.le.0.d0) then
              write(6,'(2A)') progrer,
     .      ' Negative or null masses found !'
              qmasse=.false.
          endif
      endif
      if (.not.qmasse) then
          write(6,'(2A)') program,' Masses are all set to one.'
          do i=1,natom
             massat(i)=1.d0
          enddo
      endif
c     Test/identification des residus.
      if (qpdb) then
      nunknown=0
      do i=1,nres
         ires=fatres(i)
         idres(i)=index(residus_standards,resnam(ires))/4
         if (idres(i).le.0) then
             nunknown=nunknown+1
             if (nunknown.lt.10) then
                 write(6,'(4A)') progrwn," residue:'",
     .           resnam(ires),"' is not a well known amino-acid."
                 idres(i)=idmax
             else if (nunknown.eq.10) then
                 write(6,'(2A)') progrwn,' ........'
                 idres(i)=idmax
             endif
         endif
      enddo
      if (nunknown.gt.0)
     .write(6,'(/A,I6,A)') progrwn,nunknown,' residue(s) not known.'
      endif
c   Covalent bond detection: through a distance-cutoff criterium.
c     Bonds i-j and j-i are stored, because below the matrix is
c     calculated and saved three lines at a time.
      nbonds=0
      if (typbond.ne.' ALL'.and.typbond.ne.'CONSECUTIF') goto 200
      if (typbond.eq.' ALL') then
         nbmax=0
         nbmin=999
         imax=-1
         imin=-1
         k=1
         nbond(1)=1
         do i=1,natom
            do j=1,natom
               if (i.ne.j) then
                   rx=xat(i)-xat(j)
                   ry=yat(i)-yat(j)
                   rz=zat(i)-zat(j)
                   dist=dsqrt(rx*rx + ry*ry + rz*rz)
                   if (dist.le.cutbnd) then
                       jbond(k)=j
                       k=k+1
                   endif
               endif
            enddo
            nbond(i+1)=k
            if (nbond(i+1)-nbond(i).gt.nbmax) then
                nbmax=nbond(i+1)-nbond(i)
                imax=i
            endif
            if (nbond(i+1)-nbond(i).lt.nbmin) then
                nbmin=nbond(i+1)-nbond(i)
                imin=i
            endif
            if (k-1.gt.ntopmax) then
                write(6,'(/2A,I12)') progrer,
     .        ' Too many bonds. Maximum is: ',ntopmax
                stop
            endif
         enddo
         nbonds=k-1
      else if (typbond.eq.'CONSECUTIF') then
c      Even when the CONsecutif keyword is used.
c        On fait attention aux distances...
c        Il peut y avoir plusieurs molecules,
c        plusieurs chaines, dans le systeme.
         k=1
         do i=1,natom
            nbond(i)=k
            if (i.gt.1) then
            j=i-1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            endif
            if (i.lt.natom) then
            j=i+1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            endif
            if (k.gt.ntopmax) then
                write(6,'(/2A,I12)') progrer,
     .        ' Too many bonds. Maximum is: ',ntopmax
                stop
            endif
         enddo
         nbond(natom+1)=k
         imax=2
         imin=1
         nbmin=1
         nbmax=2
         nbonds=k
      endif
      if (nbonds.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond found.'
          goto 200
      endif
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nbonds/2,' covalent bonds, i.e.,',
     .  float(nbonds)/float(2*natom),' per atom.'
      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)')
     .'         Maximum number found =',nbmax,' for atom ',imax,
     .'         Minimum number found =',nbmin,' for atom ',imin
c   Bond-ANGLes are allowed only if there are BONDs.
      nangles=0
      if (kangle.le.0.d0) then
          write(6,'(/2A)') program,' BOND but no ANGLe energy term.'
          goto 200
      endif
      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      nangle(1)=1
      do i=1,natom
         if (nbond(i+1).gt.nbond(i)) then
         do jj=nbond(i),nbond(i+1)-1
            j=jbond(jj)
            if (nbond(j+1).gt.nbond(j)) then
            do kk=nbond(j),nbond(j+1)-1
               k=jbond(kk)
               if (k.ne.i) then
                   jangle(ii)=k
                   jsomang(ii)=j
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         nangle(i+1)=ii
         if (nangle(i+1)-nangle(i).gt.namax) then
             namax=nangle(i+1)-nangle(i)
             imax=i
         endif
         if (nangle(i+1)-nangle(i).lt.namin) then
             namin=nangle(i+1)-nangle(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/2A,I12)') progrer,
     .     ' Too many angles. Maximum is: ',ntopmax
             stop
         endif
      enddo
      nangles=ii-1
      if (nangles.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond-angle found.'
          goto 200
      endif
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nangles/2,' valence angles, i.e.,',
     .  float(nangles)/float(2*natom),' per atom.'
      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)')
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin
c   Bond-DIHEdrals are allowed only if there are BONDs and ANGLes.
      ndihs=0
      if (kdihe.le.0.d0) then
         write(6,'(/2A)') program,' ANGLes but no DIHEdral energy term.'
         goto 200
      endif
      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      ndihe(1)=1
c     For each atom:
      do i=1,natom
c        All angles where it is the first atom:
         if (nangle(i+1).gt.nangle(i)) then
         do jj=nangle(i),nangle(i+1)-1
c           For the "top-atom" of this angle:
            j=jsomang(jj)
c           All angles where it is the first atom:
            if (nangle(j+1).gt.nangle(j)) then
            do kk=nangle(j),nangle(j+1)-1
               k=jangle(kk)
               if (k.ne.i.and.j.ne.jsomang(kk).and.i.ne.jsomang(kk))then
                   jdihe(ii)=k
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         ndihe(i+1)=ii
         if (ndihe(i+1)-ndihe(i).gt.namax) then
             namax=ndihe(i+1)-ndihe(i)
             imax=i
         endif
         if (ndihe(i+1)-ndihe(i).lt.namin) then
             namin=ndihe(i+1)-ndihe(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/2A,I12)') progrer,
     .     ' Too many dihes. Maximum is: ',ntopmax
             stop
         endif
      enddo
      ndihs=ii-1
      if (ndihs.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond-dihedral found.'
          goto 200
      endif
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  ndihs/2,' dihedral angles, i.e.,',
     .  float(ndihs)/float(2*natom),' per atom.'
      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)')
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin
c     ==============================
c     Matrice des derivees secondes:
c     ==============================
 200  continue
c     Lecture eventuelle de la liste des voisins:
c     -------------------------------------------
c     Formats possibles:
c     atom-number atom-number
c     atom-number atom-number force-constant
      if (qlist) then
      write(6,'(/2A)') program,' Neighbor list to be read.'
      unlst=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam2(nomlst,cformat,cstatus,unlst,.true.,
     .     qinter,qexist)
      nvois=0
      ndat=0
 210  continue
      read(unlst,'(A)',end=230,err=220) lign80
c     Les commentaires ne sont pas pris en compte:
      kcom=index(lign80,'!')
      if (kcom.le.0) kcom=index(lign80,'#')
      if (kcom.gt.0) then
          if (kcom.eq.1) then
              write(6,'(2A)') lign80(1:50),' ...'
              goto 210
          else
              lign80=lign80(1:kcom-1)
          endif
      endif
c     Plusieurs separateurs sont possibles: , ; ou blanc
      csep=","
      k=index(lign80,csep)
      if (k.le.0) then
          csep=";"
          k=index(lign80,csep)
      endif
      if (k.le.0) csep=" "
      call string_split(lign80,80,csep,
     .     mots,nmotsmax,nmots)
      if (nmots.lt.2) goto 225
      if (ndat.eq.0) then
          ndat=min(nmots,3)
          if (ndat.eq.3) qfread=.true.
      endif
      if (nmots.lt.ndat) then
          write(6,'(/A,I3,A/A)') progrer,ndat,
     .  ' data per ligne until ligne: ',lign80
          goto 220
      endif
c     Lecture de: i, j (kij le cas echeant).
      read(mots(1),*) ii
      read(mots(2),*) jj
      if (qfread) read(mots(3),*) kfce
      if (ii.le.0.or.jj.le.0) then
          write(6,'(/2A)') progrer,
     .        ' Null or negative atom number found.'
          goto 225
      endif
      if (ii.gt.natom.or.jj.gt.natom) then
          write(6,'(/2A,I6,A,I6,A)') progrer,
     .  ' Atom number: ',ii,' or ',jj,
     .  ' larger than the number of atoms.'
          stop '*Wrong file*'
      endif
      nvois=nvois+1
      if (nvois.gt.nvoismx) then
          write(6,'(/2A,I6,A)') progrer,' More than ',nvoismx,
     .  ' pairs of neighbors, the maximum allowed. Sorry.'
          stop '*Recompile with larger array*'
      endif
      ivois(nvois)=ii
      jvois(nvois)=jj
      if (qfread) fvois(nvois)=kfce
c     Ligne suivante:
      goto 210
c     Probleme de lecture:
 220  continue
      write(6,'(2A,I6)') progrer,
     .    ' While reading neigbors pair: ',nvois+1
      stop '*Wrong or corrupted file*'
 225  continue
      write(6,'(2A/A)') progrer,
     .' No (or wrong) pair of atom numbers found in ligne: ',lign80
      stop '*Wrong or corrupted file*'
c     Fin du fichier:
 230  continue
      write(6,'(/A,I6,A)') program,nvois,' pairs of neighbors.'
      endif
c     Coordonnees et masses utilisees, sauvegardees:
c     ----------------------------------------------
      unout=nunit
      nunit=nunit+1
      nomfich="pdbmat.xyzm"
      cformat="FORMATTED"
      cstatus="ove"
      call openam2(nomfich,cformat,cstatus,unout,.true.,
     .     qinter,qexist)
      do i=1,natom
         write(unout,'(4(1PG20.12),I9)')
     .   xat(i), yat(i), zat(i), massat(i), iresat(i)
      enddo
      close(unout)
      if (prtlev.gt.0)
     .write(6,'(2A)') program,
     .    ' Coordinates and masses considered are saved.'
c     Fichier de commandes pour VMD:
c     ------------------------------
      unvmd=-1
      if (nomvmd.ne.'NONE') then
      unvmd=nunit
      nunit=nunit+1
      nomfich=nomvmd
      cformat="FORMATTED"
      cstatus="ove"
      call openam2(nomfich,cformat,cstatus,unvmd,.true.,
     .     qinter,qexist)
      write(unvmd,'(A)') '#!/usr/local/bin/vmd'
      write(unvmd,'(A)') '# script for VMD (Visual Molecular Dynamics)'
      write(unvmd,'(A)') '# Goal: visualizing the elastic network'
      write(unvmd,'(A)') '# Type: vmd -e this-file'
      write(unvmd,'(A)') 'color Display {Background} white'
      write(unvmd,'(A)') 'mol new'
      write(unvmd,'(A)') 'draw color black'
      endif
c     Fichier de commandes pour Molscript:
c     ------------------------------------
      unmol=-1
      if (nommls.ne.'NONE') then
      unmol=nunit
      nunit=nunit+1
      nomfich=nommls
      cformat="FORMATTED"
      cstatus="ove"
      call openam2(nomfich,cformat,cstatus,unmol,.true.,
     .     qinter,qexist)
      write(unmol,'(A)') '! Script for Molscript (Kraulis, 1993)'
      write(unmol,'(A)') '! Goal: visualizing the elastic network'
      write(unmol,'(A)')  ' set bonddistance 99.0 ;'
      endif
c     Matrice:
c     --------
      if (qbinary) then
        if (.not.qmtx) nommtx="pdbmat.sdijb"
        cformat="UNFORMATTED"
      else
        if (.not.qmtx) nommtx="pdbmat.sdijf"
        cformat="FORMATTED"
      endif
      unout=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam2(nommtx,cformat,cstatus,unout,.true.,
     .     qinter,qexist)
c     ========================================
c     Les atomes sont tous lies deux a deux,
c     par un potentiel "universel" (M.Tirion).
c     ========================================
      elmax=0.d0
      trace=0.d0
      dmin=0.d0
      dmax=0.d0
      dmoy=0.d0
      drms=0.d0
      nnzero=0
      nntr=0
      nbig=0
      ll=0
      do i=1,natom
         ii=3*i-2
         nvoisat(i)=0
c        Liste eventuelle des voisins de i:
c        ----------------------------------
         if (qlist) then
             do j=1,natom
                qvois(j)=.false.
                kvois(j)=0.d0
             enddo
             do j=1,nvois
                if (ivois(j).eq.i) then
                    qvois(jvois(j))=.true.
                    if (qfread) kvois(jvois(j))=fvois(j)
                endif
                if (jvois(j).eq.i) then
                    qvois(ivois(j))=.true.
                    if (qfread) kvois(ivois(j))=fvois(j)
                endif
             enddo
         endif
c        On calcule trois lignes de la matrice a la fois:
c        -----------------------------------------------
         do j=1,3*natom
            der2(1,j)=0.d0
            der2(2,j)=0.d0
            der2(3,j)=0.d0
         enddo
         do j=1,natom
            if (.not.qlist.or.(qlist.and.qvois(j))) then
            if (i.ne.j) then
            jj=3*j-2
            kij=knonb
            if (qfread) kij=kvois(j)
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist2=rx*rx + ry*ry + rz*rz
            dist=dsqrt(dist2)
            if (dist.lt.rsmall) then
                write(6,'(/2A,1PG10.4,A/2(I6,2A,I6,A,1X,2A))')
     .          progrer,' Too small distance = ',dist,
     .        ' between following atoms.',
     .          i,': ',resnam(i),iresat(i),ssunam(i),atonam(i),' and ',
     .          j,': ',resnam(j),iresat(j),ssunam(j),atonam(j)
                stop '*Wrong coordinates*'
            endif
c         Hinsen's version: force constant * exp(-(dist/rkh)**2.d0)
c        (except for topological bond force constants, e.g. ANGLes).
            if (rkh.gt.0.d0) kij=kij*exp(-(dist/rkh)**2.d0)
c           Constantes de force topologiques:
            if (nbonds.gt.0) then
            if (nbond(i+1).gt.nbond(i)) then
                do k=nbond(i),nbond(i+1)-1
                   if (jbond(k).eq.j) then
                       kij=kbond
                       goto 300
                   endif
                enddo
            endif
            else
            goto 300
            endif
            if (nangles.gt.0) then
            if (nangle(i+1).gt.nangle(i)) then
                do k=nangle(i),nangle(i+1)-1
                   if (jangle(k).eq.j) then
                       kij=kangle
                       goto 300
                   endif
                enddo
            endif
            else
            goto 300
            endif
            if (ndihs.gt.0) then
            if (ndihe(i+1).gt.ndihe(i)) then
                do k=ndihe(i),ndihe(i+1)-1
                   if (jdihe(k).eq.j) then
                       kij=kdihe
                       goto 300
                   endif
                enddo
            endif
            else
            goto 300
            endif
 300        continue
c           Calcul des elements: (potentiel harmonique)
c           -------------------------------------------
            if (dist.le.cutoff.or.qlist.or.
     .         (cutoff.le.0.d0.and.rkh.gt.0.d0)) then
                ll=ll+1
                nvoisat(i)=nvoisat(i)+1
                if (j.gt.i) then
                   if (unvmd.gt.0) then
                   write(unvmd,'(A,3F12.4,A,3F12.4,A)') 'draw line {',
     .             xat(i),yat(i),zat(i),'} {',xat(j),yat(j),zat(j),'}'
                   endif
                   if (unmol.gt.0) then
                   write(unmol,'(A,I6,3A)')
     .           ' bonds require in residue ',iresat(i),', atom ',
     .             atonam(i),' and in molecule mol1 ',
     .           ' require in residue ',iresat(j),', atom ',
     .             atonam(j),' and in molecule mol1 ; '
                   endif
c                  Potentiel harmonique: 1/eval*knonb*(d - rval)**eval
                endif
                if (ll.eq.1.or.dist.lt.dmin) dmin=dist
                if (ll.eq.1.or.dist.gt.dmax) dmax=dist
                dmoy=dmoy+dist
                drms=drms+dist2
c               Elements diagonaux des blocs i et j:
c               -----------------------------------
                ddf=kij/dist2
                elemnt=rx*rx*ddf
                der2(1,ii)=der2(1,ii)+elemnt
                der2(1,jj)=der2(1,jj)-elemnt
                elemnt=ry*ry*ddf
                der2(2,ii+1)=der2(2,ii+1)+elemnt
                der2(2,jj+1)=der2(2,jj+1)-elemnt
                elemnt=rz*rz*ddf
                der2(3,ii+2)=der2(3,ii+2)+elemnt
                der2(3,jj+2)=der2(3,jj+2)-elemnt
c               Elements extra-diagonaux des deux blocs:
c               ---------------------------------------
                elemnt=rx*ry*ddf
                der2(1,ii+1)=der2(1,ii+1)+elemnt
                der2(2,ii)=der2(2,ii)+elemnt
                der2(1,jj+1)=der2(1,jj+1)-elemnt
                der2(2,jj)=der2(2,jj)-elemnt
                elemnt=rx*rz*ddf
                der2(1,ii+2)=der2(1,ii+2)+elemnt
                der2(3,ii)=der2(3,ii)+elemnt
                der2(1,jj+2)=der2(1,jj+2)-elemnt
                der2(3,jj)=der2(3,jj)-elemnt
                elemnt=ry*rz*ddf
                der2(2,ii+2)=der2(2,ii+2)+elemnt
                der2(3,ii+1)=der2(3,ii+1)+elemnt
                der2(2,jj+2)=der2(2,jj+2)-elemnt
                der2(3,jj+1)=der2(3,jj+1)-elemnt
            endif
            endif
            endif
         enddo
c        Sortie de la matrice-bande calculee:
c        -----------------------------------
c       (Uniquement la demi-matrice superieure)
c        Level-shift, pour eviter les zeros numeriques,
c        lors de la diagonalisation a venir
c       (la minimisation est parfaite, par definition).
c        Le hasard est la pour lever la degenerescence
c        des six valeurs propres nulles, et differencier
c        rotations et translations.
         der2(1,ii)  =der2(1,ii)   + levelshft*random(iseed)
         der2(2,ii+1)=der2(2,ii+1) + levelshft*random(iseed)
         der2(3,ii+2)=der2(3,ii+2) + levelshft*random(iseed)
         do j=ii,3*natom
            jat=(j-1)/3+1
            if (der2(1,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout)
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)')
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(1,j)).gt.rbig)  nbig=nbig+1
                if (dabs(der2(1,j)).gt.elmax) elmax=dabs(der2(1,j))
            endif
         enddo
         do j=ii+1,3*natom
            jat=(j-1)/3+1
            if (der2(2,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout)
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)')
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(2,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(2,j)).gt.elmax) elmax=dabs(der2(2,j))
            endif
         enddo
         do j=ii+2,3*natom
            jat=(j-1)/3+1
            if (der2(3,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout)
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)')
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(3,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(3,j)).gt.elmax) elmax=dabs(der2(3,j))
            endif
         enddo
         elemnt=(der2(1,ii)+der2(2,ii+1)+der2(3,ii+2))/massat(i)
         if (elemnt.eq.0.d0) then
             write(6,'(2A,I6,A)') progrwn,
     .     ' Atom ',i,' has a null second derivatives...'
         else
             nntr=nntr+1
         endif
         trace=trace+elemnt
      enddo
      close(unout)
      if (unvmd.gt.0) then
          write(unvmd,'(2A)') 'mol load pdb ',nompdb(1:lnompdb)
          close(unvmd)
      endif
      if (unmol.gt.0) close(unmol)
      if (unrsd.gt.0) close(unrsd)
      nmoy=0.d0
      nrms=0.d0
      nmin=natom
      nmax=0
      do i=1,natom
         if (nvoisat(i).gt.nmax) nmax=nvoisat(i)
         if (nvoisat(i).lt.nmin) nmin=nvoisat(i)
         nmoy=nmoy+nvoisat(i)
         nrms=nrms+nvoisat(i)**2.d0
      enddo
      nmoy=nmoy/float(natom)
      nrms=nrms/float(natom)-nmoy**2.d0
      if (nrms.gt.0.d0) nrms=dsqrt(nrms)
      if (ll.eq.0) then
          write(6,'(/2A,I12,A)') progrer,
     .  ' No atom-atom interaction found. Too short cutoff ?'
          stop '*Empty matrix*'
      else
          dmoy=dmoy/float(ll)
          drms=drms/float(ll)-dmoy**2.d0
          if (drms.gt.0.d0) drms=dsqrt(drms)
      endif
      if (prtlev.gt.0)
     .write(6,'(/2A)') program,' Matrix statistics:'
      write(6,'(/2A,F8.4,A)') program,' The matrix is ',
     .  100.d0*dfloat(nnzero)/dfloat(3*natom*(3*natom+1)/2),' % Filled.'
      write(6,'(A,I12,A)') program,nnzero,'  non-zero elements.'
      if (prtlev.gt.0) then
      write(6,'(A,I12,A)') program,ll/2,' atom-atom interactions.'
      write(6,'(/2A,F9.2,A,F9.2/(A,I6))') program,
     .        ' Number per atom= ',nmoy,' +/- ',nrms,
     .'         Maximum number = ',nmax,
     .'         Minimum number = ',nmin
      write(6,'(/2A,F9.2,A,F9.2/(A,F9.2))') program,
     .        ' Average dist.  = ',dmoy,' +/- ',drms,
     .'         Maximum dist.  = ',dmax,
     .'         Minimum dist.  = ',dmin
      endif
      write(6,'(/2A,1PG12.6)') program,' Matrix trace   = ',trace
      if (prtlev.gt.0) then
      write(6,'(2A,1PG12.6)') program,' Larger element = ',elmax
      write(6,'(A,I6,A,1PG8.1)') program,
     .      nbig,' elements larger than +/- ',rbig
      endif
      write(6,'(/2A)') program,' Hessian matrix ready.'
      write(6,'(2A)') program,
     .' To diagonalize it (and get the modes),'//
     .' you may use diagstd, blzpack, diagrtb...'
      write(6,'(/2A)') program,' Normal end.'
      stop
      end
"ENDOF pdbmat.f"
#=======================================================================
echo projmod.f
cat >projmod.f <<"ENDOF projmod.f"
      program projmod
c     Projection of a difference-vector onto a set of vectors,
c     and/or (qdepl=T) displacement of the structure along one of them.
c     Difference-vector: difference of two coordinates files (pdb format)
c     Vectors: CERFACS format.
c     YHS-Fev-1997: Premiere version (Toulouse).
c     YHS-Avr-2003: Derniere version (Lyon).
c     A ajouter:
c     La collectivite de chaque mode; celle du changement de conf.
      implicit none
      integer nvecmx, natmax, nresmx, nmotsmax
      parameter(natmax=10000,nresmx=10000,nvecmx=26,nmotsmax=132)
      integer fatres(nresmx+1), i, ii, iresat(natmax),
     .        iresatc(natmax), ivec, j, k, lmot, lnomeig,
     .        lnompdb, lnompdbc, lnompdbr, mdmax, modnum, natdif, natom,
     .        natomc, adm,
     .        natomeff, nddl, nddleff, nmots, nres, nresc, nrsdif,
     .        numvec(natmax), nunit, nvec, nzero, atequiv(natmax),
     .        uneig, unmor, unout, unpdb, unpdbc, unpdbq, unpdbr
      double precision avemass, bigzero, dq, dr2,
     .       freq(3*natmax), freqmin, lowfreq, massat(natmax),
     .       matvec(3*natmax,nvecmx), norme, qdiff, qnorm, qp, qproj,
     .       qtot(nvecmx),recmax, rmsat, rmsmass, small, tot, w(natmax),
     .       xat(natmax), xref(natmax), xconf(natmax), xcurr(natmax),
     .       yat(natmax), yref(natmax), yconf(natmax), ycurr(natmax),
     .       zat(natmax), zref(natmax), zconf(natmax), zcurr(natmax),
     .       wmain(natmax),x,y,z,dist,dmax,xc,yc,zc
      logical qdepl, qdir, qerror, qexist,  qinter, qmasse, qok
      character atonam(natmax)*4, atonamc(natmax)*4, cformat*10,
     .       cstatus*10, mots(nmotsmax)*132, namfil*64, nom*7,
     .       nomeig*64, nompdb*64, nompdbc*64, nompdbr*64, program*9,
     .       progrer*12, progrwn*12, resnam(natmax)*4,
     .       resnamc(natmax)*4, segid(natmax)*4, ssunam(natmax)*4,
     .       ssusel*1, version*40
cDefault:
      version=' Version 1.36, April 2003.'
      small=1e-4
      nom='Projmod'
      program=' '//nom//'>'
      progrwn='%'//nom//'-Wn>'
      progrer='%'//nom//'-Er>'
cBegin:
      write(6,'(2A)') program,version
      write(6,'(2A)') program,
     .    ' Projection of a difference vector'//
     .    ' on a set of eigenvectors.'
c     Ouverture des fichiers:
c     ----------------------
c     En lecture:
      call getnam('CERFACS file with the eigenvectors ?',
     .     nomeig,lnomeig,qok)
      if (.not.qok) stop
      cformat="FORMATTED"
      cstatus="old"
      nunit=10
      uneig=nunit
      nunit=nunit+1
      call openam(nomeig,cformat,cstatus,uneig,.true.,
     .            qinter,qexist)
      if (qinter.or..not.qexist) stop
c     On recherche l'information/sous-unite:
c     Structure de reference:
      call getnam('Pdb file with the reference structure ?',
     .     nompdbr,lnompdbr,qok)
      if (.not.qok) stop
      nompdb=nompdbr
      lnompdb=lnompdbr
      call string_split(nompdb,lnompdb,':',
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lmot)
      if (nmots.gt.1) then
          call stringcl(mots(2),lmot)
          ssusel=mots(nmots)(1:1)
          write(6,*) 'Pdbsel> Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,
     .      ' The end of pdb name, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)(1:64)
      unpdb=nunit
      nunit=nunit+1
      call openam(nompdb,cformat,cstatus,unpdb,.true.,
     .            qinter,qexist)
      if (qinter.or..not.qexist) stop
      unpdbr=unpdb
c     Autre conformere (eventuel):
      call getnam('Pdb file with the other conformer ?',
     .     nompdbc,lnompdbc,qok)
      if (.not.qok) stop
      nompdb=nompdbc
      lnompdb=lnompdbc
      call string_split(nompdb,lnompdb,':',
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lmot)
      if (nmots.gt.1) then
          call stringcl(mots(2),lmot)
          ssusel=mots(nmots)(1:1)
          write(6,*) 'Pdbsel> Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,
     .      ' The end of pdb name, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)(1:64)
      unpdbc=-1
      inquire(file=nompdb,exist=qexist)
      if (qexist) then
      unpdb=nunit
      nunit=nunit+1
      call openam(nompdb,cformat,cstatus,unpdb,.true.,
     .            qinter,qexist)
      if (qinter) stop
      unpdbc=unpdb
      endif
      call getrep(
     .   ' Are the masses given in the pdb file ? (y/n)',
     .     qmasse,qok)
      if (.not.qok) qmasse=.false.
      if (qmasse) then
          write(6,'(/2A)') program,
     .  ' Masses will be picked in the pdb files.'
      else
          write(6,'(/2A)') program,
     .  ' All masses will all be assumed to be of 1.'
      endif
      call getrep(
     .    'Displacement along one mode ? (y/n)',
     .     qdepl,qok)
      if (.not.qok) qdepl=.false.
      modnum=-1
      qdir=.true.
      if (qdepl) then
          call getchi('Dq ?',dq,qok)
          if (.not.qok) then
              write(6,'(/2A)') progrer,
     .      ' While reading dq. No displacement.'
              qdepl=.false.
          endif
          call getnum('Mode number ?',modnum,1,-1,qok)
          if (.not.qok.or.modnum.le.0) then
            if (unpdbc.gt.0) then
              write(6,'(/2A)') progrwn,
     .      ' Wrong mode number. Displacement along best mode.'
            else
              write(6,'(/2A)') progrer,
     .      ' Wrong mode number. No displacement.'
              stop
            endif
            modnum=-1
          endif
          if (modnum.lt.0) then
          call getrep(
     .        'Displacement along difference-vector direction ? (y/n)',
     .         qdir,qok)
          if (.not.qok) qdir=.true.
          endif
      endif
c     En ecriture:
      if (unpdbc.gt.0) then
      cformat="FORMATTED"
      cstatus="ove"
      namfil='projmod.res'
      unout=nunit
      nunit=nunit+1
      call openam(namfil,cformat,cstatus,unout,.true.,
     .            qinter,qexist)
      namfil='dr.res'
      unmor=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam(namfil,cformat,cstatus,unmor,.true.,
     .            qinter,qexist)
      endif
      if (qdepl) then
cg      namfil='projmod_dq.pdb'
      namfil='fort.12'
      unpdbq=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam(namfil,cformat,cstatus,unpdbq,.true.,
     .            qinter,qexist)
      endif
c     Lecture des fichiers:
c     --------------------
      call rdmodfacs(uneig,3*natmax,nvecmx,numvec,freq,
     .               matvec,nddl,nvec)
      write(6,'(/A,I5,A,I6,A)')
     .' Projmod> ',nvec,' vectors, ',nddl,' coordinates in file.'
      call rdatompdb(unpdbr,ssusel,xref,yref,zref,massat,
     .     atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .     fatres,nresmx,nres,qerror)
      if (nddl.eq.3*natom) then
          write(6,'(/2A)') program,
     .  ' Cartesian (eigen)vectors will be studied.'
      else if (nddl.gt.0) then
          write(6,'(/2A)') progrer,
     .  ' Vector and reference pdb file are not consistent.'
          stop
      else if (nddl.le.0.or.nvec.le.0) then
          write(6,'(/2A)') progrer,
     .  ' Nothing can be done.'
          stop
      endif
      if (.not.qmasse) then
      do i=1,natom
          massat(i)=1.d0
      enddo
      endif
      if (unpdbc.le.0) then
          write(6,'(/2A)') program,' One conformer only.'
          if (qdepl) then
              goto 500
          else
              write(6,'(/2A)') progrer,' Nothing done.'
              stop
          endif
      endif
      call rdatompdb(unpdbc,ssusel,xconf,yconf,zconf,wmain,
     .     atonamc,iresatc,resnamc,ssunam,segid,natmax,natomc,
     .     fatres,nresmx,nresc,qerror)
c     Conformite des deux conformeres consideres:
      if (natomc.ne.natom.or.nresc.ne.nres) then
      if (natomc.ne.natom) then
          write(6,'(/2A)') progrwn,
     .  ' Different number of atoms for the other conformer.'
      endif
      if (nresc.ne.nres) then
          write(6,'(/2A)') progrwn,
     .  ' Different number of residues for the other conformer.'
      endif
      if (natomc.lt.natom.or.nresc.lt.nres) then
          write(6,'(2A)') progrer,' Not enough of them.'
          stop
      endif
      do i=1,natom
         do j=1,natomc
         if (iresatc(j).eq.iresat(i)) then
         if (atonamc(j).eq.atonam(i).and.resnamc(j).eq.resnam(i)) then
             atequiv(i)=j
             xcurr(i)=xconf(j)
             ycurr(i)=yconf(j)
             zcurr(i)=zconf(j)
             goto 100
         endif
         endif
         enddo
         write(6,'(2A,I6,2A,I6,1X,2A)') progrer,
     . ' Atom ',i,' of first conformer: ',
     .   resnam(i),iresat(i),atonam(i),' not found in second one.'
         stop
 100  continue
      enddo
      do i=1,natom
         xconf(i)=xcurr(i)
         yconf(i)=ycurr(i)
         zconf(i)=zcurr(i)
         atonamc(i)=atonam(i)
         resnamc(i)=resnam(i)
         iresatc(i)=iresat(i)
      enddo
      write(6,'(2A)') program,
     .  ' All atoms of first conformer were found in second one.'
      endif
c     Vecteur difference:
c     ------------------
      natdif=0
      nrsdif=0
      do i=1,natom
c        Un des atomes est inconnu ?
         if ((xref(i).gt.9998.and.yref(i).gt.9998.and.
     .        zref(i).gt.9998).or.
     .       (xconf(i).gt.9998.and.yconf(i).gt.9998.and.
     .        zconf(i).gt.9998)) then
              massat(i)=0.d0
         endif
         xat(i)=xconf(i)-xref(i)
         yat(i)=yconf(i)-yref(i)
         zat(i)=zconf(i)-zref(i)
         if (atonamc(i).ne.atonam(i)) then
             natdif=natdif+1
             if (natdif.lt.5) then
                 write(6,'(2A,I6,5A)') progrwn,' Atom ',i,
     .         ' has name ',atonamc(i),' in a file and ',
     .           atonam(i),' in the other.'
             elseif (natdif.eq.5) then
                 write(6,'(2A)') progrwn,' ... '
             endif
         endif
         if (resnamc(i).ne.resnam(i)) then
             nrsdif=nrsdif+1
             if (nrsdif.lt.5) then
                 write(6,'(2A,I6,5A)') progrwn,' Atom ',i,
     .         ' belongs to residue ',resnamc(i),' in a file and ',
     .           resnam(i),' in the other.'
             elseif (nrsdif.eq.5) then
                 write(6,'(2A)') progrwn,' ... '
             endif
         endif
      enddo
      if (natdif.gt.0)
     .   write(6,'(/A,I6,A)') progrwn,
     .   natdif,' atoms have different names in the two pdb files.'
      if (nrsdif.gt.0)
     .   write(6,'(/A,I6,A)') progrwn,nrsdif,
     . ' atoms belong to different residues in the two pdb files.'
      write(6,'(/2A/)') program,
     . ' File dr.res: displacement=f(atom number).'
      avemass=0.d0
      rmsmass=0.d0
      rmsat=0.d0
      natomeff=0
      do i=1,natom
         dr2=0.d0
         if (massat(i).gt.small) then
             avemass=avemass+massat(i)
             rmsmass=rmsmass+massat(i)**2.d0
             dr2=xat(i)**2.d0+yat(i)**2.d0+zat(i)**2.d0
             rmsat=rmsat+dr2
             natomeff=natomeff+1
         endif
         write(unmor,'(I6,F12.4)') i,dsqrt(dr2)
      enddo
      write(6,'(A,I6,A)') program,natomeff,' atoms are considered.'
      nddleff=3*natomeff
      if (natomeff.eq.0) then
          write(6,'(/2A)') progrer,
     .  ' All atoms have unknown coordinate(s) or have zero mass.'
          stop
      endif
      avemass=avemass/float(natomeff)
      rmsmass=dsqrt(rmsmass/float(natomeff)-avemass**2.d0)
      rmsat=dsqrt(rmsat)
      write(6,'(/2A,F8.2)') program,
     .    ' Atomic r.m.s. displacements=  ',
     .      rmsat/sqrt(float(natomeff))
      write(6,'(A,2(A,F8.2))') program,
     .    ' Atomic average masses      =  ',avemass,' +/- ',rmsmass
      if (dabs(avemass-1.).gt.small.or.rmsmass.gt.small) then
      norme=0.d0
      do i=1,natom
         if (massat(i).gt.small)
     .   norme=norme+
     .  (xat(i)**2.d0+yat(i)**2.d0+zat(i)**2.d0)*massat(i)
      enddo
      if (norme.le.small) then
         write(6,'(2A,I6,A)') progrer,
     . ' Difference-vector has null norm. Projection skipped.'
         goto 500
      endif
      norme=dsqrt(norme)
      write(6,'(2A,F8.2)') program,
     .    ' Atomic mass-weighted rmsd  =  ',
     .      norme/sqrt(float(natomeff))
      else
      norme=rmsat
      endif
      write(6,'(/2A/)') program,' File projmod.res:'//
     .    ' dr.vector=f(fqcy), and cumulative square sum.'
      mdmax=-1
      recmax=0.d0
      tot=0.d0
      do ivec=1,nvec
         qnorm=0.d0
         qproj=0.d0
         do i=1,natom
          if (massat(i).gt.small) then
            ii=3*i-2
            qp=matvec(ii,ivec)*xat(i)*dsqrt(massat(i))+
     .         matvec(ii+1,ivec)*yat(i)*dsqrt(massat(i))+
     .         matvec(ii+2,ivec)*zat(i)*dsqrt(massat(i))
            qproj=qproj+qp
            qnorm=qnorm+
     .         matvec(ii,ivec)**2.d0+
     .         matvec(ii+1,ivec)**2.d0+
     .         matvec(ii+2,ivec)**2.d0
          endif
         enddo
c        Cosinus:
c        -------
         if (qnorm.le.small) then
             write(6,'(2A,I6,A)') progrwn,' Eigenvector ',ivec,
     .           ' has null norm. It is skipped.'
             goto 400
         endif
         if ((natom.eq.natomeff.and.dabs(qnorm-1.).gt.small).or.
     .       qnorm-1.gt.small) then
             write(6,'(2A,I6,A,F12.4)') progrwn,' Eigenvector ',ivec,
     .           ' Norm= ',qnorm
         endif
         qnorm=dsqrt(qnorm)
         qdiff=qproj/qnorm
         qproj=qdiff/norme
         if (dabs(qproj).gt.dabs(recmax)) then
             recmax=qproj
             mdmax=ivec
         endif
c        Somme cumulee des cosinus-carres:
         qtot(ivec)=qproj**2.d0
         tot=tot+qtot(ivec)
         write(6,'(A,I6,A,F8.2,2(4X,A,F6.3),A,F12.4)')
     .   ' Vector: ',ivec,' F= ',freq(ivec),
     .   ' Cos= ',qproj,' Sum= ',tot,
     .   ' q= ',qdiff
         write(unout,'(I6,3F12.4)') ivec,freq(ivec),qproj,tot
c        Vecteur suivant:
 400     continue
      enddo
c     Compter les modes a frequence "nulle":
c     a) Frequence la plus proche de la nullite.
c     b) Frequences "proches" de cette valeur (a un facteur 100 pres).
      freqmin=2.e+38
      do i=1,nvec
         if (dabs(freq(i)).lt.freqmin) freqmin=dabs(freq(i))
      enddo
      write(6,'(/2A,F13.6)') program,
     .' Best zero-frequency found  :',freqmin
      nzero=0
      bigzero=freqmin
      lowfreq=2.e+38
      do i=1,nvec
         if (dabs(freq(i)).lt.100*freqmin) then
             if (dabs(freq(i)).gt.bigzero) bigzero=dabs(freq(i))
             nzero=nzero+1
         else
             if (dabs(freq(i)).lt.lowfreq) lowfreq=dabs(freq(i))
         endif
      enddo
      if (nzero.eq.nvec) then
          write(6,'(2A/)') progrwn,
     .    ' ...but this is not a null frequency.'
          nzero=0
          bigzero=freqmin
          lowfreq=freqmin
      endif
      write(6,'(A,I6,A,F12.6)') program,nzero,
     . ' frequencies less than: ',bigzero
      write(6,'(2A,F13.6)') program,
     .' Lowest non-zero frequency  :',lowfreq
      write(6,'(/2A,F6.2,A,I6,A,F8.2,A)') program,
     .  ' Best overlap with diff.vect. = ',recmax,' for mode ',mdmax,
     .'   with F= ',freq(mdmax),' cm-1.'
c     Tri des cosinus-carres par ordre decroissant:
      call trier(qtot,nvec,nvecmx,w,.false.)
c     Sommation:
      tot=0.d0
      do i=1,nvec
         tot=tot+w(i)
         qtot(i)=tot
      enddo
      if (nvec.gt.nzero+60) then
         write(6,'(/2A,6(F5.3,2X))')  program,
     . ' 1-3-6-9-12-all-best contrb. = ',
     .   qtot(1),qtot(3),qtot(6),qtot(9),
     .   qtot(12),qtot(nvec)
      else if (nvec.gt.nzero+9) then
         write(6,'(/2A,5(F5.3,2X))')  program,
     . ' 1-3-6-9-all-best contrb. = ',
     .   qtot(1),qtot(3),qtot(6),qtot(9),
     .   qtot(nvec)
      else if (nvec.gt.nzero+3) then
         write(6,'(/2A,5(F5.3,2X))')  program,
     . ' 1-3-all-best contributions = ',
     .   qtot(1),qtot(3),qtot(nvec)
      else
         write(6,'(/2A,1(F5.3,2X))')  program,
     . ' All known contributions = ',qtot(nvec)
      endif
c     Displacement along one mode:
c     ----------------------------
 500  continue
      if (qdepl) then
      if (modnum.le.0) modnum=mdmax
      if (modnum.gt.nvec) then
          write(6,'(/2A,I6)') progrwn,
     .  ' Displacement can not be performed along mode ',modnum
          modnum=nvec
      endif
c     Direction:
      if ((recmax.gt.0.and..not.qdir).or.(recmax.lt.0.and.qdir)) then
          write(6,'(/2A)') progrwn,' Displacement reversed: '//
     .   'along vector-diff. direction, as required.'
          dq=-dq
      endif
      write(6,'(/2A,F12.6,A,I6)') program,
     .  ' Displacement of dq= ',dq,' to be performed along mode ',modnum
      k=0
      do i=1,natom
         if (massat(i).gt.small) then
            ii=3*i-2
            xat(i)=xref(i)+(dq*matvec(ii,modnum))/dsqrt(massat(i))
            yat(i)=yref(i)+(dq*matvec(ii+1,modnum))/dsqrt(massat(i))
            zat(i)=zref(i)+(dq*matvec(ii+2,modnum))/dsqrt(massat(i))
         else
            k=k+1
            if (k.le.5)
     .      write(6,'(2A,I6,A)') progrer,' Atom ',k,
     .    ' has unknown coordinates or zero mass.'
            if (k.eq.5)
     .      write(6,'(2A)') progrer,' ... '
         endif
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(6,'(/a)') 'Col1 : Mode number'
      write(6,'(a)') 'Col2 : Maximal displacement'
      write(6,'(a/)') 'Col3 : Most implied atome'
      do modnum=7,nvecmx
      dmax=0.
      xc=0.
      yc=0.
      zc=0.
      adm=0
      do i=1,natom
            ii=3*(i-1)
            x=matvec(ii+1,modnum)
            y=matvec(ii+2,modnum)
            z=matvec(ii+3,modnum)
            xc=xc+x
            yc=yc+y
            zc=zc+z
            dist=x**2+y**2+z**2
            if (dist.ge.dmax) then
            dmax = dist
            adm=i
            endif
      enddo
      write(6,'(i5,f10.5,i5)')
     .  modnum,sqrt(dmax),adm
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (k.eq.0) then
      call writpdb(unpdbq,xat,yat,zat,massat,
     .     atonam,iresat,resnam,segid,natom)
      else
      write(6,'(/A,I6,A)') progrer,k,
     .   ' atoms have unknown coordinates or zero mass.'
      stop
      endif
      endif
      write(6,'(/2A)') program,' Normal end.'
      stop
      end
"ENDOF projmod.f"
#=======================================================================
echo random.f
cat >random.f <<"ENDOF random.f"
C-----------------------------------------------------------------------
      REAL*8 FUNCTION RANDOM(ISEED)
C-----------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR: UNIFORM DISTRIBUTION (0,1)
C     ISEED: SEED FOR GENERATOR. ON THE FIRST CALL THIS HAS TO
C     HAVE A VALUE IN THE EXCLUSIVE RANGE (1, 2147483647)
C     AND WILL BE REPLACED BY A NEW VALUE TO BE USED IN
C     FOLLOWING CALL.
C
C     REF: Lewis, P.A.W., Goodman, A.S. & Miller, J.M. (1969)
C     "Pseudo-random number generator for the System/360", IBM
C     Systems Journal 8, 136.
C
C     This is a "high-quality" machine independent generator.
C     INTEGERS are supposed to be 32 bits or more.
C     The same algorithm is used as the basic IMSL generator.
C
C     Author: Lennart Nilsson
C
      implicit none
      INTEGER ISEED
      REAL*8 DSEED,DIVIS,DENOM,MULTIP
      DATA  DIVIS/2147483647.D0/
      DATA  DENOM /2147483711.D0/
      DATA  MULTIP/16807.D0/
C
      IF(ISEED.LE.1) ISEED=314159
      DSEED=MULTIP*ISEED
      DSEED=MOD(DSEED,DIVIS)
      RANDOM=DSEED/DENOM
      ISEED=DSEED
C
      RETURN
      END
"ENDOF random.f"
#=======================================================================
echo rdatompdb.f
cat >rdatompdb.f <<"ENDOF rdatompdb.f"
c-----------------------------------------------------------------------
      subroutine rdatompdb(unpdb,ssusel,xat,yat,zat,binfo,
     .           atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .           fatres,nresmx,nres,qerror)
c
c     Lecture ligne a ligne d'un fichier pdb.
c     Uniquement les lignes commencant par 'ATOM'.
c     Uniquement ceux de la sous-unite selectionnee.
c
c     fatres(i): numero du premier atome du residu i.
c
c     YHS-nov-1996: premiere version (Toulouse).
c     YHS-mar-2003: derniere version (Lyon).
c
      implicit none
cI/O:
      integer unpdb, natmax, iresat(*), natom,  lnom,
     .        nresmx, nres, fatres(*)
      double precision xat(*), yat(*), zat(*), binfo(*)
      logical qerror
      character*4 atonam(*), resnam(*), segid(*)
      character*1 ssusel, ssunam(*)
cLocal:
      integer nerr, iatom, irs, irsprev
      double precision x, y, z, bfact
      character*1  ssu
      character*4  ren, segat
      character*5  atncur
      character*80 lign80
cBegin:
      write(6,'(/A)') ' Rdatompdb> Reading pdb file.'
c
      qerror=.false.
      nerr=0
c
      irsprev=-1
      nres=0
      iatom=0
 105  continue
      read(unpdb,'(A)',end=200,err=110) lign80
c
      goto 120
 110  continue
      nerr=nerr+1
c
 120  continue
      if (lign80(1:4).eq.'ATOM') then
      read(lign80,'(12X,A4,1X,A4,A1,I4,4X,3F8.3,6X,F6.2,6X,A4)')
     .            atncur, ren, ssu, irs, x, y, z,
     .            bfact, segat
      if (iatom.lt.natmax) then
          if (ssu.eq.ssusel.or.ssusel.eq.' ') then
          iatom=iatom+1
          xat(iatom)=x
          yat(iatom)=y
          zat(iatom)=z
          binfo(iatom)=bfact
c
          call stringcl(atncur,lnom)
          atonam(iatom)=atncur(1:4)
          call stringcl(ren,lnom)
          resnam(iatom)=ren
          iresat(iatom)=irs
          ssunam(iatom)=ssu
          segid(iatom)=segat
c
          if (irs.ne.irsprev) then
              nres=nres+1
              if (nres.gt.nresmx) then
                  write(6,'(A/A,I6)')
     .          '%Rdatompdb-Er> Too many residues in this file.',
     .          ' Maximum allowed is = ',nresmx
                  stop
              endif
              irsprev=irs
              fatres(nres)=iatom
          endif
          endif
      else
          write(6,'(A/A,I6)')
     .      '%Rdatompdb-Er> Too many atoms in this file.',
     .      ' Maximum allowed is = ',natmax
          stop
      endif
      endif
c
c     2) Ligne suivante du fichier pdb :
c
      goto 105
c
c     3) Fin de la lecture du fichier pdb :
c
 200  continue
      write(6,*) 'Rdatompdb> End of file reached.'
      write(6,*) 'Rdatompdb> Number of I/O errors: ',
     .            nerr
c
      natom=iatom
      fatres(nres+1)=natom+1
      irs=0
      if (natom.gt.0) irs=iresat(natom)
c
      write(6,'(/(A,I6))')
     .' Rdatompdb> Number of residues found = ',nres,
     .'            First residue number     = ',iresat(1),
     .'            Last  residue number     = ',irs,
     .'            Number of atoms found    = ',natom
      write(6,'(A,F8.1)')
     .'            Mean number per residue  = ',float(natom)/float(nres)
c
      if (natom.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No atom found in file.'
          qerror=.true.
      endif
      if (nres.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No residue found in file.'
          qerror=.true.
      endif
c
      return
      end
"ENDOF rdatompdb.f"
#=======================================================================
echo rdatompdbp.f
cat >rdatompdbp.f <<"ENDOF rdatompdbp.f"
c-----------------------------------------------------------------------
      subroutine rdatompdbp(unpdb,ssusel,xat,yat,zat,binfo,
     .           atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .           fatres,nresmx,nres,qerror,prtlev)
c     Lecture ligne a ligne d'un fichier pdb.
c     Uniquement les lignes commencant par 'ATOM'.
c     Uniquement ceux de la sous-unite selectionnee.
c     fatres(i): numero du premier atome du residu i.
c     YHS-nov-1996: version 1.0 (Toulouse).
c     YHS-sep-2004: version 1.1 (Lyon).
      implicit none
cI/O:
      integer unpdb, natmax, iresat(*), natom,  lnom,
     .        nresmx, nres, fatres(*), prtlev
      double precision xat(*), yat(*), zat(*), binfo(*)
      logical qerror
      character*4 atonam(*), resnam(*), segid(*)
      character*1 ssusel, ssunam(*)
cLocal:
      integer iatom, irs, irsprev, nerr, ntit
      double precision bfact, x, y, z
      character*1  ssu
      character*4  ren, segat
      character*5  atncur
      character*80 lign80
cBegin:
      if (prtlev.gt.0)
     .write(6,'(/A)') ' Rdatompdb> Reading pdb file.'
      qerror=.false.
      nerr=0
      irsprev=-1
      ntit=0
      nres=0
      iatom=0
 105  continue
      read(unpdb,'(A)',end=200,err=110) lign80
      goto 120
 110  continue
      nerr=nerr+1
 120  continue
      if (lign80(1:4).eq.'ATOM') then
      read(lign80,'(12X,A4,1X,A4,A1,I4,4X,3F8.3,6X,F6.2,6X,A4)',
     .            end=130,err=130)
     .            atncur, ren, ssu, irs, x, y, z,
     .            bfact, segat
 130  continue
      if (iatom.lt.natmax) then
          if (ssu.eq.ssusel.or.ssusel.eq.' ') then
          iatom=iatom+1
          xat(iatom)=x
          yat(iatom)=y
          zat(iatom)=z
          binfo(iatom)=bfact
          call stringcl2(atncur,lnom)
          atonam(iatom)=atncur(1:4)
          call stringcl2(ren,lnom)
          resnam(iatom)=ren
          iresat(iatom)=irs
          ssunam(iatom)=ssu
          segid(iatom)=segat
          if (irs.ne.irsprev) then
              nres=nres+1
              if (nres.gt.nresmx) then
                  write(6,'(A/A,I6)')
     .          '%Rdatompdb-Er> Too many residues in this file.',
     .          ' Maximum allowed is = ',nresmx
                  stop
              endif
              irsprev=irs
              fatres(nres)=iatom
          endif
          endif
      else
          write(6,'(A/A,I6)')
     .      '%Rdatompdb-Er> Too many atoms in this file.',
     .      ' Maximum allowed is = ',natmax
          stop
      endif
      else if (lign80(1:6).eq.'REMARK'.and.prtlev.gt.0) then
          ntit=ntit+1
          if (ntit.le.10) then
              write(6,'(A)') lign80
          else if (ntit.eq.11) then
              write(6,'(A)') ' .... '
          endif
      endif
c     2) Ligne suivante du fichier pdb :
      goto 105
c     3) Fin de la lecture du fichier pdb :
 200  continue
      if (prtlev.gt.1) then
      write(6,*) 'Rdatompdb> End of file reached.'
      write(6,*) 'Rdatompdb> Number of I/O errors: ',
     .            nerr
      endif
      natom=iatom
      fatres(nres+1)=natom+1
      irs=0
      if (natom.gt.0) irs=iresat(natom)
      write(6,'(/(A,I6))')
     .' Rdatompdb> Number of residues found = ',nres,
     .'            First residue number     = ',iresat(1),
     .'            Last  residue number     = ',irs,
     .'            Number of atoms found    = ',natom
      if (prtlev.gt.0)
     .write(6,'(A,F8.1)')
     .'            Mean number per residue  = ',float(natom)/float(nres)
      if (natom.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No atom found in file.'
          qerror=.true.
      endif
      if (nres.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No residue found in file.'
          qerror=.true.
      endif
      return
      end
"ENDOF rdatompdbp.f"
#=======================================================================
echo rdmodfacs.f
cat >rdmodfacs.f <<"ENDOF rdmodfacs.f"
c
      subroutine rdmodfacs(uneig,nddlmax,nvecmx,numvec,freq,
     .           matvec,nddl,nvec)
c
c     Lecture de modes CERFACS.
c     Devra remplacer rdcerfacs.
c     Difference: comptage de l'ordre de la matrice.
c    (et pas des atomes)
c
c     Premieres versions (rdcerfacs):
c     YHS-Nov-1996.
c     Dernieres modifications:
c     YHS-Jan-2001.
cI/O:
      integer numvec(*), nvecmx, nddlmax, nvec, nddl, uneig
      double precision freq(*), matvec(nddlmax,*)
cLocal:
      integer nmotsmax
      parameter(nmotsmax=100)
      integer nerr, ivec, indnm_cfacs, nmots,
     .        i, ii, k
      double precision wtofreq
      logical qfound, qold, qfirst
      character*1 carnum
      character*132 lign132, mots(nmotsmax)
cDefaut:
c     Facteur de conversion (2*pi*f)**2 -> f (cm-1):
      wtofreq=108.586698
c
      nerr=0
      qfirst=.true.
      qold=.false.
      qfound=.false.
 100  continue
      read (uneig,'(A)',end=300,err=110) lign132
      goto 120
 110  continue
      nerr=nerr+1
 120  continue
c
      qfound=qfound.or.
     .      (index(lign132,' value ').gt.0.and.
     .       index(lign132,' vector ').gt.0.and.
     .       index(lign132,' residual ').le.0)
      qold=qold.or.
     .      (index(lign132,' VALUE ').gt.0.and.
     .       index(lign132,' VECTOR ').gt.0)
c
      if (.not.qfound.and..not.qold) goto 100
c________________________________________
c
c     Lecture des frequences des modes :
c________________________________________
c
      if (qfirst) then
          if (qold) then
          write(6,'(/A)')
     .  ' Rdmodfacs> Old Blzpack file format detected.'
          else
          write(6,'(/A)')
     .  ' Rdmodfacs> Blzpack file format detected.'
          endif
          qfirst=.false.
      endif
c
      ivec=0
      nvec=0
 250  continue
      ivec=ivec+1
      if (ivec.gt.nvecmx) then
          write(6,'(/A,I5,A)')
     .  '%Rdmodfacs-Warning> More than ',nvecmx,' vectors in file.'
          return
      endif
c
      read(lign132,'(7X,I5,12X,G12.4)',end=240,err=240)
     .     numvec(ivec), freq(ivec)
      freq(ivec)=wtofreq*dsqrt(abs(freq(ivec)))
c
      goto 255
 240  continue
      write(6,'(/3A)')
     .    '%Rdmodfacs-W> Pb with ligne: ',lign132(1:36),'...'
 255  continue
c
      nvec=ivec
c--------nettoyage de l affichage-----------
c      write(6,'(/A,I6)')
c     .    ' Rdmodfacs> Numero du vecteur CERFACS en lecture:',
c     .      numvec(ivec)
c      write(6,'(A,1PG12.4)')
c     .    ' Rdmodfacs> Frequence du vecteur en lecture:',
c     .      freq(ivec)
c
      if (numvec(ivec).le.0)
     .    write(6,'(/A/A)')
     .  '%Rdmodfacs-W> Vector number was expected in:',
     .    lign132
c
      read(uneig,'(A)',end=230,err=230) lign132
 230  continue
      read(lign132,'(1X,A1)',end=232,err=232) carnum
 232  continue
      if ((qfound.and.carnum.ne.'=').or.
     .    (qold.and.carnum.ne.'-')) then
          write(6,'(2A/A)')
     .       ' %Rdmodfacs-Warning> Unexpected character ',
     .       ' in second column of line:',
     .    lign132
      endif
c____________________________________________________
c
c     2) Lecture des coordonnees des modes CERFACS :
c        Format libre.
c____________________________________________________
c
      k=0
 257  continue
      if (k.ge.nddlmax) then
          write(6,'(/A,I6,A,I5)')
     .  '%Rdmodfacs-Err> More than ',nddlmax,
     .  ' coordinates for vector ',ivec
          return
      endif
c
      read(uneig,'(A)',end=300,err=270) lign132
c
c     Nombre de coordonnees par ligne:
      call string_split(lign132,132,' ',
     .                  mots,nmotsmax,nmots)
c
      if (lign132.eq.' ') then
          read(uneig,'(A)',end=300,err=260) lign132
      else if (.not.qold.or.index(lign132,' VALUE ').le.0) then
          read(lign132,*,end=258)
     .   (matvec(k+ii,ivec),ii=1,nmots)
          k=k+nmots
          goto 257
 258      continue
      endif
      nddl=k
c
 260  continue
      indnm_cfacs=index(lign132,'       VALUE')
      if (indnm_cfacs.le.0)
     .    indnm_cfacs=index(lign132,'       value')
      if (indnm_cfacs.gt.0) then
          goto 250
      else
          write(6,'(A,A/A/A)')
     .  ' Rdmodfacs: Lecture des modes CERFACS terminee.',
     .  ' Item VALUE non trouve dans la ligne:',lign132
          goto 300
      endif
c
 270  continue
      write(6,'(A,I6,A)')
     .   ' %Rdmodfacs-Error: durant la lecture de la coordonnee ',
     .      i,' du mode.'
      stop
c
      continue
c*****Ligne suivante de la lecture du fichier des modes en cours :
c
      goto 100
c
c     Fin de la lecture du fichier des modes :
c
 300  continue
      return
      end
"ENDOF rdmodfacs.f"
#=======================================================================
echo readxyz.f
cat >readxyz.f <<"ENDOF readxyz.f"
c-----------------------------------------------------------------------
      subroutine readxyz(uninp,x,y,z,w,ic,nmax,ncoor,ndat,qerror,prtlev)
c     Reads at most NMAX coordinates in free format.
c     Either:
c     x, y, z
c     or:
c     x, y, z, w
c     or:
c     x, y, z, w, ic
c     If first word in ligne is not a number, the whole ligne is
c     assumed to be a title or a commentary.
c     YHS-Sep-03: First version (Lyon).
cI/O:
      logical qerror
      integer ic(*), ncoor, ndat, nmax, prtlev, uninp
      double precision w(*), x(*), xc, y(*), yc, wc, z(*), zc
cLocal:
      integer nmotsmax
      parameter(nmotsmax=255)
      integer i, lmot, nchi, nlmax, nl, nlu, nmots, stats(0:nmotsmax)
      double precision rlu
      character chiffres*15, lignlg*(nmotsmax),
     .        mots(nmotsmax)*(nmotsmax), program*9, progrer*12,
     .        progrwn*12
cBegin:
      program=' Readxyz>'
      progrer='%Readxyz-Er>'
      progrwn='%Readxyz-Wn>'
      chiffres='1234567890.eE+-'
      wc=0
      xc=0
      yc=0
      zc=0
      do i=1,nmotsmax
         stats(i)=0
      enddo
c     Lecture ligne a ligne:
c     ----------------------
      if (prtlev.gt.1) write(6,'(/2A)') program,
     .  ' Comments, or lignes with less than three numbers: '
      qerror=.false.
      ncoor=0
      nl=0
 100  continue
      read(uninp,'(A)',end=200) lignlg
      nl=nl+1
      call string_split(lignlg,nmotsmax," ",mots,nmotsmax,nmots)
      nfound=0
      do i=1,nmots
         call stringcl2(mots(i),lmot)
         if (lmot.le.0) goto 150
c        Commentaire ?
         if (mots(i)(1:1).eq.'!') goto 150
c        Chiffre ?
         do k=1,lmot
            if (index(chiffres,mots(i)(k:k)).le.0) goto 110
         enddo
         nfound=nfound+1
         if (nfound.le.4) then
             read(mots(i)(1:lmot),*,err=110) rlu
             if (nfound.eq.1) xc=rlu
             if (nfound.eq.2) yc=rlu
             if (nfound.eq.3) zc=rlu
             if (nfound.eq.4) wc=rlu
         else if (nfound.eq.5) then
             read(mots(i)(1:lmot),*,err=110) nlu
         endif
c        Mot suivant:
 110     continue
c        Le premier mot n'est pas un chiffre => ligne de commentaires
         if (nfound.eq.0) goto 150
      enddo
 150  continue
c     Stockage des coordonnees:
c     -------------------------
      stats(nfound)=stats(nfound)+1
      if (nfound.ge.3) then
          ncoor=ncoor+1
          if (ncoor.le.nmax) then
              x(ncoor)=xc
              y(ncoor)=yc
              z(ncoor)=zc
              if (nfound.eq.4) then
                  w(ncoor)=wc
              else
                  w(ncoor)=wc
                  ic(ncoor)=nlu
              endif
          else
              write(6,'(/2A,I9,A)') progrer,' More than ',
     .        nmax,' particles in file.'
              write(6,'(2A)') progrer,
     .      ' Please increase program memory limits (Sorry for that).'
              stop
          endif
      else
          if (prtlev.gt.1) then
              write(6,'(2A)') lignlg(1:72),'...'
          endif
      endif
c     Ligne suivante:
      goto 100
 200  continue
      if (prtlev.gt.1.and.nl.eq.ncoor) write(6,'(2A)') program,' None.'
      write(6,'(/2A,I7)') program,
     .' Number of particles in file (with x,y,z coordinates): ',ncoor
      if (ncoor.eq.0) then
          write(6,'(/2A)') progrer,' No coordinate found in file.'
          qerror=.true.
      endif
      nchi=0
      ndat=0
      nlmax=0
      do i=1,nmotsmax
         if (stats(i).gt.0) nchi=nchi+1
         if (stats(i).gt.nlmax) then
             nlmax=stats(i)
             ndat=i
         endif
      enddo
      do i=0,nmotsmax
         if (stats(i).gt.0.and.(prtlev.gt.1.or.nchi.gt.1)) then
            write(6,'(A,I6,A,I7,A)') program,i,
     .    ' numbers found in ',stats(i),' lignes.'
         endif
      enddo
      return
      end
"ENDOF readxyz.f"
#=======================================================================
echo string_split.f
cat >string_split.f <<"ENDOF string_split.f"
c
      subroutine string_split(chaine,taille,delimiteur,
     .                        souschaine,nbremax,nbre)
c
c     "Chaine" est coupee en "nbre" "souschaine" de part et d'autre du
c     "delimiteur"
c      YHS-Sep-93, Uppsala
c I/O:
      integer taille, nbremax, nbre
      character*(*) chaine, souschaine(*), delimiteur
c Local:
      integer icar, iprev
c
      nbre=1
      iprev=1
      souschaine(1)=chaine
      do icar=1,taille
         if (chaine(icar:icar).eq.delimiteur) then
            if (icar-1.ge.iprev) then
               souschaine(nbre)=chaine(iprev:icar-1)
               nbre=nbre+1
               if (nbre.le.nbremax) then
                  if (icar+1.le.taille.and.
     .               chaine(icar+1:taille).ne.' ') then
                     souschaine(nbre)=chaine(icar+1:taille)
                  else
                     nbre=nbre-1
                     return
                  endif
               else
                  write(6,'(A,I6,A/A)')
     .               ' %String_split-Err: more than ',nbremax,
     .               ' substrings in : ',chaine
                  return
               endif
            endif
            iprev=icar+1
         endif
      enddo
c
      return
      end
"ENDOF string_split.f"
#=======================================================================
echo stringcl.f
cat >stringcl.f <<"ENDOF stringcl.f"
c
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-95, Toulouse.
c     YHS-Oct-00, Bordeaux.
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      nonblancs=0
      taille=len(chaine)
      if (taille.le.0) return
c
      if (index(chaine(1:taille),' ').le.0) then
          nonblancs=taille
          return
      endif
c
c*****Nettoyage des blancs a gauche.
c     Premier non-blanc:
c
      do icar=1,taille
         if (chaine(icar:icar).ne.' ') goto 150
      enddo
      icar=taille
 150  continue
      chaine=chaine(icar:taille)
c
c*****Nettoyage des blancs au milieu.
c
          icar=1
          ncar=1
 170      continue
          icar=icar+1
          ncar=ncar+1
          if (chaine(icar:icar).eq.' ') then
              chaine=chaine(1:icar-1)//chaine(icar+1:taille)
              icar=icar-1
          endif
          if (ncar.lt.taille-1) goto 170
c
      nonblancs=index(chaine,' ')-1
c
      return
      end
"ENDOF stringcl.f"
#=======================================================================
echo stringcl2.f
cat >stringcl2.f <<"ENDOF stringcl2.f"
c-----------------------------------------------------------------------
      subroutine stringcl2(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires
c    (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-95, Toulouse.
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      taille=len(chaine)
c
      if (index(chaine(1:taille),' ').le.0) then
          nonblancs=taille
          return
      endif
c
c*****Nettoyage des blancs a gauche.
c     Premier non-blanc:
c
      do icar=1,taille
         if (chaine(icar:icar).ne.' ') goto 150
      enddo
 150  continue
      chaine=chaine(icar:taille)
c
c*****Nettoyage des blancs au milieu.
c
          icar=1
          ncar=1
 170      continue
          icar=icar+1
          ncar=ncar+1
          if (chaine(icar:icar).eq.' ') then
              chaine=chaine(1:icar-1)//chaine(icar+1:taille)
              icar=icar-1
          endif
          if (ncar.lt.taille-1) goto 170
c
      nonblancs=index(chaine,' ')-1
c
      return
      end
"ENDOF stringcl2.f"
#=======================================================================
echo tqli.f
cat >tqli.f <<"ENDOF tqli.f"
c-----------------------------------------------------------------------
      SUBROUTINE TQLI(D,E,N,NP,Z)
c     Finds the eigenvalues and eigenvectors of a tridiagonal matrix:
      integer i, iter, k, l, m, n, np
      double precision b, c, D(NP), dd, E(NP), f, g, p, r, s,
     .       Z(NP,NP)
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) STOP 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
"ENDOF tqli.f"
#=======================================================================
echo tred2.f
cat >tred2.f <<"ENDOF tred2.f"
c-----------------------------------------------------------------------
      SUBROUTINE TRED2(A,N,NP,D,E)
c     Reduce the matrix to tridiagonal form.
      integer i, j, k, l, n, np
      double precision A(NP,NP), D(NP), E(NP), f, g, h, hh,
     .  scale
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END
"ENDOF tred2.f"
#=======================================================================
echo trier.f
cat >trier.f <<"ENDOF trier.f"
c----------------------------------------------------------------------
      subroutine trier(y,npoint,nmax,ysort,qcrois)
c
c     Tri par ordre croissant ou decroissant.
c     YHS-Jun-2002: premiere version (Bordeaux).
      implicit none
      integer i, j, nmax, npoint
      logical qcrois
      double precision y(*), ycur, ysort(*)
      character progrer*10
      progrer='%Trier-Er>'
      if (npoint.gt.nmax) then
          write(6,'(A,I9,A,I9,A)') progrer,npoint,
     .  ' points to be sorted, i.e., more than ',nmax,' Sorry.'
          stop
      endif
      do i=1,npoint
         ysort(i)=y(i)
      enddo
      if (qcrois) then
      do i=1,npoint
         do j=1,npoint
            if (ysort(i).lt.ysort(j)) then
                ycur=ysort(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
            endif
         enddo
      enddo
      else
      do i=1,npoint
         do j=1,npoint
            if (ysort(i).gt.ysort(j)) then
                ycur=ysort(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
            endif
         enddo
      enddo
      endif
      return
      end
"ENDOF trier.f"
#=======================================================================
echo trieri.f
cat >trieri.f <<"ENDOF trieri.f"
c-----------------------------------------------------------------------
      subroutine trieri(y,npoint,nmax,ysort,iord,qcrois)
c
c     Tri par ordre croissant (qcrois=T) ou non.
c     Version triviale...
c     YHS-Jun-2002: Premiere version (Bordeaux).
c     YHS-Sep-2002: Derniere version (Bordeaux).
      implicit none
      logical qcrois
      integer i, icur, iord(*), j, nmax, npoint
      double precision y(*), ycur, ysort(*)
      character progrer*10
      progrer='%Trier-Er>'
      if (npoint.gt.nmax) then
          write(6,'(A,I9,A,I9,A)') progrer,npoint,
     .  ' points to be sorted, i.e., more than ',nmax,' Sorry.'
          stop
      endif
      do i=1,npoint
         ysort(i)=y(i)
         iord(i)=i
      enddo
      do i=1,npoint
        do j=1,npoint
          if (qcrois) then
            if (ysort(i).lt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          else
            if (ysort(i).gt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          endif
        enddo
      enddo
      return
      end
"ENDOF trieri.f"
#=======================================================================
echo vecstat.f
cat >vecstat.f <<"ENDOF vecstat.f"
c-----------------------------------------------------------------------
      subroutine vecstat(vect,nmax,rmin,rmax,rave,rdev)
c     Statistics for vector Vect(NMAX):
c     minimum, maximum, average (rave) and standard deviation (rdev).
c     YHS-Sep-03: First version (Lyon).
cI/O:
      integer nmax
      double precision rave, rdev, rmax, rmin, vect(*)
cLocal:
      integer i
      character program*9, progrer*12, progrwn*12
cBegin:
      program=' Vecstat>'
      progrer='%Vecstat-Er>'
      progrwn='%Vecstat-Wn>'
      rave=0.d0
      rdev=0.d0
      rmin=-9999.d0
      rmax=9999.d0
      if (nmax.le.0) then
          write(6,'(2A)') progrer,' Zero-length vector.'
          return
      endif
      do i=1,nmax
         if (vect(i).gt.rmax.or.i.eq.1) rmax=vect(i)
         if (vect(i).lt.rmin.or.i.eq.1) rmin=vect(i)
         rave=rave+vect(i)
         rdev=rdev+vect(i)**2.0
      enddo
      rave=rave/dfloat(nmax)
      rdev=rdev/dfloat(nmax)-rave*rave
      if (rdev.gt.0.d0) rdev=dsqrt(rdev)
      return
      end
"ENDOF vecstat.f"
#=======================================================================
echo writpdb.f
cat >writpdb.f <<"ENDOF writpdb.f"
c-----------------------------------------------------------------------
      subroutine writpdb(unpdb,xat,yat,zat,binfo,
     .           atonam,ires,resnam,segid,natom)
c
c     Ecriture d'un fichier en format pdb.
c     YHS-Mars-1996.
c
      implicit none
cI/O:
      integer unpdb, ires(*), natom
      double precision xat(*), yat(*), zat(*), binfo(*)
      character*4 atonam(*), resnam(*), segid(*)
cLocal:
      integer i
      character*5 atom
cBegin:
      write(6,'(/A,I5)') ' Writpdb> Writing pdb file.',unpdb
c
      do i=1,natom
c
c     Un petit probleme avec le cadrage des noms des atomes:
c    'iHjj' ou ' NNN'
c
      atom=atonam(i)
      if (index(atonam(i),' ').gt.1) atom=' '//atonam(i)
c
c     write(unpdb,'(A,I7,1X,A4,A4,2X,I4,4X,3F8.3,6X,F6.2,6X,A)')
      write(unpdb,'(A,I7,1X,A4,1X,A4,1X,I4,4X,3F8.3,6X,F6.2,6X,A)')
     .      'ATOM', i, atom(1:4), resnam(i), ires(i),
     .      xat(i), yat(i), zat(i), binfo(i), segid(i)
      enddo
      write(unpdb,'(A)') 'ENDMDL'
c
      write(6,'(A,I6,A)') ' Writpdb> ',natom,' atoms saved.'
      return
      end
"ENDOF writpdb.f"
