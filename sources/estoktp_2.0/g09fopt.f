c g09fopt a program to start and read output from g09
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.

      SUBROUTINE g09fopt(ilevcode,tau,ntau,natom,natomt,nshared,gmemo,
     $ coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $ comline2,icharge, ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ired)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c  lc: interface to case insensitive function
      interface
            function Up_to_low(vlstring)
               implicit none
            character(len=*),intent(in) :: vlstring 
            character(len=len_trim(vlstring)) :: Up_to_low
            integer ::  i             ! COUNTER 
            intrinsic :: achar
            end function Up_to_low
      end interface
c  lc END
      

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c     parameter (ntaumx=10, nmdmx=300, natommx=100, ndim=3)

      integer natom,iatom,idim
      dimension freq(nmdmx),tau(ntaumx),tauopt(ntaumx)
      dimension coord(natommx,ndim),xint(3*natommx),abcrot(ndim)
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension idummy(natommx)
      dimension grad(3*natommx)
      character*30 gkeyword,igkey
      character*(*) comline1,comline2
      character*30 gmemo
      character*60 atomlabel(natommx)
      character*60 irconslabel
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*2 aname
      character*4 cgroup
      character*120 command1
      character*30 cjunk,cjunk1,iqclab,cjunk2,cjunk3
      character*250 longcommand1
      character*512 ExpGauexedir
 
c     character*160 cname

      LOGICAL leof,lsec,ltit
c  lc mod , single point or optimization flag
c  lc mod, BasPath is a string containing the path to
c          the custom basis set 
      integer lc_flag,BasFlag,Pt2Flag,AnhFlag,AnhCounter
      character*300 BasPath
c  lc mod,  Pt2Model contains information about the pt2model
c           to use 
      character*50  Pt2Model, lc_word 
      character*1024 vlstring
      character*1024 hindfilecheck

c  lc mod,  COUNTERPOISE CORRECTION, INSERTION OF TWO FRAGMENTS
C           VARIABLES 
      integer icharge_frag1,ispin_frag1,icharge_frag2,ispin_frag2
      integer BSSE_flag
      integer EndPoi
      integer Frag_index
      logical OptFreq
      logical HindCheck
c
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7
      character*20 ct1,ct2,ct3
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)
c  lc mod: hanlde EigenvalueFollow gaussian algorithm
      logical EFalgo
      integer NormTermCounter
c  lc mod: Logical flag to be activated if skip of computation of geometric potential
c          must be done for infinite computation
      logical SkInfGeomCorr
c      COMMON /keyword/ leof,lsec,ltit
c      COMMON /key/ line,sename,string,word,word2,word3
c     $ ,word4,word5,title,title1

c      common /tsdat/ aabs1,babs1,aabs2,babs2,babs3,rts
c      common /itsdat/ iabs,iadd,ivar,isite,ji,ki
c      common /strucword/ stoich
      include 'filcomm.f'

c111   continue
c  lc mod, setting lc_flag locally to -1:  
      lc_flag=-1 
c  lc mod, setting Pt2Flag locally to -1
c          if Pt2flag is not -1 this will activate
c           the Pt2model search 
      Pt2Flag=-1
      AnhFlag=-1
      AnhCounter=1   
      BSSE_flag=-1
      FrGueF=-1
      NormTermCounter=0
      EFalgo=.FALSE.
      OptFreq=.FALSE.
      HindCheck=.FALSE.
      SkInfGeomCorr=.FALSE.
c  lc mod,  setting BasFlag to -1 at the start of the procedure
c           if this will be set to a different value is because
c           the gen keyword is read in the actual command line 
c           for the actual computation. A secondary file will be
c           searched. This file will contain the path to the custom 
c           basis set and it will be inserted in the .com file.
      BasFlag=-1 
      BasPath='Path_not_setted'
      Pt2Model='Pt2Model not set'


c  lc Check hindered file procedure
      open(unit=127,file='HindNormalCheck.log',iostat=iEoFfl,
     & status='old' )
      if (iEoFfl.eq.0) then 
         read(127,'(A)',iostat=IfileCheck) hindfilecheck
         hindfilecheck=Up_to_low(hindfilecheck)
         if (index(hindfilecheck,'normal') .gt. 0) then
            read(127,'(A)') hindfilecheck
            longcommand1='cp '//trim(hindfilecheck)//' geom.log'
            call commrun(longcommand1)
            hindfilecheck=
     &      trim(hindfilecheck(1:len_trim(hindfilecheck)-3))//'chk'
            longcommand1='cp '//trim(hindfilecheck)//' tmp.chk'
            call commrun(longcommand1)
            longcommand1=' '
            close(127,status='delete')
            HindCheck=.TRUE.
            longcommand1='echo HINDERED ROTORS READING ACTIVATED'
            call commrun(longcommand1)
         end if 
      endif

c  lc Check for skip geom corr file and chk. If one does not find 
c  lc potocorrgeom and also is chk file it assumes that a normal
c  lc procedure must be performed

      open(unit=127,file='skipgeominf.dat',iostat=iEoFfl,
     & status='old')
      if (iEoFfl.eq.0) then
         read(127,*) iSkipGeo
         if (iSkipGeo.eq.6) then 
            open(unit=128,file='geoms/potcorr_geom.log',iostat=iEoFfl,
     & status='old')
            open(unit=129,file='geoms/potcorr_geom.chk',iostat=iEoFfl2,
     & status='old')
            if ( 
     &         (iEoFFl.eq.0)
     &            .and. 
     &         (iEoFFl2.eq.0)
     &          ) then
                  if (
     &                index(Up_to_low(comline1),'opt') .gt.0 
     &                            .or. 
     &                index(Up_to_low(comline2),'opt') .gt.0 
     &                ) then 
                     SkInfGeomCorr=.TRUE.
                     longcommand1='cp geoms/potcorr_geom.log geom.log'
                     call commrun(longcommand1)
                     longcommand1='cp geoms/potcorr_geom.chk tmp.chk'
                     call commrun(longcommand1)
                     longcommand1=' '
                     write(6,*) 'optimization skipped for inf'
                     write(6,*) comline1
                     write(6,*) comline2
                     write(6,*) '-----------------------'
                  end if
            end if
            close(128,status='keep')
            close(129,status='keep')
         end if 
         close(127,status='delete')
      end if  

      command1="rm -f temp_route.dat"
      call commrun(command1)

cccccccccccccc print out gaussian com file  ccccccccccccc
c initialize word, word2, word3, word4, word5
      call LineRead (0)
      ichecken=0
      irepeat=0
      iqc=0
      icsy=0
      ilin_fr=0
      ifilu=0
      ifreq_save=ifreq
      if(natom.eq.2)ilin=0
      if(ilin.eq.1) ilin_fr=1
 999  continue
      ifreq=ifreq_save
cc determine igkey and gkeyword
      call en_key_g09(comline1,igkey,gkeyword,ispin)
      write(*,*)'gkeyword is ',gkeyword
      write(*,*)'igkey is ',igkey
      write(*,*)'ired is ',ired

c      stop

cc if using redundant coordinates then read z-matrix

c disable redundant coord if using dummy atoms in reactants

      if(natomt.ne.natom)then
         if(ibeta.eq.1.or.iiso.eq.1)then
            ired=0
         else if (iadd.eq.1)then
            if((natomt-1).ne.natom)then
c               ired=0
            endif
         endif
      endif
c
cc determine of the optimization is to be done in internal coords
      ired_check1=0
      ired_check2=0
      iredadd=0
      call int_key_g09(comline1,ired_check1)
      call int_key_g09(comline2,ired_check2)

      if( ired_check1.eq.1.or.ired_check2.eq.1)then
         ired=1
         iredadd=1
      endif

c      write(*,*)'ok up to here 1'
c      write(*,*)'ired is',ired
c      write(*,*)'ired_check1 is',ired_check1
c      write(*,*)'ired_check2 is',ired_check2
c
      if(ired.eq.1)then
c         write(*,*)'ok up to here 2' 
        call read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)
         open (unit=99,file='command.dat',status='unknown')
         write (99,'(A300)') comline1
         write (99,'(A300)') comline2
         close(99)
         open (unit=108,status='unknown')
         write (108,1200)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,1201)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,1202)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,1203)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         

 1200 format (" sed -ie 's/Internal,/modredundant,/g' command.dat")      
 1201 format (" sed -ie 's/internal,/modredundant,/g' command.dat")      
 1202 format (" sed -ie 's/Internal/modredundant/g' command.dat")      
 1203 format (" sed -ie 's/internal/modredundant/g' command.dat")      

         open (unit=99,file='command.dat',status='unknown')
         read (99,'(A300)') comline1
         read (99,'(A300)') comline2
         close(99)
      else if (ired.eq.2)then
cc         write(*,*)'pass from here'
         call read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)
         open (unit=99,file='command.dat',status='unknown')
         write (99,'(A300)') comline1
         write (99,'(A300)') comline2
         close(99)
         open (unit=108,status='unknown')
         write (108,2200)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,2201)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,2202)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,2203)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         
         open (unit=108,status='unknown')
         write (108,2204)
         rewind (108)
         read (108,'(A300)') command1
         close (108)
         call commrun(command1)         

 2200 format (" sed -ie 's/Internal,/cartesian,/g' command.dat")      
 2201 format (" sed -ie 's/internal,/cartesian,/g' command.dat")      
 2202 format (" sed -ie 's/Internal/cartesian/g' command.dat")      
 2203 format (" sed -ie 's/internal/cartesian/g' command.dat")      
 2204 format (" sed -ie 's/modredundant/cartesian/g' command.dat")      

         open (unit=99,file='command.dat',status='unknown')
         read (99,'(A300)') comline1
         read (99,'(A300)') comline2
         close(99)
         ired=1
      endif

      write(*,*)'ok up to here 2'
      command1='echo $PWD'
      call commrun(command1)         

      OPEN(UNIT=10,STATUS='unknown',FILE='geom.com')
      REWIND (10)
c  lc mod: open file to write route section 
      open(999,file='temp_route.dat',status='unknown')
      write(999,*) comline1,comline2 
      close(999)
c  lc : Mods, introduce a control section to handle opt+freq gaussian bug
c       the bug consists in the fact that in an opt+freq computation the route section
c       is taken without the iop. 
c       check if the frequency computation is required
c   1) Obtain one route section to scan
      write(vlstring,'(A,2x,A)') trim(adjustl(comline1)),
     & trim(adjustl(comline2))
c   2) Convert in lower case, case insensitive research
      vlstring=Up_to_low(vlstring)
c   3) Check if the route contains freq and opt togheter
      if (
     &    (index(vlstring,'freq') .gt. 0) 
     &     .and.
     &    (index(vlstring,'opt') .gt. 0) 
     &    .and.
     &    (index(vlstring,'maxcyc=1') .eq. 0))then
         OptFreq=.True.
         write(10,'(A)') '--link1--'
      end if
      if (
     &   index(vlstring,',ef,') .gt. 0
     &      .or.
     &   index(vlstring,'eigenvaluefollow') .gt. 0   
     &   ) then 
       EFalgo=.true. 
       write(*,*) 'Eigenvalue algorithm checked'
      end if

c  lc : Input gaussian 
      write (10,*) '%mem=',gmemo
      write (10,*) '%chk=tmp'
      write (10,*) '%NProcShared=',nshared
cc      write (10,*) '#',comline1
cc      write (10,*) '#',comline2
cc modified to cut length
      write (10,*) '#',comline1(1:len_trim(comline1))
      write (10,*) '#',comline2(1:len_trim(comline2))

c  lc mod:  search in command.dat file the presence of gen keyword
c           to handle custom basis set 
      command1="grep -icw 'gen' temp_route.dat > LC_custom_basis.dat"
      call commrun(command1)
c lc  mod: from this file, read the basis flag 
      open(unit=903,file='LC_custom_basis.dat',status='old')
      read(903,*) BasFlag 
      close(903)
      command1="rm -f LC_custom_basis.dat"
      call commrun(command1)

c lc  if the basflag is set to a value different from 0, then 
c     the path to the basis set file will be stored 
      if (BasFlag .ne. 0) then 
         open(unit=977,file='data/Custom_Basis_Path.dat',action='read'
     $    ,status='old')
         read(977,'(A)') BasPath
         close(977)
      end if 
c  lc check for anharmonic computation
      command1="grep -ic -P '(?<!read)anh' temp_route.dat > 
     & LC_ahn_flag.dat"
      call commrun(command1)
      open(unit=977,file='LC_ahn_flag.dat') 
      read(977,*) AnhFlag
      close(977)
      command1="rm -f LC_ahn_flag.dat"
      call commrun(command1)

c  lc check for ptmodel to use 
      command1="grep -ic 'readanh' temp_route.dat > 
     & LC_pt2model_flag.dat"
      call commrun(command1)
      open(977,file='LC_pt2model_flag.dat')
      read(977,*) Pt2Flag
      write(977,*) 'Pt2flag val=',Pt2Flag
      close(977) 
      if (Pt2Flag .ne. 0) then 
         command1="grep 'pt2model' data/theory.dat | awk '{print $2}' 
     &> LC_pt2Model_name.dat"
      call commrun(command1)
      open(977,file='LC_pt2Model_name.dat')
      read(977,*) Pt2Model 
      close(977) 
      end if 
      command1="rm -f LC_pt2Model_name.dat"
      call commrun(command1)
      command1="rm -f LC_pt2model_flag.dat"
      call commrun(command1)

c  lc :  CHECK FOR OPTIMIZATION FLAG START
      command1="grep -ic 'opt' temp_route.dat > LC_Opt_index.dat"
      open(unit=902,file='LC_Opt_index.dat',status='unknown')
      call commrun(command1)
      rewind(902)
      read(902,*) lc_flag  
      close(902,status='keep') 
      command1="rm -f LC_Opt_index.dat"
      call commrun(command1)
c  lc :  CHECK FOR OPTIMIZATION FLAG END

      if(iqc.eq.1.and.ired.eq.0)then  
         write (10,*) '# scf=qc'
      else if (iqc.eq.1.and.ired.eq.1)then
         write (10,*) '# scf=qc iop(5/6=7)'
cc         write (10,*) '# scf=qc iop(4/5=1)'
cc         write (10,*) '# scf=qc '
         ires=ires+1
      endif
c     write (10,*) '# nosym'
c     if (ifreq.eq.1) then 
c        write (10,*) '# freq'
c     endif
      write (10,*) 
c  lc :  title card
      write (10,*) 'geom ', ismp
      write (10,*) 
c      WRITE (10,*) icharge,ispin
c  lc :  case insensitive conversion of counterpoise correction
c        in geom.com file
      command1="grep -ic 'counterpoise' temp_route.dat > BSSE_flag.dat"
      call commrun(command1)
         open(977,file='BSSE_flag.dat') 
            read(977,*) BSSE_flag
         close(977)
         command1="rm -f BSSE_flag.dat"
         call commrun(command1)
      command1="grep -ic 'fragment' temp_route.dat > FragmentGuess.dat"
      call commrun(command1)
         open(977,file='FragmentGuess.dat') 
            read(977,*) FrGueF
         close(977)
         command1="rm -f FragmentGuess.dat"
         call commrun(command1)

c  lc :  search for substring counterpoise in comline1
      if ((BSSE_flag .gt. 0) .or. (FrGueF .gt. 0)) then
         open(977,file='data/counter_fragments.dat')
         read(977,*) icharge_frag1,ispin_frag1,
     &               icharge_frag2,ispin_frag2
         WRITE (10,*) icharge,ispin,
     &                icharge_frag1,ispin_frag1,
     &                icharge_frag2,ispin_frag2
      else 
         WRITE (10,*) icharge,ispin
      end if 
c lc  todo COUNTERPOISE : COMMENT UNTIL CLOSE(977) TO RETURN TO STANDARD VERSION
C     IF COUNTERPOISE FLAG IS ACTIVATED
C     CHECK IF THE ATOM LABEL IS NOT PRESENT
C     THE SUBSTRING "FRAGMENT" 
C     IF IT IS NOT PRESENT
C     OPEN FILE COUNTERPOISE .DAT
C     WRITE(ATOMLABEL(IATOM(1:3)),*) TRIM(ATOMLABEL(IATOM(1:3)))//"FRAGMENT",FRAGMENT INDEX_nf//"ATOMLABEL(IATOM(3:))
C     WRITE(10,*) AS IT IS.      
C     F77 Index_nf function for check of substring in a string
C     IS OVERWRITTEN BY INDEX_nf VARIABLE... USE A CHECK FOR ) 
C     AT THE END OF THE STRING (I WILL USE F90 INT FUNC)
C     Defining a pointer to the end of the atomlabel 
C     In this case Fragment= is assumed to be present
C     The writing will be as the standard version 
C
      if ((BSSE_flag .gt. 0) .or. (FrGueF .gt. 0)) then 
         if(ixyz.eq.0)then
            DO iatom = 1 , natomt
               read(atomlabel(iatom),*) lc_word 
               EndPoi=len_trim(lc_word)
               if ( lc_word(EndPoi:EndPoi) == ')' ) then 
                  write (10,*) atomlabel(iatom)
               else
                  atomlabel(iatom)=adjustl(trim(atomlabel(iatom)))
                  lc_word=adjustl(trim(lc_word))
                  read(977,*) Frag_index
                  write(10,'(A,A,I0,A,3x,A)') 
     &             trim(lc_word),
     &             '(Fragment=',Frag_index,')',
     &             atomlabel(iatom)(1+len_trim(lc_word):)
               end if
            ENDDO
         else
            DO iatom = 1 , natom
               read(atomlabel(iatom),*) lc_word 
               EndPoi=len_trim(lc_word)
               if ( lc_word(EndPoi:EndPoi) == ')' ) then 
                  write (10,*) atomlabel(iatom)
               else
                  atomlabel(iatom)=adjustl(trim(atomlabel(iatom)))
                  lc_word=adjustl(trim(lc_word))
                  read(977,*) Frag_index
                  write(10,'(A,A,I0,A,3x,A)') 
     &             trim(lc_word),
     &             '(Fragment=',Frag_index,')',
     &             atomlabel(iatom)(1+len_trim(lc_word):)
               end if
            ENDDO
         endif
      else 
         if(ixyz.eq.0)then
            DO iatom = 1 , natomt
               write (10,*) atomlabel(iatom)
            ENDDO
         else
            DO iatom = 1 , natom
               write (10,*) atomlabel(iatom)
            ENDDO
         endif
      end if 
 
      close(977)

C UNCOMMENT TO RETURN TO STANDARD VERSION
c      if(ixyz.eq.0)then
c         DO iatom = 1 , natomt
c            write (10,*) atomlabel(iatom)
c         ENDDO
c      else
c         DO iatom = 1 , natom
c            write (10,*) atomlabel(iatom)
c         ENDDO
c      endif
c  lc:   since the if algorithm does not print out the freezed coordinates
c        a file will be used to storage these         
      if (EFalgo) then 
         open(unit=157,file='freezecoord.dat',status='unknown')
      end if 
      write (10,*)
      if(ixyz.eq.0)then
         ncoord = natom*3-6-ntau-ircons
         if (natom.eq.2) ncoord = 1
         do icoord = 1 , ncoord
            write (10,*) intcoor(icoord),xint(icoord)
         enddo
         do itau = 1 , ntau
            write (10,*) bislab(itau), tau(itau)
         enddo
         icoord = natom*3-6
         if (ircons.eq.1) then
            write (10,*)
            write (10,*) intcoor(icoord),xint(icoord)
         else if (ircons.eq.2) then
            write (10,*)
            write (10,*) intcoor(icoord-1),xint(icoord-1)
            write (10,*) intcoor(icoord),xint(icoord)
         else if (ircons.eq.3) then
            write (10,*)
            write (10,*) intcoor(icoord-2),xint(icoord-2)
            write (10,*) intcoor(icoord-1),xint(icoord-1)
            write (10,*) intcoor(icoord),xint(icoord)
         else if (ircons.eq.4) then
            write (10,*)
            write (10,*) intcoor(icoord-3),xint(icoord-3)
            write (10,*) intcoor(icoord-2),xint(icoord-2)
            write (10,*) intcoor(icoord-1),xint(icoord-1)
            write (10,*) intcoor(icoord),xint(icoord)
         else if (ircons.eq.5) then
            write (10,*)
            write (10,*) intcoor(icoord-4),xint(icoord-4)
            write (10,*) intcoor(icoord-3),xint(icoord-3)
            write (10,*) intcoor(icoord-2),xint(icoord-2)
            write (10,*) intcoor(icoord-1),xint(icoord-1)
            write (10,*) intcoor(icoord),xint(icoord)
         else if (ircons.gt.5) then
            write (10,*)
            do j=1,ircons
               index_nf=ircons-j
               write (10,*) intcoor(icoord-index_nf),
     &          xint(icoord-index_nf)
            enddo
         endif
c         write (10,*)
      endif
c
      icheckdummy=0
      if(natomt.ne.natom)then
c         if(ibeta.eq.1.or.iiso.eq.1.or.iadd.eq.1)then
         if(ibeta.eq.1.or.iiso.eq.1)then
            icheckdummy=1
c         else if (iadd.eq.1)then
c            if((natomt-1).ne.natom)then
c               icheckdummy=1
c            endif
         endif
      endif

      if(ired.eq.1.and.icheckdummy.eq.0.and.iredadd.ne.1)then
         if(ibeta.eq.1.or.iiso.eq.1)then
            if(ireact.ne.0)then
               if(isite.ne.1)then
                  write(10,*)
                  write(10,*)isite-1,ireact
               else if(ireact.ne.isite+1)then
                  write(10,*)
                  write(10,*)isite+1,ireact
               endif
            endif
            write(10,*)
         endif
c         if(iadd.eq.1)then
c            natom2=natom-ireact+1
c            write(10,*)ireact,isited
c            write(10,*)ireact,jsited
c            write(10,*)ireact,ksited
c            if(ireact.ne.isite+1)then
c               write(10,*)isite-1,ireact
c               if(isite.ne.1) then
c                  if(ireact.lt.natom)then
c                     write(10,*) isite-1,ireact+1
c                     if(isited.ne.1.and.natom2.ne.1) then 
c                        write(10,*) ireact+1,isited-1
c                     endif
c                     if(isited.ne.1.and.natom2.gt.2) then 
c                        write(10,*)ireact+2,isited-1
c                     endif
c                  else
c                     write(10,*) isite-1,ireact
c                     write(10,*)ireact,isited
c                     write(10,*)ireact,jsited
c                     write(10,*)ireact,ksited
c                     if(isited.ne.1.and.natom2.ne.1) then 
c                        write(10,*) ireact+1,isited-1
c                     endif
c                     if(isited.ne.1.and.natom2.gt.2) then 
c                        write(10,*)ireact+2,isited-1
c                     endif
c                  endif
c               else if(ireact.ne.isite+1)then
c                  write(10,*)isite+1,ireact
c
c               endif
c            else if(ireact.eq.isite+1)then
c               write(10,*)isite-1,ireact
c               write(10,*)isite-2,ireact
c            endif
c         endif
         if(iabs.eq.1.or.iadd.eq.1)then
            natom2=natom-ireact+1
            write(10,*)
            write(10,*)ireact,isited
            write(10,*)ireact,jsited
            write(10,*)ireact,ksited
            if(isited.ne.1.and.natom2.ne.1) then 
               write(10,*)
               write(10,*) ireact+1,isited-1
            endif
            if(isited.ne.1.and.natom2.gt.2) then 
               write(10,*)
               write(10,*)ireact+2,isited-1
            endif
            if(natom2.gt.2) then 
               write(10,*)
               write(10,*)ireact+2,jsited
            endif
         endif
      endif
      
      if (BasFlag .ne. 0) then
         if (natom.gt.1) write(10,*) 
         write(10,'(A)') '@'//BasPath
      end if 
      write(10,*)
      if (OptFreq .eqv. .false.) then 
         if (Pt2Flag .ne. 0) then 
            write(10,*) 'pt2model='//Pt2Model 
         endif 
      end if
      
c
c

c  lc: mod link1 
      if (OptFreq .eqv. .true.) then
         write (10,'(A)') '--link1--'
         write (10,*) '%mem=',gmemo
         write (10,*) '%chk=tmp'
         write (10,*) '%NProcShared=',nshared
cc         write (10,*) '#',comline1
cc         write (10,*) '#',comline2
cc modified to cut length
         write (10,*) '#',comline1(1:len_trim(comline1))
         write (10,*) '#',comline2(1:len_trim(comline2))

c         write (10,*) '#','guess=read geom=allcheck ChkBasis'
         write (10,*) '#','guess=read geom=allcheck '
         write (10,*)
         write (10,*)
         write (10,*) 
         if (Pt2Flag .ne. 0) write(10,*) 'pt2model='//Pt2Model  
         write (10,*)
         write (10,*)
         write (10,*) 
         command1=
     &    "sed -i -E ':a;N;$!ba;s/opt=?\(([^()]*)\)/ /2' geom.com"
         call commrun(command1)
         command1=
     &    "sed -i -E ':a;N;$!ba;s/freq(=?\(([^()]*)\))?/ /1' geom.com"
         call commrun(command1)
         command1=
     &    "sed -i -E ':a;N;$!ba;s/,?always/ /2I' geom.com"
         call commrun(command1)
         if (BasFlag .ne. 0) then
         command1=
     &    "sed -i -E ':a;N;$!ba;s/\bgen\b/ /2I' geom.com"
         call commrun(command1)
         end if
         
      end if

      close (unit=10,status='keep')

      if(ilin.eq.1.and.ired.eq.0) then
         write(command1,919)
         call commrun(command1)         
      endif
 919  format(" sed -ie 's/nosym/symm(loose,follow)/g' geom.com")

c            write (*,*)'intcoord is ',intcoor(1)
c      stop

ccccccccccccc run g09  cccccccccccccccc

c      stop
c  lc :  new feature
c        grep the path to gaussian directory
c        from data/GPaths.dat and then export
c        to env variable GAUSS_EXEDIR
      IF ((HindCheck .eqv. .FALSE.)
     &             .and.
     &    (SkInfGeomCorr.eqv..FALSE.)) then 
      if(irecov.ne.1)then
         if(ilevcode.eq.1)then
            ExpGauexedir=' '
            longcommand1='export GAUSS_EXEDIR='
            open(991,file='data/GPaths.dat',status='old',action='read',
     &      iostat=iOpstat)
c  lc First check if the file with gaussian paths exists. 
            if (iOpstat.eq.0) then 
c  lc Search for g09 keyword
               do 
                  read(991,*,iostat=istatus) cjunk
                  cjunk=Up_to_low(cjunk)
                  if (index(cjunk,'g09').gt.0) then 
                     rewind(991)
                     read(991,'(A)') ExpGauexedir
                     ExpGauexedir=
     &               ExpGauexedir(index(ExpGauexedir,'path')+5:)
                     exit 
                  end if
               end do 
               close(991,status='keep')
               longcommand1=trim(longcommand1)//trim(ExpGauexedir)
c  lc problem to fix:
c                    Up to now if I use the croutine to set
c                    The environment variable to set the new 
c                    Gaussian directory path the computation starts
c                    But does not write anything to the .log file
c                    And I cannot change the environment variable from
c                    A subprocess by fortran command or executing
c                    A shell script 
c-----------------------------------------------------------------------
c               open(991,file='setenv.sh')
c               write(991,'(A)') '#!/bin/bash'
c               write(991,'(A)') longcommand1
c               close(991)
c               call setgvar(ExpGauexedir)
c               longcommand1='chmod +x setenv.sh'
c               call commrun(longcommand1)
c               call execute_command_line(longcommand1)
c               longcommand1='. ./setenv.sh'
c               call commrun(longcommand1)
c               call commrun(longcommand1)
c              longcommand1=' '
c              ExpGauexedir=' '
            else 
               write(*,*) 'WARNING, NO GPATHS FILE DETECTED IN DATA' 
            end if
            call g09run()
         else if (ilevcode.eq.3) then
            ExpGauexedir=' '
            longcommand1='export GAUSS_EXEDIR='
            open(991,file='data/GPaths.dat',status='old',action='read',
     &      iostat=iOpstat)
c  lc First check if the file with gaussian paths exists. 
            if (iOpstat.eq.0) then 
c  lc Search for g09 keyword
               do 
                  read(991,*,iostat=istatus) cjunk
                  cjunk=Up_to_low(cjunk)
                  if (index(cjunk,'g16').gt.0) then 
                     rewind(991)
                     read(991,'(A)') ExpGauexedir
                     ExpGauexedir=
     &               ExpGauexedir(index(ExpGauexedir,'path')+5:)
                     exit 
                  end if
               end do 
               close(991,status='keep')
               longcommand1=trim(longcommand1)//trim(ExpGauexedir)
               open(991,file='setenv.sh')
               write(991,'(A)') '#!/bin/bash'
               write(991,'(A)') longcommand1
               close(991)
               longcommand1='chmod +x setenv.sh'
               call commrun(longcommand1)
               longcommand1='./setenv.sh'
               call commrun(longcommand1)
               longcommand1='rm -f setenv.sh'
               call commrun(longcommand1)
              longcommand1=' '
              ExpGauexedir=' '
            else 
               write(*,*) 'WARNING, NO GPATHS FILE DETECTED IN DATA' 
            end if
            call g16run()
         endif
      endif
      ENDIF 

ccccccccccccccccccccccccccccccccccccccc



cc first check if the framework group is linear
      open(unit=920,file='temp0.dat',status='unknown')
      write(920,*)' Framework group  C1 '
      close(920)
      command1='egrep Frame geom.log > temp1.dat'
      call commrun(command1) 
      command1='cat temp0.dat temp1.dat > temp.dat'
      call commrun(command1)         
      command1='tail -n 1 temp.dat > temp1.dat'
      call commrun(command1)         
      write(command1,921)
      call commrun(command1)         

      open(unit=920,file='temp1.dat',status='unknown')
      read(920,*)cjunk,cjunk,cgroup
      if((cgroup.ne.'C*V').and.(cgroup.ne.'D*H')) then 
         ilin_fr=0
         ilin=0
      endif
      close(920)
      write(*,*)'cgroup ',cgroup
      write(*,*)'ilin_fr ',ilin_fr
c      stop
 921  format(" sed -ie 's/\[/ /g' temp1.dat")

cc now check if frequencies were calculated        
      if (index(vlstring,'freq') .eq. 0) ifreq=0
c      command1='egrep Freq geom.log > temp.dat'
c      call commrun(command1) 
c      command1='wc temp.dat > temp1.dat'
c      call commrun(command1)
c      nlines=0
c      open(unit=920,file='temp1.dat',status='unknown')
c      read(920,*)nlines
c      close(920)
c      if(nlines.eq.0) then 
c         ifreq=0
c         write(*,*) 'FREQUENCIES WILL NOT BE PROCESSED'
c      end if 

ccc now read output file
c

c  lc:   Debug for EF optimization algorithm. 
c        First read the geom.log file to search 
c        Constants internal variables and put in
c        a separate file 

      if (EFalgo) then
         write(*,*) 'eigenvalue are you checked?'
         open(unit=11,status='old',file='geom.log')
         do 
            read(11,'(A)') vlstring
            vlstring=Up_to_low(vlstring)
            if (index(vlstring,'constants:').gt.0) then
               do
                  read(11,'(A)') vlstring
                  if (len(trim(vlstring)).eq. 0) then 
                     rewind(11)
                     rewind(157)
                     close(11)
                     go to 777
                  else 
                     read(vlstring,*) cjunk2,cjunk3
                     write(157,*) cjunk3,cjunk2  
                  endif
               end do
            endif
         enddo 
      end if 

777   continue
      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
c      open(unit=107,status='unknown',file='LC_output.out')

c start with reference potential
c  lc debug: here read energy from fchk file and not from 
c           the .log because the igkey point to double hybrid
c           functionals
      ie = 0
      if ((igkey.eq.'notfound').or.(igkey.eq.'4'))then
         write(*,*)'the gkeyword was not found'
         write(*,*)'the energy will be determined from the chk'
         write(*,*)'using the formchk.exe code'
         write(*,*)'to speed up data recovery insert the gkeyword'
         write(*,*)'for this method in the en_key_09 subroutine'
         if(ilevcode.eq.1)then
            command1='formchk09 tmp.chk'
         else if (ilevcode.eq.3)then
            command1='formchk09 tmp.chk'
         endif
         call commrun(command1)       

         open(unit=778,file='tmp.fchk',iostat=iStatus,status='old')
            if (iStatus .ne. 0 ) then 
               command1='formchk09 tmp.chk'
               call commrun(command1)  
            end if 
            close(778,status='keep')
         command1="egrep 'Total Energy' tmp.fchk > temp.log"
         call commrun(command1)         
         open(unit=65,file='temp.log',status='unknown')
         read(65,*)cjunk,cjunk,cjunk,vtot
         if (ie.eq.0) vtot_0 = vtot
         ie = ie + 1
         ichecken=1
         close(65)
      endif
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)
c      write (6,*) 'gkey test',gkeyword,word3,word4
c      write (6,*) 'line test',line
      if(igkey.eq.'3') then
         IF (WORD3.EQ.gkeyword) THEN
            OPEN (unit=65,status='unknown')
            REWIND (65)
            WRITE (65,1000) WORD4
            REWIND (65)
            READ (65,*) vtot
            close(unit=65,status='keep')
            if (ie.eq.0) vtot_0 = vtot
            ie = ie + 1
            ichecken=1
         ENDIF
      endif
c  lc :  This is the read energy even for a single point
c        write(107,*) 'The read vtot is',vtot

 1000       FORMAT (A30)

c  lc:   comment this for handle strange behaviour with 
c        combined computation
c check for early end of job
      if (WORD.eq.'JOB') then
         write(*,*) 'lc debug, line 727'
         goto 9000
      endif

cc read and write gradient
      IF(
     &    ((WORD.EQ.'VARIABLE').AND.(WORD2.EQ.'OLD')) .or.
     &    ((word3 .eq. 'GRADIENT').and.(WORD4.EQ.'DISPLACEMENT'))
     &    )
     &    then
         gradval=0.
         OPEN (unit=64,file='grad.dat',status='unknown')
         if (.not.EFalgo) CALL LineRead2 (11)
         nint = natom*3-6
         write(64,*)nint
         if (natom.eq.2) nint=1
c  lc:   EF algorithm handling, since the print is changed for 
c        this opt algo, now the grad reading is changed to 
c        the right number of coordinates displayed that are 
c        equal to nint-ircons that are the number of frozen
c        coordinates
         if (EFalgo) nint=nint-ircons
         do iread = 1 , nint
            read(11,*)cjunk,cjunk1,gradval
            if (.not.EFalgo) then 
               grad(iread)=-gradval
               write(64,*)cjunk,grad(iread)
            else 
               grad(iread)=gradval
               write(64,*)intcoor(iread),grad(iread)
            endif
         enddo
         if (EFalgo) then 
            do iwrite=0,ircons-1
               write(64,*) intcoor(nint+ircons-iwrite),0.00000000
            end do
         end if
         close(64)
      ENDIF
      IF(WORD.EQ.'CENTER'.AND.WORD2.EQ.'ATOMIC'.AND.
     $ WORD3.eq.'FORCES')then
         gradv1=0.
         gradv2=0.
         gradv3=0.
         OPEN (unit=64,file='grad_xyz.dat',status='unknown')
         CALL LineRead2 (11)
         CALL LineRead2 (11)
         nint = natom
         write(64,*)nint
         do iread = 1 , nint
            read(11,*)cjunk,cjunk1,gradv1,gradv2,gradv3
            write(64,*)cjunk1,gradv1,gradv2,gradv3
         enddo
         close(64)
      ENDIF

c lc : if the lc_flag is greater than 0 it means that
c      opt in route section of gauss.com is detected
c      so these check on the optimization failure are 
c      neeeded
c  lc control, write lc flag in LC_output.dat
c      write(107,*) 'LC flag is=',lc_flag
      if (lc_flag.gt.0) then
      IF(
     &   (
     &   (WORD.EQ.'OPTIMIZATION').AND.(WORD2.EQ.'STOPPED.').AND.
     $   (IRES.NE.0).AND.(ixyz.ne.1)
     &   )
     &   .OR.
     &   (
     &   (WORD2.EQ.'OPTIMIZATION').AND.(WORD3.EQ.'STOPPED').AND.
     &   (IRES.NE.0).AND.(ixyz.ne.1)
     &    )
     &  ) THEN
 412     CONTINUE
         write(*,*) 'Optimization stopped?'
         CALL LineRead2 (11)
         if (WORD.eq.'JOB') then
            write(*,*) 'lc debug, line 775'
            go to 9000
         endif
c        write (6,*) 'g09test2',word,natom
         if ((WORD.eq.'NORMAL').and.(natom.eq.1)) then
            write(*,*) 'Debug pring line 854'
            go to 9100
         endif
         if (word2.eq.'NAME'.and.word3.eq.'VALUE')then
            CALL LineRead2 (11)
            nint = natom*3-6
            if (natom.eq.2) nint=1
            do iread = 1 , nint
               CALL LineRead2 (11)
c              write (6,*) 'iread in g09 test',iread,nint,ntau,word,
c    $ word2 
               OPEN (unit=64,status='unknown')
               REWIND (64)
               WRITE (64,910) word3,word2,word
c              write (6,*) 'word test',word3,word2,word
               REWIND (64)
c              write (6,*) 'iread test',iread,natom,ntau
cc modified index_nf if (iread.le.nint) then
               if (iread.le.(nint-ntau)) then
                  READ (64,*) xint(iread)
               else
                  itau = iread-(natom*3-6-ntau)
c                 write (6,*) 'itau test',itau
                  READ (64,*) tauopt(itau)
                  if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $             tauopt(itau) - 360.0d0
                  if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $             tauopt(itau) + 360.0d0
               endif
               close (unit=64,status='keep')
            enddo
            write (6,*) 'tauopt',(tauopt(itau),itau=1,ntau)
            go to 414
          endif
         if (word.eq.'CENTER'.and.word2.eq.'ATOMIC'.and.word3.eq.
     $       'ATOMIC'.and.natom.ne.1) then
c            write(*,*)'passing from here'
c            write(*,*)'natomt is',natomt
            READ (11,*)
            READ (11,*)
            do j=1,natom
               READ (11,*)ijunk,iatype,ijunk,coox(j),cooy(j),cooz(j)
            enddo
         endif
         go to 412
      endif
      endif
 414  continue

cc check for failure due to small interatomic distances
cc the code will not stop as this can happen for large molecules

      if (WORD.EQ.'SMALL'.AND.WORD2.EQ.'INTERATOMIC') then
         ichecken=1
         write(*,*) 'lc debug, line 833'
         goto 9000
      endif

c check for optimization and then read coordinates and frequencies
cc     modified to  consider also job 'completed on the basis of 
cc     negligible forces'
cc lc:   for ef algo, the keyword is **  CONVERGENCE 
      IF (
     &(WORD2.EQ.'COMPLETED.'.OR.WORD2.EQ.'COMPLETED') .OR.
     &(WORD2.EQ.'CONVERGENCE'.AND. WORD5.EQ.'SATISFIED')
     & ) THEN
c  This part relies on the completed optimization, the problem 
c  relies on single point calculations
         write(*,*) 'Optimization convergence criterion met'
212      continue
         CALL LineRead2 (11)
c        write (6,*) 'word in g09 test',word,word2,word3
c        write (6,*) 'g09test2',word,natom
         if ((WORD.eq.'NORMAL').and.(natom.eq.1)) then
            write(*,*) 'Debug print line 928'
            go to 9100
         endif
         if (WORD.eq.'JOB'.and.natom.ne.1) then
            write(*,*) 'lc debug, line 851'
            go to 9000
         endif
c        write (6,*) 'line2 test',line
c        write (6,*) 'words test',word2,word3
         if (word2.eq.'NAME'.and.word3.eq.'VALUE'.and.ixyz.ne.1) then
            CALL LineRead2 (11)
            nint = natom*3-6
            if (natom.eq.2) nint=1
c  lc:   Debug for EigenvalueFollow optimization algorithm start
            if (.not. EFalgo) then
               do iread = 1 , nint
                  CALL LineRead2 (11)
                  OPEN (unit=64,status='unknown')
                  REWIND (64)
                  WRITE (64,910) word3,word2,word
                  WRITE (*,910) word3,word2,word
910              FORMAT (3A20)
                  REWIND (64)
                  if (iread.le.(nint-ntau)) then
                     READ (64,*) xint(iread)
                  else
                     itau = iread-(natom*3-6-ntau)
                     READ (64,*) tauopt(itau)
                     if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $             tauopt(itau) - 360.0d0
                     if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $             tauopt(itau) + 360.0d0
                  endif
                  close (unit=64,status='keep')
               enddo
c  lc :  debugging for eigenfollow algorithm, it must write
c        also the freezed internal coordinates. So if
c        iread is equal to nint that is setted to the number
c        of internal coordinates minus the number of freezed coordinates
c        so in this case it must write 
            else 
               do iread = 1 , nint
                  CALL LineRead2 (11)
                  OPEN (unit=64,status='unknown')
                  REWIND (64)
                  if (iread .le. (nint-ircons)) then 
                     WRITE (64,'(3(A20))') word3,word2,word
                     WRITE (*,'(3(A20))') word3,word2,word
                  else 
                     read(157,'(A)') vlstring
                     write(64,*) trim(adjustl(vlstring))//'  !'
                     write(*,*)  trim(vlstring)//'  !'
                  end if 
                  rewind(64)
                  if (iread.le.(nint-ntau)) then
                     READ (64,*) xint(iread)
                  else
                     itau = iread-(natom*3-6-ntau)
                     READ (64,*) tauopt(itau)
                     if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $             tauopt(itau) - 360.0d0
                     if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $             tauopt(itau) + 360.0d0
                  endif
                  close (unit=64,status='keep')
               enddo
               close(157)
            end if 
c  lc:   Debug for EigenvalueFollow optimization algorithm end
            write (6,*) 'tauopt',(tauopt(itau),itau=1,ntau)
         endif
         if (word.eq.'CENTER'.and.word2.eq.'ATOMIC'.and.word3.eq.
     $       'ATOMIC') then
c     $       'ATOMIC'.and.natom.ne.1) then
c            write(*,*)'passing from here'
c            write(*,*)'natomt is',natomt
c            stop
            open(unit=64,FILE='geom.xyz',status='unknown')
            READ (11,*)
            READ (11,*)
            index_nf=0

c            if(iabs.eq.1.and.ired.eq.0) then
            if(ired.eq.0) then
               natomread=natomt
            else
               natomread=natom
            endif
            if(natom.gt.50)then
               natomread=natom
            endif
c            do j=1,natomt
c            write(*,*)'check natom is ',natomread
            do j=1,natomread
               READ (11,*)ijunk,iatype,ijunk,coox(j),cooy(j),cooz(j)
c            write(*,*)'check ',ijunk,iatype,ijunk
               if (iatype.eq.1)  aname='H'
               if (iatype.eq.2)  aname='He'
               if (iatype.eq.3)  aname='Li'
               if (iatype.eq.4)  aname='Be'
               if (iatype.eq.5)  aname='B'
               if (iatype.eq.6)  aname='C'
               if (iatype.eq.7)  aname='N'
               if (iatype.eq.8)  aname='O'
               if (iatype.eq.9)  aname='F'
               if (iatype.eq.10) aname='Ne'
               if (iatype.eq.11) aname='Na'
               if (iatype.eq.12) aname='Mg'
               if (iatype.eq.13) aname='Al'
               if (iatype.eq.14) aname='Si'
               if (iatype.eq.15) aname='P'
               if (iatype.eq.16) aname='S'
               if (iatype.eq.17) aname='Cl'
               if (iatype.eq.18) aname='Ar'
               if (iatype.eq.19) aname='K'
               if (iatype.eq.20) aname='Ca'
               if (iatype.eq.31) aname='Ga'
               if (iatype.eq.32) aname='Ge'
               if (iatype.eq.33) aname='As'
               if (iatype.eq.34) aname='Se'
               if (iatype.eq.35) aname='Br'
               if (iatype.eq.36) aname='Kr'
               if (iatype.eq.45) aname='Rh'
               if (iatype.eq.49) aname='In'
               if (iatype.eq.50) aname='Sn'
               if (iatype.eq.51) aname='Sb'
               if (iatype.eq.52) aname='Te'
               if (iatype.eq.53) aname='I'
               if (iatype.eq.54) aname='Xe'
               if(iatype.ne.-1)then
                  index_nf=index_nf+1
                  write(64,1204)aname, coox(j),cooy(j),cooz(j)
                  coord(j,1)=coox(j)
                  coord(j,2)=cooy(j)
                  coord(j,3)=cooz(j)
                  if(index_nf.eq.natom) goto 1222
               endif
            enddo
 1222       continue
            write(64,*)
            close(64)
            goto 264
         endif

 1204    format(1X,A2,1X,F9.5,1X,F9.5,1X,F9.5)

         go to 212
         if(ixyz.eq.1) goto 264


c 214     continue
c        write (6,*) 'starting coordinate determination'

c         CALL LineRead (11)
c         if ((word.eq.'NUMBER').and.(word2.eq.'NUMBER')) then
c            call LineRead (11)
c            do iatom = 1 , natomt
c               read (11,*) it1,it2,it3,(coord(iatom,idim),idim=1,3)
c            enddo

 264     continue
         if(natom.eq.1) goto 9100

         call LineRead (11)
         if (word.eq.'ROTATIONAL') then
            write(*,*) 'Rotational constants check'
            OPEN (unit=64,status='unknown')
            REWIND (64)
            WRITE (64,*) line
            close (64)
         command1="sed -ie 's/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/ 999999/g'
     +        fort.64"
            call commrun(command1)
            OPEN (unit=64,status='unknown') 
            if(ilin.eq.0) then
               read (64,*) ct1,ct2,ct3,(abcrot(idim),idim=1,3)
            else if (ilin.eq.1) then
               read (64,*) ct1,ct2,ct3,ct4,abcrot(1)
            endif
            close (unit=64,status='keep')
            write(*,*) 'ifreq is   ',ifreq
            if ((ifreq.eq.1) .or. (OptFreq.eqv..true.)) then
               jfreqtot=3*natom-6
               if (natom.eq.2) jfreqtot=1
               jfreq0=0
c  lc : mod start anharmonic frequencies reading from here START
c  lc:   in this case the anharmonic fundamental bands are taken
c        with a simple bash command that is composed by a format
c        specification  
c        the command is 
c        grep -i -m1 -A(jfreqtot+2) fundamental bands geom.log | 
c        tail -n(jfreqtot) | sed 's/[H,L]//g' | awk '{print $4}' 

               if (AnhFlag .ne. 0) then
                  write(*,*) 'Anharmonic frequencies reading activated'
                  write(longcommand1,987) jfreqtot+2,jfreqtot
                  call commrun(longcommand1) 
                  open(977,file='Anharmonic_frequencies.dat')
                  do AnhCounter=1,jfreqtot 
                     read(977,*) freq(AnhCounter)
                  end do 
                  close(977)
                  goto 9100                  
 987        FORMAT(                                                    
     &   "grep -i -m1 -A",I0,1x," 'fundamental bands' geom.log |",1x,
     &   "tail -n",I0,1x,"| sed 's/[H,L]//g' |",1x, 
     &   " awk '{print $4}' >",1x,
     &   "Anharmonic_frequencies.dat"
     &       ) 
               end if

C  lc: mod start anharmonic frequencies reading from here END
               write(*,*) 'line debugging 1138'
 314           continue
               CALL LineRead (11)
               if (WORD.eq.'FREQUENCIES') then
                  write(*,*) 'Debug print, frequencies checked 1'
                  OPEN (unit=65,status='unknown')
                  REWIND (65)
                  if (natom.eq.2) then 
                     WRITE (65,1009) word3
 1009                FORMAT (A20)
                     REWIND (65)
                     READ (65,*) freq(1)
                     close (unit=65,status='keep')
                     jfreq0 = jfreq0+1
                  else
                     WRITE (65,1010) word3,WORD4,word5
 1010                FORMAT (3A20)
                     REWIND (65)
                     READ (65,*) (freq(jfreq0+jfreq),jfreq=1,3)
                     close (unit=65,status='keep')
                     jfreq0 = jfreq0+3
                  endif
               ENDIF
               if (jfreq0.lt.jfreqtot) go to 314
               if (jfreq0.eq.jfreqtot.and.ilin_fr.eq.1) then
 315              continue
                  CALL LineRead (11)
                  if (WORD.eq.'FREQUENCIES') then
                  write(*,*) 'Debug print, frequencies checked 2'
                     jfreq0 = jfreq0+1
                     OPEN (unit=65,status='unknown')
                     REWIND (65)
                     WRITE (65,1011) word3
 1011                FORMAT (A20)
                     REWIND (65)
                     READ (65,*) freq(jfreq0)
                     close (unit=65,status='keep')
                     go to 9100
                  endif
                  go to 315
               endif
               goto 9100
            else
               goto 9100
            endif
         else
c  lc : this goto 264 sets the cycle to go up to 
c        rotational keyword search. The frequencies
c        parsing is activated if rotational is checked
            go to 264
         endif
      endif 
      go to 114
      
 9000 continue
      write (*,*) 'ERROR: Complete convergence failure V set to 1.0d20'
c      stop
c  lc :  arrivi a questa flag se la opt fallisce 
      if (lc_flag .gt. 0) then
         vtot = 2.0d20
         irepeat=irepeat+1
         write(6,*)'ires = ', ires
         write(6,*)'irepeat = ', irepeat
         if(irepeat.le.ires) then
cc       check if the code died for scf convergence

            iqc_term=0
            open (unit=99,status='unknown')
            rewind (99)
            write (99,1919)
            rewind (99)
            read (99,1920) command1
            close (99)
            call commrun(command1)
            open (unit=99,file='tmp2.dat',status='unknown')
            write(99,*)'qc file end'
            close (99)
            command1='cat tmp1.dat tmp2.dat > tmp.dat'
            call commrun(command1)
            open (unit=99,file='tmp.dat',status='unknown')
            read(99,*)cjunk,cjunk,iqclab
            if(iqclab.ne.'end')then
               if(iqclab.eq.'criterion') then 
                  iqc=1
                  iqc_term=1
                  write(6,*)'setting iqc  = 1 '
               endif
c            read(99,*)cjunk,cjunk,iqclab
            endif
            close (99)

 1919    format("egrep 'Convergence criterion not met' geom.log > 
     +   tmp1.dat")
 1920 format (A300)

cc   check if error in internal coord system and switch to cartesian
cc   if vdw or mmodred activated.

            icsy_term=0
            open (unit=99,status='unknown')
            rewind (99)
            write (99,2919)
            rewind (99)
            read (99,1920) command1
            close (99)
            call commrun(command1)
            open (unit=99,file='tmp2.dat',status='unknown')
            write(99,*)'csy file end'
            close (99)
            command1='cat tmp1.dat tmp2.dat > tmp.dat'
            call commrun(command1)
            open (unit=99,file='tmp.dat',status='unknown')
            read(99,*)cjunk,cjunk,iqclab
            if(iqclab.ne.'end')then
               if(iqclab.eq.'internal'.and.ired.eq.1) then 
                  icsy=1
                  icsy_term=1
                  ired=2
                  if(irecov.ne.1)irepeat=irepeat-1
                  write(6,*)'setting icsy  = 1 '
                  write(6,*)'setting ired  = 2 '
               endif
c            read(99,*)cjunk,cjunk,iqclab
            endif

            close (99)
            open (unit=99,status='unknown')
            rewind (99)
            write (99,2920)
            rewind (99)
            read (99,1920) command1
            close (99)
            call commrun(command1)
            open (unit=99,file='tmp2.dat',status='unknown')
            write(99,*)'csy file end'
            close (99)
            command1='cat tmp1.dat tmp2.dat > tmp.dat'
            call commrun(command1)
            open (unit=99,file='tmp.dat',status='unknown')
            read(99,*)cjunk,cjunk,iqclab
c         write(*,*)'iqclab is ',iqclab
            if(iqclab.ne.'end')then
               if(iqclab.eq.'Z-matrix'.and.ired.eq.1) then 
                  icsy=1
                  icsy_term=1
                  ired=2
                  if(irecov.ne.1)irepeat=irepeat-1
                  write(6,*)'setting icsy  = 1 '
                  write(6,*)'setting ired  = 2 '
               endif
c            read(99,*)cjunk,cjunk,iqclab
            endif

            close (99)

 2919    format("egrep 'Error in internal coordinate' geom.log > 
     +   tmp1.dat")
 2920    format("egrep 'Error on Z-matrix card' geom.log > 
     +   tmp1.dat")

            if(ired.eq.1.and.iqc_term.ne.1)then
c         if(ired.eq.1)then
            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,coox,cooy,cooz,xint,tauopt,
     $       ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)
            if(xint(1).eq.0)then
               write(*,*)'failed in updating geometry'
               write(*,*)'probably error of g09, check output'
               write(*,*)'stopping now'
               stop
            endif
         endif
cc         
            goto 999
         endif
      end if
 9100 continue
      if(ifreq.eq.1)then
         write(*,'(A,/,50("-"))') 'Frequencies'
         write(*,'(f0.8,/)')  freq(1:jfreqtot)
      endif
      close(unit=11,status='keep')
cc      if(ichecken.eq.0.and.gkeyword.ne.'level1')then
      if(ired.eq.1)then
         call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $ ,idconn,bname,anname,dname,atname,coox,cooy,cooz,xint,tauopt,
     $  ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)
         if(xint(1).eq.0)then
            write(*,*)'failed in updating geometry'
            write(*,*)'probably error of g09, check output'
            write(*,*)'stopping now'
            stop
         endif
      endif
      

      if(ichecken.eq.0)then
         write(*,*)'did not found energy in the output'
         write(*,*)'or error of g09'
         write(*,*)'check you have indicated the proper gkey'
         write(*,*)'or if g09 has failed'
         stop
      endif
      write(*,*)'out of g09opt' 
      RETURN
      END 


      SUBROUTINE read_g09_ircout(force_con,
     $  natom,iread,numpointsf,numpointsb,atgeom_me,
     $  ifcread,rc_ene,rc_coord,grad,ionlyfor)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c     parameter (ntaumx=10, nmdmx=300, natommx=100, ndim=3)

c
      dimension rc_ene(nircmx)
      dimension rc_coord(nircmx)
 
      character*180 command1

      character*80 atgeom_me(natommx,nircmx)
      character*80 force_con(3*natommx*3*natommx/10,nircmx)
      character*80 grad(natommx,nircmx)
      character*2 aname
      character*20 filename
      character*20 step(nircmx)
      character*30 cjunk

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      character*20 ct1,ct2,ct3
      include 'filcomm.f'
 
c
      call LineRead(0)

      iread=0
      imd=imdtunn
      write(*,*)'entering g09 IRC reading sub'

 100  continue
      if(iread.eq.0)then
         command1='cp -f ./irc_files/irc_g09f.log 
     $            ./irc_g09.log'
         call commrun(command1)
      else if (iread.eq.1.and.ionlyfor.ne.1)then
         command1='cp -f ./irc_files/irc_g09b.log 
     $            ./irc_g09.log'
         call commrun(command1)
      endif

c         write(*,*)'ok up to here cc3'
c         stop


      open (unit=108,status='unknown')
      write (108,1200)natom+6
      rewind (108)
      read (108,1030) command1
      close (108)
 1200 format (" egrep -A"I0.3 " 'Input orientation|CHANGE' 
     +  ./irc_g09.log > ./temp.tmp")
      call commrun(command1)
      
      open (unit=108,status='unknown')
      write (108,1201)natom+6
      rewind (108)
      read (108,1030) command1
      close (108)
 1201 format (" egrep -B"I0.3 " 'CHANGE'
     +  ./temp.tmp > ./temp1.tmp")
      call commrun(command1)
      
      open (unit=108,status='unknown')
      write (108,1202)natom+3
      rewind (108)
      read (108,1030) command1
      close (108)      
 1202 format (" egrep -A"I0.3 " 'Number'
     +  ./temp1.tmp > ./temp2.tmp")
      call commrun(command1)

      open (unit=108,status='unknown')
      write (108,1206)
      rewind (108)
      read (108,1030) command1
      close (108)
 1206 format (" egrep -v 'Number|Distance|--' 
     +  ./temp2.tmp > ./str1.out")
      call commrun(command1)

      command1=' rm -f ./*tmp'
      call commrun(command1)

cc now determine number of IRC points 
      call LineRead (0)
      open (unit=108,file='./irc_g09.log',status='unknown')
      do while ((word.NE.'TOTAL').and.(word4.NE.'POINTS:'))
         call LineRead (108)
         if ((word.EQ.'TOTAL').and.(word4.EQ.'POINTS:')) then
         endif
      enddo

      read(word5,'(i5.2)') numpoints
cc      if(iread.eq.1)stop
cc      print *,  numpoints
      
      write(*,*)'the number of IRC points to read is ',numpoints
      if(iread.eq.0)numpointsf=numpoints
      if(iread.eq.1)numpointsb=numpoints

      close (108)

cc now read  str1 and convert it to structure.out, determining 

      aname='X'
      open (unit=108,file='./str1.out',status='unknown')
      open (unit=109,file='./structure.out',status='unknown')
      do i=1,numpoints
         do j=1,natom
            read(108,*)ijunk,iatype,kjunk,xcoo,ycoo,zcoo
            if (iatype.eq.1)  aname='H'
            if (iatype.eq.2)  aname='He'
            if (iatype.eq.3)  aname='Li'
            if (iatype.eq.4)  aname='Be'
            if (iatype.eq.5)  aname='B'
            if (iatype.eq.6)  aname='C'
            if (iatype.eq.7)  aname='N'
            if (iatype.eq.8)  aname='O'
            if (iatype.eq.9)  aname='F'
            if (iatype.eq.10) aname='Ne'
            if (iatype.eq.11) aname='Na'
            if (iatype.eq.12) aname='Mg'
            if (iatype.eq.13) aname='Al'
            if (iatype.eq.14) aname='Si'
            if (iatype.eq.15) aname='P'
            if (iatype.eq.16) aname='S'
            if (iatype.eq.17) aname='Cl'
            if (iatype.eq.18) aname='Ar'
            if (iatype.eq.19) aname='K'
            if (iatype.eq.20) aname='Ca'
            if (iatype.eq.31) aname='Ga'
            if (iatype.eq.32) aname='Ge'
            if (iatype.eq.33) aname='As'
            if (iatype.eq.34) aname='Se'
            if (iatype.eq.35) aname='Br'
            if (iatype.eq.36) aname='Kr'
            if (iatype.eq.45) aname='Rh'
            if (iatype.eq.49) aname='In'
            if (iatype.eq.50) aname='Sn'
            if (iatype.eq.51) aname='Sb'
            if (iatype.eq.52) aname='Te'
            if (iatype.eq.53) aname='I'
            if (iatype.eq.54) aname='Xe'
            if (aname.eq.'X') then
               write(*,*)'undefined atom name in IRC'
               write(*,*)'modify the code'
               stop
            endif
            write(109,1203)aname,xcoo,ycoo,zcoo
         enddo
      enddo

 1203 format(1X,A2,1X,F9.5,1X,F9.5,1X,F9.5)

      close (108)
      close (109)

cc now read struct.out 

c vector of geometries for me input
      open (unit=108,file='./structure.out',status='unknown')
      do inumpoints = 1, numpoints 
         if(iread.eq.0) then
            index_nf=numpoints+1-inumpoints
         else if(iread.eq.1) then
            index_nf=numpointsf+inumpoints+1
         endif
         do iatom = 1, natom
            read (108,'(A80)') atgeom_me(iatom,index_nf)
c            write (*,*)iatom, atgeom_me(iatom,index_nf)
         enddo 
      enddo 
      close(108)
c
cc now get geometry of TS from me files
      open (unit=108,file='./me_files/ts_ge.me',status='unknown')
      read(108,*)
      read(108,*)
      read(108,*)
      do iatom = 1, natom
         read (108,'(A80)') atgeom_me(iatom,numpointsf+1)
      enddo 
      close(108)
c
c extract gradients and force constants
 
      igrad=natom+2
      open (unit=108,status='unknown')
      write (108,1208)igrad
      rewind (108)
      read (108,1030) command1
      close (108)      
 1208 format (" egrep -A"I0.3" '  Forces |CHANGE'
     +  ./irc_g09.log > grad_tmp1.dat ")
      call commrun(command1)

      open (unit=108,status='unknown')
      write (108,1303)igrad
      rewind (108)
      read (108,1030) command1
      close (108)
 1303 format (" egrep -B"I0.3 " 'CHANGE'
     +  grad_tmp1.dat > grad1.dat")
      call commrun(command1)

      ifcread=0
      ifcmax=natom*3/5
      do j=0,ifcmax
         ifcread=ifcread+natom*3-5*j+1
      enddo
      open (unit=108,status='unknown')
      write (108,1302)ifcread
      rewind (108)
      read (108,1401) command1
      close (108)      
 1302 format (" egrep -A"I0.3"
     $  'Force constants in Cartesian coordinates'
     $ ./irc_g09.log > force_const1.dat")
      call commrun(command1)
 1401 format (A300)

c force_constants for proj input
      open (unit=108,status='unknown')
      write (108,1211)
      rewind (108)
      read (108,1035) command1
      close (108)
 1211 format (" sed -e 's/D/E/g' force_const1.dat
     + > force_constants.dat")      
      call commrun(command1) 
 1035 format (A300)

c vector of gradients for proj input
      open (unit=108,file='grad1.dat',status='unknown')
      do inumpoints =1, numpoints
         read (108,*)
         if(iread.eq.0) then
            index_nf=numpoints+1-inumpoints
            do iatom = 1, natom
               read (108,'(A70)') grad(iatom,index_nf)
            enddo
         else if(iread.eq.1) then
            index_nf=numpointsf+inumpoints+1
            do iatom = 1, natom
               read (108,'(A70)') grad(iatom,index_nf)
            enddo
         endif
         read (108,*)
         read (108,*)
         if(inumpoints.ne.numpoints)read (108,*)
      enddo
      close(108)
c      write(*,*)'ok 1'

c initialize ts gradient to two previous points
cc it is ok for projection, be careful for tunneling

      if(numpointsf.ne.1)then
         if(imd.ne.1) then
            do iatom = 1, natom
               grad(iatom,numpointsf)=grad(iatom,numpointsf-1)
               grad(iatom,numpointsf+1)=grad(iatom,numpointsf-1)
               grad(iatom,numpointsf+2)=grad(iatom,numpointsf+3)
cc         grad(iatom,numpointsf+1)='1  1  0.0 0.0 0.0'
cc         grad(iatom,numpointsf+1)=grad(iatom,numpointsf)
            enddo
         else if(imd.eq.1) then
            do iatom = 1, natom
c               grad(iatom,numpointsf+1)='1  1  0.0 0.0 0.0'
               grad(iatom,numpointsf+1)=grad(iatom,numpointsf)
            enddo
         endif
      endif

c vector of force constants
      open (unit=108,file='force_constants.dat',status='unknown')  

      do inumpoints = 1, numpoints+1 
c         write(*,*)'inumpoints ',inumpoints
         read (108,*)      
         if(iread.eq.0) then
            index_nf=numpoints+2-inumpoints           
            do inumlines = 1, ifcread
               read (108,'(A77)') force_con(inumlines,index_nf)
c               write (*,*) force_con(inumlines,index_nf)
c               write (*,*) 'line ',inumlines
c               write(*,*)'inumpoints ',inumpoints
            enddo
         else if (iread.eq.1) then
            index_nf=numpointsf+inumpoints
            do inumlines = 1, ifcread
               read (108,'(A77)') force_con(inumlines,index_nf)
            enddo
         endif
         if (inumpoints.ne.numpoints+1)then
            read (108,*)
         endif
      enddo
      close(108)

c now read energies

cc first get energy of TS            
      open (unit=108,status='unknown')
      write (108,1209)
      rewind (108)
      read (108,1035) command1
      close (108)
 1209 format (" egrep 'Energies reported' 
     +  ./irc_g09.log > ./tsen_irc_level.tmp")
      call commrun(command1)

      open (unit=108,status='unknown')
      write (108,1210)
      rewind (108)
      read (108,1035) command1
      close (108)
 1210 format (" awk '{print $9}' ./tsen_irc_level.tmp 
     +   > ./tsen.out")
      call commrun(command1)

c      write(*,*)'ok 3'

cc now get energies along reaction coordinate

      open (unit=108,status='unknown')
      write (108,1230)numpoints+1
      rewind (108)
      read (108,1030) command1
      close (108)      
 1230 format (" egrep -A"I0.3 " 'Energy    RxCoord' 
     +  ./irc_g09.log > ./energy1a.tmp")
      call commrun(command1)

c      write(*,*)'ok 4'

      open (unit=108,status='unknown')
      write (108,1231)numpoints+1
      rewind (108)
      read (108,1030) command1
      close (108)      
 1231 format (" egrep -A"I0.3 " 'Energy   Rx Coord' 
     +  ./irc_g09.log > ./energy1b.tmp")
      call commrun(command1)

      command1='cat ./energy1a.tmp ./energy1b.tmp 
     +        > ./energy1.tmp'
      call commrun(command1)

      write(*,*)'numpoint is',numpoints
      open (unit=108,file='./energy1.tmp',status='unknown')
      read (108,*)cjunk
      write (*,*)cjunk
      do inumpoints = 1, numpoints+1 
        call LineRead(108)
        read(word2,'(f9.6)') en_irc
        read(word3,'(f9.6)') coord_irc
        write(*,*)'en irc', en_irc
        write(*,*)'word2 is ', word2
        write(*,*)'coord irc', coord_irc
        write(*,*)'word3 is ', word2
         if(iread.eq.0) then
            index_nf=numpoints+2-inumpoints           
         else if (iread.eq.1) then
            index_nf=numpointsf+numpointsb-inumpoints+2
         endif
         write(*,*)'index_nf is ',index_nf
         write(*,*)'nircmx is  ',nircmx
         rc_ene(index_nf) = en_irc
         rc_coord(index_nf) = coord_irc
      enddo
      close(108)

      if(iread.eq.0.and.ionlyfor.ne.1) then
         iread=1
         goto 100
      endif
 1030 format (A300)

      write(*,*)'out of g09 IRC reading sub'

      return
      end

C     *****************************
      subroutine en_key_g09(comline1,igkey,gkeyword,ispin)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*30 gkeyword,igkey
      character*(*) comline1
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      open(unit=99,file='temp.tmp',status='unknown')
      write(99,*)comline1
      close(99)

      command1="sed -ie 's/\// /' temp.tmp"
      call commrun(command1)

      open(unit=99,file='temp.tmp',status='unknown')
      call LineRead(99)
      close(99)

      gkeyword='notfound'
      igkey='notfound'

      if(word.eq.'WB97XD')then
         if(ispin.eq.1)then
            gkeyword='E(RWB97XD)'
         else
            gkeyword='E(UWB97XD)'
         endif
         igkey='3'
      else if(word.eq.'UWB97XD')then
         gkeyword='E(UWB97XD)'
         igkey='3'
      else if(word.eq.'RWB97XD')then
         gkeyword='E(RWB97XD)'
         igkey='3'
      endif

      if(word.eq.'M062X')then
         if(ispin.eq.1)then
            gkeyword='E(RM062X)'
         else
            gkeyword='E(UM062X)'
         endif
         igkey='3'
      else if(word.eq.'UM062X')then
         gkeyword='E(UM062X)'
         igkey='3'
      else if(word.eq.'RM062X')then
         gkeyword='E(RM062X)'
         igkey='3'
      endif

      if(word.eq.'B3LYP')then
         if(ispin.eq.1)then
            gkeyword='E(RB3LYP)'
         else
            gkeyword='E(UB3LYP)'
         endif
         igkey='3'
      else if(word.eq.'UB3LYP')then
         gkeyword='E(UB3LYP)'
         igkey='3'
      else if(word.eq.'RB3LYP')then
         gkeyword='E(RB3LYP)'
         igkey='3'
      endif

      if(word.eq.'B2PLYP')then
         if(ispin.eq.1)then
            gkeyword='E(B2PLYP)'
         else
            gkeyword='E(B2PLYP)'
         endif
         igkey='4'
      else if(word.eq.'UB2PLYP')then
         gkeyword='E(B2PLYP)'
         igkey='4'
      else if(word.eq.'RB3LYP')then
         gkeyword='E(B2PLYP)'
         igkey='4'
      endif

      if(word.eq.'B2PLYPD3')then
         if(ispin.eq.1)then
            gkeyword='E(B2PLYPD3)'
         else
            gkeyword='E(B2PLYPD3)'
         endif
         igkey='4'
      else if(word.eq.'UB2PLYPD3')then
         gkeyword='E(B2PLYPD3)'
         igkey='4'
      else if(word.eq.'RB3LYPD3')then
         gkeyword='E(B2PLYPD3)'
         igkey='4'
      endif

      return
      end


C     *****************************

      subroutine int_key_g09(comline1,ired_check)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c      character*30 gkeyword,igkey
      character*(*) comline1
c      character*300 comline1
      character*300 comline1_up
c     character*120 command1
c      LOGICAL leof,lsec,ltit

c      CHARACTER*1000 line,string
c      CHARACTER*160 sename,word,word2,word3
c     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'

      ired_check=0

c      command1="sed -ie 's/,/ /g' temp.tmp"
c      call commrun(command1)
c      command1="sed -ie 's/(/ /g' temp.tmp"
c      call commrun(command1)
c      open(unit=99,file='temp.tmp',status='unknown')
c      call LineRead(99)

c   1) Convert to upper case, case insensitive research
      comline1_up=comline1
      call upcase2(comline1_up)

c   2) Check if the route contains intern and opt togheter
      if (
     &    (index(comline1_up,'INTERN') .eq. 0)
     &     .and.
     &    (index(comline1_up,'GEOM=CHECK') .eq. 0)
     &     .and.
     &    (index(comline1_up,'OPT') .gt. 0)) then
         ired_check=1
      endif

c      write(*,*)'ired is ',ired_check
c      stop

      return
      end

C     *****************************
      subroutine comline34_g09(ispecies,comline1,comline2,
     $                         comline3,comline4)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*30 gkeyword,igkey
      character*(*) comline1
      character*(*) comline2
      character*(*) comline3
      character*(*) comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      open(unit=99,file='temp.tmp',status='unknown')
      write(99,*)comline1
      rewind(99)
      call Lineread(99)
      close(99)
      comline3=word
      write(*,*)'comline 3 is ',comline3


      open(unit=99,file='temp.tmp',status='unknown')
      write(99,*)comline1
      write(99,*)comline2
      close(99)
      command1="sed -ie 's/=/ /g' temp.tmp"
      call commrun(command1)

      ifine=0
      open(unit=99,file='temp.tmp',status='unknown')
      call Lineread(99)
      if(word.eq.'ULTRAFINE'.or.word.eq.'ULTRA')ifine=1
      if(word2.eq.'ULTRAFINE'.or.word2.eq.'ULTRA')ifine=1
      if(word3.eq.'ULTRAFINE'.or.word3.eq.'ULTRA')ifine=1
      if(word4.eq.'ULTRAFINE'.or.word4.eq.'ULTRA')ifine=1
      if(word5.eq.'ULTRAFINE'.or.word5.eq.'ULTRA')ifine=1
      if(word6.eq.'ULTRAFINE'.or.word6.eq.'ULTRA')ifine=1
      if(word7.eq.'ULTRAFINE'.or.word7.eq.'ULTRA')ifine=1

      call Lineread(99)
      if(word.eq.'ULTRAFINE'.or.word.eq.'ULTRA')ifine=1
      if(word2.eq.'ULTRAFINE'.or.word2.eq.'ULTRA')ifine=1
      if(word3.eq.'ULTRAFINE'.or.word3.eq.'ULTRA')ifine=1
      if(word4.eq.'ULTRAFINE'.or.word4.eq.'ULTRA')ifine=1
      if(word5.eq.'ULTRAFINE'.or.word5.eq.'ULTRA')ifine=1
      if(word6.eq.'ULTRAFINE'.or.word6.eq.'ULTRA')ifine=1
      if(word7.eq.'ULTRAFINE'.or.word7.eq.'ULTRA')ifine=1

      close(99)

c      stop
      if(ispecies.eq.0) then
         if(ifine.eq.0) then
            comline4='opt(calcfc,ts,maxcyc=1,noeigentest)  
     $ iop(7/33=1) guess=read geom=check'
         else
          comline4='opt(calcfc,ts,maxcyc=1,noeigentest)  
     $ iop(7/33=1) guess=read geom=check int=ultra 
     $'
         endif
      else
         if(ifine.eq.0) then
            comline4='opt(calcfc,maxcyc=1,noeigentest) iop(7/33=1) 
     $ guess=read geom=check'
         else
      comline4='opt(calcfc,maxcyc=1,noeigentest) iop(7/33=1)  
     $ guess=read geom=check int=ultra'
         endif
      endif

      return
      end

C     *****************************
      subroutine comline56_g09(ispecies,comline1,comline2,comline5,
     +   comline6)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*30 gkeyword,igkey
      character*(*) comline1
      character*(*) comline2
      character*(*) comline5
      character*(*) comline6
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      open(unit=99,file='temp.tmp',status='unknown')
      write(99,*)comline1
      write(99,*)comline2
      close(99)
      command1="sed -ie 's/=/ /g' temp.tmp"
      call commrun(command1)

      ifine=0
      open(unit=99,file='temp.tmp',status='unknown')
      call Lineread(99)
      if(word.eq.'ULTRAFINE'.or.word.eq.'ULTRA')ifine=1
      if(word2.eq.'ULTRAFINE'.or.word2.eq.'ULTRA')ifine=1
      if(word3.eq.'ULTRAFINE'.or.word3.eq.'ULTRA')ifine=1
      if(word4.eq.'ULTRAFINE'.or.word4.eq.'ULTRA')ifine=1
      if(word5.eq.'ULTRAFINE'.or.word5.eq.'ULTRA')ifine=1
      if(word6.eq.'ULTRAFINE'.or.word6.eq.'ULTRA')ifine=1
      if(word7.eq.'ULTRAFINE'.or.word7.eq.'ULTRA')ifine=1

      call Lineread(99)
      if(word.eq.'ULTRAFINE'.or.word.eq.'ULTRA')ifine=1
      if(word2.eq.'ULTRAFINE'.or.word2.eq.'ULTRA')ifine=1
      if(word3.eq.'ULTRAFINE'.or.word3.eq.'ULTRA')ifine=1
      if(word4.eq.'ULTRAFINE'.or.word4.eq.'ULTRA')ifine=1
      if(word5.eq.'ULTRAFINE'.or.word5.eq.'ULTRA')ifine=1
      if(word6.eq.'ULTRAFINE'.or.word6.eq.'ULTRA')ifine=1
      if(word7.eq.'ULTRAFINE'.or.word7.eq.'ULTRA')ifine=1
      close(99)

      open(unit=99,file='temp.tmp',status='unknown')
      if(ifine.eq.0)then
         write(99,*)comline1
         rewind(99)
         call Lineread(99)
         comline5=word
      else if (ifine.eq.1)then
c      rewind (99)
         write(99,*)comline1
         rewind(99)
         call Lineread(99)
         rewind(99)
         write(99,100)word,' int=ultra'
         rewind(99)
         read(99,101)comline5
      endif
      close(99)
      write(*,*)'comline 5 is ',comline5
c      stop
      if(ispecies.eq.0) then
         comline6='opt(calcfc,ts,maxcycle=1,noeigentest) iop(7/33=1) 
     $ guess=read geom=check'
      else
         comline6='opt(calcfc,maxcycle=1) iop(7/33=1) guess=read 
     $ geom=check'
      endif

 100  format(A30,A20)
 101  format(A300)

      return
      end

C     *****************************
      subroutine readgrad_g09(nvar,grad)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension grad(3*natommx)

      character*30 gkeyword,igkey
      character*70 cjunk,cjunk1
      character*300 comline1
      character*300 comline2
      character*300 comline3
      character*300 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)

cc read gradient in internal coordinates

      IF((WORD.EQ.'VARIABLE').AND.(WORD2.EQ.'OLD'))then
         gradval=0.
c         OPEN (unit=64,file='grad.dat',status='unknown')
         CALL LineRead2 (11)
c         nint = natom*3-6
c         write(64,*)nint
c         if (natom.eq.2) nint=1
         do iread = 1 , nvar
            read(11,*)cjunk,cjunk1,gradval
            grad(iread)=-gradval
c            write(*,*)iread,grad(iread)
         enddo
c         close(64)
      ENDIF
c      write(*,*)'word is ',word
      if (WORD.eq.'JOB') then 
      write(*,*) 'lc debug, line 1941'
      goto 9000
      endif 
      goto 114
      
 9000 continue
      close(11)
c      stop

      return
      end

C     *****************************
      subroutine readxyzgrad_g09(natom,ianum,grad_xyz)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension ianum(3*natommx)
      dimension grad_xyz(3*natommx)

      character*30 gkeyword,igkey
      character*70 cjunk
      character*300 comline1
      character*300 comline2
      character*300 comline3
      character*300 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)
cc read and write gradient
      IF(WORD.EQ.'CENTER'.AND.WORD2.EQ.'ATOMIC'.AND.
     $ WORD3.eq.'FORCES')then
         gradv1=0.
         gradv2=0.
         gradv3=0.
c         OPEN (unit=64,file='grad_xyz.dat',status='unknown')
         CALL LineRead2 (11)
         CALL LineRead2 (11)
         ip=1
         do iread = 1 , natom
            read(11,*)cjunk,ianum(ip),grad_xyz(ip),grad_xyz(ip+1)
     + ,grad_xyz(ip+2)
            ianum(ip+1)=ianum(ip)
            ianum(ip+2)=ianum(ip)
 1          ip=ip+3
         enddo
c         stop
      ENDIF
      if (WORD.eq.'JOB') then 
         write(*,*) 'lc debug, line 1941'
         goto 9000
      endif
      goto 114
      
 9000 continue
      close(11)

      return
      end

C     *****************************
      subroutine readen_g09(energy,gkeyword,igkey)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension grad_xyz(3*natommx)

      character*30 gkeyword,igkey
      character*70 cjunk
      character*300 comline1
      character*300 comline2
      character*300 comline3
      character*300 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
c  lc : goto loop untile job reading
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)
cc read energy
c  lc debug : insert energy reading 
      
      if ((igkey.eq.'notfound').or.(igkey.eq.'4'))then
         if(ilevcode.eq.1)then
            command1='formchk09 tmp.chk'
         else if (ilevcode.eq.3)then
            command1='formchk09 tmp.chk'
         endif
         call commrun(command1) 
         open(unit=778,file='tmp.fchk',iostat=iStatus,status='old')
            if (iStatus .ne. 0 ) then 
               command1='formchk09 tmp.chk'
               call commrun(command1)  
            end if 
            close(778,status='keep')        
         command1="egrep 'Total Energy' tmp.fchk > temp.log"
         call commrun(command1)         
         open(unit=65,file='temp.log',status='unknown')
         read(65,*)cjunk,cjunk,cjunk,vtot
         energy = vtot
         close(65)

      end if

      if(igkey.eq.'3') then
         IF (WORD3.EQ.gkeyword) THEN
            OPEN (unit=65,status='unknown')
            REWIND (65)
            WRITE (65,1000) WORD4
            REWIND (65)
            READ (65,*) vtot
            close(unit=65,status='keep')
            energy = vtot
         ENDIF
      endif 
c         else if (igkey.eq.'4') then
c            IF (WORD4.EQ.gkeyword) THEN
c               OPEN (unit=65,status='unknown')
c               REWIND (65)
c               WRITE (65,1000) WORD5
c               REWIND (65)
c               READ (65,*) vtot
c               close(unit=65,status='keep')
c               energy = vtot
c            ENDIF
         if (WORD.eq.'JOB') then 
            write(*,*) 'lc debug, line 2089'
            goto 9000
         endif
         goto 114
      
 9000 continue

      close(11)

 1000 FORMAT (A30)

      return
      end
C     *****************************
      subroutine readhess_g09(natom,hess_xyz)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c      dimension ianum(3*natommx)
      dimension hess_xyz(3*natommx,3*natommx)

      character*30 gkeyword,igkey
      character*70 cj
      character*300 comline1
      character*300 comline2
      character*300 comline3
      character*300 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      numbl_col5=int(3*natom/5)
      numlbl_col=3*natom-numbl_col5*5
      if(numlbl_col.ne.0)numbl_col5=numbl_col5+1

c      write(*,*)'numbl1 ',numbl_col5
c      write(*,*)'numbl2 ',numlbl_col
      do i=1,3*natom
         do j=1,3*natom
            hess_xyz(i,j)=0.
         enddo
      enddo
c      stop

      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)
cc read and write gradient
      IF(WORD.EQ.'FORCE'.AND.WORD2.EQ.'CONSTANTS'.AND.
     $ WORD4.eq.'CARTESIAN')then

         irow=0
         jcol=0
         do j=1,numbl_col5
c         do j=1,1
            CALL LineRead2 (11)
            do i = 1+irow, natom*3
               if(i-irow.eq.1)then
                  read(11,*)cj,hess_xyz(i,j+jcol)
               else if(i-irow.eq.2)then
                  read(11,*)cj,hess_xyz(i,j+jcol),hess_xyz(i,j+1+jcol)
               else if(i-irow.eq.3)then
                  read(11,*)cj,hess_xyz(i,j+jcol),hess_xyz(i,j+1+jcol),
     + hess_xyz(i,j+2+jcol)
               else if(i-irow.eq.4)then
                  read(11,*)cj,hess_xyz(i,j+jcol),hess_xyz(i,j+1+jcol),
     + hess_xyz(i,j+2+jcol),hess_xyz(i,j+3+jcol)
               else
                  read(11,*)cj,hess_xyz(i,j+jcol),hess_xyz(i,j+1+jcol),
     + hess_xyz(i,j+2+jcol),hess_xyz(i,j+3+jcol),hess_xyz(i,j+4+jcol)
               endif
            enddo
            irow=irow+5
            jcol=jcol+4
         enddo
c         if(numlbl_col.ne.0)then
c            do i = 1+irow, natom*3
c
c            enddo
c         endif
      ENDIF
      if (WORD.eq.'JOB') then 
         write(*,*) 'lc debug, line 2182'
         goto 9000
      endif
      goto 114
      
 9000 continue
      close(11)

c      stop

      return
      end

C     *****************************
      subroutine readxyzgeom_g09(natom,idummy,coox,cooy,cooz)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension coox(natommx)
      dimension cooy(natommx)
      dimension cooz(natommx)
      dimension idummy(natommx)

      character*30 gkeyword,igkey
      character*70 cjunk
      character*300 comline1
      character*300 comline2
      character*300 comline3
      character*300 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      iread=0
      atype=0.
      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)
cc read and write gradient
      IF(WORD.EQ.'Z-MATRIX'.AND.WORD2.EQ.'ORIENTATION:')then
         iread=1
         CALL LineRead2 (11)
         CALL LineRead2 (11)
         CALL LineRead2 (11)
         CALL LineRead2 (11)
         do i = 1 , natom
            read(11,*)cjunk,atype,cjunk,coox(i),cooy(i),cooz(i)
c            write(*,*)'atype is',atype
            if(atype.eq.-1.)idummy(i)=1
c            write(*,*)'idummy is is',idummy(i)
         enddo
c         stop
      ENDIF
      if (WORD.eq.'JOB') then 
         write(*,*) 'lc debug, line 2247'
         goto 9000
      endif
      goto 114
      
 9000 continue
      if(iread.eq.0)then
         rewind(11)
 115     continue
         CALL LineRead (0)
         CALL LineRead (11)
cc read geometry
         IF(WORD.EQ.'INPUT'.AND.WORD2.EQ.'ORIENTATION:')then
            iread=1
            CALL LineRead2 (11)
            CALL LineRead2 (11)
            CALL LineRead2 (11)
            CALL LineRead2 (11)
            do i = 1 , natom
               read(11,*)cjunk,cjunk,cjunk,coox(i),cooy(i),cooz(i)
            enddo
c         stop
         ENDIF
         if (WORD.eq.'JOB') then 
            write(*,*) 'lc debug, line 2271'
            goto 9000
         endif
         goto 115
      endif
      close(11)

      return
      end

         function Up_to_low(vlstring)
         implicit none
         character(len=*),intent(in) :: vlstring 
         character(len=len_trim(vlstring)) :: Up_to_low
         integer ::  i             ! COUNTER 
         intrinsic :: achar
      
         Up_to_low=vlstring
         do i=1,len_trim(vlstring)      
            select case (vlstring(i:i))
            case ('A':'Z')
                     Up_to_low(i:i)=achar(  ichar(vlstring(i:i))+32  )
                  case default
                     cycle
            end select
         end do
      
      end function Up_to_low

C     *****************************
