c Electronic Structure to k(T,P)
c code to build TS from estoktp files
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.

c
      program tsbuild

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension dhindmn(nhindmx),dhindmx(nhindmx),
     $     nhindsteps(nhindmx),nhindsymm(nhindmx)

      LOGICAL leof,lsec,ltit
      CHARACTER*100 line,sename,string,word,word2,word3,title,title1,
     $ word4,word5,word6,word7
      character*100 command1
      character*40 label
      character*30 gmemlow,gmemhigh
      character*30 gmem
      character*10 hindlabread
      character*10 hindlab(nhindmx)
      include 'filcomm.f'

      open (unit=10,file='ts.dat',status='UNKNOWN')
      open (unit=11,file='tsr1.dat',status='UNKNOWN')
      open (unit=12,file='tsr2.dat',status='UNKNOWN')
      open (unit=13,file='reac1.dat',status='UNKNOWN')
      open (unit=14,file='reac2.dat',status='UNKNOWN')
      open (unit=15,file='estoktp.dat',status='UNKNOWN')
      open (unit=16,file='wellr.dat',status='UNKNOWN')
      open (unit=17,file='wellp.dat',status='UNKNOWN')
      open (unit=21,file='tsbuild.res',status='UNKNOWN')


cc first check reaction type
      rewind(15)
      do while (WORD.NE.'REACTIONTYPE')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (*,*) 'reaction type must be defined in estoktp.dat'
            write (21,*) 'reaction type must be defined in estoktp.dat'
            stop
         endif
      enddo
      if (WORD2.EQ.'ABSTRACTION') iabs = 1
      if (WORD2.EQ.'ADDITION') iadd = 1
      if (WORD2.EQ.'ISOMERIZATION') iiso = 1
      if (WORD2.EQ.'BETASCISSION') ibeta = 1

cc read if reac1 aor reac2 are linear

cc now read charge and spin from reac1 and reac2

      do while (WORD.NE.'NATOM')
         call LineRead (13)
         if (WORD.EQ.'END') then
            write (*,*) 'natom must be defined in reac1.dat'
            write (21,*) 'natom must be defined in reac1.dat'
            stop
         endif
      enddo
      read(13,*)natom1,natomt1,ilin1
      rewind(13)
      word=''

      do while (WORD.NE.'NATOM')
         call LineRead (14)
         if (WORD.EQ.'END') then
            write (*,*) 'natom must be defined in reac2.dat'
            write (21,*) 'natom must be defined in reac2.dat'
            stop
         endif
      enddo
      read(14,*)natom2,natomt2,ilin2
      rewind(14)
      word=''


c     read dihed infos

      do while (WORD.NE.'NHIND')
         call LineRead (11)
c         write(*,*)'word is ',word
         if (WORD.EQ.'END') then
            write (*,*) 'hind rotors must be defined in tsr1.dat'
            write (21,*) 'hind rotors must be defined in tsr1.dat'
            stop
         endif
      enddo
      read (11,*) nhind
c      write(*,*)'nhind is ', nhind
      nhind1=nhind
      read (11,*)
      do ihind = 1 , nhind
         read (11,*) hindlab(ihind),dhindmn(ihind),dhindmx(ihind),
     +              nhindsteps(ihind),nhindsymm(ihind)
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) hindlab(ihind)
         rewind (99)
         call LineRead (99)
         hindlab(ihind)=WORD
         close (unit=99,status='keep')
      enddo
      rewind(11)

      word=''
      do while (WORD.NE.'NHIND')
         call LineRead (12)
         if (WORD.EQ.'END') then
            write (*,*) 'hind rotors must be defined in tsr2.dat'
            write (21,*) 'hind rotors must be defined in tsr2.dat'
            stop
         endif
      enddo
      read (12,*) nhind
      nhindtot=nhind1+nhind
      ijump=0
      iprog=nhind1+1
      read (12,*)
      do ihind = nhind1+1, nhindtot
         read (12,*) hindlabread,dhindmnread,dhindmxread,
     +              nhindstepsread,nhindsymmread
            write(*,*)'hindlabread ',hindlabread
            write(*,*)'ijump ',ijump
         if(hindlabread.eq.'BABS2'.and.ilin1.eq.1.and.iabs.eq.1)ijump=1
         if(hindlabread.eq.'babs2'.and.ilin1.eq.1.and.iabs.eq.1)ijump=1
         if(ijump.eq.0)then
            write(*,*)'hindlabread ',hindlabread
            hindlab(iprog)=hindlabread
            dhindmn(iprog)=dhindmnread
            dhindmx(iprog)=dhindmxread
            nhindsteps(iprog)=nhindstepsread
            nhindsymm(iprog)=nhindsymmread
c         read (12,*) hindlab(ihind),dhindmn(ihind),dhindmx(ihind),
c     +              nhindsteps(ihind),nhindsymm(ihind)
            open (unit=99,status='unknown')
            rewind (99)
            write (99,*) hindlab(iprog)
            rewind (99)
            call LineRead (99)
            hindlab(iprog)=WORD
            close (unit=99,status='keep')
            iprog=iprog+1
         endif
         if(ijump.eq.1)nhindtot=nhindtot-1
         ijump=0
      enddo
      rewind(12)

cc now write dihedral input for stochastic scan

      nsmp=nhindtot*nhindtot+3

      if(nhindtot.eq.0)nsmp=1
      
      write(10,*)' nosmp dthresh,ethresh '
      write(10,*)nsmp,' 1.0 0.00001'
      write(10,*)
      write(10,*)'ntau'
      write(10,*)nhindtot
      write(10,*)'-->scanned dihedral coordinates, min and max values'
      do j=1,nhindtot
         write(10,100)hindlab(j),dhindmn(j),dhindmx(j)
      enddo
      write(10,*)

 100  format(A8,F5.1,3X,F5.1)
cc now write dihedral input for HR scan

      write(10,*)'nhind'
      write(10,*)nhindtot
      write(10,*)'-->namehind,hindmn,hindmx,hindstep periodicity'
      do j=1,nhindtot
         write(10,101)hindlab(j),dhindmn(j),dhindmx(j),
     +              nhindsteps(j),nhindsymm(j)
      enddo
      write(10,*)
 101  format(A8,F5.1,3X,F5.1,3x,I3,3x,I2)

cc for wellr and wellp write default parameters

      write(16,*)' nosmp dthresh,ethresh '
      write(16,*)'1 1.0 0.00001'
      write(16,*)
      write(16,*)'ntau'
      write(16,*)'0'
      write(16,*)'-->scanned dihedral coordinates, min and max values'
      write(16,*)
      write(16,*)'nhind'
      write(16,*)'0'
      write(16,*)'-->namehind,hindmn,hindmx,hindstep periodicity'
      write(16,*)

      write(17,*)' nosmp dthresh,ethresh '
      write(17,*)'1 1.0 0.00001'
      write(17,*)
      write(17,*)'ntau'
      write(17,*)'0'
      write(17,*)'-->scanned dihedral coordinates, min and max values'
      write(17,*)
      write(17,*)'nhind'
      write(17,*)'0'
      write(17,*)'-->namehind,hindmn,hindmx,hindstep periodicity'
      write(17,*)

cc now read and write abstraction/addition site

      do while (WORD.NE.'ISITE')
         call LineRead (11)
         if (WORD.EQ.'END') then
            write (*,*) 'abs/add site must be defined in tsr1.dat'
            write (21,*) 'abs/add site must be defined in tsr1.dat'
            stop
         endif
      enddo
      call LineRead (11)
      write(10,*)'isite ji ki'
      write(10,*)line
      write(10,*)
      rewind(11)

cc check if reaction is exhothermic or endothermic
      
      prod2_en=0.
      prod2_zpe=0.

cc now read and write grid parameters

cc for abstraction reactions check if it is exothermic
      if(iabs.eq.1)then
         open (unit=22,file='../me_files/reac1_en.me',status='UNKNOWN')
         read(22,*)reac1_en
         close(22)
         open (unit=22,file='../me_files/reac1_zpe.me',status='UNKNOWN')
         read(22,*)reac1_zpe
         close(22)
         open (unit=22,file='../me_files/reac2_en.me',status='UNKNOWN')
         read(22,*)reac2_en
         close(22)
         open (unit=22,file='../me_files/reac2_zpe.me',status='UNKNOWN')
         read(22,*)reac2_zpe
         close(22)
         open (unit=22,file='../me_files/prod1_en.me',status='UNKNOWN')
         read(22,*)prod1_en
         close(22)
         open (unit=22,file='../me_files/prod1_zpe.me',status='UNKNOWN')
         read(22,*)prod1_zpe
         close(22)
         open (unit=22,file='../me_files/prod2_en.me',status='UNKNOWN')
         read(22,*)prod2_en
         close(22)
         open (unit=22,file='../me_files/prod2_zpe.me',status='UNKNOWN')
         read(22,*)prod2_zpe
         close(22)
         reac_dhr=prod2_en+prod2_zpe+prod1_en+prod1_zpe-reac2_en
     $         -reac2_zpe
     $         -reac1_en-reac1_zpe
         reac_dhr=reac_dhr*627.5

         write(*,*)'the enthalpy change is ',reac_dhr,' kcal/mol'
         write(21,*)'the enthalpy change is ',reac_dhr,' kcal/mol'
      endif

      iback=0
      if(iabs.eq.1.and.reac_dhr.gt.5.0)iback=1

      do while (WORD.NE.'RMIN')
         call LineRead (12)
         if (WORD.EQ.'END') then
            write (*,*) 'grid must be defined in tsr2.dat'
            write (21,*) 'grid must be defined in tsr2.dat'
            stop
         endif
      enddo
      if(iback.eq.0)then
         write(10,*)'rmin rmax nr '
      else
         write(10,*)'rmin reverse '
      endif
      call LineRead (12)
      write(10,*)line
      call LineRead (12)
      write(10,*)line
      call LineRead (12)
      write(10,*)line

c      read(12,*)label
c      write(10,*)label
c      read(12,*)label
c      write(10,*)label
c      read(12,*)label
c      write(10,*)label
      write(10,*)
      rewind(12)

cc now read charge and spin from reac1 and reac2

      do while (WORD.NE.'CHARGE')
         call LineRead (13)
         if (WORD.EQ.'END') then
            write (*,*) 'charge must be defined in reac1.dat'
            write (21,*) 'charge must be defined in reac1.dat'
            stop
         endif
      enddo
      read(13,*)nch1,nsp1
      rewind(13)
      word=''

      do while (WORD.NE.'CHARGE')
         call LineRead (14)
         if (WORD.EQ.'END') then
            write (*,*) 'charge must be defined in reac2.dat'
            write (21,*) 'charge must be defined in reac2.dat'
            stop
         endif
      enddo
      read(14,*)nch2,nsp2
      rewind(14)

      nch=nch1+nch2
      nsp=nsp1+nsp2-1
      write(10,*)'charge spin'
      write(10,*)nch,nsp
      write(10,*)
cc do the same for wellr and wellp
      write(16,*)'charge spin'
      write(16,*)nch,nsp
      write(16,*)

      write(17,*)'charge spin'
      write(17,*)nch,nsp
      write(17,*)

cc write symmetry, not decided at this step 

      write(10,*)'SymmetryFactor'
      write(10,*)'1'
      write(10,*)
cc do the same for wellr and wellp
      write(16,*)'SymmetryFactor'
      write(16,*)'1'
      write(16,*)

      write(17,*)'SymmetryFactor'
      write(17,*)'1'
      write(17,*)

cc decide electronic states, assume only ground state for now

      write(10,*)'nelec '
      write(10,*)'1 '
      write(10,*)'0. ',nsp

      write(10,*)
      write(10,*)'End'
      write(10,*)

cc do the same for wellr and wellp

      write(16,*)'nelec '
      write(16,*)'1 '
      write(16,*)'0. ',nsp
      write(16,*)
      write(16,*)'End'
      write(16,*)

      write(17,*)'nelec '
      write(17,*)'1 '
      write(17,*)'0. ',nsp
      write(17,*)
      write(17,*)'End'
      write(17,*)


      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(21)
      stop

      end

      SUBROUTINE LineRead(IUnit)

C     Subroutine reads the next non-blank non-comment line from IUNIT.
C     WORD contains the first keyword (uppercase).
C     WORD2 contains the second keyword, if any.
C     WORD3 contains the third keyword, if any.
C     WORD4 contains the fourth keyword, if any.
C     WORD5 contains the fifth keyword, if any.
C     NWORD is 1 if one keyword, 2 if two keywords, 3 if three keywords.
C     If the first non-blank character is *, LSEC is set to
C     .TRUE. and WORD is the section name.
C     If the first non-blank character is &, the line is printed. 
C     Comments begin with ! and may begin anywhere in a line.
C     Comments following the keywords are detected and ignored.
C     If end of file is detected LEOF is set to .TRUE.

      implicit double precision(a-h,o-z)
      implicit integer (i-n)


      logical leof,lsec,ltit
      character*100 line,sename,string,word,word2,word3,title,title1
     #,word4,word5,word6,word7

      include 'filcomm.f'
cc

      if(iunit.eq.0) then
        word=''
        word2=''
        word3=''
        word4=''
        word5=''
        string=''
        return
      endif

 5    Continue
      LSEC = .FALSE.
      LEOF = .FALSE.
      LTIT=.FALSE.

 10   READ(iunit,'(A100)',END=9000) STRING
c     write (6,*) 'string test',string

      IP = 0
      NWORD=1

C     Find the first nonblank character in the line

 20   IP = IP + 1
      IF (((STRING(IP:IP).EQ.' ').OR.(STRING(IP:IP).EQ.'        '))
     &     .AND.(IP .Lt. 100)) GO TO 20           

C     If it  is all blank or is a comment line then read the next line

      IF (IP .ge. 100 .OR. STRING(IP:IP) .EQ. '!')  GO TO 10

C     If it is a title line write it and read next line.

      IF (STRING(IP:IP) .EQ. '&') THEN
         IPP=IP+1 
         TITLE1=STRING(IPP:)
         Write (6,27) TITLE1
 27      Format (a)
         Go to 5
      ENDIF

C     Convert to uppercase

      call upcase

C     Find the position of a comment character, if any,
C     in a nonblank noncomment line

      ICOM = IP 
 50   ICOM = ICOM + 1
      IF (LINE(ICOM:ICOM) .NE. '!' .AND. ICOM .LT. 100)  GOTO 50      

 55   IF (IP .LT. ICOM) THEN
         IBEG = IP
 60      IP = IP + 1
         IF (LINE(IP:IP) .NE. ' ' .AND. IP .LT. ICOM) GOTO 60 
         WORD = LINE(IBEG:IP-1)
      ENDIF
c     write (6,*) 'word test',word,ip

c     Find second keyword

 70   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 70
      NWORD=2
      IBEG2=IP

 80   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 80 
      WORD2 = LINE(IBEG2:IP-1)
c     write (6,*) 'word2 test',word2,ip

C     Find third keyword

 75   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 75
      NWORD=3
      IBEG3=IP

 85   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 85 
      WORD3 = LINE(IBEG3:IP-1)
c     write (6,*) 'word3 test',word3,ip

C     Find fourth keyword

175   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 175
      NWORD=4
      IBEG4=IP

185   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 185 
      WORD4 = LINE(IBEG4:IP-1)
      IF (WORD4 .EQ. '=') GOTO 175
c     write (6,*) 'word4 test',word4,ip

C     Find fifth keyword

275   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 275
      NWORD=5
      IBEG5=IP

285   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 285 
      WORD5 = LINE(IBEG5:IP-1)
      IF (WORD5 .EQ. '=') GOTO 275
c     write (6,*) 'word5 test',word5,ip

C     See if WORD is a section header

 90   IF (WORD(1:1) .EQ. '*') THEN
         LSEC=.TRUE.
         SENAME = WORD(2:)
         WORD = SENAME
         If (WORD.eq.'END') Then
            LEOF = .TRUE.
            Return
         Endif
         LSEC = .TRUE.
      ENDIF

      RETURN

 9000 Continue
      Write (6,*) 'Unexpected end of input file upon calling LineRead'
      Write (6,*) 'Last keyword read was ',WORD
      go to 9900

 9900 continue

      Stop
      end

      subroutine upcase

C     Function which takes a string of 80 characters and converts the 
C     lower case letters in the string to upper case letters
C     This function is a modified version of CASE which was written
C     by Rozeanne Steckler

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      CHARACTER * 1 XLETT

c     keyword common block

      logical leof,lsec,ltit
      character*100 line,sename,string,word,word2,word3,title,title1,
     #word4,word5,word6,word7
      common /keyword/ leof,lsec,ltit
      common /key/ line,sename,string,word,word2,word3,word4,word5,
     $ word6,word7,title,title1
      common /nkey/ nword

      LINE = STRING
      DO 10 I = 1, 100
         XLETT = LINE(I:I)
         ITRY = ICHAR (XLETT)
         IF (XLETT .GE. 'a' .AND. XLETT .LE. 'z') THEN 
            ITRY = ITRY - 32
            LINE(I:I) = CHAR (ITRY)
         ENDIF
 10   CONTINUE

      RETURN
      END


