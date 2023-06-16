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
      character*1000 line,string
      character*160 sename,word,word2,word3,title,title1
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

 10   READ(iunit,'(A1000)',END=9000) STRING
c     write (6,*) 'string test',string

      IP = 0
      NWORD=1

C     Find the first nonblank character in the line

 20   IP = IP + 1
      IF (((STRING(IP:IP).EQ.' ').OR.(STRING(IP:IP).EQ.'	'))
     &     .AND.(IP .Lt. 1000)) GO TO 20           

C     If it  is all blank or is a comment line then read the next line

      IF (IP .ge. 1000 .OR. STRING(IP:IP) .EQ. '!')  GO TO 10

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
      IF (LINE(ICOM:ICOM) .NE. '!' .AND. ICOM .LT. 1000)  GOTO 50      

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

c     **********************************************************************

      SUBROUTINE LineRead2(IUnit)

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

c     keyword common block

      logical leof,lsec,ltit
      character*1000 line,string
      character*160 sename,word,word2,word3,title,title1
     #,word4,word5,word6,word7

      include 'filcomm.f'

 5    Continue
      LSEC = .FALSE.
      LEOF = .FALSE.
      LTIT=.FALSE.

 10   READ(iunit,'(A1000)',END=9000) STRING
c     write (6,*) 'string2 test',string

      IP = 0
      NWORD=1

C     Find the first nonblank character in the line

 20   IP = IP + 1
      IF (((STRING(IP:IP).EQ.' ').OR.(STRING(IP:IP).EQ.'        '))
     &     .AND.(IP .Lt. 1000)) GO TO 20           

C     If it  is all blank or is a comment line then read the next line

      IF (IP .ge. 1000)  GO TO 10

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
      IF (ICOM .LT. 1000)  GOTO 50      

 55   IF (IP .LT. ICOM) THEN
         IBEG = IP
 60      IP = IP + 1
         IF (LINE(IP:IP) .NE. ' ' .AND. IP .LT. ICOM) GOTO 60 
         WORD = LINE(IBEG:IP-1)
      ENDIF

c     Find second keyword

 70   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 70
      NWORD=2
      IBEG2=IP

 80   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 80 
      WORD2 = LINE(IBEG2:IP-1)

C     Find third keyword

 75   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 75
      NWORD=3
      IBEG3=IP

 85   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 85 
      WORD3 = LINE(IBEG3:IP-1)

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


c     **********************************************************************

      SUBROUTINE LineRead3(IUnit)

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

c     keyword common block

      logical leof,lsec,ltit
      character*1000 line,string
      character*160 sename,word,word2,word3
      character*160 title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      if(iunit.eq.0) then
        word=''
        word2=''
        word3=''
        word4=''
        word5=''
        word6=''
        word7=''

        string=''
        return
      endif


 5    Continue
      LSEC = .FALSE.
      LEOF = .FALSE.
      LTIT=.FALSE.

 10   READ(iunit,'(A1000)',END=9000) STRING
c     write (6,*) 'string2 test',string

      IP = 0
      NWORD=1

C     Find the first nonblank character in the line

 20   IP = IP + 1
      IF (((STRING(IP:IP).EQ.' ').OR.(STRING(IP:IP).EQ.'        '))
     &     .AND.(IP .Lt. 1000)) GO TO 20           

C     If it  is all blank or is a comment line then read the next line

      IF (IP .ge. 1000)  GO TO 10

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
      IF (ICOM .LT. 1000)  GOTO 50      

 55   IF (IP .LT. ICOM) THEN
         IBEG = IP
 60      IP = IP + 1
         IF (LINE(IP:IP) .NE. ' ' .AND. IP .LT. ICOM) GOTO 60 
         WORD = LINE(IBEG:IP-1)
      ENDIF

c     Find second keyword

 70   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 70
      NWORD=2
      IBEG2=IP

 80   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 80 
      WORD2 = LINE(IBEG2:IP-1)

C     Find third keyword

 75   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 75
      NWORD=3
      IBEG3=IP

 85   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 85 
      WORD3 = LINE(IBEG3:IP-1)

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

C     Find sixth keyword

375   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 375
      NWORD=6
      IBEG6=IP

385   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 385 
      WORD6 = LINE(IBEG6:IP-1)
      IF (WORD6 .EQ. '=') GOTO 375

C     Find seventh keyword

475   IP=IP+1
      IF (IP.EQ.ICOM) GO TO 90
      IF (LINE(IP:IP).EQ.' ') GO TO 475
      NWORD=7
      IBEG7=IP

485   IP = IP + 1
      IF (LINE(IP:IP).NE. ' ' .AND.IP.LT. ICOM) GOTO 485 
      WORD7 = LINE(IBEG7:IP-1)
      IF (WORD7 .EQ. '=') GOTO 475


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



c     **********************************************************************
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
      character*1000 line,string
      character*160 sename,word,word2,word3,title,title1,
     #word4,word5,word6,word7
      common /keyword/ leof,lsec,ltit
      common /key/ line,sename,string,word,word2,word3,word4,word5,
     $ word6,word7,title,title1
      common /nkey/ nword

      LINE = STRING
      DO 10 I = 1, 1000
         XLETT = LINE(I:I)
         ITRY = ICHAR (XLETT)
         IF (XLETT .GE. 'a' .AND. XLETT .LE. 'z') THEN 
            ITRY = ITRY - 32
            LINE(I:I) = CHAR (ITRY)
         ENDIF
 10   CONTINUE

      RETURN
      END


c     **********************************************************************
      subroutine upcase2(convline)

C     Function which takes a string of 80 characters and converts the 
C     lower case letters in the string to upper case letters
C     This function is a modified version of CASE which was written
C     by Rozeanne Steckler

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      CHARACTER * 1 XLETT

c     keyword common block

      logical leof,lsec,ltit
      character*1000 line,string
      character*1000 convline
      character*160 sename,word,word2,word3,title,title1,
     #word4,word5,word6,word7
      common /keyword/ leof,lsec,ltit
      common /key/ line,sename,string,word,word2,word3,word4,word5,
     $ word6,word7,title,title1
      common /nkey/ nword

      LINE = CONVLINE
      DO 10 I = 1, 1000
         XLETT = LINE(I:I)
         ITRY = ICHAR (XLETT)
         IF (XLETT .GE. 'a' .AND. XLETT .LE. 'z') THEN 
            ITRY = ITRY - 32
            LINE(I:I) = CHAR (ITRY)
         ENDIF
 10   CONTINUE
      convline=line

      RETURN
      END



