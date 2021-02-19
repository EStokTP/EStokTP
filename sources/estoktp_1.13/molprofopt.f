c molprofopt a program to start and read output from molpro 2010-2015
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.

      SUBROUTINE molprofopt(tau,ntau,natom,natomt,nshared,gmemo,
     $ coord,vtot,freq,ifreq,ilin,ismp,icharge, ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies,iaspace)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c     parameter (ntaumx=10, nmdmx=300, natommx=100, ndim=3)

      integer natom,iatom,idim
      dimension freq(nmdmx),tau(ntaumx),tauopt(ntaumx)
      dimension coord(natommx,ndim),xint(3*natommx),abcrot(ndim)
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension gradx(natommx),grady(natommx),gradz(natommx)

      logical leof,lsec,ltit

      character*160 stringr
      character*100 comline1
      character*100 command1
      character*30 gmemo
      character*30 cjunk
      character*60 atomlabel(natommx)
      character*80 irconslabel
      character*60 label
      character*30 intcoor(3*natommx)
      character*30 cootoread
      character*20 bislab(ntaumx)
      character*2 aname
      character*2 atomtype(natommx)
      character*2 atomtype_new
      
c     character*160 cname

 
      character*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

c      character*160 sename,word,word2,word3
c     $ ,title,title1,word4,word5,word6,word7

      character*20 ct1,ct2,ct3

c      COMMON /keyword/ leof,lsec,ltit
c      COMMON /key/ line,sename,string,word,word2,word3
c     $ ,word4,word5,title,title1

c      common /tsdat/ aabs1,babs1,aabs2,babs2,babs3,rts
c      common /itsdat/ iabs,iadd,ivar,isite,ji,ki
c      common /strucword/ stoich
      include 'filcomm.f'

111   continue

cccccccccccccc print out molpro input file  ccccccccccccc
c initialize word, word2, word3, word4, word5
      call LineRead (0)
      ichecken=0
      irepeat=0
c      ncoord = natom*3-6-ntau-ircons
      ncoord = natom*3-6-ntau
      if (natom.eq.2) ncoord = 1
c      write(*,*)'natomt is',natomt
c      stop
c      write(*,*)'ircons is ',ircons

 999  continue

      open(unit=10,status='unknown',file='molpro.inp')
      if(iaspace.eq.1)then
         open(unit=98,status='unknown',file='molpro_step2.inp')
      endif
      if(ilev.eq.0) then
         open(unit=11,status='unknown',file='level0_molpro.dat')
      else if(ilev.eq.1) then
         open(unit=11,status='unknown',file='level1_molpro.dat')
      else if (ilev.eq.2) then
         open(unit=11,status='unknown',file='hl_molpro.dat')
      else if (ilev.eq.3.or.ilev.eq.4) then
         open(unit=11,status='unknown',file='natst_molpro.dat')
      else if (ilev.eq.10) then
         open(unit=11,status='unknown',file='onedtau_molpro.dat')
      else
         write(*,*)'molpro call not implemented yet '
         write(*,*)'for this type of calculation'
         stop
      endif

cc assign memory for calculation
cc it is assumed the memory in the input is given as MW

      open(unit=99,status='unknown')

      igmem=100
      open(unit=99,file='temp.tmp',status='unknown')
      write(99,*)gmemo
      close(99)
      command1="sed -ie 's/MW/ /' temp.tmp"
      call commrun(command1)
      open(unit=99,file='temp.tmp',status='unknown')
      read(99,*)igmem
      close(99)
      write(10,*)'memory,',igmem,',m'

 100  continue
      read (11,'(A100)') comline1
      if (comline1.EQ.'End1'.or.comline1.eq.' End1') go to 200
      write (10,*) comline1
      if(iaspace.eq.1)write(98,*) comline1
      goto 100
 200  continue
      write (10,*)'geometry={angstrom'
      if(iaspace.eq.1)write(98,*)'geometry={angstrom'
      if(ixyz.ne.1)then
         do iatom = 1 , natomt
            write (10,*)atomlabel(iatom)
            if(iaspace.eq.1)write (98,*)atomlabel(iatom)
         enddo
      else
         do iatom = 1 , natom
            open(unit=99,status='unknown')
            write(99,*)atomlabel(iatom)
            rewind(99)
            read(99,*)aname,tcoox,tcooy,tcooz
            write(10,*)aname,', ',tcoox,', ',tcooy,', ',tcooz
            if(iaspace.eq.1)then
               write(98,*)aname,', ',tcoox,', ',tcooy,', ',tcooz
            endif
            close(99)
         enddo
      endif
      write (10,*)'}'
      if(iaspace.eq.1) write (98,*)'}'
      
      if(ixyz.ne.1)then
         do icoord = 1 , ncoord-1
            write (10,*) intcoor(icoord),' = ',xint(icoord)
            if(iaspace.eq.1)then
               write (98,*) intcoor(icoord),' = ',xint(icoord)
            endif
         enddo
         if(ilev.eq.10.and.iaspace.eq.1)then
            open(unit=99,file='hrcondcoord.dat',status='unknown')
            read(99,'(A100)')comline1
            close(99)
            write(10,*)comline1
            if(ispin.eq.1)then
               open(unit=15,file='activespace.dat',status='unknown')
               read(15,*)nlines
               do j=1,nlines
                  read(15,'(A100)')comline1
               enddo
               write (10,*) comline1
               write(10,*)'gthresh,step=5.d-3'
               close(15)
            endif
            write(10,*)
         endif
         write (10,*) intcoor(ncoord),' = ',xint(ncoord)
         if(iaspace.eq.1)then
            write (98,*) intcoor(ncoord),' = ',xint(ncoord)
         endif
         do itau = 1 , ntau
            write (10,*) bislab(itau),' = ',tau(itau)
            if(iaspace.eq.1)then
               write (98,*) bislab(itau),' = ',tau(itau)
            endif
         enddo
      endif
      write (10,*)
      write (10,*) 'SET,SPIN=',ispin-1
      if(iaspace.eq.1)write (98,*)
      if(iaspace.eq.1)write (98,*)'SET,SPIN=',ispin-1

 101  continue
      word2=''
      read (11,'(A100)') comline1
      write(99,*)'line ',comline1
      rewind(99)
      call lineread(99)
      rewind(99)
c      write(7,*)word2
c      write(7,*)comline1
      if(iaspace.eq.1.and.word2.eq.'ACTIVE_SPACE1')then
         if(ispin.eq.1)then
            open(unit=15,file='activespace.dat',status='unknown')
            read(15,*)nlines
            do j=1,nlines
               read(15,'(A100)')comline1
               write (10,*) comline1
               write (98,*) comline1
            enddo
            close(15)
            write(98,*)'End1'
            write(98,*)
            write (98,*) comline1
         endif
         goto 101
      endif
      if(iaspace.eq.1.and.word2.eq.'ACTIVE_SPACE2')then
         if(ispin.eq.1)then
            open(unit=15,file='activespace.dat',status='unknown')
            read(15,*)nlines
            do j=1,nlines
               read(15,'(A90)')comline1
            enddo
            write (10,*) comline1
            write (98,*) comline1
            close(15)
         endif
         goto 101
      endif
      if (comline1.EQ.'End2'.or.comline1.eq.' End2') go to 201
      if (ispin.eq.1) write (10,*) comline1
      if (ispin.eq.1)then 
         write (98,*) comline1
      endif
      goto 101
 201  continue
      close(99)
      if (ispin.eq.1)then 
         write (98,*) 'End2'
      endif
c      close(7)
c      stop
      if (ispin.eq.1) go to 401
cc else read input for open shell
 301  continue

c      read (11,'(A70)') comline1
      word2=''
      read (11,'(A100)') comline1
      write(99,*)'line ',comline1
      rewind(99)
      call lineread(99)
      rewind(99)
c      write(7,*)word2
c      write(7,*)comline1
      if(iaspace.eq.1.and.word2.eq.'ACTIVE_SPACE1')then
c         if(ispin.eq.1)then
            open(unit=15,file='activespace.dat',status='unknown')
            read(15,*)nlines
            do j=1,nlines
               read(15,'(A100)')comline1
               write (10,*) comline1
               write (98,*) comline1
            enddo
            close(15)
            write(98,*)'End1'
            write(98,*)'End2'
            write(98,*)
            write (98,*) comline1
c         endif
         goto 301
      endif
      if(iaspace.eq.1.and.word2.eq.'ACTIVE_SPACE2')then
c        if(ispin.eq.1)then
            open(unit=15,file='activespace.dat',status='unknown')
            read(15,*)nlines
            do j=1,nlines
               read(15,'(A100)')comline1
            enddo
            write (10,*) comline1
            write (98,*) comline1
            close(15)
c         endif
         goto 301
      endif

      if (ispin.gt.1)then 
         write (98,*) comline1
      endif
      if (comline1.EQ.'End3'.or.comline1.eq.' End3') go to 401
      write (10,*) comline1
      goto 301
 401  continue
 
      close(11)
      write(10,*)'put,molden,molpro.molden'
      write(10,*)'---'
      if (ispin.eq.1)then 
         write(98,*)'put,molden,molpro.molden'
         write(98,*)'---'
         write (98,*) 'End3'
      endif
      close(98)
      close (unit=10,status='keep')

      ntotcoord = natom*3-6
      if(ircons.ne.0)then
         open(unit=99,file='temp.tmp',status='unknown')
         if (ircons.eq.1) then
            write (99,*) 'inactive,',intcoor(ntotcoord)
         else if (ircons.eq.2) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1)
         else if (ircons.eq.3) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1),
     $           ',',intcoor(ntotcoord-2)
         else if (ircons.eq.4) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1),
     $           ',',intcoor(ntotcoord-2),
     $           ',',intcoor(ntotcoord-3)
         else if (ircons.eq.5) then
            write (99,*)'inactive,',intcoor(ntotcoord),
     $           ',',intcoor(ntotcoord-1),
     $           ',',intcoor(ntotcoord-2),
     $           ',',intcoor(ntotcoord-3),
     $           ',',intcoor(ntotcoord-4)
         endif
         close(99)
         command1="sed -ie 's/ //g' temp.tmp"
         call commrun(command1)

         open(unit=99,file='temp.tmp',status='unknown')
         read(99,'(A80)')irconslabel
         rewind(99)
         write(99,1200)irconslabel
         rewind(99)
         read(99,'(A80)')irconslabel
         close(99)
         command1=irconslabel
         call commrun(command1)

 1200    format (" sed -ie 's/inactive/"A35"/g' molpro.inp")
         write(*,*)'irconslab is ',irconslabel
      endif
c
      if(ispecies.eq.0.and.ircons.eq.3.and.ilev.eq.10)then
         open(unit=99,status='unknown')
         write(99,1202)
         rewind(99)
         read(99,'(A60)')label
         close(99)
         command1=label
         call commrun(command1)
 1202    format (" sed -ie 's/root=2/root=1/g' molpro.inp")
         open(unit=99,status='unknown')
         write(99,1203)
         rewind(99)
         read(99,'(A60)')label
         close(99)
         command1=label
         call commrun(command1)
 1203    format (" sed -ie 's/Root=2/root=1/g' molpro.inp")
      endif

c            write (*,*)'intcoord is ',intcoor(1)
c            stop

ccccccccccccc run molpro  cccccccccccccccc
      command1='rm -f molpro.xml'       
      call commrun(command1)
      command1='rm -f molpro.log_*'       
      call commrun(command1)

      if(irecov.ne.1)then
         command1='rm -f molpro.out'       
         call commrun(command1)
         call molprorun(numproc)
      endif
      call check_molpro_exit(ilev)

cc first copy output file to geom.log, for saving

      command1='cp -f molpro.out geom.log'       
      call commrun(command1)

cc now read results

      call LineRead (0)

cc first get energy

      command1='egrep CBSEN molpro.out > temp1.log'
      call commrun(command1)
      command1="egrep 'CBSEN\(1' molpro.out >> temp1.log"
      call commrun(command1)
      command1='tail -n 1 temp1.log > temp.log'
      call commrun(command1)

      open (unit=99,status='unknown',file='temp.log')
      read(99,*)cjunk,cjunk,cjunk,vtot
      close(99)
      write(*,*)'vtot is',vtot

cc then get coordinates if not an xyz calculations
      
c      open (unit=100,status='unknown',file='molpro.molden')
      if(ixyz.ne.1)then
         do j=1,ncoord+ntau
c         rewind(100)
            open (unit=99,status='unknown',file='temp.log')
            if(j.le.ncoord) write(99,*) intcoor(j)
            if(j.gt.ncoord) write(99,*) bislab(j-ncoord)
            close(99)
            call trimtemp
            open (unit=99,status='unknown',file='temp1.log')
            call lineread(99)
            cootoread=word
            write(*,*)'cootoread is ',cootoread
            rewind(99)
            write(99,1201)cootoread
 1201       format (" egrep "A20" molpro.molden > temp.log")
            rewind(99)
            read(99,'(A60)')command1
            call commrun(command1)
            close(99)         
            command1="sed -ie 's/=/ /g' temp.log"
            call commrun(command1)
            open (unit=99,status='unknown',file='temp.log')
            read(99,*)cjunk,xint(j)         
            close(99)
c         stop
c         do while (WORD.ne.cootoread)
c            CALL LineRead (100)
c            write(*,*)'word is ',word
c         enddo
c         open (unit=99,status='unknown')
c         write(99,*) word2
c         rewind(99)
c         read(99,*)xint(j)
c         close(99)
         enddo
c      close(100)
         iread=1
         itau=0
         do j=1,ncoord+ntau
            write(*,*)'opt coord is',xint(j)
            itau = iread-(natom*3-6-ntau)
            if (iread.gt.ncoord) then
               tauopt(itau)=xint(j)
               if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $              tauopt(itau) - 360.0d0
               if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $              tauopt(itau) + 360.0d0
            endif
            iread=iread+1
         enddo
      endif
c      stop
cc then get and write xyz coordinates

      open (unit=100,status='unknown',file='molpro.molden')
      open (unit=101,status='unknown',file='geom.xyz')

      do while (WORD.ne.'[ATOMS]')
         CALL LineRead (100)
      enddo
      nread=natomt
      if (ixyz.eq.1) nread=natom
      do j=1,nread
         read(100,*)atomtype(j),djunk,djunk,coox(j),cooy(j),cooz(j)
         if(atomtype(j).ne.'X')then
            call update_atomtype(atomtype(j),atomtype_new)
            write(101,1204)atomtype_new,coox(j),cooy(j),cooz(j)
            coord(j,1)=coox(j)
            coord(j,2)=cooy(j)
            coord(j,3)=cooz(j)
         endif
      enddo
      write(101,*)
      close(100)
      close(101)

 1204    format(1X,A2,1X,F9.5,1X,F9.5,1X,F9.5)

cc now read frequencies, if computed

cc first do a check in the molden file to see if the frequencies key was there

      command1='grep -q "\[FREQ" molpro.molden && echo 1 > temp.log || 
     $  echo 0 > temp.log'
      call commrun(command1)

      open (unit=99,status='unknown',file='temp.log')
      read(99,*)ifreq
      close(99)

cc then read fcmat.log file, even if no freqs must be written

      open (unit=101,status='unknown',file='fcmat.log')
      write(101,*)'frequencies not found in molden file'
      close(101)

cc then read frequencies

      if(ifreq.eq.1)then
         open (unit=100,status='unknown',file='molpro.molden')

         do while (WORD.ne.'[FREQ]')
            CALL LineRead (100)
         enddo
         index=0.
         nread=3*natom
         do j=1,nread
            read(100,*)freqread
c            write(*,*)'freq ',j,' is ',freqread
            if(freqread.gt.0.5)then
               index=index+1
               freq(index)=freqread
cc now assign negative value to imaginary frequency
               if(freq(index).lt.freq(index-1))then
                  freq(index-1)=-freq(index-1)
c                  index=index-1
               endif
            endif
         enddo
         nfreq=3*natom-6
         if(ilin.eq.1)  nfreq=3*natom-5
         if(natom.eq.2)  nfreq=1
c         if(ispecies.eq.0)nfreq=nfreq-1
         if(nfreq.ne.index)then
            write(*,*)' there is disagreement between '
            write(*,*)' expected and read frequencies'
            write(*,*)' the program will be stopped'
            write(*,*)' read frequencies ',index
            write(*,*)' expected frequencies ',nfreq
            stop
         endif
         do j=1,nfreq
            write(*,*)'freq ',j,' is ',freq(j)
         enddo
         close(100)

cc now proceed to saving the hessian
         open (unit=100,status='unknown',file='molpro.out')
         open (unit=101,status='unknown',file='fcmat.log')
         write(101,*)'Force constants in Cartesian coordinates'
         numlines=0
         nummax=natom*3/5
         do j=0,nummax
            numlines=numlines+natom*3-5*j+1
         enddo
c         numlines=numlines-1
         write(*,*)'numlines =',numlines

         do while (stringr.ne.
     $  ' Force Constants (Second Derivatives of the Energy) in [a.u.]')
            read(100,'(A160)')stringr
            if(stringr.eq.' Variable memory released')then
               write(*,*)'did not found Hessian infos in the output'
               goto 500
            endif
         enddo
         do j=1,numlines
            read(100,'(A160)')stringr
            write(101,'(A100)')stringr
         enddo
         write(101,*)
         write(101,*)
 500     continue
         write(101,*)'Input orientation'
         write(101,*)'---'
         write(101,*)'Center'
         write(101,*)'Number'
         write(101,*)'---'

         do j=1,natom
            if(atomtype(j).eq.'H')ianumb=1
            if(atomtype(j).eq.'B')ianumb=5
            if(atomtype(j).eq.'C')ianumb=6
            if(atomtype(j).eq.'N')ianumb=7
            if(atomtype(j).eq.'O')ianumb=8
            if(atomtype(j).eq.'F')ianumb=9
            if(atomtype(j).eq.'Ne')ianumb=10
            if(atomtype(j).eq.'Al')ianumb=13
            if(atomtype(j).eq.'Si')ianumb=14
            if(atomtype(j).eq.'P')ianumb=15
            if(atomtype(j).eq.'S')ianumb=16
            if(atomtype(j).eq.'Cl')ianumb=17
            if(atomtype(j).eq.'CL')ianumb=17
            if(atomtype(j).eq.'Ar')ianumb=18
            if(atomtype(j).eq.'Ga')ianumb=31
            if(atomtype(j).eq.'Ge')ianumb=32
            if(atomtype(j).eq.'As')ianumb=33
            if(atomtype(j).eq.'Se')ianumb=34
            if(atomtype(j).eq.'Br')ianumb=35
            if(atomtype(j).eq.'Kr')ianumb=36
            if(atomtype(j).eq.'Rh')ianumb=45
            if(atomtype(j).eq.'In')ianumb=49
            if(atomtype(j).eq.'Sn')ianumb=50
            if(atomtype(j).eq.'Sb')ianumb=51
            if(atomtype(j).eq.'Te')ianumb=52
            if(atomtype(j).eq.'I')ianumb=53
            if(atomtype(j).eq.'Xe')ianumb=54
            write(101,1002)j,ianumb,0,coox(j),cooy(j),cooz(j)
         enddo

cc now save the gradient
         open(unit=99,status='unknown') 
         write(99,900)natom+3
         rewind(99)
         read(99,'(A100)')command1
         call commrun(command1)

         rewind(99)
         write(99,901)natom+1
         rewind(99)
         read(99,'(A100)')command1
         call commrun(command1)
         close(99)

         open(unit=100,file='end.log',status='unknown')
         write(100,*)'End'
         close(100)
         
         command1='cat grad2.txt end.log > grad.txt'
         call commrun(command1)

         command1='rm -f grad2.txt'
         call commrun(command1)
         command1='rm -f end.log'
         call commrun(command1)

cc count lines

         icount=0
         open(unit=100,file='grad.txt',status='unknown')
         do while (WORD.NE.'END')
            call LineRead(100)
            icount=icount+1
         enddo
         write(*,*)'icount= ',icount
         if(icount.eq.natom+1)then
            rewind(100)
            do j=1,natom
               read(100,*)cjunk,gradx(j),grady(j),gradz(j)
            enddo
         else
            write(101,*)'unable to read gradient'
            write(101,*)'gradient set to 0.'
            do j=1,natom
               gradx(j)=0.
               grady(j)=0.
               gradz(j)=0.
            enddo
         endif

 900     format("egrep -A"I0.2 " 'GRADIENT FOR' molpro.log > 
     +       grad1.txt")
 901     format("tail -n"I0.2 " grad1.txt > grad2.txt")

         write(101,*)
         write(101,*)
         write(101,*)'Forces '
         write(101,*)'---'
         write(101,*)'Center'

         do j=1,natom
            if(atomtype(j).eq.'H')ianumb=1
            if(atomtype(j).eq.'B')ianumb=5
            if(atomtype(j).eq.'C')ianumb=6
            if(atomtype(j).eq.'N')ianumb=7
            if(atomtype(j).eq.'O')ianumb=8
            if(atomtype(j).eq.'F')ianumb=9
            if(atomtype(j).eq.'Ne')ianumb=10
            if(atomtype(j).eq.'Al')ianumb=13
            if(atomtype(j).eq.'Si')ianumb=14
            if(atomtype(j).eq.'P')ianumb=15
            if(atomtype(j).eq.'S')ianumb=16
            if(atomtype(j).eq.'Cl')ianumb=17
            if(atomtype(j).eq.'CL')ianumb=17
            if(atomtype(j).eq.'Ar')ianumb=18
            if(atomtype(j).eq.'Ga')ianumb=31
            if(atomtype(j).eq.'Ge')ianumb=32
            if(atomtype(j).eq.'As')ianumb=33
            if(atomtype(j).eq.'Se')ianumb=34
            if(atomtype(j).eq.'Br')ianumb=35
            if(atomtype(j).eq.'Kr')ianumb=36
            if(atomtype(j).eq.'Rh')ianumb=45
            if(atomtype(j).eq.'In')ianumb=49
            if(atomtype(j).eq.'Sn')ianumb=50
            if(atomtype(j).eq.'Sb')ianumb=51
            if(atomtype(j).eq.'Te')ianumb=52
            if(atomtype(j).eq.'I')ianumb=53
            if(atomtype(j).eq.'Xe')ianumb=54
            write(101,1001)j,ianumb,gradx(j),grady(j),gradz(j)
         enddo
         write(101,*)
         write(101,*)
         close(100)
         close(101)
      endif

1001  format(1X,I2,1X,I2,1X,3(1x,F11.6))
1002  format(1X,I2,1X,I2,1X,I2,1X,3(1x,F9.4))

      write(6,*)'out of molprofopt' 

      RETURN
      END 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_molpro_exit(ilev)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      character*100 command1
      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7
  
      include 'filcomm.f'

      command1='tail -n 1 molpro.out > temp.log'
      call commrun(command1)

      open(unit=99,file='temp.log',status='unknown')
      call LineRead (99)
      close(99)
      if(word2.eq.'ERROR')then
        open(unit=99,file='failed',status='unknown')
        write(99,*)'molpro job failed'
        if(ilev.eq.0) then
           write(99,*)'doing level 0 calculations'        
        else if(ilev.eq.1) then
           write(99,*)'doing level 1 calculations'                
        else if (ilev.eq.2) then
           write(99,*)'doing high level calculations'        
        else if (ilev.eq.10) then
           write(99,*)'doing hind rot calculations'        
        else
           write(99,*)'at undefined level '
        endif
        close(99)
        stop
      endif

      command1="egrep 'PROGRAM \* OPT|No conv' molpro.out > temp.log"
      call commrun(command1)

      command1='tail -n 1 molpro.out > temp.log'
      call commrun(command1)

      open(unit=99,file='temp.log',status='unknown')
      call LineRead (99)
      close(99)
      if(word.eq.'NO'.and.word2.eq.'CONVERGENCE')then
        open(unit=99,file='failed',status='unknown')
        write(99,*)'molpro job failed'
        if(ilev.eq.0) then
           write(99,*)'doing level 0 calculations'        
        else if(ilev.eq.1) then
           write(99,*)'doing level 1 calculations'                
        else if (ilev.eq.2) then
           write(99,*)'doing high level calculations'        
        else if (ilev.eq.10) then
           write(99,*)'doing hind rot calculations'        
        else
           write(99,*)'at undefined level '
        endif
        write(99,*)'in the optimization stage'
        close(99)
        stop
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_atomtype(atomtype_in,atomtype_out)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      character*100 command1
      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7
      character*2 atomtype_in,atomtype_out

      include 'filcomm.f'

      atomtype_out=atomtype_in
      if(atomtype_in.eq.'CL')atomtype_out='Cl'
      if(atomtype_in.eq.'NE')atomtype_out='Ne'
      if(atomtype_in.eq.'AL')atomtype_out='Al'
      if(atomtype_in.eq.'SI')atomtype_out='Si'
      if(atomtype_in.eq.'AR')atomtype_out='Ar'
      if(atomtype_in.eq.'GA')atomtype_out='Ga'
      if(atomtype_in.eq.'GE')atomtype_out='Ge'
      if(atomtype_in.eq.'AS')atomtype_out='As'
      if(atomtype_in.eq.'SE')atomtype_out='Se'
      if(atomtype_in.eq.'BR')atomtype_out='Br'
      if(atomtype_in.eq.'KR')atomtype_out='Kr'
      if(atomtype_in.eq.'RH')atomtype_out='Rh'
      if(atomtype_in.eq.'SN')atomtype_out='Sn'
      if(atomtype_in.eq.'SB')atomtype_out='Sb'
      if(atomtype_in.eq.'TE')atomtype_out='Te'
      if(atomtype_in.eq.'XE')atomtype_out='Xe'

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readgrad_molpro(nvar,grad)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension grad(3*natommx)

      character*30 gkeyword,igkey
      character*70 cjunk,cjunk1,cjunk2
      character*70 comline1
      character*70 comline2
      character*70 comline3
      character*70 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      command1='sed -i "s/\///g" geom.log'
      call commrun(command1)

      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)

cc read gradient in internal coordinates

      if((WORD.EQ.'OPTIMIZATION').AND.(WORD2.EQ.'POINT'))then
         gradval=0.
         read(11,*)cjunk
         read(11,*)cjunk

         do iread = 1 , nvar
            read(11,*)cjunk,cjunk,cjunk,cjunk,cjunk,gradval
            grad(iread)=gradval
c            grad(iread)=0.1
c            write(*,*)iread,gradval
         enddo
      endif
c      write(*,*)'word is ',word
      if (WORD.eq.'ENDFILE')goto 9000
      goto 114
      
 9000 continue
      close(11)
c      stop

      return
      end

C     *****************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readen_molpro(En1)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      character*30 gkeyword,igkey
      character*70 cjunk,cjunk1,cjunk2
      character*70 comline1
      character*70 comline2
      character*70 comline3
      character*70 comline4
      character*120 command1
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

      include 'filcomm.f'

      En1=0.
      OPEN (unit=11,status='old',file='geom.log')
      rewind (11)
114   CONTINUE
      CALL LineRead (0)
      CALL LineRead (11)

cc read gradient in internal coordinates

      if((WORD.EQ.'SETTING').AND.(WORD2.EQ.'CBSEN'))then
c         gradval=0.

         read(word4,100)En1
 100     format (F14.7)

         read(11,*)cjunk
         read(11,*)cjunk

      endif
c      write(*,*)'word is ',word
      if (WORD.eq.'ENDFILE')goto 9000
      goto 114
      
 9000 continue
      close(11)
c      stop

      return
      end

C     *****************************
