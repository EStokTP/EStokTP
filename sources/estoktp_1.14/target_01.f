C      Copyright (C) 2019 Carlo de Falco
C      
C      This program is free software: you can redistribute it and/or modify it
C      under the terms of the GNU General Public License as published by
C      the Free Software Foundation, either version 3 of the License, or
C      (at your option) any later version.
C      
C      This program is distributed in the hope that it will be useful, but
C      WITHOUT ANY WARRANTY; without even the implied warranty of
C      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C      GNU General Public License for more details.
C      
C      You should have received a copy of the GNU General Public License
C      along with this program.  If not, see <https://www.gnu.org/licenses/>.

      
c      block data init
c      implicit none
c      integer n
c      parameter (n = 2)
c      double precision  qinit(n)
c      common /general/ qinit
c      data qinit / 0.0D0, 0.0D0/
c      end
      
C     *****************************
      subroutine problem_size(nres)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      integer n, nres
      character*30 cjunk

      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,n,nconst
      close(15)
      nres = n-nconst
c      write(*,*)'nvar is',nres
c      stop
      
      end subroutine

C     *****************************
      subroutine initial_guess(q)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx),qinit(nmdmx),ibond(nmdmx)
      character*30 cjunk

      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,np,nconst

      do j=1,natomt
         read(15,*)cjunk
      enddo
      do j=1,np
         read(15,*)cjunk,qinit(j),ibond(j)
      enddo
      close(15)
     
      do j=1,np
         q(j) = qinit(j)
         if(ibond(j).eq.1)then
            q(j)=q(j)/CAUTOANG
         else
            q(j)=q(j)/180.*pigr
         endif
      enddo
      
      end subroutine

C     *****************************
      subroutine V1(q, En1)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx)
      dimension ibond(nmdmx)
      dimension freq(nmdmx),tau(ntaumx),tauopt(ntaumx)
      dimension coord(natommx,ndim),xint(3*natommx),
     $ abcrot(ndim) 

      LOGICAL leof,lsec,ltit
      
      character*70 comline1,comline2
      character*70 comline3,comline4
      character*30 intcoor(3*natommx)
      character*60 atomlabel(natommx)
      character*20 bislab(ntaumx)
      character*30 cjunk
      character*30 gmem
      character*30 namepes1,namepes2
      character*100 commandcopy

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

c      integer n
c      parameter (n = 2)

      En1=0
      iguess=0
c      write(*,*)'En 1 is', En1
c      write(*,*)'q 1 is', q(1)
c      stop
      write(*,*)'entering V1 '

      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,np,nconst

      do j=1,natomt
         read(15,'(A70)')atomlabel(j)
      enddo
      do j=1,np
         read(15,*)intcoor(j),xint(j),ibond(j)
      enddo

      do j=1,np-nconst
         if(ibond(j).eq.1)then
            xint(j)=q(j)*CAUTOANG
         else
            xint(j)=q(j)*180./pigr
         endif
      enddo

      read(15,*)cjunk
      read(15,*)ilevcode
      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         read(15,'(A70)')comline1
         read(15,'(A70)')comline2
         read(15,'(A70)')comline3
         read(15,'(A70)')comline4
      else if (ilevcode.eq.2)then
         read(15,'(A30)')namepes1
         read(15,'(A30)')namepes2
         read(15,*)cjunk
         read(15,*)cjunk
         write(commandcopy,100)namepes1
 100     format ('cp -f ',A30,' natst_molpro.dat')
c         write(7,*)commandcopy
         call commrun(commandcopy)
c         write(7,*)'ok 1'
c         close(7)
c         stop
      endif
      read(15,*)cjunk
      read(15,*)icharge,ispin,iguess
      close(15)

c      write(7,*)'ok 1'
c      stop
c      write(*,*)atomlabel(5)
c      stop
      ismp=0
      ifreq=0
      ircons=0
      ixyz=0
      ired=0
      ntau=0
      ilin=0
      ires=0
      ispecies=0
      iaspace=0
      ilev=0

      if(ilevcode.eq.1) then
         if(iguess.eq.1)then
            commandcopy='cp -f tmp_na1.chk tmp.chk'
            call commrun(commandcopy)
         endif
         call g09fopt(ilevcode,tau,ntau,natom,natomt,numproc,gmem,
     $        coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $        comline2,icharge,ispin,ircons,
     $        atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $        ires,ixyz,ired)

      else if (ilevcode.eq.2) then
         ilev=3
         if(iguess.eq.1.and.ilevcode.eq.3)then
            commandcopy='cp -f tmp_na1.chk tmp.chk'
            call commrun(commandcopy)
         endif
         call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot_0,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies,iaspace)
      else
         call elstructopt(ilevcode,tau,ntau,natom,natomt,numproc,
     $        gmem,
     $        coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $        comline2,icharge,ispin,ircons,
     $        atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $        ires,ixyz,ired,ispecies,iaspace)

      endif
   
      open (unit=15,file='geom_pes1.xyz',status='unknown')
      write(15,*)natomt
      write(15,*)'point on PES1'
      do j=1,natomt
         write(15,*)j,(coord(j,idim),idim=1,3)
      enddo
      close(15)

c      stop
      En1 = vtot_0

      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         commandcopy='cp -f geom.log natst_V1.log'
      else if (ilevcode.eq.2)then
         commandcopy='cat molpro.out molpro.log > geom.log'
         call commrun(commandcopy)
         commandcopy='echo endfile >> geom.log'
         call commrun(commandcopy)
         commandcopy='cp -f geom.log natst_V1.log'
      endif

      call commrun(commandcopy)

c      commandcopy='cp -f grad_xyz.dat grad1_xyz.dat'
c      call commrun(commandcopy)

      write(*,*)'En1 is ',En1
c      write(7,*)'En1 is ',En1
c      stop
      write(*,*)'out of V1 '

      end subroutine

C     *****************************
      subroutine V2(q, En2)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx)
      dimension freq(nmdmx),tau(ntaumx),tauopt(ntaumx)
      dimension coord(natommx,ndim),xint(3*natommx),
     $ abcrot(ndim) 
      dimension ibond(nmdmx)

      LOGICAL leof,lsec,ltit
      
      character*70 comline1,comline2
      character*70 comline3,comline4
      character*30 intcoor(3*natommx)
      character*60 atomlabel(natommx)
      character*20 bislab(ntaumx)
      character*30 cjunk
      character*30 gmem
      character*30 namepes1,namepes2
      character*100 commandcopy

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

c      integer n
c      parameter (n = 2)
      write(*,*)'entering V2 '

      En2=0
c      write(*,*)'En 1 is', En1
c      write(*,*)'q 1 is', q(1)
c      stop

      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,np,nconst

      do j=1,natomt
         read(15,'(A70)')atomlabel(j)
      enddo
      do j=1,np
         read(15,*)intcoor(j),xint(j),ibond(j)
      enddo

      do j=1,np-nconst
         if(ibond(j).eq.1)then
            xint(j)=q(j)*CAUTOANG
         else
            xint(j)=q(j)*180./pigr
         endif
      enddo

      read(15,*)cjunk
      read(15,*)ilevcode

      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         read(15,'(A70)')comline1
         read(15,'(A70)')comline2
         read(15,'(A70)')comline3
         read(15,'(A70)')comline4
      else if (ilevcode.eq.2)then
         read(15,'(A30)')namepes1
         read(15,'(A30)')namepes2
         read(15,*)cjunk
         read(15,*)cjunk
         write(commandcopy,100)namepes2
 100     format ('cp -f ',A30,' natst_molpro.dat')
c         write(7,*)commandcopy
         call commrun(commandcopy)
      endif
      read(15,*)cjunk
      read(15,*)cjunk,cjunk
      read(15,*)icharge,ispin,iguess
      close(15)

c      write(*,*)atomlabel(5)
c      stop
      ismp=0
      ifreq=0
      ircons=0
      ixyz=0
      ired=0
      ntau=0
      ilin=0
c      ilevcode=1
      ired=0
      ires=0
      ilev=0

      if(ilevcode.eq.1) then
         if(iguess.eq.1)then
            commandcopy='cp -f tmp_na2.chk tmp.chk'
            call commrun(commandcopy)
         endif

         call g09fopt(ilevcode,tau,ntau,natom,natomt,numproc,gmem,
     $        coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline3,
     $        comline4,icharge,ispin,ircons,
     $        atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $        ires,ixyz,ired)

      else if (ilevcode.eq.2) then
         ilev=3
         call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot_0,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies,iaspace)
      else
         if(iguess.eq.1.and.ilevcode.eq.3)then
            commandcopy='cp -f tmp_na2.chk tmp.chk'
            call commrun(commandcopy)
         endif
         call elstructopt(ilevcode,tau,ntau,natom,natomt,numproc,gmem,
     $        coord,vtot_0,vtot,freq,ifreq,ilin,ismp,comline1,
     $        comline2,icharge,ispin,ircons,
     $        atomlabel,intcoor,bislab,tauopt,xint,abcrot,
     $        ires,ixyz,ired,ispecies,iaspace)

      endif
   
      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         commandcopy='cp -f geom.log natst_V2.log'
      else if (ilevcode.eq.2)then
         do j=1,np-nconst
            if(ibond(j).eq.1)then
               xint(j)=q(j)*CAUTOANG
            else
               xint(j)=q(j)*180./pigr
            endif
         enddo
         commandcopy='cat molpro.out molpro.log > geom.log'
         call commrun(commandcopy)
         commandcopy='echo endfile >> geom.log'
         call commrun(commandcopy)
         commandcopy='cp -f geom.log natst_V2.log'
      endif

      call commrun(commandcopy)

      En2 = vtot_0
      write(*,*)'En2 is ',En2
c
cc this call overwrites geom.log
c

      call readV1(En1)
      write(7,*)'En1: ',En1,'En2: ',En2, 'En1-En2: ',En1-En2 
      write(31,*)'En1: ',En1,'En2: ',En2, 'En1-En2: ',En1-En2 

c      stop

c      commandcopy='cp -f grad.dat grad2.dat'
c      call commrun(commandcopy)
c      commandcopy='cp -f grad_xyz.dat grad2_xyz.dat'
c      call commrun(commandcopy)

      open (unit=15,file='na_mecp.out',status='unknown')
      write(15,*)'MECP structure '
c      write(15,*)natomt,np
      do j=1,natomt
         write(15,*)atomlabel(j)
      enddo
      do j=1,np
         write(15,*)intcoor(j),xint(j)
      enddo
      close(15)

c
c      integer n
c      parameter (n = 2)
c      double precision q(n), En2
c      En2 = q(1)**2 + q(2)**2 - (q(1) + q(2) -1.0D0)

      end subroutine


C     *****************************
      subroutine dV1(q, dz)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx)
      dimension dz(nmdmx)
      dimension grad(3*natommx)
      dimension ibond(nmdmx)

      character*30 cjunk
      character*100 commandcopy

c      implicit none
c      integer n
c      double precision q(n), dz(n)

      write(*,*)'entering dv1 '
      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,np,nconst

      do j=1,natomt
         read(15,*)cjunk
      enddo
      do j=1,np
         read(15,*)cjunk,cjunk,ibond(j)
      enddo

      read(15,*)cjunk
      read(15,*)ilevcode
      close(15)

c      write(*,*)'ilevcode is ',ilevcode
      commandcopy='cp -f natst_V1.log geom.log'
      call commrun(commandcopy)

      call problem_size(nint)
      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         call readgrad_g09(nint,grad)
      else if (ilevcode.eq.2) then
         call readgrad_molpro(nint,grad)
         do j=1,nint
            if(ibond(j).eq.1)then
               grad(j)=grad(j)
            else
               grad(j)=grad(j)*180./pigr
            endif
         enddo
      endif

      do j=1,nint
         dz(j)=grad(j)
      enddo

c      En1=0.
c      call V1(q,En1)
c      write(*,*)'last grad value is ',dz(nint)
c      stop

c      dz(1) = 2.0D0*q(1)
c      dz(2) = 2.0D0*q(2)
      write(*,*)'out of dv1'

      end subroutine

C     *****************************
      subroutine dV2(q, dz)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx)
      dimension dz(nmdmx)
      dimension grad(3*natommx)
      dimension ibond(nmdmx)

      character*30 cjunk
      character*100 commandcopy

c      implicit none
c      integer n
c      double precision q(n), dz(n)

c      write(*,*)'entering in grad'
c      En2=0.
c      call V2(q,En2)

      write(*,*)'entering dv2 '

      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,np,nconst

      do j=1,natomt
         read(15,*)cjunk
      enddo
      do j=1,np
         read(15,*)cjunk,cjunk,ibond(j)
      enddo

      read(15,*)cjunk
      read(15,*)ilevcode
      close(15)

      commandcopy='cp -f natst_V2.log geom.log'
      call commrun(commandcopy)

      call problem_size(nint)

      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         call readgrad_g09(nint,grad)
      else if (ilevcode.eq.2) then
         call readgrad_molpro(nint,grad)
         do j=1,nint
            if(ibond(j).eq.1)then
               grad(j)=grad(j)
            else
               grad(j)=grad(j)*180./pigr
            endif
         enddo
      endif

c      call readgrad_g09(nint,grad)

      do j=1,nint
         dz(j)=grad(j)
      enddo
c      write(*,*)'last grad value is ',dz(nint)
c      stop
c      dz(1) = 2.0D0*q(1)
c      dz(2) = 2.0D0*q(2)
      write(*,*)'out of dv2'

      end subroutine

C     *****************************
      subroutine readV1(En1)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx)
      dimension dz(nmdmx)
      dimension grad(3*natommx)

      character*30 cjunk
      character*30 gkeyword,igkey
      character*70 comline1,comline2
      character*100 commandcopy

c      implicit none
c      integer n
c      double precision q(n), dz(n)

c      call problem_size(nint)
c      write(*,*)'entering in grad'
c      write(*,*)'nvar is ',nint

      write(*,*)'entering readV1'

      commandcopy='cp -f natst_V1.log geom.log'
      call commrun(commandcopy)

      open (unit=15,file='na_input.dat',status='unknown')
      read(15,*)cjunk
      read(15,*)natom,natomt,np,nconst

      do j=1,natomt
         read(15,*)cjunk
      enddo
      do j=1,np
         read(15,*)cjunk
      enddo

      read(15,*)cjunk
      read(15,*)ilevcode
      read(15,'(A70)')comline1
      read(15,'(A70)')comline2
      read(15,*)cjunk
      read(15,*)cjunk
      read(15,*)cjunk
      read(15,*)icharge,ispin,iguess
      close(15)

c      write(*,*)'comline ',comline1
c      write(*,*)'spin ',ispin

      En1=0.
      if(ilevcode.eq.1.or.ilevcode.eq.3)then
         call en_key_g09(comline1,igkey,gkeyword,ispin)
         call readen_g09(En1,gkeyword,igkey)
      else if (ilevcode.eq.2)then
         call readen_molpro(En1)
      endif

      write(*,*)'energy is ',En1
c      stop
      write(*,*)'out of readV1'

      end subroutine

C     *****************************
C     *****************************
      subroutine calcxyzgrad_molpro(inapes,natom,ianum,grad_xyz)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension q(nmdmx)
      dimension ibond(nmdmx)
      dimension freq(nmdmx),tau(ntaumx),tauopt(ntaumx)
      dimension coord(natommx,ndim),xint(3*natommx),
     $ abcrot(ndim) 
      dimension ianum(3*natommx)
      dimension grad_xyz(3*natommx)

      LOGICAL leof,lsec,ltit
      
      character*30 gmem
      character*70 comline1,comline2
      character*70 comline3,comline4
      character*30 intcoor(3*natommx)
      character*60 atomlabel(natommx)
      character*20 bislab(ntaumx)
      character*30 cjunk
      character*30 namepes1,namepes2
      character*100 commandcopy,comm1

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      include 'filcomm.f'

c      integer n
c      parameter (n = 2)

c      stop
      write(*,*)'entering calcxyzgrad_molpro'

      if(inapes.eq.1.or.inapes.eq.2)then
         open (unit=15,file='na_input.dat',status='unknown')
         read(15,*)cjunk
         read(15,*)natom,natomt,np,nconst
         do j=1,natomt
            read(15,*)cjunk
         enddo
         do j=1,np
            read(15,*)cjunk
         enddo
         read(15,*)cjunk
         read(15,*)ilevcode
         read(15,'(A30)')namepes1
         read(15,'(A30)')namepes2
         read(15,*)cjunk
         read(15,*)cjunk
         read(15,*)cjunk
         if(inapes.eq.1)then
            read(15,*)icharge,ispin,iguess
         else if (inapes.eq.2)then
            read(15,*)cjunk
            read(15,*)icharge,ispin,iguess
         endif
         close(15)
         if(inapes.eq.1)then
            write(commandcopy,100)namepes1
         else if (inapes.eq.2)then
            write(commandcopy,100)namepes2
         endif
 100     format ('cp -f ',A30,' natst_molpro.dat')
c         write(7,*)commandcopy
         call commrun(commandcopy)
         comm1="sed -ie 's/optg,maxit=1/optg,varsave,maxit=1;coord,3N'/ 
     $ natst_molpro.dat"
         call commrun(comm1)
      endif

      open (unit=15,file='./output/na_mecp.out',status='unknown')
      read(15,*)cjunk
      do j=1,natomt
         read(15,'(A70)')atomlabel(j)
      enddo
      do j=1,np
         read(15,*)intcoor(j),xint(j)
      enddo
      close(15)


c      write(7,*)'ok 1'
c      stop
c      write(*,*)atomlabel(5)
c      stop
      ismp=0
      ifreq=0
      ircons=0
      ixyz=0
      ired=0
      ntau=0
      ilin=0
      ires=0
      ispecies=0
      iaspace=0
      ilev=4

      commandcopy='cp -f ./geoms/natst_mecp.molden molpro.molden'
      call commrun(commandcopy)

      call molprofopt(tau,ntau,natom,natomt,numproc,gmem,
     $ coord,vtot_0,freq,ifreq,ilin,ismp,icharge,ispin,ircons,
     $ atomlabel,intcoor,bislab,tauopt,xint,abcrot,ires,ixyz,ilev
     $ ,ispecies,iaspace)

      commandcopy='cat molpro.out molpro.log > geom.log'
      call commrun(commandcopy)
      commandcopy='echo endfile >> geom.log'
      call commrun(commandcopy)
      comm1='sed -i "s/\///g" geom.log'
      call commrun(comm1)

      commandcopy='echo endfile >> molpro.molden'
      call commrun(commandcopy)
   
      if(inapes.eq.1)then
         commandcopy='cp -f geom.log ./na_tst/natst_V1_xyz.log'
      else if (inapes.eq.2)then
         commandcopy='cp -f geom.log ./na_tst/natst_V2_xyz.log'
      endif
      call commrun(commandcopy)

      open (unit=11,status='old',file='molpro.molden')
      rewind (11)
114   continue
      call LineRead (0)
      call LineRead (11)

cc read atom numbers

      if((WORD.EQ.'[ATOMS]').AND.(WORD2.EQ.'ANGS'))then
c         gradval=0.
         iprog=0
         do ip = 1 , natom
            read(11,*)cjunk,cjunk,ianum(iprog+1)
            write(*,*)cjunk,cjunk,ianum(iprog+1)
            ianum(iprog+2)=ianum(iprog+1)
            ianum(iprog+3)=ianum(iprog+1)
            iprog=iprog+3
         enddo
      endif
c      write(*,*)'word is ',word
      if (WORD.eq.'ENDFILE')goto 9000
      goto 114
      
 9000 continue
      close(11)

      open (unit=11,status='old',file='geom.log')
c      rewind (11)
115   continue
      call LineRead (0)
      call LineRead (11)

cc read gradient 

      if((WORD.EQ.'OPTIMIZATION').AND.(WORD2.EQ.'POINT'))then
         gradval=0.
         read(11,*)cjunk
         read(11,*)cjunk
         nvar=natom*3

         do iread = 1 , nvar
            read(11,*)cjunk,cjunk,cjunk,cjunk,cjunk,gradval
            grad_xyz(iread)=gradval
            write(*,*)'grad is ',grad_xyz(iread)
         enddo
      endif
      if (WORD.eq.'ENDFILE')goto 9001
      goto 115
      
 9001 continue
      close(11)

c      commandcopy='cp -f grad_xyz.dat grad1_xyz.dat'
c      call commrun(commandcopy)
c      stop

      return
      end 
ccccccccccccccccc
