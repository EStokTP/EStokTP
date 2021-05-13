c a set of programs to read and convert z-matrices
c Copyright (C) 2018  by Carlo Cavallotti and Matteo Pelucchi, Politecnico di MIlano, Italy
c and Stephen J. Klippenstein, Argonne National Laboratories, USA.
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, version 3 of the License
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

      subroutine read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
 
c      integer ibconn
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension idummy(natommx)
      character*60 atomlabel(natommx)
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)
      character*1 xread
      character*5 ccheck
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7
      character*20 ct1,ct2,ct3
      include 'filcomm.f'

c      write(*,*)'natomt is ', natomt

      do j=1,natomt
         atname(j)=''
         bconnt(j)=''
         bname(j)=''
         aconnt(j)=''
         anname(j)=''
         dconnt(j)=''
         dname(j)=''
         ibconn(j)=0
         iaconn(j)=0
         idconn(j)=0
      enddo
c      write(*,*)'ok 1'
c      stop
c
      do j=1,natomt
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) atomlabel(j)
         rewind (99)
         call LineRead3 (99)
         close(99)
         atname(j)=word
         bconnt(j)=word2
         bname(j)=word3
         if(j.ge.2)then
            aconnt(j)=word4
            anname(j)=word5
         endif
         if(j.ge.3)then
            dconnt(j)=word6
            dname(j)=word7
         endif
      enddo

      do k=1,natomt
c         do j=1,natommx
         do j=1,19
            if(j.lt.10)write(ccheck,"(I1)")j
            if(j.lt.100.and.j.gt.9)write(ccheck,"(I2)")j
            if(j.lt.1000.and.j.gt.99)write(ccheck,"(I3)")j
c            write(*,*)'ccheck is ',ccheck
c            write(*,*)'bconnt is ',bconnt(k)
            if(bconnt(k).eq.ccheck)then
               bconnt(k)=atname(j)
            endif
            if(aconnt(k).eq.ccheck)then
               aconnt(k)=atname(j)
            endif
            if(dconnt(k).eq.ccheck)then
               dconnt(k)=atname(j)
            endif
         enddo
c         stop
      enddo
cc now determine connectivity in terms of progressive atom numbering

      do j=1,natomt
         do i=1,natomt
            if(bconnt(j).eq.atname(i))then
               ibconn(j)=i
            endif
            if(aconnt(j).eq.atname(i))then
               iaconn(j)=i
            endif
            if(dconnt(j).eq.atname(i))then
               idconn(j)=i
            endif
         enddo
      enddo

c      do j=1,natomt
c         write(7,*)'at ',j,' conn is: ',ibconn(j),iaconn(j),idconn(j)
c      enddo

c      do j=1,natomt
c         write(7,*)'at ',j,'bond and dihed are: ',bname(j),anname(j)
c     +        ,dname(j)
c      enddo
c      return

cc now check for dummy atoms
      ndummy=0
      do j=1,natomt
         open (unit=99,status='unknown')
         rewind (99)
         write (99,*) atname(j)
         rewind (99)
         read(99,*)xread
         close(99)
         if(xread.eq.'X') then
            idummy(j)=1
         else
            idummy(j)=0
         endif
         write(*,*)'idummy is j: ',j, idummy(j)
      enddo
      isited=0
      jsited=0
      ksited=0
      do j=1,isite
         if(idummy(j).eq.0)isited=isited+1
      enddo
      do j=1,jsite
         if(idummy(j).eq.0)jsited=jsited+1
      enddo
      do j=1,ksite
         if(idummy(j).eq.0)ksited=ksited+1
      enddo

c      write(*,*)'isite isited ',isite,isited
c      write(*,*)'jsite jsited ',jsite,jsited
c      write(*,*)'ksite ksited ',ksite,ksited

      write(*,*)'out of read_zmat'

      return
      end

c ******************************************************

      subroutine update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $ ,idconn,bname,anname,dname,atname,coox,cooy,cooz,xint,tauopt
     $ ,ntau,idummy,ilin,aconnt,bconnt,dconnt,atomlabel,ifilu)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
 
c     
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension xint(3*natommx)
      dimension xintt(3*natommx)
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension cooxs(natommx),cooys(natommx),coozs(natommx)
      dimension cooxt(natommx),cooyt(natommx),coozt(natommx)
      dimension tauopt(ntaumx)
      dimension idummy(natommx)

      character*1 aname1(natommx)
      character*30 intcoor(3*natommx)
      character*20 bislab(ntaumx)
      character*20 cname
      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)
      LOGICAL leof,lsec,ltit
      character*60 atomlabel(natommx)

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7


c      include 'filcomm.f'

c      write(*,*)'ifile is ',1
c      stop

      nint = natom*3-6-ntau
      if (natom.eq.2) nint=1
      
      do j=1,nint
         xintt(j)=xint(j)
         xint(j)=0.
      enddo
 
      do j=1,ntau
         tauopt(j)=0.
      enddo

      do j=1,natom
         cooxs(j)=coox(j)
         cooys(j)=cooy(j)
         coozs(j)=cooz(j)
      enddo

c      open (unit=98,status='unknown')
c      do j=1,natom
c         write(98,*)atname(j),coox(j),cooy(j),cooz(j)
c      enddo
c      close(98)
c      stop


c rototraslate the xyz matrix
c      write(*,*)'ilin = ',ilin
c      do j=1,natom
c         write(*,*)'in config ',coox(j),cooy(j),cooz(j)
c      enddo


      call rototrasl(natom,coox,cooy,cooz,ilin)

c      do j=1,natom
c         write(*,*)'out config ',coox(j),cooy(j),cooz(j)
c      enddo

cc
cc write xyz matrix   

c      open (unit=97,status='unknown')
c      do j=1,natom
c         write(97,*)atname(j),coox(j),cooy(j),cooz(j)
c      enddo
c      write(97,*)'ilin is',ilin
c      close(97)
c      stop

cc first introduce fake coordinates for dummy atoms
cc to be consistent with numering in input z-mat

      ind=0
      iangname=0
      do j=1,natomt
c         write(*,*)'idummy is ', j,idummy(j)
         if(idummy(j).eq.0) then
            ind=ind+1
            cooxt(j)=coox(ind)
            cooyt(j)=cooy(ind)
            coozt(j)=cooz(ind)
         else if (j.eq.3.and.idummy(j).eq.2) then
cc it is assumed a dummy atom is on the z axis with y=0 perp to atom2 or 1 
cc this has been deactivated
            read(bname(j),'(f10.0)')bd
            read(anname(j),'(f10.0)')ang

c            write (99,*) bname(j)
c            write (99,*) anname(j)
c            rewind (99)
c            read(99,*)bd
c            read(99,*)ang
c            close(99)
            if(ibconn(j).eq.1)then
               cooxt(j)=0.
               cooyt(j)=0.
               coozt(j)=bd
            else
               cooxt(j)=cooxt(2)
               cooyt(j)=0.
               coozt(j)=bd
            endif
c            stop
         else if (j.eq.3.and.idummy(j).eq.1) then
cc it is assumed a dummy atom is on the z axis with y=0 perp to atom2 or 1 
c            open (unit=99,status='unknown')
            read(bname(j),'(f10.0)')bd
            read(anname(j),'(f10.0)')ang
c            open (unit=99,status='unknown')
c            rewind (99)
c            write (99,*) bname(j)
c            write (99,*) anname(j)
c            rewind (99)
c            read(99,*)bd
c            read(99,*)ang
c            close(99)
c            if(iabs.eq.1.and.ireact.gt.3) write (99,*) dname(j)
            iangindex=0
            do ik=4,natomt
               if(idummy(ik).ne.1)then
                  iangname2=0
                  do ij=1,nint
                     if(intcoor(ij).eq.dname(ik))then
                        iangname2=1
                     endif
                  enddo
                  if (iangname2.eq.0)then
                     iangindex=ik
                  endif
               endif
            enddo
c
            if (iangindex.ne.0) then
               adist2=bd**2
               iat1=ibconn(j)
               iat2=iaconn(j)
               da=bd
               db=sqrt((coox(iat1)-coox(iat2))**2+
     $                    (cooy(iat1)-cooy(iat2))**2+
     $                    (cooz(iat1)-cooz(iat2))**2)
               alfa=ang
               dc=0.
               call distang(da,db,dc,alfa)
               bdist2=dc**2
               iat3=iangindex
               itestind=0
               if(bconnt(iat3).eq.atname(j))itestind=1
               if(aconnt(iat3).eq.atname(j))itestind=1
               if(itestind.eq.1)then
                  write(*,*)'found reference atom for dummy atom3'
               endif

               iatemp1=iangindex
               icorr=0
               do jk=1,iatemp1-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iatemp1=iatemp1-icorr
               iai=iatemp1
               iatemp2=ibconn(iangindex)
               icorr=0
               do jk=1,iatemp2-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iatemp2=iatemp2-icorr
               iab=iatemp2
               iatemp3=idconn(iangindex)
               icorr=0
               do jk=1,iatemp3-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iatemp3=iatemp3-icorr
               iac=iatemp3

               da=sqrt((coox(iai)-coox(iab))**2+
     $                    (cooy(iai)-cooy(iab))**2+
     $                    (cooz(iai)-cooz(iab))**2)

               db=sqrt((coox(iac)-coox(iatemp2))**2+
     $                    (cooy(iac)-cooy(iatemp2))**2+
     $                    (cooz(iac)-cooz(iatemp2))**2)

c              open (unit=99,status='unknown')
c               rewind (99)
c               write (99,*) dname(iangindex)
c               rewind (99)
c               read(99,*)dih_ia
c               rewind (99)
c               write (99,*) anname(j)
c               rewind (99)
c               read(99,*) ang_abc
c               close(99)

               read(dname(iangindex),'(f10.0)')dih_ia
               read(anname(j),'(f10.0)')ang_abc

c               write(*,*)'dihed is',dih_ia
c               write(*,*)'ang_abc is',ang_abc
               if(dih_ia.eq.180.)then
                  write(*,*)'recognized opposite planar configuration'
               else
c                  write(*,*)'failed z-mat conv. in routine update_zmat' 
               endif
               dac=sqrt((coox(iai)-coox(iac))**2+
     $                    (cooy(iai)-cooy(iac))**2+
     $                    (cooz(iai)-cooz(iac))**2)

               alfa2=0.
               if(abs(dac-da-db).gt.1.0e-6)then
                  call distang2(da,db,dac,alfa2)
               else
                  alfa2=180.
               endif
               alfa=alfa2-ang_abc
               
               db=bd
               call distang(da,db,dc,alfa)
               cdist2=dc**2
               icorr=0
               do jk=1,iat3-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iat3=iat3-icorr

               call twoan_to_xyz(adist2,bdist2,cdist2,coox(iat1)
     $  ,cooy(iat1),cooz(iat1),coox(iat2),cooy(iat2),cooz(iat2),
     $  coox(iat3),cooy(iat3),cooz(iat3),xa,ya,za)

               cooxt(j)=xa
               cooyt(j)=ya
               coozt(j)=za
            else
cc it is assumed a dummy atom is on the z axis with y=0 perp to atom2 or 1 
c               open (unit=99,status='unknown')
c               rewind (99)
c               write (99,*) bname(j)
c               write (99,*) anname(j)
c               rewind (99)
c               read(99,*)bd
c               read(99,*)ang

               read(bname(j),'(f10.0)')bd
               read(anname(j),'(f10.0)')ang

c               close(99)
               if(ibconn(j).eq.1)then
                  cooxt(j)=0.
                  cooyt(j)=0.
                  coozt(j)=bd
               else
                  cooxt(j)=cooxt(2)
                  cooyt(j)=0.
                  coozt(j)=bd
               endif
            endif
         else
c            rewind (99)
c            write (99,*) bname(j)
            read(bname(j),'(f10.0)')bd
c            write (99,*) anname(j)
            read(anname(j),'(f10.0)')ang
            
            if(j.gt.3) then
               do ij=1,nint
                  if(intcoor(ij).eq.dname(j))then
                     iangname=1
c                     write(99,*)'iname is ',iname
c                     write (99,*) xintt(ij)
                     dihed=xintt(ij)
c                     stop
                  endif
               enddo
               if(iangname.eq.0)then
c                  write (99,*) dname(j)
                  read(dname(j),'(f10.0)')dihed
c                  iangname=0
               endif
            endif
c            rewind (99)
c            read(99,*)bd
c            read(99,*)ang
c            if(j.gt.3)read(99,*)dihed
c            close(99)

c            write(*,*)'bd is ',bd
c            write(*,*)'ang is ', ang
c            write(*,*)'dihed is ', dihed
c            stop
c            write(*,*)'check dummy iangname',iangname
c            write(*,*)'check dummy ireact',ireact

            if(iangname.eq.0)then
c               xa=0.
c               ya=0.
c               za=0.
            call zmat_to_xyz(xa,ya,za,cooxt(ibconn(j)),cooyt(ibconn(j)),
     $ coozt(ibconn(j)),cooxt(iaconn(j)),cooyt(iaconn(j)),
     $ coozt(iaconn(j)),
     $ cooxt(idconn(j)),cooyt(idconn(j)),coozt(idconn(j)),
     $ bd,ang,dihed)
            else if (iangname.eq.1) then
               dist1=bd
               iat2=ibconn(j)
               iat1=iaconn(j)
               ang1=ang
c               write(*,*) 'working on that'
c               write(*,*)'dist1 ', dist1
c               write(*,*)'iat1 ', iat1
               icorr=0
               do jk=1,iat1-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iat1=iat1-icorr
c               write(*,*)'iat1 ', iat1
               icorr=0
               do jk=1,iat2-1
                  if(idummy(jk).eq.1)icorr=icorr+1
               enddo
               iat2=iat2-icorr
c               write(*,*)'iat2 ', iat2
c               write(*,*)'ang1 is ',ang
c               write(*,*)'j is ',j
               jind=j
               icheckang=0
               do jind=j+1,natomt
                 if(ibconn(jind).eq.ibconn(j).and.iaconn(jind).eq.j)then
                  write(*,*)'recognized 2 successive ang dihedral def'
                  icheckang=1
c                  iat3=j+1
                  iat3=jind
c                  write(*,*)'iat3 ', iat3
                  icorr=0
                  do jk=1,iat3-1
                     if(idummy(jk).eq.1)icorr=icorr+1
                  enddo
                  iat3=iat3-icorr
c                  write(*,*)'iat3 ', iat3
                  do ij=1,nint
c                     if(intcoor(ij).eq.anname(j+1))then
                     if(intcoor(ij).eq.anname(jind))then
                        write (*,*)'failed to recog dummy 2 ang conf'
                        write (*,*)'stopping now'
                        stop
                     endif
                  enddo
c                  open (unit=99,status='unknown')

c                  rewind (99)
c                  write (99,*) anname(j+1)
c                  write (99,*) anname(jind)
c                  rewind (99) 
c                  read(99,*) ang2
                  read(anname(jind),'(f10.0)')ang2

c                  close(99)
c                  write(*,*)'ang2 is ',ang2
                 endif
               enddo

               if(icheckang.eq.0)then
                  write (*,*)'failed2 to recog dummy 2 ang conf'
                  write (*,*)'stopping now'
                  stop
               endif

               dist2=sqrt((coox(iat1)-coox(iat2))**2+
     $                    (cooy(iat1)-cooy(iat2))**2+
     $                    (cooz(iat1)-cooz(iat2))**2)
               dist3=sqrt((coox(iat2)-coox(iat3))**2+
     $                    (cooy(iat2)-cooy(iat3))**2+
     $                    (cooz(iat2)-cooz(iat3))**2)
c               write(*,*)'dist2 ', dist2
c               write(*,*)'dist3 ', dist3
               adist2=dist2**2+(dist1+dist2*
     $                cos(3.1415/180*(180-ang1)))**2
               bdist2=dist1*dist1
               cdist2=dist3**2+(dist1+dist3*
     $                cos(3.1415/180*(180-ang1)))**2
c               write(*,*)'adist2 ', adist2
c               write(*,*)'bdist2 ', bdist2
c               write(*,*)'cdist2 ', cdist2
               call twoan_to_xyz(adist2,bdist2,cdist2,coox(iat1)
     $  ,cooy(iat1),cooz(iat1),coox(iat2),cooy(iat2),cooz(iat2),
     $  coox(iat3),cooy(iat3),cooz(iat3),xa,ya,za)
c               stop
            endif

            cooxt(j)=xa
            cooyt(j)=ya
            coozt(j)=za
c            write(*,*)'xa is ,',xa
c            write(*,*)'ya is ,',ya
c            write(*,*)'za is ,',za

         endif
      enddo

c      close(99)
c      write(*,*)'ok 1'
c      stop
c      ifilu=0
      if(ifilu.eq.1)then         
         open (unit=99,file='struct.xyz',status='unknown')
         write(99,*)natomt
         write(99,*)'xyz zmat from update z.mat with dummies'
         do j=1,natomt
            aname1(j)=atname(j)
            write(99,*)aname1(j),cooxt(j),cooyt(j),coozt(j)
         enddo
         close(99)
      endif
c      stop


cc  update  nint intcoor variables

      do j=1,nint
c         write(*,*)'updating var j ',j,intcoor(j)
         ind=0
         do i=1,natomt
            if(intcoor(j).eq.bname(i)) then
               xint(j)=sqrt((cooxt(i)-cooxt(ibconn(i)))**2
     $        +(cooyt(i)-cooyt(ibconn(i)))**2
     $        +(coozt(i)-coozt(ibconn(i)))**2)      
            endif
            if(intcoor(j).eq.anname(i)) then
c               open (unit=99,file='dihed.dat',status='unknown')
c               write(99,*)cooxt(i),cooyt(i),coozt(i)
c               ianum=ibconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               ianum=iaconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               ianum=idconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               close(99)
c               call dihedral
c               open(unit=23,file='dihed.res',status='unknown')
c               read(23,*)xint(j)
c               close(23)
               ib=ibconn(i)
               ia=iaconn(i)
               id=idconn(i)
               if(i.le.3)then
c use fake coordinates to avoid overflow. Impacts dihedral, but only angle matters
                  call xyz_to_zmat(cooxt(i),cooyt(i),coozt(i),cooxt(ib),
     +                cooyt(ib),coozt(ib),cooxt(ia),cooyt(ia),coozt(ia),
     +                cooxt(ia)+0.1,cooyt(ia)+0.1,coozt(ia)+0.1,ang,
     +                dihed)                   
               else
                  call xyz_to_zmat(cooxt(i),cooyt(i),coozt(i),cooxt(ib),
     +                cooyt(ib),coozt(ib),cooxt(ia),cooyt(ia),coozt(ia),
     +                cooxt(id),cooyt(id),coozt(id),ang,dihed)                   
               endif
c               write(*,*)'ib is ',ib
c               write(*,*)'ia is ',ia
c               write(*,*)'id is ',id
c               write(*,*)'i is ',i
c               stop

               xint(j)=ang

            endif
            if(intcoor(j).eq.dname(i)) then
c               open (unit=99,file='dihed.dat',status='unknown')
c               write(99,*)j,i
c               write(99,*)'dname ',dname(i)
c               write(99,*)cooxt(i),cooyt(i),coozt(i)
c               ianum=ibconn(i)
               ib=ibconn(i)
c               write(99,*)cooxt(ib),cooyt(ib),coozt(ib)
c               write(99,*)ib
c               ianum=iaconn(i)
               ia=iaconn(i)
c               write(99,*)cooxt(ia),cooyt(ia),coozt(ia)
c               ianum=idconn(i)
               id=idconn(i)
c               write(99,*)cooxt(id),cooyt(id),coozt(id)
c               close(99)
c               call dihedral
c               open(unit=23,file='dihed.res',status='unknown')
c               read(23,*)cjunk
c               read(23,*)xint(j)
c               close(23)
               call xyz_to_zmat(cooxt(i),cooyt(i),coozt(i),cooxt(ib),
     +              cooyt(ib),coozt(ib),cooxt(ia),cooyt(ia),coozt(ia),
     +              cooxt(id),cooyt(id),coozt(id),ang,dihed) 
               xint(j)=dihed
            endif
         enddo
c         write(*,*)'xint is',j,xint(j)
      enddo
c      write(*,*)'ok 3'
c      stop

      do j=1,ntau
         nprog=natom*3-6-ntau+j
c         write(*,*)'nprog is ',nprog
         do i=1,natomt
            if(bislab(j).eq.bname(i)) then
               tauopt(j)=sqrt((cooxt(i)-cooxt(ibconn(i)))**2+
     $       (cooyt(i)-cooyt(ibconn(i)))**2+(coozt(i)-coozt(ibconn(i))))      
            endif
            if(intcoor(j).eq.anname(i)) then
c               open (unit=99,file='dihed.dat',status='unknown')
c               write(99,*)cooxt(i),cooyt(i),coozt(i)
c               ianum=ibconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               ianum=iaconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               ianum=idconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               close(99)
c               call dihedral
c               open(unit=23,file='dihed.res',status='unknown')
c               read(23,*)tauopt(j)
c               close(23)
               ib=ibconn(i)
               ia=iaconn(i)
               id=idconn(i)
               call xyz_to_zmat(cooxt(i),cooyt(i),coozt(i),cooxt(ib),
     +              cooyt(ib),coozt(ib),cooxt(ia),cooyt(ia),coozt(ia),
     +              cooxt(id),cooyt(id),coozt(id),ang,dihed) 
               tauopt(j)=ang
            endif
            if(intcoor(j).eq.dname(i)) then
c               open (unit=99,file='dihed.dat',status='unknown')
c               write(99,*)cooxt(i),cooyt(i),coozt(i)
c               ianum=ibconn(i)
c               ib=ibconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               ianum=iaconn(i)
c               ia=iaconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               ianum=idconn(i)
c               id=idconn(i)
c               write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c               close(99)
c               call dihedral
c               open(unit=23,file='dihed.res',status='unknown')
c               read(23,*)cjunk
c               read(23,*)tauopt(j)
c               close(23)
               ib=ibconn(i)
               ia=iaconn(i)
               id=idconn(i)
               call xyz_to_zmat(cooxt(i),cooyt(i),coozt(i),cooxt(ib),
     +              cooyt(ib),coozt(ib),cooxt(ia),cooyt(ia),coozt(ia),
     +              cooxt(id),cooyt(id),coozt(id),ang,dihed) 
               tauopt(j)=dihed
            endif
         enddo
c      write(*,*)'ok 4'


         if (tauopt(itau).gt.360.0d0) tauopt(itau) = 
     $        tauopt(itau) - 360.0d0
         if (tauopt(itau).lt.0.0d0) tauopt(itau) = 
     $        tauopt(itau) + 360.0d0
         write(*,*)'tauopt is',j,tauopt(j)
      enddo
cc now update atomlabel if iangindex not equal to zero (dummy atom at position
cc 3)
      
c      if(iangindex.ne.0)then
c      write(*,*)'iangindex is',iangindex
c      write(*,*)'atomlabel ',atomlabel(iangindex)
c      write(*,*)'atname(iangindex) ',atname(iangindex)
c      write(*,*)'bconnt ',bconnt(iangindex)
c      write(*,*)'bname ',bname(iangindex)
c      write(*,*)'aconnt ',aconnt(iangindex)
c      write(*,*)'aname ',anname(iangindex)
c      write(*,*)'dconnt ',dconnt(iangindex)
c      write(*,*)'dname ',dname(iangindex)
c         ia=iangindex
c         open (unit=99,file='dihed.dat',status='unknown')
c         write(99,*)cooxt(ia),cooyt(ia),coozt(ia)
c         ianum=ibconn(ia)
c         write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c         ianum=iaconn(ia)
c         write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c         ianum=idconn(ia)
c         write(99,*)cooxt(ianum),cooyt(ianum),coozt(ianum)
c         close(99)
c         call dihedral
c         open(unit=23,file='dihed.res',status='unknown')
c         read(23,*)cjunk
c         read(23,*)cname
c         close(23)
c         write(*,*)'cname is ',cname
c      stop
c
c         open (unit=99,status='unknown')
c         write (99,2604) atname(ia),bconnt(ia),bname(ia),
c     $        aconnt(ia),anname(ia),dconnt(ia),cname
c         rewind(99)
c         read (99,'(A60)') atomlabel(ia)         
c         close(99)
c         write(*,*)'atomlabel ',atomlabel(iangindex)
c      endif

c 2604 format (1x,3a6,1x,a6,1x,a6,1x,a6,1x,a6)

c      stop

      inda=0
      do j=1,natom
c         if(idummy(j).ne.1)then
c            inda=inda+1
            coox(j)=cooxs(j)
            cooy(j)=cooys(j)
            cooz(j)=coozs(j)
c         endif
      enddo

      return
      end

c************************************************************
      
      subroutine zmat_to_xyz(xa,ya,za,xb,yb,zb,xc,yc,zc,
     $ xd,yd,zd,bd,ang,dihed)

cc returns coordinates for given coords of b,c,d and ab distance, abc angle 
cc and abcd dihedral angle

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'

c this is a geometric approach for calculating the xyz coordinates of atom A
c when the xyz coordinates of the B C and D are known and A position is defined with respect 
c to B C and D in Z-matrix notation.
c the adopted approach consists in translating to B the frame and then rotating 
c so that B C and D are on the xy plane, with C on the y axis
c A coordinates are immediately determined in this frame, which is then
c rototranslated back to the BCD original frame of reference   
c _rt variables are in the rototranslated reference system
c _t variables are in the translated ref system 
c this is similar to the NERF algorithm (see for example 
c Parsons et al. J.Comp.Chem. 26(10) 1063, 2005.)

c first determines coordinates of ABCD in the RT frame of reference
      
c      bc_dist=sqrt((xc-xb)**2+(yc-yb)**2+(zc-zb)**2)
c      bd_dist=sqrt((xb-xd)**2+(yb-yd)**2+(zb-zd)**2)
c      cd_dist=sqrt((xc-xd)**2+(yc-yd)**2+(zc-zd)**2)

c      write(*,*)'bd ',bd
c      write(*,*)'ang ',ang
c      write(*,*)'dihed ',dihed
c      write(*,*)'bc_dist ',bc_dist
c      write(*,*)'bd_dist ',bd_dist
c      write(*,*)'cd_dist ',cd_dist
c      if(cd_dist.gt.bd_dist)dihed=-dihed
c      if(dihed.gt.180.)dihed=dihed-360.

      ang1=ang*pigr/180.
      dihed1=dihed*pigr/180.
c      write(*,*)'cos dihed ',cos(dihed1)

      bdist=bd*cos(ang1)
      ddist=bd*sin(ang1)

      xb_rt=0.
      yb_rt=0.
      zb_rt=0.

      xa_rt=ddist*cos(dihed1)
      ya_rt=bdist
      za_rt=ddist*sin(dihed1)


      bc_dist=sqrt((xc-xb)**2+(yc-yb)**2+(zc-zb)**2)
      bd_dist=sqrt((xb-xd)**2+(yb-yd)**2+(zb-zd)**2)
      cd_dist=sqrt((xc-xd)**2+(yc-yd)**2+(zc-zd)**2)

      xc_rt=0.
      yc_rt=bc_dist
      zc_rt=0.

      yd_rt=(bc_dist**2+bd_dist**2-cd_dist**2)/2./bc_dist
      xd_rt=sqrt(bd_dist**2-yd_rt**2)
      zd_rt=0. 

c      check1=sqrt((za_rt-zd_rt)**2+(ya_rt-yd_rt)**2+(xa_rt-xd_rt)**2)
c      write(*,*)'check1 is ',check1
c      write(*,*)'xa rt = ',xa_rt
c      write(*,*)'ya rt = ',ya_rt
c      write(*,*)'za rt = ',za_rt
c      write(*,*)'xc_rt  = ',xc_rt
c      write(*,*)'yc_rt  = ',yc_rt
c      write(*,*)'zc_rt  = ',zc_rt
c      write(*,*)'xd_rt  = ',xd_rt
c      write(*,*)'yd_rt  = ',yd_rt
c      write(*,*)'zd_rt  = ',zd_rt
c      write(*,*)'xd  = ',xd
c      write(*,*)'yd  = ',yd
c      write(*,*)'zd  = ',zd
      za_rt=-za_rt

cc translate original frame of ref coords so that B is at (0,0,0)

      xb_t=0.
      yb_t=0.
      zb_t=0.
      xc_t=xc-xb
      yc_t=yc-yb
      zc_t=zc-zb
      xd_t=xd-xb
      yd_t=yd-yb
      zd_t=zd-zb


      if(yd_t.eq.0)yd_t=0.00000001
      if(xd_t.eq.0)xd_t=0.00000001
      if(xd_rt.eq.0)xd_rt=0.00000001

cc now determine the rotation matrix to rotate back to the original ref system
cc this is based on determinig the position of a (0,0,1) in the orginal frame
cc and determining the consistent rotation matrix

c      write(*,*)'yc_rt ',yc_rt
c      write(*,*)'xd_rt ',xd_rt
c      check1b=sqrt((za_rt-zd_rt)**2+(ya_rt-yd_rt)**2+(xa_rt-xd_rt)**2)
c      write(*,*)'check1 is ',check1

      r12=(xc-xb)/yc_rt
      r22=(yc-yb)/yc_rt
      r32=(zc-zb)/yc_rt

      r11=(xd-xb-yd_rt*r12)/xd_rt
      r21=(yd-yb-yd_rt*r22)/xd_rt
      r31=(zd-zb-yd_rt*r32)/xd_rt
c
      anum_aconst=yc_t-yd_t/xd_t*xc_t
      den_aconst=zc_t-zd_t/xd_t*xc_t

cc this is to avoid overflows for planar systems

      if(abs(anum_aconst).lt.1.0e-6.and.abs(den_aconst).lt.1.0e-6)then
         if(anum_aconst.lt.0)aconst=-1.0E20
         if(anum_aconst.gt.0)aconst=1.0E20
      else if(abs(den_aconst).lt.1.0e-6)then
         if(anum_aconst.lt.0)aconst=-1.0E20
         if(anum_aconst.gt.0)aconst=1.0E20
      else
         aconst=(yc_t-yd_t/xd_t*xc_t)/(zc_t-zd_t/xd_t*xc_t)
      endif

      den1=(yd_t/xd_t-aconst*zd_t/xd_t)
      if(den1.eq.0.)den1=1.0e-20
      bconst=1./den1

      dvect=1.0

      xe_t=-(dvect)/sqrt(1+(bconst**2)*(1+aconst**2))
      ye_t=-xe_t*bconst
      ze_t=-ye_t*aconst

c      write(*,*)'num_aconst ',anum_aconst
c      write(*,*)'den_aconst ',den_aconst
c      write(*,*)'bconst ',bconst

c      write(*,*)'xc_t ',(xc_t)
c      write(*,*)'xd_t ',(xd_t)
c      write(*,*)'yc_t ',(yc_t)
c      write(*,*)'yd_t ',(yd_t)
c      write(*,*)'zc_t ',(zc_t)
c      write(*,*)'zd_t ',(zd_t)
c      write(*,*)'den ',(yd_t/xd_t-aconst*zd_t/xd_t)
c      write(*,*)'den1 ',(yd_t/xd_t)
c      write(*,*)'den2 ',(aconst*zd_t/xd_t)
c      write(*,*)'den21 ',(aconst)
c      write(*,*)'den22 ',(zd_t/xd_t)
c      write(*,*)'xd ',xd_t
c      write(*,*)'yd ',yd_t
c      write(*,*)'zd ',zd_t

      r13=xe_t/dvect
      r23=ye_t/dvect
      r33=ze_t/dvect

      r13n=-xe_t/dvect
      r23n=-ye_t/dvect
      r33n=-ze_t/dvect

c      xe_rt=0.
c      ye_rt=0.
c      ze_rt=1.

cc now rotate and translate back

cc here I check  the (001) vector direction to decide whether to take the positive of negative results of the 
cc square root taken above
c      write(*,*)'r12 is ',r12
c      write(*,*)'r22 is ',r22
c      write(*,*)'r32 is ',r32
c      write(*,*)'r13 is ',r13
c      write(*,*)'r23 is ',r23
c      write(*,*)'r33 is ',r33

      xap=xb+r11*xa_rt+r12*ya_rt+r13*za_rt
      yap=yb+r21*xa_rt+r22*ya_rt+r23*za_rt
      zap=zb+r31*xa_rt+r32*ya_rt+r33*za_rt

      xan=xb+r11*xa_rt+r12*ya_rt+r13n*za_rt
      yan=yb+r21*xa_rt+r22*ya_rt+r23n*za_rt
      zan=zb+r31*xa_rt+r32*ya_rt+r33n*za_rt

      b1v=xb-xc
      b2v=yb-yc
      b3v=zb-zc
      c1v=xc-xd
      c2v=yc-yd
      c3v=zc-zd
      vec1=b2v*c3v-b3v*c2v
      vec1c=xe_t
      vec2=b3v*c1v-b1v*c3v
      vec2c=ye_t
      vec3=b1v*c2v-b2v*c1v
      vec3c=ze_t

      if(abs(vec1c).gt.1.0e-5)then
         checkv=vec1/vec1c
      else if(abs(vec2c).gt.1.0e-5)then
         checkv=vec2/vec2c
      else
         checkv=vec3/vec3c
      endif
cc check
cc      if(checkv.lt.0)then
      if(checkv.gt.0)then
         xa=xap
         ya=yap
         za=zap
cc      else if (checkv.gt.0)then
      else if (checkv.lt.0)then
         xa=xan
         ya=yan
         za=zan
      else
         xa=xap
         ya=yap
         za=zap
      endif

c      check2=sqrt((za-zd)**2+(ya-yd)**2+(xa-xd)**2)
c      check3=sqrt((za-zd)**2+(ya-yd)**2+(xa-xd)**2)
c      check2n=sqrt((zan-zd)**2+(yan-yd)**2+(xan-xd)**2)
c      write(*,*)'vec1 is ',vec1
c      write(*,*)'vec1c is ',vec1c
c      write(*,*)'vec2c is ',vec2c
c      write(*,*)'vec2 is ',vec2
c      write(*,*)'checkv is ',checkv


c      write(*,*)'check2 is ',check2
c      write(*,*)'check2n is ',check2n
c      write(*,*)'xan,p is ',xan,xap
c      write(*,*)'yan,p is ',yan,yap
c      write(*,*)'zan,p is ',zan,zap

      return
      end

c************************************************************
      subroutine rototrasl(natom,coox,cooy,cooz,ilin)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
 
c     
      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension cooxt(natommx),cooyt(natommx),coozt(natommx)
      dimension cooxt1(natommx),cooyt1(natommx),coozt1(natommx)

      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'


cc performs traslation to center of axis
cc and zy rotation to get second atom on x axis
      if(ilin.eq.0)then
         do j=1,natom
            cooxt(j)=coox(j)-coox(1)
            cooyt(j)=cooy(j)-cooy(1)
            coozt(j)=cooz(j)-cooz(1)
         enddo

c         do j=1,natom
c            write(*,*)'coox ',cooxt(j)
c            write(*,*)'cooy ',cooyt(j)
c            write(*,*)'cooz ',coozt(j)
c         enddo
c         stop

cc first rotate around z
         ialfajump=0
         dist=sqrt(cooxt(2)*cooxt(2)+cooyt(2)*cooyt(2))
         if(dist.gt.1.0e-10) then
            alfa=dacos((cooxt(2))/dist)
         else
            alfa=0.
            ialfajump=1
         endif
c         if(cooyt(2).gt.0)then
c            if(alfa.lt.0)alfa=3.1415-alfa
c         endif
c         if(cooyt(2).lt.0)then
c            if(alfa.gt.0)alfa=3.1415*2-alfa
c            if(alfa.lt.0)alfa=3.1415-alfa
c         endif

c         write(*,*)'dist is = ',dist
c         write(*,*)'alfa is = ',alfa

         if(cooyt(2).gt.0)then
            alfa=pigr*2-alfa
c            if(cooxt(2).gt.0)alfa=3.1415-alfa+3.1415
c            if(cooxt(2).lt.0)alfa=3.1415-alfa+3.1415
c            if(alfa.lt.0)alfa=3.1415-abs(alfa)
         endif

c      sinpar=(cooyt(2))/dist

c         write(*,*)'alfa z = ',alfa
c      write(*,*)'alfa y = ',asin(sinpar)
         if(ialfajump.eq.0)then
            do j=1,natom
               cooxt1(j)=cooxt(j)*dcos(alfa)-cooyt(j)*dsin(alfa)
               cooyt1(j)=cooxt(j)*dsin(alfa)+cooyt(j)*dcos(alfa)
               coozt1(j)=coozt(j)
            enddo
         else
            do j=1,natom
               cooxt1(j)=cooxt(j)
               cooyt1(j)=cooyt(j)
               coozt1(j)=coozt(j)
            enddo
         endif
c         do j=1,natom
c            write(*,*)'coox1rt ',cooxt1(j)
c            write(*,*)'cooy1rt ',cooyt1(j)
c            write(*,*)'cooz1rt ',coozt1(j)
c         enddo


cc now rotate around y
         ialfajump=0
         dist=sqrt(cooxt1(2)*cooxt1(2)+coozt1(2)*coozt1(2))
         if(dist.gt.1.0e-10) then
            alfa=dacos((cooxt1(2))/dist)
         else
            ialfajump=1
            alfa=0
         endif
c      write(*,*)'alfa y = ',alfa
c      write(*,*)'dist = ',dist
c      write(*,*)'alfajump = ',alfajump
c      write(*,*)'ialfa y = ',ialfajump
c      do j=1,natom
c         write(*,*)cooxt1(j),cooyt1(j),coozt1(j)
c      enddo

         if(coozt1(2).lt.0)then
            alfa=pigr*2.-alfa
c            if(cooxt(2).gt.0)alfa=3.1415-alfa+3.1415
c            if(cooxt(2).lt.0)alfa=3.1415-alfa+3.1415
c            if(alfa.lt.0)alfa=3.1415-abs(alfa)
         endif

         if(ialfajump.eq.0)then
            do j=1,natom
               coox(j)=cooxt1(j)*dcos(alfa)+coozt1(j)*dsin(alfa)
               cooy(j)=cooyt1(j)
               cooz(j)=-cooxt1(j)*dsin(alfa)+coozt1(j)*dcos(alfa)
            enddo
         else
            do j=1,natom
               coox(j)=cooxt1(j)
               cooy(j)=cooyt1(j)
               cooz(j)=coozt1(j)
            enddo
         endif
c         do j=1,natom
c            write(*,*)'coox2 ',coox(j)
c            write(*,*)'cooy2 ',cooy(j)
c            write(*,*)'cooz2 ',cooz(j)
c         enddo


c      do j=1,natom
c         write(*,*)coox(j),cooy(j),cooz(j)
c      enddo
c      write(*,*)'ok 2'
      else
cc check molecule axis: x (1), y(2) or z(3)
c         write(*,*)'arrived here'
c         stop

         cooxtot=0.
         cooytot=0.
         cooztot=0.
         do j=1,natom
            cooxtot=cooxtot+abs(coox(j))
            cooytot=cooytot+abs(cooy(j))
            cooztot=cooztot+abs(cooz(j))
         enddo
         idir=0
         if(cooxtot.gt.cooytot.and.cooxtot.gt.cooztot)idir=1
         if(cooytot.gt.cooxtot.and.cooytot.gt.cooztot)idir=2
         if(cooztot.gt.cooxtot.and.cooztot.gt.cooytot)idir=3
         if(idir.eq.2)then
            do j=1,natom
               coox(j)=cooy(j)
               cooy(j)=0.
            enddo
            do j=1,natom
               coox(j)=coox(j)-coox(1)
            enddo
         endif
         if(idir.eq.3)then
            do j=1,natom
               coox(j)=cooz(j)
               cooz(j)=0.
            enddo
            do j=1,natom
               coox(j)=coox(j)-coox(1)
            enddo
         endif
      endif

      return
      end

c************************************************************
         subroutine twoan_to_xyz(da2,db2,dc2,x1,y1,z1,x2,y2,z2,x3,y3,z3
     $  ,xa,ya,za)

cc returns coordinates xa ya za for an atom defined with respect to
cc three atoms whose coordinates are known and for which the distances 
cc are known

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'

cc first check if reference atoms are co-linear:

      distcheck1=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
      ref1=sqrt(abs(dc2-da2))+sqrt(abs(db2-da2))
c         write(*,*)'dist1',distcheck1
c         write(*,*)'dist2',ref1
c         write(*,*)'da2 ',da2
c         write(*,*)'db2 ',db2
c         write(*,*)'dc2 ',dc2
c         stop
      if(abs(distcheck1-ref1).lt.0.001)then
c         write(*,*)'i am here'
cc if system colinear place x on z plane
c         write(*,*)'3-1-2 colinear atoms'
cc translate ref system
         xt1=0.
         yt1=0.
         zt1=0.
         xt2=x2-x1
         yt2=y2-y1
         zt2=z2-z1
         xt3=x3-x1
         yt3=y3-y1
         zt3=z3-z1
         if(yt2.ne.0)then
            alfa=-xt2/yt2
            beta=(db2-da2-xt2*xt2-yt2*yt2-zt2*zt2)/2/yt2
            xa=-alfa*beta/(1+alfa*alfa)+sqrt((alfa*beta/
     +            (1+alfa*alfa))**2-
     +           (beta*beta-da2)/(1+alfa*alfa))
            ya=alfa*xa-beta
            za=0.
         else
            xa=-(db2-da2-xt2*xt2-yt2*yt2-zt2*zt2)/2./xt2
            ya=sqrt(abs(da2-xa*xa))
            za=0.
         endif
         xa=xa+x1
         ya=ya+y1
         za=za+z1
      else
         itrasl=0   
         if(x2.eq.0.and.y2.eq.0.and.z2.eq.0)then
            xt1=0.
            yt1=0.
            zt1=0.
            xt2=x2-x1
            yt2=y2-y1
            zt2=z2-z1
            xt3=x3-x1
            yt3=y3-y1
            zt3=z3-z1
         else if(x1.eq.0.and.y1.eq.0.and.z1.eq.0)then
            xt1=x1
            yt1=y1
            zt1=z1
            xt2=x2
            yt2=y2
            zt2=z2
            xt3=x3
            yt3=y3
            zt3=z3
         else
            itrasl=1
            xt1=0.
            yt1=0.
            zt1=0.
            xt2=x2-x1
            yt2=y2-y1
            zt2=z2-z1
            xt3=x3-x1
            yt3=y3-y1
            zt3=z3-z1
c            write(*,*)'this algorithm is not going to work'
c            write(*,*)'this is going to need some debugging'
c            stop
         endif
         if(itrasl.ne.1)then
            xa=(da2-db2+xt2*xt2)/2/xt2
            if(yt3.ne.0.and.zt3.ne.0)then
               Ac=-zt3/yt3
               Bc=-(dc2-da2+2*xa*xt3-xt3*xt3-yt3*yt3-zt3*zt3)/2/yt3
               za=(-Ac*Bc+sqrt(abs(Ac*Ac*Bc*Bc-(1+Ac*Ac)*
     +           (xta*xta-da2+Bc*Bc))))/(1+Ac*Ac)
               ya=Ac*za+Bc
            else if (yt3.ne.0.and.zt3.eq.0)then
               ya=-(dc2-da2+2*xa*xt3-xt3*xt3-yt3*yt3-zt3*zt3)/2/yt3
               za=sqrt(abs(da2-xa*xa-ya*ya))
            else
               za=-(dc2-da2+2*xa*xt3-xt3*xt3-yt3*yt3-zt3*zt3)/2/zt3
               ya=sqrt(abs(da2-xa*xa-za*za))
            endif
            if(x2.eq.0.and.y2.eq.0.and.z2.eq.0)then
               xa=xa+x1
               ya=ya+y1
               za=za+z1
            endif
         else
            a1=-yt2/xt2
            a2=-zt2/xt2
            a3=(da2+xt2**2+yt2**2+zt2**2-db2)/2./xt2
            b1=(2*a2*xt3+2*zt3)/(-2*a1*xt3-2*yt3)
            b2=(dc2-xt3**2-yt3**2-zt3**2-da2+2*a3*xt3)/(-2*a1*xt3-2*yt3)
            c1=a1*b1+a2
            c2=a1*b2+a3
            aa=c1**2+b1**2+1.
            bb=2*c1*c2+2*b1*b2
            cc=c2**2+b2**2-da2
            delta=bb**2-4*aa*cc
c            write(*,*)'bb is',bb
c            write(*,*)'cc is',cc
c            write(*,*)'delta is ',delta
cc
cc taking the positive value of the quadratic should not impact
cc the reconstruction of the z-matrix
c
            za=(-bb+sqrt(bb**2-4*aa*cc))/2/aa
            ya=za*b1+b2
            xa=za*c1+c2
            xa=xa+x1
            ya=ya+y1
            za=za+z1
         endif

c         xa=(da2-db2+x2*x2)/2/x2
c         Ac=-z3/y3
c         Bc=-(dc2-da2+2*xa*x3-x3*x3-y3*y3-z3*z3)/2/y3
c         za=(-Ac*Bc+sqrt(abs(Ac*Ac*Bc*Bc-(1+Ac*Ac)*(xa*xa-da2+Bc*Bc))))
c     +        /(1+Ac*Ac)
c         ya=Ac*za+Bc
      endif

c      write(*,*)'xa is ',xa
c      write(*,*)'ya is ',ya
c      write(*,*)'za is ',za
c      write(*,*)'sqrt is ',Ac*Ac*Bc*Bc-(1+Ac*Ac)*(xa*xa-da2+Bc*Bc)

c      xa=0.
c      ya=0.
c      za=0.
c      if(abs(xa-x1).lt.0.1)xa=0.5
c      if(abs(ya-y1).lt.0.1)ya=0.5
c      if(abs(za-z1).lt.0.1)za=0.5
c      ya=1.5
c      za=1.5

c      do j=1,100
c
c         d1=-((xa-x1)**2+(ya-y1)**2+(za-z1)**2-da2)/2
c         d2=-((xa-x2)**2+(ya-y2)**2+(za-z2)**2-db2)/2
c         d3=-((xa-x3)**2+(ya-y3)**2+(za-z3)**2-dc2)/2
c         a1=xa-x1
c         a2=xa-x2
c         a3=xa-x3
c         b1=ya-y1
c         b2=ya-y2
c         b3=ya-y3
c         c1=za-z1
c         c2=za-z2
c         c3=za-z3

c         write(*,*)'d1 is',d1
c         write(*,*)'d2 is',d2
c         write(*,*)'d3 is',d3
      
c         gr_i=(d2-d1/a1*a2)/(b2-b1/a1*a2)
c         gr_h=(c1/a1*a2-c2)/(b2-b1/a1*a2)
c         gr_m=-gr_h*b1/a1*a3+b3*gr_h-c1/a1*a3+c3
c         gr_n=d3-d1/a1*a3-gr_i*(-b1/a1*a3+b3)
c         delta3=gr_n/gr_m
c         delta2=delta3*gr_h+gr_i
c         delta1=d1/a1-b1/a1*delta2-c1/a1*delta3

c         write(*,*)'delta1 is',delta1
c         write(*,*)'delta2 is',delta2
c         write(*,*)'delta3 is',delta3
         
c         xa=xa+delta1
c         ya=ya+delta2
c         za=za+delta3
c         error_tot=abs(d1)+abs(d2)+abs(d3)
c         if (error_tot.lt.1.0e-10) goto 100
         
c         write(*,*)'xa is',xa
c         write(*,*)'ya is',ya
c         write(*,*)'za is',za

c         write(*,*)' error 1 is ',d1
c         write(*,*)' error 2 is ',d2
c         write(*,*)' error 2 is ',d3
c      enddo
c 100  continue
c      if (error_tot.gt.1.0) then
c         write(*,*)'failed xyz translation for dummy atoms '
c         write(*,*)'error in zmat-update'
c         stop
c      endif


c      write(*,*)'convergence reached in dummy atom xyz coo search'
c      stop

      return
      end
c************************************************************
      subroutine distang(da,db,dc,alfa)
c
c returns dc given da,db and alfa
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'
c      pi=3.14159265359
      alfa=pigr*alfa/180.

      da1=db*cos(pigr-alfa)
      db1=db*sin(pigr-alfa)

      dc=sqrt((da+da1)**2+db1**2)

      return
      end
c************************************************************
      subroutine distang2(da,db,dc,alfa)
c
c returns alfa given da, db, dc
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'
c      write(*,*)'da is',da
c      write(*,*)'db is',db
c      write(*,*)'dc is',dc
c      pi=3.14159265359
      beta=acos((dc*dc-da*da-db*db)/2/da/db)
c     write(*,*)'acos is ',(dc*dc-da*da-db*db)/2/da/db
c     write(*,*)'beta is ',beta
c      stop
      alfa=(pigr-beta)*180./pigr

      return
      end
c************************************************************
      subroutine xyz_to_zmat(xa,ya,za,xb,yb,zb,xc,yc,zc,
     $ xd,yd,zd,ang,dihed)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'
c     
      dimension a_b(3),b_c(3),c_d(3)
      dimension vn1(3),vn2(3),vm(3)


      LOGICAL leof,lsec,ltit

      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3
     $ ,title,title1,word4,word5,word6,word7

c      include 'filcomm.f'

      xbt=xb-xa
      ybt=yb-ya
      zbt=zb-za
      xct=xc-xa
      yct=yc-ya
      zct=zc-za
      xdt=xd-xa
      ydt=yd-ya
      zdt=zd-za
      xat=0.
      yat=0.
      zat=0.

      a_b(1)=xat-xbt
      a_b(2)=yat-ybt
      a_b(3)=zat-zbt

      b_c(1)=xbt-xct
      b_c(2)=ybt-yct
      b_c(3)=zbt-zct

      c_d(1)=xct-xdt
      c_d(2)=yct-ydt
      c_d(3)=zct-zdt

      axb=a_b(1)*b_c(1)+a_b(2)*b_c(2)+a_b(3)*b_c(3)
      anorm=sqrt(a_b(1)*a_b(1)+a_b(2)*a_b(2)+a_b(3)*a_b(3))
      bnorm=sqrt(b_c(1)*b_c(1)+b_c(2)*b_c(2)+b_c(3)*b_c(3))
      cnorm=sqrt(c_d(1)*c_d(1)+c_d(2)*c_d(2)+c_d(3)*c_d(3))

c      write(*,*)'xa ya za ',xat,yat,zat
c      write(*,*)'xb yb zb ',xbt,ybt,zbt
c      write(*,*)'xc yc zc ',xct,yct,zct
c      write(*,*)'xd yd zd ',xdt,ydt,zdt
      if(anorm.lt.1.0e-20)then
         anorm=1.0e-20
      endif
      if(bnorm.lt.1.0e-20)then
         bnorm=1.0e-20
      endif
c      write(*,*)'cnorm is ',cnorm
c      write(*,*)'axb is ',axb
c      write(*,*)'acos is ',acos(axb/anorm/bnorm)

      ang=180./pigr*(pigr-acos(axb/anorm/bnorm))

c      write(*,*)'ang is ',ang

      a_b(1)=a_b(1)/anorm
      a_b(2)=a_b(2)/anorm
      a_b(3)=a_b(3)/anorm

      b_c(1)=b_c(1)/bnorm
      b_c(2)=b_c(2)/bnorm
      b_c(3)=b_c(3)/bnorm

      c_d(1)=c_d(1)/cnorm
      c_d(2)=c_d(2)/cnorm
      c_d(3)=c_d(3)/cnorm

      vn1(1)= a_b(2)*b_c(3)-a_b(3)*b_c(2)
      vn1(2)=-a_b(1)*b_c(3)+a_b(3)*b_c(1)
      vn1(3)= a_b(1)*b_c(2)-a_b(2)*b_c(1)

      vn2(1)= b_c(2)*c_d(3)-b_c(3)*c_d(2)
      vn2(2)=-b_c(1)*c_d(3)+b_c(3)*c_d(1)
      vn2(3)= b_c(1)*c_d(2)-b_c(2)*c_d(1)

      vm(1)=  b_c(2)*vn1(3)-b_c(3)*vn1(2)
      vm(2)= -b_c(1)*vn1(3)+b_c(3)*vn1(1)
      vm(3)=  b_c(1)*vn1(2)-b_c(2)*vn1(1)

      vx=vn1(1)*vn2(1)+vn1(2)*vn2(2)+vn1(3)*vn2(3)
      vy=vm(1)*vn2(1)+vm(2)*vn2(2)+vm(3)*vn2(3)

c      write(*,*)'atan is ',atan2(vy,vx)
c      write(*,*)'acos1 is ',acos(vx)
c      write(*,*)'acos2 is ',acos(vy)
c      stop

      dihed=-atan2(vy,vx)*180./pigr
      dihed2=acos(vx)*180./pigr
c      dihed=atan2(vx,vy)*180./pigr

c      write(*,*)'ang is ',ang
c      write(*,*)'dihed is ',dihed
c      write(*,*)'dihed 2 is',dihed2

c      stop

      return
      end
