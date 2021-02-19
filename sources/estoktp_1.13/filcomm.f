c
c common file for estoktop.f
c
      character*20 stoich

      common /keyword/ leof,lsec,ltit
      common /key/line,sename,string,word,word2,word3,word4,word5,
     $ word6,word7,title,title1
      common /nodes/ nodenum,numproc,numprocll,numprochl
      common /mem/gmem
      common /debug/ idebug
      common /tsdat/ aabs1,babs1,aabs2,babs2,babs3,rts
      common /itsdat/ iabs,iadd,iiso,ibeta,ibarr,ifrozrts,isite,jsite
     $ ,ksite,ireact,irecov,iprojrcoo,ibstep
      common /ircdat/ivar,iresirc,imdtunn,inotunn,intfreq,
     $ irotd_lr,irotd_sr
      common /irtype/ irw,ipw,nts,ip1,ip2,ipr1
c      common /strucword/ stoich
      common /igeom/ frozcoo,igeom_wellp,igeom_wellr

