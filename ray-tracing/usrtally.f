************************************************************************
*                                                                      *
      subroutine usrtally(ncol)
*                                                                      *
*        sample subroutine for user defined tally.                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       ncol = 1 : start of calculation                                *
*              2 : end of calculation                                  *
*              3 : end of a batch                                      *
*              4 : source                                              *
*              5 : detection of geometry error                         *
*              6 : recovery of geometry error                          *
*              7 : termination by geometry error                       *
*              8 : termination by weight cut-off                       *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : (n,x) reaction                                      *
*             14 : (n,n'x) reaction                                    *
*             15 : sequential transport only for tally                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        In the distributed-memory parallel computing,                 *
*          npe : total number of used Processor Elements               *
*          me : ID number of each processor                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        In the shared memory parallel computing,                      *
*          ipomp : ID number of each core                              *
*          npomp : total number of used cores                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nocas : current event number in this batch                    :
*        nobch : current batch number                                  *
*        rcasc : real number of NOCAS+maxcas*(NOBCH-1)                 *
*        rsouin : sum of the weight of source particle                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        no : cascade id in this event                                 *
*        idmn(mat) : material id                                       *
*        ityp : particle type                                          *
*        ktyp : particle kf-code                                       *
*        jtyp : charge number of the particle                          *
*        mtyp : baryon number of the particle                          *
*        rtyp : rest mass of the particle (MeV)                        *
*        oldwt : weight of the particle at (x,y,z)                     *
*        qs : dE/dx of electron at (x,y,z)                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        iblz1 : cell id at (x,y,z)                                    *
*        iblz2 : cell id after crossing                                *
*        ilev1 : level structure id of the cell at (x,y,z)             *
*          ilat1(i,j)                                                  *
*        ilev2 : level structure id of the cell after crossing         *
*          ilat2(i,j)                                                  *
*        costha : cosine of theta on surface crossing                  *
*        uang(1) : x of position(?) at surface crossing                *
*        uang(2) : y of position(?) at surface crossing                *
*        uang(3) : z of position(?) at surface crossing                *
*        nsurf : surface number                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        name(no,ipomp+1) : collision number of the particle                   *
*        ncnt(1,no,ipomp+1) : values of counter 1                              *
*        ncnt(2,no,ipomp+1) : values of counter 2                              *
*        ncnt(3,no,ipomp+1) : values of counter 3                              *
*        wt(no,ipomp+1) : weight of the particle at (xc,yc,zc)                 *
*        u(no,ipomp+1) : x, y, z-components of unit vector of                  *
*        v(no,ipomp+1) :momentum of the particle                               *
*        w(no,ipomp+1) :                                                       *
*        e(no,ipomp+1) : energy of the particle at (x,y,z) (MeV)               *
*        t(no,ipomp+1) : time of the particle at (x,y,z) (nsec)                *
*        x(no,ipomp+1) : x, y, z-position coordinates of                       *
*        y(no,ipomp+1) :the preceding event point (cm)                         *
*        z(no,ipomp+1) :                                                       *
*        ec(no,ipomp+1) : energy of the particle at (xc,yc,zc) (MeV)           *
*        tc(no,ipomp+1) : time of the particle at (xc,yc,zc) (nsec)            *
*        xc(no,ipomp+1) : x, y, z-position coordinates of                      *
*        yc(no,ipomp+1) :the particle (cm)                                     *
*        zc(no,ipomp+1) :                                                      *
*        spx(no,ipomp+1) : x, y, z-components of unit vector of                *
*        spy(no,ipomp+1) :spin direction of the particle                       *
*        spz(no,ipomp+1) :                                                     *
*        nzst(no,ipomp+1)                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts : the number of produced particle and nucleus          *
*        mathz : Z number of the mother nucleus                        *
*        mathn : N number of the mother nucleus                        *
*        jcoll : reaction type id1                                     *
*        kcoll : reaction type id2                                     *
*                                                                      *
*----------------------------------------------------------------------*
*        jcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, High Energy Nuclear collisions                  *
*                =  5, Heavy Ion reactions                             *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
*                =  8, Electron reactions by data                      *
*                =  9, Proton reactions by data                        *
*                = 10, Neutron event mode                              *
*                = 11, Delta Ray production                            *
*                = 12, Muon interaction                                *
*                = 13, Photon by EGS5                                  *
*                = 14, Electron by EGS5                                *
*                                                                      *
*----------------------------------------------------------------------*
*        kcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight                                         *
*                  = 9, time                                           *
*                  = 10, x                                             *
*                  = 11, y                                             *
*                  = 12, z                                             *
*                                                                      *
*        jcount(i,nclsts)                                              *
*                                                                      *
************************************************************************
*                                                                      *
*        kf code table                                                 *
*                                                                      *
*           kf-code: ityp :  description                               *
*                                                                      *
*             2212 :   1  :  proton                                    *
*             2112 :   2  :  neutron                                   *
*              211 :   3  :  pion (+)                                  *
*              111 :   4  :  pion (0)                                  *
*             -211 :   5  :  pion (-)                                  *
*              -13 :   6  :  muon (+)                                  *
*               13 :   7  :  muon (-)                                  *
*              321 :   8  :  kaon (+)                                  *
*              311 :   9  :  kaon (0)                                  *
*             -321 :  10  :  kaon (-)                                  *
*                                                                      *
*               11 :  12  :  electron                                  *
*              -11 :  13  :  positron                                  *
*               22 :  14  :  photon                                    *
*          1000002 :  15  :  deuteron                                  *
*          1000003 :  16  :  triton                                    *
*          2000003 :  17  :  3he                                       *
*          2000004 :  18  :  alpha                                     *
*      Z*1000000+A :  19  :  nucleus                                   *
*                                                                      *
*           kf-code of the other transport particles (ityp=11)         *
*               12 :         nu_e                                      *
*               14 :         nu_mu                                     *
*              221 :         eta                                       *
*              331 :         eta'                                      *
*             -311 :         k0bar                                     *
*            -2112 :         nbar                                      *
*            -2212 :         pbar                                      *
*             3122 :         Lanbda0                                   *
*             3222 :         Sigma+                                    *
*             3212 :         Sigma0                                    *
*             3112 :         Sigma-                                    *
*             3322 :         Xi0                                       *
*             3312 :         Xi-                                       *
*             3334 :         Omega-                                    *
*                                                                      *
************************************************************************

      use MMBANKMOD !FURUTA
      use usrtalmod

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /cusrtally/ iusrtally, iudtf(50)
      common /cudtpara/ udtpara(0:9)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

cKN 2012/10/31
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb,iransb64 
!$OMP THREADPRIVATE(/randtp/)

      common /taliin/ rsouin, nzztin, nrgnin
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

      integer nudtvar
      common /cnudtvar/ nudtvar

      real*8 xold,yold,zold,dist,totalflux
      save xold,yold,zold,dist,totalflux
!$OMP THREADPRIVATE(xold,yold,zold,dist,totalflux)
!      save totalflux ! total flux

      io = iudtf(1)

*-----------------------------------------------------------------------
*     similar to [t-cross] tally
*     This tally returns the distances traveled in each history 
*     through each distinct material, as a large list.
*     If using a neutrino source, this allows drawing straight lines.
*     Minimally, just the history count, material ID number, 
*     and distance traveled from the previous point in cm is listed.
*     Full output with more columns can be printed if uncommented.
*-----------------------------------------------------------------------

      if( ncol .eq. 1 ) then ! Start of calculation
       totalflux=0.0
! minimized output (#history, matID, dist)
       write(io,'(''#history   matID    dist(cm)'')')
       return
! below is more full output  
!       write(io,'(''  #history    #batch      ncol   kf-code'',
!     & ''  cell_old  cell_new surfaceID     matID'',
!     & ''      E(MeV)       X(cm)       Y(cm)       Z(cm)'',
!     & ''      Weight  Cos(theta)    dist(cm)'')')
       return
      end if

!      if( ncol .eq. 2 ) then ! End of calculation
!       write(io,'(''Total Flux (1/cm^2/source)='',es12.4)')
!     & totalflux/rsouin
!       return
!      endif
      if( ncol .eq. 4 ) then ! source spawn
       xold=xc(no,ipomp+1)
       yold=yc(no,ipomp+1)
       zold=zc(no,ipomp+1)
       dist=0.0
      endif 
      if(ncol.eq.10) then ! timing for crossing surface
       erg=ec(no,ipomp+1)
!       lmat = 1  ! water should be always material 1 in this tally
!       if(ityp.eq.12.or.ityp.eq.13) then 
!        dedx=getdEdxH2O(erg) ! for electron and positron, determined from database
!       else
!        call dedxas(erg,dedx,lmat,ityp,ktyp,jtyp,rtyp)  ! dedx in MeV/cm
!       endif
!       dedx=dedx/10.0  ! LET in keV/um
       if(iblz1.ne.iblz2) then ! only update crossing between regions of differing materials
        dist=sqrt((xc(no,ipomp+1)-xold)**2+ (yc(no,ipomp+1)-yold)**2+
     &            (zc(no,ipomp+1)-zold)**2)
! minimized output (#history, matID, dist)
         write(io,'(2i8,1es12.4)') int(rcasc),idmn(mat),dist
! below is more full output     
!         write(io,'(8i10,20es12.4)') nocas,nobch,ncol,ktyp,iblz1,iblz2,
!     &   idsn(nsurf),idmn(mat),ec(no,ipomp+1),
!     &   xc(no,ipomp+1),yc(no,ipomp+1),zc(no,ipomp+1),oldwt,costha,dist
!         totalflux=totalflux+oldwt/costha
        xold=xc(no,ipomp+1)
        yold=yc(no,ipomp+1)
        zold=zc(no,ipomp+1)
       endif  
      endif

      return
      end
