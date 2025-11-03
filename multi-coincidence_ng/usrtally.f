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

      parameter(maxregion=200) ! maximum number of region to calculate deposition energy
      dimension depeventn(maxregion) ! deposition energy per neutron event
      save depeventn
!$OMP THREADPRIVATE(depeventn)
      dimension depeventg(maxregion) ! deposition energy per gamma event
      save depeventg
!$OMP THREADPRIVATE(depeventg)
      dimension xfirstcol(maxregion) ! x (cm) of energy deposit collision
      save xfirstcol
!$OMP THREADPRIVATE(xfirstcol)
      dimension yfirstcol(maxregion) ! y (cm) of energy deposit collision
      save yfirstcol
!$OMP THREADPRIVATE(yfirstcol)
      dimension zfirstcol(maxregion) ! z (cm) of energy deposit collision
      save zfirstcol
!$OMP THREADPRIVATE(zfirstcol)
      dimension tfirstcol(maxregion) ! t (ns) of energy deposit collision
      save tfirstcol
!$OMP THREADPRIVATE(tfirstcol)

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

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout
      real*8 edep
      common/EPCONTedep/ edep

      integer nudtvar
      common /cnudtvar/ nudtvar

      integer nobchold,nocasold,nameold,noold
      save nobchold,nocasold,nameold,noold
!$OMP THREADPRIVATE(nobchold,nocasold,nameold,noold)

      integer nctoneold,ncttwoold,nctonereq,ncttworeq,iwrtrxns
      save nctoneold,ncttwoold,nctonereq,ncttworeq,iwrtrxns
!$OMP THREADPRIVATE(nctoneold,ncttwoold,nctonereq,ncttworeq,iwrtrxns)

      integer nregnedep,nreggedep
      save nregnedep,nreggedep
!$OMP THREADPRIVATE(nregnedep,nreggedep)

      CHARACTER(LEN=10000) :: outstr
      save outstr
!$OMP THREADPRIVATE(outstr)
      CHARACTER(LEN=500) :: tmpstr,tmpstrb
      save tmpstr,tmpstrb
!$OMP THREADPRIVATE(tmpstr)
      character(:), allocatable :: outstrfin
      save outstrfin
!$OMP THREADPRIVATE(outstrfin)

      real*8 edepminn,edepming
      save edepminn,edepming
!$OMP THREADPRIVATE(edepminn,edepming)

      io = iudtf(1)

*-----------------------------------------------------------------------
*     similar to [t-deposit2] and [t-product] tallies
*
* The thought is to only score when counter 1 or counter 2 end a history
* with some final value or greater, recording the energy deposited in each 
* specified region in the udtvar(:) variables, recording the cell number,
* energy deposited, and (x,y,z) coordinates and time of each deposition.
*
* The initial use of this was to, in a <list-mode> fashion, list the 
* histories/events where multiple cells/detectors are hit in a single
* history, with the multiplicity requirements enforced with the counters,
* and also being able to distinguish between neutron and gamma-ray events.
* 
* The last six udtvar(n) variables have different meanings
* udtvar(nudtvar-5) = required chmax(1) / #regions w/ above threshold n Edep
* udtvar(nudtvar-4) = required chmax(2) / #regions w/ above threshold g Edep
* udtvar(nudtvar-3) = minimum energy deposit for neutrons to be recorded (MeV)
* udtvar(nudtvar-2) = minimum energy deposit for gammas to be recorded (MeV)
* udtvar(nudtvar-1) = toggle writing of reaction lines too (1=yes, 0=no)
* udtvar(nudtvar) = 1 / X, specifies if above are specified and used, 
*                         otherwise means they do not exist
* To be clear, the FINAL udtvar (n) specifies whether the five before it, 
*   (n-1) to (n-5) mean to control the max counter values, energy
*   deposition thresholds, and writing of extra reaction infor or are 
*   just  cell numbers (including this final variable).
* If this final udtvar is not equal to 1, then all udtvars are 
*   interpreted to be cell numbers.
* If equal to one, then the preceeding five are for the desired values of 
*   the specified counters and control for whether the reaction lines,
*   which describe secondary reaction products and incident particle eng,
*   are written in addition to the standard energy deposition lines.
*-----------------------------------------------------------------------

      if( ncol .eq. 2 .or. ncol .eq. 4 ) then ! Source generation or End of calculation
       if(nobch.eq.1.and.nocas.eq.1) then ! 1st history of 1st batch
!        write(io,'(''        #iomp    #batch  #history'',200i12)') 
!     &  ,(idnint(udtvar(ir)),ir=1,nudtvar)
        if(idnint(udtvar(nudtvar)).eq.1) then ! final utdvars are used for setting counter vals
         nctonereq = idnint(udtvar(nudtvar-5))
         ncttworeq = idnint(udtvar(nudtvar-4))
         edepminn = udtvar(nudtvar-3)
         edepming = udtvar(nudtvar-2)
         iwrtrxns = idnint(udtvar(nudtvar-1))
         nudtvar = nudtvar - 6
        else  
         nctonereq = 1
         ncttworeq = 1
         iwrtrxns = 1
         edepminn = 0.0
         edepming = 0.0
        endif  
        write(io,'(A,i8,A,i8)') "!Required final counter 1 value =", 
     &  nctonereq, " ; Required final counter 2 value =", ncttworeq
        write(io,'(A)') "!       #iomp    #batch  #history" //
     &   "       #no     #name"//
     &   "        #reg  EdepA(MeV)      xA(cm)      yA(cm)" // 
     &   "      zA(cm)      tA(ns)" // 
     &   "        #reg  EdepB(MeV)      xB(cm)      yB(cm)" // 
     &   "      zB(cm)      tB(ns)" //
     &   "        #reg  EdepC(MeV)      xC(cm)      yC(cm)" // 
     &   "      zC(cm)      tC(ns)" 
        write(io,'(A)') "!ncol   Z   N jcl kcl nclsts " 
        write(io,'(A)') "!In/Out kf-code     E(MeV)      weight " 
       else
!        if(maxval(depeventn(:)).gt.0.0) then 
        nregnedep = 0
        nreggedep = 0
        if(nctoneold.gt.0 .or. ncttwoold.ge.0) then ! only run calc if one of the counters is nonzero
         do ir=1,nudtvar
          if(depeventn(ir).ge.edepminn)then 
           nregnedep = nregnedep + 1
          else 
           depeventn(ir) = 0.0
          endif 
          if(depeventg(ir).ge.edepming)then 
           nreggedep = nreggedep + 1
          else 
           depeventg(ir) = 0.0
          endif 
         enddo 
!         nregnedep = count(depeventn(:nudtvar)>0.0)
!         nreggedep = count(depeventg(:nudtvar)>0.0)
        endif
        if((nregnedep.ge.nctonereq).or. ! nctoneold.ge.nctonereq.or.
     &     (nreggedep.ge.ncttworeq)) then ! ncttwoold.ge.ncttworeq.or.
!         write(io,'(A,3i5,4x,3i5)')"dbg",nctoneold,nregnedep,nctonereq,
!     &     ncttwoold,nreggedep,ncttworeq
!         if(nctoneold.ge.nctonereq) then ! if counter 1 (# neutron cols) >= 2
         if(nregnedep.ge.nctonereq) then ! if num regs exceeding energy limit >= required #
!!         write(io,'(''ne '',3i10,200es12.4)') ipomp,nobchold,nocasold,
!!     &   (depeventn(ir),ir=1,nudtvar)
          write(tmpstr,'(5i10)') ipomp,nobchold,nocasold,noold,nameold
          tmpstrb = "ne " // trim(tmpstr) // " ; "
          do ir=1,nudtvar
           if(depeventn(ir).gt.0.0)then 
            write(tmpstr,'(i10,10es12.4)')
     &       idnint(udtvar(ir)),depeventn(ir),
     &       xfirstcol(ir),yfirstcol(ir),zfirstcol(ir),tfirstcol(ir)
            tmpstrb = trim(tmpstrb) // trim(tmpstr) // " , "
           endif 
          enddo
          outstr = trim(tmpstrb) // trim(outstr)
         endif 
!         if(maxval(depeventg(:)).gt.0.0) then
!         if(ncttwoold.ge.ncttworeq) then ! if counter 2 (# gamma cols) >= 3
         if(nreggedep.ge.ncttworeq) then ! if num regs exceeding energy limit >= required #
!!          write(io,'(''ge '',3i10,200es12.4)') ipomp,nobchold,nocasold,
!!     &    (depeventg(ir),ir=1,nudtvar)
          write(tmpstr,'(5i10)') ipomp,nobchold,nocasold,noold,nameold
          tmpstrb = "ge " // trim(tmpstr) // " ; "
          do ir=1,nudtvar
           if(depeventg(ir).gt.0.0)then 
            write(tmpstr,'(i10,10es12.4)') 
     &       idnint(udtvar(ir)),depeventg(ir),
     &       xfirstcol(ir),yfirstcol(ir),zfirstcol(ir),tfirstcol(ir)
            tmpstrb = trim(tmpstrb) // trim(tmpstr) // " , "
           endif 
          enddo
          outstr = trim(tmpstrb) // trim(outstr) 
         endif
         outstrfin = CHAR(13)//CHAR(10)// trim(ADJUSTL(outstr))
         write(io,'(A)') outstrfin
        endif 
       endif
       depeventn(:)=0.0 ! initialization
       depeventg(:)=0.0 ! initialization
       xfirstcol(:)=0.0 ! initialization
       yfirstcol(:)=0.0 ! initialization
       zfirstcol(:)=0.0 ! initialization
       tfirstcol(:)=0.0 ! initialization
       outstr=""
       tmpstr=""
       tmpstrb=""
       nobchold=nobch
       nocasold=nocas
       if(ncol .ne. 2) then 
        noold=no
        nameold=name(no,ipomp+1)
       endif 
       nctoneold=0
       ncttwoold=0
       return
      endif

      if(ncol.ge.5) then ! timing for depositing energy
       do ir=1,nudtvar
        if(iblz1.eq.idnint(udtvar(ir))) then ! region check
         ein=e(no,ipomp+1)   ! current particle energy 
         ecc=ec(no,ipomp+1)  ! particle energy before this event
         if((ityp.eq.12.or.ityp.eq.13.or.ityp.eq.14).and.iegsemi.ne.0)
     &    then  ! EGS mode
          tlngth = sqrt( ( xc(no,ipomp+1) - x(no,ipomp+1) )**2
     &                 + ( yc(no,ipomp+1) - y(no,ipomp+1) )**2
     &                 + ( zc(no,ipomp+1) - z(no,ipomp+1) )**2 )
          call egs5edxde(ecc,tlngth,ein,mat,ityp)
          if(ncol.eq.9 .or. ncol.eq.11) then
           ecc = 0d0
           if( ityp.eq.12 .or. ityp.eq.13 ) ein = ein + edep
          endif
         else ! original PHITS mode
          ecc = ein - dedxfd * ( ein - ecc ) ! modification factor for electron
          if((ncol.eq.9.or.ncol.eq.11).and.jtyp.ne.0) ecc = 0.d0
         endif
         if(ein.ne.ecc) then
          edepo=ein-ecc ! energy changed
          nctoneold=ncnt(1,no,ipomp+1)
          ncttwoold=ncnt(2,no,ipomp+1)
          if(depeventg(ir).eq.0.0.and.depeventn(ir).eq.0.0) then ! first edep in reg
           xfirstcol(ir)=x(no,ipomp+1)
           yfirstcol(ir)=y(no,ipomp+1)
           zfirstcol(ir)=z(no,ipomp+1)
           tfirstcol(ir)=t(no,ipomp+1)
           if(noold.eq.1.and.nameold.eq.1) then ! overwrite initial value
            noold=no
            nameold=name(no,ipomp+1)
           endif 
          endif 
          if(ktyp.eq.22.or.abs(ktyp).eq.11) then ! gamma event
!           if(depeventg(ir).eq.0.0) then ! this is the first energy dep in cell
!            write(tmpstr,'(''gx '',2i10,2i4,8es12.4)') ktyp,iblz1,
!     &      ncnt(1,no,ipomp+1),ncnt(2,no,ipomp+1),
!     &      x(no,ipomp+1),y(no,ipomp+1),z(no,ipomp+1),t(no,ipomp+1)
!            outstr = trim(outstr) //CHAR(13)//CHAR(10)// trim(tmpstr)
!           endif
           depeventg(ir)=depeventg(ir)+edepo*oldwt
!            write(io,'(''g '',5i10,24es12.4)') nocas,nobch,ncol,ktyp,
!     &      iblz1,edepo,oldwt,e(no,ipomp+1),x(no,ipomp+1),y(no,ipomp+1),
!     &      z(no,ipomp+1),t(no,ipomp+1),ecc,
!     &      xc(no,ipomp+1),yc(no,ipomp+1),zc(no,ipomp+1),tc(no,ipomp+1)
          else ! neutron event
!           if(depeventn(ir).eq.0.0) then ! this is the first energy dep in cell
!            write(tmpstr,'(''nx '',2i10,2i4,8es12.4)') ktyp,iblz1,
!     &      ncnt(1,no,ipomp+1),ncnt(2,no,ipomp+1),
!     &      x(no,ipomp+1),y(no,ipomp+1),z(no,ipomp+1),t(no,ipomp+1)
!            outstr = trim(outstr) //CHAR(13)//CHAR(10)// trim(tmpstr)
!           endif 
           depeventn(ir)=depeventn(ir)+edepo*oldwt
!           write(io,'(''n '',5i10,24es12.4)') nocas,nobch,ncol,ktyp,
!     &     iblz1,edepo,oldwt,e(no,ipomp+1),x(no,ipomp+1),y(no,ipomp+1),
!     &     z(no,ipomp+1),t(no,ipomp+1),ecc,
!     &     xc(no,ipomp+1),yc(no,ipomp+1),zc(no,ipomp+1),tc(no,ipomp+1)
          endif 
!          write(io,'(5i10,20es12.4)') nocas,nobch,ncol,ktyp,iblz1,edepo,
!     &    oldwt,e(no,ipomp+1),x(no,ipomp+1),y(no,ipomp+1),z(no,ipomp+1),
!     &    ecc,xc(no,ipomp+1),yc(no,ipomp+1),zc(no,ipomp+1)
         endif
        endif
       enddo
      endif

!      if((ncol.eq.13.or.ncol.eq.14).and.mat.gt.0.and.nclsts.gt.0) then ! dead particle is produced by a reaction
      if(ncol.eq.13.or.ncol.eq.14) then ! reaction timing
       if(mat.gt.0.and.nclsts.gt.0) then ! dead particle is produced by a reaction
        do ir=1,nudtvar
         if(iblz1.eq.idnint(udtvar(ir))) then ! region check
          do j = 1, nclsts
           if( jclusts(4,j) .lt. 0 .and. jclusts(5,j) .ne. 0 ) then ! dead particle and !=0 charge
            ein = qclusts(7,j)
            ecc = 0.0d0
            edepo=ein
            kfd = jclusts(7,j)
            xd  = qclusts(10,j)
            yd  = qclusts(11,j)
            zd  = qclusts(12,j)
            td  = qclusts(9,j)
            nctoneold=ncnt(1,no,ipomp+1)
            ncttwoold=ncnt(2,no,ipomp+1)
            if(depeventg(ir).eq.0.0.and.depeventn(ir).eq.0.0) then ! first edep in reg
             xfirstcol(ir)=xd
             yfirstcol(ir)=yd
             zfirstcol(ir)=zd
             tfirstcol(ir)=td
             if(noold.eq.1.and.nameold.eq.1) then ! overwrite initial value
              noold=no
              nameold=name(no,ipomp+1)
             endif 
            endif
!            write(io,'(5i10,20es12.4)') nocas,nobch,ncol,kfd,iblz1,edepo,
!     &      oldwt,ein,xd,yd,zd,ecc,xd,yd,zd
            if(ktyp.eq.22.or.abs(ktyp).eq.11) then ! gamma event
!             if(depeventg(ir).eq.0.0) then ! this is the first energy dep in cell
!              write(tmpstr,'(''gx '',2i10,2i4,8es12.4)') kfd,iblz1,
!     &        ncnt(1,no,ipomp+1),ncnt(2,no,ipomp+1), 
!     &        xd,yd,zd,td
!              outstr = trim(outstr)//CHAR(13)//CHAR(10)//trim(tmpstr)
!             endif
             depeventg(ir)=depeventg(ir)+edepo*oldwt
!             write(io,'(''g '',5i10,24es12.4)') nocas,nobch,ncol,kfd,
!     &       iblz1,edepo,oldwt,ein,xd,yd,zd,td,ecc,xd,yd,zd,td
            else ! neutron event
!             if(depeventn(ir).eq.0.0) then ! this is the first energy dep in cell
!              write(tmpstr,'(''nx '',2i10,2i4,8es12.4)') kfd,iblz1,
!     &        ncnt(1,no,ipomp+1),ncnt(2,no,ipomp+1), 
!     &        xd,yd,zd,td
!              outstr = trim(outstr)//CHAR(13)//CHAR(10)//trim(tmpstr)
!             endif 
             depeventn(ir)=depeventn(ir)+edepo*oldwt
!             write(io,'(''n '',5i10,24es12.4)') nocas,nobch,ncol,kfd,
!     &       iblz1,edepo,oldwt,ein,xd,yd,zd,td,ecc,xd,yd,zd,td
            endif 
           endif
          enddo
         endif
        enddo
       endif 
       if(iwrtrxns.eq.1) then 
        if(jcoll.eq.6.or.jcoll.eq.8.or.jcoll.eq.9) return ! neutron, electron, proton reaction by nuclear or atomic data
        if(jcoll.eq.7.and.ipnint.eq.0) return ! photon reaction by atomic data
        do ir=1,nudtvar
         if(iblz1.eq.idnint(udtvar(ir)).and.abs(ktyp).ne.11) then ! region check, exclude electrons/positrons
!    Output reaction information
!         write(tmpstr,'(''  #history    #batch ncol   cell-ID'',
!     &   '' Z_mother N_mother jcoll kcoll nclsts'',
!     &   ''       x(cm)       y(cm)       z(cm)'')') 
!         outstr = trim(outstr) // CHAR(13)//CHAR(10) // trim(tmpstr)
          write(tmpstr,'(A2,i3,2i4,2i3,i4,3i10,2x,i10,12x,3es12.4)')  
     &    "rx",ncol,mathz,mathn,jcoll,kcoll,nclsts,nocas,
     &    no,name(no,ipomp+1),iblz1,
     &    xc(no,ipomp+1),yc(no,ipomp+1),zc(no,ipomp+1)
          outstr = trim(outstr) // CHAR(13)//CHAR(10) // trim(tmpstr)
          if(noold.eq.1.and.nameold.eq.1) then ! overwrite initial value
           noold=no
           nameold=name(no,ipomp+1)
          endif 
!    Output Primary particle information
!         write(tmpstr,'(''   Pri/Sec   kf-code      E(MeV)'',
!     &   ''      weight'')')  
!         outstr = trim(outstr) // CHAR(13)//CHAR(10) // trim(tmpstr)
          write(tmpstr,'(''  In'',i10,4es12.4)') 
     &    ktyp,ec(no,ipomp+1),oldwt
!     &    u(no,ipomp+1),v(no,ipomp+1),w(no,ipomp+1)
          outstr = trim(outstr) // CHAR(13)//CHAR(10) // trim(tmpstr)
!    Output Secondary particle information
          do j=1,nclsts
           px=qclusts(1,j)
           py=qclusts(2,j)
           pz=qclusts(3,j)
           pt=sqrt(px**2+py**2+pz**2)
           if(j.eq.1) then 
            write(tmpstr,'('' Out'')') 
            outstr = trim(outstr) //CHAR(13)//CHAR(10)// trim(tmpstr)
           endif 
           write(tmpstr,'(i10,8es12.4)') 
     &     jclusts(7,j),qclusts(7,j),qclusts(8,j) ! px/pt,py/pt,pz/pt = u, v, w
           if(j.eq.nclsts) then
            outstr = trim(outstr) // trim(tmpstr)
           else
            outstr = trim(outstr) // trim(tmpstr) // " , "
           endif 
!           totalproduct=totalproduct+qclusts(8,j)
          enddo
         endif
        enddo 
       endif
       
      endif

      return
      end
