MODULE inp_out_solv

! contains subroutines:
!
! read_parm
! read_sol_parm
! xtc2fortran_trj                !not used in versions > 7.0
! top_gro
! top_gro_sol

!... Solution Routines

! wat_bio_histo            - inplemented 2023
! wat_env_orient_tcfs
! wat_HB_tcfs              - not implemented
! wat_tetra
! wat_tetra_sol
! wat_LSI                  - development     
! prot_contact
! potent_map
! pot_map_umbrella
! prot_water_umbrella
! wat_env_diffusion
! wat_env_diff_stat        - discontinued
! vmd_trj

!... Pure Water Routines

! pure_wat_orient_tcfs
! pure_wat_LDL_HDL_otcf
! pure_wat_LDL_HDL_otcf_HB     - Lars HB strong/weak deconvolution
! pure_wat_HB_tcfs         
! pure_wat_HB_LDL_HDL_tcfs  
! pure_wat_LSI_tcfs            - development
! pure_wat_tetra
! pure_wat_diffus_msd
! pure_wat_LDL_HDL_msd
! pure_wat_diffus_msd_Env
! pure_wat_tetra_Env
! pure_wat_Voronoi

USE water_solv

implicit none

 contains
 
SUBROUTINE read_parm(ngrmx,ffname,inputformat,trjfile,unfoldtrjfile,unfold_trj,topfile,LJfile,ksLJ,boxfile,nsyst,nens,cube,nstep,dt,dfr,nmolsol,natmsol,&
                     nmolwat,nions,nwatsites,kbiowater,biowsamp,biowRadius,kwtcfcrsh,mtdelay,m_time,mtdelay_hb,m_time_hb,rwtcfbulk,khbtcf,hbtdelay,hbs_time,kwtetra,rwtetbulk,&
                     kwsoltetra,rwsoltetbulk,tetsampw,tetsampws,Voronsampw,khbpureW,hbtframe,hbs_samp,kscontact,intsf,rcont,kspotential,kumbrella_trj,kntrj,&
                     ncontsamp,kMIC,kMIC_,intsf_,rcont_,Coul_th_min,Coul_th_max,LJ126_th_min,LJ126_th_max,nattp,ncrossLJ,mcontsamp,ksprotwat,kprot_wat,kwat_wat,&
                     kdifwsol,msd_delay,difsampw,difbulk,krmvmd,vmd_samp,ktcfpureW,ktcfhbW,ktetpureW,kVoronoipureW,rcut,Nenvt,R5Tet,kdifpureW,QMMM_pureW,Nenv,Nenvrm,R5Dif,Nenv_hb,&
                     kLDL_HDL,kLDL_HDL_,LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2,Nenv_hb_HB,LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_,&
                     klsipureW,lsitframe,lsi_samp,Nenv_lsi,LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II,KWLSI,LSI_sample,rwlsibulk,WATMODEL)
!
    integer,intent(in)            :: ngrmx
!    
    character(len=7),intent(out)  :: inputformat ! type of MD input: TINKER, GROMACS, BOMD
    character(len=9),intent(out)  :: ffname      ! ffname: GROMOS, AMBER, CHARMM, AMOEBA, OPLS, BLYP etc
    character(len=15),intent(out) :: trjfile,topfile,LJfile,boxfile        ! GROMACS trj file (.xtc); GROMACS topology file (.top or .itp)
    character(len=14),intent(out) :: unfoldtrjfile                         ! internal - unfolded trajectory file name for Einstein diffusion trj_unfold.xtc
    integer,intent(out)           :: ksLJ                                  !0: do not read LJ1 2-6 parameters from file 1: read LJ 12-6 parm
    integer,intent(out)           :: nsyst       ! 0: pure water; 1: solute and water system
    integer,intent(out)           :: nens        ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(out)      :: cube(3)     ! box side length NVE and NVT; TINKER and BOMD type       
    integer,intent(out)           :: nstep       ! number of frames
    real(kind=8),intent(out)      :: dt          ! time-step (fs)
    integer,intent(out)           :: dfr         ! number of timesteps between frames  
    integer,intent(out)           :: nmolsol     ! number of molecules of solute     
    integer,intent(out)           :: natmsol     ! number of atoms per molecule of solute
    character(len=7),intent(out)  :: WATMODEL
    integer,intent(out)           :: nmolwat     ! number of water molecules
    integer,intent(out)           :: nions       ! number of ions
    integer,intent(out)           :: nwatsites   ! number of sites per water molecule (SPC/E = 3; TIP4P-EW = 4)
    integer,intent(out)           :: kbiowater   ! 0: do not calculate; 1: calculate bio-bulk water location
    integer,intent(out)           :: biowsamp    ! sampling time
    real(kind=8),intent(out)      :: biowRadius  ! bio-water cut-off
    integer,intent(out)           :: kwtcfcrsh   ! 0: do not calculate; 1: calculate - radial reor. and freq. tcfs
    integer,intent(out)           :: mtdelay     ! delay time
    integer,intent(out)           :: m_time      ! time freq. to sample water OH groups in radial shells for reor. tcf calculation
    integer,intent(out)           :: mtdelay_hb,m_time_hb    ! delay time; time freq. to sample water OH groups in radial shells for reor. tcf calculation
    real(kind=8)                  :: rwtcfbulk               ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk reorientational tcf
    integer,intent(out)           :: khbtcf                  ! 0: do not calculate; 1: calculate - hydrogen-bond tcfs (continuous and intermittent) 
    integer,intent(out)           :: hbtdelay,hbs_time       ! delay time and sampling time for HB tcfs calculation
    integer,intent(out)           :: kwtetra                 ! Flag for tetrahedrality calculation
    real(kind=8),intent(out)      :: rwtetbulk               ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk tetrahedrality
    integer,intent(out)           :: kwsoltetra              ! Flag for tetrahedrality calculation - include solute heavy atoms 
    real(kind=8),intent(out)      :: rwsoltetbulk            ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk tetrahedrality
    integer,intent(out)           :: KWLSI                   ! 0: do not calculate; 1: calculate - LSI for solutions    
    real(kind=8),intent(out)      :: rwlsibulk               ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk LSI
    integer,intent(out)           :: tetsampw,tetsampws      ! interval for samplig for the tetrahedrality - time-steps
    integer,intent(out)           :: Voronsampw              ! interval for samplig for Voronoi - time-steps
    real(kind=4),intent(out)      :: rcut
    integer,intent(out)           :: LSI_sample              ! LSI sampling interval - time-steps
    integer,intent(out)           :: kscontact,intsf         ! Flag for contact map between proteins; contact map trajectory sample frequency (steps) [routine prot_contact]
    real(kind=8),intent(out)      :: rcont                   ! distance cut-off for contact map - residues at a r < rcont are accounted [routine prot_contact]
    integer,intent(out)           :: kMIC,kMIC_              ! apply Minimum Image Convention or not to proteins
    integer,intent(out)           :: intsf_                  ! contact map trajectory sample frequency (steps) [routine potent_map]
    real(kind=8),intent(out)      :: rcont_                  ! distance cut-off for contact map - residues at a r < rcont are accounted [routine potent_map]
    integer,intent(out)           :: kspotential,ncontsamp   ! Flag for contact map electrostatics and van der Waals; sampling frequency
    integer,intent(out)           :: mcontsamp               ! sampling frequency protein-water and water-water
    integer,intent(out)           :: kumbrella_trj           ! 0: Trjs are not from umbrella sampling 1: Trjs are from umbrella sampling [no averages over trjs will be calculated]
    integer,intent(out)           :: kntrj                   ! Number of MD trajectories to build a mean electrostatic contact map
!V15    real(kind=8),intent(out)      :: Coul_th_min,LJ126_th_min        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!V15    real(kind=8),intent(out)      :: Coul_th_max,LJ126_th_max        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!Version 15 - High energy RES-RES contacts will be chosen based on a single criterion so that Coulombic and van der Waals pairs are the same
    real(kind=8),intent(out)      :: Coul_th_min,Coul_th_max ! Potential energy threshold for printing res-res interactions [Min,Max]
    real(kind=8),intent(out)      :: LJ126_th_min,LJ126_th_max     ! Lennard-Jones 12-6 thresholds for printing res-res interactions
    integer,intent(out)           :: nattp                   ! Number of distinct atomtypes in the force field file (ffnonbonded.itp - Lennard-Jones parm)
    integer,intent(out)           :: ksprotwat,kprot_wat,kwat_wat !protein-water options
    integer,intent(out)           :: ncrossLJ                ! Number of LJ126 crossed IJ terms - for reading non-geometric C12
    integer,intent(out)           :: kdifwsol,difsampw       ! Flag for diffusion calculation - msd Einstein
    real(kind=8)                  :: difbulk                 ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk diffusion
    integer,intent(out)           :: msd_delay
    integer,intent(out)           :: krmvmd,vmd_samp
!Pure water    
!    integer,intent(out)           :: nbwmol,nwatpoints      ! number of water molecules, number of water sites - pure water
!    real(kind=8),intent(out)      :: wside                  ! box side length (Ang) - pure water
    integer,intent(out)           :: ktcfpureW               ! 0: do not calculate; 1: calculate reorientational tcfs - pure water
    integer,intent(out)           :: ktcfhbW                 ! 0: do not calculate; 1: calculate reorientational tcfs - pure water - use HB definition for OH tcf deconvolution
!    integer,intent(out)           :: tcfframe,tcf_samp      ! delay time and sampling time for reorientational tcf calculation - pure water                                 
    integer,intent(out)           :: khbpureW                ! 0: do not calculate; 1: calculate - hydrogen-bond tcfs (continuous and intermittent) - pure water
    integer,intent(out)           :: hbtframe,hbs_samp       ! delay time and sampling time for HB tcfs calculation - pure water
    integer,intent(out)           :: ktetpureW               ! 0: do not calculate; 1: calculate - tetrahedrality of pure water
    integer,intent(out)           :: kVoronoipureW           ! 0: do not calculate; 1: calculate - Voronoi polyhedra
    real(kind=8),intent(out)      :: R5Tet                   ! cut-off to study tetrahedrality in different environments in pure water
    integer,intent(out)           :: kdifpureW               ! 0: do not calculate; 1: calculate - diffusion of pure water
    integer,intent(out)           :: unfold_trj              ! 0: unfold trj for msd calculation; 1: do not unfold trajectory
    integer,intent(out)           :: QMMM_pureW              ! 0: do not print; 1: print QMMM conf. for pure water
    integer,intent(out)           :: Nenv,Nenvrm,Nenvt       ! Number of environments to calculate the orientational tcfs/diffusion/tetrahedrality for pure water (LDL/HDL)
    real(kind=8),intent(out)      :: R5Dif                   ! cut-off to study diffusion in different environments in pure water
!    integer,intent(out)           :: nbulk                   ! 0 skip tetrahedrality calculation for bulk water; 1: include bulk tetrahedrality
    integer,intent(out)           :: Nenv_hb,Nenv_hb_HB      ! Number of environments to calculate the orientational tcfs - [HB deconvolution version 24/25]
    integer,intent(out)           :: kLDL_HDL,kLDL_HDL_      ! 0:read LDL/HDL list; 1:calculate LDL/HDL on the fly
    real(kind=8),intent(out)      :: LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2
    real(kind=8),intent(out)      :: LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_ 
    integer,intent(out)           :: klsipureW                ! 0: do not calculate; 1: calculate - LSI tcfs (continuous and intermittent) - pure water
    integer,intent(out)           :: lsitframe,lsi_samp,Nenv_lsi
    real(kind=8),intent(out)      :: LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II
!    
    unfoldtrjfile='trj_unfold.xtc'
    read(ngrmx,*)inputformat
    read(ngrmx,*)ffname
    read(ngrmx,*)trjfile
    read(ngrmx,*)unfold_trj
    read(ngrmx,*)topfile
    read(ngrmx,*)LJfile
    read(ngrmx,*)nattp,ncrossLJ
    read(ngrmx,*)boxfile
    read(ngrmx,*)nsyst
    read(ngrmx,*)nens
    read(ngrmx,*)cube(1),cube(2),cube(3)
    read(ngrmx,*)nstep
    read(ngrmx,*)dt
    read(ngrmx,*)dfr
    read(ngrmx,*)                                      !input comment line
    read(ngrmx,*)nmolsol,natmsol
    read(ngrmx,*)WATMODEL
    read(ngrmx,*)nmolwat,nwatsites
    read(ngrmx,*)nions
    read(ngrmx,*)kbiowater,biowsamp,biowRadius
    read(ngrmx,*)kwtcfcrsh,mtdelay,m_time,rwtcfbulk
    read(ngrmx,*)khbtcf,hbtdelay,hbs_time
    read(ngrmx,*)kwtetra,tetsampw,rwtetbulk
    read(ngrmx,*)kwsoltetra,tetsampws,rwsoltetbulk
    read(ngrmx,*)KWLSI,LSI_sample,rwlsibulk
    read(ngrmx,*)kscontact,intsf,rcont
    read(ngrmx,*)kspotential,kumbrella_trj,kntrj,ncontsamp
    read(ngrmx,*)intsf_,rcont_,kMIC
!    read(ngrmx,*)Coul_th_min,Coul_th_max,LJ126_th_min,LJ126_th_max
    read(ngrmx,*)Coul_th_min,Coul_th_max
    read(ngrmx,*)LJ126_th_min,LJ126_th_max
    read(ngrmx,*)ksprotwat,kMIC_,mcontsamp
    read(ngrmx,*)kprot_wat,kwat_wat
    read(ngrmx,*)kdifwsol,msd_delay,difsampw,difbulk
    read(ngrmx,*)krmvmd,vmd_samp
        
    if(nsyst==0)then                                     !pure water
       read(ngrmx,*)                                     !input comment line
       read(ngrmx,*)nmolwat,nwatsites 
       read(ngrmx,*)cube(1),cube(2),cube(3)
!Version 24 Nenv - number of environments - to calculate the otcf for bulk, LDL, HDL, and other       
       read(ngrmx,*)ktcfpureW,mtdelay,m_time,Nenv,kLDL_HDL
       read(ngrmx,*)LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2
       read(ngrmx,*)ktcfhbW,mtdelay_hb,m_time_hb,Nenv_hb,kLDL_HDL_
       read(ngrmx,*)khbpureW,hbtframe,hbs_samp,Nenv_hb_HB
       read(ngrmx,*)LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_
       read(ngrmx,*)klsipureW,lsitframe,lsi_samp,Nenv_lsi
       read(ngrmx,*)LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II
       read(ngrmx,*)ktetpureW,tetsampw,Nenvt,R5Tet
       read(ngrmx,*)kVoronoipureW,Voronsampw,rcut
       read(ngrmx,*)kdifpureW,msd_delay,difsampw,Nenvrm,R5Dif
       read(ngrmx,*)QMMM_pureW
    endif
    
!Version_12   
! FIX_ME: The program cannot open multiple .xtc trj files if it has already ran another routine...
     if(kspotential==1.and.kntrj>1)then
        write(*,*)
        write(*,*)'WARNING - POTENTIAL MAP with multiple trjs TURNED ON'
        write(*,*)'No additional analysis will be performed...'
        write(*,*)
        kwtcfcrsh = 0
        khbtcf = 0
        kwtetra = 0
        kwsoltetra = 0
        kscontact = 0
      endif  
!      
      ksLJ = 0
      if(kspotential==1)ksLJ = 1
      if(ksprotwat==1)ksLJ = 1

      !    if(kspotential==1.and.kscontact==0)then
!       write(*,*)
!       write(*,*)'WARNING - POTENTIAL MAP is on while CONTACT MAP is off'
!       write(*,*)
!       write(*,*)'CONTACT MAP will be turned on'
!       write(*,*)'CONTACT MAP sampling frequency (steps) = ',intsf
!       write(*,*)'CONTACT MAP distance cut-off (Ang) = ',rcont
!       write(*,*)
!       kscontact = 1
!    endif   
! End Version 12
    
       
    
    close(ngrmx)
     
END subroutine read_parm



SUBROUTINE read_sol_parm(nsolprm,nsoltypes,Zattype,mattype,niontypes,niontype,miontype,Ziontype,nchains,&
                         nchain_nat,nchain_lab,nb_res_ch,chain_label_res,natHSh,natsol_id,oatsol_hs,ratsol_hs,nb_res)
!
    integer,intent(in)                                        :: nsolprm
!    
    integer,intent(out)                                       :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(:),allocatable,intent(out)     :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(:),allocatable,intent(out)         :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(out)                                       :: niontypes   ! Number of ions distinct species (Na+, Cl- etc)
    integer,dimension(:),allocatable,intent(out)              :: niontype    ! Number of ions of type (Na+, Cl- etc)
    integer,intent(out)                                       :: nchains     ! Number of protein chains (e.g., HgbS = 8)
    integer,dimension(:),allocatable,intent(out)              :: nchain_nat  ! Number of atoms in each chain
    character(len=1),dimension(:),allocatable,intent(out)     :: nchain_lab  ! Chain label (e.g., A, B, C etc) read
    character(len=1),dimension(:),allocatable,intent(out)     :: chain_label_res ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(:),allocatable,intent(out)              :: nb_res_ch   ! number of residues of each chain
    real(kind=8),dimension(:),allocatable,intent(out)         :: miontype    ! Mass of ions of type (Na+, Cl- etc)
    character(len=2),dimension(:),allocatable,intent(out)     :: Ziontype    ! Atomic symbol of ion types (Na+, Cl- etc)
    integer,intent(out)                                       :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(:),allocatable,intent(out)              :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
    real(kind=8),dimension(:),allocatable,intent(out)         :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(:),allocatable,intent(out)         :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    integer,intent(out)                                       :: nb_res      ! Number of residues [GROMACS]
!   Local
    integer                                                   :: i,j,k
!
    read(nsolprm,*)                     !comment line
    read(nsolprm,*)nb_res               !number of residues
!atom types    
    read(nsolprm,*)nsoltypes
    allocate(Zattype(nsoltypes),mattype(nsoltypes))
    do i = 1,nsoltypes
       read(nsolprm,*)Zattype(i),mattype(i)
    end do   
!ions    
    read(nsolprm,*)niontypes
    allocate(Ziontype(niontypes),niontype(niontypes),miontype(niontypes))
    do i = 1,niontypes
       read(nsolprm,*)Ziontype(i),niontype(i),miontype(i)
    end do 
!number of protein chains
    read(nsolprm,*)nchains
    allocate(nchain_nat(nchains),nchain_lab(nchains),nb_res_ch(nchains),chain_label_res(nb_res))
    do i = 1,nchains
       read(nsolprm,*)nchain_nat(i),nchain_lab(i),nb_res_ch(i)
    end do
    k = 1                                                !residue counter
    do i = 1,nchains
       do j=1,nb_res_ch(i)
          chain_label_res(k) = nchain_lab(i)
!check          write(*,*)k,chain_label_res(k)            
          k = k + 1
       end do  
    end do   
!solute atoms solvation    
    read(nsolprm,*)natHSh 
    allocate(natsol_id(natHSh),ratsol_hs(natHSh),oatsol_hs(natHSh))
    do i = 1,natHSh
       read(nsolprm,*)natsol_id(i),oatsol_hs(i),ratsol_hs(i) 
    end do 
 
    close(nsolprm)
    
END subroutine read_sol_parm
 
 
 
SUBROUTINE xtc2fortran_trj(ntrjunf,ncell,natms,nstep,trjfile)

    character(len=15),intent(in)          :: trjfile
    integer,intent(in)                    :: ntrjunf,ncell
    integer,intent(in)                    :: natms 
    integer,intent(in)                    :: nstep
 
!    real(kind=8),intent(out)       :: x(natms),y(natms),z(natms) 
!    real(kind=4),intent(out)       :: xx(natms),yy(natms),zz(natms) 
!    real(kind=8),intent(out)       :: cell(3)

!    vmd plugin variables

!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4 npart, maxatom, handle(4), status
    real(kind=4)                          :: box(6)
    real(kind=4)                          :: xyz(natms*3) 
    character infile*200, intype*10
    integer i,j

!    infile = 'gro.xtc'
    infile = TRIM(trjfile)
    intype = 'auto'
    npart  = -1
    handle(1) = -1
    handle(2) = -1
    handle(3) = -1
    handle(4) = -1
    
!    print*,'filename: ', infile
!    print*,'type:     ', intype      

!     set up everything and register all static plugins
!NG    call f77_molfile_init
    call f77_molfile_open_read(handle(1),npart,infile,intype)
    if (handle(1).lt.0) then
       print*,'file type unknown or not registered'
       stop
    else
       print*,'file successfully opened:'
       print*,'handle:',handle(1)
       print*,'npart: ',npart
       print*,'nsteps:',nstep
    end if
      
    if(npart /= natms)then
      write(*,*)'error - number of atoms is wrong in solv_inp.dat'
      stop
    endif
    do i=1,nstep
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
!      print*,'read ',i,'  status:',status
!      print*,'atom(1)', (xyz(j),j=1,3)
!      print*,'atom(10)',(xyz(j),j=31,33)
!      print*,'atom(100)',(xyz(j),j=301,303)
!      print*,'box',box

    write(ncell,*)box(1),box(2),box(3)
    write(ntrjunf)dble(box(1)),dble(box(2)),dble(box(3))
    do j = 1,npart*3,3
!check       write(*,'(3(f14.5))')xyz(j),xyz(j+1),xyz(j+2)                   
!       write(ntrjunf)dble(xyz(j)),dble(xyz(j+1)),dble(xyz(j+2))        !double precision
       write(ntrjunf)xyz(j),xyz(j+1),xyz(j+2)                           !single precision
!       write(*,*)xyz(j),xyz(j+1),xyz(j+2)     
    end do
    if(status.eq.0) then
       write(*,*)'error on reading trajectory file - step',i
       stop
    endif
    end do
    
!    write(*,*)'MD box side length in Ang:', box(1),box(2),box(3)
 
    call f77_molfile_finish
!    close(ntrjunf)

    write(*,*)'finished writing fortran unformatted trj...'

    return

END SUBROUTINE xtc2fortran_trj



SUBROUTINE top_gro(ngroitp,ffname,nLJitp,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,niontypes,niontype,nb_res,nb_res_ch,nchains,segname,resindex,resname,&
                   chrg,atomname,atomname_d,atomtype,nsoltypes,Zattype,mattype,ZMBIO,Ziontype,chain_label_res,Natresid,name_res,numb_res,ksLJ,nattp,ncrossLJ,&
                   C6kjmol,C12kjmol,C6_IJ,C12_IJ,C6_IW,C12_IW)

! System topology for GROMACS
    integer,intent(in)  :: ngroitp
    character(len=9),intent(in)  :: ffname                           ! ffname: GROMOS, AMBER, CHARMM, AMOEBA, OPLS, BLYP etc
    integer,intent(in)  :: nLJitp
    integer,intent(in)  :: natms
    integer,intent(in)  :: nmolsol
    integer,intent(in)  :: natmsol
    integer,intent(in)  :: nmolwat,nwatsites,nions
    integer,intent(in)                                :: nb_res
    integer,intent(in)                                :: nchains     ! Number of protein chains (e.g., HgbS = 8)
    integer,intent(in)                                :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in)  :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)      :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                                :: niontypes   ! Number of ions distinct species (Na+, Cl- etc)
    integer,dimension(niontypes),intent(in)           :: niontype    ! Number of ions of type (Na+, Cl- etc)
    character(len=2),dimension(niontypes),intent(in)  :: Ziontype    ! Atomic symbol of ion types (Na+, Cl- etc)
    character(len=1),dimension(nb_res),intent(in)     :: chain_label_res  ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nchains),intent(in)             :: nb_res_ch        ! number of residues of each chain
    integer,intent(in)                                :: ksLJ             !0: do not read LJ1 2-6 parameters from file 1: read LJ 12-6 parm
    integer,intent(in)                                :: nattp            !number of distinct atomtypes in the force field file (ffnonbonded.itp - Lennard-Jones parm)
    integer,intent(in)                                :: ncrossLJ         !number of LJ126 crossed IJ terms - for reading non-geometric C12
                                                      
    character(len=3),dimension(natms),intent(out)     :: segname
    integer,dimension(natms),intent(out)              :: resindex 
    character(len=4),dimension(natms),intent(out)     :: resname
    real(kind=4),dimension(natms),intent(out)         :: chrg
    character(len=4),dimension(natms),intent(out)     :: atomname
    character(len=4),dimension(natms),intent(out)     :: atomname_d
    character(len=4),dimension(natms),intent(out)     :: atomtype 
    real(kind=8),dimension(:),allocatable,intent(out) :: ZMBIO
    integer,dimension(:),allocatable,intent(out)      :: Natresid                           ! Number of atoms of each residue
    character(len=4),dimension(:),allocatable,intent(out)         :: name_res
    integer,dimension(:),allocatable,intent(out)      :: numb_res                           ! Number of the residue: from 1 to the total number of residues of each chain
    real(kind=8),dimension(:),allocatable,intent(out) :: C6kjmol,C12kjmol                   ! Solute Lennard-Jones parameters
    real(kind=8),dimension(:,:),allocatable,intent(out)           :: C6_IJ,C12_IJ           ! Solute Lennard-Jones C6 and C12 crossed parameters
    real(kind=8),dimension(:),allocatable,intent(out)             :: C6_IW,C12_IW           ! Solute-water (oxygen) Lennard-Jones C6 and C12 crossed parameters 
    
    
!Local
    character(len=4),dimension(natms)       :: atomtype_gro
    character(len=4),dimension(natms)       :: atomname_gro
    character(len=4),dimension(natms)       :: resname_gro
    
    character(len=1)                        :: res_label
    character(len=1)                        :: label_res
    character(len=1),dimension(nsoltypes)   :: Zattype_1
    character(len=1),dimension(natms)       :: atomname_1
    
    character(len=7)                        :: resid_label
    character(len=7)                        :: label_resid
    
    integer :: i,j,k,iw,kw,res_i,iat,nit,is
    integer :: cgnr                                                   !Charge Group NumbeR
    
    character(len=5),dimension(natmsol)     :: atomtype_LJ            ! input atomtype [read]
    character(len=4),dimension(natmsol)     :: atomtype_LJ_04         !atom type reduced to length 4 to compare with atomtype [length 4]
    real(kind=4)                            :: sm,sc
    character(len=1)                        :: sp
    real(kind=8),dimension(:),allocatable   :: C6LJ,C12LJ             ! Solute Lennard-Jones parameters
    integer                                 :: nat_LJ_01
    
    integer                                 :: k1,k2,ic1,ic2,iats1,iats2
    
    character(len=5),dimension(ncrossLJ)     :: atomtype_LJ_i            ! input atomtype [read]
    character(len=5),dimension(ncrossLJ)     :: atomtype_LJ_j            ! input atomtype [read]
    character(len=4),dimension(ncrossLJ)     :: atomtype_LJ_i_04         !atom type reduced to length 4 to compare with atomtype [length 4]
    character(len=4),dimension(ncrossLJ)     :: atomtype_LJ_j_04         !atom type reduced to length 4 to compare with atomtype [length 4]
    integer                                  :: nf
    integer,dimension(:,:),allocatable       :: nC12rep
    integer,dimension(:),allocatable         :: nC12rep_W    
    real(kind=8)                             :: C6WOkjmol,C12WOkjmol
    

allocate(ZMBIO(natmsol)) 
allocate(Natresid(nb_res))
allocate(name_res(nb_res))
allocate(numb_res(nb_res))
allocate(C6LJ(natmsol),C12LJ(natmsol))
allocate(C6kjmol(natmsol),C12kjmol(natmsol))
allocate(C6_IJ(natmsol,natmsol),C12_IJ(natmsol,natmsol))
allocate(C6_IW(natmsol),C12_IW(natmsol))
allocate(nC12rep(natmsol,natmsol))
allocate(nC12rep_W(natmsol))

ZMBIO = 0.0d0
Natresid = 0
numb_res = 0
   C6LJ = 0.0d0
   C12LJ = 0.0d0
   C6kjmol = 0.0d0
   C12kjmol = 0.0d0
nat_LJ_01 = 0 

res_label = ';'
res_i = 0
iat = 0

resid_label='residue'
label_resid='frkafka'

do i=1,nsoltypes 
   Zattype_1(i) = trim(adjustl(Zattype(i)))
!NG**   write(*,*)Zattype_1(i)
end do

write(*,*)'Reading topology...'
    
!read GROMACS top header - the number of lines changes with the system - use the label "residue" to find the begining of topology     
!do j = 1,26           
!do j = 1,25
!   read(ngroitp,'()')
!end do

do while(label_resid/=resid_label)
   read(ngroitp,'(a,1x,a7)')label_res,label_resid
end do 
backspace(ngroitp)                 !Rewind a single line

!read GROMACS topology file
do j = 1,nmolsol*natmsol + nb_res    
   read(ngroitp,'(a)')label_res
   if(label_res==res_label)then
      res_i = res_i + 1
!NG**      write(*,*)'Reading residue',res_i
   else
      iat = iat + 1
      backspace(ngroitp)                 !Rewind a single line
!Version 10 contact map
      Natresid(res_i) = Natresid(res_i) + 1 
!End version 10
!      read(ngroitp,'(1x,i5,7x,a4,4x,i5,2x,a4,2x,a4)')i,atomtype_gro(iat),resindex(iat),resname_gro(iat),atomname_gro(iat)
      read(ngroitp,'(1x,i5,7x,a4,4x,i5,2x,a4,2x,a4,3x,i4,4x,f7.3)')is,atomtype_gro(iat),resindex(iat),resname_gro(iat),atomname_gro(iat),&
      cgnr,chrg(iat)
      write(atomtype(iat),'(a)') trim(adjustl(atomtype_gro(iat)))
      write(resname(iat) ,'(a)') trim(adjustl(resname_gro(iat)))
      write(atomname(iat),'(a)') trim(adjustl(atomname_gro(iat)))  
      write(atomname_d(iat),'(a)') trim(adjustl(atomname_gro(iat))) 
      segname(iat) = 'MOL'
      atomname_1(iat) = trim(adjustl(atomname(iat)))                                                 !re-define atomic symbol by a single character 
!NG      write(*,500)iat,segname(iat),resindex(iat),resname(iat),atomname(iat),atomtype(iat)
!      write(*,*)iat,atomtype(iat),chrg(iat)
      do i=1,nsoltypes                                              !Assign atomic masses
         if(atomname_1(iat)==Zattype_1(i)) ZMBIO(iat)=mattype(i)
      end do
      name_res(res_i) = resname(iat)
!check      
!      write(*,*) iat,atomname(iat),ZMBIO(iat)                   
   endif   
end do

if(res_i/=nb_res)then
   write(*,*)'Error - number of residues in gro.itp different from input'
   write(*,*)res_i,' /= ',nb_res
   stop
endif
   
!Version 10
!
k = 0
do i=1,nchains
   do j=1,nb_res_ch(i)
      k=k+1
      numb_res(k) = j
   end do
end do

!check do j=1,nb_res                                                                       !Loop over total number of residues
!check    write(*,*)j,Natresid(j),name_res(j),numb_res(j),chain_label_res(j)
!check end do
!End version 10
  
do j = 1,nmolsol*natmsol
   if(ZMBIO(j)==0.0)then
      write(*,*)'Error!!! atomic mass undefined'
      write(*,*)'index = ',j,' atom name = ',atomname_1(j)
   endif
end do   

write(*,*)
write(*,*)'Ended reading topology'

!Water topology
iw = nb_res + 1
kw = nmolsol*natmsol + 1 
do j = 1,nmolwat
   segname(kw)  = 'WAT' 
   resindex(kw) =  iw 
   if(nwatsites==3)resname(kw)  = 'SPCE'
   if(nwatsites==4)resname(kw)  = 'TIP4'
   atomname(kw) = 'OH2' 
   atomtype(kw) = 'OW'
!NG**   write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
   kw = kw + 1
   segname(kw)  = 'WAT' 
   resindex(kw) =  iw 
   if(nwatsites==3)resname(kw)  = 'SPCE'
   if(nwatsites==4)resname(kw)  = 'TIP4'  
   atomname(kw) = 'H1' 
   atomtype(kw) = 'HW'
!NG**   write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
   kw = kw + 1
   segname(kw)  = 'WAT' 
   resindex(kw) =  iw 
   if(nwatsites==3)resname(kw)  = 'SPCE'
   if(nwatsites==4)resname(kw)  = 'TIP4'   
   atomname(kw) = 'H2' 
   atomtype(kw) = 'HW'
!NG**   write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
   kw = kw + 1
   if(nwatsites==4)then
      segname(kw)  = 'WAT' 
      resindex(kw) =  iw 
!      if(nwatsites==3)resname(kw)  = 'SPCE'
!      if(nwatsites==4)resname(kw)  = 'TIP4'
      resname(kw)  = 'TIP4'
      atomname(kw) = 'M' 
      atomtype(kw) = 'MW'
!NG**      write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
      kw = kw + 1
   endif 
end do


if(nions > 0)then
   iw = nb_res + 2
   do i = 1,niontypes
      nit = niontype(i) 
      do j = 1,nit
         segname(kw)  = 'ION' 
         resindex(kw) =  iw
         resname(kw)  = 'ION'
         write(atomname(kw),*)Ziontype(i) 
         write(atomtype(kw),*)Ziontype(i)
   !NG**      write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
         kw = kw + 1
      end do
      iw = iw + 1
   end do   
endif 

write(*,*)
write(*,*)'Total number of atoms = ',natms
write(*,*)



!Read Lennard-Jones 12-6 parameters - to be used in prot_cont_map calculation

! PROBLEMS: GROMOS54A7 does not use geometric combining rules for all LJ C12 coefficients - dx.doi.org/10.1002/jcc.20090 [J Comput Chem 25: 1656–1676, 2004]
! SOLUTION: First apply geometric rule to all C6 and C12 parameters and store them in a matrix. Second read [ nonbond_params ] from ffnonbonded.itp
!           and replace the C12 parameters which are different from those obtained through the geometric rule.

if(ksLJ==1)then
!   read(nLJitp,*)nattp            [input parameter]
   read(nLJitp,*)
   read(nLJitp,*)
   write(*,*)'Reading Lennard-Jones parameters...'
   write(*,*)'Number of Lennard-Jones atom types = ',nattp
   do i=1,nattp
      read(nLJitp,'(a5,3x,i2,2(4x,f7.3)5x,a1,2(2x,f12.10))')atomtype_LJ(i),is,sm,sc,sp,C6LJ(i),C12LJ(i)
      atomtype_LJ_04(i) = trim(adjustl(atomtype_LJ(i)))
!check      write(*,*)atomtype_LJ(i),atomtype_LJ_04(i),is,sm,sc,sp,C6LJ(i),C12LJ(i)
   end do
   do j = 1,nmolsol*natmsol
!    do j = 1,1                         !testing - 1 atom
      do i=1,nattp 
!check         write(*,*)atomtype(j),atomtype_LJ_04(i)
         if(atomtype(j)==atomtype_LJ_04(i))then
            nat_LJ_01 = nat_LJ_01 + 1
            C6kjmol(j)= C6LJ(i)
            C12kjmol(j)= C12LJ(i)
!check            write(*,*)'Found LJ atom type...'
            goto 5
         endif
      end do
      5 continue
   end do
   
   write(*,*)'Found Lennard-Jones parameters for ',nat_LJ_01,' atoms out of ',nmolsol*natmsol
   write(*,*)
   if(nat_LJ_01.ne.nmolsol*natmsol)then
      write(*,*)'Some solute Lennard-Jones parameters were not found'
      stop
   endif   
    
write(*,*)

!Version22
!SPC/E model - oxygen
    C6WOkjmol  = 0.0026173456d0 
    C12WOkjmol = 2.634129d-6 
!End Version 22

! Calculate interaction IJ Lennard-Jones parameters from geometric rules for C6 and C12 

   k1 = 0                               !Chains ABCD atomic index
   do ic1=1,nb_res/2                    !residues of ABCD
      do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
         k1 = k1 + 1
         k2 = natmsol/2                 !Chains EFGH atomic index                  
         do ic2=nb_res/2+1,nb_res       !residues of EFGH
            do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
               k2 = k2 + 1
               C6_IJ(k1,k2) = dsqrt(C6kjmol(k1)*C6kjmol(k2))                  !kJ/mol nm**6
               C12_IJ(k1,k2) = dsqrt(C12kjmol(k1)*C12kjmol(k2))                  !kJ/mol nm**12
             end do
         end do   
      end do
   end do
   
!Version_22
!Protein-water LJ parameters from geometric rules for C6 and C12
   k1 = 0                               !Chains ABCD and EFGH atomic index
   do ic1=1,nb_res                      !residues of ABCD and EFGH
      do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
         k1 = k1 + 1
         C6_IW(k1)  =  dsqrt(C6kjmol(k1)*C6WOkjmol)                    !kJ/mol nm**6
         C12_IW(k1) =  dsqrt(C12kjmol(k1)*C12WOkjmol)                  !kJ/mol nm**12
      end do
   end do   
!End Version_22
   
! Replace C6 and C12 IJ Lennard-Jones parameters for some atom pairs according to [ nonbond_params ] from ffnonbonded.itp 
! For C6 this is not necessary as they always obey the geometric rule - values are changed nevertheless to use the exact same values in ffnonbonded.itp
   if(ffname=='GROMOS   ')then 
      write(*,'(a71)')'Force Field is GROMOS - read modified cross C12(i,j) parms from input'
      write(*,*)
      read(nLJitp,*)
      read(nLJitp,*)
      read(nLJitp,*)
      write(*,*)'Reading Lennard-Jones C6 and C12 ij [nonbond_params] ...'
      write(*,*)'Number of Lennard-Jones input nonbond_params = ',ncrossLJ
      write(*,*)
      do i=1,ncrossLJ
         read(nLJitp,'(a5,3x,a5,3x,i4,4x,e12.6,4x,e12.6)')atomtype_LJ_i(i),atomtype_LJ_j(i),nf,C6LJ(i),C12LJ(i)
         atomtype_LJ_i_04(i) = trim(adjustl(atomtype_LJ_i(i)))
         atomtype_LJ_j_04(i) = trim(adjustl(atomtype_LJ_j(i)))
!check      write(*,'(i6,1x,a4,1x,a4,3x,e12.6,3x,e12.6)')i,atomtype_LJ_i_04(i),atomtype_LJ_j_04(i),C6LJ(i),C12LJ(i)
      end do
!Replace parameters in the respective arrays  
      write(*,*)'Start checking C6 and C12 ij [nonbond_params] ...'
      nC12rep = 0
      k1 = 0                               !Chains ABCD atomic index
      do ic1=1,nb_res/2                    !residues of ABCD
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            k2 = natmsol/2                 !Chains EFGH atomic index
            if(mod(k1,1000)==0)write(*,*)'checking C12 for atom ',k1
            do ic2=nb_res/2+1,nb_res       !residues of EFGH
               do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                  k2 = k2 + 1
                  do i=1,ncrossLJ
                     if(atomtype(k1)/=atomtype(k2))then                !For equal atoms C12 is = C12 of the atom
                        if(atomtype_LJ_i_04(i)==atomtype(k1).and.atomtype_LJ_j_04(i)==atomtype(k2))then
                           C6_IJ(k1,k2)  = C6LJ(i)
                           C12_IJ(k1,k2) = C12LJ(i)
                           nC12rep(k1,k2)=1
                           goto 10
                        elseif(atomtype_LJ_j_04(i)==atomtype(k1).and.atomtype_LJ_i_04(i)==atomtype(k2))then
                           C6_IJ(k1,k2)  = C6LJ(i)
                           C12_IJ(k1,k2) = C12LJ(i)
                           nC12rep(k1,k2)=1
                           goto 10
                        endif
                     endif   
                  end do 
                  10 continue
               end do
            end do   
         end do
      end do 
      write(*,*)
!Version_22
!Replace parameters in the respective arrays  
      write(*,*)'Start checking C6 and C12 iw protein-water [nonbond_params] ...'
      nC12rep_W = 0                        !keeps track of replacements
      k1 = 0                               !Chains ABCD and EFGH atomic index
      do ic1=1,nb_res                      !residues of ABCD and EFGH
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            if(mod(k1,1000)==0)write(*,*)'checking C12 for atom ',k1
            do i=1,ncrossLJ                           
               if(atomtype_LJ_i_04(i)==atomtype(k1).and.atomtype_LJ_j_04(i)=='OW')then
                  C6_IW(k1)  = C6LJ(i)
                  C12_IW(k1) = C12LJ(i)
                  nC12rep_W(k1)=1
                  goto 15  
               elseif(atomtype_LJ_j_04(i)==atomtype(k1).and.atomtype_LJ_i_04(i)=='OW')then
                  C6_IW(k1)  = C6LJ(i)
                  C12_IW(k1) = C12LJ(i)
                  nC12rep_W(k1)=1
                  goto 15     
               endif   
            end do 
            15 continue   
         end do
      end do  
!End Version_22           
   else
      write(*,*)'Geometric rule will be applied to C6 and C12 LJ parm'
      write(*,*)
   endif   
   write(*,*)
   write(*,*)'Ended checking C12 parameters'
   write(*,*)

!end C12 replacements



!check missing replacements
   if(ffname=='GROMOS   ')then 
      write(*,*)'Check for protein-protein vdW missing parameters'
      write(*,*)
      k1 = 0                               !Chains ABCD atomic index
      do ic1=1,nb_res/2                    !residues of ABCD
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            k2 = natmsol/2                 !Chains EFGH atomic index                  
            do ic2=nb_res/2+1,nb_res       !residues of EFGH
               do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                  k2 = k2 + 1
                  if((atomtype(k1)/=atomtype(k2)).and.(nC12rep(k1,k2)==0))then
                     write(*,*)atomtype(k1),atomtype(k2),' - C12(i,j) not found' 
                     write(*,*)'Check for identation problems in the input file'
                     write(*,*)'Increase the nb of atom pairs to read from ffnonbonded.itp'
                     write(*,*)
                     stop
                  endif
               end do
            end do   
         end do
      end do 
   endif
   write(*,*)
   
!Version_22
!check missing replacements
   if(ffname=='GROMOS   ')then 
      write(*,*)'Check for protein-water vdW missing parameters'
      write(*,*)
      k1 = 0                               !Chains ABCD and EFGH atomic index
      do ic1=1,nb_res                      !residues of ABCD and EFGH
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            if(nC12rep_W(k1)==0)then
               write(*,*)atomtype(k1),'-OW',' - C12(i,w) not found' 
               write(*,*)'Check for identation problems in the input file'
               write(*,*)'Increase the nb of atom pairs to read from ffnonbonded.itp'
               write(*,*)
               stop
            endif   
         end do
      end do 
   endif
!End Version_22



!Warn about zero C6 and or C12 parameters
   write(*,*)'Check for zero parameters'
   write(*,*)
   do j = 1,nmolsol*natmsol
      if((C6kjmol(j)==0.or.C12kjmol(j)==0).and.(atomtype(j).ne.'H'))then
         write(*,'(a36,1x,i7,1x,a4,1x,a4,1x,a23)')'WARNING!!! - SOLUTE ATOM =/ Hydrogen',j,'type',atomtype(j),'- ZERO LJ C6 and/or C12'
      endif   
   end do 
   write(*,*)


endif    !End LJ parameters

500  format(3x,i6,2x,a3,2x,i4,4x,a4,1x,a4,1x,a4)

return 

END SUBROUTINE top_gro

SUBROUTINE top_gro_sol(ngroitp,ffname,nLJitp,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,niontypes,niontype,nb_res,nb_res_ch,nchains,segname,resindex,resname,&
                   chrg,atomname,atomtype,nsoltypes,Zattype,mattype,ZMBIO,Ziontype,chain_label_res,Natresid,name_res,numb_res,ksLJ,nattp,ncrossLJ,&
                   C6kjmol,C12kjmol,C6_IJ,C12_IJ,C6_IW,C12_IW)

! System topology for GROMACS
    integer,intent(in)  :: ngroitp
    character(len=9),intent(in)  :: ffname                           ! ffname: GROMOS, AMBER, CHARMM, AMOEBA, OPLS, BLYP etc
    integer,intent(in)  :: nLJitp
    integer,intent(in)  :: natms
    integer,intent(in)  :: nmolsol
    integer,intent(in)  :: natmsol
    integer,intent(in)  :: nmolwat,nwatsites,nions
    integer,intent(in)                                :: nb_res
    integer,intent(in)                                :: nchains     ! Number of protein chains (e.g., HgbS = 8)
    integer,intent(in)                                :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in)  :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)      :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                                :: niontypes   ! Number of ions distinct species (Na+, Cl- etc)
    integer,dimension(niontypes),intent(in)           :: niontype    ! Number of ions of type (Na+, Cl- etc)
    character(len=2),dimension(niontypes),intent(in)  :: Ziontype    ! Atomic symbol of ion types (Na+, Cl- etc)
    character(len=1),dimension(nb_res),intent(in)     :: chain_label_res  ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nchains),intent(in)             :: nb_res_ch        ! number of residues of each chain
    integer,intent(in)                                :: ksLJ             !0: do not read LJ1 2-6 parameters from file 1: read LJ 12-6 parm
    integer,intent(in)                                :: nattp            !number of distinct atomtypes in the force field file (ffnonbonded.itp - Lennard-Jones parm)
    integer,intent(in)                                :: ncrossLJ         !number of LJ126 crossed IJ terms - for reading non-geometric C12
                                                      
    character(len=3),dimension(natms),intent(out)     :: segname
    integer,dimension(natms),intent(out)              :: resindex 
    character(len=4),dimension(natms),intent(out)     :: resname
    real(kind=4),dimension(natms),intent(out)         :: chrg
    character(len=4),dimension(natms),intent(out)     :: atomname
    character(len=4),dimension(natms),intent(out)     :: atomtype 
    real(kind=8),dimension(:),allocatable,intent(out) :: ZMBIO
    integer,dimension(:),allocatable,intent(out)      :: Natresid                           ! Number of atoms of each residue
    character(len=4),dimension(:),allocatable,intent(out)         :: name_res
    integer,dimension(:),allocatable,intent(out)      :: numb_res                           ! Number of the residue: from 1 to the total number of residues of each chain
    real(kind=8),dimension(:),allocatable,intent(out) :: C6kjmol,C12kjmol                   ! Solute Lennard-Jones parameters
    real(kind=8),dimension(:,:),allocatable,intent(out)           :: C6_IJ,C12_IJ           ! Solute Lennard-Jones C6 and C12 crossed parameters
    real(kind=8),dimension(:),allocatable,intent(out)             :: C6_IW,C12_IW           ! Solute-water (oxygen) Lennard-Jones C6 and C12 crossed parameters 
    
    
!Local
    character(len=4),dimension(natms)       :: atomtype_gro
    character(len=4),dimension(natms)       :: atomname_gro
    character(len=4),dimension(natms)       :: resname_gro
    
    character(len=1)                        :: res_label
    character(len=1)                        :: label_res
    character(len=1),dimension(nsoltypes)   :: Zattype_1
    character(len=1),dimension(natms)       :: atomname_1
    
!v29    character(len=7)                        :: resid_label
    character(len=2)                        :: resid_label
!v29    character(len=7)                        :: label_resid
    character(len=2)                        :: label_resid
    
    integer :: i,j,k,iw,kw,res_i,iat,nit,is
    integer :: cgnr                                                   !Charge Group NumbeR
    
    character(len=5),dimension(natmsol)     :: atomtype_LJ            ! input atomtype [read]
    character(len=4),dimension(natmsol)     :: atomtype_LJ_04         !atom type reduced to length 4 to compare with atomtype [length 4]
    real(kind=4)                            :: sm,sc
    character(len=1)                        :: sp
    real(kind=8),dimension(:),allocatable   :: C6LJ,C12LJ             ! Solute Lennard-Jones parameters
    integer                                 :: nat_LJ_01
    
    integer                                 :: k1,k2,ic1,ic2,iats1,iats2
    
    character(len=5),dimension(ncrossLJ)     :: atomtype_LJ_i            ! input atomtype [read]
    character(len=5),dimension(ncrossLJ)     :: atomtype_LJ_j            ! input atomtype [read]
    character(len=4),dimension(ncrossLJ)     :: atomtype_LJ_i_04         !atom type reduced to length 4 to compare with atomtype [length 4]
    character(len=4),dimension(ncrossLJ)     :: atomtype_LJ_j_04         !atom type reduced to length 4 to compare with atomtype [length 4]
    integer                                  :: nf
    integer,dimension(:,:),allocatable       :: nC12rep
    integer,dimension(:),allocatable         :: nC12rep_W    
    real(kind=8)                             :: C6WOkjmol,C12WOkjmol
    

allocate(ZMBIO(natmsol)) 
allocate(Natresid(nb_res))
allocate(name_res(nb_res))
allocate(numb_res(nb_res))
allocate(C6LJ(natmsol),C12LJ(natmsol))
allocate(C6kjmol(natmsol),C12kjmol(natmsol))
allocate(C6_IJ(natmsol,natmsol),C12_IJ(natmsol,natmsol))
allocate(C6_IW(natmsol),C12_IW(natmsol))
allocate(nC12rep(natmsol,natmsol))
allocate(nC12rep_W(natmsol))

ZMBIO = 0.0d0
Natresid = 0
numb_res = 0
   C6LJ = 0.0d0
   C12LJ = 0.0d0
   C6kjmol = 0.0d0
   C12kjmol = 0.0d0
nat_LJ_01 = 0 

res_label = ';'
res_i = 0
iat = 0

!V29 resid_label='residue'
resid_label='nr'
!v29 label_resid='frkafka'
label_resid='fk'

do i=1,nsoltypes 
   Zattype_1(i) = trim(adjustl(Zattype(i)))
   write(*,*)Zattype_1(i)
end do

write(*,*)'Reading topology...'
    
!read GROMACS top header - the number of lines changes with the system - use the label "residue" to find the begining of topology     
!do j = 1,26           
!do j = 1,25
!   read(ngroitp,'()')
!end do

do while(label_resid/=resid_label)
   read(ngroitp,'(a,3x,a2)')label_res,label_resid
!check   write(*,*)label_res,label_resid
end do 
!V29 backspace(ngroitp)                 !Rewind a single line

!read GROMACS topology file
!V29 do j = 1,nmolsol*natmsol + nb_res
do j = 1,nmolsol*natmsol
!v29    read(ngroitp,'(a)')label_res
!v29    write(*,*)label_res
!V29   if(label_res==res_label)then
!V29     res_i = res_i + 1
!NG**      write(*,*)'Reading residue',res_i
!V29   else
      iat = iat + 1
!v29      backspace(ngroitp)                 !Rewind a single line
!Version 10 contact map
!V29      Natresid(res_i) = Natresid(res_i) + 1 
!End version 10
!      read(ngroitp,'(1x,i5,7x,a4,4x,i5,2x,a4,2x,a4)')i,atomtype_gro(iat),resindex(iat),resname_gro(iat),atomname_gro(iat)
!v29      read(ngroitp,'(1x,i5,7x,a4,4x,i5,2x,a4,2x,a4,3x,i4,4x,f7.3)')is,atomtype_gro(iat),resindex(iat),resname_gro(iat),atomname_gro(iat),&
!v29      cgnr,chrg(iat)
          read(ngroitp,'(1x,i5,3x,a2,3x,i3,3x,a3,2x,a4,2x,i3,4x,f9.6)')is,atomtype_gro(iat),resindex(iat),resname_gro(iat),atomname_gro(iat),&
      cgnr,chrg(iat)
      write(atomtype(iat),'(a)') trim(adjustl(atomtype_gro(iat)))
      write(resname(iat) ,'(a)') trim(adjustl(resname_gro(iat)))
      write(atomname(iat),'(a)') trim(adjustl(atomname_gro(iat)))  
      segname(iat) = 'MOL'
      atomname_1(iat) = trim(adjustl(atomname(iat)))                                                 !re-define atomic symbol by a single character 
!NG      write(*,500)iat,segname(iat),resindex(iat),resname(iat),atomname(iat),atomtype(iat)
!      write(*,*)iat,atomtype(iat),chrg(iat)
      do i=1,nsoltypes                                              !Assign atomic masses
         if(atomname_1(iat)==Zattype_1(i)) ZMBIO(iat)=mattype(i)
      end do
!v29      name_res(res_i) = resname(iat)
!      write(*,*) iat,atomname(iat),ZMBIO(iat)                   0k
!V29   endif
      write(*,*)j,is,atomtype_gro(iat),resindex(iat),resname_gro(iat),atomname_gro(iat),&
      cgnr,chrg(iat)
end do

if(res_i/=nb_res)then
   write(*,*)'Error - number of residues in gro.itp different from input'
   write(*,*)res_i,' /= ',nb_res
   stop
endif
   
!Version 10
!
!v29   k = 0
!v29   do i=1,nchains
!v29      do j=1,nb_res_ch(i)
!v29         k=k+1
!v29         numb_res(k) = j
!v29      end do
!v29   end do

!check do j=1,nb_res                                                                       !Loop over total number of residues
!check    write(*,*)j,Natresid(j),name_res(j),numb_res(j),chain_label_res(j)
!check end do
!End version 10
  
do j = 1,nmolsol*natmsol
   if(ZMBIO(j)==0.0)then
      write(*,*)'Error!!! atomic mass undefined'
      write(*,*)'index = ',j,' atom name = ',atomname_1(j)
   endif
end do   

write(*,*)
write(*,*)'Ended reading topology'

!Water topology
iw = nb_res + 1
kw = nmolsol*natmsol + 1 
do j = 1,nmolwat
   segname(kw)  = 'WAT' 
   resindex(kw) =  iw 
   if(nwatsites==3)resname(kw)  = 'SPCE'
   if(nwatsites==4)resname(kw)  = 'TIP4'
   atomname(kw) = 'OH2' 
   atomtype(kw) = 'OW'
!NG**   write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
   kw = kw + 1
   segname(kw)  = 'WAT' 
   resindex(kw) =  iw 
   if(nwatsites==3)resname(kw)  = 'SPCE'
   if(nwatsites==4)resname(kw)  = 'TIP4'  
   atomname(kw) = 'H1' 
   atomtype(kw) = 'HW'
!NG**   write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
   kw = kw + 1
   segname(kw)  = 'WAT' 
   resindex(kw) =  iw 
   if(nwatsites==3)resname(kw)  = 'SPCE'
   if(nwatsites==4)resname(kw)  = 'TIP4'   
   atomname(kw) = 'H2' 
   atomtype(kw) = 'HW'
!NG**   write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
   kw = kw + 1
   if(nwatsites==4)then
      segname(kw)  = 'WAT' 
      resindex(kw) =  iw 
!      if(nwatsites==3)resname(kw)  = 'SPCE'
!      if(nwatsites==4)resname(kw)  = 'TIP4'
      resname(kw)  = 'TIP4'
      atomname(kw) = 'M' 
      atomtype(kw) = 'MW'
!NG**      write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
      kw = kw + 1
   endif 
end do


if(nions > 0)then
   iw = nb_res + 2
   do i = 1,niontypes
      nit = niontype(i) 
      do j = 1,nit
         segname(kw)  = 'ION' 
         resindex(kw) =  iw
         resname(kw)  = 'ION'
         write(atomname(kw),*)Ziontype(i) 
         write(atomtype(kw),*)Ziontype(i)
   !NG**      write(*,500)kw,segname(kw),resindex(kw),resname(kw),atomname(kw),atomtype(kw)
         kw = kw + 1
      end do
      iw = iw + 1
   end do   
endif 

write(*,*)
write(*,*)'Total number of atoms = ',natms
write(*,*)



!Read Lennard-Jones 12-6 parameters - to be used in prot_cont_map calculation

! PROBLEMS: GROMOS54A7 does not use geometric combining rules for all LJ C12 coefficients - dx.doi.org/10.1002/jcc.20090 [J Comput Chem 25: 1656–1676, 2004]
! SOLUTION: First apply geometric rule to all C6 and C12 parameters and store them in a matrix. Second read [ nonbond_params ] from ffnonbonded.itp
!           and replace the C12 parameters which are different from those obtained through the geometric rule.

if(ksLJ==1)then
!   read(nLJitp,*)nattp            [input parameter]
   read(nLJitp,*)
   read(nLJitp,*)
   write(*,*)'Reading Lennard-Jones parameters...'
   write(*,*)'Number of Lennard-Jones atom types = ',nattp
   do i=1,nattp
      read(nLJitp,'(a5,3x,i2,2(4x,f7.3)5x,a1,2(2x,f12.10))')atomtype_LJ(i),is,sm,sc,sp,C6LJ(i),C12LJ(i)
      atomtype_LJ_04(i) = trim(adjustl(atomtype_LJ(i)))
!check      write(*,*)atomtype_LJ(i),atomtype_LJ_04(i),is,sm,sc,sp,C6LJ(i),C12LJ(i)
   end do
   do j = 1,nmolsol*natmsol
!    do j = 1,1                         !testing - 1 atom
      do i=1,nattp 
!check         write(*,*)atomtype(j),atomtype_LJ_04(i)
         if(atomtype(j)==atomtype_LJ_04(i))then
            nat_LJ_01 = nat_LJ_01 + 1
            C6kjmol(j)= C6LJ(i)
            C12kjmol(j)= C12LJ(i)
!check            write(*,*)'Found LJ atom type...'
            goto 5
         endif
      end do
      5 continue
   end do
   
   write(*,*)'Found Lennard-Jones parameters for ',nat_LJ_01,' atoms out of ',nmolsol*natmsol
   write(*,*)
   if(nat_LJ_01.ne.nmolsol*natmsol)then
      write(*,*)'Some solute Lennard-Jones parameters were not found'
      stop
   endif   
    
write(*,*)

!Version22
!SPC/E model - oxygen
    C6WOkjmol  = 0.0026173456d0 
    C12WOkjmol = 2.634129d-6 
!End Version 22

! Calculate interaction IJ Lennard-Jones parameters from geometric rules for C6 and C12 

   k1 = 0                               !Chains ABCD atomic index
   do ic1=1,nb_res/2                    !residues of ABCD
      do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
         k1 = k1 + 1
         k2 = natmsol/2                 !Chains EFGH atomic index                  
         do ic2=nb_res/2+1,nb_res       !residues of EFGH
            do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
               k2 = k2 + 1
               C6_IJ(k1,k2) = dsqrt(C6kjmol(k1)*C6kjmol(k2))                  !kJ/mol nm**6
               C12_IJ(k1,k2) = dsqrt(C12kjmol(k1)*C12kjmol(k2))                  !kJ/mol nm**12
             end do
         end do   
      end do
   end do
   
!Version_22
!Protein-water LJ parameters from geometric rules for C6 and C12
   k1 = 0                               !Chains ABCD and EFGH atomic index
   do ic1=1,nb_res                      !residues of ABCD and EFGH
      do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
         k1 = k1 + 1
         C6_IW(k1)  =  dsqrt(C6kjmol(k1)*C6WOkjmol)                    !kJ/mol nm**6
         C12_IW(k1) =  dsqrt(C12kjmol(k1)*C12WOkjmol)                  !kJ/mol nm**12
      end do
   end do   
!End Version_22
   
! Replace C6 and C12 IJ Lennard-Jones parameters for some atom pairs according to [ nonbond_params ] from ffnonbonded.itp 
! For C6 this is not necessary as they always obey the geometric rule - values are changed nevertheless to use the exact same values in ffnonbonded.itp
   if(ffname=='GROMOS   ')then 
      write(*,'(a71)')'Force Field is GROMOS - read modified cross C12(i,j) parms from input'
      write(*,*)
      read(nLJitp,*)
      read(nLJitp,*)
      read(nLJitp,*)
      write(*,*)'Reading Lennard-Jones C6 and C12 ij [nonbond_params] ...'
      write(*,*)'Number of Lennard-Jones input nonbond_params = ',ncrossLJ
      write(*,*)
      do i=1,ncrossLJ
         read(nLJitp,'(a5,3x,a5,3x,i4,4x,e12.6,4x,e12.6)')atomtype_LJ_i(i),atomtype_LJ_j(i),nf,C6LJ(i),C12LJ(i)
         atomtype_LJ_i_04(i) = trim(adjustl(atomtype_LJ_i(i)))
         atomtype_LJ_j_04(i) = trim(adjustl(atomtype_LJ_j(i)))
!check      write(*,'(i6,1x,a4,1x,a4,3x,e12.6,3x,e12.6)')i,atomtype_LJ_i_04(i),atomtype_LJ_j_04(i),C6LJ(i),C12LJ(i)
      end do
!Replace parameters in the respective arrays  
      write(*,*)'Start checking C6 and C12 ij [nonbond_params] ...'
      nC12rep = 0
      k1 = 0                               !Chains ABCD atomic index
      do ic1=1,nb_res/2                    !residues of ABCD
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            k2 = natmsol/2                 !Chains EFGH atomic index
            if(mod(k1,1000)==0)write(*,*)'checking C12 for atom ',k1
            do ic2=nb_res/2+1,nb_res       !residues of EFGH
               do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                  k2 = k2 + 1
                  do i=1,ncrossLJ
                     if(atomtype(k1)/=atomtype(k2))then                !For equal atoms C12 is = C12 of the atom
                        if(atomtype_LJ_i_04(i)==atomtype(k1).and.atomtype_LJ_j_04(i)==atomtype(k2))then
                           C6_IJ(k1,k2)  = C6LJ(i)
                           C12_IJ(k1,k2) = C12LJ(i)
                           nC12rep(k1,k2)=1
                           goto 10
                        elseif(atomtype_LJ_j_04(i)==atomtype(k1).and.atomtype_LJ_i_04(i)==atomtype(k2))then
                           C6_IJ(k1,k2)  = C6LJ(i)
                           C12_IJ(k1,k2) = C12LJ(i)
                           nC12rep(k1,k2)=1
                           goto 10
                        endif
                     endif   
                  end do 
                  10 continue
               end do
            end do   
         end do
      end do 
      write(*,*)
!Version_22
!Replace parameters in the respective arrays  
      write(*,*)'Start checking C6 and C12 iw protein-water [nonbond_params] ...'
      nC12rep_W = 0                        !keeps track of replacements
      k1 = 0                               !Chains ABCD and EFGH atomic index
      do ic1=1,nb_res                      !residues of ABCD and EFGH
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            if(mod(k1,1000)==0)write(*,*)'checking C12 for atom ',k1
            do i=1,ncrossLJ                           
               if(atomtype_LJ_i_04(i)==atomtype(k1).and.atomtype_LJ_j_04(i)=='OW')then
                  C6_IW(k1)  = C6LJ(i)
                  C12_IW(k1) = C12LJ(i)
                  nC12rep_W(k1)=1
                  goto 15  
               elseif(atomtype_LJ_j_04(i)==atomtype(k1).and.atomtype_LJ_i_04(i)=='OW')then
                  C6_IW(k1)  = C6LJ(i)
                  C12_IW(k1) = C12LJ(i)
                  nC12rep_W(k1)=1
                  goto 15     
               endif   
            end do 
            15 continue   
         end do
      end do  
!End Version_22           
   else
      write(*,*)'Geometric rule will be applied to C6 and C12 LJ parm'
      write(*,*)
   endif   
   write(*,*)
   write(*,*)'Ended checking C12 parameters'
   write(*,*)

!end C12 replacements



!check missing replacements
   if(ffname=='GROMOS   ')then 
      write(*,*)'Check for protein-protein vdW missing parameters'
      write(*,*)
      k1 = 0                               !Chains ABCD atomic index
      do ic1=1,nb_res/2                    !residues of ABCD
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            k2 = natmsol/2                 !Chains EFGH atomic index                  
            do ic2=nb_res/2+1,nb_res       !residues of EFGH
               do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                  k2 = k2 + 1
                  if((atomtype(k1)/=atomtype(k2)).and.(nC12rep(k1,k2)==0))then
                     write(*,*)atomtype(k1),atomtype(k2),' - C12(i,j) not found' 
                     write(*,*)'Check for identation problems in the input file'
                     write(*,*)'Increase the nb of atom pairs to read from ffnonbonded.itp'
                     write(*,*)
                     stop
                  endif
               end do
            end do   
         end do
      end do 
   endif
   write(*,*)
   
!Version_22
!check missing replacements
   if(ffname=='GROMOS   ')then 
      write(*,*)'Check for protein-water vdW missing parameters'
      write(*,*)
      k1 = 0                               !Chains ABCD and EFGH atomic index
      do ic1=1,nb_res                      !residues of ABCD and EFGH
         do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
            k1 = k1 + 1
            if(nC12rep_W(k1)==0)then
               write(*,*)atomtype(k1),'-OW',' - C12(i,w) not found' 
               write(*,*)'Check for identation problems in the input file'
               write(*,*)'Increase the nb of atom pairs to read from ffnonbonded.itp'
               write(*,*)
               stop
            endif   
         end do
      end do 
   endif
!End Version_22



!Warn about zero C6 and or C12 parameters
   write(*,*)'Check for zero parameters'
   write(*,*)
   do j = 1,nmolsol*natmsol
      if((C6kjmol(j)==0.or.C12kjmol(j)==0).and.(atomtype(j).ne.'H'))then
         write(*,'(a36,1x,i7,1x,a4,1x,a4,1x,a23)')'WARNING!!! - SOLUTE ATOM =/ Hydrogen',j,'type',atomtype(j),'- ZERO LJ C6 and/or C12'
      endif   
   end do 
   write(*,*)


endif    !End LJ parameters

500  format(3x,i6,2x,a3,2x,i4,4x,a4,1x,a4,1x,a4)

return 

END SUBROUTINE top_gro_sol





!LAB 2023 BIO WATER

SUBROUTINE wat_bio_histo(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                     biowsamp,biowRadius,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO,WATMODEL)               
!Calculate the distance of each water molecule to the heavy atoms of the protein;Water molecules at a distance below biowRadius are considered biological water
!Print the time average and stdev for biological and bulk water
    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
    character(len=7),intent(in)                      :: WATMODEL
    character(len=4),dimension(natms)                :: atomname    
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)  
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment   
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze    
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),intent(in)                          :: biowRadius  ! Distance relative to solute heavy atoms for onset of bulk water
    integer,intent(in)                               :: biowsamp
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
   
! Local variables    
    character(len=1),dimension(natms)                 :: atomname_1
    integer                                           :: nw_bulk                    ! Bulk water counter
    integer                                           :: NBIO_wat                   ! Biological Water counter
    integer,dimension(:),allocatable                  :: kwlabel                    ! 
    integer,dimension(:),allocatable                  :: kwbio_time                 ! 
    real(kind=4)                                      :: AV_NBIO_wat,av_nw_bulk
    real(kind=4)                                      :: stdev_bio,stdev_bio_2
    real(kind=4)                                      :: stdev_bulk,stdev_bulk_2    
    integer      :: n0, n1
    integer      :: n21
    integer      :: i, j, k
    integer      :: k_s
    integer      :: io, ilw
    integer,parameter                                 :: INTMAX = 2147483647
    integer                                           :: natcheck
    integer                                           :: at_id                    ! solute atomic index to analyse the respective solvation environment
    integer                                           :: kstps
    logical                                           :: filex
    logical                                           :: lprnt_biow
    integer                                           :: biow_vmd
!End Local variables    
    
! 
!NG real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!   real(kind=4)                                     :: x(natms),y(natms),z(natms) 
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8)                                     :: dx,dy,dz,dr,dr2


    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
allocate(kwlabel(nmolwat))          
allocate(kwbio_time(nstep))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj

inquire(file='Bio_Water_out/log_Bio_Wat.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Bio_Water_out')

n0 = 100
n1 = 110

n21 = 310

open(n0,file='Bio_Water_out/log_Bio_Wat.dat',status='unknown',action='write')
open(n1,file='Bio_Water_out/Bio_Wat_hist.dat',status='unknown',action='write')
open(n21,file='Bio_Water_out/Bio_Wat_vmd.gro',status='unknown',action='write')

lprnt_biow = .true.
biow_vmd = biowsamp*50000

!VMD configuration .gro with solute and biological water
if (lprnt_biow) then
   write(*,*)
   write(*,*)'WARNING!!! VISUAL'
   write(*,*)'PRINT WATER MOLECULES FOR A SINGLE CONFIGURATION - VISUAL'
   write(*,*)'VMD gro cfg will be printed at time-step ',biow_vmd
   write(*,*)'SOLUTE SHOULD BE CENTERED IN THE BOX'
   write(*,*)'gmx trjconv -f npt.xtc -s npt.tpr -center -pbc mol -o npt_pbc_mol.xtc'
   write(*,*)
   write(n0,*)
   write(n0,*)'WARNING!!! VISUAL'
   write(n0,*)'PRINT WATER MOLECULES FOR A SINGLE CONFIGURATION - VISUAL'
   write(n0,*)'VMD gro cfg will be printed at time-step ',biow_vmd
   write(n0,*)'SOLUTE SHOULD BE CENTERED IN THE BOX'
   write(n0,*)'gmx trjconv -f npt.xtc -s npt.tpr -center -pbc mol -o npt_pbc_mol.xtc'
   write(n0,*)
!   
   write(n21,'(a52)')'Protein and Bio-Water - add nb atoms of water below'                               !gro file comment line
   write(n21,'(i4)')nmolsol*natmsol
   ilw = nmolsol*natmsol
endif

kstps = 0
NBIO_wat = 0
nw_bulk = 0 
kwlabel = 0                                                               !array with waters ID for bulk water counting
kwbio_time = 0 
                                                           !array with number of bio water found at each time-step for stdev calculation

AV_NBIO_wat = 0.0 
av_nw_bulk = 0.0
stdev_bio= 0.0 
stdev_bio_2 = 0.0      
stdev_bulk = 0.0
stdev_bulk_2 = 0.0

! Start biological water calculation
write(*,*)
write(*,*)'Starting biological water calculation...'
write(*,*)'Routine wat_bio_histo'
write(*,*)'Results are printed to Bio_Water_out'
write(*,*)'water sampling frequency (steps) = ',biowsamp
write(*,*)'Biological Water cut-off (Ang) =',biowRadius 
write(*,*)'total number of time-steps =',nstep
write(*,*)

write(n0,*)
write(n0,*)'Starting biological water calculation...'
write(n0,*)'Routine wat_bio_histo'
write(n0,*)'Results are printed to Bio_Water_out'
write(n0,*)'water sampling frequency (steps) = ',biowsamp
write(n0,*)'Biological Water cut-off (Ang) =',biowRadius 
write(n0,*)'total number of time-steps =',nstep
write(n0,*)

do j = 1,nstep
   kwlabel = 0
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!Print protein atoms for visual gro
!WARNING!!! write protein only for time-step biow_vmd - NOT EVERY biowsamp 
!NG   if (j==biowsamp .and. lprnt_biow)then 
   if (j==biow_vmd .and. lprnt_biow)then  
      do i=1,nmolsol*natmsol       
         write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'1MOL',trim(adjustl(atomname(i))),i,0.1*x(i),0.1*y(i),0.1*z(i)
      end do
   endif   
!   if (j==biow_vmd .and. lprnt_biow)write(n21,'(3(F9.4))')0.1*cell

! sample frequency - assess biological water every biowsamp time-steps  

   if(j==1.or.mod(j,biowsamp)==0)then
      write(*,*)'time-step',j,' sampling water molecules'
      kwlabel = 0                                                                  !array with waters ID for bulk water counting
      if (j==biow_vmd .and. lprnt_biow)write(*,*)j,' visual bio water'
      kstps = kstps + 1                                                            !count number of origins sampled
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
        
!Loop over protein heavy atoms 
         do at_id=1,nmolsol*natmsol                           !Loop over protein atomic species 
            atomname_1(at_id) = trim(adjustl(atomname(at_id)))  
!check      if (atomname_1(at_id) /=  'H')write(*,*)atomname_1(at_id)
            if (atomname_1(at_id) /=  'H')then                !use protein heavy atoms only
               k= at_id               
               dx = x(i) - x(k)
               dy = y(i) - y(k)
               dz = z(i) - z(k)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               
               IF (dr <= biowRadius)THEN
                  NBIO_wat = NBIO_wat + 1                                                  !accumulate number of water found for average
                  kwlabel(io) = 1
                  kwbio_time(kstps)= kwbio_time(kstps) + 1                                   !store number of bio waters found for each cfg - stdev                  
                  if(j==biow_vmd .and. lprnt_biow)then
!                     write(*,*)ilw
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i),ilw,0.1*x(i),0.1*y(i),0.1*z(i)
                     ilw = ilw + 1                                                                                             
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i+1),ilw,0.1*x(i+1),0.1*y(i+1),0.1*z(i+1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i+2),ilw,0.1*x(i+2),0.1*y(i+2),0.1*z(i+2)                          
                  endif  
                  goto 10                                                                  !jump to next water molecule
               ENDIF 
            endif

         end do                                                                     !end loop over atomic species
         if(kwlabel(io) == 0) nw_bulk = nw_bulk + 1                                 !water molecule not labeled bio water 
         10 continue 
      end do                                                                        !end loop over water oxygens
      write(*,*)kstps,kwbio_time(kstps),nmolwat-kwbio_time(kstps)                   !bio water and bulk water found every time-step
      write(*,*)kstps,NBIO_wat,nw_bulk 
   endif                                                                            !end sampling of waters
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj
   if (j==biow_vmd .and. lprnt_biow)write(n21,'(3(F9.4))')0.1*cell
   if (j==biow_vmd .and. lprnt_biow)write(n0,*)
   if (j==biow_vmd .and. lprnt_biow)write(n0,*)kwbio_time(kstps),' bio water molecules'
   if (j==biow_vmd .and. lprnt_biow)write(n0,*)kwbio_time(kstps)*3,' bio water atoms - add to visual nb atoms'
   if (j==biow_vmd .and. lprnt_biow)write(n0,*)
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                                                                                                                                        

!Average number of waters on each environment

write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps

!Print mean values and stdev
if(NBIO_wat > INTMAX)write(*,*)'WARNING!!! Number of accumulated bio waters too large'
if(NBIO_wat > INTMAX)write(n0,*)'WARNING!!! Number of accumulated bio waters too large'
if(nw_bulk > INTMAX)write(*,*)'WARNING!!! Number of accumulated bulk waters too large'
if(nw_bulk > INTMAX)write(n0,*)'WARNING!!! Number of accumulated bulk waters too large'
AV_NBIO_wat = FLOAT(NBIO_wat)/FLOAT(kstps)                                              
av_nw_bulk  = FLOAT(nw_bulk)/FLOAT(kstps)   

do k_s=1,kstps
   stdev_bio_2 = stdev_bio_2 + (FLOAT(kwbio_time(k_s)) - AV_NBIO_wat)**2.0                
   stdev_bulk_2 = stdev_bulk_2 + (FLOAT(nmolwat-kwbio_time(k_s)) - av_nw_bulk)**2.0   
end do
stdev_bio = sqrt(stdev_bio_2/float(kstps-1))
stdev_bulk = sqrt(stdev_bulk_2/float(kstps-1))
WRITE(n1,*)
WRITE(n1,*)'AVERAGE <Nbio_wat> +/- stdev = ',AV_NBIO_wat,stdev_bio
WRITE(n1,*)'AVERAGE <Nbulk_wat> +/- stdev = ',av_nw_bulk,stdev_bulk 
WRITE(n1,*)
WRITE(n1,*)'CONTROL SUM of WATER MOLECULES = ',AV_NBIO_wat+av_nw_bulk,nmolwat

close(n21) 
close(n0) 
close(n1) 


deallocate (kwlabel)
deallocate(kwbio_time)
deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT(27X,F14.4,5X,F17.7)
    return

END SUBROUTINE wat_bio_histo



!END LAB 2023 BIO WATER


SUBROUTINE wat_env_orient_tcfs(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nions,nwatsites,atomname,&
                                mtdelay,m_time,rwtcfbulk,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO,&
                                resindex,resname,atomtype,chrg,nb_res,chain_label_res)                 
! Calculate water's OH reorientational tcfs - 1st, 2nd and 3rd Legendre Polynomials - Solutions
! Every m_time time-steps the code samples waters on different environments 
    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
!    character(len=4),dimension(natms),intent(out)    :: atomname
    character(len=4),dimension(natms)                :: atomname
!Version 9.0    
    integer,dimension(natms),intent(in)              :: resindex 
    character(len=4),dimension(natms),intent(in)     :: resname
    character(len=4),dimension(natms),intent(in)     :: atomtype 
    real(kind=4),dimension(natms),intent(in)         :: chrg
!    
    integer,intent(in)                               :: mtdelay     ! reorientational tcf delay time (time-windows/ps)
    integer,intent(in)                               :: m_time      ! time freq. to sample water OH groups in radial shells for reor. tcf calculation
    real(kind=8)                                     :: rwtcfbulk   ! onset (Ang) for bulk water relative to the solute COM - for bulk reorientational tcf
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
    integer,intent(in)                               :: nb_res      ! Number of residues [GROMCAS]
    character(len=1),dimension(nb_res),intent(in)    :: chain_label_res ! Chain label (e.g., A, B, C etc) internal code use
   
! Local variables
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8),dimension(:,:,:),allocatable        :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable        :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                                     :: TSTEP,PSECS
    real(kind=8)                                     :: dx,dy,dz,dr,dr2
    real(kind=8)                                     :: dx_ref,dy_ref,dz_ref,dr_ref,dr2_ref
    real(kind=8)                                     :: xoh,yoh,zoh,sqoh
    integer                                          :: NDELS,KINTVL,NDELS1
    integer,dimension(:,:,:),allocatable             :: NOHidO,NOHidH
    integer,dimension(:,:),allocatable               :: ntime
    integer,dimension(:,:),allocatable               :: nOHorg
    integer,dimension(:),allocatable                 :: kv
    integer,dimension(:,:),allocatable               :: kwlabel      !index of waters in each environment
    integer      :: n0, n1, n2, n3, n11, n12, n13, n14
    integer      :: i, j, k, L, LR
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer,dimension(:),allocatable           :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable           :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                               :: dr_bulk, dr4nb
    real(kind=8)                               :: RADMAX,RADMIN
    integer                                    :: norgmax
    integer                                    :: NRSHMAX   
    integer,dimension(:),allocatable           :: norg_0HB  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!Orientational tcf
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b         !bulk
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs        !hsh
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: tcf_mean_1   ,tcf_mean_2   ,tcf_mean_3   
    real(kind=8)                                 :: X0TLEG,ARGLEG1,ARGLEG2,ARGLEG3
    integer,dimension(:),allocatable             :: nOHmean
    real(kind=8),dimension(:),allocatable        :: stdevnOH_sum,stdevnOH

    integer                                      :: bulk_flag,kr_start
    integer                                      :: at_id                ! solute atomic index to analyse the respective solvation environment
    integer                                      :: maxwat_solv
    integer                                      :: maxwat_bulk                !maximum number of waters in the bulk
    integer,dimension(:),allocatable             :: koc                  !Nb of origins corrected for origins where no waters were found in a given environment

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj
allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    
!end xtc2fortran_trj
     
inquire(file='Leg_env_out/log_env_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Leg_env_out')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
n11 = 210
n12 = 220
n13 = 230
n14 = 240

open(n0,file='Leg_env_out/log_env_tcf.dat',status='unknown',action='write')
open(n1,file='Leg_env_out/wat_reor_tcf_1.dat',status='unknown',action='write')
open(n2,file='Leg_env_out/wat_reor_tcf_2.dat',status='unknown',action='write')
open(n3,file='Leg_env_out/wat_reor_tcf_3.dat',status='unknown',action='write')
open(n11,file='Leg_env_out/tcf_1_rot_tau.dat',status='unknown',action='write')
open(n12,file='Leg_env_out/tcf_2_rot_tau.dat',status='unknown',action='write')
open(n13,file='Leg_env_out/tcf_3_rot_tau.dat',status='unknown',action='write') 
open(n14,file='Leg_env_out/solute_id_nth.dat',status='unknown',action='write') 

maxwat_solv = 10000                           !for allocation purposes - maximum number of waters in solvation shells (move to input if needed)
if(nmolwat<maxwat_solv)maxwat_solv = nmolwat
maxwat_bulk = 500*2                     !number of OH groups 2*number of molecules
if(nmolwat<maxwat_bulk)maxwat_bulk = nmolwat

!time-window for calculation of each tcf
! IDINT - convert to integer - truncate

NDELS = IDINT(mtdelay*1000d0/dt)                  !sets the delay time (in frames)

dr_bulk = rwtcfbulk                               !set distance from the solute com to assume bulk water
bulk_flag = 1                                     !default - calculate tcf for bulk water
kr_start = 1                                      !environments loop starter - if bulk_flag = 0; kr_start = 2 

if(dr_bulk==0.0)then
   bulk_flag = 0
   write(*,*)
   write(*,*)'Warning!!! No bulk region defined for reorientational tcf calculation'
   write(*,*)'Bulk tcf will not be calculated accordingly'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! No bulk region defined for reorientational tcf calculation'
   write(n0,*)'Bulk tcf will not be calculated accordingly'
   write(n0,*)
else 
   write(*,*)
   write(*,*)'Warning!!! Maximum bulk waters for reorientational tcf  = ',maxwat_bulk/2
   write(*,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! Maximum bulk waters for reorientational tcf = ',maxwat_bulk/2
   write(n0,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(n0,*)   
endif   

NRSHMAX = 4                                      !maximum number of environments
if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/m_time                                  !for memory allocation purposes

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))  
NTDW4O = 0; NONTDW = 0 

!orient. tcf

!allocate(RXt0(norgmax,2*nmolwat,NRSHMAX),RYt0(norgmax,2*nmolwat,NRSHMAX),RZt0(norgmax,2*nmolwat,NRSHMAX))
!allocate(RXt(norgmax,2*nmolwat,NRSHMAX),RYt(norgmax,2*nmolwat,NRSHMAX),RZt(norgmax,2*nmolwat,NRSHMAX))
!allocate(NOHidO(norgmax,2*nmolwat,NRSHMAX),NOHidH(norgmax,2*nmolwat,NRSHMAX))

allocate(RXt0(norgmax,2*maxwat_solv,NRSHMAX),RYt0(norgmax,2*maxwat_solv,NRSHMAX),RZt0(norgmax,2*maxwat_solv,NRSHMAX))
allocate(RXt(norgmax,2*maxwat_solv,NRSHMAX),RYt(norgmax,2*maxwat_solv,NRSHMAX),RZt(norgmax,2*maxwat_solv,NRSHMAX))
allocate(NOHidO(norgmax,2*maxwat_solv,NRSHMAX),NOHidH(norgmax,2*maxwat_solv,NRSHMAX))

allocate(nOHorg(norgmax,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(tcf_leg_1_b(NDELS),tcf_leg_1_hs(NDELS),tcf_leg_1_th(NDELS),tcf_leg_1_nth(NDELS))
allocate(tcf_leg_2_b(NDELS),tcf_leg_2_hs(NDELS),tcf_leg_2_th(NDELS),tcf_leg_2_nth(NDELS))
allocate(tcf_leg_3_b(NDELS),tcf_leg_3_hs(NDELS),tcf_leg_3_th(NDELS),tcf_leg_3_nth(NDELS))
allocate(tcf_mean_1(NDELS,NRSHMAX) ,tcf_mean_2(NDELS,NRSHMAX)  ,tcf_mean_3(NDELS,NRSHMAX))
tcf_leg_1_b = 0.0d0;tcf_leg_1_hs = 0.0d0;tcf_leg_1_th = 0.0d0;tcf_leg_1_nth = 0.0d0
tcf_leg_2_b = 0.0d0;tcf_leg_2_hs = 0.0d0;tcf_leg_2_th = 0.0d0;tcf_leg_2_nth = 0.0d0
tcf_leg_3_b = 0.0d0;tcf_leg_3_hs = 0.0d0;tcf_leg_3_th = 0.0d0;tcf_leg_3_nth = 0.0d0
tcf_mean_1  = 0.0d0;tcf_mean_2   = 0.0d0;tcf_mean_3   = 0.0d0 
ntime = 1; nOHorg = 0

allocate(norg_0HB(NRSHMAX))
norg_0HB = 0
allocate(kv(NRSHMAX))
allocate(kwlabel(NRSHMAX,nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled) 
allocate(nOHmean(NRSHMAX))
allocate(stdevnOH_sum(NRSHMAX),stdevnOH(NRSHMAX))
allocate(koc(NRSHMAX))
nOHmean = 0
stdevnOH_sum = 0
stdevnOH = 0
kwlabel = 0
koc = 0

 cmxsol = 0.0d0; cmysol=0.0d0; cmzsol=0.0d0

write(*,*)'Solute atoms id and HSh onset/cut-off (Ang)' 
write(n0,*)'Solute atoms id and HSh onset/cut-off (Ang)'
do i = 1,natHSh
   write(*,'(i6,1x,f7.2,a1,f7.2)')natsol_id(i),oatsol_hs(i),'/',ratsol_hs(i)
   write(n0,'(i6,1x,f7.2,a1,f7.2)')natsol_id(i),oatsol_hs(i),'/',ratsol_hs(i)
end do
write(*,*)
write(n0,*)

! Start tcfs calculation

write(*,*)'Starting water environment reorientation tcfs calculation...'
write(*,*)'Routine wat_env_tcfs'
write(*,*)'Results are printed to Leg_env_out'
write(*,*)'water tcfs delay time-window (ps) =',mtdelay
write(*,*)'tcfs OH sampling frequency (steps) = ',m_time
write(*,*)'total number of time-steps =',nstep
write(*,*)'tcfs delay-time (steps) =',NDELS
write(*,*)
write(n0,*)'Starting water environment reorientation tcfs calculation...'
write(n0,*)'Routine wat_env_tcfs'
write(n0,*)'Results are printed to Leg_env_out'
write(n0,*)'water tcfs delay time-window (ps) =',mtdelay
write(n0,*)'tcfs OH sampling frequency (steps) = ',m_time
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'tcfs delay-time (steps) =',NDELS
write(n0,*)
write(n14,*)'Solute-ref. atom; ','Solute-vertex-atom; ','Sol-vertex-RESindex; ','Sol-vertex-CHAIN ', 'Sol-vertex-RESname; ','Sol-vertex-ATOMname; ','Sol-vertex-charge; ','dist_to_res/Ang'

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!check      if(j==1)write(*,'(i5,1x,A4,1x,3(F14.5))')natmid(i),atomname(i),x(i),y(i),z(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN  
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
      
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules in specific environments
   if((j==1.or.mod(j,m_time)==0).and.(j.le.nstep-NDELS))then               !time origin: sample OH groups in specific environments
      write(n14,*)
      write(n14,*)'time-step = ',j,' out of ',nstep
      write(n14,*)
      NTDW4O = 0; NONTDW = 0
      write(*,*)'time-step',j,' sampling for OH groups'
      kv = 0                                                               !OH groups in each environment 
      ko = ko + 1                                                          !number of origins = number of blocks used to compute the tcfs mean and stdev
      kwlabel = 0
   
!Find solute centre of mass for bulk water tcfs calculation     

      if(bulk_flag==1)call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                              !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water      
         dx = x(i) - cmxsol
         dy = y(i) - cmysol
         dz = z(i) - cmzsol
!         dx = dx - anint(dx/cell(1))*cell(1)
!         dy = dy - anint(dy/cell(2))*cell(2)
!         dz = dz - anint(dz/cell(3))*cell(3)
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
!check         write(*,*)'distance to the centre of mass = ',dr
         if(bulk_flag==0)dr_bulk = 100.d0*cell(1)                            !define a distance larger than the box "radius"
!Version 13         if(dr>=dr_bulk)then
         if(dr>=dr_bulk.and.kv(1)<maxwat_bulk)then
            if(bulk_flag==0)write(*,*)'Warning!!! Found a water molecule in the bulk'
!bulk tcf(t0)             
            do iH = 1,2
               nradsh = 1                                     !bulk environment
               kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
               kvi = kv(nradsh)                               !OH number
!check               write(*,*)'bulk',j,i,iH,kvi
               nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
               NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
               NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
               dx = x(i+iH) - x(i)
               dy = y(i+iH) - y(i)
               dz = z(i+iH) - z(i)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               xoh = dx
               yoh = dy
               zoh = dz
               sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
               RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = yoh/sqoh
               RZt0(ko,kvi,nradsh) = zoh/sqoh 
               if(kvi >=2 *maxwat_solv)then
                  write(*,*)'error - too many waters in the bulk'
                  write(*,*)'action - increase maxwat_solv and re-compile'
                  write(*,*)'or avoid bulk calculation - make cut-off = 0 '
                  stop
               endif   
            end do         
!bulk tcf(t0) end                     
         endif 
         
!Hydration shell - tetrahedral and non tetrahedral waters in the first HSh         
         do at_id=1,natHSh                                      !Loop over atomic species to analyze the environment         
            k = natsol_id(at_id)
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2.0d0 + dy**2.0d0 + dz**2.0d0
            dr  = dsqrt(dr2)
!Hydration shell tcf(t0)                        
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then                                 !water within hydration shell radius 
!Avoid water molecules repetition
               nradsh = 2                                               !hydration shell environment
               if(kwlabel(nradsh,io)==0)then
                  kwlabel(nradsh,io)=1                                                       !Oxygen id for correction of the nb of waters
               else
                  goto 10                                                                    !sample another water molecule
               endif   
!End avoid water molecules repitition 
               do iH = 1,2
                  kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                  kvi = kv(nradsh)                               !OH number
!check                  write(*,*)'HSh',j,i,iH,kvi
                  nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                  NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  xoh = dx
                  yoh = dy
                  zoh = dz
                  sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                  RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = yoh/sqoh
                  RZt0(ko,kvi,nradsh) = zoh/sqoh
                  if(kvi >= 2*maxwat_solv)then
                     write(*,*)'error - too many waters in the HSh'
                     write(*,*)'action - increase maxwat_solv and re-compile'
                     stop
                  endif   
               end do         
!Hydration shell tcf(t0) end                              

!Find tetrahedral and non-tetrahedral waters in the hydration shell for sub-ensemble tetrahedral and non-tetrahedral tcfs        
!FIND THE FOUR NEAREST O ATOMS J TO EACH OXYGEN I
               RADMAX  = 25.0D0
               RADMIN  =  0.0D0        
               JT      = 1
               DO WHILE(JT.LE.4) 
                  do kw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
                     if(kw.ne.i)then
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOMS i AND OXYGEN ATOMS kw  
                        dx = x(i) - x(kw)
                        dy = y(i) - y(kw)
                        dz = z(i) - z(kw)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        IF((dr>RADMIN).AND.(dr<RADMAX))THEN
                           RADMAX = dr
!NTDW4O STORES THE O ATOMS NEIGHBORS OF i: JT = 1 FIRST NEIGHBOR; JT = 2 SECOND NEIGHBOR; JT = 3 THIRD NEIGHBOR; JT = 4 FOURTH NEIGHBOR
                           NTDW4O(JT) = kw
                        ENDIF
                     endif
                  end do
!check                  write(*,*)j,i,io,JT,NTDW4O(JT)
                  JT = JT + 1
                  RADMIN = RADMAX
                  RADMAX = 25.0D0
               END DO
!CHECK IF SOLUTE ATOMS ARE CLOSER THAN THE 4th NEAREST O ATOM          
               KI=NTDW4O(4)                                       !Oxygen ID of the 4th nearest water neighbor 
               dx = x(i) - x(KI)
               dy = y(i) - y(KI)
               dz = z(i) - z(KI)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               dr4nb = dr                                           !distance to the fourth nearest neighbor
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY SOLUTE ATOM IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = 1                           !run over all solute atoms - heavy and non-heavy atoms
               DO WHILE(L.LE.nmolsol*natmsol)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN             
                     NONTDW(io) = 1                                  !water i is non-tetrahedral
!Version 9 identify the solute atom to which the water molecule is nearby 
!Find distance of the solute atom near the water to the ref. solute atom
                     dx_ref = x(k)-x(L)
                     dy_ref = y(k)-y(L)
                     dz_ref = z(k)-z(L)
                     dx_ref = dx_ref - dnint(dx_ref/cell(1))*cell(1)
                     dy_ref = dy_ref - dnint(dy_ref/cell(2))*cell(2)
                     dz_ref = dz_ref - dnint(dz_ref/cell(3))*cell(3)
                     dr2_ref = dx_ref**2 + dy_ref**2 + dz_ref**2
                     dr_ref  = dsqrt(dr2_ref)
                     LR = resindex(L)
                     write(n14,'(2(1x,I5),1x,i4,1x,A1,1x,a4,1x,a5,1x,4(1x,f7.3))')&
                     natsol_id(at_id),L,resindex(L),chain_label_res(LR),resname(L),atomname(L),chrg(L),dr,dr_ref,ratsol_hs(at_id)                    
                     goto 5
                  ENDIF
                  L = L+1
               END DO     
               5 continue           
               if(NONTDW(io)==1)then                                 !Non-tetrahedral water               
!Hydration shell - non tetrahedral tcf(t0)             
                  do iH = 1,2
                     nradsh = 4                                      !hydration shell (non-tetrahedral) environment
                     kv(nradsh) = kv(nradsh) + 1                     !number of water OH groups (index)
                     kvi = kv(nradsh)                                !OH number
!check                     write(*,*)'non-tetrah',j,i,iH,kvi
                     nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1       !number of OH groups in environment nradsh at origin ko
                     NOHidO(ko,kvi,nradsh) = i                       ! Oxygen id
                     NOHidH(ko,kvi,nradsh) = iH                      ! Hydrogen id 
                     dx = x(i+iH) - x(i)
                     dy = y(i+iH) - y(i)
                     dz = z(i+iH) - z(i)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     xoh = dx
                     yoh = dy
                     zoh = dz
                     sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                     RXt0(ko,kvi,nradsh) = xoh/sqoh                  !origin = ko, OH number = kvi, environment type = nradsh
                     RYt0(ko,kvi,nradsh) = yoh/sqoh
                     RZt0(ko,kvi,nradsh) = zoh/sqoh     
                  end do         
!Hydration shell - non tetrahedral tcf(t0) end                                           
               elseif(NONTDW(io)==0)then                             !Tetrahedral water 
!Hydration shell - tetrahedral tcf(t0)             
                  do iH = 1,2
                     nradsh = 3                                      !hydration shell (tetrahedral) environment
                     kv(nradsh) = kv(nradsh) + 1                     !number of water OH groups (index)
                     kvi = kv(nradsh)                                !OH number
!check                     write(*,*)'tetrah',j,i,iH,kvi
                     nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1       !number of OH groups in environment nradsh at origin ko
                     NOHidO(ko,kvi,nradsh) = i                       !Oxygen id
                     NOHidH(ko,kvi,nradsh) = iH                      !Hydrogen id 
                     dx = x(i+iH) - x(i)
                     dy = y(i+iH) - y(i)
                     dz = z(i+iH) - z(i)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     xoh = dx
                     yoh = dy
                     zoh = dz
                     sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                     RXt0(ko,kvi,nradsh) = xoh/sqoh                  !origin = ko, OH number = kvi, environment type = nradsh
                     RYt0(ko,kvi,nradsh) = yoh/sqoh
                     RZt0(ko,kvi,nradsh) = zoh/sqoh     
                  end do         
!Hydration shell - tetrahedral tcf(t0) end  
               endif
!            
            endif                                                    !end water is in the coordination sphere      
         end do                                                      !end loop over atomic species
         10 continue
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nOHorg(ko,1),nOHorg(ko,2),nOHorg(ko,3),nOHorg(ko,4)
   endif                                                             !end sampling new time-origin      
   
!calculate tcfs - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nOHorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then                !ntime counts the delay times for which the tcf has already been calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nOHorg(kt,kr)                                       !loop over OH groups in environment kr: nOHorg can be zero for a given time-origin            
!Check if water molecule is still in the specific environment and update number of OH groups for tcfs for waters at every time in a given environment             

!End check environment           
               i=NOHidO(kt,iw,kr)
               ihat = NOHidH(kt,iw,kr)                                    !Hydrogen atom (ihat = 1 or 2)        
               dx = x(i+ihat) - x(i)
               dy = y(i+ihat) - y(i)
               dz = z(i+ihat) - z(i)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               xoh = dx
               yoh = dy
               zoh = dz
               sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
               RXt(kt,iw,kr) = xoh/sqoh
               RYt(kt,iw,kr) = yoh/sqoh
               RZt(kt,iw,kr) = zoh/sqoh
               X0TLEG = RXt0(kt,iw,kr)*RXt(kt,iw,kr)+RYt0(kt,iw,kr)*RYt(kt,iw,kr)+RZt0(kt,iw,kr)*RZt(kt,iw,kr)
               ARGLEG1= X0TLEG
               ARGLEG2= 0.5d0*(3.0d0*X0TLEG**2.0-1.0d0)
               ARGLEG3= 0.5d0*(5.0d0*X0TLEG**3.0-3.0d0*X0TLEG)
               nt = ntime(kt,kr)
               if(kr==1)then                                            !bulk
                  tcf_leg_1_b(nt)=tcf_leg_1_b(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_b(nt)=tcf_leg_2_b(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_b(nt)=tcf_leg_3_b(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==2)then                                        !hshell
                  tcf_leg_1_hs(nt)=tcf_leg_1_hs(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_hs(nt)=tcf_leg_2_hs(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_hs(nt)=tcf_leg_3_hs(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==3)then                                        !hshell-tetrahedral
                  tcf_leg_1_th(nt)=tcf_leg_1_th(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_th(nt)=tcf_leg_2_th(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_th(nt)=tcf_leg_3_th(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==4)then                                        !hShell-non-tetrahedral
                  tcf_leg_1_nth(nt)=tcf_leg_1_nth(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_nth(nt)=tcf_leg_2_nth(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_nth(nt)=tcf_leg_3_nth(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean tcfs over different origins
do kr = 1,NRSHMAX         !loop over environments
   do it = 1,NDELS
!nOHorg(kt,kr) - number of OH groups found in origin (block) kt and in the environment kr               
      if(kr==1)then
         tcf_mean_1(it,kr) = tcf_leg_1_b(it)
         tcf_mean_2(it,kr) = tcf_leg_2_b(it)
         tcf_mean_3(it,kr) = tcf_leg_3_b(it)
      elseif(kr==2)then   
         tcf_mean_1(it,kr) = tcf_leg_1_hs(it)
         tcf_mean_2(it,kr) = tcf_leg_2_hs(it)
         tcf_mean_3(it,kr) = tcf_leg_3_hs(it)
      elseif(kr==3)then
         tcf_mean_1(it,kr) = tcf_leg_1_th(it)
         tcf_mean_2(it,kr) = tcf_leg_2_th(it)
         tcf_mean_3(it,kr) = tcf_leg_3_th(it)
      elseif(kr==4)then
         tcf_mean_1(it,kr) = tcf_leg_1_nth(it)
         tcf_mean_2(it,kr) = tcf_leg_2_nth(it)
         tcf_mean_3(it,kr) = tcf_leg_3_nth(it)
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of OH groups on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nOHmean(kr) = nOHmean(kr) + nOHorg(kt,kr)
      if(nOHorg(kt,kr)>0)then
         koc(kr) = koc(kr) + 1                                    !number of origins where at least a single water was found in a given environment
      endif
   end do
end do

!Version 21
!standard deviation for the number of waters on each environment
do kr = 1,NRSHMAX
   do kt = 1, ko
      stdevnOH_sum(kr) = stdevnOH_sum(kr) + ( dfloat(nOHorg(kt,kr))/2.0 - dfloat(nOHmean(kr))/dfloat(ko*2) )**2.0 
   end do      
end do
!End Version 21


if(bulk_flag==0)kr_start=2                                        !do not print the bulk

do kr = kr_start,NRSHMAX
   if(ko > 1)stdevnOH(kr) = dsqrt(stdevnOH_sum(kr)/dfloat(ko-1))
   if(ko == 1)stdevnOH(kr) = dsqrt(stdevnOH_sum(kr)/dfloat(ko))
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a45,1x,i7)')'#Number of origins where waters were found = ',koc(kr)
   write(*,'(a28,1x,F14.2,1x,a7,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko),'stdev =',stdevnOH(kr) 
   write(*,'(a25,1x,F14.2,1x,a7,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2),'stdev =',stdevnOH(kr)  
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a45,1x,i7)')'#Number of origins where waters were found = ',koc(kr)
   write(n0,'(a28,1x,F14.2,1x,a7,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko),'stdev =',stdevnOH(kr) 
   write(n0,'(a25,1x,F14.2,1x,a7,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2),'stdev =',stdevnOH(kr)  
end do   
      
deallocate(tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b)
deallocate(tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs)
deallocate(tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th)
deallocate(tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth)
      
!Print mean orient. tcfs      
write(n2,9)
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a5)')'#Bulk'
   if(kr==2)write(n1,'(a4)')'#HSh'
   if(kr==3)write(n1,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n1,'(a15)')'#HSh-non-tetrah'
   if(kr==1)write(n2,'(a5)')'#Bulk'
   if(kr==2)write(n2,'(a4)')'#HSh'
   if(kr==3)write(n2,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n2,'(a15)')'#HSh-non-tetrah'
   if(kr==1)write(n3,'(a5)')'#Bulk'
   if(kr==2)write(n3,'(a4)')'#HSh'
   if(kr==3)write(n3,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n3,'(a15)')'#HSh-non-tetrah'
   do it = 1,NDELS
      jt = it -1                                                     !print tcf starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
!NG      WRITE(n1,19)PSECS,tcf_mean_1(it,kr)/dfloat(ko)     !mean orientational tcf rank 1 over different time-origins 
!NG      WRITE(n2,19)PSECS,tcf_mean_2(it,kr)/dfloat(ko)     !mean orientational tcf rank 2 over different time-origins   
!NG      WRITE(n3,19)PSECS,tcf_mean_3(it,kr)/dfloat(ko)     !mean orientational tcf rank 3 over different time-origins 
      WRITE(n1,19)PSECS,tcf_mean_1(it,kr)/dfloat(koc(kr))     !mean orientational tcf rank 1 over different time-origins 
      WRITE(n2,19)PSECS,tcf_mean_2(it,kr)/dfloat(koc(kr))     !mean orientational tcf rank 2 over different time-origins   
      WRITE(n3,19)PSECS,tcf_mean_3(it,kr)/dfloat(koc(kr))     !mean orientational tcf rank 3 over different time-origins 
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
end do

!Calculate the relaxation times - tau - trapezoid method
do kr = kr_start,NRSHMAX
   if(kr==1)then
      write(n11,'(a5)')'#Bulk'
      write(n12,'(a5)')'#Bulk'
      write(n13,'(a5)')'#Bulk'
   elseif(kr==2)then
      write(n11,'(a4)')'#HSh'
      write(n12,'(a4)')'#HSh'
      write(n13,'(a4)')'#HSh'
   elseif(kr==3)then
      write(n11,'(a11)')'#HSh-tetrah'
      write(n12,'(a11)')'#HSh-tetrah'
      write(n13,'(a11)')'#HSh-tetrah'
   elseif(kr==4)then
      write(n11,'(a15)')'#HSh-non-tetrah'
      write(n12,'(a15)')'#HSh-non-tetrah'
      write(n13,'(a15)')'#HSh-non-tetrah'
   endif  
   
!NG   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!NG   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!NG   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)   
         
   if(kr==1)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==2)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==3)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==4)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   endif   
end do 

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n11)
close(n12)
close(n13)

deallocate(NTDW4O,NONTDW)
deallocate(nOHorg,NOHidO,NOHidH)
deallocate(ntime)
deallocate(tcf_mean_1,tcf_mean_2,tcf_mean_3)
deallocate(norg_0HB,kv,nOHmean)
deallocate(stdevnOH_sum,stdevnOH)
deallocate(koc)

deallocate(xyz)
deallocate(x,y,z)


  9   FORMAT('#',17X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
      return

END SUBROUTINE wat_env_orient_tcfs



SUBROUTINE wat_tetra(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                     nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,rwtetbulk,tetsampw,ZMBIO,WATMODEL)               
! Calculate water's tetrahedrality in different environments - do not include solute heavy atoms as possible tetrahedron vertices

    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
    character(len=7),intent(in)                      :: WATMODEL
    character(len=4),dimension(natms)                :: atomname    
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment   
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),intent(in)                          :: rwtetbulk   ! Distance relative to solute com for onset of bulk water
    integer,intent(in)                               :: tetsampw
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
   
! Local variables
! Tetrahedrality and distance distributions
    real(kind=8)                                     :: SG,CHHASG,RO4
    integer                                          :: NWTDNB(4)
    REAL(KIND=4)                                     :: RNTDWH(2)
    real(kind=8)                                     :: RO_nb(5)
    real(kind=8)                                     :: RH2,RHO1
    integer                                          :: nsh,nsh_
    real(kind=8)                                     :: SGDEL,RDEL,ENRDEL
    integer,dimension(:,:),allocatable               :: NPDSG,NDPRO4,NDPRH2,NDPRHO1
    integer,dimension(:,:,:),allocatable             :: NDPRO
    integer,dimension(:,:,:),allocatable             :: NDPEO
    real(kind=8),dimension(:),allocatable            :: AVWSG,AVCHHASG,AVRO4,AVRH2,AVRHO1
    real(kind=8),dimension(:,:),allocatable          :: AVRO
    real(kind=8),dimension(:,:),allocatable          :: AVEPO
    integer,parameter                                :: nshdist = 50000                      !max. number of bins in tetrahedrality and distance dist.
    integer,parameter                                :: nshener = 50000                     !max. number of bins in pair potential energy
    integer                                          :: NTDELS,NRDELS
    integer                                          :: NENRDELS
    real(kind=8)                                     :: TETHQ,RD_ANG
    real(kind=8)                                     :: TETHDIST,DISTR_R
    real(kind=8)                                     :: E_kjmol,DISTR_E
! 
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms) 
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8)                                     :: dx,dy,dz,dr,dr2
    integer,dimension(:),allocatable                 :: kv
    integer,dimension(:,:),allocatable               :: kwlabel      !index of waters in each environment 
    integer      :: n0, n1, n2, n3,n4, n5, n6, n7, n8, n9
    integer      :: n00, n10, n11, n12, n13, n14, n20, n21
    integer      :: i, j, k
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: iH, ihat, nt, ko, jt
    integer      :: io, ilw
    integer      :: iv1, iv2, iv3, iv4
    integer      :: L, KN, IT, NKR
    integer,parameter                                 :: INTMAX = 2147483647
    integer,dimension(:),allocatable                  :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable                  :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                                      :: dr_bulk,dr4nb
    real(kind=8)                                      :: RADMAX,RADMIN
    integer                                           :: NRSHMAX   
    integer                                           :: nradsh    
    integer                                           :: natcheck
!    integer                                           :: kor
    integer                                           :: bulk_flag,kr_start
    integer                                           :: at_id                    ! solute atomic index to analyse the respective solvation environment
    integer                                           :: kstps
    logical                                           :: filex
    
    integer                                           :: maxwat_bulk                !maximum number of waters in the bulk
    integer                                           :: nw_bulk                    !Bulk water counter
!version 31
    real(kind=8)                                      :: ENvdw_nb(5),ENcoul_nb(5),ENpot_nb(5)        !water-water interaction energy with the 5 nearest nb
    logical                                           :: lprnt_teth
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Tetra_env_out/log_env_tetra.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Tetra_env_out')

inquire(file='Pair_En_env_out/log_env_pair_en.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Pair_En_env_out')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140
n5 = 150
n6 = 160
n7 = 170
n8 = 180
n9 = 190
n00 = 195
n10 = 200
n11 = 210
n12 = 220
n13 = 230
n14 = 240
n20 = 300
n21 = 310

open(n0,file='Tetra_env_out/log_env_tetra.dat',status='unknown',action='write')
open(n1,file='Tetra_env_out/wat_tetra_dist.dat',status='unknown',action='write')
open(n2,file='Tetra_env_out/wat_rO4_dist.dat',status='unknown',action='write')
open(n3,file='Tetra_env_out/wat_rOH2_dist.dat',status='unknown',action='write')
open(n4,file='Tetra_env_out/wat_rHO1_dist.dat',status='unknown',action='write')
open(n5,file='Tetra_env_out/wat_rO1_nb_dist.dat',status='unknown',action='write')
open(n6,file='Tetra_env_out/wat_rO2_nb_dist.dat',status='unknown',action='write')
open(n7,file='Tetra_env_out/wat_rO3_nb_dist.dat',status='unknown',action='write')
open(n8,file='Tetra_env_out/wat_rO4_nb_dist.dat',status='unknown',action='write')
open(n9,file='Tetra_env_out/wat_rO5_nb_dist.dat',status='unknown',action='write')
open(n20,file='Tetra_env_out/Label_tetrah.dat',status='unknown',action='write')
open(n21,file='Tetra_env_out/tetrahedrons.gro',status='unknown',action='write')
!
open(n00,file='Pair_En_env_out/log_env_pair_en.dat',status='unknown',action='write')
open(n10,file='Pair_En_env_out/wat_EpO1_nb_dist.dat',status='unknown',action='write')
open(n11,file='Pair_En_env_out/wat_EpO2_nb_dist.dat',status='unknown',action='write')
open(n12,file='Pair_En_env_out/wat_EpO3_nb_dist.dat',status='unknown',action='write')
open(n13,file='Pair_En_env_out/wat_EpO4_nb_dist.dat',status='unknown',action='write')
open(n14,file='Pair_En_env_out/wat_EpO5_nb_dist.dat',status='unknown',action='write')    
    
    
dr_bulk = rwtetbulk                              !set distance from the solute com to assume bulk water
bulk_flag = 1                                    !default - calculate tetrahedrality for bulk water
kr_start = 1                                     !environments loop started - if bulk_flag = 0; kr_start = 2 (skip environment bulk) 
lprnt_teth = .true.

!VMD configuration .gro with solute and water tetrahedrons
if (lprnt_teth) then
   write(*,*)
   write(*,*)'WARNING!!!'
   write(*,*)'PRINT WATER SHELL TETRAHEDRONS FOR A SINGLE CONF - VISUAL PURPOSES'
   write(*,*)'SOLUTE SHOULD BE CENTERED IN THE BOX'
   write(*,*)'gmx trjconv -f npt.xtc -s npt.tpr -center -pbc mol -o npt_pbc_mol.xtc'
   write(*,*)
   write(n0,*)
   write(n0,*)'WARNING!!!'
   write(n0,*)'PRINT WATER SHELL TETRAHEDRONS FOR A SINGLE CONF - VISUAL PURPOSES'
   write(n0,*)'SOLUTE SHOULD BE CENTERED IN THE BOX'
   write(n0,*)'gmx trjconv -f npt.xtc -s npt.tpr -center -pbc mol -o npt_pbc_mol.xtc'
   write(n0,*)
   write(n20,*)
   write(n20,'(8(2x,A6))')'STEP','CENT_W','VERT1','VERT2','VERT3','VERT4','INDEX','TetraH'
   write(n21,'(a39)')'Sol+Water - add nb atoms of water below'
   write(n21,'(i4)')nmolsol*natmsol
   ilw = nmolsol*natmsol
endif

maxwat_bulk = 500
if(nmolwat<maxwat_bulk)maxwat_bulk = nmolwat

if(dr_bulk==0.0)then
   bulk_flag = 0
   write(*,*)
   write(*,*)'Warning!!! No bulk region defined for tetrahedrality calculation'
   write(*,*)'Bulk tetrahedrality will not be calculated accordingly'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! No bulk region defined for tetrahedrality calculation'
   write(n0,*)'No bulk tetrahedrality will be calculated accordingly'
   write(n0,*)
   write(n00,*)
   write(n00,*)'Warning!!! No bulk region defined for tetrahedrality calculation'
   write(n00,*)'No bulk tetrahedrality will be calculated accordingly'
   write(n00,*)
else 
   write(*,*)
   write(*,*)'Warning!!! Maximum bulk waters for tetrahedrality calculation = ',maxwat_bulk
   write(*,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! Maximum bulk waters for tetrahedrality calculation = ',maxwat_bulk
   write(n0,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(n0,*)
   write(n00,*)
   write(n00,*)'Warning!!! Maximum bulk waters for tetrahedrality calculation = ',maxwat_bulk
   write(n00,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(n00,*)   
endif   

NRSHMAX = 4                                      !maximum number of environments
SGDEL   = 0.01D0                                 !tetrahedrality bin size 
RDEL    = 0.01D0                                 !bin size in distance distributions Ang
NRDELS  = 10.0d0/RDEL                            !maximum number of bins in distance distributions - set to 10 Ang
ENRDEL   = 0.4D0                                 !bin size in energy distributions kJ/mol
NENRDELS = 1000/ENRDEL                           !maximum number of bins in energy distributions

if(NRDELS>nshdist)then
   write(*,*)'error - routine wat_tetra - too many bins'
   stop
endif
if(NENRDELS>nshener)then
   write(*,*)'error - routine wat_tetra - too many bins'
   stop
endif
if(-NENRDELS<-nshener)then
   write(*,*)'error - routine wat_tetra - too many bins'
   stop
endif 

!Tetrahedrality and distances
allocate(NPDSG(NRSHMAX,nshdist),NDPRO4(NRSHMAX,nshdist),NDPRH2(NRSHMAX,nshdist),NDPRHO1(NRSHMAX,nshdist))
allocate(NDPRO(NRSHMAX,5,nshdist),NDPEO(NRSHMAX,5,-nshener:nshener))
allocate(AVWSG(NRSHMAX),AVCHHASG(NRSHMAX),AVRO4(NRSHMAX),AVRH2(NRSHMAX),AVRHO1(NRSHMAX))
allocate(AVRO(NRSHMAX,5),AVEPO(NRSHMAX,5))
NPDSG=0;NDPRO4=0;NDPRH2=0;NDPRHO1=0
NDPRO=0
NDPEO=0
AVWSG=0.0d0;AVCHHASG=0.0d0;AVRO4=0.0d0;AVRH2=0.0d0;AVRHO1=0.0d0
AVRO=0.0d0
AVEPO=0.0d0

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))   
allocate(kv(NRSHMAX))
allocate(kwlabel(NRSHMAX,nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled)    
NTDW4O = 0; NONTDW = 0
kv = 0
kstps = 0
kwlabel = 0

 cmxsol = 0.0d0; cmysol=0.0d0; cmzsol=0.0d0

! Start tetrahedrality calculation

write(*,*)
write(*,*)'Starting water tetrahedrality calculation...'
write(*,*)'Routine wat_tetra'
write(*,*)'Results are printed to Tetra_env_out'
write(*,*)'water sampling frequency (steps) = ',tetsampw
write(*,*)'total number of time-steps =',nstep
write(*,*)

write(n0,*)
write(n0,*)'Starting water tetrahedrality calculation...'
write(n0,*)'Routine wat_tetra'
write(n0,*)'Results are printed to Tetra_env_out'
write(n0,*)'water sampling frequency (steps) = ',tetsampw
write(n0,*)'total number of time-steps =',nstep
write(n0,*)

write(n00,*)
write(n00,*)'Starting water pair interaction energy calculation...'
write(n00,*)'Routine wat_tetra'
write(n00,*)'Results are printed to Pair_En_env_out'
write(n00,*)'water sampling frequency (steps) = ',tetsampw
write(n00,*)'total number of time-steps =',nstep
write(n00,*)

do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
   NTDW4O = 0; NONTDW = 0
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

   if (j==tetsampw .and. lprnt_teth)then
      do i=1,nmolsol*natmsol       
         write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'1MOL',trim(adjustl(atomname(i))),i,0.1*x(i),0.1*y(i),0.1*z(i)
      end do
   endif   
!   if (j==tetsampw .and. lprnt_teth)write(n21,'(3(F9.4))')0.1*cell

!sample for water molecules in specific environments


! sample frequency - calculate the tetrahedrality every tetsampw time-step  

   if(j==1.or.mod(j,tetsampw)==0)then
      nw_bulk = 0                                                               !bulk water counter
      write(*,*)'time-step',j,' sampling water molecules'
      kstps = kstps + 1                                                         !count number of origins sampled
      kwlabel = 0                                                                 !array with waters ID to avoid repetition
      
!Find solute centre of mass for bulk water tetrahedrality         
      if(bulk_flag==1)call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water      
         dx = x(i) - cmxsol
         dy = y(i) - cmysol
         dz = z(i) - cmzsol
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
         if(bulk_flag==0)dr_bulk = 100.d0*cell(1)                                        !define a distance larger than the box "radius"
!Version 13         if(dr>=dr_bulk)then                                                             !bulk water
         if(dr>=dr_bulk.and.nw_bulk<maxwat_bulk)then
            nw_bulk = nw_bulk + 1
            if(bulk_flag==0)write(*,*)'Warning!!! Found a water molecule in the bulk'    !bulk water should not be found if bulk_flag=0                                                 
            nradsh = 1                                                                   !hydration shell environment
            kv(nradsh) = kv(nradsh) + 1                                                  !number of waters
            call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                             SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH) 
!version 31
            call water_water_En(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,WATMODEL,ENvdw_nb,ENcoul_nb,ENpot_nb) 
! distributions amd means
!            nsh = IDNINT(3.0D0/SGDEL+SG/SGDEL+0.5D0)
            nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0
            NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh)+1
            AVWSG(nradsh) = AVWSG(nradsh) + SG
            AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
            SG = 0.0D0
! CALCULATE THE RO4 DISTRIBUTION AND AVERAGE VALUE
            nsh = RO4/RDEL+0.5D0
            NDPRO4(nradsh,nsh) = NDPRO4(nradsh,nsh) + 1
            AVRO4(nradsh) = AVRO4(nradsh) + RO4
            RO4 = 0.0D0
! CALCULATE THE RH2 DISTRIBUTION AND AVERAGE VALUE
            nsh = RH2/RDEL+0.5D0
            NDPRH2(nradsh,nsh) = NDPRH2(nradsh,nsh) + 1
            AVRH2(nradsh) = AVRH2(nradsh) + RH2
            RH2 = 0.0D0
! CALCULATE THE RHO1 DISTRIBUTION AND AVERAGE VALUE
            nsh = RHO1/RDEL+0.5D0
            NDPRHO1(nradsh,nsh) = NDPRHO1(nradsh,nsh) + 1
            AVRHO1(nradsh) = AVRHO1(nradsh) + RHO1
            RHO1 = 0.0D0
! CALCULATE THE RO DISTRIBUTIONS AND AVERAGE VALUES
            do KN = 1,5                                                          !run over five distinct oxygen-oxygen distances
               nsh = RO_nb(KN)/RDEL+0.5D0
               NDPRO(nradsh,KN,nsh) = NDPRO(nradsh,KN,nsh) + 1
               AVRO(nradsh,KN) = AVRO(nradsh,KN) + RO_nb(KN)
               RO_nb(KN) = 0.0D0       
            end do       
! CALCULATE THE POTENTIAL ENERGY DISTRIBUTIONS AND AVERAGE VALUES
            do KN = 1,5                                                          !run over five distinct neighbors
!              nsh_ = ENpot_nb(KN)/ENRDEL+0.5D0                Warning !!! this is wrong for negative and positive bins
               nsh = IDNINT(ENpot_nb(KN)/ENRDEL)
!check               write(*,*)nsh,nsh_
               NDPEO(nradsh,KN,nsh) = NDPEO(nradsh,KN,nsh) + 1
               AVEPO(nradsh,KN) = AVEPO(nradsh,KN) + ENpot_nb(KN)
!check               write(*,*)nsh,KN,ENpot_nb(KN)
               ENpot_nb(KN) = 0.0D0
            end do
!Electrostatic            
!Coul            do KN = 1,5                                                          !run over five distinct neighbors
!Coul               nsh = ENcoul_nb(KN)/ENRDEL+0.5D0
!Coul               NDPEO_coul(nradsh,KN,nsh) = NDPEO_coul(nradsh,KN,nsh) + 1
!Coul               AVEPO_coul(nradsh,KN) = AVEPO_coul(nradsh,KN) + ENcoul_nb(KN)
!Coul               ENcoul_nb(KN) = 0.0D0
!Coul            end do            
!van der Waals
!vdW            do KN = 1,5                                                          !run over five distinct neighbors
!vdW               nsh = ENvdw_nb(KN)/ENRDEL+0.5D0
!vdW               NDPEO_vdW(nradsh,KN,nsh) = NDPEO_vdW(nradsh,KN,nsh) + 1
!vdW               AVEPO_vdW(nradsh,KN) = AVEPO_vdW(nradsh,KN) + ENvdw_nb(KN)
!vdW               ENvdw_nb(KN) = 0.0D0
!vdW            end do                               
! End distributions and means                     
         endif                                                                   !end bulk

!Loop over atomic species whose solvation environment will be studied      

         do at_id=1,natHSh                                   !Loop over protein atomic species to analyze

!Hydration shell - tetrahedral and non tetrahedral water in the first HSh         
            k = natsol_id(at_id)
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)              
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then                            !water within hydration shell radius
               nradsh = 2                                                                    !hydration shell environment
!Avoid water molecules repetition
               if(kwlabel(nradsh,io)==0)then
                  kwlabel(nradsh,io)=1                                                       !Oxygen id for correction of the nb of waters
                  kv(nradsh) = kv(nradsh) + 1                                                !number of waters in the HSh (nradsh = 2) 
               else
!NG                  write(*,*)'water molecule already sampled ',io,' = ',i
                  goto 10                                                                    !sample another water molecule
               endif   
!End avoid water molecules repitition                             
               call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)
!version 31
               call water_water_En(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,WATMODEL,ENvdw_nb,ENcoul_nb,ENpot_nb)                               
! distributions amd means
               nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0                                        
               NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh) + 1
               AVWSG(nradsh) = AVWSG(nradsh) + SG
               AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
               SG = 0.0D0
! CALCULATE THE RO4 DISTRIBUTION AND AVERAGE VALUE
               nsh = RO4/RDEL+0.5D0
               NDPRO4(nradsh,nsh) = NDPRO4(nradsh,nsh) + 1
               AVRO4(nradsh) = AVRO4(nradsh) + RO4
               RO4 = 0.0D0
! CALCULATE THE RH2 DISTRIBUTION AND AVERAGE VALUE
               nsh = RH2/RDEL+0.5D0
               NDPRH2(nradsh,nsh) = NDPRH2(nradsh,nsh) + 1
               AVRH2(nradsh) = AVRH2(nradsh) + RH2
               RH2 = 0.0D0
! CALCULATE THE RHO1 DISTRIBUTION AND AVERAGE VALUE
               nsh = RHO1/RDEL+0.5D0
               NDPRHO1(nradsh,nsh) = NDPRHO1(nradsh,nsh) + 1
               AVRHO1(nradsh) = AVRHO1(nradsh) + RHO1
               RHO1 = 0.0D0
! CALCULATE THE RO DISTRIBUTIONS AND AVERAGE VALUES
               do KN = 1,5
                  nsh = RO_nb(KN)/RDEL+0.5D0
                  NDPRO(nradsh,KN,nsh) = NDPRO(nradsh,KN,nsh) + 1
                  AVRO(nradsh,KN) = AVRO(nradsh,KN) + RO_nb(KN)
                  RO_nb(KN) = 0.0D0       
               end do
! CALCULATE THE POTENTIAL ENERGY DISTRIBUTIONS AND AVERAGE VALUES
               do KN = 1,5                                                          !run over five distinct neighbors
                  nsh = IDNINT(ENpot_nb(KN)/ENRDEL)
                  NDPEO(nradsh,KN,nsh) = NDPEO(nradsh,KN,nsh) + 1
                  AVEPO(nradsh,KN) = AVEPO(nradsh,KN) + ENpot_nb(KN)
                  ENpot_nb(KN) = 0.0D0
               end do                                  
! End distributions and means          

!tetrahedral and non-tetrahedral waters in the hydration shell for sub-ensemble tetrahedral and non-tetrahedral distributions      
! FIND THE FOUR NEAREST O ATOMS J TO EACH OXYGEN I
               RADMAX  = 25.0D0
               RADMIN  =  0.0D0        
               JT      = 1
               DO WHILE(JT.LE.4) 
                  do kw = nmolsol*natmsol+1,natms-nions,nwatsites                         !Loop over water oxygens
                     if(kw.ne.i)then
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOMS i AND OXYGEN ATOMS kw  
                        dx = x(i) - x(kw)
                        dy = y(i) - y(kw)
                        dz = z(i) - z(kw)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        IF((dr>RADMIN).AND.(dr<RADMAX))THEN
                           RADMAX = dr
!NTDW4O STORES THE O ATOMS NEIGHBORS OF i: JT = 1 FIRST NEIGHBOR; JT = 2 SECOND NEIGHBOR; JT = 3 THIRD NEIGHBOR; JT = 4 FOURTH NEIGHBOR
                           NTDW4O(JT) = kw
                        ENDIF
                     endif
                  end do
                  JT = JT + 1
                  RADMIN = RADMAX
                  RADMAX = 25.0D0
               END DO
!CHECK IF SOLUTE ATOMS ARE CLOSER THAN THE 4th NEAREST O ATOM          
               KI=NTDW4O(4)                                                     !Oxygen ID of the 4th nearest water neighbor 
               dx = x(i) - x(KI)
               dy = y(i) - y(KI)
               dz = z(i) - z(KI)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               dr4nb = dr                                                           !distance to the fourth nearest neighbor
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY SOLUTE ATOM IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = 1                                                                !run over all solute atoms - heavy and non-heavy atoms
               DO WHILE(L.LE.nmolsol*natmsol)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN             
                     NONTDW(io) = 1                                               !water i is non-tetrahedral
                     goto 5       
                  ENDIF
                  L = L+1
               END DO     
!               5 continue
               
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY ION IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = nmolsol*natmsol + nmolwat*nwatsites + 1                                                                !run over all ions
               DO WHILE(L.LE.natms)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN   
                     write(*,*)'Water molecule next to ionic species'
                     NONTDW(io) = 1                                               !water i is non-tetrahedral
                     goto 5       
                  ENDIF
                  L = L+1
               END DO     
               5 continue  
               
!Non-tetrahedral waters (less than 4 water neighbors)               
               if(NONTDW(io)==1)then                                                            !Non-tetrahedral water    
                  nradsh = 4                                                                    !hydration shell (non-tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                                                   !number of waters
                  call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                   SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)
!version 31
                  call water_water_En(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,WATMODEL,ENvdw_nb,ENcoul_nb,ENpot_nb)                                
! distributions amd means
                  nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0                                        
                  NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh) + 1
                  AVWSG(nradsh) = AVWSG(nradsh) + SG
                  AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
                  SG = 0.0D0
! CALCULATE THE RO4 DISTRIBUTION AND AVERAGE VALUE
                  nsh = RO4/RDEL+0.5D0
                  NDPRO4(nradsh,nsh) = NDPRO4(nradsh,nsh) + 1
                  AVRO4(nradsh) = AVRO4(nradsh) + RO4
                  RO4 = 0.0D0
! CALCULATE THE RH2 DISTRIBUTION AND AVERAGE VALUE
                  nsh = RH2/RDEL+0.5D0
                  NDPRH2(nradsh,nsh) = NDPRH2(nradsh,nsh) + 1
                  AVRH2(nradsh) = AVRH2(nradsh) + RH2
                  RH2 = 0.0D0
! CALCULATE THE RHO1 DISTRIBUTION AND AVERAGE VALUE
                  nsh = RHO1/RDEL+0.5D0
                  NDPRHO1(nradsh,nsh) = NDPRHO1(nradsh,nsh) + 1
                  AVRHO1(nradsh) = AVRHO1(nradsh) + RHO1
                  RHO1 = 0.0D0
! CALCULATE THE RO DISTRIBUTIONS AND AVERAGE VALUES
                  do KN = 1,5
                     nsh = RO_nb(KN)/RDEL+0.5D0
                     NDPRO(nradsh,KN,nsh) = NDPRO(nradsh,KN,nsh) + 1
                     AVRO(nradsh,KN) = AVRO(nradsh,KN) + RO_nb(KN)
                     RO_nb(KN) = 0.0D0       
                  end do
! CALCULATE THE POTENTIAL ENERGY DISTRIBUTIONS AND AVERAGE VALUES
                  do KN = 1,5                                                          !run over five distinct neighbors
                     nsh = IDNINT(ENpot_nb(KN)/ENRDEL)
                     NDPEO(nradsh,KN,nsh) = NDPEO(nradsh,KN,nsh) + 1
                     AVEPO(nradsh,KN) = AVEPO(nradsh,KN) + ENpot_nb(KN)
                     ENpot_nb(KN) = 0.0D0
                  end do                                                    
! End distributions and means          

               elseif(NONTDW(io)==0)then                                                        !Tetrahedral water 
                  nradsh = 3                                                                    !hydration shell (tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                                                   !number of waters
                  call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                   SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)
!version 31
                  call water_water_En(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,WATMODEL,ENvdw_nb,ENcoul_nb,ENpot_nb)                               
!check 0k            write(*,*)j,i,sg 
!Version 35  
                  if(j == tetsampw .and. lprnt_teth)then
!                     write(*,*)'Printing Tetrahedrons...'
                     ilw = ilw + 1
                     write(n20,'(7(2x,i5),3x,F7.3)')j,i,NWTDNB(1),NWTDNB(2),NWTDNB(3),NWTDNB(4),ilw,SG
! Print central water molecules
! FOR SEQUENCIAL NUMBERS PRINT ilw INSTEAD OF i - HOWEVER THE FILE n20 ALLOWS IDENTIFYING TETRAHEDRONS BY i
!                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i),ilw,0.1*x(i),0.1*y(i),0.1*z(i)
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i),ilw,0.1*x(i),0.1*y(i),0.1*z(i)
                     ilw = ilw + 1                                                                                             
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i+1),ilw,0.1*x(i+1),0.1*y(i+1),0.1*z(i+1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'2SOL',atomname(i+2),ilw,0.1*x(i+2),0.1*y(i+2),0.1*z(i+2)       
!Print tetrahedron vertices 
!WARNING!!! THERE WILL BE SOME REPETITION - SOME CENTRAL WATERS ARE ALSO VERTICES
                     iv1 = NWTDNB(1)
                     iv2 = NWTDNB(2)
                     iv3 = NWTDNB(3)
                     iv4 = NWTDNB(4)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'3SOL',atomname(iv1),ilw,0.1*x(iv1),0.1*y(iv1),0.1*z(iv1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'3SOL',atomname(iv1+1),ilw,0.1*x(iv1+1),0.1*y(iv1+1),0.1*z(iv1+1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'3SOL',atomname(iv1+2),ilw,0.1*x(iv1+2),0.1*y(iv1+2),0.1*z(iv1+2) 
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'4SOL',atomname(iv2),ilw,0.1*x(iv2),0.1*y(iv2),0.1*z(iv2)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'4SOL',atomname(iv2+1),ilw,0.1*x(iv2+1),0.1*y(iv2+1),0.1*z(iv2+1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'4SOL',atomname(iv2+2),ilw,0.1*x(iv2+2),0.1*y(iv2+2),0.1*z(iv2+2) 
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'5SOL',atomname(iv3),ilw,0.1*x(iv3),0.1*y(iv3),0.1*z(iv3)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'5SOL',atomname(iv3+1),ilw,0.1*x(iv3+1),0.1*y(iv3+1),0.1*z(iv3+1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'5SOL',atomname(iv3+2),ilw,0.1*x(iv3+2),0.1*y(iv3+2),0.1*z(iv3+2) 
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'6SOL',atomname(iv4),ilw,0.1*x(iv4),0.1*y(iv4),0.1*z(iv4)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'6SOL',atomname(iv4+1),ilw,0.1*x(iv4+1),0.1*y(iv4+1),0.1*z(iv4+1)
                     ilw = ilw + 1
                     write(n21,'(4x,A4,4x,A3,I5,3(1x,F7.3))')'6SOL',atomname(iv4+2),ilw,0.1*x(iv4+2),0.1*y(iv4+2),0.1*z(iv4+2) 
                  endif   
!End version 35

! distributions amd means
                  nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0                                        
                  NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh) + 1
                  AVWSG(nradsh) = AVWSG(nradsh) + SG
                  AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
                  SG = 0.0D0
! CALCULATE THE RO4 DISTRIBUTION AND AVERAGE VALUE
                  nsh = RO4/RDEL+0.5D0
                  NDPRO4(nradsh,nsh) = NDPRO4(nradsh,nsh) + 1
                  AVRO4(nradsh) = AVRO4(nradsh) + RO4
                  RO4 = 0.0D0
! CALCULATE THE RH2 DISTRIBUTION AND AVERAGE VALUE
                  nsh = RH2/RDEL+0.5D0
                  NDPRH2(nradsh,nsh) = NDPRH2(nradsh,nsh) + 1
                  AVRH2(nradsh) = AVRH2(nradsh) + RH2
                  RH2 = 0.0D0
! CALCULATE THE RHO1 DISTRIBUTION AND AVERAGE VALUE
                  nsh = RHO1/RDEL+0.5D0
                  NDPRHO1(nradsh,nsh) = NDPRHO1(nradsh,nsh) + 1
                  AVRHO1(nradsh) = AVRHO1(nradsh) + RHO1
                  RHO1 = 0.0D0
! CALCULATE THE RO DISTRIBUTIONS AND AVERAGE VALUES
                  do KN = 1,5
                     nsh = RO_nb(KN)/RDEL+0.5D0
                     NDPRO(nradsh,KN,nsh) = NDPRO(nradsh,KN,nsh) + 1
                     AVRO(nradsh,KN) = AVRO(nradsh,KN) + RO_nb(KN)
                     RO_nb(KN) = 0.0D0       
                  end do
! CALCULATE THE POTENTIAL ENERGY DISTRIBUTIONS AND AVERAGE VALUES
                  do KN = 1,5                                                          !run over five distinct neighbors
                     nsh = IDNINT(ENpot_nb(KN)/ENRDEL)
                     NDPEO(nradsh,KN,nsh) = NDPEO(nradsh,KN,nsh) + 1
                     AVEPO(nradsh,KN) = AVEPO(nradsh,KN) + ENpot_nb(KN)
                     ENpot_nb(KN) = 0.0D0
                  end do                                                                      
! End distributions and means          

               endif                                                                !end water non-tetrahedral/tetrahedral            
            endif                                                                   !end water is in the coordination sphere 
         end do                                                                     !end loop over atomic species
         10 continue
      end do                                                                        !end loop over water oxygens
   endif                                                                            !end sampling of waters
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
!   lprnt_teth = .false.                                                             !version 35
if (j==tetsampw .and. lprnt_teth)write(n21,'(3(F9.4))')0.1*cell
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                                                                                                                                        

!Average number of waters on each environment

write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps

do kr = 1,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
!   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep) 
   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
!   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep)
   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps)   
   write(n00,*)
   write(n00,'(a15,1x,i4)')'#Environment = ',kr 
   write(n00,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps)
end do   

!Print distributions and mean values 

if(bulk_flag==0)kr_start=2                                        !do not print the bulk

write(n0,*)
write(n0,*)'TETRAHEDRALITY AND DISTANCE DISTRIBUTIONS'
write(n0,*)'========================================='
write(n0,*)

write(n00,*)
write(n00,*)'PAIR INTERACTION ENERGY DISTRIBUTIONS'
write(n00,*)'====================================='
write(n00,*)
     
do kr = kr_start,NRSHMAX                                                      !loop over environments
   IF(kr==1)WRITE(n0,*) 
   IF(kr==1)WRITE(n0,*)'Bulk'
   IF(kr==1)WRITE(n0,*)'****'
   IF(kr==2)WRITE(n0,*)
   IF(kr==2)WRITE(n0,*)'Hydration Shell'
   IF(kr==2)WRITE(n0,*)'***************'
   IF(kr==3)WRITE(n0,*)
   IF(kr==3)WRITE(n0,*)'HSh-tetrahedral'
   IF(kr==3)WRITE(n0,*)'***************'
   IF(kr==4)WRITE(n0,*)
   IF(kr==4)WRITE(n0,*)'HSh-non-tetrahedral' 
   IF(kr==4)WRITE(n0,*)'*******************' 
!
   IF(kr==1)WRITE(n00,*) 
   IF(kr==1)WRITE(n00,*)'Bulk'
   IF(kr==1)WRITE(n00,*)'****'
   IF(kr==2)WRITE(n00,*)
   IF(kr==2)WRITE(n00,*)'Hydration Shell'
   IF(kr==2)WRITE(n00,*)'***************'
   IF(kr==3)WRITE(n00,*)
   IF(kr==3)WRITE(n00,*)'HSh-tetrahedral'
   IF(kr==3)WRITE(n00,*)'***************'
   IF(kr==4)WRITE(n00,*)
   IF(kr==4)WRITE(n00,*)'HSh-non-tetrahedral' 
   IF(kr==4)WRITE(n00,*)'*******************'    
   
! q tetrahedrality
   IF(kr==1)write(n1,'(a5)')'#Bulk'
   IF(kr==2)write(n1,'(a4)')'#HSh'
   IF(kr==3)write(n1,'(a11)')'#HSh-tetrah'
   IF(kr==4)write(n1,'(a15)')'#HSh-non-tetrah'   
! q VARIES BETWEEN -3 AND 1 FOR THE DEBENEDETTI NORMALIZATION IMPLEMENTED HERE
   NTDELS = 4.0D0/SGDEL
   DO IT = 1,NTDELS
      IF(NPDSG(kr,IT).GT.0)THEN
         TETHQ    = -3.01D0+SGDEL*DFLOAT(IT)
         TETHDIST = DFLOAT(NPDSG(kr,IT))/DFLOAT(kv(kr))                    !kv is the sum of the number of waters over the number of origins 
         WRITE(n1,9)TETHQ,TETHDIST
      ENDIF
   END DO
   WRITE(n1,*)
!
   AVWSG(kr)    = AVWSG(kr)/DFLOAT(kv(kr))
   AVCHHASG(kr) = AVCHHASG(kr)/DFLOAT(kv(kr))
   WRITE(n0,*)
   WRITE(n0,*)
   WRITE(n0,*)'Errington & Debenedetti perfect tetrahedron <q> = 1'
   WRITE(n0,*)'Average Errington & Debenedetti <q> = ',AVWSG(kr)
   WRITE(n0,*)
   WRITE(n0,*)'Chau & Hardwick perfect tetrahedron <q> = 0 '
   WRITE(n0,*)'Average Chau & Hardwick <Sg> =  <q> = ',AVCHHASG(kr)
   WRITE(n0,*)
   
! rO4 - mean over the 4 nearest oxygens 
   IF(kr==1)write(n2,'(a5)')'#Bulk'
   IF(kr==2)write(n2,'(a4)')'#HSh'
   IF(kr==3)write(n2,'(a11)')'#HSh-tetrah'
   IF(kr==4)write(n2,'(a15)')'#HSh-non-tetrah'    
   DO NKR =1,NRDELS
      IF(NDPRO4(kr,NKR).GT.0)THEN
         RD_ANG = RDEL*DFLOAT(NKR)          
         DISTR_R=DFLOAT(NDPRO4(kr,NKR))/DFLOAT(kv(kr))
         WRITE(n2,9)RD_ANG,DISTR_R
      ENDIF
   END DO
   WRITE(n2,*)
!
   AVRO4(kr) = AVRO4(kr)/DFLOAT(kv(kr))
   WRITE(n0,*)
   WRITE(n0,*)
   WRITE(n0,*)'AVERAGE <RO_O4> Ang = ',AVRO4(kr)
   WRITE(n0,*)
   
! rOH2 = rH2 - mean over the 2 nearest hydrogens 
   IF(kr==1)write(n3,'(a5)')'#Bulk'
   IF(kr==2)write(n3,'(a4)')'#HSh'
   IF(kr==3)write(n3,'(a11)')'#HSh-tetrah'
   IF(kr==4)write(n3,'(a15)')'#HSh-non-tetrah'    
   DO NKR =1,NRDELS
      IF(NDPRH2(kr,NKR).GT.0)THEN
         RD_ANG = RDEL*DFLOAT(NKR)          
         DISTR_R=DFLOAT(NDPRH2(kr,NKR))/DFLOAT(kv(kr))
         WRITE(n3,9)RD_ANG,DISTR_R
      ENDIF
   END DO
   WRITE(n3,*)
!
   AVRH2(kr) = AVRH2(kr)/DFLOAT(kv(kr))
   WRITE(n0,*)
   WRITE(n0,*)
   WRITE(n0,*)'AVERAGE <RO_H2> Ang = ',AVRH2(kr)
   WRITE(n0,*)

! rHO1 - mean over the nearest oxygen of each hydrogen 
   IF(kr==1)write(n4,'(a5)')'#Bulk'
   IF(kr==2)write(n4,'(a4)')'#HSh'
   IF(kr==3)write(n4,'(a11)')'#HSh-tetrah'
   IF(kr==4)write(n4,'(a15)')'#HSh-non-tetrah'    
   DO NKR =1,NRDELS
      IF(NDPRHO1(kr,NKR).GT.0)THEN
         RD_ANG = RDEL*DFLOAT(NKR)          
         DISTR_R=DFLOAT(NDPRHO1(kr,NKR))/DFLOAT(kv(kr))
         WRITE(n4,9)RD_ANG,DISTR_R
      ENDIF
   END DO
   WRITE(n4,*)
!
   AVRHO1(kr) = AVRHO1(kr)/DFLOAT(kv(kr))
   WRITE(n0,*)
   WRITE(n0,*)
   WRITE(n0,*)'AVERAGE <RH_O1> Ang = ',AVRHO1(kr)
   WRITE(n0,*)

!rOi [i = 1,2,3,4,5] - mean over the nearest, the second nearest, the third nearest...etc oxygen
   do KN = 1,5
      IF(kr==1)write(n4+10*KN,'(a5)')'#Bulk'
      IF(kr==2)write(n4+10*KN,'(a4)')'#HSh'
      IF(kr==3)write(n4+10*KN,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n4+10*KN,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRO(kr,KN,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRO(kr,KN,NKR))/DFLOAT(kv(kr))
            WRITE(n4+10*KN,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n4+10*KN,*)
!
      AVRO(kr,KN) = AVRO(kr,KN)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RO_O',KN,'> Ang = ',AVRO(kr,KN)
      WRITE(n0,*)                    
   end do
!
!EpOi [i = 1,2,3,4,5] - mean over the nearest, the second nearest, the third nearest...etc oxygen
   do KN = 1,5
      IF(kr==1)write(n9+10*KN,'(a5)')'#Bulk'
      IF(kr==2)write(n9+10*KN,'(a4)')'#HSh'
      IF(kr==3)write(n9+10*KN,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n9+10*KN,'(a15)')'#HSh-non-tetrah'    
      DO NKR =-NENRDELS,NENRDELS
!         write(*,*)NDPEO(kr,KN,NKR)
         IF(NDPEO(kr,KN,NKR).GT.0)THEN
            E_kjmol = ENRDEL*DFLOAT(NKR)
            DISTR_E=DFLOAT(NDPEO(kr,KN,NKR))/(DFLOAT(kv(kr))*ENRDEL)                   !normalize by the bin size
            WRITE(n9+10*KN,9)E_kjmol,DISTR_E
!check            WRITE(*,*)NKR,E_kjmol,DISTR_E
         ENDIF
      END DO
      WRITE(n9+10*KN,*)
!
      AVEPO(kr,KN) = AVEPO(kr,KN)/DFLOAT(kv(kr))
      WRITE(n00,*)
      WRITE(n00,*)
      WRITE(n00,*)'AVERAGE <EW_W',KN,'> kJ/mol = ',AVEPO(kr,KN)
      WRITE(n00,*)                    
   end do   

end do

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n4) 
close(n5) 
close(n6) 
close(n7) 
close(n8) 
close(n9) 

deallocate(NPDSG,NDPRO4,NDPRH2,NDPRHO1,NDPRO)
deallocate(AVWSG,AVCHHASG,AVRO4,AVRH2,AVRHO1)
deallocate(AVRO)
deallocate(NTDW4O)
deallocate(NONTDW)
deallocate(kv)

deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT(27X,F14.4,5X,F17.7)
    return

END SUBROUTINE wat_tetra



SUBROUTINE wat_tetra_sol(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                         nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,rwsoltetbulk,tetsampws,ZMBIO)          
                         
! Calculate water's tetrahedrality in different environments - include solute heavy atoms as possible tetrahedral vertices 

    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
    character(len=4),dimension(natms)                :: atomname    
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),intent(in)                          :: rwsoltetbulk   ! Distance relative to solute com for onset of bulk water
    integer,intent(in)                               :: tetsampws
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
   
! Local variables
! Tetrahedrality and distance distributions
    real(kind=8)                                     :: SG,CHHASG
    integer                                          :: nsh
    real(kind=8)                                     :: SGDEL
    integer,dimension(:,:),allocatable               :: NPDSG
    real(kind=8),dimension(:),allocatable            :: AVWSG,AVCHHASG
    integer,parameter                                :: nshdist = 50000                      !max. number of bins in tetrahedrality and distance dist.

    integer                                          :: NTDELS
    real(kind=8)                                     :: TETHQ
    real(kind=8)                                     :: TETHDIST
! 
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8)                                     :: dx,dy,dz,dr,dr2
    integer,dimension(:),allocatable                 :: kv
    integer,dimension(:,:),allocatable               :: kwlabel      !index of waters in each environment
    integer      :: n0, n1
    integer      :: i, j, k
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: iH, ihat, nt, ko, jt
    integer      :: io
    integer      :: L, KN, IT, NKR
    integer,parameter                                 :: INTMAX = 2147483647
    integer,dimension(:),allocatable                  :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable                  :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                                      :: dr_bulk,dr4nb
    real(kind=8)                                      :: RADMAX,RADMIN
    integer                                           :: NRSHMAX   
    integer                                           :: nradsh    
    integer                                           :: natcheck
!    integer                                           :: kor
!    character(len=1),dimension(natms)                 :: atom_symb
    integer                                           :: bulk_flag,kr_start
    integer                                           :: at_id       ! solute atomic index to analyse the respective solvation environment
    integer                                           :: kstps
    logical                                           :: filex
    
    integer                                           :: maxwat_bulk                !maximum number of waters in the bulk
    integer                                           :: nw_bulk                    !Bulk water counter

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Tetra_env_out/log_env_sol_tetra.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Tetra_env_out')

n0 = 100
n1 = 110

open(n0,file='Tetra_env_out/log_env_sol_tetra.dat',status='unknown',action='write')
open(n1,file='Tetra_env_out/wat_sol_tetra_dist.dat',status='unknown',action='write')    
    
    
dr_bulk = rwsoltetbulk                           !set distance from solute com to assume bulk water
bulk_flag = 1                                    !default - calculate tetrahedrality for bulk water
kr_start = 1                                     !environments loop started - if bulk_flag = 0; kr_start = 2 

maxwat_bulk = 500
if(nmolwat<maxwat_bulk)maxwat_bulk = nmolwat

if(dr_bulk==0.0)then
   bulk_flag = 0
   write(*,*)
   write(*,*)'Warning!!! No bulk region defined for water-solute tetrahedrality calculation'
   write(*,*)'Bulk tetrahedrality will not be calculated accordingly'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! No bulk region defined for water-solute tetrahedrality calculation'
   write(n0,*)'No bulk tetrahedrality will be calculated accordingly'
   write(n0,*)
else 
   write(*,*)
   write(*,*)'Warning!!! Maximum bulk waters for tetrahedrality calculation = ',maxwat_bulk
   write(*,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! Maximum bulk waters for tetrahedrality calculation = ',maxwat_bulk
   write(n0,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(n0,*)   
endif   

NRSHMAX = 4                                      !maximum number of environments
SGDEL   = 0.01D0                                 !tetrahedrality bin size 

 cmxsol = 0.0d0; cmysol=0.0d0; cmzsol=0.0d0

!Tetrahedrality
allocate(NPDSG(NRSHMAX,nshdist))
allocate(AVWSG(NRSHMAX),AVCHHASG(NRSHMAX))
NPDSG=0
AVWSG=0.0d0;AVCHHASG=0.0d0

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))   
allocate(kv(NRSHMAX))
allocate(kwlabel(NRSHMAX,nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled) 
NTDW4O = 0; NONTDW = 0
kv = 0
kstps = 0
kwlabel = 0

! Start tetrahedrality calculation

write(*,*)
write(*,*)'Starting water-solute tetrahedrality calculation...'
write(*,*)'Routine wat_tetra_sol'
write(*,*)'Results are printed to Tetra_env_out'
write(*,*)'water sampling frequency (steps) = ',tetsampws
write(*,*)'total number of time-steps =',nstep
write(*,*)

write(n0,*)
write(n0,*)'Starting water-solute tetrahedrality calculation...'
write(n0,*)'Routine wat_tetra_sol'
write(n0,*)'Results are printed to Tetra_env_out'
write(n0,*)'water sampling frequency (steps) = ',tetsampws
write(n0,*)'total number of time-steps =',nstep
write(n0,*)

do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
   NTDW4O = 0; NONTDW = 0
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!NG         if(j==1)atom_symb(i)=TRIM(atomname(i))
!check         if(j==1.and.atom_symb(i)/='H')write(*,*)atom_symb(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample for water molecules in specific environments
   if(j==1.or.mod(j,tetsampws)==0)then
      nw_bulk = 0
      write(*,*)'time-step',j,' sampling water molecules'
      kstps = kstps + 1                                                         !count number of origins sampled
      kwlabel = 0
      
!Find solute centre of mass for bulk water tetrahedrality         
      if(bulk_flag==1)call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water      
         dx = x(i) - cmxsol
         dy = y(i) - cmysol
         dz = z(i) - cmzsol
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
         if(bulk_flag==0)dr_bulk = 100.d0*cell(1)                                !define a distance larger than the box "radius"
!Version 13         if(dr>=dr_bulk)then                                                     !bulk water
         if(dr>=dr_bulk.and.nw_bulk<maxwat_bulk)then
            nw_bulk = nw_bulk + 1
            if(bulk_flag==0)write(*,*)'Warning!!! Found a water molecule in the bulk'    !bulk water should not be found if bulk_flag=0 
            nradsh = 1                                                           !hydration shell environment
            kv(nradsh) = kv(nradsh) + 1                                          !number of waters
            call sol_tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                 SG,CHHASG,atomname)
! distributions amd means
!            nsh = IDNINT(3.0D0/SGDEL+SG/SGDEL+0.5D0)
            nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0
            NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh)+1
            AVWSG(nradsh) = AVWSG(nradsh) + SG
            AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
            SG = 0.0D0
! End distributions and means                     
         endif
  
!Loop over atomic species whose solvation environment will be studied      

         do at_id=1,natHSh                                   !Loop over protein atomic species to analyze

!Hydration shell - tetrahedral and non tetrahedral water in the first HSh         
            k = natsol_id(at_id)
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then                   !water within hydration shell radius
               nradsh = 2                                                           !hydration shell environment              
!Avoid water molecules repetition
               if(kwlabel(nradsh,io)==0)then
                  kwlabel(nradsh,io)=1                                                       !Oxygen id for correction of the nb of waters
                  kv(nradsh) = kv(nradsh) + 1                                                !number of waters in the HSh (nradsh = 2) 
               else
                  goto 10                                                                    !sample another water molecule
               endif   
!End avoid water molecules repitition                             
               call sol_tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                    SG,CHHASG,atomname)
! distributions amd means
               nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0                                        
               NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh) + 1
               AVWSG(nradsh) = AVWSG(nradsh) + SG
               AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
               SG = 0.0D0
! End distributions and means          

!tetrahedral and non-tetrahedral waters in the hydration shell for sub-ensemble tetrahedral and non-tetrahedral distributions      
! FIND THE FOUR NEAREST O ATOMS J TO EACH OXYGEN I
               RADMAX  = 25.0D0
               RADMIN  =  0.0D0        
               JT      = 1
               DO WHILE(JT.LE.4) 
                  do kw = nmolsol*natmsol+1,natms,nwatsites                         !Loop over water oxygens
                     if(kw.ne.i)then
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOMS i AND OXYGEN ATOMS kw  
                        dx = x(i) - x(kw)
                        dy = y(i) - y(kw)
                        dz = z(i) - z(kw)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        IF((dr>RADMIN).AND.(dr<RADMAX))THEN
                           RADMAX = dr
!NTDW4O STORES THE O ATOMS NEIGHBORS OF i: JT = 1 FIRST NEIGHBOR; JT = 2 SECOND NEIGHBOR; JT = 3 THIRD NEIGHBOR; JT = 4 FOURTH NEIGHBOR
                           NTDW4O(JT) = kw
                        ENDIF
                     endif
                  end do
                  JT = JT + 1
                  RADMIN = RADMAX
                  RADMAX = 25.0D0
               END DO
!CHECK IF SOLUTE ATOMS ARE CLOSER THAN THE 4th NEAREST O ATOM          
               KI=NTDW4O(4)                                                     !Oxygen ID of the 4th nearest water neighbor 
               dx = x(i) - x(KI)
               dy = y(i) - y(KI)
               dz = z(i) - z(KI)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               dr4nb = dr                                                           !distance to the fourth nearest neighbor
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY SOLUTE ATOM IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = 1                                                                !run over all solute atoms - heavy and non-heavy atoms
               DO WHILE(L.LE.nmolsol*natmsol)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN             
                     NONTDW(io) = 1                                               !water i is non-tetrahedral
                     goto 5       
                  ENDIF
                  L = L+1
               END DO     
               5 continue           
               if(NONTDW(io)==1)then                                               !Non-tetrahedral water    
                  nradsh = 4                                                        !hydration shell (non-tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                                       !number of waters
                  call sol_tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                      SG,CHHASG,atomname)
! distributions amd means
                  nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0                                        
                  NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh) + 1
                  AVWSG(nradsh) = AVWSG(nradsh) + SG
                  AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
                  SG = 0.0D0
! End distributions and means          

               elseif(NONTDW(io)==0)then                                            !Tetrahedral water 
                  nradsh = 3                                                         !hydration shell (tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                                        !number of waters
                  call sol_tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                                       SG,CHHASG,atomname)
!check 0k            write(*,*)j,i,sg                 
! distributions amd means
                  nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0                                        
                  NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh) + 1
                  AVWSG(nradsh) = AVWSG(nradsh) + SG
                  AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
                  SG = 0.0D0
! End distributions and means          

               endif                                                                !end water non-tetrahedral/tetrahedral            
            endif                                                                   !end water is in the coordination sphere
         end do                                                                     !end solute atoms loop
         10 continue
      end do                                                                        !end loop over water oxygens
   endif                                                                         !end sampling of waters
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!Average number of waters on each environment

write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps

do kr = 1,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
!   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep) 
   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
!   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep)
   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps)
end do   

!Print distributions and mean values 

if(bulk_flag==0)kr_start=2   

write(n0,*)
write(n0,*)'TETRAHEDRALITY DISTRIBUTION'
write(n0,*)'==========================='
write(n0,*)
     
do kr = kr_start,NRSHMAX                                                      !loop over environments
   IF(kr==1)WRITE(n0,*) 
   IF(kr==1)WRITE(n0,*)'Bulk'
   IF(kr==1)WRITE(n0,*)'****'
   IF(kr==2)WRITE(n0,*)
   IF(kr==2)WRITE(n0,*)'Hydration Shell'
   IF(kr==2)WRITE(n0,*)'***************'
   IF(kr==3)WRITE(n0,*)
   IF(kr==3)WRITE(n0,*)'HSh-tetrahedral'
   IF(kr==3)WRITE(n0,*)'***************'
   IF(kr==4)WRITE(n0,*)
   IF(kr==4)WRITE(n0,*)'HSh-non-tetrahedral' 
   IF(kr==4)WRITE(n0,*)'*******************' 
   
! q tetrahedrality
   IF(kr==1)write(n1,'(a5)')'#Bulk'
   IF(kr==2)write(n1,'(a4)')'#HSh'
   IF(kr==3)write(n1,'(a11)')'#HSh-tetrah'
   IF(kr==4)write(n1,'(a15)')'#HSh-non-tetrah'   
! q VARIES BETWEEN -3 AND 1 FOR THE DEBENEDETTI NORMALIZATION IMPLEMENTED HERE
   NTDELS = 4.0D0/SGDEL
   DO IT = 1,NTDELS
      IF(NPDSG(kr,IT).GT.0)THEN
         TETHQ    = -3.01D0+SGDEL*DFLOAT(IT)
         TETHDIST = DFLOAT(NPDSG(kr,IT))/DFLOAT(kv(kr))                    !kv is the sum of the number of waters over the number of origins 
         WRITE(n1,9)TETHQ,TETHDIST
      ENDIF
   END DO
   WRITE(n1,*)
!
   AVWSG(kr)    = AVWSG(kr)/DFLOAT(kv(kr))
   AVCHHASG(kr) = AVCHHASG(kr)/DFLOAT(kv(kr))
   WRITE(n0,*)
   WRITE(n0,*)
   WRITE(n0,*)'Errington & Debenedetti perfect tetrahedron <q> = 1'
   WRITE(n0,*)'Average Errington & Debenedetti <q> = ',AVWSG(kr)
   WRITE(n0,*)
   WRITE(n0,*)'Chau & Hardwick perfect tetrahedron <q> = 0 '
   WRITE(n0,*)'Average Chau & Hardwick <Sg> =  <q> = ',AVCHHASG(kr)
   WRITE(n0,*)
   
end do

close(n0)
close(n1)

deallocate(NPDSG)
deallocate(AVWSG,AVCHHASG)
deallocate(NTDW4O)
deallocate(NONTDW)
deallocate(kv)

deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT(27X,F9.4,5X,F14.5)
  
    return

END SUBROUTINE wat_tetra_sol


!Start development


SUBROUTINE wat_LSI(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                   nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,LSI_sample,rwlsibulk,ZMBIO)               
! Calculate water's Local Structure Index (LSI) in different environments

    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
    character(len=4),dimension(natms)                :: atomname    
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment   
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),intent(in)                          :: rwlsibulk   ! Distance relative to solute com for onset of bulk water
    integer,intent(in)                               :: LSI_sample
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
   
! Local variables
! LSI distributions
    real(kind=8),dimension(:),allocatable        :: LSI
    integer,dimension(:,:),allocatable           :: NLSI 
    integer,dimension(:),allocatable             :: NLSI_count
    real(kind=8),dimension(:),allocatable        :: AVLSI,DISTR_LSI
    real(kind=8)                                 :: r_lsi
    real(kind=8)                                 :: LSI_del,LSI_ang2
    integer                                      :: nsh,LSI_NDELS,nkr
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms) 
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8)                                     :: dx,dy,dz,dr,dr2
    integer,dimension(:),allocatable                 :: kv
    integer,dimension(:,:),allocatable               :: kwlabel      !index of waters in each environment 
    integer      :: n0, n1
    integer      :: i, j, k
    integer      :: kw, iw, kt, kr, ki
    integer      :: nt, ko, jt
    integer      :: io
    integer      :: L, KN, IT
    integer,parameter                                 :: INTMAX = 2147483647
    integer,dimension(:),allocatable                  :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable                  :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                                      :: dr_bulk,dr4nb
    real(kind=8)                                      :: RADMAX,RADMIN
    integer                                           :: NRSHMAX   
    integer                                           :: nradsh    
    integer                                           :: natcheck
    integer                                           :: bulk_flag,kr_start
    integer                                           :: at_id                    ! solute atomic index to analyse the respective solvation environment
    integer                                           :: kstps
    logical                                           :: filex
    
    integer                                           :: maxwat_bulk                !maximum number of waters in the bulk
    integer                                           :: nw_bulk                    !Bulk water counter

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!LSI values
r_lsi  = 3.70d0                                 !original value of the LSI
LSI_del = 0.0005d0                              !bin for LSI distribution
LSI_NDELS= 1.0d0/LSI_del


!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='LSI_env_out/log_env_tetra.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p LSI_env_out')

n0 = 100
n1 = 110

open(n0,file='LSI_env_out/log_env_LSI.dat',status='unknown',action='write')
open(n1,file='LSI_env_out/wat_LSI_dist.dat',status='unknown',action='write')
    
    
dr_bulk = rwlsibulk                              !set distance from the solute COM to assume bulk water
bulk_flag = 1                                    !default - calculate tetrahedrality for bulk water
kr_start = 1                                     !environments loop started - if bulk_flag = 0; kr_start = 2 (skip environment bulk) - for printing purposes

maxwat_bulk = 500
if(nmolwat<maxwat_bulk)maxwat_bulk = nmolwat

if(dr_bulk==0.0)then
   bulk_flag = 0
   write(*,*)
   write(*,*)'Warning!!! No bulk region defined for LSI calculation'
   write(*,*)'Bulk LSI will not be calculated accordingly'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! No bulk region defined for LSI calculation'
   write(n0,*)'No bulk LSI will be calculated accordingly'
   write(n0,*)
else 
   write(*,*)
   write(*,*)'Warning!!! Maximum bulk waters for LSI calculation = ',maxwat_bulk
   write(*,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! Maximum bulk waters for LSI calculation = ',maxwat_bulk
   write(n0,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(n0,*)   
endif   

NRSHMAX = 4                                               !maximum number of environments

allocate(LSI(nmolwat))
allocate(NLSI(NRSHMAX,LSI_NDELS))
allocate(NLSI_count(NRSHMAX))
allocate(AVLSI(NRSHMAX))
allocate(DISTR_LSI(NRSHMAX))
LSI = 0
NLSI = 0
NLSI_count = 0
AVLSI = 0.0d0
DISTR_LSI = 0.0d0

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))   
allocate(kv(NRSHMAX))
allocate(kwlabel(NRSHMAX,nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled)    
NTDW4O = 0; NONTDW = 0
kv = 0
kstps = 0
kwlabel = 0

 cmxsol = 0.0d0; cmysol=0.0d0; cmzsol=0.0d0

! Start LSI calculation

write(*,*)
write(*,*)'Starting water LSI calculation...'
write(*,*)'Routine wat_LSI'
write(*,*)'Results are printed to LSI_env_out'
write(*,*)'water sampling frequency (steps) = ',LSI_sample
write(*,*)'total number of time-steps =',nstep
write(*,*)

write(n0,*)
write(n0,*)'Starting water LSI calculation...'
write(n0,*)'Routine wat_LSI'
write(n0,*)'Results are printed to LSI_env_out'
write(n0,*)'water sampling frequency (steps) = ',LSI_sample
write(n0,*)'total number of time-steps =',nstep
write(n0,*)

do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
   NTDW4O = 0; NONTDW = 0
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample for water molecules in specific environments


! sample frequency - calculate the LSI every LSI_sample time-step  

   if(j==1.or.mod(j,LSI_sample)==0)then
      nw_bulk = 0                                             !bulk water counter - re-initialized every time-step
      write(*,*)'time-step',j,' sampling water molecules'
      kstps = kstps + 1                                       !count number of origins sampled
      kwlabel = 0                                             !array with waters ID to avoid repetition
      
!Find solute centre of mass for bulk water tetrahedrality         
      if(bulk_flag==1)call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
         iw = i
!bulk water      
         dx = x(i) - cmxsol
         dy = y(i) - cmysol
         dz = z(i) - cmzsol
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
         if(bulk_flag==0)dr_bulk = 100.d0*cell(1)                                        !define a distance larger than the box "radius"
!Version 13         if(dr>=dr_bulk)then                                                             !bulk water
         if(dr>=dr_bulk.and.nw_bulk<maxwat_bulk)then
            nw_bulk = nw_bulk + 1
            if(bulk_flag==0)write(*,*)'Warning!!! Found a water molecule in the bulk'    !bulk water should not be found if bulk_flag=0                                                 
            nradsh = 1                                                                   !hydration shell environment
            kv(nradsh) = kv(nradsh) + 1                                                  !number of waters
!            call LSI_water(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
            call LSI_water_fast(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
!Compute LSI distribution and average value
            nsh = LSI(io)/LSI_del+0.5D0            
!check         if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check         write(*,*)LSI(io),nsh
            if(nsh==0)then
               write(*,*) 'Warning!!! LSI bin = 0; transformed to 1'
               write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
               nsh = 1
            endif   
!            nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!            write(*,*)nsh,nsh_
!            nsh = nsh_
            NLSI(nradsh,nsh) = NLSI(nradsh,nsh) + 1
            AVLSI(nradsh) = AVLSI(nradsh) + LSI(io)
            NLSI_count(nradsh) = NLSI_count(nradsh) + 1       !count number of values used for distribution and average value calculation           
!            write(*,*)nradsh,nsh,NLSI(nradsh,nsh) 
         endif                                                                   !end bulk

!Loop over atomic species whose solvation environment will be studied      

         do at_id=1,natHSh                                   !Loop over protein atomic species to analyze

!Hydration shell - tetrahedral and non tetrahedral water in the first HSh         
            k = natsol_id(at_id)
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)              
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then                            !water within hydration shell radius
               nradsh = 2                                                                    !hydration shell environment
!Avoid water molecules repetition
               if(kwlabel(nradsh,io)==0)then
                  kwlabel(nradsh,io)=1                                                       !Oxygen id for correction of the nb of waters
                  kv(nradsh) = kv(nradsh) + 1                                                !number of waters in the HSh (nradsh = 2) 
               else
!NG                  write(*,*)'water molecule already sampled ',io,' = ',i
                  goto 10                                                                    !sample another water molecule
               endif   
!End avoid water molecules repitition  
!            call LSI_water(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
               call LSI_water_fast(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
!Compute LSI distribution and average value
            nsh = LSI(io)/LSI_del+0.5D0            
!check         if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check         write(*,*)LSI(io),nsh
            if(nsh==0)then
               write(*,*) 'Warning!!! LSI bin = 0; transformed to 1'
               write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
               nsh = 1
            endif   
!            nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!            write(*,*)nsh,nsh_
!            nsh = nsh_
            NLSI(nradsh,nsh) = NLSI(nradsh,nsh) + 1
            AVLSI(nradsh) = AVLSI(nradsh) + LSI(io)
            NLSI_count(nradsh) = NLSI_count(nradsh) + 1       !count number of values used for distribution and average value calculation           

! End distributions and means          

!tetrahedral and non-tetrahedral waters in the hydration shell for sub-ensemble tetrahedral and non-tetrahedral distributions      
! FIND THE FOUR NEAREST O ATOMS J TO EACH OXYGEN I
               RADMAX  = 25.0D0
               RADMIN  =  0.0D0        
               JT      = 1
               DO WHILE(JT.LE.4) 
                  do kw = nmolsol*natmsol+1,natms-nions,nwatsites                         !Loop over water oxygens
                     if(kw.ne.i)then
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOMS i AND OXYGEN ATOMS kw  
                        dx = x(i) - x(kw)
                        dy = y(i) - y(kw)
                        dz = z(i) - z(kw)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        IF((dr>RADMIN).AND.(dr<RADMAX))THEN
                           RADMAX = dr
!NTDW4O STORES THE O ATOMS NEIGHBORS OF i: JT = 1 FIRST NEIGHBOR; JT = 2 SECOND NEIGHBOR; JT = 3 THIRD NEIGHBOR; JT = 4 FOURTH NEIGHBOR
                           NTDW4O(JT) = kw
                        ENDIF
                     endif
                  end do
                  JT = JT + 1
                  RADMIN = RADMAX
                  RADMAX = 25.0D0
               END DO
!CHECK IF SOLUTE ATOMS ARE CLOSER THAN THE 4th NEAREST O ATOM          
               KI=NTDW4O(4)                                                     !Oxygen ID of the 4th nearest water neighbor 
               dx = x(i) - x(KI)
               dy = y(i) - y(KI)
               dz = z(i) - z(KI)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               dr4nb = dr                                                           !distance to the fourth nearest neighbor
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY SOLUTE ATOM IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = 1                                                                !run over all solute atoms - heavy and non-heavy atoms
               DO WHILE(L.LE.nmolsol*natmsol)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN             
                     NONTDW(io) = 1                                               !water i is non-tetrahedral
                     goto 5       
                  ENDIF
                  L = L+1
               END DO     
!               5 continue
               
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY ION IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = nmolsol*natmsol + nmolwat*nwatsites + 1                                                                !run over all ions
               DO WHILE(L.LE.natms)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN   
                     write(*,*)'Water molecule next to ionic species'
                     NONTDW(io) = 1                                               !water i is non-tetrahedral
                     goto 5       
                  ENDIF
                  L = L+1
               END DO     
               5 continue  
               
!Non-tetrahedral waters (less than 4 water neighbors)               
               if(NONTDW(io)==1)then                                                            !Non-tetrahedral water    
                  nradsh = 4                                                                    !hydration shell (non-tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                                                   !number of waters
!                  call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
!                                   SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)
!                 call LSI_water(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
                  call LSI_water_fast(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
!Compute LSI distribution and average value
            nsh = LSI(io)/LSI_del+0.5D0            
!check         if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check         write(*,*)LSI(io),nsh
            if(nsh==0)then
               write(*,*) 'Warning!!! LSI bin = 0; transformed to 1'
               write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
               nsh = 1
            endif   
!            nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!            write(*,*)nsh,nsh_
!            nsh = nsh_
            NLSI(nradsh,nsh) = NLSI(nradsh,nsh) + 1
            AVLSI(nradsh) = AVLSI(nradsh) + LSI(io)
            NLSI_count(nradsh) = NLSI_count(nradsh) + 1       !count number of values used for distribution and average value calculation           

! End distributions and means          

               elseif(NONTDW(io)==0)then                                                        !Tetrahedral water 
                  nradsh = 3                                                                    !hydration shell (tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                                                   !number of waters
!                  call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
!                                   SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)
!                 call LSI_water(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
                  call LSI_water_fast(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
!Compute LSI distribution and average value
            nsh = LSI(io)/LSI_del+0.5D0            
!check         if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check         write(*,*)LSI(io),nsh
            if(nsh==0)then
               write(*,*) 'Warning!!! LSI bin = 0; transformed to 1'
               write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
               nsh = 1
            endif   
!            nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!            write(*,*)nsh,nsh_
!            nsh = nsh_
            NLSI(nradsh,nsh) = NLSI(nradsh,nsh) + 1
            AVLSI(nradsh) = AVLSI(nradsh) + LSI(io)
            NLSI_count(nradsh) = NLSI_count(nradsh) + 1   !count number of values used for distribution and average value calculation           

! End distributions and means          

               endif                                                                !end water non-tetrahedral/tetrahedral            
            endif                                                                   !end water is in the coordination sphere 
         end do                                                                     !end loop over atomic species
         10 continue
      end do                                                                        !end loop over water oxygens
   endif                                                                            !end sampling of waters
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj  


!Print LSI distribution

write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps

do kr = 1,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr
   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps)
end do   

if(bulk_flag==0)kr_start=2                                        !do not print the bulk

write(n0,*)
write(n0,*)'LSI DISTRIBUTIONS'
write(n0,*)'================='
write(n0,*)

do kr = kr_start,NRSHMAX                                          !loop over environments sampled
   IF(kr==1)WRITE(n0,*) 
   IF(kr==1)WRITE(n0,*)'Bulk'
   IF(kr==1)WRITE(n0,*)'****'
   IF(kr==2)WRITE(n0,*)
   IF(kr==2)WRITE(n0,*)'Hydration Shell'
   IF(kr==2)WRITE(n0,*)'***************'
   IF(kr==3)WRITE(n0,*)
   IF(kr==3)WRITE(n0,*)'HSh-tetrahedral'
   IF(kr==3)WRITE(n0,*)'***************'
   IF(kr==4)WRITE(n0,*)
   IF(kr==4)WRITE(n0,*)'HSh-non-tetrahedral' 
   IF(kr==4)WRITE(n0,*)'*******************' 
!
   IF(kr==1)write(n1,'(a5)')'#Bulk'
   IF(kr==2)write(n1,'(a4)')'#HSh'
   IF(kr==3)write(n1,'(a11)')'#HSh-tetrah'
   IF(kr==4)write(n1,'(a15)')'#HSh-non-tetrah'   
   do nkr = 1,LSI_NDELS
!check      write(*,*)LSI_NDELS,kr,nkr,NLSI(kr,nkr),NLSI_count(kr)  
      if(NLSI(kr,nkr) > 0)then
         LSI_ang2  = LSI_del*DFLOAT(nkr)
         DISTR_LSI(kr) = DFLOAT(NLSI(kr,nkr))/(LSI_del*DFLOAT(NLSI_count(kr)))                             !Normalize
         WRITE(n1,39)LSI_ang2,DISTR_LSI(kr)
      endif
   end do 
   WRITE(n1,*)
   AVLSI(kr) = AVLSI(kr)/DFLOAT(NLSI_count(kr))
   write(n0,*)
   write(n0,*)'LSI mean value (Ang) = ',AVLSI(kr)
   write(n0,*)
end do   



close(n0) 
close(n1) 


deallocate(NTDW4O)
deallocate(NONTDW)
deallocate(kv)

deallocate(LSI)
deallocate(NLSI)
deallocate(NLSI_count)
deallocate(AVLSI)
deallocate(DISTR_LSI)

deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT(27X,F9.4,5X,F14.5)
 39   FORMAT(27X,F9.4,5X,F14.5)
    return

END SUBROUTINE wat_LSI

!End developement


SUBROUTINE prot_contact(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,&
                        nsoltypes,intsf,rcont,atomtype,nb_res,chain_label_res,Natresid,name_res,numb_res,nres_cont_pair_t0,res_cont_pair_t0)  

! Calculate protein-protein contact map 
! Finds all residue-residue contacts during a trajectory based on a treshold distance "rcont"; sample freq. is "intsf" steps

    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname    
    integer,intent(in)                               :: nsoltypes               ! Number of solute distinct atomic species (C, H, O, N etc)
    integer,intent(in)                               :: intsf                   ! Protein contact map sampling frequency
    real(kind=8),intent(in)                          :: rcont                   ! distance for contact map - residues at a r < rcont are accounted
    character(len=4),dimension(natms),intent(in)     :: atomtype 
    integer,intent(in)                               :: nb_res                 ! Number of residues [GROMCAS]
    character(len=1),dimension(nb_res),intent(in)    :: chain_label_res        ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nb_res),intent(in)             :: Natresid               ! Number of atoms of each residue
    character(len=4),dimension(nb_res),intent(in)    :: name_res
    integer,dimension(nb_res),intent(in)             :: numb_res               ! Number of the residue: from 1 to the total number of residues of each chain

    integer,dimension(:,:),allocatable,intent(out)   :: res_cont_pair_t0
    integer,intent(out)                              :: nres_cont_pair_t0
    
!Local    
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                 !atomic index (solute and solvent index)
    integer      :: n0, n1, n3, n4
    integer      :: i, j, k
    integer      :: ic1,ic2,iats1,iats2,k1,k2,it,jt,jats1,jats2,kj1,kj2                                         !chain indexes for contact map
    real(kind=8)                                     :: dx,dy,dz,dr,dr2,rhi,rhf
    integer                                          :: natcheck
    logical                                          :: filex
    
    integer,dimension(:),allocatable                 :: res_cont_solv    
    character(len=1),dimension(natms)                :: atom_n          !first atom of atomname
    
    
!xtc2fortran_trj
!    character(len=15),intent(in)            :: trjfile    
    character(len=15)                       :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

allocate(res_cont_pair_t0(nb_res,nb_res))
allocate(res_cont_solv(nb_res))



!xtc2fortran
allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
!end xtc2fortran

inquire(file='Prot_cont_out/log_prot_contact.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Prot_cont_out')

n0 = 100
n1 = 110
n3 = 130
n4 = 140

open(n0,file='Prot_cont_out/log_prot_contact.dat',status='unknown',action='write')
open(n1,file='Prot_cont_out/prot_cont_map.dat',status='unknown',action='write')
open(n3,file='Prot_cont_out/List_RES_Solv_1.dat',status='unknown',action='write') 
open(n4,file='Prot_cont_out/List_RES_Solv_2.dat',status='unknown',action='write')

!Print headers

write(n1,*)'contact pair (FIRST) : RES_in_system; atom; RES_name; RES_in_chain; chain'
write(n1,*)'contact pair (SECOND): RES_in_system; atom; RES_name; RES_in_chain; chain'
write(n1,*)'total nb of residues = ',nb_res,' total nb of protein atoms = ',natmsol
write(n1,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(n1,*)

write(n3,*)'List of atom sites to be used as input for solvation analysis - sol_sites.dat'
write(n3,*)
write(n4,*)'List of atom sites to be used as input for solvation analysis - sol_sites.dat'
write(n4,*)

write(*,*)
write(*,*)'Starting protein contact map calculation...'
write(*,*)'Routine prot_contact'
write(*,*)'Results are printed to Prot_cont_out'
write(*,*)'traj sampling freq for contacts (steps) = ',intsf
write(*,*)'max. residue-residue contact dist. (Ang) = ',rcont
write(*,*)'total number of time-steps = ',nstep
write(*,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(*,*)
write(n0,*)
write(n0,*)'Starting protein contact map calculation...'
write(n0,*)'Routine prot_contact'
write(n0,*)'Results are printed to Prot_cont_out'
write(n0,*)'traj sampling freq for contact map (steps) = ',intsf
write(n0,*)'max. residue-residue contact dist. (Ang) = ',rcont
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(n0,*)

!Hydration shell onset and cut-off used to print out lists of residues (carbons-beta) to study solvation
rhi =  0.0
rhf = 10.0
             
kj1 = 0
kj2 = 0

res_cont_pair_t0 = 0
nres_cont_pair_t0 = 0


!first atom of atomname - used to choose all carbon atoms of HEM for solvation analysis
atom_n=atomname
!check write(*,*)atom_n


!xtc2fortran_trj
   
IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG      call f77_molfile_init

   call f77_molfile_open_read(handle(1),npart,infile,intype)

   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    
   
!end xtc2fortran_trj  

!Start trajectory analysis

do j = 1,nstep                                               !Loop over number of steps of trajectory
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
     read(ninput,*)natcheck
     if(natcheck.ne.natms)then
        write(*,*)'error - number of atoms from solv_inp is wrong'
        stop
     endif   
     do i = 1,natms
        read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!NG         if(j==1)atom_symb(i)=TRIM(atomname(i))
!check         if(j==1.and.atom_symb(i)/='H')write(*,*)atom_symb(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
      
!xtc2fortran_trj 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
   status = 1   ! status=1 on entry means read
   call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
   
   cell(1)=dble(box(1))
   cell(2)=dble(box(2))
   cell(3)=dble(box(3))
   kp=1
   do i = 1,npart*3,3
      x(kp)=xyz(i)
      y(kp)=xyz(i+1)
      z(kp)=xyz(i+2)
      kp = kp + 1
   end do   
   
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!Calculate protein contacts between chains ABCD and EFGH EVERY "intsf" steps
!The contacts will change along time and these will be accumulated in the array res_cont_pair_t0(ic1,ic2)...

   if(j==1.or.mod(j,intsf)==0)then                     !sample new contacts every intsf steps 
      write(n1,*)'Sample step = ',j
      write(n1,*)
      k1 = 0                               !Chains ABCD atomic index
      do ic1=1,nb_res/2
         do iats1 = 1,Natresid(ic1)
            k1 = k1 + 1
            k2 = natmsol/2                 !Chains EFGH atomic index
            do ic2=nb_res/2+1,nb_res
               do iats2 = 1,Natresid(ic2)
                  k2 = k2 + 1
                  dx = x(k1)-x(k2)
                  dy = y(k1)-y(k2)
                  dz = z(k1)-z(k2)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
!check                  write(*,*)dr
                  IF((dr.LE.rcont).and.(res_cont_pair_t0(ic1,ic2)==0))THEN             
                     res_cont_pair_t0(ic1,ic2)=1 
                     nres_cont_pair_t0 = nres_cont_pair_t0 + 1
                     write(n1,'(2(1x,I8),1x,a3,1x,i5,1x,a2,5x,2(1x,I8),1x,a3,1x,i5,1x,a2,3x,a7,i7)')&
                     ic1,k1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,k2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'step = ',j
                  ENDIF
               end do       !End number of atoms of residue ic2
            end do          !End residues ic2
         end do             !End number of atoms of residue ic1
      end do                !End residued ic1
      write(n1,*)
      write(n0,*)
      if(j==1)write(*,*)
      write(n0,'(a30,i7,1x,a1,1x,i9)')'nb of contact pairs after step',j,'=',nres_cont_pair_t0
      write(*,'(a30,i7,1x,a1,1x,i9)')'nb of contact pairs after step',j,'=',nres_cont_pair_t0
   
!check purposes only
!NG      do ic1 = 1,natmsol/2
!NG         do ic2 = natmsol/2+1,natmsol
!NG            dx = x(ic1)-x(ic2)
!NG            dy = y(ic1)-y(ic2)
!NG            dz = z(ic1)-z(ic2)
!NG            dx = dx - dnint(dx/cell(1))*cell(1)
!NG            dy = dy - dnint(dy/cell(2))*cell(2)
!NG            dz = dz - dnint(dz/cell(3))*cell(3)
!NG            dr2 = dx**2 + dy**2 + dz**2
!NG            dr  = dsqrt(dr2)
!NG            write(*,*)dr
!NG         end do
!NG      end do   
!end check         
      
   endif                                                                           !end sampling for res-res contacts 
          
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0)then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
    
end do                                                      !end time-step main loop
   
!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj      

write(*,*)
write(*,*)'finished sampling protein-protein contacts'
write(*,*)

!Solvation analysis input file - list of sites to analyse solvation based on the contact pairs 
!Prints a list of atomic indexes (carbons-beta) of contact residues that can then be used as input for solvation analysis
!For GLY the carbon alpha is printed and for HEM all carbons are printed
k1 = 0                               !Chains ABCD atomic index
do ic1=1,nb_res/2                    !residues of ABCD
   do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
      k1 = k1 + 1
      k2 = natmsol/2                 !Chains EFGH atomic index                  
      do ic2=nb_res/2+1,nb_res       !residues of EFGH
         do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
            k2 = k2 + 1
            if(res_cont_pair_t0(ic1,ic2)==1.and.res_cont_solv(ic1)==0)then
               res_cont_solv(ic1) = 1
               do jats1 = 1,Natresid(ic1)        !atoms of residue ic1
                  kj1 = k1 + jats1 - 1
                  if(atomname(kj1)=='CB')then
                     write(n3,'(i7,2(1x,f5.1),10x,a1,a3,1x,i5,1x,a2,1x,a3,1x,a4,i7)')&
                     kj1,rhi,rhf,'!',name_res(ic1),numb_res(ic1),chain_label_res(ic1),atomtype(kj1),atomname(kj1),ic1
                  elseif(name_res(ic1)=='GLY'.and.atomname(kj1)=='CA')then
                     write(n3,'(i7,2(1x,f5.1),10x,a1,a3,1x,i5,1x,a2,1x,a3,1x,a4,i7)')&
                     kj1,rhi,rhf,'!',name_res(ic1),numb_res(ic1),chain_label_res(ic1),atomtype(kj1),atomname(kj1),ic1
                  elseif(name_res(ic1)=='HEM'.and.atom_n(kj1)=='C')then
                     write(n3,'(i7,2(1x,f5.1),10x,a1,a3,1x,i5,1x,a2,1x,a3,1x,a4,i7)')&
                     kj1,rhi,rhf,'!',name_res(ic1),numb_res(ic1),chain_label_res(ic1),atomtype(kj1),atomname(kj1),ic1
                  endif              
               end do
            endif
            if(res_cont_pair_t0(ic1,ic2)==1.and.res_cont_solv(ic2)==0)then
               res_cont_solv(ic2) = 1
               do jats2 = 1,Natresid(ic2)        !atoms of residue ic2
                  kj2 = k2 + jats2 - 1
                  if(atomname(kj2)=='CB')then
                     write(n4,'(i7,2(1x,f5.1),10x,a1,a3,1x,i5,1x,a2,1x,a3,1x,a4,i7)')&
                     kj2,rhi,rhf,'!',name_res(ic2),numb_res(ic2),chain_label_res(ic2),atomtype(kj2),atomname(kj2),ic2
                  elseif(name_res(ic2)=='GLY'.and.atomname(kj2)=='CA')then
                     write(n4,'(i7,2(1x,f5.1),10x,a1,a3,1x,i5,1x,a2,1x,a3,1x,a4,i7)')&
                     kj2,rhi,rhf,'!',name_res(ic2),numb_res(ic2),chain_label_res(ic2),atomtype(kj2),atomname(kj2),ic2
                  elseif(name_res(ic2)=='HEM'.and.atom_n(kj2)=='C')then
                     write(n4,'(i7,2(1x,f5.1),10x,a1,a3,1x,i5,1x,a2,1x,a3,1x,a4,i7)')&
                     kj2,rhi,rhf,'!',name_res(ic2),numb_res(ic2),chain_label_res(ic2),atomtype(kj2),atomname(kj2),ic2
                  endif
               end do   
            endif  
            
         end do
      end do                   !end residues chain EFGH
   end do
end do                         !end residues chain ABDC


deallocate(xyz)
deallocate(x,y,z)
deallocate(res_cont_pair_t0)
deallocate(res_cont_solv)

return

END SUBROUTINE prot_contact


                   
SUBROUTINE potent_map(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,atomname_d,Zattype,mattype,ZMBIO,&
                      nsoltypes,kumbrella_trj,ncontsamp,kntrj,Coul_th_min,Coul_th_max,LJ126_th_min,LJ126_th_max,atomtype,chrg,nb_res,chain_label_res,&
                      Natresid,name_res,numb_res,kMIC,intsf_,rcont_,C6kjmol,C12kjmol,C6_IJ,C12_IJ)                      

!VERSION 12 - call a routine inside potent_map() to compute protein contacts instead of running the routine prot_contact() first...   
! Computes electrostatic and van der Waals interactions between the contact residues along time 
! The routine allows averaging over different trajectories, however, the contacts list is found from the first trajectory alone in the routine prot_contact
! This is convenient for averaging over equivalent trajectories e.g., SMD trajectories

!Version 15 - a single energetic criterion started beeing used to have the same resiudes in Coul, LJ, and the pot time profiles
!             the drawback is that the LJ interactions are not necessarily the strongest and in fact they will never be

    integer,intent(in)                                    :: ninput,nbox
    character(len=7),intent(in)                           :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                               :: dt                   ! time between frames (fs)
    integer,intent(in)                                    :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                               :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                                    :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites
    character(len=4),dimension(natms)                     :: atomname  
    character(len=4),dimension(natms),intent(in)          :: atomname_d  
    integer,intent(in)                                    :: nsoltypes               ! Number of solute distinct atomic species (C, H, O, N etc)
    integer,intent(in)                                    :: ncontsamp               ! Protein contact map electrostatics sampling frequency
    integer,intent(in)                                    :: kntrj                   ! Number of trajectories to average - SMD 
    integer,intent(in)                                    :: kumbrella_trj           ! 0: Trjs are not from umbrella sampling 1: Trjs are from umbrella sampling
!    real(kind=8),intent(in)                               :: Coul_th,LJ126_th        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!V15    real(kind=8),intent(in)                               :: Coul_th_min,LJ126_th_min        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!V15    real(kind=8),intent(in)                               :: Coul_th_max,LJ126_th_max        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
    real(kind=8),intent(in)                               :: Coul_th_min,Coul_th_max ! Potential energy threshold for printing res-res interactions [Min,Max]
    real(kind=8),intent(in)                               :: LJ126_th_min,LJ126_th_max
    character(len=4),dimension(natms),intent(in)          :: atomtype 
    real(kind=4),dimension(natms),intent(in)              :: chrg     
    integer,intent(in)                                    :: nb_res                 ! Number of residues [GROMCAS]
    character(len=1),dimension(nb_res),intent(in)         :: chain_label_res        ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nb_res),intent(in)                  :: Natresid               ! Number of atoms of each residue
    character(len=4),dimension(nb_res),intent(in)         :: name_res
    integer,dimension(nb_res),intent(in)                  :: numb_res               ! Number of the residue: from 1 to the total number of residues of each chain
    real(kind=8),dimension(natmsol),intent(in)            :: C6kjmol,C12kjmol       ! Solute Lennard-Jones parameters
!Version_12    integer,dimension(nb_res,nb_res),intent(in)           :: res_cont_pair_t0
!Version_12    integer,intent(in)                                    :: nres_cont_pair_t0
!Version_12
    integer,dimension(:,:),allocatable                    :: res_cont_pair_t0         !Res-Res identity that are below some input distance
    integer                                               :: nres_cont_pair_t0        
    integer,intent(in)                                    :: kMIC                     ! apply Minimum Image Convention or not to proteins
    integer,intent(in)                                    :: intsf_                   ! Protein contact map sampling frequency
    real(kind=8),intent(in)                               :: rcont_                   ! distance for contact map - residues at a r < rcont are accounted
!End !Version_12    
       
    real(kind=8),dimension(natmsol,natmsol),intent(in)    :: C6_IJ,C12_IJ           ! Solute Lennard-Jones C6 and C12 crossed parameters kj/mol nm**6; kj/mol nm**12
    character(len=2),dimension(nsoltypes),intent(in)      :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)          :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    real(kind=8),dimension(natmsol),intent(in)            :: ZMBIO
    
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                 !atomic index (solute and solvent index)
    integer      :: n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17
    integer      :: n19, n20, n21, n22, n25, n26, n27, n28, n29
    integer      :: i, j, k
    integer      :: ic1,ic2,iats1,iats2,k1,k2,it,jt,kt,jats1,jats2,kj1,kj2                                         !chain indexes for contact map
    integer      :: ktr,kl  
    integer      :: icm
    integer      :: ipC,ipthC,ipLJ,ipthLJ,ipp,ipthp
    integer      :: kstr                                              !correct a bug when sampling is every time-step
   
    real(kind=8)                                     :: dx,dy,dz,dr,dr2,rhi,rhf
    integer                                          :: natcheck
    logical                                          :: filex
    
    integer                                          :: ndist
    real(kind=8),dimension(:,:,:),allocatable        :: PCoul_pair
    real(kind=8),dimension(:,:),allocatable          :: PCoul_pair_mean
    real(kind=8),dimension(:,:),allocatable          :: PCoul_pair_stdev
    
    real(kind=8),dimension(:,:,:),allocatable        :: PLJ_pair
    real(kind=8),dimension(:,:),allocatable          :: PLJ_pair_mean
    real(kind=8),dimension(:,:),allocatable          :: PLJ_pair_stdev
    
    real(kind=8),dimension(:,:,:),allocatable        :: pot_pair
    real(kind=8),dimension(:,:),allocatable          :: pot_pair_mean
    real(kind=8),dimension(:,:),allocatable          :: pot_pair_stdev
!Version 32
    real(kind=8),dimension(:,:,:),allocatable        :: rdist_pair
    real(kind=8),dimension(:,:),allocatable          :: rdist_pair_mean
    real(kind=8),dimension(:,:),allocatable          :: rdist_pair_stdev
    
    
    real(kind=8)                                     :: C6,C12
    
    real(kind=8)                                     :: TSTEP,PSECS,NSECS
    character(len=2)                                 :: klab            ! trj file number label 01, 02, etc, 99
    character(len=1),dimension(natms)                :: atom_n          !first atom of atomname
    
!Electrostatics   
    real(kind=8),parameter                              :: AVSNO = 6.02214199D+23
! 1/(4*PI*EPS0) 8.98750D9 Nm2C-2 = JAngC-2 8.9875D19
    real(kind=8),parameter                              :: PMTTV=8.987551760D19
! electron (C) 
    real(kind=8),parameter                              :: ELCH=1.602177D-019
    real(kind=8)                                        :: enconv_kjmol
    real(kind=8)                                        :: Q1,Q2
    real(kind=8)                                        :: coul_pot,coul_pot_t0
    real(kind=8)                                        :: lj_pot,lj_pot_t0
    real(kind=8)                                        :: pe_pot
!Version 15 - turned pe_pot_t0 an array     
    real(kind=8),dimension(:,:),allocatable             :: pe_pot_t0
    
    real(kind=8),dimension(:),allocatable               :: sum_Coul_en_t       !sum of the Coul en for all contacts with pe outside threshold window
    real(kind=8),dimension(:),allocatable               :: sum_LJ_en_t         !sum of the LJ en for all contacts with pe outside threshold window
    real(kind=8),dimension(:),allocatable               :: sum_pot_en_t        !sum of the pe for all contacts with pe outside threshold window
    real(kind=8)                                        :: coul_pot_t,lj_pot_t,pe_pot_t
    
    real(kind=8),dimension(:),allocatable               :: sum_Coul_en_all       !sum of the Coul en for all contacts 
    real(kind=8),dimension(:),allocatable               :: sum_LJ_en_all         !sum of the LJ en for all contacts 
    real(kind=8),dimension(:),allocatable               :: sum_pot_en_all        !sum of the pe for all contacts 
    
    real(kind=8)                                        :: sum_Coul_en_t_mean,sum_LJ_en_t_mean,sum_pot_en_t_mean
    real(kind=8)                                        :: sum_Coul_en_all_mean,sum_LJ_en_all_mean,sum_pot_en_all_mean
    
    integer                                             :: maxres              !maximum number of contacts allowed for printing full time-energy profile
    integer                                             :: ncont_th,ncont_Coul_th,ncont_LJ_th    !number of contacts with energy outside threshold window  
    integer                                             :: ncont_th_,ncont_Coul_th_,ncont_LJ_th_
    real(kind=8)                                        :: cmxsol_I,cmysol_I,cmzsol_I         !solute centre of mass
    real(kind=8)                                        :: cmxsol_II,cmysol_II,cmzsol_II      !solute centre of mass
    real(kind=8)                                        :: dcmx,dcmy,dcmz,dcmr,dcmr2
    real(kind=8),dimension(:),allocatable               :: dcm_trj
    integer                                             :: nstart,nend
!stdeviations mean energies
    real(kind=8)                                        :: sum_Coul_en_thr_stdev 
    real(kind=8)                                        :: sum_LJ_en_thr_stdev   
    real(kind=8)                                        :: sum_pot_en_thr_stdev     
    real(kind=8)                                        :: sum_Coul_en_all_stdev 
    real(kind=8)                                        :: sum_LJ_en_all_stdev   
    real(kind=8)                                        :: sum_pot_en_all_stdev
!
!SMD variables
    real(kind=8),dimension(:,:),allocatable             :: dcm_trj_SMD
    real(kind=8),dimension(:),allocatable               :: stdev_dcm_SMD    
    real(kind=8),dimension(:,:,:,:),allocatable         :: PCoul_SMD
    real(kind=8),dimension(:,:,:),allocatable           :: stdev_PCoul_SMD
    real(kind=8),dimension(:,:,:,:),allocatable         :: PLJ_SMD
    real(kind=8),dimension(:,:,:),allocatable           :: stdev_PLJ_SMD
    real(kind=8),dimension(:,:,:,:),allocatable         :: pot_SMD
    real(kind=8),dimension(:,:,:),allocatable           :: stdev_pot_SMD
    
    real(kind=8),dimension(:,:),allocatable             :: sum_Coul_en_all_SMD
    real(kind=8),dimension(:,:),allocatable             :: sum_LJ_en_all_SMD
    real(kind=8),dimension(:,:),allocatable             :: sum_pot_en_all_SMD
    real(kind=8),dimension(:,:),allocatable             :: sum_Coul_en_thr_SMD
    real(kind=8),dimension(:,:),allocatable             :: sum_LJ_en_thr_SMD
    real(kind=8),dimension(:,:),allocatable             :: sum_pot_en_thr_SMD
    
    real(kind=8),dimension(:),allocatable               :: stdev_Coul_all_SMD
    real(kind=8),dimension(:),allocatable               :: stdev_LJ_all_SMD  
    real(kind=8),dimension(:),allocatable               :: stdev_pot_all_SMD 
    real(kind=8),dimension(:),allocatable               :: stdev_Coul_t_SMD  
    real(kind=8),dimension(:),allocatable               :: stdev_LJ_t_SMD    
    real(kind=8),dimension(:),allocatable               :: stdev_pot_t_SMD   
    
    real(kind=8),dimension(:),allocatable               :: sum_check_Coul
!END SMD variables
    
    real(kind=8),dimension(:,:,:,:),allocatable         :: pot_thr_ftrj
    real(kind=8),dimension(:,:),allocatable             :: pot_thr_ftrj_mean
    real(kind=8),dimension(:,:,:,:),allocatable         :: PLJ_thr_ftrj
    real(kind=8),dimension(:,:),allocatable             :: PLJ_thr_ftrj_mean
    
    
!xtc2fortran_trj
!    character(len=15),intent(in)            :: trjfile    
    character(len=15)                       :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp, nrc, nrc_C, nrc_LJ, nrc_p
    
    integer                                 :: nrc_C_, nrc_LJ_, nrc_p_
!end xtc2fortran_trj

kstr = 1
!Version 18   if(ncontsamp==1)kstr=0

!Version 12
allocate(res_cont_pair_t0(nb_res,nb_res))

allocate(dcm_trj(nstep/ncontsamp+kstr))
!Version 15
allocate(pe_pot_t0(nb_res,nb_res))
allocate(PCoul_pair(nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(PCoul_pair_mean(nb_res,nb_res),PCoul_pair_stdev(nb_res,nb_res))
allocate(PLJ_pair(nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(PLJ_pair_mean(nb_res,nb_res),PLJ_pair_stdev(nb_res,nb_res))
allocate(pot_pair(nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(pot_pair_mean(nb_res,nb_res),pot_pair_stdev(nb_res,nb_res))
!Version 32
allocate(rdist_pair(nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(rdist_pair_mean(nb_res,nb_res),rdist_pair_stdev(nb_res,nb_res))
!end v32
allocate(sum_Coul_en_t(nstep/ncontsamp+kstr),sum_LJ_en_t(nstep/ncontsamp+kstr),sum_pot_en_t(nstep/ncontsamp+kstr))
allocate(sum_Coul_en_all(nstep/ncontsamp+kstr),sum_LJ_en_all(nstep/ncontsamp+kstr),sum_pot_en_all(nstep/ncontsamp+kstr))

!****************************************************************************************************************SMD
!Multiple trajectories (SMD) - calculate average and stdev over trajectories
if(kntrj > 1)then
   allocate(dcm_trj_SMD(kntrj,nstep/ncontsamp+kstr))
   allocate(stdev_dcm_SMD(nstep/ncontsamp+kstr))   
   allocate(PCoul_SMD(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
   allocate(stdev_PCoul_SMD(nb_res,nb_res,nstep/ncontsamp+kstr))
   allocate(PLJ_SMD(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
   allocate(stdev_PLJ_SMD(nb_res,nb_res,nstep/ncontsamp+kstr))
   allocate(pot_SMD(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
   allocate(stdev_pot_SMD(nb_res,nb_res,nstep/ncontsamp+kstr))
!instantaneous values for the full and threshold sums for each trj
   allocate(sum_Coul_en_all_SMD(kntrj,nstep/ncontsamp+kstr)) 
   allocate(sum_LJ_en_all_SMD(kntrj,nstep/ncontsamp+kstr))    
   allocate(sum_pot_en_all_SMD(kntrj,nstep/ncontsamp+kstr))   
   allocate(sum_Coul_en_thr_SMD(kntrj,nstep/ncontsamp+kstr))
   allocate(sum_LJ_en_thr_SMD(kntrj,nstep/ncontsamp+kstr))  
   allocate(sum_pot_en_thr_SMD(kntrj,nstep/ncontsamp+kstr)) 
   
   allocate(stdev_Coul_all_SMD(nstep/ncontsamp+kstr))
   allocate(stdev_LJ_all_SMD  (nstep/ncontsamp+kstr))
   allocate(stdev_pot_all_SMD (nstep/ncontsamp+kstr))  
   allocate(stdev_Coul_t_SMD  (nstep/ncontsamp+kstr))   
   allocate(stdev_LJ_t_SMD    (nstep/ncontsamp+kstr))     
   allocate(stdev_pot_t_SMD   (nstep/ncontsamp+kstr))  
   
   allocate(sum_check_Coul(nstep/ncontsamp+kstr))
endif
!*************************************************************************************************************end SMD

allocate(pot_thr_ftrj(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(pot_thr_ftrj_mean(nb_res,nb_res)) 

!Version 20
allocate(PLJ_thr_ftrj(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(PLJ_thr_ftrj_mean(nb_res,nb_res)) 
! End V20

!xtc2fortran
allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
!end xtc2fortran

!call SYSTEM('mkdir -p Prot_potent_out')
inquire(file='Prot_potent_out/log_prot_potent.dat',exist=filex)
if(.not.filex) call SYSTEM('mkdir -p Prot_potent_out')

!call SYSTEM('mkdir -p Prot_potent_out')
call SYSTEM('mkdir -p Prot_potent_out/Coul_e')
call SYSTEM('mkdir -p Prot_potent_out/vdW_e')
call SYSTEM('mkdir -p Prot_potent_out/pot_e')
call SYSTEM('mkdir -p Prot_potent_out/com_d')
call SYSTEM('mkdir -p Prot_potent_out/VDW_threshold')

call SYSTEM('mkdir -p Prot_potent_out/SMD_trj_averages')

!log
n0 = 100
!com-com
n1 = 110
!en
n2 = 120
n3 = 130
n4 = 140
n5 = 150
n6 = 160
n7 = 170
n8 = 180
n9 = 190
n10 = 200
n11 = 210
n12 = 220
n13 = 230
n15 = 250
n16 = 260
n17 = 270

!time-com-com
n19 = 290

!SMD - com-com en
n20 = 300
n21 = 310
n22 = 320

n25 = 350
n26 = 360
n27 = 370

!Version 32
n28 = 380
n29 = 390
!end v32

open(n0,file='Prot_potent_out/log_prot_potent.dat',status='unknown',action='write')

!centre of mass
open(n1,file='Prot_potent_out/com_d/prot_com_first_trj.dat',status='unknown',action='write')

!Full res-res interactions - electrostatics, vdW, pot - time (ns) PROFILES
open(n2,file='Prot_potent_out/Coul_e/cont_map_Coul_inter.dat',status='unknown',action='write')
open(n3,file='Prot_potent_out/vdW_e/cont_map_LJ_inter.dat',status='unknown',action='write')
open(n4,file='Prot_potent_out/pot_e/cont_map_POTEn_inter.dat',status='unknown',action='write')
!Threshold res-res interactions - electrostatics, vdW, pot - time (ns) PROFILES en cut-off = potential energy
open(n5,file='Prot_potent_out/Coul_e/cont_map_Coul_thrsh.dat',status='unknown',action='write')
open(n6,file='Prot_potent_out/vdW_e/cont_map_LJ_thrsh.dat',status='unknown',action='write')
open(n7,file='Prot_potent_out/pot_e/cont_map_POTEn_thrsh.dat',status='unknown',action='write')
!Full time average res-res interactions - electrostatics, vdW, pot [TIME STANDARD DEVIATIONS]
open(n8,file='Prot_potent_out/Coul_e/cont_map_Coul_time_mean.dat',status='unknown',action='write')
open(n9,file='Prot_potent_out/vdW_e/cont_map_LJ_time_mean.dat',status='unknown',action='write')
open(n10,file='Prot_potent_out/pot_e/cont_map_POTEn_time_mean.dat',status='unknown',action='write')
!Threshold time averageres-res interactions - electrostatics, vdW, pot [TIME STANDARD DEVIATIONS]
open(n11,file='Prot_potent_out/Coul_e/cont_map_Coul_time_mean.agr',status='unknown',action='write')
open(n12,file='Prot_potent_out/vdW_e/cont_map_LJ_time_mean.agr',status='unknown',action='write')
open(n13,file='Prot_potent_out/pot_e/cont_map_POTEn_time_mean.agr',status='unknown',action='write')
!Threshold time averageres-res interactions - electrostatics, vdW, pot [TIME STANDARD DEVIATIONS] - INCLUDES THE LABELS of EACH res-res interaction
open(n15,file='Prot_potent_out/Coul_e/cont_map_Coul_time_mean_thrsh.dat',status='unknown',action='write')
open(n16,file='Prot_potent_out/vdW_e/cont_map_LJ_time_mean_thrsh.dat',status='unknown',action='write')
open(n17,file='Prot_potent_out/pot_e/cont_map_POTEn_time_mean_thrsh.dat',status='unknown',action='write')

!Version 20 - vdW thrshold energies
!Threshold res-res interactions - electrostatics, vdW, pot - time (ns) PROFILES - en cut-off = vdW
open(n25,file='Prot_potent_out/VDW_threshold/cont_map_Coul_thrsh_vdW.dat',status='unknown',action='write')
open(n26,file='Prot_potent_out/VDW_threshold/cont_map_LJ_thrsh_vdW.dat',status='unknown',action='write')
open(n27,file='Prot_potent_out/VDW_threshold/cont_map_POTEn_thrsh_vdW.dat',status='unknown',action='write')
!End v20

!Averages and stdev over multiple SMD trajectories - time versus mean values +/- stdev

!time versus com
open(n19,file='Prot_potent_out/SMD_trj_averages/SMD_trjs_mean_com.dat',status='unknown',action='write')

!com distance versus energy [TRJ STANDARD DEVIATIONS stdev(com-com) and stdev(en)] 
open(n20,file='Prot_potent_out/SMD_trj_averages/SMD_trjs_COM_EN_COUL.dat',status='unknown',action='write')
open(n21,file='Prot_potent_out/SMD_trj_averages/SMD_trjs_COM_EN_LJ.dat',status='unknown',action='write')
open(n22,file='Prot_potent_out/SMD_trj_averages/SMD_trjs_COM_EN_POT.dat',status='unknown',action='write')

!Version 32
open(n28,file='Prot_potent_out/pot_e/cont_map_rdist_time_mean.agr',status='unknown',action='write')
open(n29,file='Prot_potent_out/pot_e/contour_plot_inter_dist_POTen.dat',status='unknown',action='write')
!End v32           
           
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3                      !

enconv_kjmol = 1.0d-3*AVSNO                   !conversion factor J ---> kJ/mol

dcm_trj = 0.0d0

PCoul_pair = 0.0d0
PCoul_pair_mean = 0.0d0 
PCoul_pair_stdev = 0.0d0

PLJ_pair = 0.0d0
PLJ_pair_mean = 0.0d0 
PLJ_pair_stdev = 0.0d0

pot_pair = 0.0d0
pot_pair_mean = 0.0d0 
pot_pair_stdev = 0.0d0

!Version 32
rdist_pair = 0.0d0
rdist_pair_mean = 0.0d0 
rdist_pair_stdev = 0.0d0
!End v32

sum_Coul_en_t = 0.0d0
sum_LJ_en_t = 0.0d0
sum_pot_en_t = 0.0d0

sum_Coul_en_all = 0.0d0
sum_LJ_en_all = 0.0d0
sum_pot_en_all = 0.0d0

sum_Coul_en_t_mean = 0.0d0 
sum_LJ_en_t_mean = 0.0d0 
sum_pot_en_t_mean = 0.0d0

sum_Coul_en_all_mean = 0.0d0 
sum_LJ_en_all_mean = 0.0d0 
sum_pot_en_all_mean = 0.0d0

sum_Coul_en_thr_stdev = 0.0d0
sum_LJ_en_thr_stdev   = 0.0d0 
sum_pot_en_thr_stdev  = 0.0d0

sum_Coul_en_all_stdev = 0.0d0
sum_LJ_en_all_stdev   = 0.0d0
sum_pot_en_all_stdev  = 0.0d0

kj1 = 0
kj2 = 0

maxres = 5000           !Maximum number of residue contacts allowed for printing energy-time profiles

ncont_th = 0            !Number of contacts with potential energy outside threshold
ncont_Coul_th = 0       !Number of contacts with Coulomb energy outside threshold
ncont_LJ_th = 0         !Number of contacts with vdW energy outside threshold
ncont_th_ = 0            !Number of contacts with potential energy outside threshold
ncont_Coul_th_ = 0       !Number of contacts with Coulomb energy outside threshold
ncont_LJ_th_ = 0         !Number of contacts with vdW energy outside threshold

nrc = 0

nrc_C = 0
nrc_LJ = 0
nrc_p = 0
nrc_C_ = 0
nrc_LJ_ = 0
nrc_p_ = 0

 cmxsol_I = 0.0d0;  cmysol_I = 0.0d0;  cmzsol_I = 0.0d0
 cmxsol_II = 0.0d0; cmysol_II = 0.0d0; cmzsol_II = 0.0d0
 
ndist=0 

!****************************************************************************************************************SMD
if(kntrj > 1)then 
    PCoul_SMD = 0.0d0
    PLJ_SMD = 0.0d0
    pot_SMD = 0.0d0
   
   stdev_dcm_SMD   = 0.0d0
   stdev_PCoul_SMD = 0.0d0
   stdev_PLJ_SMD   = 0.0d0
   stdev_pot_SMD   = 0.0d0
   
   sum_Coul_en_all_SMD = 0.0d0
   sum_LJ_en_all_SMD   = 0.0d0
   sum_pot_en_all_SMD  = 0.0d0
   sum_Coul_en_thr_SMD = 0.0d0
   sum_LJ_en_thr_SMD   = 0.0d0
   sum_pot_en_thr_SMD  = 0.0d0
   
   stdev_Coul_all_SMD = 0.0d0
   stdev_LJ_all_SMD   = 0.0d0
   stdev_pot_all_SMD  = 0.0d0
   stdev_Coul_t_SMD   = 0.0d0
   stdev_LJ_t_SMD     = 0.0d0
   stdev_pot_t_SMD    = 0.0d0  
   
   sum_check_Coul = 0.0d0
endif   
!****************************************************************************************************************end SMD   
   
pot_thr_ftrj = 0.0d0
pot_thr_ftrj_mean = 0.0d0   

!Version 20
PLJ_thr_ftrj = 0.00D0
PLJ_thr_ftrj_mean = 0.00D0
!End V20
   
!Version_12
res_cont_pair_t0 = 0
nres_cont_pair_t0 = 0


!VERSION_12
!
write(*,*)
write(*,'(10X,A30,10X)')'*****************************'
write(*,'(10X,A30,10X)')'CALLING ROUTINE prot_res_cont'
write(*,'(10X,A30,10X)')'*****************************'
write(*,*)
write(*,'(a72)')'Calculating protein contacts along the first trajectory [md_trj.xtc]...'
write(*,'(a131)')'These will be used in every traj [md_trj01.xtc = md_trj.xtc; md_trj02.xtc; md_trj03.xtc; ... ] to compute mean values and stdev...' 
call f77_molfile_init     ![called in main before potent_map call]
call prot_res_cont(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,&
                   nsoltypes,kumbrella_trj,kMIC,intsf_,rcont_,atomtype,nb_res,chain_label_res,Natresid,name_res,numb_res,nres_cont_pair_t0,res_cont_pair_t0) 
!
!END VERSION_12


!Start potential map calculation

write(*,*)
write(*,*)'************************************************************'
write(*,*)'*                        rEaD mE!!!                        *'
write(*,*)'* . For number of trajectories > 1 these must be named as: *'
write(*,*)'* md_trj01.xtc, md_trj02.xtc,...,md_trjmn.xtc              *'
write(*,*)'* with mn = kntrj = number of trajectories                 *'
write(*,*)'* . Trajectories must have the exact same time-lenght      *'
write(*,*)'* . The contact map used is found from md_trj.xtc          *'
write(*,*)'* where md_trj.xtc may be one of the md_trjmn.xtc trajs    *'
write(*,*)'* IMPORTANT!!! this must be named md_trj.xtc               *'
write(*,*)'*                                                          *'
write(*,*)'* . Trajectories from different SMD trajectories may be    *'
write(*,*)'* used - contacts found from md_trj.xtc will be used for   *'
write(*,*)'* all trajectories                                         *'
write(*,*)'* HeRe AvErAgEs OvEr TrJs WiLl Be CaLcUlAtEd               *'
write(*,*)'*                                                          *'
write(*,*)'* ENERGY THRESHOLD FILES:                                  *'
write(*,*)'*                                                          *'
write(*,*)'* residue pairs with intermolecular potential en outside   *'
write(*,*)'* input threshold (kJ/mol) ',Coul_th_min,'-',Coul_th_max,' *'
write(*,*)'*                                                          *'
write(*,*)'* cont_map_Coul_thrsh.dat:            Coulomb time profile *'
write(*,*)'*                                                          *'
write(*,*)'* cont_map_LJ_thrsh.dat:        van der Waals time profile *'
write(*,*)'*                                                          *'
write(*,*)'* cont_map_POTEn_thrsh.dat:       intermol pe time profile *'
write(*,*)'*                                                          *'
write(*,*)'************************************************************'

write(*,*)
write(*,*)'Starting protein potential map calculation...'
write(*,*)'Routine potent_map'
write(*,*)'Results are printed to Prot_potent_out'
write(*,*)'Traj sampling freq for Coul and vdW potential (steps) = ',ncontsamp
write(*,*)'number of trajectories = ',kntrj
write(*,*)'total number of time-steps = ',nstep
write(*,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(*,*)'Number of contact pairs to be analysed = ',nres_cont_pair_t0
write(*,'(a48,1x,f12.4)')'Fraction of contact pairs to be analysed (%) = ',100.0*dfloat(nres_cont_pair_t0)/dfloat((nb_res/2)**2)
if(nres_cont_pair_t0 > maxres)then
   write(*,*)
   write(*,*)'WARNING!!! - too many contacts'
   write(*,*)'protein contacts time-energy profiles will not be printed'
   write(*,*)'IF YOU STILL NEED TO PRINT INCREASE maxres = ',maxres 
endif   
write(*,*)
write(n0,*)
write(n0,*)'Starting protein potential map calculation...'
write(n0,*)'Routine potent_map'
write(n0,*)'Results are printed to Prot_potent_out'
write(n0,*)'Traj sampling freq for Coul and vdW potential (steps) = ',ncontsamp
write(n0,*)'number of trajectories = ',kntrj
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(n0,*)'Number of contact pairs to be analysed = ',nres_cont_pair_t0
write(n0,'(a48,1x,f12.4)')'Fraction of contact pairs to be analysed (%) = ',100.0*dfloat(nres_cont_pair_t0)/dfloat((nb_res/2)**2)
if(nres_cont_pair_t0 > maxres)then
   write(n0,*)
   write(n0,*)'WARNING!!! - too many contacts '
   write(n0,*)'protein contacts time-energy profiles will not be printed'
   write(n0,*)'IF YOU STILL NEED TO PRINT INCREASE maxres = ',maxres 
endif   
write(n0,*)

write(n1,'(a1,3x a7,3x,a12)')'#','time/ns','com dist/Ang'

write(n2,'(a62)')'# time (ns) versus AVERAGE [over trajectories] Coul en(kJ/mol)'
write(n3,'(a62)')'# time (ns) versus AVERAGE [over trajectories] LJ   en(kJ/mol)'
write(n4,'(a62)')'# time (ns) versus AVERAGE [over trajectories] pot  en(kJ/mol)'

write(n5,'(a62)')'# time (ns) versus AVERAGE [over trajectories] Coul en(kJ/mol)'
write(n6,'(a62)')'# time (ns) versus AVERAGE [over trajectories] LJ   en(kJ/mol)'
write(n7,'(a62)')'# time (ns) versus AVERAGE [over trajectories] pot  en(kJ/mol)'

!Version 20
write(n25,'(a62)')'# time (ns) versus AVERAGE [over trajectories] Coul en(kJ/mol)'
write(n26,'(a62)')'# time (ns) versus AVERAGE [over trajectories] LJ   en(kJ/mol)'
write(n27,'(a62)')'# time (ns) versus AVERAGE [over trajectories] pot  en(kJ/mol)'
!End v20

write(n8,'(a76)')'# STANDARD DEVIATION RELATIVE TO TIME AVERAGE -  NOT NUMBER OF TRAJECTORIES'
write(n8,'(a66)')'# EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES'

write(n9,'(a76)')'# STANDARD DEVIATION RELATIVE TO TIME AVERAGE -  NOT NUMBER OF TRAJECTORIES'
write(n9,'(a66)')'# EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES'

write(n10,'(a76)')'# STANDARD DEVIATION RELATIVE TO TIME AVERAGE -  NOT NUMBER OF TRAJECTORIES'
write(n10,'(a66)')'# EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES'

write(n15,'(a76)')'# STANDARD DEVIATION RELATIVE TO TIME AVERAGE -  NOT NUMBER OF TRAJECTORIES'
write(n15,'(a66)')'# EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES'

write(n16,'(a76)')'# STANDARD DEVIATION RELATIVE TO TIME AVERAGE -  NOT NUMBER OF TRAJECTORIES'
write(n16,'(a66)')'# EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES'

write(n17,'(a76)')'# STANDARD DEVIATION RELATIVE TO TIME AVERAGE -  NOT NUMBER OF TRAJECTORIES'
write(n17,'(a66)')'# EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES'

!SMD headers
!time - d<com-com>
write(n19,'(a24,1x,i6,1x,a12)')'# AVERAGE and STDEV over',kntrj,'trajectories'
write(n19,'(a1,3x a7,3x,a14,3x,a9)')'#','time/ns','<com dist>/Ang','stdev/Ang'
write(n19,*)

!d<com-com> - <en>
write(n20,'(a24,1x,i6,1x,a12)')'# AVERAGE and STDEV over',kntrj,'trajectories'
write(n20,'(a1,3x a7,3x,a17,3x,a18,3x,a18)')'#','com/Ang','<en Coul>/kJ/mol','stdev_com/Ang','stdev_Coul/kJ/mol'
write(n20,*)

write(n21,'(a24,1x,i6,1x,a12)')'# AVERAGE and STDEV over',kntrj,'trajectories'
write(n21,'(a1,3x a7,3x,a17,3x,a18,3x,a16)')'#','com/Ang','<en LJ  >/kJ/mol','stdev_com/Ang','stdev_LJ/kJ/mol'
write(n21,*)

write(n22,'(a24,1x,i6,1x,a12)')'# AVERAGE and STDEV over',kntrj,'trajectories'
write(n22,'(a1,3x a7,3x,a17,3x,a18,3x,a17)')'#','com/Ang','<en pot >/kJ/mol','stdev_com/Ang','stdev_pot/kJ/mol'
write(n22,*)

!Loop over number of trajectories for calculating mean potential map [electrostatic and vdW interactions]
do ktr = 1,kntrj                              !Loop over number of trajectories
   IF(inputformat.eq.'GROMACS'.and.ktr > 1)call f77_molfile_init                          !This corrected the problem of running this routine after other analyses
   write(*,*)'Starting analysis of trajectory ', ktr,' out of ',kntrj
   write(*,*)

!Label different trajectories input files: md_trj01.xtc, md_trj02.xtc, etc
   IF(kntrj>1)THEN
!   if(ktr==1)trjfile='md_trj01.xtc'
!   if(ktr==2)trjfile='md_trj02.xtc'
!   etc

      write (klab,'(I2.2)') ktr                             ! converting integer (ktr) to string (klab) using an 'internal file'
      trjfile='md_trj'//trim(klab)//'.xtc'
!check      write(*,*) trjfile 
   ENDIF

!time counter of electrostatic interactions (1,2,3,...) with 1,2,3 separated by ncontsamp
   it = 0
   jt = 0
   kt = 0
   icm = 0
   
!xtc2fortran_trj
   
   IF(inputformat.eq.'GROMACS')THEN
      infile = TRIM(trjfile)
      intype = 'auto'
      npart  = -1
      handle(1) = -1
      handle(2) = -1
      handle(3) = -1
      handle(4) = -1
!set up everything and register all static plugins
!NG      call f77_molfile_init

      call f77_molfile_open_read(handle(1),npart,infile,intype)

      if (handle(1).lt.0) then
         print*,'file type unknown or not registered'
         stop
      else
         print*,'file successfully opened:'
         print*,'handle:',handle(1)
         print*,'npart: ',npart
         print*,'nsteps:',nstep
      end if
        
      if(npart /= natms)then
        write(*,*)'error - number of atoms is wrong in solv_inp.dat'
        stop
      endif  
   ENDIF    
   
!end xtc2fortran_trj  

!Start trajectory analysis
   
!Version 18   do j = 1,nstep                                               !Loop over number of steps of trajectory
   do j = 0,nstep
      if(mod(j,500)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/    
      IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
         if(nens==1)then
            read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
         elseif(nens==0)then
            cell(1) = cube(1)
            cell(2) = cube(2)
            cell(3) = cube(3)
         endif   
!read solute and water coordinates
         read(ninput,*)natcheck
         if(natcheck.ne.natms)then
            write(*,*)'error - number of atoms from solv_inp is wrong'
            stop
         endif   
         do i = 1,natms
            read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!NG         if(j==1)atom_symb(i)=TRIM(atomname(i))
!check         if(j==1.and.atom_symb(i)/='H')write(*,*)atom_symb(i)
         end do
      ELSEIF(inputformat.eq.'GROMACS')THEN
      
!xtc2fortran_trj 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
       status = 1   ! status=1 on entry means read
       call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
       cell(1)=dble(box(1))
       cell(2)=dble(box(2))
       cell(3)=dble(box(3))
       kp=1
!       write(*,*)cell
       do i = 1,npart*3,3
          x(kp)=xyz(i)
          y(kp)=xyz(i+1)
          z(kp)=xyz(i+2)
          kp = kp + 1
       end do   
   
!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
      ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 


!Calculate proteins centre of mass (COM) distance along time for each trajectory
!Version 18   if(j==1.or.mod(j,ncontsamp)==0)then 
      if(mod(j,ncontsamp)==0)then 
!check         write(*,*)'calculating com distance...trj ',ktr,' step ',j
         icm = icm + 1
         kt = kt + 1                                                                              !pseudo-time counter
         nstart = 1
         nend   = natmsol/2
         call sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol_I,cmysol_I,cmzsol_I)
         nstart = natmsol/2 + 1
         nend   = natmsol
         call sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol_II,cmysol_II,cmzsol_II)
         dcmx = cmxsol_I - cmxsol_II
         dcmy = cmysol_I - cmysol_II
         dcmz = cmzsol_I - cmzsol_II
         if(kMIC==1)then
            dcmx = dcmx - dnint(dcmx/cell(1))*cell(1)
            dcmy = dcmy - dnint(dcmy/cell(2))*cell(2)
            dcmz = dcmz - dnint(dcmz/cell(3))*cell(3)
         endif
         dcmr2 = dcmx**2 + dcmy**2 + dcmz**2
         dcmr  = dsqrt(dcmr2)
         dcm_trj(kt) = dcm_trj(kt) + dcmr                     !accumulate com distance at each time for all trajectories
!Multiple trjs         
         if(kntrj > 1)then
            dcm_trj_SMD(ktr,kt) = dcmr                        !accumulate com distance at each time for each trajectory
         endif   
!check for a single trajectory         
         PSECS = TSTEP*dfloat(icm-1)*dfloat(ncontsamp)
         NSECS = 1e-3*PSECS                                   !nsecs
         if(ktr==1)write(n1,'(e14.5,3x,f12.4)')NSECS,dcmr                     !print check for a single trajectory - to compare with the mean over trajectories 
!NG         if(ktr==1)write(n1,'(e19.3,3x,f12.4,2x,a6,1x,i6)')NSECS,dcmr,'trj = ',ktr  
      endif 

!End COM calculation


!Calculate electrostatic interactions of residue contact pairs along time
      if(j==0)write(n0,*)
      if(j==0)write(n0,*)'trajectory number = ',ktr
      if(j==0)write(n0,*)
      if(j==0)write(n0,*)'Electrostatic interactions at time ZERO'
      if(j==0)write(n0,*)
      if(ktr==1.and.j==0)write(n0,'(2(a5,2x,a4,2x,a12,2x,a9,7x),3x,2(2x,a4),1x,a3,2x,a11)')&
      'atom1','RES1','N_atoms_RES1','name_RES1','atom2','RES2','N_atoms_RES2','name_RES2','Q1/e','Q2/e','r/A','COUL/kJ/mol'
      if(j==0)write(n0,*)     
!Version 18   if(j==1.or.mod(j,ncontsamp)==0)then       
      if(mod(j,ncontsamp)==0)then
         ndist = 0                            !number of times the distances are accumulated
         it = it +1                           !pseudo-time counter
         k1 = 0                               !Chains ABCD atomic index
         do ic1=1,nb_res/2                    !residues of ABCD
            do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
               k1 = k1 + 1
               k2 = natmsol/2                 !Chains EFGH atomic index                  
               do ic2=nb_res/2+1,nb_res       !residues of EFGH
                  do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                     k2 = k2 + 1
                     if(res_cont_pair_t0(ic1,ic2)==1)then                  !calculate electrostatic potential of the pair
                        dx = x(k1)-x(k2)
                        dy = y(k1)-y(k2)
                        dz = z(k1)-z(k2)
                        if(kMIC==1)then
                           dx = dx - dnint(dx/cell(1))*cell(1)
                           dy = dy - dnint(dy/cell(2))*cell(2)
                           dz = dz - dnint(dz/cell(3))*cell(3)
                        endif 
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        Q1 = chrg(k1)*ELCH
                        Q2 = chrg(k2)*ELCH
                        PCoul_pair(ic1,ic2,it) = PCoul_pair(ic1,ic2,it) + enconv_kjmol*PMTTV*Q1*Q2/dr                  !accumulate over different trajectories
                        pot_pair(ic1,ic2,it) = pot_pair(ic1,ic2,it) + enconv_kjmol*PMTTV*Q1*Q2/dr                      !accumulate over different trajectories
                        
!Version 32              
                        if(atomname_d(k1)=='CA'.and.atomname_d(k2)=='CA')then
                           rdist_pair(ic1,ic2,it) = rdist_pair(ic1,ic2,it) + dr
!                            write(*,*)j,k1,k2,atomname_d(k1),atomname_d(k2)
                           ndist = ndist + 1
                        elseif(atomname_d(k1)=='FE'.and.atomname_d(k2)=='CA')then
                           rdist_pair(ic1,ic2,it) = rdist_pair(ic1,ic2,it) + dr
!                            write(*,*)j,k1,k2,atomname_d(k1),atomname_d(k2)
                            ndist = ndist + 1
                        elseif(atomname_d(k1)=='CA'.and.atomname_d(k2)=='FE')then
                           rdist_pair(ic1,ic2,it) = rdist_pair(ic1,ic2,it) + dr
!                           write(*,*)j,k1,k2,atomname_d(k1),atomname_d(k2) 
                           ndist = ndist + 1
                        elseif(atomname_d(k1)=='FE'.and.atomname_d(k2)=='FE')then
                           rdist_pair(ic1,ic2,it) = rdist_pair(ic1,ic2,it) + dr
!                           write(*,*)j,k1,k2,atomname_d(k1),atomname_d(k2)
                           ndist = ndist + 1
                        endif   
!
                        
!End v32        
                        pot_thr_ftrj(ktr,ic1,ic2,it) = pot_thr_ftrj(ktr,ic1,ic2,it) + enconv_kjmol*PMTTV*Q1*Q2/dr    !en used to probe threshold
!Multiple trjs - SMD                        
                        if(kntrj > 1)then
                           PCoul_SMD(ktr,ic1,ic2,it) = PCoul_SMD(ktr,ic1,ic2,it) + enconv_kjmol*PMTTV*Q1*Q2/dr
                           pot_SMD(ktr,ic1,ic2,it)   = pot_SMD(ktr,ic1,ic2,it)   + enconv_kjmol*PMTTV*Q1*Q2/dr
                        endif
                        
                        if(ktr==1.and.j==0)write(n0,'(2(i9,1x,i7,1x,i2,1x,a3,14x),5x,3(F9.3,1x),e12.3)') &
                        k1,ic1,Natresid(ic1),name_res(ic1),k2,ic2,Natresid(ic2),name_res(ic2),chrg(k1),chrg(k2),dr,PCoul_pair(ic1,ic2,it)         
                     endif
                  end do
               end do   
            end do
         end do
      endif                          !end electrostatics sampling
      write(*,*)'step',j,'Number of distances calculated',ndist
   
!End electrostatics   


!Calculate van der Waals interactions of residue contact pairs along time
      if(j==0)write(n0,*)
      if(j==0)write(n0,*)'Van der Waals interactions at time ZERO'
      if(j==0)write(n0,*)
      if(ktr==1.and.j==0)write(n0,'(2(a5,2x,a4,2x,a12,2x,a9,4x),2(1x,a15),2(1x,a16),1x,a3,1x,a10)')&
      'atom1','RES1','N_atoms_RES1','name_RES1','atom2','RES2','N_atoms_RES2','name_RES2','LJC6/kJnm6/mol','LJC6/kJnm12/mol',&
      'LJC12/kJnm6/mol','LJC12/kJnm12/mol','r/A','vdW/kJ/mol'
      if(j==0)write(n0,*)    
!Version 18   if(j==1.or.mod(j,ncontsamp)==0)then       
      if(mod(j,ncontsamp)==0)then
         jt = jt +1                           !pseudo-time counter
         k1 = 0                               !Chains ABCD atomic index
         do ic1=1,nb_res/2                    !residues of ABCD
            do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
               k1 = k1 + 1
               k2 = natmsol/2                 !Chains EFGH atomic index                  
               do ic2=nb_res/2+1,nb_res       !residues of EFGH
                  do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                     k2 = k2 + 1
                     if(res_cont_pair_t0(ic1,ic2)==1)then                  !calculate Lennard-Jones potential of the pair
                        dx = x(k1)-x(k2)
                        dy = y(k1)-y(k2)
                        dz = z(k1)-z(k2)
                        if(kMIC==1)then
                           dx = dx - dnint(dx/cell(1))*cell(1)
                           dy = dy - dnint(dy/cell(2))*cell(2)
                           dz = dz - dnint(dz/cell(3))*cell(3)
                        endif
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
!Apply geometric rule - units [kJ/mol nm**6]  and [kJ/mol nm**12] - convert to [kJ/mol Ang**6]  and [kJ/mol Ang**12]                   
!kJ/mol Ang**6                        C6  = 1.0d6*dsqrt(C6kjmol(k1)*C6kjmol(k2))  
!kJ/mol Ang**12                       C12 = 1.0d12*dsqrt(C12kjmol(k1)*C12kjmol(k2)) 
!r Ang                        PLJ_pair(ic1,ic2,jt) = PLJ_pair(ic1,ic2,jt) + (C12/dr**12.0 - C6/dr**6.0)                      !accumulate over different trajectories
!r Ang                        pot_pair(ic1,ic2,jt) = pot_pair(ic1,ic2,jt) + (C12/dr**12.0 - C6/dr**6.0)                      !accumulate over different trajectories
                        C6  = dsqrt(C6kjmol(k1)*C6kjmol(k2))                    !kJ/mol nm**6  
                        C12 = dsqrt(C12kjmol(k1)*C12kjmol(k2))                  !kJ/mol nm**12
!For GROMOS force field C12 was modified according to GROMOS interactions rules for C12                        
                        PLJ_pair(ic1,ic2,jt) = PLJ_pair(ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)  !r nm !accumulate over =/ trajs
                        pot_pair(ic1,ic2,jt) = pot_pair(ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)  !r nm !accumulate over =/ trajs
!en used to probe threshold
                        pot_thr_ftrj(ktr,ic1,ic2,jt) = pot_thr_ftrj(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0) 
!Version 20                        
                        PLJ_thr_ftrj(ktr,ic1,ic2,jt) = PLJ_thr_ftrj(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)
!End V20                        
!Multiple trjs                       
                        if(kntrj > 1)then
                           PLJ_SMD(ktr,ic1,ic2,jt) = PLJ_SMD(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)
                           pot_SMD(ktr,ic1,ic2,jt) = pot_SMD(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)
                        endif
                        
                        if(ktr==1.and.j==0)write(n0,'(2(i9,1x,i7,1x,i2,1x,a3,14x),5x,4(e12.4,1x),1x,F9.3,1x,e12.3)') &
                        k1,ic1,Natresid(ic1),name_res(ic1),k2,ic2,Natresid(ic2),name_res(ic2),C6,C6_IJ(k1,k2),C12,C12_IJ(k1,k2),dr,PLJ_pair(ic1,ic2,jt)               !check
                     endif
                  end do
               end do   
            end do
         end do
      endif                          !end vdW sampling   
   
!End van der Waals   

     
!xtc2fortran_trj
      if(inputformat.eq.'GROMACS'.and.status.eq.0) then
         write(*,*)'error on reading trajectory file - step',i
         stop
      endif
!end xtc2fortran_trj    
    
   end do                !end time-step main loop

!xtc2fortran_trj
   IF(inputformat.eq.'GROMACS')THEN
      call f77_molfile_finish
      write(*,*)'finished reading xtc trj...'
   ENDIF   
!end xtc2fortran_trj
   
end do                   !end trajectories loop   

!=======================================================================================================
!VERSION 16
!Find mean energy to probe threshold - first trj
!A mean value is obtained for each residue pair...

do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then
! +1
         do it=1,nstep/ncontsamp + kstr
            pot_thr_ftrj_mean(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2) + pot_thr_ftrj(1,ic1,ic2,it)
            PLJ_thr_ftrj_mean(ic1,ic2) = PLJ_thr_ftrj_mean(ic1,ic2) + PLJ_thr_ftrj(1,ic1,ic2,it)
            
         end do
! +1         
!Version 18         pot_thr_ftrj_mean(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2)/(dfloat(nstep)/dfloat(ncontsamp)+dfloat(kstr))
         pot_thr_ftrj_mean(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
         PLJ_thr_ftrj_mean(ic1,ic2) = PLJ_thr_ftrj_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
      endif
   end do
end do
!=======================================================================================================


!Print electrostatic, van der Waals, and total intermolecular potential energy, time profiles 
!Mean [over number of trajectories] electrostatic interactions along time
do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then
         nrc = nrc + 1
         if(nres_cont_pair_t0 < maxres)write(n2,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)')   &
         '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc
         if(nres_cont_pair_t0 < maxres)write(n3,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)')   &
         '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc
         if(nres_cont_pair_t0 < maxres)write(n4,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)')   &
         '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc
!
!=======================================================================================================instantaneous en
!VERSION 16         coul_pot_t0 = PCoul_pair(ic1,ic2,1)/dfloat(kntrj)                          !value at time zero  [no longer used]
!VERSION 16         lj_pot_t0 = PLJ_pair(ic1,ic2,1)/dfloat(kntrj)                              !value at time zero  [no longer used]
!VERSION 16         pe_pot_t0(ic1,ic2) = pot_pair(ic1,ic2,1)/dfloat(kntrj)                     !value at time zero - use for Coul, LJ and Coul+vdW
         
! Version 14 - account for res-res interactions above or below threshold at times > t0
! The largest value is used to choose residues to print - if the interaction between any 2 residues is above or below a threshold 
! at any time during the trj this will be printed

!VERSION 16         do it=2,nstep/ncontsamp + 1
!VERSION 16            coul_pot_t = PCoul_pair(ic1,ic2,it)/dfloat(kntrj)
!VERSION 16            if((coul_pot_t>0.0d0).and.(coul_pot_t>coul_pot_t0))then
!VERSION 16               coul_pot_t0 = coul_pot_t
!VERSION 16            elseif((coul_pot_t<0.0d0).and.(coul_pot_t<coul_pot_t0))then
!VERSION 16               coul_pot_t0 = coul_pot_t
!VERSION 16            endif
!VERSION 16            
!VERSION 16            lj_pot_t = PLJ_pair(ic1,ic2,it)/dfloat(kntrj)           
!VERSION 16            if((lj_pot_t>0.0d0).and.(lj_pot_t>lj_pot_t0))then
!VERSION 16               lj_pot_t0 = lj_pot_t
!VERSION 16            elseif((lj_pot_t<0.0d0).and.(lj_pot_t<lj_pot_t0))then
!VERSION 16               lj_pot_t0 = lj_pot_t   
!VERSION 16            endif
!VERSION 16            
!VERSION 16            pe_pot_t = pot_pair(ic1,ic2,it)/dfloat(kntrj)
!VERSION 16            if((pe_pot_t>0.0d0).and.(pe_pot_t>pe_pot_t0(ic1,ic2)))then
!VERSION 16               pe_pot_t0(ic1,ic2) = pe_pot_t
!VERSION 16            elseif((pe_pot_t<0.0d0).and.(pe_pot_t<pe_pot_t0(ic1,ic2)))then   
!VERSION 16               pe_pot_t0(ic1,ic2) = pe_pot_t
!VERSION 16            endif
!VERSION 16         end do
         
!Version 16 - use the mean en instead of the instantaneous energies

         pe_pot_t0(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2)
!=======================================================================================================end instantaneous en         
         

!PRINT HEADERS OF THRESHOLD FILES         
!Coulomb - threshold header
!Version 15 -print the same pairs for Coul. LJ and potential - use a single threshold en
!Version 15        if((coul_pot_t0 < Coul_th_min).or.(coul_pot_t0 > Coul_th_max))then
         if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then
            nrc_C = nrc_C + 1
            write(n5,*)
            write(n5,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_C
!van der Waals - threshold header       
!Version 15        if((lj_pot_t0 < LJ126_th_min).or.(lj_pot_t0 > LJ126_th_max))then
            nrc_LJ = nrc_LJ + 1
            write(n6,*)
            write(n6,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_LJ
!Intermolecular potential - sum of Coulomb and van der Waals - threshold header         
            nrc_p = nrc_p + 1
            write(n7,*)
            write(n7,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_p         
         endif
         
         
!Number of contacts with energies within the threshold window [Coul_th_min,Coul_th_max]         
!Version 15         if((coul_pot_t0 > Coul_th_min).and.(coul_pot_t0 < Coul_th_max))then
         if((pe_pot_t0(ic1,ic2) > Coul_th_min).and.(pe_pot_t0(ic1,ic2) < Coul_th_max))then
               ncont_Coul_th = ncont_Coul_th + 1
!Version 15         if((lj_pot_t0 > LJ126_th_min).and.(lj_pot_t0 < LJ126_th_max))then
               ncont_LJ_th = ncont_LJ_th + 1
               ncont_th = ncont_th + 1
         endif
         
!================================================================================================================================         
!Version 20

!vDW - threshold
         if((PLJ_thr_ftrj_mean(ic1,ic2) < LJ126_th_min).or.(PLJ_thr_ftrj_mean(ic1,ic2) > LJ126_th_max))then
            nrc_C_ = nrc_C_ + 1
            write(n25,*)
            write(n25,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_C_
!van der Waals - threshold header       
!Version 15        if((lj_pot_t0 < LJ126_th_min).or.(lj_pot_t0 > LJ126_th_max))then
            nrc_LJ_ = nrc_LJ_ + 1
            write(n26,*)
            write(n26,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_LJ_
!Intermolecular potential - sum of Coulomb and van der Waals - threshold header         
            nrc_p_ = nrc_p_ + 1
            write(n27,*)
            write(n27,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_p_
         endif
         
         
!Number of contacts with energies within the threshold window [Coul_th_min,Coul_th_max]         
!Version 15         if((coul_pot_t0 > Coul_th_min).and.(coul_pot_t0 < Coul_th_max))then
         if((PLJ_thr_ftrj_mean(ic1,ic2) > LJ126_th_min).and.(PLJ_thr_ftrj_mean(ic1,ic2) < LJ126_th_max))then
               ncont_Coul_th_ = ncont_Coul_th_ + 1
!Version 15         if((lj_pot_t0 > LJ126_th_min).and.(lj_pot_t0 < LJ126_th_max))then
               ncont_LJ_th_ = ncont_LJ_th_ + 1
               ncont_th_ = ncont_th_ + 1
         endif

!End Version 20
!================================================================================================================================ 
         
!============================================================ SMD           
!SMD sum of energies for each trj  - for stdev calculation           
         if(kntrj > 1)then        
            do ktr = 1,kntrj
               do kt=1,nstep/ncontsamp + kstr
                  sum_Coul_en_all_SMD(ktr,kt) = sum_Coul_en_all_SMD(ktr,kt) + PCoul_SMD(ktr,ic1,ic2,kt)
                  sum_LJ_en_all_SMD(ktr,kt)   = sum_LJ_en_all_SMD(ktr,kt)   + PLJ_SMD(ktr,ic1,ic2,kt)
                  sum_pot_en_all_SMD(ktr,kt)  = sum_pot_en_all_SMD(ktr,kt)  + pot_SMD(ktr,ic1,ic2,kt)
                  if((pe_pot_t0(ic1,ic2) > Coul_th_min).and.(pe_pot_t0(ic1,ic2) < Coul_th_max))then
                     sum_Coul_en_thr_SMD(ktr,kt) = sum_Coul_en_thr_SMD(ktr,kt) + PCoul_SMD(ktr,ic1,ic2,kt)
                     sum_LJ_en_thr_SMD(ktr,kt)   = sum_LJ_en_thr_SMD(ktr,kt)   + PLJ_SMD(ktr,ic1,ic2,kt)
                     sum_pot_en_thr_SMD(ktr,kt)  = sum_pot_en_thr_SMD(ktr,kt)  + pot_SMD(ktr,ic1,ic2,kt)
                  endif
               end do
            end do
         endif
!============================================================ END SMD                    
         
!Write mean potential [over trajectories] along time for each residue pair         
!NG         do it=1,nstep/ncontsamp + 1                        !This is not time - positions in the array - each position is separated in time by ncontsamp 
! +1
         do it=1,nstep/ncontsamp + kstr
            PSECS = TSTEP*dfloat(it-1)*dfloat(ncontsamp)                 !psecs
            NSECS = 1e-3*PSECS                                           !nanosecs
            coul_pot = PCoul_pair(ic1,ic2,it)/dfloat(kntrj)
            lj_pot   = PLJ_pair(ic1,ic2,it)/dfloat(kntrj)
            pe_pot   = pot_pair(ic1,ic2,it)/dfloat(kntrj)
            
            if(nres_cont_pair_t0 < maxres)write(n2,'(e14.5,3x,f12.4)')NSECS,PCoul_pair(ic1,ic2,it)/dfloat(kntrj)
            if(nres_cont_pair_t0 < maxres)write(n3,'(e14.5,3x,f12.4)')NSECS,PLJ_pair(ic1,ic2,it)/dfloat(kntrj)
            if(nres_cont_pair_t0 < maxres)write(n4,'(e14.5,3x,f12.4)')NSECS,pot_pair(ic1,ic2,it)/dfloat(kntrj)
            
!calculate the sum of the potential energy of all contacts with an interaction energy within the threshold window [Coul_th_min,Coul_th_max]

!Version 15            if((coul_pot_t0 > Coul_th_min).and.(coul_pot_t0 < Coul_th_max))then
            if((pe_pot_t0(ic1,ic2) > Coul_th_min).and.(pe_pot_t0(ic1,ic2) < Coul_th_max))then
               sum_Coul_en_t(it) = sum_Coul_en_t(it) + PCoul_pair(ic1,ic2,it)/dfloat(kntrj) 
!Version 15            if((lj_pot_t0 > LJ126_th_min).and.(lj_pot_t0 < LJ126_th_max))then
               sum_LJ_en_t(it) = sum_LJ_en_t(it) + PLJ_pair(ic1,ic2,it)/dfloat(kntrj)
               sum_pot_en_t(it) = sum_pot_en_t(it) + pot_pair(ic1,ic2,it)/dfloat(kntrj)
            endif
            
!calculate the sum of the potential energy of all contacts

            sum_Coul_en_all(it) = sum_Coul_en_all(it) + PCoul_pair(ic1,ic2,it)/dfloat(kntrj)
            sum_LJ_en_all(it)   = sum_LJ_en_all(it)   + PLJ_pair(ic1,ic2,it)/dfloat(kntrj)
            sum_pot_en_all(it)  = sum_pot_en_all(it)  + pot_pair(ic1,ic2,it)/dfloat(kntrj)
                     
!Print threshold contact energies            
!Version 15            if((coul_pot_t0 < Coul_th_min).or.(coul_pot_t0 > Coul_th_max))then
            if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then
               write(n5,'(e14.5,3x,f12.4)')NSECS,coul_pot
!Version 15            if((lj_pot_t0 < LJ126_th_min).or.(lj_pot_t0 > LJ126_th_max))then
               write(n6,'(e14.5,3x,f12.4)')NSECS,lj_pot
               write(n7,'(e14.5,3x,f12.4)')NSECS,pe_pot
            endif
!Version 20 ========================================================           
            if((PLJ_thr_ftrj_mean(ic1,ic2) < LJ126_th_min).or.(PLJ_thr_ftrj_mean(ic1,ic2) > LJ126_th_max))then
               write(n25,'(e14.5,3x,f12.4)')NSECS,coul_pot
               write(n26,'(e14.5,3x,f12.4)')NSECS,lj_pot
               write(n27,'(e14.5,3x,f12.4)')NSECS,pe_pot
            endif  
            
!End Version 20 ====================================================
!
         end do  
         write(n2,*)
         write(n3,*)
         write(n4,*)
      endif
   end do                            !End loop over residues of EFGH
end do                               !End loop over residues of ABCD

!======================================================
!check sums
!NG_check     write(*,*)
!NG_check     write(*,*)'check sum'
!NG_check     write(*,*)
!NG_check     if(kntrj > 1)then
!NG_check        do kt=1,nstep/ncontsamp + 1
!NG_check           do ktr = 1,kntrj
!NG_check               sum_check_Coul(kt) = sum_check_Coul(kt) + sum_Coul_en_all_SMD(ktr,kt)
!NG_check           end do
!NG_check        end do   
!NG_check        sum_check_Coul = sum_check_Coul/dfloat(kntrj)                                   !this must be equal to sum_Coul_en_all
!NG_check        do it=1,nstep/ncontsamp + 1
!NG_check           write(*,*)sum_Coul_en_all(it),sum_check_Coul(it)
!NG_check        end do
!NG_check     endif
!NG_check     write(*,*)
!NG_check     write(*,*)'Ended check sum'
!NG_check     write(*,*)
!======================================================

!Print mean [over number of trajectories] energies of contacts with energy within threshold [Coul_th_min - Coul_th]
write(n5,*)
write(n5,'(a59,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean Coul energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n5,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n5,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_Coul_th
write(n5,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_Coul_th)/dfloat(nres_cont_pair_t0)*100.0

write(n6,*)
write(n6,'(a57,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean LJ energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n6,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n6,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_LJ_th
write(n6,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_LJ_th)/dfloat(nres_cont_pair_t0)*100.0

write(n7,*)
write(n7,'(a58,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean pot energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n7,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n7,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_th
write(n7,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_th)/dfloat(nres_cont_pair_t0)*100.0

!Version 20
write(n25,*)
write(n25,'(a59,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean Coul energy of contacts with vdW energy in the range','[',LJ126_th_min,'-',LJ126_th_max,']'
write(n25,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n25,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_Coul_th
write(n25,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_Coul_th)/dfloat(nres_cont_pair_t0)*100.0

write(n26,*)
write(n26,'(a57,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean LJ energy of contacts with vdW energy in the range','[',LJ126_th_min,'-',LJ126_th_max,']'
write(n26,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n26,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_LJ_th
write(n26,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_LJ_th)/dfloat(nres_cont_pair_t0)*100.0

write(n27,*)
write(n27,'(a58,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean pot energy of contacts with vdW energy in the range','[',LJ126_th_min,'-',LJ126_th_max,']'
write(n27,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n27,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_th
write(n27,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_th)/dfloat(nres_cont_pair_t0)*100.0

!End Version 20

!NG do it=1,nstep/ncontsamp + 1                        !This is not time - positions in the array - each position is separated in time by ncontsamp 
! +1
do it=1,nstep/ncontsamp + kstr   
   PSECS = TSTEP*dfloat(it-1)*dfloat(ncontsamp)                 !psecs
   NSECS = 1e-3*PSECS                                           !nanosecs
   write(n5,'(e14.5,3x,f12.4)')NSECS,sum_Coul_en_t(it)   
   write(n6,'(e14.5,3x,f12.4)')NSECS,sum_LJ_en_t(it)
   write(n7,'(e14.5,3x,f12.4)')NSECS,sum_pot_en_t(it)
   write(n25,'(e14.5,3x,f12.4)')NSECS,sum_Coul_en_t(it)   
   write(n26,'(e14.5,3x,f12.4)')NSECS,sum_LJ_en_t(it)
   write(n27,'(e14.5,3x,f12.4)')NSECS,sum_pot_en_t(it)
end do 


!Print mean (over number of trajectories) energies of all contacts 
write(n5,*)
write(n5,'(a34)')'# Mean Coul energy of all contacts'
write(n5,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n6,*)
write(n6,'(a32)')'# Mean LJ energy of all contacts' 
write(n6,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n7,*)
write(n7,'(a33)')'# Mean pot energy of all contacts'
write(n7,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

!Version 20
write(n25,*)
write(n25,'(a34)')'# Mean Coul energy of all contacts'
write(n25,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n26,*)
write(n26,'(a32)')'# Mean LJ energy of all contacts' 
write(n26,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n27,*)
write(n27,'(a33)')'# Mean pot energy of all contacts'
write(n27,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0



!End Version 20

!NG do it=1,nstep/ncontsamp + 1                        !This is not time - positions in the array - each position is separated in time by ncontsamp 
! +1
do it=1,nstep/ncontsamp + kstr   
   PSECS = TSTEP*dfloat(it-1)*dfloat(ncontsamp)                 !psecs
   NSECS = 1e-3*PSECS                                           !nanosecs   
   write(n5,'(e14.5,3x,f12.4)')NSECS,sum_Coul_en_all(it)
   write(n6,'(e14.5,3x,f12.4)')NSECS,sum_LJ_en_all(it)
   write(n7,'(e14.5,3x,f12.4)')NSECS,sum_pot_en_all(it)
   write(n25,'(e14.5,3x,f12.4)')NSECS,sum_Coul_en_all(it)
   write(n26,'(e14.5,3x,f12.4)')NSECS,sum_LJ_en_all(it)
   write(n27,'(e14.5,3x,f12.4)')NSECS,sum_pot_en_all(it) 
end do 

!Calculate the MEAN electrostatic and van der Waals potential [average over time and trajectories] for each residue pair
!A mean value is obtained for each residue pair...

do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then                  
!         write(n8,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!NG         do it=1,nstep/ncontsamp+1                        !This is not time - positions in the array - each position is separated in time by ncontsamp 
! +1
         do it=1,nstep/ncontsamp + kstr
            PCoul_pair_mean(ic1,ic2) = PCoul_pair_mean(ic1,ic2) + PCoul_pair(ic1,ic2,it)/dfloat(kntrj)
            PLJ_pair_mean(ic1,ic2)   = PLJ_pair_mean(ic1,ic2)   + PLJ_pair(ic1,ic2,it)/dfloat(kntrj)
            pot_pair_mean(ic1,ic2)   = pot_pair_mean(ic1,ic2)   + pot_pair(ic1,ic2,it)/dfloat(kntrj)
!Version 32
! Mean distance between residues
            rdist_pair_mean(ic1,ic2) = rdist_pair_mean(ic1,ic2) + rdist_pair(ic1,ic2,it)/dfloat(kntrj)
!end v32
         end do
! +1
!Version 18
         PCoul_pair_mean(ic1,ic2) = PCoul_pair_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
         PLJ_pair_mean(ic1,ic2)   = PLJ_pair_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
         pot_pair_mean(ic1,ic2)   = pot_pair_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
!Version 32
         rdist_pair_mean(ic1,ic2) = rdist_pair_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
!end v32
      endif
   end do
end do

!Calculate the MEAN [over time] electrostatic and van der Waals potential for the sum of the energy of all residues and residues within the threshold window
! +1
do it=1,nstep/ncontsamp + kstr
   sum_Coul_en_t_mean   = sum_Coul_en_t_mean + sum_Coul_en_t(it)
   sum_LJ_en_t_mean     = sum_LJ_en_t_mean   + sum_LJ_en_t(it)
   sum_pot_en_t_mean    = sum_pot_en_t_mean  + sum_pot_en_t(it)
   
   sum_Coul_en_all_mean = sum_Coul_en_all_mean + sum_Coul_en_all(it)
   sum_LJ_en_all_mean   = sum_LJ_en_all_mean   + sum_LJ_en_all(it)
   sum_pot_en_all_mean  = sum_pot_en_all_mean  + sum_pot_en_all(it)
end do
!Version 18
sum_Coul_en_t_mean   = sum_Coul_en_t_mean /dfloat(nstep/ncontsamp + kstr)  
sum_LJ_en_t_mean     = sum_LJ_en_t_mean   /dfloat(nstep/ncontsamp + kstr)
sum_pot_en_t_mean    = sum_pot_en_t_mean  /dfloat(nstep/ncontsamp + kstr)
                                          
sum_Coul_en_all_mean = sum_Coul_en_all_mean/dfloat(nstep/ncontsamp + kstr)
sum_LJ_en_all_mean   = sum_LJ_en_all_mean  /dfloat(nstep/ncontsamp + kstr)
sum_pot_en_all_mean  = sum_pot_en_all_mean /dfloat(nstep/ncontsamp + kstr)

!Calculate the STANDARD DEVIATION [over time] of the mean electrostatic and van der Waals potential for each residue pair
!The STANDARD DEVIATION IS FOR TIME NOT NUMBER OF TRAJECTORIES - EACH VALUE AT EACH TIME IS HOWEVER AN AVERAGE OVER TRAJECTORIES
!A mean value and standard deviation is obtained for each residue pair 
do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then                  
!NG         do it=1,nstep/ncontsamp+1                        !This is not time - positions in the array - each position is separated in time by ncontsamp 
! +1
         do it=1,nstep/ncontsamp + kstr
            PCoul_pair_stdev(ic1,ic2) = PCoul_pair_stdev(ic1,ic2) + (PCoul_pair(ic1,ic2,it)/dfloat(kntrj)-PCoul_pair_mean(ic1,ic2))**2.0
            PLJ_pair_stdev(ic1,ic2)   = PLJ_pair_stdev(ic1,ic2)   + (PLJ_pair(ic1,ic2,it)/dfloat(kntrj)  -PLJ_pair_mean(ic1,ic2))**2.0
            pot_pair_stdev(ic1,ic2)   = pot_pair_stdev(ic1,ic2)   + (pot_pair(ic1,ic2,it)/dfloat(kntrj)  -pot_pair_mean(ic1,ic2))**2.0
            rdist_pair_stdev(ic1,ic2) = rdist_pair_stdev(ic1,ic2) + (rdist_pair(ic1,ic2,it)/dfloat(kntrj)-rdist_pair_mean(ic1,ic2))**2.0
         end do  
!Version 18         
         PCoul_pair_stdev(ic1,ic2) = dsqrt(PCoul_pair_stdev(ic1,ic2)/dfloat(nstep/ncontsamp)) 
         PLJ_pair_stdev(ic1,ic2)   = dsqrt(PLJ_pair_stdev(ic1,ic2)/dfloat(nstep/ncontsamp))
         pot_pair_stdev(ic1,ic2)   = dsqrt(pot_pair_stdev(ic1,ic2)/dfloat(nstep/ncontsamp))                     !n-1 (standard deviation finite sample +1 -1 = 0) 
         rdist_pair_stdev(ic1,ic2) = dsqrt(rdist_pair_stdev(ic1,ic2)/dfloat(nstep/ncontsamp))
      endif
   end do
end do

!Calculate the STANDARD DEVIATION [over time] of the mean electrostatic and van der Waals potential for the sum of the energy of all residues and residues within threshold window
! +1
do it=1,nstep/ncontsamp + kstr
   sum_Coul_en_thr_stdev = sum_Coul_en_thr_stdev + (sum_Coul_en_t(it)-sum_Coul_en_t_mean)**2.0
   sum_LJ_en_thr_stdev   = sum_LJ_en_thr_stdev   + (sum_LJ_en_t(it)-sum_LJ_en_t_mean)**2.0
   sum_pot_en_thr_stdev  = sum_pot_en_thr_stdev  + (sum_pot_en_t(it)-sum_pot_en_t_mean)**2.0
   sum_Coul_en_all_stdev = sum_Coul_en_all_stdev + (sum_Coul_en_all(it)-sum_Coul_en_all_mean)**2.0
   sum_LJ_en_all_stdev   = sum_LJ_en_all_stdev   + (sum_LJ_en_all(it)-sum_LJ_en_all_mean)**2.0
   sum_pot_en_all_stdev  = sum_pot_en_all_stdev  + (sum_pot_en_all(it)-sum_pot_en_all_mean)**2.0
end do
!Version 18
sum_Coul_en_thr_stdev = dsqrt(sum_Coul_en_thr_stdev/dfloat(nstep/ncontsamp))
sum_LJ_en_thr_stdev   = dsqrt(sum_LJ_en_thr_stdev  /dfloat(nstep/ncontsamp))
sum_pot_en_thr_stdev  = dsqrt(sum_pot_en_thr_stdev /dfloat(nstep/ncontsamp))
sum_Coul_en_all_stdev = dsqrt(sum_Coul_en_all_stdev/dfloat(nstep/ncontsamp))
sum_LJ_en_all_stdev   = dsqrt(sum_LJ_en_all_stdev  /dfloat(nstep/ncontsamp))
sum_pot_en_all_stdev  = dsqrt(sum_pot_en_all_stdev /dfloat(nstep/ncontsamp))

ipC = 0
ipthC = 0
ipLJ = 0
ipthLJ = 0
ipp = 0
ipthp = 0

!Print mean electrostatic potential and respective standard deviations
do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then
!Coulomb
         ipC = ipC + 1
         write(n8,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!NG         write(n8,'(e12.4,3x,a3,3x,e12.4)')PCoul_pair_mean(ic1,ic2),'+/-',PCoul_pair_stdev(ic1,ic2)
         write(n8,'(i9,3x,e12.4,3x,e12.4)')ipC,PCoul_pair_mean(ic1,ic2),PCoul_pair_stdev(ic1,ic2)
         write(n8,*)
         write(n11,'(i9,3x,e12.4,3x,e12.4)')ipC,PCoul_pair_mean(ic1,ic2),PCoul_pair_stdev(ic1,ic2)
!         write(n11,*)
!Version 15         if((PCoul_pair_mean(ic1,ic2)< Coul_th_min).or.(PCoul_pair_mean(ic1,ic2) > Coul_th_max))then
!Version 15         if((pot_pair_mean(ic1,ic2)< Coul_th_min).or.(pot_pair_mean(ic1,ic2) > Coul_th_max))then                !use mean en                    
         if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then                                  !use instantaneous en
            ipthC = ipthC + 1
            write(n15,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!NG            write(n15,'(e12.4,3x,a3,3x,e12.4)')PCoul_pair_mean(ic1,ic2),'+/-',PCoul_pair_stdev(ic1,ic2)
!            write(n15,'(i9,3x,e12.4,3x,e12.4)')ipthC,PCoul_pair_mean(ic1,ic2),PCoul_pair_stdev(ic1,ic2)
            write(n15,'(i9,3x,e15.7,3x,e15.7)')ipthC,PCoul_pair_mean(ic1,ic2),PCoul_pair_stdev(ic1,ic2)
            write(n15,*)
         endif   
         
!van der Waals       
         ipLJ = ipLJ + 1
         write(n9,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!         write(n9,'(e12.4,3x,a3,3x,e12.4)')PLJ_pair_mean(ic1,ic2),'+/-',PLJ_pair_stdev(ic1,ic2)
         write(n9,'(i9,3x,e12.4,3x,e12.4)')ipLJ,PLJ_pair_mean(ic1,ic2),PLJ_pair_stdev(ic1,ic2)
         write(n9,*)
         write(n12,'(i9,3x,e12.4,3x,e12.4)')ipLJ,PLJ_pair_mean(ic1,ic2),PLJ_pair_stdev(ic1,ic2)
!         write(n12,*)
         
!Version 15         if((PLJ_pair_mean(ic1,ic2)< LJ126_th_min).or.(PLJ_pair_mean(ic1,ic2) > LJ126_th_max))then
!Version 15         if((pot_pair_mean(ic1,ic2)< Coul_th_min).or.(pot_pair_mean(ic1,ic2) > Coul_th_max))then
         if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then                                  !use instantaneous en         
            ipthLJ = ipthLJ + 1
            write(n16,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!NG            write(n16,'(e12.4,3x,a3,3x,e12.4)')PLJ_pair_mean(ic1,ic2),'+/-',PLJ_pair_stdev(ic1,ic2)
            write(n16,'(i9,3x,e12.4,3x,e12.4)')ipthLJ,PLJ_pair_mean(ic1,ic2),PLJ_pair_stdev(ic1,ic2)
            write(n16,*)
         endif 
         
!pot energy - sum(Coul+vdW)      
         ipp = ipp + 1
         write(n10,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!NG         write(n10,'(e12.4,3x,a3,3x,e12.4)')pot_pair_mean(ic1,ic2),'+/-',pot_pair_stdev(ic1,ic2)
         write(n10,'(i9,3x,e12.4,3x,e12.4)')ipp,pot_pair_mean(ic1,ic2),pot_pair_stdev(ic1,ic2)
         write(n10,*)
         write(n13,'(i9,3x,e12.4,3x,e12.4)')ipp,pot_pair_mean(ic1,ic2),pot_pair_stdev(ic1,ic2)
!         write(n13,*)
!Version 32
         write(n28,'(i9,3x,e12.4,3x,e12.4)')ipp,rdist_pair_mean(ic1,ic2),rdist_pair_stdev(ic1,ic2)
         write(n29,'(i9,3x,e12.4,3x,e12.4)')ipp,rdist_pair_mean(ic1,ic2),pot_pair_mean(ic1,ic2)
!end v32
                  
!Version 15         if((pot_pair_mean(ic1,ic2)< Coul_th_min).or.(pot_pair_mean(ic1,ic2) > Coul_th_max))then
         if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then                                  !use instantaneous en
            ipthp = ipthp + 1
            write(n17,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x))')'#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2)
!NG            write(n17,'(e12.4,3x,a3,3x,e12.4)')pot_pair_mean(ic1,ic2),'+/-',pot_pair_stdev(ic1,ic2)
            write(n17,'(i9,3x,e12.4,3x,e12.4)')ipthp,pot_pair_mean(ic1,ic2),pot_pair_stdev(ic1,ic2)
            write(n17,*)
         endif   
      endif
   end do
end do

write(n15,*)
write(n15,'(a59,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean Coul energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n15,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n15,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_Coul_th
write(n15,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_Coul_th)/dfloat(nres_cont_pair_t0)*100.0

write(n16,*)
write(n16,'(a57,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean LJ energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n16,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n16,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_LJ_th
write(n16,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_LJ_th)/dfloat(nres_cont_pair_t0)*100.0

write(n17,*)
write(n17,'(a58,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean pot energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n17,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n17,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_th
write(n17,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_th)/dfloat(nres_cont_pair_t0)*100.0

write(n15,'(e12.4,3x,e12.4)')sum_Coul_en_t_mean,sum_Coul_en_thr_stdev 
write(n16,'(e12.4,3x,e12.4)')sum_LJ_en_t_mean,sum_LJ_en_thr_stdev     
write(n17,'(e12.4,3x,e12.4)')sum_pot_en_t_mean,sum_pot_en_thr_stdev    

write(n15,*)
write(n15,'(a34)')'# Mean Coul energy of all contacts'
write(n15,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n16,*)
write(n16,'(a32)')'# Mean LJ energy of all contacts' 
write(n16,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n17,*)
write(n17,'(a33)')'# Mean pot energy of all contacts'
write(n17,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n15,'(e12.4,3x,e12.4)')sum_Coul_en_all_mean,sum_Coul_en_all_stdev
write(n16,'(e12.4,3x,e12.4)')sum_LJ_en_all_mean,sum_LJ_en_all_stdev 
write(n17,'(e12.4,3x,e12.4)')sum_pot_en_all_mean,sum_pot_en_all_stdev 

!=========================================================================================================================== SMD 
!SMD trj averages and stdev [over trajectories] - only if more than one trj is read
nrc = 0
if(kntrj > 1)then
! +1
   do kt = 1,nstep/ncontsamp + kstr
!centre of mass   
      do ktr = 1,kntrj                                                                 !Loop over trajectories   
         stdev_dcm_SMD(kt) = stdev_dcm_SMD(kt) + (dcm_trj_SMD(ktr,kt) - dcm_trj(kt)/dfloat(kntrj))**2.0  
      end do 
      stdev_dcm_SMD(kt) = dsqrt(stdev_dcm_SMD(kt)/(dfloat(kntrj)-1.0))
   end do                                                                              !End pseudo-time loop
   
!res-res pair energies - trj stdevs
! +1
   do it = 1,nstep/ncontsamp + kstr
      do ic1=1,nb_res/2                                           !residues of ABCD
         do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
            if(res_cont_pair_t0(ic1,ic2)==1)then
               do ktr = 1,kntrj                                   !sum energies at time kt, for residues ic1 and ic2 interaction, for kntrj trjs
                  stdev_PCoul_SMD(ic1,ic2,it) = stdev_PCoul_SMD(ic1,ic2,it) + (PCoul_SMD(ktr,ic1,ic2,it) - PCoul_pair(ic1,ic2,it)/dfloat(kntrj))**2.0
                  stdev_PLJ_SMD(ic1,ic2,it)   = stdev_PLJ_SMD(ic1,ic2,it)   + (PLJ_SMD(ktr,ic1,ic2,it) - PLJ_pair(ic1,ic2,it)/dfloat(kntrj))**2.0
                  stdev_pot_SMD(ic1,ic2,it)   = stdev_pot_SMD(ic1,ic2,it)   + (pot_SMD(ktr,ic1,ic2,it) - pot_pair(ic1,ic2,it)/dfloat(kntrj))**2.0
               end do
               stdev_PCoul_SMD(ic1,ic2,it) = dsqrt(stdev_PCoul_SMD(ic1,ic2,it)/(dfloat(kntrj)-1.0))
               stdev_PLJ_SMD(ic1,ic2,it)   = dsqrt(stdev_PLJ_SMD(ic1,ic2,it)/(dfloat(kntrj)-1.0)) 
               stdev_pot_SMD(ic1,ic2,it)   = dsqrt(stdev_pot_SMD(ic1,ic2,it)/(dfloat(kntrj)-1.0))
            endif
         end do   
      end do
   end do
   
!sum of energies - threshold and all interactions - trj stdevs
!====================================================================  
   do it = 1,nstep/ncontsamp + kstr
      do ktr = 1,kntrj                                   
         stdev_Coul_all_SMD(it) = stdev_Coul_all_SMD(it) + (sum_Coul_en_all_SMD(ktr,it) - sum_Coul_en_all(it))**2.0
!         write(*,*)it,ktr,sum_Coul_en_all_SMD(ktr,it),sum_Coul_en_all(it)
         stdev_LJ_all_SMD(it)   = stdev_LJ_all_SMD(it)   + (sum_LJ_en_all_SMD(ktr,it)   - sum_LJ_en_all(it))**2.0
         stdev_pot_all_SMD(it)  = stdev_pot_all_SMD(it)  + (sum_pot_en_all_SMD(ktr,it)  - sum_pot_en_all(it))**2.0
         stdev_Coul_t_SMD(it)   = stdev_Coul_t_SMD(it)   + (sum_Coul_en_thr_SMD(ktr,it) - sum_Coul_en_t(it))**2.0
         stdev_LJ_t_SMD(it)     = stdev_LJ_t_SMD(it)     + (sum_LJ_en_thr_SMD(ktr,it)   - sum_LJ_en_t(it))**2.0
         stdev_pot_t_SMD(it)    = stdev_pot_t_SMD(it)    + (sum_pot_en_thr_SMD(ktr,it)  - sum_pot_en_t(it))**2.0
      end do
      stdev_Coul_all_SMD(it)    = dsqrt(stdev_Coul_all_SMD(it)/(dfloat(kntrj)-1.0))
      stdev_LJ_all_SMD(it)      = dsqrt(stdev_LJ_all_SMD(it)/(dfloat(kntrj)-1.0)) 
      stdev_pot_all_SMD(it)     = dsqrt(stdev_pot_all_SMD(it)/(dfloat(kntrj)-1.0))
      stdev_Coul_t_SMD(it)      = dsqrt(stdev_Coul_t_SMD(it)/(dfloat(kntrj)-1.0))
      stdev_LJ_t_SMD(it)        = dsqrt(stdev_LJ_t_SMD(it)/(dfloat(kntrj)-1.0)) 
      stdev_pot_t_SMD(it)       = dsqrt(stdev_pot_t_SMD(it)/(dfloat(kntrj)-1.0))
   end do
!====================================================================  

!Print means and stdev [over trajectories] 
!com distance
   do kt=1,nstep/ncontsamp + kstr
      PSECS = TSTEP*dfloat(kt-1)*dfloat(ncontsamp)                   !psecs
      NSECS = 1e-3*PSECS  
      write(n19,'(e14.5,3x,f12.4,5x,e14.5)')NSECS,dcm_trj(kt)/dfloat(kntrj),stdev_dcm_SMD(kt) 
   end do           
   
!energies   
!====================================================================  
   do ic1=1,nb_res/2                                           !residues of ABCD
      do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
         if(res_cont_pair_t0(ic1,ic2)==1)then            
            if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then
               nrc = nrc + 1
               write(n20,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)')   &
                  '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc
               write(n21,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)')   &
                  '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc
               write(n22,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)')   &
                  '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc                    
            endif
            do it=1,nstep/ncontsamp + kstr
               PSECS = TSTEP*dfloat(it-1)*dfloat(ncontsamp)                   !psecs
               NSECS = 1e-3*PSECS 
!threshold contacts             
               if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then                   
                  write(n20,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),PCoul_pair(ic1,ic2,it)/dfloat(kntrj),stdev_dcm_SMD(it),stdev_PCoul_SMD(ic1,ic2,it)                   
                  write(n21,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),PLJ_pair(ic1,ic2,it)/dfloat(kntrj),stdev_dcm_SMD(it),stdev_PLJ_SMD(ic1,ic2,it) 
                  write(n22,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),pot_pair(ic1,ic2,it)/dfloat(kntrj),stdev_dcm_SMD(it),stdev_pot_SMD(ic1,ic2,it)                  
               endif                          
            end do
         endif
      end do
   end do  
 
   write(n20,*)
   write(n20,'(a59,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean Coul energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
   write(n20,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
   write(n20,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_Coul_th
   write(n20,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_Coul_th)/dfloat(nres_cont_pair_t0)*100.0
   
   write(n21,*)
   write(n21,'(a57,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean LJ energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
   write(n21,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
   write(n21,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_LJ_th
   write(n21,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_LJ_th)/dfloat(nres_cont_pair_t0)*100.0
  
   write(n22,*)
   write(n22,'(a58,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean pot energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
   write(n22,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
   write(n22,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_th
   write(n22,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_th)/dfloat(nres_cont_pair_t0)*100.0
  
   do it=1,nstep/ncontsamp + kstr
      PSECS = TSTEP*dfloat(it-1)*dfloat(ncontsamp)                   !psecs
      NSECS = 1e-3*PSECS 
!en sums       
      write(n20,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),sum_Coul_en_t(it),stdev_dcm_SMD(it),stdev_Coul_t_SMD(it)
      write(n21,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),sum_LJ_en_t(it),stdev_dcm_SMD(it),stdev_LJ_t_SMD(it) 
      write(n22,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),sum_pot_en_t(it),stdev_dcm_SMD(it),stdev_pot_t_SMD(it) 
   end do
  

   write(n20,*)
   write(n20,'(a34)')'# Mean Coul energy of all contacts'
   write(n20,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
   
   write(n21,*)
   write(n21,'(a32)')'# Mean LJ energy of all contacts' 
   write(n21,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
     
   write(n22,*)
   write(n22,'(a33)')'# Mean pot energy of all contacts'
   write(n22,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
   
   do it=1,nstep/ncontsamp + 1
     PSECS = TSTEP*dfloat(it-1)*dfloat(ncontsamp)                   !psecs
     NSECS = 1e-3*PSECS 
!en sums       
     write(n20,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),sum_Coul_en_all(it),stdev_dcm_SMD(it),stdev_Coul_all_SMD(it) 
     write(n21,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),sum_LJ_en_all(it),stdev_dcm_SMD(it),stdev_LJ_all_SMD(it)
     write(n22,'(f12.4,3x,f12.4,2(4x,e14.5))')dcm_trj(it)/dfloat(kntrj),sum_pot_en_all(it),stdev_dcm_SMD(it),stdev_pot_all_SMD(it)
   end do           

!====================================================================  

endif
!End mean over trjs print

!end SMD averages   
   
deallocate(PCoul_pair)
deallocate(PCoul_pair_mean,PCoul_pair_stdev)

deallocate(PLJ_pair)
deallocate(PLJ_pair_mean,PLJ_pair_stdev)

deallocate(pot_pair)
deallocate(pot_pair_mean,pot_pair_stdev)

deallocate(sum_Coul_en_t,sum_LJ_en_t,sum_pot_en_t)
deallocate(sum_Coul_en_all,sum_LJ_en_all,sum_pot_en_all)

deallocate(dcm_trj)

if(kntrj > 1)then
   deallocate(dcm_trj_SMD)
   deallocate(PCoul_SMD)
   deallocate(PLJ_SMD)
   deallocate(pot_SMD)
   deallocate(stdev_dcm_SMD)
   deallocate(stdev_PCoul_SMD)
   deallocate(stdev_PLJ_SMD)
   deallocate(stdev_pot_SMD)
   
   deallocate(sum_Coul_en_all_SMD) 
   deallocate(sum_LJ_en_all_SMD)    
   deallocate(sum_pot_en_all_SMD)   
   deallocate(sum_Coul_en_thr_SMD)
   deallocate(sum_LJ_en_thr_SMD)  
   deallocate(sum_pot_en_thr_SMD) 
   
   deallocate(stdev_Coul_all_SMD)
   deallocate(stdev_LJ_all_SMD)
   deallocate(stdev_pot_all_SMD)  
   deallocate(stdev_Coul_t_SMD)   
   deallocate(stdev_LJ_t_SMD)     
   deallocate(stdev_pot_t_SMD)  
   
   deallocate(sum_check_Coul) 
   deallocate(pot_thr_ftrj)
   deallocate(PLJ_thr_ftrj)
   deallocate(pot_thr_ftrj_mean)
   deallocate(PLJ_thr_ftrj_mean)
endif


deallocate(xyz)
deallocate(x,y,z)

return

END SUBROUTINE potent_map



SUBROUTINE pot_map_umbrella(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,Zattype,mattype,ZMBIO,&
                            nsoltypes,kumbrella_trj,ncontsamp,kntrj,Coul_th_min,Coul_th_max,atomtype,chrg,nb_res,chain_label_res,&
                            Natresid,name_res,numb_res,kMIC,intsf_,rcont_,C6kjmol,C12kjmol,C6_IJ,C12_IJ)                      

! Computes mean electrostatic and van der Waals interactions between the contact residues for different umbrella sampling COM-COM dist. trajectories 
! The routine allows analysing umbrella sampling trajectories; the contacts list is found from the first trajectory alone in the routine prot_contact
! The routine prints the COM-COM distance and the mean electrostatic, van der Waals, and pot en (Coul+vdW) for each umbrella sampling trj corresponding to a 
! different COM-COM distance; the COM-COM and en standard deviations are also printed

    integer,intent(in)                                    :: ninput,nbox
    character(len=7),intent(in)                           :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                               :: dt                   ! time between frames (fs)
    integer,intent(in)                                    :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                               :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                                    :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites
    character(len=4),dimension(natms)                     :: atomname    
    integer,intent(in)                                    :: nsoltypes               ! Number of solute distinct atomic species (C, H, O, N etc)
    integer,intent(in)                                    :: ncontsamp               ! Protein contact map electrostatics sampling frequency
    integer,intent(in)                                    :: kntrj                   ! Number of trajectories to average - SMD 
    integer,intent(in)                                    :: kumbrella_trj           ! 0: Trjs are not from umbrella sampling 1: Trjs are from umbrella sampling
!    real(kind=8),intent(in)                              :: Coul_th,LJ126_th                ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!V15    real(kind=8),intent(in)                           :: Coul_th_min,LJ126_th_min        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!V15    real(kind=8),intent(in)                           :: Coul_th_max,LJ126_th_max        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
    real(kind=8),intent(in)                               :: Coul_th_min,Coul_th_max ! Potential energy threshold for printing res-res interactions [Min,Max]
    character(len=4),dimension(natms),intent(in)          :: atomtype 
    real(kind=4),dimension(natms),intent(in)              :: chrg     
    integer,intent(in)                                    :: nb_res                 ! Number of residues [GROMCAS]
    character(len=1),dimension(nb_res),intent(in)         :: chain_label_res        ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nb_res),intent(in)                  :: Natresid               ! Number of atoms of each residue
    character(len=4),dimension(nb_res),intent(in)         :: name_res
    integer,dimension(nb_res),intent(in)                  :: numb_res               ! Number of the residue: from 1 to the total number of residues of each chain
    real(kind=8),dimension(natmsol),intent(in)            :: C6kjmol,C12kjmol       ! Solute Lennard-Jones parameters
!Version_12    integer,dimension(nb_res,nb_res),intent(in)           :: res_cont_pair_t0
!Version_12    integer,intent(in)                                    :: nres_cont_pair_t0
!Version_12
    integer,dimension(:,:),allocatable                    :: res_cont_pair_t0         !Res-Res identity that are below some input distance
    integer                                               :: nres_cont_pair_t0        
    integer,intent(in)                                    :: kMIC                     ! apply Minimum Image Convention or not to proteins
    integer,intent(in)                                    :: intsf_                   ! Protein contact map sampling frequency
    real(kind=8),intent(in)                               :: rcont_                   ! distance for contact map - residues at a r < rcont are accounted
!End !Version_12    
       
    real(kind=8),dimension(natmsol,natmsol),intent(in)    :: C6_IJ,C12_IJ           ! Solute Lennard-Jones C6 and C12 crossed parameters kj/mol nm**6; kj/mol nm**12
    character(len=2),dimension(nsoltypes),intent(in)      :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)          :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    real(kind=8),dimension(natmsol),intent(in)            :: ZMBIO
    
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                 !atomic index (solute and solvent index)
    integer      :: n0, n1, n2, n3, n4, n5, n6, n7, n8, n9
    integer      :: i, j, k
    integer      :: ic1,ic2,iats1,iats2,k1,k2,it,jt,kt,jats1,jats2,kj1,kj2                                         !chain indexes for contact map
    integer      :: ktr,kl  
    integer      :: icm
    integer      :: ipC,ipthC,ipLJ,ipthLJ,ipp,ipthp
        
    real(kind=8)                                     :: dx,dy,dz,dr,dr2,rhi,rhf
    integer                                          :: natcheck
    logical                                          :: filex
!Version 15   
    real(kind=8),dimension(:),allocatable            :: dcm_trj_US,drcoord_trj_US
    real(kind=8),dimension(:,:),allocatable          :: dcm_trj_SMD
    real(kind=8),dimension(:,:,:),allocatable        :: PCoul_pair_US
    real(kind=8),dimension(:,:,:,:),allocatable      :: PCoul_SMD
    real(kind=8),dimension(:,:,:),allocatable        :: PLJ_pair_US
    real(kind=8),dimension(:,:,:,:),allocatable      :: PLJ_SMD
    real(kind=8),dimension(:,:,:),allocatable        :: pot_pair_US
    real(kind=8),dimension(:,:,:,:),allocatable      :: pot_SMD
       
    real(kind=8),dimension(:,:),allocatable             :: coul_pot_t0
    real(kind=8),dimension(:,:),allocatable             :: lj_pot_t0
    real(kind=8),dimension(:,:),allocatable             :: pe_pot_t0
    real(kind=8)                                        :: coul_pot_t,lj_pot_t,pe_pot_t
    real(kind=8)                                        :: dcm_mean,coul_pot,lj_pot,pe_pot
    real(kind=8)                                        :: drcoord_mean
    
    real(kind=8),dimension(:),allocatable               :: sum_Coul_en_t       !sum of the Coul en for all contacts with pe outside threshold window
    real(kind=8),dimension(:),allocatable               :: sum_LJ_en_t         !sum of the LJ en for all contacts with pe outside threshold window
    real(kind=8),dimension(:),allocatable               :: sum_pot_en_t        !sum of the pe for all contacts with pe outside threshold window
    
    real(kind=8),dimension(:),allocatable               :: sum_Coul_en_all       !sum of the Coul en for all contacts 
    real(kind=8),dimension(:),allocatable               :: sum_LJ_en_all         !sum of the LJ en for all contacts 
    real(kind=8),dimension(:),allocatable               :: sum_pot_en_all        !sum of the pe for all contacts 
!umbrella sampling trjs instantaneous en
    real(kind=8),dimension(:,:),allocatable             :: sum_Coul_en_US_thr 
    real(kind=8),dimension(:,:),allocatable             :: sum_LJ_en_US_thr 
    real(kind=8),dimension(:,:),allocatable             :: sum_pot_en_US_thr    
    real(kind=8),dimension(:,:),allocatable             :: sum_Coul_en_US_all 
    real(kind=8),dimension(:,:),allocatable             :: sum_LJ_en_US_all  
    real(kind=8),dimension(:,:),allocatable             :: sum_pot_en_US_all 
!umbrella sampling trjs standard deviations
    real(kind=8),dimension(:),allocatable                  :: stdev_dcm_US
    real(kind=8),dimension(:,:,:),allocatable              :: stdev_PCoul_US
    real(kind=8),dimension(:,:,:),allocatable              :: stdev_PLJ_US
    real(kind=8),dimension(:,:,:),allocatable              :: stdev_pot_US
    real(kind=8),dimension(:),allocatable                  :: stdev_PCoul_US_all
    real(kind=8),dimension(:),allocatable                  :: stdev_PLJ_US_all
    real(kind=8),dimension(:),allocatable                  :: stdev_pot_US_all
    real(kind=8),dimension(:),allocatable                  :: stdev_PCoul_US_thr
    real(kind=8),dimension(:),allocatable                  :: stdev_PLJ_US_thr
    real(kind=8),dimension(:),allocatable                  :: stdev_pot_US_thr 
    
    real(kind=8),dimension(:,:,:,:),allocatable         :: pot_thr_ftrj
    real(kind=8),dimension(:,:),allocatable             :: pot_thr_ftrj_mean
       
    real(kind=8)                                     :: C6,C12   
    real(kind=8)                                     :: TSTEP,PSECS,NSECS
    character(len=2)                                 :: klab            ! trj file number label 01, 02, etc, 99
    character(len=1),dimension(natms)                :: atom_n          !first atom of atomname
    
!Electrostatics   
    real(kind=8),parameter                              :: AVSNO = 6.02214199D+23
! 1/(4*PI*EPS0) 8.98750D9 Nm2C-2 = JAngC-2 8.9875D19
    real(kind=8),parameter                              :: PMTTV=8.987551760D19
! electron (C) 
    real(kind=8),parameter                              :: ELCH=1.602177D-019
    real(kind=8)                                        :: enconv_kjmol
    real(kind=8)                                        :: Q1,Q2
    integer                                             :: ncont_th,ncont_Coul_th,ncont_LJ_th    !number of contacts with energy outside threshold window  
    real(kind=8)                                        :: cmxsol_I,cmysol_I,cmzsol_I         !solute centre of mass
    real(kind=8)                                        :: cmxsol_II,cmysol_II,cmzsol_II      !solute centre of mass
    real(kind=8)                                        :: dcmx,dcmy,dcmz,dcmr,dcmr2
        
    integer                                             :: nstart,nend  
    integer                                             :: kstr
    
!xtc2fortran_trj
!    character(len=15),intent(in)            :: trjfile    
    character(len=15)                       :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp, nrc, nrc_C, nrc_LJ, nrc_p
!end xtc2fortran_trj

kstr = 1
!Version 18   if(ncontsamp==1)kstr=0 

!NEW
allocate(res_cont_pair_t0(nb_res,nb_res))

allocate(coul_pot_t0(nb_res,nb_res))
allocate(lj_pot_t0(nb_res,nb_res))
allocate(pe_pot_t0(nb_res,nb_res))

allocate(dcm_trj_US(kntrj))
allocate(dcm_trj_SMD(kntrj,nstep/ncontsamp+kstr))
allocate(drcoord_trj_US(kntrj))
allocate(PCoul_pair_US(nb_res,nb_res,kntrj))
allocate(PCoul_SMD(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(PLJ_pair_US(nb_res,nb_res,kntrj))
allocate(PLJ_SMD(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(pot_pair_US(nb_res,nb_res,kntrj))
allocate(pot_SMD(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))

allocate(sum_Coul_en_t(kntrj),sum_LJ_en_t(kntrj),sum_pot_en_t(kntrj))
allocate(sum_Coul_en_all(kntrj),sum_LJ_en_all(kntrj),sum_pot_en_all(kntrj))

!Umbrella sampling instantaneous values
allocate(sum_Coul_en_US_thr(kntrj,nstep/ncontsamp+kstr)) 
allocate(sum_LJ_en_US_thr(kntrj,nstep/ncontsamp+kstr))
allocate(sum_pot_en_US_thr(kntrj,nstep/ncontsamp+kstr))
allocate(sum_Coul_en_US_all(kntrj,nstep/ncontsamp+kstr))
allocate(sum_LJ_en_US_all(kntrj,nstep/ncontsamp+kstr))
allocate(sum_pot_en_US_all(kntrj,nstep/ncontsamp+kstr))

!Umbrella sampling stdevations   
allocate(stdev_dcm_US(kntrj))
allocate(stdev_PCoul_US(nb_res,nb_res,kntrj))
allocate(stdev_PLJ_US(nb_res,nb_res,kntrj))
allocate(stdev_pot_US(nb_res,nb_res,kntrj))
allocate(stdev_PCoul_US_all(kntrj))
allocate(stdev_PLJ_US_all(kntrj))
allocate(stdev_pot_US_all(kntrj))
allocate(stdev_PCoul_US_thr(kntrj))
allocate(stdev_PLJ_US_thr(kntrj))
allocate(stdev_pot_US_thr(kntrj))

allocate(pot_thr_ftrj(kntrj,nb_res,nb_res,nstep/ncontsamp+kstr))
allocate(pot_thr_ftrj_mean(nb_res,nb_res)) 

!xtc2fortran
allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
!end xtc2fortran

!call SYSTEM('mkdir -p Prot_potent_umbrella_out')
inquire(file='Prot_potent_umbrella_out/log_prot_potent.dat',exist=filex)
if(.not.filex) call SYSTEM('mkdir -p Prot_potent_umbrella_out')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
!n4 = 140
!n5 = 150
!n6 = 160
n7 = 170
n8 = 180
n9 = 190

open(n0,file='Prot_potent_umbrella_out/log_prot_potent.dat',status='unknown',action='write')

!mean energies as a function of the COM-COM distance from umbrella sampling trajectories
open(n1,file='Prot_potent_umbrella_out/umbrella_Coul_COM_thr.dat',status='unknown',action='write')
open(n2,file='Prot_potent_umbrella_out/umbrella_LJ_COM_thr.dat',status='unknown',action='write')
open(n3,file='Prot_potent_umbrella_out/umbrella_POTEn_COM_thr.dat',status='unknown',action='write')
!open(n4,file='Prot_potent_umbrella_out/umbrella_Coul_COM_thr.agr',status='unknown',action='write')
!open(n5,file='Prot_potent_umbrella_out/umbrella_LJ_COM_thr.agr',status='unknown',action='write')
!open(n6,file='Prot_potent_umbrella_out/umbrella_POTEn_COM_thr.agr',status='unknown',action='write')
open(n7,file='Prot_potent_umbrella_out/umbrella_Coul_RC_thr.dat',status='unknown',action='write')
open(n8,file='Prot_potent_umbrella_out/umbrella_LJ_RC_thr.dat',status='unknown',action='write')
open(n9,file='Prot_potent_umbrella_out/umbrella_POTEn_RC_thr.dat',status='unknown',action='write')
            
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3

enconv_kjmol = 1.0d-3*AVSNO                   !conversion factor J ---> kJ/mol

!Version 15
dcm_trj_US = 0.0d0
drcoord_trj_US = 0.0d0
PCoul_pair_US = 0.0d0
PLJ_pair_US = 0.0d0
pot_pair_US = 0.0d0

sum_Coul_en_t = 0.0d0
sum_LJ_en_t = 0.0d0
sum_pot_en_t = 0.0d0

sum_Coul_en_all = 0.0d0
sum_LJ_en_all = 0.0d0
sum_pot_en_all = 0.0d0

sum_Coul_en_US_thr = 0.0d0
sum_LJ_en_US_thr = 0.0d0
sum_pot_en_US_thr = 0.0d0

sum_Coul_en_US_all = 0.0d0
sum_LJ_en_US_all = 0.0d0
sum_pot_en_US_all = 0.0d0

stdev_dcm_US = 0.0d0  
stdev_PCoul_US = 0.0d0
stdev_PLJ_US = 0.0d0
stdev_pot_US = 0.0d0
stdev_PCoul_US_all = 0.0d0
stdev_PLJ_US_all = 0.0d0
stdev_pot_US_all = 0.0d0
stdev_PCoul_US_thr = 0.0d0
stdev_PLJ_US_thr = 0.0d0
stdev_pot_US_thr = 0.0d0

pot_thr_ftrj = 0.00D0
pot_thr_ftrj_mean = 0.0d0

PCoul_SMD = 0.0d0
pot_SMD   = 0.0d0
PLJ_SMD   = 0.0d0

kj1 = 0
kj2 = 0

ncont_th = 0           !Number of contacts with potential energy outside threshold
ncont_Coul_th = 0      !Number of contacts with Coulomb energy outside threshold
ncont_LJ_th = 0        !Number of contacts with vdW energy outside threshold

nrc = 0
nrc_C = 0
nrc_LJ = 0
nrc_p = 0

 cmxsol_I = 0.0d0;  cmysol_I = 0.0d0;  cmzsol_I = 0.0d0
 cmxsol_II = 0.0d0; cmysol_II = 0.0d0; cmzsol_II = 0.0d0

!Version_12
res_cont_pair_t0 = 0
nres_cont_pair_t0 = 0


!*****************************************************************************************************************************************
write(*,*)
write(*,'(10X,A30,10X)')'*****************************'
write(*,'(10X,A30,10X)')'CALLING ROUTINE prot_res_cont'
write(*,'(10X,A30,10X)')'*****************************'
write(*,*)
write(*,'(a72)')'Calculating protein contacts along the first trajectory [md_trj.xtc]...'
write(*,'(a131)')'These will be used in every traj [md_trj01.xtc = md_trj.xtc; md_trj02.xtc; md_trj03.xtc; ... ] to compute mean values and stdev...' 
call f77_molfile_init     ![called in main before potent_map call]
call prot_res_cont(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,&
                   nsoltypes,kumbrella_trj,kMIC,intsf_,rcont_,atomtype,nb_res,chain_label_res,Natresid,name_res,numb_res,nres_cont_pair_t0,res_cont_pair_t0) 
!*****************************************************************************************************************************************


! Start potential map calculation

write(*,*)
write(*,*)'************************************************************'
write(*,*)'*                        rEaD mE!!!                        *'
write(*,*)'* . For number of trajectories > 1 these must be named as: *'
write(*,*)'* md_trj01.xtc, md_trj02.xtc,...,md_trjmn.xtc              *'
write(*,*)'* with mn = kntrj = number of trajectories                 *'
write(*,*)'* . Trajectories must have the exact same time-lenght      *'
write(*,*)'* . The contact map used is found from md_trj.xtc          *'
write(*,*)'* where md_trj.xtc may be one of the md_trjmn.xtc trajs    *'
write(*,*)'* IMPORTANT!!! this must be named md_trj.xtc               *'
write(*,*)'*                                                          *'
write(*,*)'* . Trajectories corresponding to different prot-prot      *'
write(*,*)'* distances may be used e.g. different umbrella sampling   *'
write(*,*)'* windows - contacts found from md_trj.xtc will be used    *'
write(*,*)'* for all trajectories                                     *'
write(*,*)'* No AvErAgEs OvEr TrJs WiLl Be CaLcUlAtEd HeRe            *'
write(*,*)'*                                                          *'
write(*,*)'* ENERGY THRESHOLD FILES:                                  *'
write(*,*)'*                                                          *'
write(*,*)'* residue pairs with intermolecular potential en outside   *'
write(*,*)'* input threshold (kJ/mol) ',Coul_th_min,'-',Coul_th_max,' *'
write(*,*)'*                                                          *'
write(*,*)'* umbrella_Coul_thr.dat:        Coulomb COM-COM profile    *'
write(*,*)'*                                                          *'
write(*,*)'* umbrella_LJ_thr.dat:       van der Waals COM-COM profile *'
write(*,*)'*                                                          *'
write(*,*)'* umbrella_POTEn_thr.dat:      intermol pe COM-COM profile *'
write(*,*)'*                                                          *'
write(*,*)'************************************************************'

write(*,*)
write(*,*)'Starting protein umbrella potential map calculation...'
write(*,*)'Routine pot_map_umbrella'
write(*,*)'Results are printed to Prot_potent_umbrella_out'
write(*,*)'Traj sampling freq for Coul and vdW potential (steps) = ',ncontsamp
write(*,*)'number of umbrella sampling trajectories = ',kntrj
write(*,*)'total number of time-steps = ',nstep
write(*,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(*,*)'Number of contact pairs to be analysed = ',nres_cont_pair_t0
write(*,'(a48,1x,f12.4)')'Fraction of contact pairs to be analysed (%) = ',100.0*dfloat(nres_cont_pair_t0)/dfloat((nb_res/2)**2)
write(*,*)

write(n0,*)
write(n0,*)'Starting protein umbrella potential map calculation...'
write(n0,*)'Routine pot_map_umbrella'
write(n0,*)'Results are printed to Prot_potent_umbrella_out'
write(n0,*)'Traj sampling freq for Coul and vdW potential (steps) = ',ncontsamp
write(n0,*)'number of umbrella sampling trajectories = ',kntrj
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(n0,*)'Number of contact pairs to be analysed = ',nres_cont_pair_t0
write(n0,'(a48,1x,f12.4)')'Fraction of contact pairs to be analysed (%) = ',100.0*dfloat(nres_cont_pair_t0)/dfloat((nb_res/2)**2)
write(n0,*)

write(n1,'(a77)')'# COM-COM (Ang) versus Coul en(kJ/mol) time average; stdev(COM-COM),stdev(en)'   
write(n2,'(a76)')'# COM-COM (Ang) versus vdW en(kJ/mol) time average; stdev(COM-COM),stdev(en)'
write(n3,'(a81)')'# COM-COM (Ang) versus Coul+vdW en(kJ/mol) time average; stdev(COM-COM),stdev(en)'
!Version 19
write(n7,'(a56)')'# RC (nm) versus Coul en(kJ/mol) time average; stdev(en)'   
write(n8,'(a53)')'# RC (nm) versus vdW en(kJ/mol) time average; stdev(en)'
write(n9,'(a60)')'# RC (nm) versus Coul+vdW en(kJ/mol) time average; stdev(en)'



!Loop over number of trajectories available for calculating mean potential map (electrostatic and vdW interactions)
do ktr = 1,kntrj                              !Loop over number of umbrella sampling trajectories
   IF(inputformat.eq.'GROMACS'.and.ktr > 1)call f77_molfile_init                          !This corrected the problem of running this routine after other analyses
   write(*,*)'Starting analysis of trajectory ', ktr,' out of ',kntrj
   write(*,*)

!Label different trajectories input files: md_trj01.xtc, md_trj02.xtc, etc
   IF(kntrj>1)THEN
!   if(ktr==1)trjfile='md_trj01.xtc'
!   if(ktr==2)trjfile='md_trj02.xtc'
!   etc

      write (klab,'(I2.2)') ktr                             ! converting integer (ktr) to string (klab) using an 'internal file'
      trjfile='md_trj'//trim(klab)//'.xtc'
!check      write(*,*) trjfile 
   ENDIF

!time counter of electrostatic interactions (1,2,3,...) with 1,2,3 separated by ncontsamp
   it = 0
   jt = 0  
   kt = 0
   icm = 0
   
!xtc2fortran_trj
   
   IF(inputformat.eq.'GROMACS')THEN
      infile = TRIM(trjfile)
      intype = 'auto'
      npart  = -1
      handle(1) = -1
      handle(2) = -1
      handle(3) = -1
      handle(4) = -1
!set up everything and register all static plugins
!NG      call f77_molfile_init

      call f77_molfile_open_read(handle(1),npart,infile,intype)

      if (handle(1).lt.0) then
         print*,'file type unknown or not registered'
         stop
      else
         print*,'file successfully opened:'
         print*,'handle:',handle(1)
         print*,'npart: ',npart
         print*,'nsteps:',nstep
      end if
        
      if(npart /= natms)then
        write(*,*)'error - number of atoms is wrong in solv_inp.dat'
        stop
      endif  
   ENDIF    
   
!end xtc2fortran_trj  

!Start trajectory analysis
   
!Version 18   do j = 1,nstep                                               !Loop over number of steps of trajectory
   do j = 0,nstep   
      if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/    
      IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
         if(nens==1)then
            read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
         elseif(nens==0)then
            cell(1) = cube(1)
            cell(2) = cube(2)
            cell(3) = cube(3)
         endif   
!read solute and water coordinates
         read(ninput,*)natcheck
         if(natcheck.ne.natms)then
            write(*,*)'error - number of atoms from solv_inp is wrong'
            stop
         endif   
         do i = 1,natms
            read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!NG         if(j==1)atom_symb(i)=TRIM(atomname(i))
!check         if(j==1.and.atom_symb(i)/='H')write(*,*)atom_symb(i)
         end do
      ELSEIF(inputformat.eq.'GROMACS')THEN
      
!xtc2fortran_trj 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
       status = 1   ! status=1 on entry means read
       call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
       cell(1)=dble(box(1))
       cell(2)=dble(box(2))
       cell(3)=dble(box(3))
       kp=1
!       write(*,*)cell
       do i = 1,npart*3,3
          x(kp)=xyz(i)
          y(kp)=xyz(i+1)
          z(kp)=xyz(i+2)
          kp = kp + 1
       end do   
   
!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
      ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 


!Calculate proteins centre of mass (COM) distance along time for each trajectory
!Version 18      if(j==1.or.mod(j,ncontsamp)==0)then 
      if(mod(j,ncontsamp)==0)then       
!check         write(*,*)'calculating com distance...trj ',ktr,' step ',j
         icm = icm + 1
         kt = kt + 1                                                                              !pseudo-time counter
         nstart = 1
         nend   = natmsol/2
         call sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol_I,cmysol_I,cmzsol_I)
         nstart = natmsol/2 + 1
         nend   = natmsol
         call sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol_II,cmysol_II,cmzsol_II)
         dcmx = cmxsol_I - cmxsol_II
         dcmy = cmysol_I - cmysol_II
         dcmz = cmzsol_I - cmzsol_II
         if(kMIC==1)then
            dcmx = dcmx - dnint(dcmx/cell(1))*cell(1)
            dcmy = dcmy - dnint(dcmy/cell(2))*cell(2)
            dcmz = dcmz - dnint(dcmz/cell(3))*cell(3)
         endif
         dcmr2 = dcmx**2 + dcmy**2 + dcmz**2
         dcmr  = dsqrt(dcmr2)
         dcm_trj_US(ktr) = dcm_trj_US(ktr) + dcmr          !accumulate com distance for each trajectory         
         dcm_trj_SMD(ktr,kt) = dcmr                        !get com distance at each time for each trajectory - for stdev calculation
!Version19
         drcoord_trj_US(ktr) = drcoord_trj_US(ktr) + dsqrt(dcmx**2.0)           !accumulate com distance for each trajectory



!check for a single trajectory         
!              PSECS = TSTEP*dfloat(j-1)                            !psecs
!check         PSECS = TSTEP*dfloat(icm-1)*dfloat(ncontsamp)
!check         NSECS = 1e-3*PSECS                                   !nsecs
!check         if(ktr==1)write(n1,'(e19.3,3x,f12.4)')NSECS,dcmr                     !print check for a single trajectory - to compare with the mean over trajectories 
!NG         if(ktr==1)write(n1,'(e19.3,3x,f12.4,2x,a6,1x,i6)')NSECS,dcmr,'trj = ',ktr  
      endif 

!End COM calculation


!Calculate electrostatic interactions of residue contact pairs

!Version 18      if(j==1.or.mod(j,ncontsamp)==0)then
      if(mod(j,ncontsamp)==0)then   
         it = it +1                           !pseudo-time counter
         k1 = 0                               !Chains ABCD atomic index
         do ic1=1,nb_res/2                    !residues of ABCD
            do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
               k1 = k1 + 1
               k2 = natmsol/2                 !Chains EFGH atomic index                  
               do ic2=nb_res/2+1,nb_res       !residues of EFGH
                  do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                     k2 = k2 + 1
                     if(res_cont_pair_t0(ic1,ic2)==1)then                  !calculate electrostatic potential of the pair
                        dx = x(k1)-x(k2)
                        dy = y(k1)-y(k2)
                        dz = z(k1)-z(k2)
                        if(kMIC==1)then
                           dx = dx - dnint(dx/cell(1))*cell(1)
                           dy = dy - dnint(dy/cell(2))*cell(2)
                           dz = dz - dnint(dz/cell(3))*cell(3)
                        endif 
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        Q1 = chrg(k1)*ELCH
                        Q2 = chrg(k2)*ELCH
!Warning!!! en is summed for [nstep/ncontsamp] +1 steps                        
                        PCoul_pair_US(ic1,ic2,ktr) = PCoul_pair_US(ic1,ic2,ktr) + enconv_kjmol*PMTTV*Q1*Q2/dr                  !accumulate en for each trj
                        pot_pair_US(ic1,ic2,ktr)   = pot_pair_US(ic1,ic2,ktr) + enconv_kjmol*PMTTV*Q1*Q2/dr                    !accumulate en for each trj

                        PCoul_SMD(ktr,ic1,ic2,it) = PCoul_SMD(ktr,ic1,ic2,it) + enconv_kjmol*PMTTV*Q1*Q2/dr
                        pot_SMD(ktr,ic1,ic2,it)   = pot_SMD(ktr,ic1,ic2,it)   + enconv_kjmol*PMTTV*Q1*Q2/dr
                        
                        pot_thr_ftrj(ktr,ic1,ic2,it) = pot_thr_ftrj(ktr,ic1,ic2,it) + enconv_kjmol*PMTTV*Q1*Q2/dr

                     endif
                  end do
               end do   
            end do
         end do
      endif                          !end electrostatics sampling   
   
!End electrostatics   


!Calculate van der Waals interactions of residue contact pairs
!Version 18      if(j==1.or.mod(j,ncontsamp)==0)then
      if(mod(j,ncontsamp)==0)then   
         jt = jt +1                           !pseudo-time counter
         k1 = 0                               !Chains ABCD atomic index
         do ic1=1,nb_res/2                    !residues of ABCD
            do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
               k1 = k1 + 1
               k2 = natmsol/2                 !Chains EFGH atomic index                  
               do ic2=nb_res/2+1,nb_res       !residues of EFGH
                  do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                     k2 = k2 + 1
                     if(res_cont_pair_t0(ic1,ic2)==1)then                  !calculate Lennard-Jones potential of the pair
                        dx = x(k1)-x(k2)
                        dy = y(k1)-y(k2)
                        dz = z(k1)-z(k2)
                        if(kMIC==1)then
                           dx = dx - dnint(dx/cell(1))*cell(1)
                           dy = dy - dnint(dy/cell(2))*cell(2)
                           dz = dz - dnint(dz/cell(3))*cell(3)
                        endif
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
!Apply geometric rule - units [kJ/mol nm**6]  and [kJ/mol nm**12] - convert to [kJ/mol Ang**6]  and [kJ/mol Ang**12]                   
!kJ/mol Ang**6                        C6  = 1.0d6*dsqrt(C6kjmol(k1)*C6kjmol(k2))  
!kJ/mol Ang**12                       C12 = 1.0d12*dsqrt(C12kjmol(k1)*C12kjmol(k2)) 
!r Ang                        PLJ_pair(ic1,ic2,jt) = PLJ_pair(ic1,ic2,jt) + (C12/dr**12.0 - C6/dr**6.0)                      !accumulate over different trajectories
!r Ang                        pot_pair(ic1,ic2,jt) = pot_pair(ic1,ic2,jt) + (C12/dr**12.0 - C6/dr**6.0)                      !accumulate over different trajectories
                        C6  = dsqrt(C6kjmol(k1)*C6kjmol(k2))                    !kJ/mol nm**6  
                        C12 = dsqrt(C12kjmol(k1)*C12kjmol(k2))                  !kJ/mol nm**12
!For GROMOS force field C12 was modified according to GROMOS interactions rules for C12                        

                        PLJ_pair_US(ic1,ic2,ktr) = PLJ_pair_US(ic1,ic2,ktr) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)   !r nm   !accumulate en for each trj
                        pot_pair_US(ic1,ic2,ktr) = pot_pair_US(ic1,ic2,ktr) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)  !r nm    !accumulate en for each trj
                        
                        PLJ_SMD(ktr,ic1,ic2,jt) = PLJ_SMD(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)
                        pot_SMD(ktr,ic1,ic2,jt) = pot_SMD(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)     
                           
                        pot_thr_ftrj(ktr,ic1,ic2,jt) = pot_thr_ftrj(ktr,ic1,ic2,jt) + (C12_IJ(k1,k2)/(0.1d0*dr)**12.0 - C6_IJ(k1,k2)/(0.1d0*dr)**6.0)
                        
                     endif
                  end do
               end do   
            end do
         end do
      endif                          !end vdW sampling   
   
!End van der Waals   

     
!xtc2fortran_trj
      if(inputformat.eq.'GROMACS'.and.status.eq.0) then
         write(*,*)'error on reading trajectory file - step',i
         stop
      endif
!end xtc2fortran_trj    
    
   end do                !end time-step main loop

!xtc2fortran_trj
   IF(inputformat.eq.'GROMACS')THEN
      call f77_molfile_finish
      write(*,*)'finished reading xtc trj...'
   ENDIF   
!end xtc2fortran_trj
   
end do                   !end trajectories loop 

!=======================================================================================================
!VERSION 16
!Find mean energy to probe threshold - first trj
!A mean value is obtained for each residue pair...

do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then
! +1
         do it=1,nstep/ncontsamp + kstr
            pot_thr_ftrj_mean(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2) + pot_thr_ftrj(1,ic1,ic2,it)
         end do
! +1         
!Version 18         pot_thr_ftrj_mean(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2)/(dfloat(nstep)/dfloat(ncontsamp)+dfloat(kstr))
         pot_thr_ftrj_mean(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2)/dfloat(nstep/ncontsamp + kstr)
      endif
   end do
end do
!=======================================================================================================


!Print electrostatic, van der Waals, and total intermolecular potential energy, COM-COM profiles 

do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then                  !Analyse only res pairs within input distance cut-off
         nrc = nrc + 1        
          
!=======================================================================================================instantaneous en
!VERSION 16
!VERSION 16         coul_pot_t0(ic1,ic2) = PCoul_SMD(1,ic1,ic2,1)                 !value for umbella sampling trj 1 at time zero [NOT USED ANYMORE]
!VERSION 16         lj_pot_t0(ic1,ic2)   = PLJ_SMD(1,ic1,ic2,1)                   !value for umbella sampling trj 1 at time zero [NOT USED ANYMORE]
!VERSION 16         pe_pot_t0(ic1,ic2)   = pot_SMD(1,ic1,ic2,1)                   !value for umbella sampling trj 1 at time zero
!VERSION 16         
!VERSION 16!Version 14 - account for res-res interactions above or below threshold at times > t0
!VERSION 16! The largest value is used to choose residues to print - if the interaction between any two residues is above or below a threshold 
!VERSION 16! at any time during the trj this will be printed - the res-res for trj 1 are used for all umbrella sampling trjs
!VERSION 16
!VERSION 16!NG check         do it=2,nstep/ncontsamp + 1
!VERSION 16         do it=2,nstep/ncontsamp + 1
!VERSION 16            coul_pot_t = PCoul_SMD(1,ic1,ic2,it)
!VERSION 16            if((coul_pot_t>0.0d0).and.(coul_pot_t>coul_pot_t0(ic1,ic2)))then                    !REPULSION
!VERSION 16               coul_pot_t0(ic1,ic2) = coul_pot_t
!VERSION 16            elseif((coul_pot_t<0.0d0).and.(coul_pot_t<coul_pot_t0(ic1,ic2)))then                !ATTRACTION
!VERSION 16               coul_pot_t0(ic1,ic2) = coul_pot_t
!VERSION 16            endif
!VERSION 16            
!VERSION 16            lj_pot_t = PLJ_SMD(1,ic1,ic2,it) 
!VERSION 16            if((lj_pot_t>0.0d0).and.(lj_pot_t>lj_pot_t0(ic1,ic2)))then
!VERSION 16               lj_pot_t0(ic1,ic2) = lj_pot_t
!VERSION 16            elseif((lj_pot_t<0.0d0).and.(lj_pot_t<lj_pot_t0(ic1,ic2)))then
!VERSION 16               lj_pot_t0(ic1,ic2) = lj_pot_t   
!VERSION 16            endif
!VERSION 16            
!VERSION 16            pe_pot_t = pot_SMD(1,ic1,ic2,it)
!VERSION 16            if((pe_pot_t>0.0d0).and.(pe_pot_t>pe_pot_t0(ic1,ic2)))then
!VERSION 16               pe_pot_t0(ic1,ic2) = pe_pot_t
!VERSION 16            elseif((pe_pot_t<0.0d0).and.(pe_pot_t<pe_pot_t0(ic1,ic2)))then   
!VERSION 16               pe_pot_t0(ic1,ic2) = pe_pot_t
!VERSION 16            endif
!VERSION 16         end do
        
!Version 16 - use the mean en instead of the instantaneous energies

         pe_pot_t0(ic1,ic2) = pot_thr_ftrj_mean(ic1,ic2)
!=======================================================================================================end instantaneous en         
         
       
!Number of contacts with energies within the threshold window   [Coul_th_min,Coul_th_max]

!Version 15         if((coul_pot_t0(ic1,ic2) > Coul_th_min).and.(coul_pot_t0(ic1,ic2) < Coul_th_max))then
         if((pe_pot_t0(ic1,ic2) > Coul_th_min).and.(pe_pot_t0(ic1,ic2) < Coul_th_max))then
            ncont_Coul_th = ncont_Coul_th + 1
!Version 15         if((lj_pot_t0(ic1,ic2) > LJ126_th_min).and.(lj_pot_t0(ic1,ic2) < LJ126_th_max))then
            ncont_LJ_th = ncont_LJ_th + 1
!Version 15         if((pe_pot_t0(ic1,ic2) > Coul_th_min).and.(pe_pot_t0(ic1,ic2) < Coul_th_max))then
            ncont_th = ncont_th + 1
         endif
 
!calculate the energy of all contacts with an interaction energy within the threshold window for each trj [Coul_th_min,Coul_th_max]
         do ktr = 1,kntrj                                            !LOOP OVER TRJs            
!Version 15            if((coul_pot_t0(ic1,ic2) > Coul_th_min).and.(coul_pot_t0(ic1,ic2) < Coul_th_max))then
            if((pe_pot_t0(ic1,ic2) > Coul_th_min).and.(pe_pot_t0(ic1,ic2) < Coul_th_max))then
               sum_Coul_en_t(ktr) = sum_Coul_en_t(ktr) + PCoul_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
!Version 15            if((lj_pot_t0(ic1,ic2) > LJ126_th_min).and.(lj_pot_t0(ic1,ic2) < LJ126_th_max))then
               sum_LJ_en_t(ktr) = sum_LJ_en_t(ktr) + PLJ_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
               sum_pot_en_t(ktr) = sum_pot_en_t(ktr) + pot_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
!               
!calculate energy at each time-step for stdev calculation
! +1
               do kt = 1,nstep/ncontsamp + kstr
                  sum_Coul_en_US_thr(ktr,kt) = sum_Coul_en_US_thr(ktr,kt) + PCoul_SMD(ktr,ic1,ic2,kt)
                  sum_LJ_en_US_thr(ktr,kt) = sum_LJ_en_US_thr(ktr,kt) + PLJ_SMD(ktr,ic1,ic2,kt)
                  sum_pot_en_US_thr(ktr,kt) = sum_pot_en_US_thr(ktr,kt) + pot_SMD(ktr,ic1,ic2,kt)
               end do   
            endif
            
!calculate the sum of the potential energy of all contacts
            sum_Coul_en_all(ktr) = sum_Coul_en_all(ktr) + PCoul_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
            sum_LJ_en_all(ktr)   = sum_LJ_en_all(ktr)   + PLJ_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
            sum_pot_en_all(ktr)  = sum_pot_en_all(ktr)  + pot_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
!            
!calculate energy at each time-step for stdev calculation    
!WARNING!!! kt should change between 1 and nstep/ncontsamp + 1
! + 1
            do kt = 1,nstep/ncontsamp + kstr
               sum_Coul_en_US_all(ktr,kt) = sum_Coul_en_US_all(ktr,kt) + PCoul_SMD(ktr,ic1,ic2,kt)
               sum_LJ_en_US_all(ktr,kt) = sum_LJ_en_US_all(ktr,kt) + PLJ_SMD(ktr,ic1,ic2,kt)
               sum_pot_en_US_all(ktr,kt) = sum_pot_en_US_all(ktr,kt) + pot_SMD(ktr,ic1,ic2,kt)
            end do
            
         end do           !End loop over trjs
         
      endif
   end do                            !End loop over residues of EFGH
end do                               !End loop over residues of ABCD

!Calculate standard deviations over time for each trajectory  
do ktr = 1,kntrj  
! +1
   do kt = 1,nstep/ncontsamp + kstr 
      stdev_dcm_US(ktr)           = stdev_dcm_US(ktr) + (dcm_trj_SMD(ktr,kt) - dcm_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr))**2.0            
      
      stdev_PCoul_US_thr(ktr)     = stdev_PCoul_US_thr(ktr) + (sum_Coul_en_US_thr(ktr,kt) - sum_Coul_en_t(ktr))**2.0
      stdev_PLJ_US_thr(ktr)       = stdev_PLJ_US_thr(ktr) + (sum_LJ_en_US_thr(ktr,kt) - sum_LJ_en_t(ktr))**2.0
      stdev_pot_US_thr(ktr)       = stdev_pot_US_thr(ktr) + (sum_pot_en_US_thr(ktr,kt) - sum_pot_en_t(ktr))**2.0       
      
      stdev_PCoul_US_all(ktr)     = stdev_PCoul_US_all(ktr) + (sum_Coul_en_US_all(ktr,kt) - sum_Coul_en_all(ktr))**2.0
      stdev_PLJ_US_all(ktr)       = stdev_PLJ_US_all(ktr) + (sum_LJ_en_US_all(ktr,kt) - sum_LJ_en_all(ktr))**2.0
      stdev_pot_US_all(ktr)       = stdev_pot_US_all(ktr) + (sum_pot_en_US_all(ktr,kt) - sum_pot_en_all(ktr))**2.0
   end do
!Version 18   
   stdev_dcm_US(ktr)           = dsqrt(stdev_dcm_US(ktr)/dfloat(nstep/ncontsamp))
   stdev_PCoul_US_all(ktr)     = dsqrt(stdev_PCoul_US_all(ktr)/dfloat(nstep/ncontsamp))
   stdev_PLJ_US_all(ktr)       = dsqrt(stdev_PLJ_US_all(ktr)/dfloat(nstep/ncontsamp))
   stdev_pot_US_all(ktr)       = dsqrt(stdev_pot_US_all(ktr)/dfloat(nstep/ncontsamp))
   stdev_PCoul_US_thr(ktr)     = dsqrt(stdev_PCoul_US_thr(ktr)/dfloat(nstep/ncontsamp))
   stdev_PLJ_US_thr(ktr)       = dsqrt(stdev_PLJ_US_thr(ktr)/dfloat(nstep/ncontsamp))
   stdev_pot_US_thr(ktr)       = dsqrt(stdev_pot_US_thr(ktr)/dfloat(nstep/ncontsamp))
end do   

do ktr = 1,kntrj
   do ic1=1,nb_res/2                                           !residues of ABCD
      do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
         if(res_cont_pair_t0(ic1,ic2)==1)then       
            if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then
!+1            
               do kt = 1,nstep/ncontsamp + kstr          
                  stdev_PCoul_US(ic1,ic2,ktr) = stdev_PCoul_US(ic1,ic2,ktr) + (PCoul_SMD(ktr,ic1,ic2,kt) - PCoul_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr))**2.0 
                  stdev_PLJ_US(ic1,ic2,ktr)   = stdev_PLJ_US(ic1,ic2,ktr) + (PLJ_SMD(ktr,ic1,ic2,kt) - PLJ_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr))**2.0 
                  stdev_pot_US(ic1,ic2,ktr)   = stdev_pot_US(ic1,ic2,ktr) + (pot_SMD(ktr,ic1,ic2,kt) - pot_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr))**2.0 
               end do
            endif
!Version 18            
            stdev_PCoul_US(ic1,ic2,ktr) = dsqrt(stdev_PCoul_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp))
            stdev_PLJ_US(ic1,ic2,ktr)   = dsqrt(stdev_PLJ_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp))
            stdev_pot_US(ic1,ic2,ktr)   = dsqrt(stdev_pot_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp))
         endif
      end do
   end do
end do   !End loop over trjs
  
!Print 
do ic1=1,nb_res/2                                           !residues of ABCD
   do ic2=nb_res/2+1,nb_res                                 !residues of EFGH
      if(res_cont_pair_t0(ic1,ic2)==1)then 
!PRINT HEADERS OF THRESHOLD FILES
!Version 15 - print the same pairs for Coul. LJ and potential - use a single threshold en
!Coulomb - threshold header
         if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then
         
            nrc_C = nrc_C + 1
            write(n1,*)
            write(n1,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_C 
            write(n7,*)
            write(n7,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_C
!van der Waals - threshold header       
!Version 15        if((lj_pot_t0(ic1,ic2) < LJ126_th_min).or.(lj_pot_t0(ic1,ic2) > LJ126_th_max))then
            nrc_LJ = nrc_LJ + 1
            write(n2,*)
            write(n2,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_LJ
            write(n8,*)
            write(n8,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_LJ
!Intermolecular potential - sum of Coulomb and van der Waals - threshold header         
!         if((pe_pot_t0(ic1,ic2) < Coul_th_min).or.(pe_pot_t0(ic1,ic2) > Coul_th_max))then
            nrc_p = nrc_p + 1
            write(n3,*)
            write(n3,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_p
             write(n9,*)
            write(n9,'(a1,1x,2(i6,1x,a4,1x,i6,a2,5x),3x,a9,i6)') &
            '#',ic1,name_res(ic1),numb_res(ic1),chain_label_res(ic1),ic2,name_res(ic2),numb_res(ic2),chain_label_res(ic2),'c-pair = ',nrc_p
         
            do ktr = 1,kntrj
               dcm_mean = dcm_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr)
               drcoord_mean = drcoord_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr)
               coul_pot = PCoul_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
               lj_pot   = PLJ_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
               pe_pot   = pot_pair_US(ic1,ic2,ktr)/dfloat(nstep/ncontsamp + kstr)
               write(n1,'(4(f12.4,3x))')dcm_mean,coul_pot,stdev_dcm_US(ktr),stdev_PCoul_US(ic1,ic2,ktr)
               write(n2,'(4(f12.4,3x))')dcm_mean,lj_pot,stdev_dcm_US(ktr),stdev_PLJ_US(ic1,ic2,ktr) 
               write(n3,'(4(f12.4,3x))')dcm_mean,pe_pot,stdev_dcm_US(ktr),stdev_pot_US(ic1,ic2,ktr)
!Version 19 - Reaction coordinate - COM-COM-x               
               write(n7,'(3(f12.4,3x))')drcoord_mean*0.1d0,coul_pot,stdev_PCoul_US(ic1,ic2,ktr)                !nm
               write(n8,'(3(f12.4,3x))')drcoord_mean*0.1d0,lj_pot,stdev_PLJ_US(ic1,ic2,ktr)                    !nm
               write(n9,'(3(f12.4,3x))')drcoord_mean*0.1d0,pe_pot,stdev_pot_US(ic1,ic2,ktr)                    !nm
            end do
            
         endif
      endif
   end do
end do
write(n1,*)
write(n2,*)
write(n3,*)
!Version 19
write(n7,*)
write(n8,*)
write(n9,*)

!Print mean energies of contacts with energy above and below threshold [Coul_th_min - Coul_th]
write(n1,*)
write(n1,'(a59,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean Coul energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n1,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n1,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_Coul_th
write(n1,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_Coul_th)/dfloat(nres_cont_pair_t0)*100.0

write(n2,*)
write(n2,'(a57,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean LJ energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n2,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n2,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_LJ_th
write(n2,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_LJ_th)/dfloat(nres_cont_pair_t0)*100.0

write(n3,*)
write(n3,'(a58,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean pot energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n3,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n3,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_th
write(n3,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_th)/dfloat(nres_cont_pair_t0)*100.0

write(n7,*)
write(n7,'(a59,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean Coul energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n7,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n7,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_Coul_th
write(n7,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_Coul_th)/dfloat(nres_cont_pair_t0)*100.0

write(n8,*)
write(n8,'(a57,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean LJ energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n8,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n8,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_LJ_th
write(n8,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_LJ_th)/dfloat(nres_cont_pair_t0)*100.0

write(n9,*)
write(n9,'(a58,1x,a1,f7.1,1x,a1,f7.1,a1)')'# Mean pot energy of contacts with pot energy in the range','[',Coul_th_min,'-',Coul_th_max,']'
write(n9,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
write(n9,'(a48,i9)')'# Number of contacts outside threshold window = ',ncont_th
write(n9,'(a48,f9.3)')'# Fraction of contacts outside threshold window (%) = ',dfloat(ncont_th)/dfloat(nres_cont_pair_t0)*100.0

do ktr = 1,kntrj
   dcm_mean = dcm_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr)
   drcoord_mean = drcoord_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr)
   write(n1,'(4(f12.4,3x))')dcm_mean,sum_Coul_en_t(ktr),stdev_dcm_US(ktr),stdev_PCoul_US_thr(ktr)
   write(n2,'(4(f12.4,3x))')dcm_mean,sum_LJ_en_t(ktr),stdev_dcm_US(ktr),stdev_PLJ_US_thr(ktr)
   write(n3,'(4(f12.4,3x))')dcm_mean,sum_pot_en_t(ktr),stdev_dcm_US(ktr),stdev_pot_US_thr(ktr)
!Version 19   
   write(n7,'(3(f12.4,3x))')drcoord_mean*0.1d0,sum_Coul_en_t(ktr),stdev_PCoul_US_thr(ktr)                  !nm
   write(n8,'(3(f12.4,3x))')drcoord_mean*0.1d0,sum_LJ_en_t(ktr),stdev_PLJ_US_thr(ktr)                      !nm
   write(n9,'(3(f12.4,3x))')drcoord_mean*0.1d0,sum_pot_en_t(ktr),stdev_pot_US_thr(ktr)                     !nm
end do
write(n1,*)
write(n2,*)
write(n3,*)
!Version 19
write(n7,*)
write(n8,*)
write(n9,*)
   
!Print mean (over number of trajectories) energies of all contacts 
write(n1,*)
write(n1,'(a34)')'# Mean Coul energy of all contacts'
write(n1,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n2,*)
write(n2,'(a32)')'# Mean LJ energy of all contacts' 
write(n2,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n3,*)
write(n3,'(a33)')'# Mean pot energy of all contacts'
write(n3,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0
!Version 19
write(n7,*)
write(n7,'(a34)')'# Mean Coul energy of all contacts'
write(n7,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n8,*)
write(n8,'(a32)')'# Mean LJ energy of all contacts' 
write(n8,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

write(n9,*)
write(n9,'(a33)')'# Mean pot energy of all contacts'
write(n9,'(a34,i9)')'# Total nb of contacts analysed = ',nres_cont_pair_t0

do ktr = 1,kntrj
   dcm_mean = dcm_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr) 
   drcoord_mean = drcoord_trj_US(ktr)/dfloat(nstep/ncontsamp + kstr)
   write(n1,'(4(f12.4,3x))')dcm_mean,sum_Coul_en_all(ktr),stdev_dcm_US(ktr),stdev_PCoul_US_all(ktr)
   write(n2,'(4(f12.4,3x))')dcm_mean,sum_LJ_en_all(ktr),stdev_dcm_US(ktr),stdev_PLJ_US_all(ktr)
   write(n3,'(4(f12.4,3x))')dcm_mean,sum_pot_en_all(ktr),stdev_dcm_US(ktr),stdev_pot_US_all(ktr)
!Version 19
   write(n7,'(3(f12.4,3x))')drcoord_mean*0.1d0,sum_Coul_en_all(ktr),stdev_PCoul_US_all(ktr)               !nm
   write(n8,'(3(f12.4,3x))')drcoord_mean*0.1d0,sum_LJ_en_all(ktr),stdev_PLJ_US_all(ktr)                   !nm
   write(n9,'(3(f12.4,3x))')drcoord_mean*0.1d0,sum_pot_en_all(ktr),stdev_pot_US_all(ktr)                  !nm
end do

deallocate(dcm_trj_US)
deallocate(drcoord_trj_US)

deallocate(coul_pot_t0)
deallocate(lj_pot_t0)
deallocate(pe_pot_t0)

deallocate(dcm_trj_US)
deallocate(dcm_trj_SMD)
deallocate(PCoul_pair_US)
deallocate(PCoul_SMD)
deallocate(PLJ_pair_US)
deallocate(PLJ_SMD)
deallocate(pot_pair_US)
deallocate(pot_SMD)
   
deallocate(sum_Coul_en_t,sum_LJ_en_t,sum_pot_en_t)
deallocate(sum_Coul_en_all,sum_LJ_en_all,sum_pot_en_all)   

deallocate(stdev_dcm_US,stdev_PCoul_US,stdev_PLJ_US,stdev_pot_US)
deallocate(stdev_PCoul_US_all,stdev_PLJ_US_all,stdev_pot_US_all)
deallocate(stdev_PCoul_US_thr,stdev_PLJ_US_thr,stdev_pot_US_thr)

deallocate(pot_thr_ftrj)
deallocate(pot_thr_ftrj_mean)

deallocate(xyz)
deallocate(x,y,z)

return

END SUBROUTINE pot_map_umbrella



!Version 22
SUBROUTINE prot_water_umbrella(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,Zattype,mattype,ZMBIO,&
                               nsoltypes,nions,mcontsamp,atomtype,chrg,nb_res,chain_label_res,Natresid,name_res,numb_res,kMIC_,C6kjmol,C12kjmol,kprot_wat,kwat_wat,&
                               C6_IW,C12_IW)                      

!   Computes the mean protein-water and water-water potential energy - no Ewald sum

    integer,intent(in)                                    :: ninput,nbox
    character(len=7),intent(in)                           :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                               :: dt                   ! time between frames (fs)
    integer,intent(in)                                    :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                               :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                                    :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
    character(len=4),dimension(natms)                     :: atomname    
    integer,intent(in)                                    :: nsoltypes            ! Number of solute distinct atomic species (C, H, O, N etc)
    integer,intent(in)                                    :: mcontsamp            ! Protein contact map electrostatics sampling frequency      
    character(len=4),dimension(natms),intent(in)          :: atomtype 
    real(kind=4),dimension(natms),intent(in)              :: chrg     
    integer,intent(in)                                    :: nb_res                 ! Number of residues [GROMCAS]
    character(len=1),dimension(nb_res),intent(in)         :: chain_label_res        ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nb_res),intent(in)                  :: Natresid               ! Number of atoms of each residue
    character(len=4),dimension(nb_res),intent(in)         :: name_res
    integer,dimension(nb_res),intent(in)                  :: numb_res               ! Number of the residue: from 1 to the total number of residues of each chain
    real(kind=8),dimension(natmsol),intent(in)            :: C6kjmol,C12kjmol       ! Solute Lennard-Jones parameters
    real(kind=8),dimension(natmsol),intent(in)            :: C6_IW,C12_IW           ! Solute Lennard-Jones C6 and C12 crossed parameters kj/mol nm**6; kj/mol nm**12

    integer,intent(in)                                    :: kMIC_                   ! apply Minimum Image Convention or not to proteins

    character(len=2),dimension(nsoltypes),intent(in)      :: Zattype                ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)          :: mattype                ! Mass of atoms of type (C, H, O, N etc)
    real(kind=8),dimension(natmsol),intent(in)            :: ZMBIO
    integer,intent(in)                                    :: kprot_wat,kwat_wat
    
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                 !atomic index (solute and solvent index)
    integer      :: n0,n1,n2,n3,n4
    integer      :: n7,n8,n9,n10
    integer      :: i, j, k
    integer      :: ic1,ic2,iats1,iats2,k1,k2,it,jt,kt,kwt,it1,it2,jt1,jt2
    integer      :: ktr,kntrj
    integer      :: icm
        
    real(kind=8)                                     :: dx,dy,dz,dr,dr2,rhi,rhf
    integer                                          :: natcheck
    logical                                          :: filex
!Version 15   
    real(kind=8)                                     :: dcm_trj_US,drcoord_trj_US
    real(kind=8)                                     :: drcoord_mean,dcm_mean      
    real(kind=8)                                     :: C6,C12   
    real(kind=8)                                     :: TSTEP,PSECS,NSECS
    character(len=2)                                 :: klab            ! trj file number label 01, 02, etc, 99
    character(len=1),dimension(natms)                :: atom_n          !first atom of atomname
    
!Electrostatics   
    real(kind=8),parameter                              :: AVSNO = 6.02214199D+23
! 1/(4*PI*EPS0) 8.98750D9 Nm2C-2 = JAngC-2 8.9875D19
    real(kind=8),parameter                              :: PMTTV=8.987551760D19
! electron (C) 
    real(kind=8),parameter                              :: ELCH=1.602177D-019
    real(kind=8)                                        :: enconv_kjmol
    real(kind=8)                                        :: Q1,Q2,QW,QW1,QW2 
    real(kind=8)                                        :: cmxsol_I,cmysol_I,cmzsol_I         !solute centre of mass
    real(kind=8)                                        :: cmxsol_II,cmysol_II,cmzsol_II      !solute centre of mass
    real(kind=8)                                        :: dcmx,dcmy,dcmz,dcmr,dcmr2       
    integer                                             :: nstart,nend  
    integer                                             :: kstr
    integer                                             :: nwm,kw,ic12                             !water model number choice
    integer                                             :: kw1,kw2,kw_
    real(kind=4)                                        :: chrg_OW,chrg_HW
    real(kind=8)                                        :: C6_OW,C12_OW
    real(kind=8)                                        :: R_Coul_cut
    character(len=4)                                    :: atom_WO,atom_WH 
    real(kind=8),dimension(:,:),allocatable             :: PCoul_RES_W                        !Coul Prot-water
    real(kind=8),dimension(:,:),allocatable             :: PLJ_RES_W                          !van der Waals Prot-water
    real(kind=8),dimension(:),allocatable               :: Pot_Prot_wat                       !Coul + van der Waals
    real(kind=8),dimension(:),allocatable               :: pot_coul_PW_m
    real(kind=8),dimension(:),allocatable               :: pot_LJ_PW_m       
    real(kind=8)                                        :: pot_prot_W_mean
    real(kind=8)                                        :: pot_coul_PW_mean
    real(kind=8)                                        :: pot_LJ_PW_mean      
    real(kind=8)                                        :: stdev_prot_W  
    real(kind=8)                                        :: stdev_coul_PW 
    real(kind=8)                                        :: stdev_LJ_PW   
    real(kind=8)                                        :: stdev_PW_pot,stdev_PW_coul,stdev_PW_LJ
    Logical                                             :: LWatWat,LProtWat
    real(kind=8),dimension(:),allocatable               :: pot_WW_coul,pot_WW_vdw,EP_WW
    real(kind=8)                                        :: pot_W_W_mean,pot_WW_coul_mean,pot_WW_vdw_mean 
    real(kind=8)                                        :: stdev_W_W,stdev_WW_pot    
    real(kind=8)                                        :: stdev_W_W_Coul,stdev_W_W_LJ,stdev_WW_Coul,stdev_WW_LJ   

!xtc2fortran_trj
!    character(len=15),intent(in)            :: trjfile    
    character(len=15)                       :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp 

!end xtc2fortran_trj

!Water models
!SPC/E nwm = 1
!TIP4P/EW nwm = 2     !not coded
nwm = 1
if(nwm ==1)then
   chrg_OW = -0.8476e0
   chrg_HW = 0.4238e0
   C6_OW   = 0.0026173456d0                    !For water-water only
   C12_OW  = 2.634129d-6                       !For water-water only
   R_Coul_cut = 10.0d0                         !Electrostatic cut-off tests only - real space sum
endif
atom_WO='OW'
atom_WH='HW'

!Calculate water-water interactions - .TRUE.
LWatWat = .FALSE.
LProtWat = .FALSE.
if(kwat_wat==1)LWatWat  = .TRUE.
if(kprot_wat==1)LProtWat = .TRUE.

kstr = 1
kntrj = 1          !Number of trajectories - not implemented for more than one trajectory - keep always == 1

!Version 18   if(mcontsamp==1)kstr=0

allocate(PCoul_RES_W(nb_res,nstep/mcontsamp+kstr))                         !Coul Prot-water
allocate(PLJ_RES_W(nb_res,nstep/mcontsamp+kstr))         !van der Waals Prot-water
allocate(Pot_Prot_wat(nstep/mcontsamp+kstr))          !Coul + van der Waals
allocate(pot_coul_PW_m(nstep/mcontsamp+kstr))
allocate(pot_LJ_PW_m(nstep/mcontsamp+kstr))

allocate (pot_WW_coul(nstep/mcontsamp+kstr),pot_WW_vdw(nstep/mcontsamp+kstr),EP_WW(nstep/mcontsamp+kstr))

!xtc2fortran
allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
!end xtc2fortran

!call SYSTEM('mkdir -p Prot_potent_umbrella_out')
inquire(file='Prot_water_potent_out/log_prot_water.dat',exist=filex)
if(.not.filex) call SYSTEM('mkdir -p Prot_water_potent_out')

n0 = 100
!protein-water
n1 = 110
n2 = 120
n3 = 130
n4 = 140
!water-water
n7 = 170
n8 = 180
n9 = 190
n10 = 200

!mean energies as a function of the COM-COM distance from umbrella sampling trajectories
open(n0,file='Prot_water_potent_out/log_prot_water.dat',status='unknown',action='write')
open(n1,file='Prot_water_potent_out/PE_prot_water.dat',status='unknown',action='write')
open(n2,file='Prot_water_potent_out/Coul_prot_water.dat',status='unknown',action='write')
open(n3,file='Prot_water_potent_out/LJ_prot_water.dat',status='unknown',action='write')
open(n4,file='Prot_water_potent_out/PE_PW_time.dat',status='unknown',action='write')

open(n7,file='Prot_water_potent_out/PE_water_water.dat',status='unknown',action='write')
open(n8,file='Prot_water_potent_out/PE_WW_time.dat',status='unknown',action='write')
open(n9,file='Prot_water_potent_out/Coul_water_water.dat',status='unknown',action='write')
open(n10,file='Prot_water_potent_out/LJ_water_water.dat',status='unknown',action='write')

write(n1,'(a1,3x,a6,3x,a24,3x,a13)')'#','RC(nm)','PE_protein-water(kJ/mol)','stdev(kJ/mol)'
write(n2,'(a1,3x,a6,3x,a24,3x,a13)')'#','RC(nm)','Coul_protein-wat(kJ/mol)','stdev(kJ/mol)'
write(n3,'(a1,3x,a6,3x,a24,3x,a13)')'#','RC(nm)','LJ_protein-water(kJ/mol)','stdev(kJ/mol)'
write(n4,'(a1,3x,a8,3x,a24)')'#','time(ns)','PE_protein-water(kJ/mol)'

write(n7,'(a1,3x,a6,3x,a22,3x,a13)')'#','RC(nm)','PE_water-water(kJ/mol)','stdev(kJ/mol)'
write(n8,'(a1,3x,a8,3x,a22)')'#','time(ns)','PE_water-water(kJ/mol)'
write(n9,'(a1,3x,a6,3x,a24,3x,a13)')'#','RC(nm)','Coul_wat-wat(kJ/mol)','stdev(kJ/mol)'
write(n10,'(a1,3x,a6,3x,a24,3x,a13)')'#','RC(nm)','LJ_wat-water(kJ/mol)','stdev(kJ/mol)'


!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3

enconv_kjmol = 1.0d-3*AVSNO                   !conversion factor J ---> kJ/mol

!Version 15
dcm_trj_US = 0.0d0
drcoord_trj_US = 0.0d0

PCoul_RES_W   = 0.0d0
PLJ_RES_W     = 0.0d0
Pot_Prot_wat  = 0.0d0
pot_coul_PW_m = 0.0d0
pot_LJ_PW_m   = 0.0d0

pot_prot_W_mean  = 0.0d0
pot_coul_PW_mean = 0.0d0
pot_LJ_PW_mean   = 0.0d0
stdev_prot_W     = 0.0d0
stdev_coul_PW    = 0.0d0
stdev_LJ_PW      = 0.0d0

!Water-Water
pot_WW_coul     = 0.0d0
pot_WW_vdw      = 0.0d0
EP_WW           = 0.0d0
pot_W_W_mean    = 0.0d0

pot_WW_coul_mean = 0.0d0
pot_WW_vdw_mean  = 0.0d0 

stdev_W_W       = 0.0d0
stdev_WW_pot    = 0.0d0
stdev_W_W_Coul  = 0.0d0
stdev_W_W_LJ    = 0.0d0
stdev_WW_Coul   = 0.0d0
stdev_WW_LJ     = 0.0d0

 cmxsol_I = 0.0d0;  cmysol_I = 0.0d0;  cmzsol_I = 0.0d0
 cmxsol_II = 0.0d0; cmysol_II = 0.0d0; cmzsol_II = 0.0d0

!Version_12

write(*,*)
write(*,*)'Starting protein-water potential energy calculation...'
write(*,*)'Routine prot_water_umbrella'
write(*,*)'Results are printed to Prot_water_potent_out'
write(*,*)'Traj sampling freq for Protein-Water PE (steps) = ',mcontsamp
write(*,*)'Total number of time-steps = ',nstep
write(*,*)'Calculate protein-Water PE ?',LProtWat
write(*,*)'Calculate Water-Water PE ?',LWatWat
write(*,*)

write(n0,*)
write(n0,*)'Starting protein-water potential energy calculation...'
write(n0,*)'Routine prot_water_umbrella'
write(n0,*)'Results are printed to Prot_water_potent_out'
write(n0,*)'Traj sampling freq for Protein-Water PE (steps) = ',mcontsamp
write(n0,*)'Total number of time-steps = ',nstep
write(n0,*)'Calculate protein-Water PE ?',LProtWat
write(n0,*)'Calculate Water-Water PE ?',LWatWat
write(n0,*)

!=======================================================================================================

!Loop over number of trajectories available for calculating mean potential map (electrostatic and vdW interactions)
do ktr = 1,kntrj                                                              !Loop over number of umbrella sampling trajectories
   IF(inputformat.eq.'GROMACS'.and.ktr > 1)call f77_molfile_init              !This corrected the problem of running this routine after other analyses
   write(*,*)'Starting analysis of trajectory ', ktr,' out of ',kntrj
   write(*,*)

!Label different trajectories input files: md_trj01.xtc, md_trj02.xtc, etc
   IF(kntrj>1)THEN
!   if(ktr==1)trjfile='md_trj01.xtc'
!   if(ktr==2)trjfile='md_trj02.xtc'
!   etc

      write (klab,'(I2.2)') ktr                             ! converting integer (ktr) to string (klab) using an 'internal file'
      trjfile='md_trj'//trim(klab)//'.xtc'
!check      write(*,*) trjfile 
   ENDIF

!time counter of electrostatic interactions (1,2,3,...) with 1,2,3 separated by mcontsamp
   it  = 0
   it1 = 0
   it2 = 0
   jt  = 0  
   jt1 = 0
   jt2 = 0
   kt  = 0
   kwt = 0
   icm = 0
   
!xtc2fortran_trj
   
   IF(inputformat.eq.'GROMACS')THEN
      infile = TRIM(trjfile)
      intype = 'auto'
      npart  = -1
      handle(1) = -1
      handle(2) = -1
      handle(3) = -1
      handle(4) = -1
!set up everything and register all static plugins
!NG      call f77_molfile_init

      call f77_molfile_open_read(handle(1),npart,infile,intype)

      if (handle(1).lt.0) then
         print*,'file type unknown or not registered'
         stop
      else
         print*,'file successfully opened:'
         print*,'handle:',handle(1)
         print*,'npart: ',npart
         print*,'nsteps:',nstep
      end if
        
      if(npart /= natms)then
        write(*,*)'error - number of atoms is wrong in solv_inp.dat'
        stop
      endif  
   ENDIF    
   
!end xtc2fortran_trj  

!Start trajectory analysis
   
!Version 18   do j = 1,nstep                                               !Loop over number of steps of trajectory
   do j = 0,nstep   
      if(mod(j,1)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/    
      IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
         if(nens==1)then
            read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
         elseif(nens==0)then
            cell(1) = cube(1)
            cell(2) = cube(2)
            cell(3) = cube(3)
         endif   
!read solute and water coordinates
         read(ninput,*)natcheck
         if(natcheck.ne.natms)then
            write(*,*)'error - number of atoms from solv_inp is wrong'
            stop
         endif   
         do i = 1,natms
            read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!NG         if(j==1)atom_symb(i)=TRIM(atomname(i))
!check         if(j==1.and.atom_symb(i)/='H')write(*,*)atom_symb(i)
         end do
      ELSEIF(inputformat.eq.'GROMACS')THEN
      
!xtc2fortran_trj 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
         status = 1   ! status=1 on entry means read
         call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
         cell(1)=dble(box(1))
         cell(2)=dble(box(2))
         cell(3)=dble(box(3))
         kp=1
!       write(*,*)cell
         do i = 1,npart*3,3
            x(kp)=xyz(i)
            y(kp)=xyz(i+1)
            z(kp)=xyz(i+2)
            kp = kp + 1
         end do   
      ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 


!Calculate proteins centre of mass (COM) distance along time for each trajectory
!Version 18      if(j==1.or.mod(j,mcontsamp)==0)then 
      if(mod(j,mcontsamp)==0)then      
         write(*,*)'Sampling trj - step',j
!check         write(*,*)'calculating com distance...trj ',ktr,' step ',j
         icm    = icm + 1
         kt     = kt + 1                                                                              !pseudo-time counter
         nstart = 1
         nend   = natmsol/2
         call sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol_I,cmysol_I,cmzsol_I)
         nstart = natmsol/2 + 1
         nend   = natmsol
         call sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol_II,cmysol_II,cmzsol_II)
         dcmx = cmxsol_I - cmxsol_II
         dcmy = cmysol_I - cmysol_II
         dcmz = cmzsol_I - cmzsol_II
         if(kMIC_==1)then
            dcmx = dcmx - dnint(dcmx/cell(1))*cell(1)
            dcmy = dcmy - dnint(dcmy/cell(2))*cell(2)
            dcmz = dcmz - dnint(dcmz/cell(3))*cell(3)
         endif
         dcmr2 = dcmx**2 + dcmy**2 + dcmz**2
         dcmr  = dsqrt(dcmr2)
         dcm_trj_US = dcm_trj_US + dcmr                               !accumulate com distance       
         drcoord_trj_US = drcoord_trj_US + dsqrt(dcmx**2.0)           !accumulate RC distance (RC = Reaction Coordinate)

!check for a single trajectory         
!              PSECS = TSTEP*dfloat(j-1)                            !psecs
!check         PSECS = TSTEP*dfloat(icm-1)*dfloat(mcontsamp)
!check         NSECS = 1e-3*PSECS                                   !nsecs
!check         if(ktr==1)write(n1,'(e19.3,3x,f12.4)')NSECS,dcmr     !print check for a single trajectory - to compare with the mean over trajectories 
!NG         if(ktr==1)write(n1,'(e19.3,3x,f12.4,2x,a6,1x,i6)')NSECS,dcmr,'trj = ',ktr  
      endif
!      write(*,*)'Ended COM-COM calculation'

!End COM calculation



!Calculate water-water interactions
      IF(LWatWat)THEN
         if(mod(j,mcontsamp)==0)then
            kwt = kwt + 1                                      !pseudo time-counter
            do kw1=nmolsol*natmsol+1,natms-nions               !loop over water atoms O and Hs 
               if(atomname(kw1)=='OH2')kw_ = kw1 + nwatsites
               if(atomname(kw1)=='H1') kw_ = kw1 + nwatsites-1
               if(atomname(kw1)=='H2') kw_ = kw1 + nwatsites-2
               do kw2=kw_,natms-nions
                  dx = x(kw1)-x(kw2)
                  dy = y(kw1)-y(kw2)
                  dz = z(kw1)-z(kw2)
                  if(kMIC_==1)then
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                  endif 
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
!check                  write(*,*)kw1,kw2                         !check compute half of the interactions 1-2 = 2-1
                  if(atomtype(kw1)==atom_WO)QW1 = chrg_OW*ELCH
                  if(atomtype(kw1)==atom_WH)QW1 = chrg_HW*ELCH
                  if(atomtype(kw2)==atom_WO)QW2 = chrg_OW*ELCH
                  if(atomtype(kw2)==atom_WH)QW2 = chrg_HW*ELCH
                  pot_WW_coul(kwt) = pot_WW_coul(kwt) + enconv_kjmol*PMTTV*QW1*QW2/dr
                  EP_WW(kwt) = EP_WW(kwt) + enconv_kjmol*PMTTV*QW1*QW2/dr
                  if((atomtype(kw1)==atom_WO).and.(atomtype(kw2)==atom_WO))then                                 !vdW - oxygens only
                     pot_WW_vdw(kwt)= pot_WW_vdw(kwt) + (C12_OW/(0.1d0*dr)**12.0 - C6_OW/(0.1d0*dr)**6.0)
                     EP_WW(kwt)     = EP_WW(kwt)      + (C12_OW/(0.1d0*dr)**12.0 - C6_OW/(0.1d0*dr)**6.0)
                  endif                 
               end do   ! End water loop 2
            end do      ! End water loop 1
         endif          ! End sampling
      ENDIF             ! End Water-Water
!End calculate water-water interactions



!Calculate protein-water interactions

!Electrostatics and van der Waals
      IF(LProtWat)THEN
         if(mod(j,mcontsamp)==0)then   
            it1 = it1 +1                         !pseudo-time counter
            k1 = 0                               !Chains ABCD atomic index
            do ic1=1,nb_res/2                    !residues of ABCD
               do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
                  k1 = k1 + 1
                  do kw=nmolsol*natmsol+1,natms-nions             !loop over water atoms O and Hs               
                     dx = x(k1)-x(kw)
                     dy = y(k1)-y(kw)
                     dz = z(k1)-z(kw)
                     if(kMIC_==1)then
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                     endif 
                     dr2 = dx**2 + dy**2 + dz**2
                     dr  = dsqrt(dr2)
                     Q1 = chrg(k1)*ELCH
!check                     if(kw==1)write(*,*)chrg(k1)
                     if(atomtype(kw)==atom_WO)QW = chrg_OW*ELCH
                     if(atomtype(kw)==atom_WH)QW = chrg_HW*ELCH
!                     write(*,*)kw,atomtype(kw),atom_WO,QW
!Warning!!! en is summed for [nstep/mcontsamp] +1 steps                   
!NG test cut-off     IF(dr <= R_Coul_cut) PCoul_RES_W(ic1,it1) = PCoul_RES_W(ic1,it1) + enconv_kjmol*PMTTV*Q1*QW/dr                          !per residue and time-step
                     PCoul_RES_W(ic1,it1) = PCoul_RES_W(ic1,it1) + enconv_kjmol*PMTTV*Q1*QW/dr
                     Pot_Prot_wat(it1)    = Pot_Prot_wat(it1)    + enconv_kjmol*PMTTV*Q1*QW/dr                                        !per time step
!van der Waals                  
                     if(atomtype(kw)==atom_WO)then
!NG                        C6  = dsqrt(C6kjmol(k1)*C6_OW)                    !kJ/mol nm**6  
!NG                        C12 = dsqrt(C12kjmol(k1)*C12_OW)                  !kJ/mol nm**12
!NG                        PLJ_RES_W(ic1,it1) = PLJ_RES_W(ic1,it1) + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0) 
!NG                        Pot_Prot_wat(it1)  = Pot_Prot_wat(it1)  + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0)
                        PLJ_RES_W(ic1,it1) = PLJ_RES_W(ic1,it1) + (C12_IW(k1)/(0.1d0*dr)**12.0 - C6_IW(k1)/(0.1d0*dr)**6.0) 
                        Pot_Prot_wat(it1)  = Pot_Prot_wat(it1)  + (C12_IW(k1)/(0.1d0*dr)**12.0 - C6_IW(k1)/(0.1d0*dr)**6.0)
                     endif
!end van der Waals                  
                  end do                  !end water atoms
               end do                     !end atoms of residue ic1
            end do                        !end number of residues of first protein   
            
            it2 = it2 +1                     !pseudo-time counter
            k2  = natmsol/2                 !Chains EFGH atomic index                  
            do ic2=nb_res/2+1,nb_res       !residues of EFGH
               do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
                  k2 = k2 + 1
                  do kw=nmolsol*natmsol+1,natms-nions             !loop over water atoms O and Hs              
                     dx = x(k2)-x(kw)
                     dy = y(k2)-y(kw)
                     dz = z(k2)-z(kw)
                     if(kMIC_==1)then
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                     endif 
                     dr2 = dx**2 + dy**2 + dz**2
                     dr  = dsqrt(dr2)
                     Q2 = chrg(k2)*ELCH
!check                      if(kw==nmolsol*natmsol+1)write(*,*)chrg(k2)
                     if(atomtype(kw)==atom_WO)QW = chrg_OW*ELCH
                     if(atomtype(kw)==atom_WH)QW = chrg_HW*ELCH
!check                     write(*,*)atomtype(kw),QW
!Warning!!! en is summed for [nstep/mcontsamp] +1 steps                   
!NG test cut-off                     IF(dr <= R_Coul_cut) PCoul_RES_W(ic2,it2) = PCoul_RES_W(ic2,it2) + enconv_kjmol*PMTTV*Q2*QW/dr
                     PCoul_RES_W(ic2,it2) = PCoul_RES_W(ic2,it2) + enconv_kjmol*PMTTV*Q2*QW/dr
                     Pot_Prot_wat(it2)    = Pot_Prot_wat(it2)    + enconv_kjmol*PMTTV*Q2*QW/dr 
!van der Waals                  
                     if(atomtype(kw)==atom_WO)then
!                        C6  = dsqrt(C6kjmol(k2)*C6_OW)                    !kJ/mol nm**6  
!                        C12 = dsqrt(C12kjmol(k2)*C12_OW)                  !kJ/mol nm**12
!                        PLJ_RES_W(ic2,it2) = PLJ_RES_W(ic2,it2) + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0)
!                        Pot_Prot_wat(it2)  = Pot_Prot_wat(it2)  + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0)
                         PLJ_RES_W(ic2,it2) = PLJ_RES_W(ic2,it2) + (C12_IW(k2)/(0.1d0*dr)**12.0 - C6_IW(k2)/(0.1d0*dr)**6.0)
                         Pot_Prot_wat(it2)  = Pot_Prot_wat(it2)  + (C12_IW(k2)/(0.1d0*dr)**12.0 - C6_IW(k2)/(0.1d0*dr)**6.0)
                     endif
!end van der Waals                                    
                  end do                  !end water atoms
               end do                     !end atoms of residue ic1
            end do                        !end number of residues
         endif                            !end electrostatics sampling  
      ENDIF
                  
!End protein-water electrostatics and van der Waals             
               

!NG   !van der Waals 
!NG   !Version 18      if(j==1.or.mod(j,mcontsamp)==0)then
!NG         if(mod(j,mcontsamp)==0)then   
!NG            jt1 = jt1 +1                           !pseudo-time counter
!NG            k1 = 0                               !Chains ABCD atomic index
!NG            do ic1=1,nb_res/2                    !residues of ABCD
!NG               do iats1 = 1,Natresid(ic1)        !atoms of residue ic1
!NG                  k1 = k1 + 1
!NG                  do kw=nmolsol*natmsol+1,natms-nions,nwatsites             !loop over water oxygens - H atoms no van der Waals 
!NG                     dx = x(k1)-x(kw)
!NG                     dy = y(k1)-y(kw)
!NG                     dz = z(k1)-z(kw)
!NG                     if(kMIC_==1)then
!NG                        dx = dx - dnint(dx/cell(1))*cell(1)
!NG                        dy = dy - dnint(dy/cell(2))*cell(2)
!NG                        dz = dz - dnint(dz/cell(3))*cell(3)
!NG                     endif
!NG                     dr2 = dx**2 + dy**2 + dz**2
!NG                     dr  = dsqrt(dr2)
!NG   !Apply geometric rule - units [kJ/mol nm**6]  and [kJ/mol nm**12] - convert to [kJ/mol Ang**6]  and [kJ/mol Ang**12]                   
!NG   !kJ/mol Ang**6                        C6  = 1.0d6*dsqrt(C6kjmol(k1)*C6kjmol(k2))  
!NG   !kJ/mol Ang**12                       C12 = 1.0d12*dsqrt(C12kjmol(k1)*C12kjmol(k2)) 
!NG   !Do not convert to [kJ/mol Ang**6]  and [kJ/mol Ang**12] - convert distances instead 
!NG                     C6  = dsqrt(C6kjmol(k1)*C6_OW)                    !kJ/mol nm**6  
!NG                     C12 = dsqrt(C12kjmol(k1)*C12_OW)                  !kJ/mol nm**12
!NG                     PLJ_RES_W(ic1,jt1) = PLJ_RES_W(ic1,jt1) + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0) 
!NG                     Pot_Prot_wat(jt1)  = Pot_Prot_wat(jt1)  + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0)
!NG                  end do        !end O loop
!NG               end do           !end atoms of residue ic1   
!NG            end do              !end number of residues
         
!NG             jt2 = jt2 +1                           !pseudo-time counter
!NG             k2 = natmsol/2                 !Chains EFGH atomic index                  
!NG             do ic2=nb_res/2+1,nb_res       !residues of EFGH
!NG                do iats2 = 1,Natresid(ic2)  !atoms of residue ic2
!NG                   k2 = k2 + 1
!NG                   do kw=nmolsol*natmsol+1,natms-nions,nwatsites             !loop over water oxygens - H atoms no van der Waals
!NG    !                  write(*,*)atomtype(kw)
!NG                      dx = x(k2)-x(kw)
!NG                      dy = y(k2)-y(kw)
!NG                      dz = z(k2)-z(kw)
!NG                      if(kMIC_==1)then
!NG                         dx = dx - dnint(dx/cell(1))*cell(1)
!NG                         dy = dy - dnint(dy/cell(2))*cell(2)
!NG                         dz = dz - dnint(dz/cell(3))*cell(3)
!NG                      endif
!NG                      dr2 = dx**2 + dy**2 + dz**2
!NG                      dr  = dsqrt(dr2)
!NG    !Apply geometric rule - units [kJ/mol nm**6]  and [kJ/mol nm**12] - convert to [kJ/mol Ang**6]  and [kJ/mol Ang**12]                   
!NG    !kJ/mol Ang**6                        C6  = 1.0d6*dsqrt(C6kjmol(k1)*C6kjmol(k2))  
!NG    !kJ/mol Ang**12                       C12 = 1.0d12*dsqrt(C12kjmol(k1)*C12kjmol(k2)) 
!NG                      C6  = dsqrt(C6kjmol(k2)*C6_OW)                    !kJ/mol nm**6  
!NG                      C12 = dsqrt(C12kjmol(k2)*C12_OW)                  !kJ/mol nm**12
!NG                      PLJ_RES_W(ic2,jt2) = PLJ_RES_W(ic2,jt2) + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0)
!NG                      Pot_Prot_wat(jt2)  = Pot_Prot_wat(jt2)  + (C12/(0.1d0*dr)**12.0 - C6/(0.1d0*dr)**6.0)
!NG                   end do     !end O loop
!NG                end do        !end atoms of residue ic2 
!NG             end do           !end number of residues
!NG          endif               !end vdW sampling 
!NG  !      write(*,*)'Ended vdW calculation' 
!NG  !End van der Waals         

     
!xtc2fortran_trj
      if(inputformat.eq.'GROMACS'.and.status.eq.0) then
         write(*,*)'error on reading trajectory file - step',i
         stop
      endif
!end xtc2fortran_trj    
    
   end do                !end trajectory time-step main loop

!xtc2fortran_trj
   IF(inputformat.eq.'GROMACS')THEN
      call f77_molfile_finish
      write(*,*)'finished reading xtc trj...'
   ENDIF   
!end xtc2fortran_trj
   
end do                   !end trajectories loop 



!Calculate mean energies and standard deviations
IF(LProtWat)THEN
!PE - protein-water
   do it=1,nstep/mcontsamp + kstr
      PSECS = TSTEP*dfloat(it-1)*dfloat(mcontsamp)                   !psecs
      NSECS = 1e-3*PSECS  
      pot_prot_W_mean = pot_prot_W_mean + Pot_Prot_wat(it) 
      write(n4,*)NSECS,Pot_Prot_wat(it) 
   end do
   pot_prot_W_mean = pot_prot_W_mean/dfloat(nstep/mcontsamp + kstr)

!Coulomb and van der Waals 
   do ic12 = 1,nb_res
!check   do ic12 = 1,nb_res/2
      do it=1,nstep/mcontsamp + kstr
         pot_coul_PW_m(it) = pot_coul_PW_m(it) + PCoul_RES_W(ic12,it)                 !total coulomb energy
         pot_LJ_PW_m(it)   =  pot_LJ_PW_m(it) + PLJ_RES_W(ic12,it)                      !total van der Waals 
      end do
   end do  
!Coulomb and van der Waals
   do it=1,nstep/mcontsamp + kstr
       pot_coul_PW_mean = pot_coul_PW_mean + pot_coul_PW_m(it)
       pot_LJ_PW_mean   =  pot_LJ_PW_mean + pot_LJ_PW_m(it)
   end do

   pot_coul_PW_mean = pot_coul_PW_mean/dfloat(nstep/mcontsamp + kstr)
   pot_LJ_PW_mean   = pot_LJ_PW_mean/dfloat(nstep/mcontsamp + kstr)
   
!Calculate standard deviations
   do kt = 1,nstep/mcontsamp + kstr 
         stdev_prot_W  = stdev_prot_W + (Pot_Prot_wat(kt) - pot_prot_W_mean)**2.0 
         
         stdev_coul_PW = stdev_coul_PW + (pot_coul_PW_m(kt) - pot_coul_PW_mean)**2.0
         stdev_LJ_PW   = stdev_LJ_PW + (pot_LJ_PW_m(kt) - pot_LJ_PW_mean)**2.0
   end do
   stdev_PW_pot  = dsqrt(stdev_prot_W/dfloat(nstep/mcontsamp))
   stdev_PW_coul = dsqrt(stdev_coul_PW/dfloat(nstep/mcontsamp))
   stdev_PW_LJ   = dsqrt(stdev_LJ_PW/dfloat(nstep/mcontsamp))
      
!mean reaction coordinate (com-com-X) and com-com distance    
   drcoord_mean = drcoord_trj_US/dfloat(nstep/mcontsamp + kstr)
   dcm_mean = dcm_trj_US/dfloat(nstep/mcontsamp + kstr)
!print results 
   write(n1,*)drcoord_mean*0.1d0,pot_prot_W_mean,stdev_PW_pot                    !nm
   write(n2,*)drcoord_mean*0.1d0,pot_coul_PW_mean,stdev_PW_coul                  !nm
   write(n3,*)drcoord_mean*0.1d0,pot_LJ_PW_mean,stdev_PW_LJ                      !nm
ENDIF
      
! Water-Water potential energy

IF(LWatWat)THEN
   do kt = 1,nstep/mcontsamp + kstr
      PSECS = TSTEP*dfloat(kt-1)*dfloat(mcontsamp)                   !psecs
      NSECS = 1e-3*PSECS  
      pot_W_W_mean = pot_W_W_mean + EP_WW(kt)
      write(n8,*)NSECS,EP_WW(kt) 
!Coulomb and van der Waals     
      pot_WW_coul_mean = pot_WW_coul_mean + pot_WW_coul(kt)
      pot_WW_vdw_mean  = pot_WW_vdw_mean  + pot_WW_vdw(kt)
      
   end do
   pot_W_W_mean = pot_W_W_mean/dfloat(nstep/mcontsamp + kstr)
   pot_WW_coul_mean = pot_WW_coul_mean/dfloat(nstep/mcontsamp + kstr)
   pot_WW_vdw_mean  = pot_WW_vdw_mean/dfloat(nstep/mcontsamp + kstr)
      
   do kt = 1,nstep/mcontsamp + kstr 
      stdev_W_W = stdev_W_W + (EP_WW(kt) - pot_W_W_mean)**2.0 
      stdev_W_W_Coul = stdev_W_W_Coul + (pot_WW_coul(kt) - pot_WW_coul_mean)**2.0 
      stdev_W_W_LJ   = stdev_W_W_LJ   + (pot_WW_vdw(kt) - pot_WW_vdw_mean)**2.0   
   end do
   stdev_WW_pot   = dsqrt(stdev_W_W/dfloat(nstep/mcontsamp))
   stdev_WW_Coul = dsqrt(stdev_W_W_Coul/dfloat(nstep/mcontsamp))
   stdev_WW_LJ = dsqrt(stdev_W_W_LJ/dfloat(nstep/mcontsamp))

!mean reaction coordinate (com-com-X) and com-com distance    
   drcoord_mean = drcoord_trj_US/dfloat(nstep/mcontsamp + kstr)
   dcm_mean = dcm_trj_US/dfloat(nstep/mcontsamp + kstr)   
!print results   
   write(n7,*)drcoord_mean*0.1d0,pot_W_W_mean,stdev_WW_pot
   write(n9,*)drcoord_mean*0.1d0,pot_WW_coul_mean,stdev_WW_Coul
   write(n10,*)drcoord_mean*0.1d0,pot_WW_vdw_mean,stdev_WW_LJ   
ENDIF

deallocate(PCoul_RES_W)       !Coul Prot-water
deallocate(PLJ_RES_W)         !van der Waals Prot-water
deallocate(Pot_Prot_wat)      !Coul + van der Waals
deallocate(pot_coul_PW_m)
deallocate(pot_LJ_PW_m) 

deallocate (pot_WW_coul)
deallocate(pot_WW_vdw)
deallocate(EP_WW)


return

END SUBROUTINE prot_water_umbrella




SUBROUTINE wat_env_diffusion(trjfile,ninput,unfoldtrjfile,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nions,nwatsites,atomname,&
                               msd_delay,difsampw,difbulk,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO)                               
! Calculate water's mean square displacement and Einstein diffusion coefficient - Solutions
! Every difsampw time-steps the code samples waters on different environments 
    integer,intent(in)                               :: ninput,nbox,nsyst
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
!    character(len=4),dimension(natms),intent(out)    :: atomname
    character(len=4),dimension(natms)                :: atomname
    
    integer,intent(in)                               :: msd_delay     ! msd delay time (time-windows/ps)
    integer,intent(in)                               :: difsampw      ! time freq. to sample water in radial shells for diffusion calculation
    real(kind=8)                                     :: difbulk       ! onset (Ang) for bulk water relative to the solute COM - for bulk diffusion       
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
   
! Local variables
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=4),dimension(:),allocatable            :: x_u,y_u,z_u
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    real(kind=8)                                     :: cell_u(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8),dimension(:,:,:),allocatable        :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable        :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                                     :: TSTEP,PSECS       
    real(kind=8)                                     :: dx,dy,dz,dr,dr2
    integer                                          :: NDELS,KINTVL,NDELS1
    integer,dimension(:,:,:),allocatable             :: NH2OidO,NOHidH
    integer,dimension(:,:),allocatable               :: ntime
    integer,dimension(:,:),allocatable               :: nH2Oorg
    integer,dimension(:),allocatable                 :: kv
    integer,dimension(:,:),allocatable               :: kwlabel      !index of waters in each environment
    integer      :: n0, n1, n2, n3,n4, n5, nf
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer,dimension(:),allocatable           :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable           :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                               :: dr_bulk, dr4nb
    real(kind=8)                               :: RADMAX,RADMIN
    integer                                    :: norgmax
    integer                                    :: NRSHMAX   
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!MSD   
    real(kind=8),dimension(:),allocatable        :: msdx_b,  msdy_b,  msdz_b                         !bulk
    real(kind=8),dimension(:),allocatable        :: msdx_hs, msdy_hs, msdz_hs                        !hsh
    real(kind=8),dimension(:),allocatable        :: msdx_th, msdy_th, msdz_th                        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: msdx_nth,msdy_nth,msdz_nth                       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: msd_mean,msd_mean_x,msd_mean_y,msd_mean_z        
    integer,dimension(:),allocatable             :: nH2Omean
    real(kind=8),dimension(:),allocatable        :: stdevnOH_sum,stdevnOH
    integer                                      :: bulk_flag,kr_start
    integer                                      :: at_id                            ! solute atomic index to analyse the respective solvation environment  
!Least Squares Fit    
    integer                                      :: LSF_begin,LSF_end
    integer                                      :: n_LSF_begin,n_LSF_end    
    real(kind=8)                                 :: a,b,d,r                    !Y = a + b*X; d = standard deviation of the fit; r = correlation coeficient
    integer                                      :: n                          !number of points used in the fit: range ]n_LSF_begin,n_LSF_end[
    integer                                      :: maxwat_solv                !maximum number of solvation waters allocated
    integer                                      :: maxwat_bulk                !maximum number of waters in the bulk
    integer,dimension(:),allocatable             :: koc                  !Nb of origins corrected for origins where no waters were found in a given environment
    integer                                      :: nroute               !integer to decide which files to open to compute the Einstein diffusion from LSF

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
    character(len=14),intent(in)            :: unfoldtrjfile                           ! internal - unfolded trajectory file name for Einstein diffusion trj_unfold.xtc
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
    
    integer*4                               :: npart_u, maxatom_u, handle_u(4), status_u
    real(kind=4)                            :: box_u(6)
    real(kind=4),dimension(:),allocatable   :: xyz_u
    character                               :: infile_u*200, intype_u*10
    integer                                 :: kp_u    
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

allocate(xyz_u(natms*3))
allocate(x_u(natms),y_u(natms),z_u(natms))

!Open PBC trajectory - used to find solvation waters 

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'PBC trj file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    


!Version 12 - MSD solvation bug fix
!Open unfolded trajectory - used to calculate the msd 

IF(inputformat.eq.'GROMACS')THEN
   infile_u = TRIM(unfoldtrjfile)
   intype_u = 'auto'
   npart_u  = -1
   handle_u(1) = -1
   handle_u(2) = -1
   handle_u(3) = -1
   handle_u(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle_u(1),npart_u,infile_u,intype_u)
   if (handle_u(1).lt.0) then
      print*,'unfolded trj file type unknown or not registered'
      stop
   else
      print*,'unfolded trj file successfully opened:'
      print*,'handle:',handle_u(1)
      print*,'npart: ',npart_u
      print*,'nsteps:',nstep
   end if
     
   if(npart_u /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF 
!End Version 12 MSD solvation bug fix


!end xtc2fortran_trj    
      
maxwat_solv = 10000                           !for allocation purposes - maximum number of waters in solvation shells (move to input if needed)
if(nmolwat<maxwat_solv)maxwat_solv = nmolwat
! maxwat_bulk = 1000
maxwat_bulk = 500                    !Water molecules - not OH groups
if(nmolwat<maxwat_bulk)maxwat_bulk = nmolwat

!define internally range for Linear Squares Fit (ps)

!LSF_begin = 400                               !Start fit after 10 ps         (move to input if needed)
!LSF_begin = 20
!LSF_begin = 500
LSF_begin = msd_delay/2
LSF_end = msd_delay - 10                     !End fit 10 ps before msd end  (move to input if needed)
nroute = 0                                   !Keep always zero (0)

n_LSF_begin = IDINT(LSF_begin*1000d0/dt) 
n_LSF_end = IDINT(LSF_end*1000d0/dt) 
    
inquire(file='MSD_env_out/log_env_msd.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p MSD_env_out')    

! Warning!!! Do not change the following numbers (n1 to n4) - used in the routine msd_LSF_Einstein
n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140
n5 = 150

open(n0,file='MSD_env_out/log_env_msd.dat',status='unknown',action='write')
open(n1,file='MSD_env_out/msd_bulk.dat'   ,status='unknown',action='write')
open(n2,file='MSD_env_out/msd_hsh.dat'    ,status='unknown',action='write')
open(n3,file='MSD_env_out/msd_th.dat'     ,status='unknown',action='write')
open(n4,file='MSD_env_out/msd_nth.dat'    ,status='unknown',action='write') 
open(n5,file='MSD_env_out/msd_env.dat'    ,status='unknown',action='write') 

!time-window for calculation of the msd
! IDINT - convert to integer - truncate

NDELS = IDINT(msd_delay*1000d0/dt)                  !sets the delay time (in frames)

dr_bulk = difbulk                                   !set distance from the solute com to assume bulk water
bulk_flag = 1                                       !default - calculate diffusion for bulk water
kr_start = 1                                        !environments loop starter - if bulk_flag = 0; kr_start = 2 

if(dr_bulk==0.0)then
   bulk_flag = 0
   write(*,*)
   write(*,*)'Warning!!! No bulk region defined for diffusion calculation'
   write(*,*)'Bulk diffusion will not be calculated accordingly'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! No bulk region defined for diffusion calculation'
   write(n0,*)'Bulk diffusion will not be calculated accordingly'
   write(n0,*)
else 
   write(*,*)
   write(*,*)'Warning!!! Maximum bulk waters for diffusion calculation = ',maxwat_bulk
   write(*,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! Maximum bulk waters for diffusion calculation = ',maxwat_bulk
   write(n0,*)'IF MORE WATERS ARE REQUIRED INCREASE VARIABLE maxwat_bulk'
   write(n0,*)
endif   

NRSHMAX = 4                                      !maximum number of environments
if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/difsampw                                  !for memory allocation purposes

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))  
NTDW4O = 0; NONTDW = 0 

!mean square displacement
allocate(nH2Oorg(norgmax,NRSHMAX))
allocate(msdx_b(NDELS),  msdy_b(NDELS),  msdz_b(NDELS))
allocate(msdx_hs(NDELS), msdy_hs(NDELS), msdz_hs(NDELS))
allocate(msdx_th(NDELS) ,msdy_th(NDELS), msdz_th(NDELS))
allocate(msdx_nth(NDELS),msdy_nth(NDELS),msdz_nth(NDELS))
allocate(msd_mean(NDELS,NRSHMAX),msd_mean_x(NDELS,NRSHMAX),msd_mean_y(NDELS,NRSHMAX),msd_mean_z(NDELS,NRSHMAX))
msdx_b=0.0d0;   msdy_b=0.0d0;   msdz_b=0.0d0
msdx_hs=0.0d0;  msdy_hs=0.0d0;  msdz_hs=0.0d0
msdx_th=0.0d0;  msdy_th=0.0d0;  msdz_th=0.0d0
msdx_nth=0.0d0; msdy_nth=0.0d0; msdz_nth=0.0d0
msd_mean=0.0d0;msd_mean_x=0.0d0;msd_mean_y=0.0d0;msd_mean_z=0.0d0


!allocate(RXt0(norgmax,nmolwat,NRSHMAX),RYt0(norgmax,nmolwat,NRSHMAX),RZt0(norgmax,nmolwat,NRSHMAX))
!allocate(RXt(norgmax,nmolwat,NRSHMAX),RYt(norgmax,nmolwat,NRSHMAX),RZt(norgmax,nmolwat,NRSHMAX))
!allocate(NH2OidO(norgmax,nmolwat,NRSHMAX))

allocate(RXt0(norgmax,maxwat_solv,NRSHMAX),RYt0(norgmax,maxwat_solv,NRSHMAX),RZt0(norgmax,maxwat_solv,NRSHMAX))
allocate(RXt(norgmax,maxwat_solv,NRSHMAX),RYt(norgmax,maxwat_solv,NRSHMAX),RZt(norgmax,maxwat_solv,NRSHMAX))
allocate(NH2OidO(norgmax,maxwat_solv,NRSHMAX))

allocate(ntime(norgmax,NRSHMAX))
ntime = 1; nH2Oorg = 0
allocate(kv(NRSHMAX))
allocate(kwlabel(NRSHMAX,nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled) 
allocate(nH2Omean(NRSHMAX))
allocate(stdevnOH_sum(NRSHMAX),stdevnOH(NRSHMAX))
allocate(koc(NRSHMAX))
nH2Omean = 0
stdevnOH_sum = 0
stdevnOH = 0
kwlabel = 0
koc = 0

 cmxsol = 0.0d0; cmysol=0.0d0; cmzsol=0.0d0

write(*,*)'Solute atoms and HSh onset/cut-off (Ang)'
write(n0,*)'Solute atoms and HSh onset/cut-off (Ang)'
do i = 1,natHSh
   write(*,'(i6,1x,f7.2,a1,f7.2)')natsol_id(i),oatsol_hs(i),'/',ratsol_hs(i)
   write(n0,'(i6,1x,f7.2,a1,f7.2)')natsol_id(i),oatsol_hs(i),'/',ratsol_hs(i)
end do
write(*,*)
write(n0,*)

! Start msd calculation
write(*,*)'Starting water environment diffusion calculation...'
write(*,*)'Routine wat_env_diffusion'
write(*,*)'PBC trajectory = ',trjfile,' used to determine solvation waters'
write(*,*)'unfolded trajectory = ',unfoldtrjfile,' used to calculated the msd'
write(*,*)'Results are printed to MSD_env_out'
write(*,*)'water diffusion delay time-window (ps) =',msd_delay
write(*,*)'msd water sampling frequency (steps) = ',difsampw
write(*,*)'total number of time-steps =',nstep
write(*,*)'msd delay-time (steps) =',NDELS
write(*,*)
write(n0,*)'Starting water environment diffusion calculation...'
write(n0,*)'Routine wat_env_diffusion'
write(n0,*)'PBC trajectory = ',trjfile,' used to determine solvation waters'
write(n0,*)'unfolded trajectory = ',unfoldtrjfile,' used to calculated the msd'
write(n0,*)'Results are printed to MSD_env_out'
write(n0,*)'water diffusion delay time-window (ps) =',msd_delay
write(n0,*)'msd water sampling frequency (steps) = ',difsampw
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'msd delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!check      if(j==1)write(*,'(i5,1x,A4,1x,3(F14.5))')natmid(i),atomname(i),x(i),y(i),z(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj

!PBC trajectory
!==============

!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      

!Unfolded trajectory
!===================

!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status_u = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle_u(1),npart_u,xyz_u(1),box_u,status_u);
    
    cell_u(1)=dble(box_u(1))
    cell_u(2)=dble(box_u(2))
    cell_u(3)=dble(box_u(3))
    kp_u=1
    do i = 1,npart_u*3,3
       x_u(kp_u)=xyz_u(i)
       y_u(kp_u)=xyz_u(i+1)
       z_u(kp_u)=xyz_u(i+2)
       kp_u = kp_u + 1
    end do   

   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules in specific environments
   if((j==1.or.mod(j,difsampw)==0).and.(j.le.nstep-NDELS))then               !time origin: sample waters in specific environments
!NG      if(j==1.or.mod(j,difsampw)==0)then                                  !time origin: sample waters in specific environments
      NTDW4O = 0; NONTDW = 0
      write(*,*)'time-step',j,' sampling water molecules'
      kv = 0                                                                 !waters in each environment 
      ko = ko + 1                                                            !number of origins = number of blocks used to compute the tcfs mean 
      kwlabel = 0
   
!Find solute centre of mass for bulk water tcfs calculation     

      if(bulk_flag==1)call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water      
         dx = x(i) - cmxsol
         dy = y(i) - cmysol
         dz = z(i) - cmzsol
!         dx = dx - anint(dx/cell(1))*cell(1)
!         dy = dy - anint(dy/cell(2))*cell(2)
!         dz = dz - anint(dz/cell(3))*cell(3)
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
!check         write(*,*)'distance to the centre of mass = ',dr
         if(bulk_flag==0)dr_bulk = 100.d0*cell(1)                            !define a distance larger than the box "radius"
!Version 13 - impose limit of waters in the bulk region         
!         if(dr>=dr_bulk)then
         if(dr>=dr_bulk.and.kv(1)<maxwat_bulk)then
            if(bulk_flag==0)write(*,*)'Warning!!! Found a water molecule in the bulk'
!bulk msd(t0)             
            nradsh = 1                                       !bulk environment
            kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
            kvi = kv(nradsh)                                 !O number
!check            write(*,*)'bulk',j,i,kvi
            nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
            NH2OidO(ko,kvi,nradsh) = i                       ! Oxygen id
            RXt0(ko,kvi,nradsh) = x_u(i)                 !origin = ko, O number = kvi, environment type = nradsh
            RYt0(ko,kvi,nradsh) = y_u(i)
            RZt0(ko,kvi,nradsh) = z_u(i)     
            if(kvi >= maxwat_solv)then
               write(*,*)'error - too many waters in the Bulk'
               write(*,*)'action - increase maxwat_solv and re-compile'
               write(*,*)'or avoid bulk calculation - make cut-off = 0 '
               stop
            endif   
!bulk msd(t0) end                     
         endif 
         
!Hydration shell - tetrahedral and non tetrahedral waters in the first HSh         
         do at_id=1,natHSh                                      !Loop over atomic species to analyze the environment         
            k = natsol_id(at_id)
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2.0d0 + dy**2.0d0 + dz**2.0d0
            dr  = dsqrt(dr2)
!Hydration shell msd(t0)            
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then       !water within hydration shell radius            
!Avoid water molecules repetition
               nradsh = 2                                               !hydration shell environment
               if(kwlabel(nradsh,io)==0)then
                  kwlabel(nradsh,io)=1                                                       !Oxygen id for correction of the nb of waters
               else
                  goto 10                                                                    !sample another water molecule
               endif   
!End avoid water molecules repitition 
!NG               nradsh = 2                                     !hydration shell environment
               kv(nradsh) = kv(nradsh) + 1                       !number of water OH groups (index)
               kvi = kv(nradsh)                                  !OH number
!check               write(*,*)'HSh',j,i,iH,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1       !number of OH groups in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                        ! Oxygen id
               RXt0(ko,kvi,nradsh) = x_u(i)                        !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y_u(i)
               RZt0(ko,kvi,nradsh) = z_u(i)     
               if(kvi >= maxwat_solv)then
                  write(*,*)'error - too many waters in the HSh'
                  write(*,*)'action - increase maxwat_solv and re-compile'
                  stop
               endif   
!Hydration shell msd(t0) end                              

!Find tetrahedral and non-tetrahedral waters in the hydration shell for sub-ensemble tetrahedral and non-tetrahedral msd        
!FIND THE FOUR NEAREST O ATOMS J TO EACH OXYGEN I
               RADMAX  = 25.0D0
               RADMIN  =  0.0D0        
               JT      = 1
               DO WHILE(JT.LE.4) 
                  do kw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
                     if(kw.ne.i)then
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOMS i AND OXYGEN ATOMS kw  
                        dx = x(i) - x(kw)
                        dy = y(i) - y(kw)
                        dz = z(i) - z(kw)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        IF((dr>RADMIN).AND.(dr<RADMAX))THEN
                           RADMAX = dr
!NTDW4O STORES THE O ATOMS NEIGHBORS OF i: JT = 1 FIRST NEIGHBOR; JT = 2 SECOND NEIGHBOR; JT = 3 THIRD NEIGHBOR; JT = 4 FOURTH NEIGHBOR
                           NTDW4O(JT) = kw
                        ENDIF
                     endif
                  end do
!check                  write(*,*)j,i,io,JT,NTDW4O(JT)
                  JT = JT + 1
                  RADMIN = RADMAX
                  RADMAX = 25.0D0
               END DO
!CHECK IF SOLUTE ATOMS ARE CLOSER THAN THE 4th NEAREST O ATOM          
               KI=NTDW4O(4)                                       !Oxygen ID of the 4th nearest water neighbor 
               dx = x(i) - x(KI)
               dy = y(i) - y(KI)
               dz = z(i) - z(KI)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               dr4nb = dr                                           !distance to the fourth nearest neighbor
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY SOLUTE ATOM IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = 1                           !run over all solute atoms - heavy and non-heavy atoms
               DO WHILE(L.LE.nmolsol*natmsol)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN             
                     NONTDW(io) = 1                               !water i is non-tetrahedral
                     goto 5
                  ENDIF
                  L = L+1
               END DO     
               5 continue           
               if(NONTDW(io)==1)then                              !Non-tetrahedral water               
!Hydration shell - non tetrahedral msd(t0)             
                  nradsh = 4                                      !hydration shell (non-tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                     !number of oxygens (index)
                  kvi = kv(nradsh)                                !O number
!check                  write(*,*)'non-tetrah',j,i,kvi
                  nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1     !number of waters in environment nradsh at origin ko
                  NH2OidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  RXt0(ko,kvi,nradsh) = x_u(i)                      !origin = ko, O number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = y_u(i)
                  RZt0(ko,kvi,nradsh) = z_u(i)     
!Hydration shell - non tetrahedral msd(t0) end                                           
               elseif(NONTDW(io)==0)then                          !Tetrahedral water 
!Hydration shell - tetrahedral msd(t0)             
                  nradsh = 3                                      !hydration shell (tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                     !number of oxygens (index)
                  kvi = kv(nradsh)                                !O number
!check                  write(*,*)'tetrah',j,i,kvi
                  nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1     !number of oxygens in environment nradsh at origin ko
                  NH2OidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  RXt0(ko,kvi,nradsh) = x_u(i)                      !origin = ko, O number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = y_u(i)
                  RZt0(ko,kvi,nradsh) = z_u(i)             
!Hydration shell - tetrahedral msd(t0) end  
               endif
!            
            endif                                                    !end water is in the coordination sphere      
         end do                                                      !end loop over atomic species
      10 continue   
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nH2Oorg(ko,1),nH2Oorg(ko,2),nH2Oorg(ko,3),nH2Oorg(ko,4)
   endif                                                             !end sampling new time-origin      
   
!calculate msd - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nH2Oorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then               !ntime counts the delay times for which the tcf has already been calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nH2Oorg(kt,kr)                                      !loop over water molecules in environment kr: nH2Oorg can be zero for a given time-origin 
               i=NH2OidO(kt,iw,kr)
               RXt(kt,iw,kr) = x_u(i)
               RYt(kt,iw,kr) = y_u(i)
               RZt(kt,iw,kr) = z_u(i)
               nt = ntime(kt,kr)
               if(kr==1)then                                                                                 !bulk
                  msdx_b(nt) = msdx_b(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_b(nt) = msdy_b(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_b(nt) = msdz_b(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==2)then                                                                             !hshell
                  msdx_hs(nt) = msdx_hs(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_hs(nt) = msdy_hs(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_hs(nt) = msdz_hs(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==3)then                                                                             !hshell-tetrahedral
                  msdx_th(nt) = msdx_th(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_th(nt) = msdy_th(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_th(nt) = msdz_th(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==4)then                                                                             !hShell-non-tetrahedral
                  msdx_nth(nt) = msdx_nth(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_nth(nt) = msdy_nth(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_nth(nt) = msdz_nth(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate average msd over different origins
do kr = 1,NRSHMAX                                      !loop over environments
   do it = 1,NDELS
      if(kr==1)then
         msd_mean_x(it,kr) = msdx_b(it)
         msd_mean_y(it,kr) = msdy_b(it)
         msd_mean_z(it,kr) = msdz_b(it)
         msd_mean(it,kr)   = (msdx_b(it)+msdy_b(it)+msdz_b(it))      
      elseif(kr==2)then
         msd_mean_x(it,kr) = msdx_hs(it)
         msd_mean_y(it,kr) = msdy_hs(it)
         msd_mean_z(it,kr) = msdz_hs(it)
         msd_mean(it,kr)   = (msdx_hs(it)+msdy_hs(it)+msdz_hs(it))
      elseif(kr==3)then
         msd_mean_x(it,kr) = msdx_th(it)
         msd_mean_y(it,kr) = msdy_th(it)
         msd_mean_z(it,kr) = msdz_th(it)
         msd_mean(it,kr)   = (msdx_th(it)+msdy_th(it)+msdz_th(it))
      elseif(kr==4)then
         msd_mean_x(it,kr) = msdx_nth(it)
         msd_mean_y(it,kr) = msdy_nth(it)
         msd_mean_z(it,kr) = msdz_nth(it)
         msd_mean(it,kr)   = (msdx_nth(it)+msdy_nth(it)+msdz_nth(it))
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of O atoms on each environment, kr, over the number of origins, ko 
!nH2Oorg(kt,kr) - number of waters found in origin (block) kt and in the environment kr               
do kr = 1,NRSHMAX
   do kt = 1, ko
      nH2Omean(kr) = nH2Omean(kr) + nH2Oorg(kt,kr)
      if(nH2Oorg(kt,kr)>0)then
         koc(kr) = koc(kr) + 1                                         !number of origins where at least a single water was found in a given environment
      endif
   end do
end do

!Version 21
!standard deviation for the number of waters on each environment
do kr = 1,NRSHMAX
   do kt = 1, ko
      stdevnOH_sum(kr) = stdevnOH_sum(kr) + ( dfloat(nH2Oorg(kt,kr)) - dfloat(nH2Omean(kr))/dfloat(ko) )**2.0 
   end do      
end do
!End Version 21

if(bulk_flag==0)kr_start=2                                        !do not print the bulk

do kr = kr_start,NRSHMAX
   if(ko > 1)stdevnOH(kr) = dsqrt(stdevnOH_sum(kr)/dfloat(ko-1))
   if(ko == 1)stdevnOH(kr) = dsqrt(stdevnOH_sum(kr)/dfloat(ko))
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a45,1x,i7)')'#Number of origins where waters were found = ',koc(kr)
   write(*,'(a25,1x,F14.2,1x,a7,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko),'stdev =',stdevnOH(kr)  
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a45,1x,i7)')'#Number of origins where waters were found = ',koc(kr)
   write(n0,'(a25,1x,F14.2,1x,a7,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko),'stdev =',stdevnOH(kr)  
end do
   write(*,*)
   write(n0,*)

deallocate(msdx_b,  msdy_b,  msdz_b)
deallocate(msdx_hs, msdy_hs, msdz_hs)
deallocate(msdx_th ,msdy_th, msdz_th)


!Print average msd     
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==1)write(n1,'(a5)')'#Bulk'
   if(kr==2)write(n2,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==2)write(n2,'(a4)')'#HSh'
   if(kr==3)write(n3,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==3)write(n3,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n4,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==4)write(n4,'(a15)')'#HSh-non-tetrah'
   
   if(kr==1)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==1)write(n5,'(a5)')'#Bulk'
   if(kr==2)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==2)write(n5,'(a4)')'#HSh'
   if(kr==3)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==3)write(n5,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==4)write(n5,'(a15)')'#HSh-non-tetrah'
   do it = 1,NDELS
      jt = it -1                                                     !print msd starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      
!check       write(*,*)kr,it,koc(kr)
      
!NG      if(kr==1)WRITE(n1,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==2)WRITE(n2,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==3)WRITE(n3,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==4)WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(ko)
     
!NG      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==2)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==3)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==4)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      
      if(kr==1)WRITE(n1,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      if(kr==2)WRITE(n2,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      if(kr==3)WRITE(n3,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      if(kr==4)WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      
      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      if(kr==2)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      if(kr==3)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      if(kr==4)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
      
   end do
   if(kr==1)write(n1,*)
   if(kr==2)write(n2,*)
   if(kr==3)write(n3,*)
   if(kr==4)write(n4,*)
!   
   if(kr==1)write(n5,*)
   if(kr==2)write(n5,*)
   if(kr==3)write(n5,*)
   if(kr==4)write(n5,*)
end do

close(n1) 
close(n2) 
close(n3) 
close(n4)
close(n5)

!Calculate the Einstein self-diffusion coefficient - MSD Linear Squares fit
write(n0,*)
write(n0,*)'Einstein self-diffusion coefficient'
write(n0,*)
do kr = kr_start,NRSHMAX
!read data points in msd_LSF_Einstein for the LSF
   if(kr==1)then
      nf = n1
      write(n0,'(a5)')'#Bulk'
   elseif(kr==2)then
      nf = n2
      write(n0,'(a4)')'#HSh'
   elseif(kr==3)then
      nf = n3
      write(n0,'(a11)')'#HSh-tetrah'
   elseif(kr==4)then
      nf = n4
      write(n0,'(a15)')'#HSh-non-tetrah'
   endif  
   
   call msd_LSF_Einstein(nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end,a,b,d,r,n)
           
   if(kr==1)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   elseif(kr==2)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   elseif(kr==3)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   elseif(kr==4)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   endif   
end do 

close(n0) 

deallocate(NTDW4O,NONTDW)
deallocate(nH2Oorg,NH2OidO)
deallocate(stdevnOH_sum,stdevnOH)
deallocate(ntime)
deallocate(msd_mean,msd_mean_x,msd_mean_y,msd_mean_z)
deallocate(kv,nH2Omean)
deallocate(koc)

deallocate(xyz)
deallocate(x,y,z)
  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
      return

END SUBROUTINE wat_env_diffusion


!WARNING This routine is no longer called and it contains a bug - the unfolded trajectory is used to determine the solvation waters using MIC - THIS IS A BUG
!ORIGINALLY THIS ROUTINE WAS DESIGNED TO ALLOW A BETTER STATISTICS AT SHORT TIMES 
SUBROUTINE wat_env_diff_stat(trjfile,ninput,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nions,nwatsites,atomname,&
                             msd_delay,difsampw,difbulk,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO)                               
! Calculate water's mean square displacement and Einstein diffusion coefficient - Solutions
! Every difsampw time-steps the code samples waters on different environments 
    integer,intent(in)                               :: ninput,nbox,nsyst
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
!    character(len=4),dimension(natms),intent(out)    :: atomname
    character(len=4),dimension(natms)                :: atomname
    
    integer,intent(in)                               :: msd_delay     ! msd delay time (time-windows/ps)
    integer,intent(in)                               :: difsampw      ! time freq. to sample water in radial shells for diffusion calculation
    real(kind=8)                                     :: difbulk       ! onset (Ang) for bulk water relative to the solute COM - for bulk diffusion       
    integer,intent(in)                               :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in) :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)     :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    integer,intent(in)                               :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)             :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
    real(kind=8),dimension(natHsh),intent(in)        :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)        :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
    real(kind=8),dimension(natmsol),intent(in)       :: ZMBIO
   
! Local variables
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                                     :: cmxsol,cmysol,cmzsol      !solute centre of mass
    real(kind=8),dimension(:,:,:),allocatable        :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable        :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                                     :: TSTEP,PSECS       
    real(kind=8)                                     :: dx,dy,dz,dr,dr2
    integer                                          :: NDELS,KINTVL,NDELS1
    integer,dimension(:,:,:),allocatable             :: NH2OidO,NOHidH
    integer,dimension(:,:),allocatable               :: ntime
    integer,dimension(:,:),allocatable               :: nH2Oorg
    integer,dimension(:),allocatable                 :: kv
    integer,dimension(:,:),allocatable               :: kwlabel      !index of waters in each environment
    integer      :: n0, n1, n2, n3,n4, n5, nf
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer,dimension(:),allocatable           :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable           :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                               :: dr_bulk, dr4nb
    real(kind=8)                               :: RADMAX,RADMIN
    integer                                    :: norgmax
    integer                                    :: NRSHMAX   
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!MSD   
    real(kind=8),dimension(:),allocatable        :: msdx_b,  msdy_b,  msdz_b                         !bulk
    real(kind=8),dimension(:),allocatable        :: msdx_hs, msdy_hs, msdz_hs                        !hsh
    real(kind=8),dimension(:),allocatable        :: msdx_th, msdy_th, msdz_th                        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: msdx_nth,msdy_nth,msdz_nth                       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: msd_mean,msd_mean_x,msd_mean_y,msd_mean_z        
    integer,dimension(:),allocatable             :: nH2Omean
    integer                                      :: bulk_flag,kr_start
    integer                                      :: at_id                            ! solute atomic index to analyse the respective solvation environment  
!Least Squares Fit    
    integer                                      :: LSF_begin,LSF_end
    integer                                      :: n_LSF_begin,n_LSF_end    
    real(kind=8)                                 :: a,b,d,r                    !Y = a + b*X; d = standard deviation of the fit; r = correlation coeficient
    integer                                      :: n                          !number of points used in the fit: range ]n_LSF_begin,n_LSF_end[
    integer                                      :: maxwat_solv                !maximum number of solvation waters allocated
    integer,dimension(:),allocatable             :: koc                        !Nb of origins corrected for origins where no waters were found in a given environment
    integer,dimension(:,:),allocatable           :: k_org                      !number of origins sampled for each delay-time
    integer                                      :: no
    integer                                      :: nroute               !integer to decide which files to open to compute the Einstein diffusion from LSF

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
      
maxwat_solv = 10000                           !for allocation purposes - maximum number of waters in solvation shells (move to input if needed)
if(nmolwat<maxwat_solv)maxwat_solv = nmolwat

!define internally range for Linear Squares Fit (ps)

LSF_begin = 400                                 !Start fit after 10 ps         (move to input if needed)
LSF_end = msd_delay - 10                       !End fit 10 ps before msd end  (move to input if needed)
nroute = 1                                     !Keep always one (1)

n_LSF_begin = IDINT(LSF_begin*1000d0/dt) 
n_LSF_end = IDINT(LSF_end*1000d0/dt) 
    
inquire(file='MSD_env_out_stat/log_env_msd.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p MSD_env_out_stat')    

! Warning!!! Do not change the following numbers (n1 to n4) - used in the routine msd_LSF_Einstein
n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140
n5 = 150

open(n0,file='MSD_env_out_stat/log_env_msd.dat',status='unknown',action='write')
open(n1,file='MSD_env_out_stat/msd_bulk.dat'   ,status='unknown',action='write')
open(n2,file='MSD_env_out_stat/msd_hsh.dat'    ,status='unknown',action='write')
open(n3,file='MSD_env_out_stat/msd_th.dat'     ,status='unknown',action='write')
open(n4,file='MSD_env_out_stat/msd_nth.dat'    ,status='unknown',action='write') 
open(n5,file='MSD_env_out_stat/msd_env.dat'    ,status='unknown',action='write') 

!time-window for calculation of the msd
! IDINT - convert to integer - truncate

NDELS = IDINT(msd_delay*1000d0/dt)                  !sets the delay time (in frames)

dr_bulk = difbulk                                   !set distance from the solute com to assume bulk water
bulk_flag = 1                                       !default - calculate diffusion for bulk water
kr_start = 1                                        !environments loop starter - if bulk_flag = 0; kr_start = 2 

if(dr_bulk==0.0)then
   bulk_flag = 0
   write(*,*)
   write(*,*)'Warning!!! No bulk region defined for diffusion calculation'
   write(*,*)'Bulk diffusion will not be calculated accordingly'
   write(*,*)
   write(n0,*)
   write(n0,*)'Warning!!! No bulk region defined for diffusion calculation'
   write(n0,*)'Bulk diffusion will not be calculated accordingly'
   write(n0,*)
endif   

NRSHMAX = 4                                      !maximum number of environments
if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/difsampw                                  !for memory allocation purposes

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))  
NTDW4O = 0; NONTDW = 0 

!mean square displacement
allocate(nH2Oorg(norgmax,NRSHMAX))
allocate(msdx_b(NDELS),  msdy_b(NDELS),  msdz_b(NDELS))
allocate(msdx_hs(NDELS), msdy_hs(NDELS), msdz_hs(NDELS))
allocate(msdx_th(NDELS) ,msdy_th(NDELS), msdz_th(NDELS))
allocate(msdx_nth(NDELS),msdy_nth(NDELS),msdz_nth(NDELS))
allocate(msd_mean(NDELS,NRSHMAX),msd_mean_x(NDELS,NRSHMAX),msd_mean_y(NDELS,NRSHMAX),msd_mean_z(NDELS,NRSHMAX))
msdx_b=0.0d0;   msdy_b=0.0d0;   msdz_b=0.0d0
msdx_hs=0.0d0;  msdy_hs=0.0d0;  msdz_hs=0.0d0
msdx_th=0.0d0;  msdy_th=0.0d0;  msdz_th=0.0d0
msdx_nth=0.0d0; msdy_nth=0.0d0; msdz_nth=0.0d0
msd_mean=0.0d0;msd_mean_x=0.0d0;msd_mean_y=0.0d0;msd_mean_z=0.0d0


!allocate(RXt0(norgmax,nmolwat,NRSHMAX),RYt0(norgmax,nmolwat,NRSHMAX),RZt0(norgmax,nmolwat,NRSHMAX))
!allocate(RXt(norgmax,nmolwat,NRSHMAX),RYt(norgmax,nmolwat,NRSHMAX),RZt(norgmax,nmolwat,NRSHMAX))
!allocate(NH2OidO(norgmax,nmolwat,NRSHMAX))

allocate(RXt0(norgmax,maxwat_solv,NRSHMAX),RYt0(norgmax,maxwat_solv,NRSHMAX),RZt0(norgmax,maxwat_solv,NRSHMAX))
allocate(RXt(norgmax,maxwat_solv,NRSHMAX),RYt(norgmax,maxwat_solv,NRSHMAX),RZt(norgmax,maxwat_solv,NRSHMAX))
allocate(NH2OidO(norgmax,maxwat_solv,NRSHMAX))

allocate(ntime(norgmax,NRSHMAX))
ntime = 1; nH2Oorg = 0
allocate(kv(NRSHMAX))
allocate(kwlabel(NRSHMAX,nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled) 
allocate(nH2Omean(NRSHMAX))
allocate(koc(NRSHMAX))
allocate(k_org(NRSHMAX,NDELS))
nH2Omean = 0
kwlabel = 0
koc = 0
k_org = 0

 cmxsol = 0.0d0; cmysol=0.0d0; cmzsol=0.0d0

write(*,*)'Solute atoms and HSh onset/cut-off (Ang)'
write(n0,*)'Solute atoms and HSh onset/cut-off (Ang)'
do i = 1,natHSh
   write(*,'(i6,1x,f7.2,a1,f7.2)')natsol_id(i),oatsol_hs(i),'/',ratsol_hs(i)
   write(n0,'(i6,1x,f7.2,a1,f7.2)')natsol_id(i),oatsol_hs(i),'/',ratsol_hs(i)
end do
write(*,*)
write(n0,*)

! Start msd calculation
write(*,*)'Starting water environment diffusion calculation...'
write(*,*)'Routine wat_env_diff_stat'
write(*,*)'Results are printed to MSD_env_out_stat'
write(*,*)'water diffusion delay time-window (ps) =',msd_delay
write(*,*)'msd water sampling frequency (steps) = ',difsampw
write(*,*)'total number of time-steps =',nstep
write(*,*)'msd delay-time (steps) =',NDELS
write(*,*)
write(n0,*)'Starting water environment diffusion calculation...'
write(n0,*)'Routine wat_env_diff_stat'
write(n0,*)'Results are printed to MSD_env_out_stat'
write(n0,*)'water diffusion delay time-window (ps) =',msd_delay
write(n0,*)'msd water sampling frequency (steps) = ',difsampw
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'msd delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
!check      if(j==1)write(*,'(i5,1x,A4,1x,3(F14.5))')natmid(i),atomname(i),x(i),y(i),z(i)
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(ninput)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(ninput)x(i),y(i),z(i)
!NG         read(ninput)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules in specific environments
!NG   if((j==1.or.mod(j,difsampw)==0).and.(j.le.nstep-NDELS))then               !time origin: sample waters in specific environments
      if(j==1.or.mod(j,difsampw)==0)then                                  !time origin: sample waters in specific environments
      NTDW4O = 0; NONTDW = 0
      write(*,*)'time-step',j,' sampling water molecules'
      kv = 0                                                                 !waters in each environment 
      ko = ko + 1                                                            !number of origins = number of blocks used to compute the tcfs mean 
      kwlabel = 0
   
!Find solute centre of mass for bulk water tcfs calculation     

      if(bulk_flag==1)call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
      io = 0
      do i=nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water      
         dx = x(i) - cmxsol
         dy = y(i) - cmysol
         dz = z(i) - cmzsol
!         dx = dx - anint(dx/cell(1))*cell(1)
!         dy = dy - anint(dy/cell(2))*cell(2)
!         dz = dz - anint(dz/cell(3))*cell(3)
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
!check         write(*,*)'distance to the centre of mass = ',dr
         if(bulk_flag==0)dr_bulk = 100.d0*cell(1)                            !define a distance larger than the box "radius"
         if(dr>=dr_bulk)then
            if(bulk_flag==0)write(*,*)'Warning!!! Found a water molecule in the bulk'
!bulk msd(t0)             
            nradsh = 1                                       !bulk environment
            kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
            kvi = kv(nradsh)                                 !O number
!check            write(*,*)'bulk',j,i,kvi
            nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
            NH2OidO(ko,kvi,nradsh) = i                       ! Oxygen id
            RXt0(ko,kvi,nradsh) = x(i)                 !origin = ko, O number = kvi, environment type = nradsh
            RYt0(ko,kvi,nradsh) = y(i)
            RZt0(ko,kvi,nradsh) = z(i)     
            if(kvi >= maxwat_solv)then
               write(*,*)'error - too many waters in the Bulk'
               write(*,*)'action - increase maxwat_solv and re-compile'
               write(*,*)'or avoid bulk calculation - make cut-off = 0 '
               stop
            endif   
!bulk msd(t0) end                     
         endif 
         
!Hydration shell - tetrahedral and non tetrahedral waters in the first HSh         
         do at_id=1,natHSh                                      !Loop over atomic species to analyze the environment         
            k = natsol_id(at_id)
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2.0d0 + dy**2.0d0 + dz**2.0d0
            dr  = dsqrt(dr2)
!Hydration shell msd(t0)            
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then       !water within hydration shell radius            
!Avoid water molecules repetition
               nradsh = 2                                               !hydration shell environment
               if(kwlabel(nradsh,io)==0)then
                  kwlabel(nradsh,io)=1                                  !oxygen id for correction of the nb of waters
               else
                  goto 10                                               !sample another water molecule
               endif   
!End avoid water molecules repitition 
!NG               nradsh = 2                                     !hydration shell environment
               kv(nradsh) = kv(nradsh) + 1                       !number of water OH groups (index)
               kvi = kv(nradsh)                                  !OH number
!check               write(*,*)'HSh',j,i,iH,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1       !number of OH groups in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                        ! Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                        !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)     
               if(kvi >= maxwat_solv)then
                  write(*,*)'error - too many waters in the HSh'
                  write(*,*)'action - increase maxwat_solv and re-compile'
                  stop
               endif   
!Hydration shell msd(t0) end                              

!Find tetrahedral and non-tetrahedral waters in the hydration shell for sub-ensemble tetrahedral and non-tetrahedral msd        
!FIND THE FOUR NEAREST O ATOMS J TO EACH OXYGEN I
               RADMAX  = 25.0D0
               RADMIN  =  0.0D0        
               JT      = 1
               DO WHILE(JT.LE.4) 
                  do kw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
                     if(kw.ne.i)then
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOMS i AND OXYGEN ATOMS kw  
                        dx = x(i) - x(kw)
                        dy = y(i) - y(kw)
                        dz = z(i) - z(kw)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        dr2 = dx**2 + dy**2 + dz**2
                        dr  = dsqrt(dr2)
                        IF((dr>RADMIN).AND.(dr<RADMAX))THEN
                           RADMAX = dr
!NTDW4O STORES THE O ATOMS NEIGHBORS OF i: JT = 1 FIRST NEIGHBOR; JT = 2 SECOND NEIGHBOR; JT = 3 THIRD NEIGHBOR; JT = 4 FOURTH NEIGHBOR
                           NTDW4O(JT) = kw
                        ENDIF
                     endif
                  end do
!check                  write(*,*)j,i,io,JT,NTDW4O(JT)
                  JT = JT + 1
                  RADMIN = RADMAX
                  RADMAX = 25.0D0
               END DO
!CHECK IF SOLUTE ATOMS ARE CLOSER THAN THE 4th NEAREST O ATOM          
               KI=NTDW4O(4)                                       !Oxygen ID of the 4th nearest water neighbor 
               dx = x(i) - x(KI)
               dy = y(i) - y(KI)
               dz = z(i) - z(KI)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               dr4nb = dr                                           !distance to the fourth nearest neighbor
!EXCLUDE WATER MOLECULES AS TETRAHEDRON ORIGINS IF ANY SOLUTE ATOM IS CLOSER THAN THE 4th VERTEX WATER OXYGEN
               L = 1                           !run over all solute atoms - heavy and non-heavy atoms
               DO WHILE(L.LE.nmolsol*natmsol)
                  dx = x(i)-x(L)
                  dy = y(i)-y(L)
                  dz = z(i)-z(L)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  IF(dr.LE.dr4nb)THEN             
                     NONTDW(io) = 1                               !water i is non-tetrahedral
                     goto 5
                  ENDIF
                  L = L+1
               END DO     
               5 continue           
               if(NONTDW(io)==1)then                              !Non-tetrahedral water               
!Hydration shell - non tetrahedral msd(t0)             
                  nradsh = 4                                      !hydration shell (non-tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                     !number of oxygens (index)
                  kvi = kv(nradsh)                                !O number
!check                  write(*,*)'non-tetrah',j,i,kvi
                  nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1     !number of waters in environment nradsh at origin ko
                  NH2OidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  RXt0(ko,kvi,nradsh) = x(i)                      !origin = ko, O number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = y(i)
                  RZt0(ko,kvi,nradsh) = z(i)     
!Hydration shell - non tetrahedral msd(t0) end                                           
               elseif(NONTDW(io)==0)then                          !Tetrahedral water 
!Hydration shell - tetrahedral msd(t0)             
                  nradsh = 3                                      !hydration shell (tetrahedral) environment
                  kv(nradsh) = kv(nradsh) + 1                     !number of oxygens (index)
                  kvi = kv(nradsh)                                !O number
!check                  write(*,*)'tetrah',j,i,kvi
                  nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1     !number of oxygens in environment nradsh at origin ko
                  NH2OidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  RXt0(ko,kvi,nradsh) = x(i)                      !origin = ko, O number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = y(i)
                  RZt0(ko,kvi,nradsh) = z(i)             
!Hydration shell - tetrahedral msd(t0) end  
               endif
!            
            endif                                                    !end water is in the coordination sphere      
         end do                                                      !end loop over atomic species
      10 continue   
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nH2Oorg(ko,1),nH2Oorg(ko,2),nH2Oorg(ko,3),nH2Oorg(ko,4)
   endif                                                             !end sampling new time-origin      
   
!calculate msd - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nH2Oorg(kt,kr)>0.and.ntime(kt,kr)<=NDELS)then  !ntime counts the delay times for which the tcf has already been calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nH2Oorg(kt,kr)                                      !loop over water molecules in environment kr: nH2Oorg can be zero for a given time-origin 
               i=NH2OidO(kt,iw,kr)
               RXt(kt,iw,kr) = x(i)
               RYt(kt,iw,kr) = y(i)
               RZt(kt,iw,kr) = z(i)
               nt = ntime(kt,kr)
               if(kr==1)then                                                                                 !bulk
                  msdx_b(nt) = msdx_b(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_b(nt) = msdy_b(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_b(nt) = msdz_b(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==2)then                                                                             !hshell
                  msdx_hs(nt) = msdx_hs(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_hs(nt) = msdy_hs(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_hs(nt) = msdz_hs(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==3)then                                                                             !hshell-tetrahedral
                  msdx_th(nt) = msdx_th(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_th(nt) = msdy_th(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_th(nt) = msdz_th(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==4)then                                                                             !hShell-non-tetrahedral
                  msdx_nth(nt) = msdx_nth(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_nth(nt) = msdy_nth(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_nth(nt) = msdz_nth(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               endif
            end do   
            no = ntime(kt,kr)                                            
            k_org(kr,no) = k_org(kr,no) + 1                              !origins at a given delay-time
            ntime(kt,kr) = ntime(kt,kr) + 1                              !delay time counter [initialized at 1]
         endif
      end do                                                             !end loop over environments
   end do                                                                !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate average msd over different origins
do kr = 1,NRSHMAX                                      !loop over environments
   do it = 1,NDELS
      if(kr==1)then
         msd_mean_x(it,kr) = msdx_b(it)
         msd_mean_y(it,kr) = msdy_b(it)
         msd_mean_z(it,kr) = msdz_b(it)
         msd_mean(it,kr)   = (msdx_b(it)+msdy_b(it)+msdz_b(it))      
      elseif(kr==2)then
         msd_mean_x(it,kr) = msdx_hs(it)
         msd_mean_y(it,kr) = msdy_hs(it)
         msd_mean_z(it,kr) = msdz_hs(it)
         msd_mean(it,kr)   = (msdx_hs(it)+msdy_hs(it)+msdz_hs(it))
      elseif(kr==3)then
         msd_mean_x(it,kr) = msdx_th(it)
         msd_mean_y(it,kr) = msdy_th(it)
         msd_mean_z(it,kr) = msdz_th(it)
         msd_mean(it,kr)   = (msdx_th(it)+msdy_th(it)+msdz_th(it))
      elseif(kr==4)then
         msd_mean_x(it,kr) = msdx_nth(it)
         msd_mean_y(it,kr) = msdy_nth(it)
         msd_mean_z(it,kr) = msdz_nth(it)
         msd_mean(it,kr)   = (msdx_nth(it)+msdy_nth(it)+msdz_nth(it))
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of O atoms on each environment, kr, over the number of origins, ko 
!nH2Oorg(kt,kr) - number of waters found in origin (block) kt and in the environment kr               
do kr = 1,NRSHMAX
   do kt = 1, ko
      nH2Omean(kr) = nH2Omean(kr) + nH2Oorg(kt,kr)
      if(nH2Oorg(kt,kr)>0)then
         koc(kr) = koc(kr) + 1                                         !number of origins where at least a single water was found in a given environment
      endif
   end do
end do

if(bulk_flag==0)kr_start=2                                        !do not print the bulk

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a45,1x,i7)')'#Number of origins where waters were found = ',koc(kr)
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a45,1x,i7)')'#Number of origins where waters were found = ',koc(kr)
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko) 
end do   

deallocate(msdx_b,  msdy_b,  msdz_b)
deallocate(msdx_hs, msdy_hs, msdz_hs)
deallocate(msdx_th ,msdy_th, msdz_th)


!Print average msd     
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==1)write(n1,'(a5)')'#Bulk'
   if(kr==2)write(n2,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==2)write(n2,'(a4)')'#HSh'
   if(kr==3)write(n3,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==3)write(n3,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n4,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==4)write(n4,'(a15)')'#HSh-non-tetrah'
   
   if(kr==1)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==1)write(n5,'(a5)')'#Bulk'
   if(kr==2)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==2)write(n5,'(a4)')'#HSh'
   if(kr==3)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==3)write(n5,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
   if(kr==4)write(n5,'(a15)')'#HSh-non-tetrah'
   do it = 1,NDELS
      jt = it -1                                                     !print msd starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)

!check      write(*,*)kr,it,k_org(kr,it)
      
!if waters were found in every environment for every origin and the number of origins for each "it" is the same       

!NG      if(kr==1)WRITE(n1,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==2)WRITE(n2,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==3)WRITE(n3,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==4)WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG
!NG      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==2)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==3)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
!NG      if(kr==4)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)


!if for some origins no waters were found on a given environment but the number of origins for each "it" is the same   

!NG      if(kr==1)WRITE(n1,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      if(kr==2)WRITE(n2,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      if(kr==3)WRITE(n3,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      if(kr==4)WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      
!NG      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      if(kr==2)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      if(kr==3)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))
!NG      if(kr==4)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(koc(kr))

!if for some origins no waters were found on a given environment and the number of origins for each "it" is not the same
!this allows to sample origins for which the msd cannot be calculated until NDELS for every origin - statistics decreases with time

      if(kr==1)WRITE(n1,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      if(kr==2)WRITE(n2,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      if(kr==3)WRITE(n3,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      if(kr==4)WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      
      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      if(kr==2)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      if(kr==3)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))
      if(kr==4)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(k_org(kr,it))      
      
   end do
   if(kr==1)write(n1,*)
   if(kr==2)write(n2,*)
   if(kr==3)write(n3,*)
   if(kr==4)write(n4,*)
!   
   if(kr==1)write(n5,*)
   if(kr==2)write(n5,*)
   if(kr==3)write(n5,*)
   if(kr==4)write(n5,*)
end do

close(n1) 
close(n2) 
close(n3) 
close(n4)
close(n5)

!Calculate the Einstein self-diffusion coefficient - MSD Linear Squares fit
write(n0,*)
write(n0,*)'Einstein self-diffusion coefficient'
write(n0,*)
do kr = kr_start,NRSHMAX
!read data points in msd_LSF_Einstein for the LSF
   if(kr==1)then
      nf = n1
      write(n0,'(a5)')'#Bulk'
   elseif(kr==2)then
      nf = n2
      write(n0,'(a4)')'#HSh'
   elseif(kr==3)then
      nf = n3
      write(n0,'(a11)')'#HSh-tetrah'
   elseif(kr==4)then
      nf = n4
      write(n0,'(a15)')'#HSh-non-tetrah'
   endif  
   
   call msd_LSF_Einstein(nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end,a,b,d,r,n)
           
   if(kr==1)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   elseif(kr==2)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   elseif(kr==3)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   elseif(kr==4)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
      write(*,*)
   endif   
end do 

close(n0) 

deallocate(NTDW4O,NONTDW)
deallocate(nH2Oorg,NH2OidO)
deallocate(ntime)
deallocate(msd_mean,msd_mean_x,msd_mean_y,msd_mean_z)
deallocate(kv,nH2Omean)
deallocate(koc)

deallocate(xyz)
deallocate(x,y,z)
  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
      return

END SUBROUTINE wat_env_diff_stat



! Pure Water routines

SUBROUTINE pure_wat_orient_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay,m_time)  

! Calculate water's OH reorientational tcfs - 1st, 2nd, and 3rd Legendre Polynomials - Pure Water
! Every m_time time-steps the code samples every water in the system 
    integer,intent(in)                               :: npurew,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: mtdelay     ! reorientational tcf delay time (time-windows/ps)
    integer,intent(in)                               :: m_time      ! time freq. to sample water OH groups in radial shells for reor. tcf calculation
   
! Local variables
!NG    real(kind=8)                               :: x(natms),y(natms),z(natms) 
!    real(kind=4)                               :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG    real(kind=4)                               :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                               :: cell(3)  
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    
    real(kind=8),dimension(:,:,:),allocatable  :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable  :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                               :: TSTEP,PSECS
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    real(kind=8)                               :: xoh,yoh,zoh,sqoh
    integer                                    :: NDELS,KINTVL,NDELS1
    integer,dimension(:,:,:),allocatable       :: NOHidO,NOHidH
    integer,dimension(:,:),allocatable         :: ntime
    integer,dimension(:,:),allocatable         :: nOHorg
    integer,dimension(:),allocatable           :: kv
    integer      :: n0, n1, n2, n3, n11, n12, n13
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt, jw
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer                                    :: norgmax
    integer                                    :: NRSHMAX   
    integer,dimension(:),allocatable           :: norg_0HB  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!Orientational tcf
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b         !bulk
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs        !hsh
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: tcf_mean_1   ,tcf_mean_2   ,tcf_mean_3   
    real(kind=8)                                 :: X0TLEG,ARGLEG1,ARGLEG2,ARGLEG3
    integer,dimension(:),allocatable             :: nOHmean
    integer                                      :: kr_start
    integer,dimension(:),allocatable             :: koc                  !Nb of origins corrected for origins where no waters were found in a given environment

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Leg_PW_out/log_PW_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Leg_PW_out')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
!
n11 = 210
n12 = 220
n13 = 230

open(n0,file='Leg_PW_out/log_PW_tcf.dat',status='unknown',action='write')
open(n1,file='Leg_PW_out/wat_reor_tcf_1.dat',status='unknown',action='write')
open(n2,file='Leg_PW_out/wat_reor_tcf_2.dat',status='unknown',action='write')
open(n3,file='Leg_PW_out/wat_reor_tcf_3.dat',status='unknown',action='write')
!
open(n11,file='Leg_PW_out/tcf_1_rot_tau.dat',status='unknown',action='write')
open(n12,file='Leg_PW_out/tcf_2_rot_tau.dat',status='unknown',action='write')
open(n13,file='Leg_PW_out/tcf_3_rot_tau.dat',status='unknown',action='write')    

!time-window for calculation of each tcf
! IDINT - convert to integer - truncate

NDELS = IDINT(mtdelay*1000d0/dt)                  !sets the delay time (in frames)

kr_start = 1                                      !environments loop starter - for pure water single environment 

write(*,*)'input delay-time (ps) =',mtdelay
write(*,*)'delay-time (input frames units) =',NDELS
write(n0,*)'input delay-time (ps) =',mtdelay
write(n0,*)'delay-time (input frames units) =',NDELS

!NRSHMAX = 4                                      !maximum number of environments
NRSHMAX = 1                                       !maximum number of environments - 1 for pure water
if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/m_time                                  !for memory allocation purposes

!orient. tcf
allocate(nOHorg(norgmax,NRSHMAX))
allocate(RXt0(norgmax,2*nmolwat,NRSHMAX),RYt0(norgmax,2*nmolwat,NRSHMAX),RZt0(norgmax,2*nmolwat,NRSHMAX))
allocate(RXt(norgmax,2*nmolwat,NRSHMAX),RYt(norgmax,2*nmolwat,NRSHMAX),RZt(norgmax,2*nmolwat,NRSHMAX))
allocate(NOHidO(norgmax,2*nmolwat,NRSHMAX),NOHidH(norgmax,2*nmolwat,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(tcf_leg_1_b(NDELS),tcf_leg_1_hs(NDELS),tcf_leg_1_th(NDELS),tcf_leg_1_nth(NDELS))
allocate(tcf_leg_2_b(NDELS),tcf_leg_2_hs(NDELS),tcf_leg_2_th(NDELS),tcf_leg_2_nth(NDELS))
allocate(tcf_leg_3_b(NDELS),tcf_leg_3_hs(NDELS),tcf_leg_3_th(NDELS),tcf_leg_3_nth(NDELS))
allocate(tcf_mean_1(NDELS,NRSHMAX) ,tcf_mean_2(NDELS,NRSHMAX)  ,tcf_mean_3(NDELS,NRSHMAX))
tcf_leg_1_b = 0.0d0;tcf_leg_1_hs = 0.0d0;tcf_leg_1_th = 0.0d0;tcf_leg_1_nth = 0.0d0
tcf_leg_2_b = 0.0d0;tcf_leg_2_hs = 0.0d0;tcf_leg_2_th = 0.0d0;tcf_leg_2_nth = 0.0d0
tcf_leg_3_b = 0.0d0;tcf_leg_3_hs = 0.0d0;tcf_leg_3_th = 0.0d0;tcf_leg_3_nth = 0.0d0
tcf_mean_1  = 0.0d0;tcf_mean_2   = 0.0d0;tcf_mean_3   = 0.0d0 
ntime = 1; nOHorg = 0

allocate(norg_0HB(NRSHMAX))
norg_0HB = 0
allocate(kv(NRSHMAX))
allocate(nOHmean(NRSHMAX))
allocate(koc(NRSHMAX))
nOHmean = 0
koc = 0

write(*,*)
write(*,*)'Pure Water reorientation tcfs calculation'
write(*,*)'Routine pure_wat_orient_tcfs'
write(*,*)'water tcfs delay time-window (ps) =',mtdelay
write(*,*)'Results are printed to Leg_PW_out'
write(*,*)
write(n0,*)
write(n0,*)'Pure Water reorientation tcfs calculation'
write(n0,*)'Routine pure_wat_orient_tcfs'
write(n0,*)'water tcfs delay time-window (ps) =',mtdelay
write(n0,*)'Results are printed to Leg_PW_out'
write(n0,*)

! Start tcfs calculation
write(*,*)
write(*,*)'Starting reor. tcf calculation...'
write(*,*)'tcfs OH sampling frequency (steps) = ',m_time
write(*,*)'total number of time-steps =',nstep
write(*,*)'tcfs delay-time (steps) =',NDELS
write(*,*)
write(n0,*)
write(n0,*)'Starting reor. tcf calculation...'
write(n0,*)'tcfs OH sampling frequency (steps) = ',m_time
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'tcfs delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif      
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!End Version 9      
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules in specific environments - pure water - single environment
   if((j==1.or.mod(j,m_time)==0).and.(j.le.nstep-NDELS))then               !time origin: sample OH groups in specific environments
      write(*,*)'time-step',j,' sampling for OH groups'
      kv = 0                                                               !OH groups in each environment 
      ko = ko + 1                                                          !number of origins = number of blocks used to compute the tcfs mean and stdev
      io = 0
      do i=1,natms,nwatsites                              !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...

!bulk tcf(t0)             
         do iH = 1,2
            nradsh = 1                                     !bulk environment
            kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
            kvi = kv(nradsh)                               !OH number
!check            write(*,*)j,i,iH,kvi
            nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
            NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
            NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
            dx = x(i+iH) - x(i)
            dy = y(i+iH) - y(i)
            dz = z(i+iH) - z(i)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            xoh = dx
            yoh = dy
            zoh = dz
            sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
            RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
            RYt0(ko,kvi,nradsh) = yoh/sqoh
            RZt0(ko,kvi,nradsh) = zoh/sqoh     
         end do         
!bulk tcf(t0) end                     
         
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nOHorg(ko,1)
   endif                                                             !end sampling new time-origin      
   
!calculate tcfs - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nOHorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then                !ntime counts the delay times for which the tcf has already been calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nOHorg(kt,kr)                                       !loop over OH groups in environment kr: nOHorg can be zero for a given time-origin
                  
               i=NOHidO(kt,iw,kr)
               ihat = NOHidH(kt,iw,kr)                                    !Hydrogen atom (ihat = 1 or 2)        
               dx = x(i+ihat) - x(i)
               dy = y(i+ihat) - y(i)
               dz = z(i+ihat) - z(i)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               xoh = dx
               yoh = dy
               zoh = dz
               sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
               RXt(kt,iw,kr) = xoh/sqoh
               RYt(kt,iw,kr) = yoh/sqoh
               RZt(kt,iw,kr) = zoh/sqoh
               X0TLEG = RXt0(kt,iw,kr)*RXt(kt,iw,kr)+RYt0(kt,iw,kr)*RYt(kt,iw,kr)+RZt0(kt,iw,kr)*RZt(kt,iw,kr)
               ARGLEG1= X0TLEG
               ARGLEG2= 0.5d0*(3.0d0*X0TLEG**2.0-1.0d0)
               ARGLEG3= 0.5d0*(5.0d0*X0TLEG**3.0-3.0d0*X0TLEG)
               nt = ntime(kt,kr)
               if(kr==1)then                                            !bulk
                  tcf_leg_1_b(nt)=tcf_leg_1_b(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_b(nt)=tcf_leg_2_b(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_b(nt)=tcf_leg_3_b(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==2)then                                        !hshell
                  tcf_leg_1_hs(nt)=tcf_leg_1_hs(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_hs(nt)=tcf_leg_2_hs(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_hs(nt)=tcf_leg_3_hs(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==3)then                                        !hshell-tetrahedral
                  tcf_leg_1_th(nt)=tcf_leg_1_th(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_th(nt)=tcf_leg_2_th(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_th(nt)=tcf_leg_3_th(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==4)then                                        !hShell-non-tetrahedral
                  tcf_leg_1_nth(nt)=tcf_leg_1_nth(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_nth(nt)=tcf_leg_2_nth(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_nth(nt)=tcf_leg_3_nth(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean tcfs over different origins
do kr = 1,NRSHMAX         !loop over environments
   do it = 1,NDELS
!nOHorg(kt,kr) - number of OH groups found in origin (block) kt and in the environment kr               
      if(kr==1)then
         tcf_mean_1(it,kr) = tcf_leg_1_b(it)
         tcf_mean_2(it,kr) = tcf_leg_2_b(it)
         tcf_mean_3(it,kr) = tcf_leg_3_b(it)
      elseif(kr==2)then   
         tcf_mean_1(it,kr) = tcf_leg_1_hs(it)
         tcf_mean_2(it,kr) = tcf_leg_2_hs(it)
         tcf_mean_3(it,kr) = tcf_leg_3_hs(it)
      elseif(kr==3)then
         tcf_mean_1(it,kr) = tcf_leg_1_th(it)
         tcf_mean_2(it,kr) = tcf_leg_2_th(it)
         tcf_mean_3(it,kr) = tcf_leg_3_th(it)
      elseif(kr==4)then
         tcf_mean_1(it,kr) = tcf_leg_1_nth(it)
         tcf_mean_2(it,kr) = tcf_leg_2_nth(it)
         tcf_mean_3(it,kr) = tcf_leg_3_nth(it)
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of OH groups on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nOHmean(kr) = nOHmean(kr) + nOHorg(kt,kr)
   end do
end do

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a28,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko) 
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a28,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko)
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2) 
end do   
      
deallocate(tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b)
deallocate(tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs)
deallocate(tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th)
deallocate(tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth)
      
!Print mean orient. tcfs      
write(n2,9)
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a5)')'#Bulk'
   if(kr==2)write(n1,'(a4)')'#HSh'
   if(kr==3)write(n1,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n1,'(a15)')'#HSh-non-tetrah'
   if(kr==1)write(n2,'(a5)')'#Bulk'
   if(kr==2)write(n2,'(a4)')'#HSh'
   if(kr==3)write(n2,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n2,'(a15)')'#HSh-non-tetrah'
   if(kr==1)write(n3,'(a5)')'#Bulk'
   if(kr==2)write(n3,'(a4)')'#HSh'
   if(kr==3)write(n3,'(a11)')'#HSh-tetrah'
   if(kr==4)write(n3,'(a15)')'#HSh-non-tetrah'
   do it = 1,NDELS
      jt = it -1                                                     !print tcf starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      WRITE(n1,19)PSECS,tcf_mean_1(it,kr)/dfloat(ko)     !mean orientational tcf rank 1 over different time-origins 
      WRITE(n2,19)PSECS,tcf_mean_2(it,kr)/dfloat(ko)     !mean orientational tcf rank 2 over different time-origins   
      WRITE(n3,19)PSECS,tcf_mean_3(it,kr)/dfloat(ko)     !mean orientational tcf rank 3 over different time-origins 
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
end do

!Calculate the relaxation times - tau - trapezoid method
do kr = kr_start,NRSHMAX
   if(kr==1)then
      write(n11,'(a5)')'#Bulk'
      write(n12,'(a5)')'#Bulk'
      write(n13,'(a5)')'#Bulk'
   elseif(kr==2)then
      write(n11,'(a4)')'#HSh'
      write(n12,'(a4)')'#HSh'
      write(n13,'(a4)')'#HSh'
   elseif(kr==3)then
      write(n11,'(a11)')'#HSh-tetrah'
      write(n12,'(a11)')'#HSh-tetrah'
      write(n13,'(a11)')'#HSh-tetrah'
   elseif(kr==4)then
      write(n11,'(a15)')'#HSh-non-tetrah'
      write(n12,'(a15)')'#HSh-non-tetrah'
      write(n13,'(a15)')'#HSh-non-tetrah'
   endif  
   koc(kr)=ko
!   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
         
   if(kr==1)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==2)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==3)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==4)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   endif   
end do 

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n11)
close(n12)
close(n13)

deallocate(nOHorg,NOHidO,NOHidH)
deallocate(ntime)
deallocate(tcf_mean_1,tcf_mean_2,tcf_mean_3)
deallocate(norg_0HB,kv,nOHmean)
deallocate(koc)

deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT('#',17X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
   return

END SUBROUTINE pure_wat_orient_tcfs



SUBROUTINE pure_wat_LDL_HDL_otcf(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay,m_time,Nenv,kLDL_HDL,&
                                 LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2)  
! Lars Pettersson and Gaia Camisasca collaboration - June 2018 - coded in version 24
! The program reads a list from an input file with information on LDL and HDL waters
! Calculate OH reorientational tcfs - 1st, 2nd, and 3rd Legendre Polynomials - pure water LDL/HDL ensembles
! Every m_time time-steps the code samples water molecules in the system that are either LDL, HDL or neither - output 3 orientational tcfs 

    integer,intent(in)                               :: npurew,nbox,Nenv,kLDL_HDL
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: mtdelay     ! reorientational tcf delay time (time-windows/ps)
    integer,intent(in)                               :: m_time      ! time freq. to sample water OH groups in radial shells for reor. tcf calculation
    real(kind=8),intent(in)                          :: LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2
   
! Local variables
!NG    real(kind=8)                               :: x(natms),y(natms),z(natms) 
!    real(kind=4)                               :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG    real(kind=4)                               :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                               :: cell(3)  
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    
    real(kind=8),dimension(:,:,:),allocatable  :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable  :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                               :: TSTEP,PSECS
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    real(kind=8)                               :: xoh,yoh,zoh,sqoh
    integer                                    :: NDELS,KINTVL,NDELS1
    integer,dimension(:,:,:),allocatable       :: NOHidO,NOHidH
    integer,dimension(:,:),allocatable         :: ntime
    integer,dimension(:,:),allocatable         :: nOHorg
    integer,dimension(:),allocatable           :: kv
    integer      :: n0, n1, n2, n3, n11, n12, n13, n15, n16, n20, n25
    integer      :: ihb,jhb
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt, jw
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer                                    :: norgmax
    integer                                    :: NRSHMAX   
    integer,dimension(:),allocatable           :: norg_0HB  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!Orientational tcf
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b         !bulk
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs        !LDL
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th        !HDL
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth       !OTHER
    real(kind=8),dimension(:,:),allocatable      :: tcf_mean_1   ,tcf_mean_2   ,tcf_mean_3   
    real(kind=8)                                 :: X0TLEG,ARGLEG1,ARGLEG2,ARGLEG3
    integer,dimension(:),allocatable             :: nOHmean
    integer                                      :: kr_start
    integer,dimension(:),allocatable             :: koc                  !Nb of origins corrected for origins where no waters were found in a given environment
!Version 24    
    integer                                      :: nw_list
    integer,dimension(:),allocatable             :: list_s
    integer,dimension(:),allocatable             :: list_LSI               !LSI index - LSI computed locally to compare with Gaia List
    integer                                      :: pwnmolsol,pwnatmsol,pwnions                !solution parameters - set to 0
    real(kind=8),dimension(:),allocatable        :: LSI
    integer,dimension(:),allocatable             :: NLSI
    real(kind=8)                                 :: r_lsi
    real(kind=8)                                 :: LSI_del,AVLSI,LSI_ang2,DISTR_LSI
    integer                                      :: nsh,LSI_NDELS,NLSI_count,nkr
!    integer                                      :: nsh_
    real(kind=8)                                 :: R_rdf_OO,RDEL,cubeh_,cube_,RADIUS,PI,VOLSHL,DSTNUM,FNOM
    integer                                      :: NSHLOO,NRDELS
    integer,dimension(:),allocatable             :: NRDFOO,NRDFOO_LDL,NRDFOO_HDL
    real(kind=8),dimension(:),allocatable        :: GROO,GROO_LDL,GROO_HDL
    

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!LSI values
r_lsi  = 3.70d0                                 !original value of the LSI

LSI_del = 0.0005d0                              !bin for LSI distribution
LSI_NDELS= 1.0d0/LSI_del

RDEL    = 0.05d0                                !bin size in rdf/distance distributions Ang
NRDELS  = 100.0d0/RDEL

PI=3.14159265d0
FNOM = dfloat(nmolwat)

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Leg_LDL_HDL/log_PW_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Leg_LDL_HDL')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
!
n11 = 210
n12 = 220
n13 = 230
!
n15 = 250
n16 = 260
!
n20 = 300

!check   n25 = 500

open(n0,file='Leg_LDL_HDL/log_PW_tcf.dat',status='unknown',action='write')
open(n1,file='Leg_LDL_HDL/wat_reor_tcf_1.dat',status='unknown',action='write')
open(n2,file='Leg_LDL_HDL/wat_reor_tcf_2.dat',status='unknown',action='write')
open(n3,file='Leg_LDL_HDL/wat_reor_tcf_3.dat',status='unknown',action='write')
!
open(n11,file='Leg_LDL_HDL/tcf_1_rot_tau.dat',status='unknown',action='write')
open(n12,file='Leg_LDL_HDL/tcf_2_rot_tau.dat',status='unknown',action='write')
open(n13,file='Leg_LDL_HDL/tcf_3_rot_tau.dat',status='unknown',action='write')

open(n15,file='Leg_LDL_HDL/LSI_dist.dat',status='unknown',action='write')
open(n16,file='Leg_LDL_HDL/rdf_OO.dat',status='unknown',action='write')

!Version 24 - read LDL/HDL water molecules label from input file
!Lars Pettersson - list
! 1 is LDL 
! 2 is HDL
! 3 skip

!kLDL_HDL == 0 read LDL/HDL from input file; kLDL_HDL == 1 calculate LSI on the fly
if(Nenv>1.and.kLDL_HDL==0)open(n20,file='List_W.dat',status='old',action='read')

!check open(n25,file='Leg_LDL_HDL/check_List_W.dat',status='unknown',action='write')

!time-window for calculation of each tcf
! IDINT - convert to integer - truncate

NDELS = IDINT(mtdelay*1000d0/dt)                  !sets the delay time (in frames)

kr_start = 1                                      !environments loop starter - for pure water single environment 

write(*,*)'input delay-time (ps) =',mtdelay
write(*,*)'delay-time (input frames units) =',NDELS
write(n0,*)'input delay-time (ps) =',mtdelay
write(n0,*)'delay-time (input frames units) =',NDELS

if(Nenv>1)then
   write(*,*)
   write(*,*)'Environments: 1 (Bulk), 2(LDL), 3(HDL), 4(other)'
   write(*,*)
   write(n0,*)
   write(n0,*)'Environments: 1 (Bulk), 2(LDL), 3(HDL), 4(other)'
   write(n0,*)
endif

if(kLDL_HDL==1)write(*,*)'LDL/HDL populations calculated on the fly'
if(kLDL_HDL==0)write(*,*)'LDL/HDL populations read from input file'
if(kLDL_HDL==1)then
   write(*,*)'LSI range for LDL = ',LSI_min_1,'-',LSI_max_1
   write(*,*)'LSI range for HDL = ',LSI_min_2,'-',LSI_max_2
endif
if(kLDL_HDL==1)write(n0,*)'LDL/HDL populations calculated on the fly'
if(kLDL_HDL==0)write(n0,*)'LDL/HDL populations read from input file'
if(kLDL_HDL==1)then
   write(n0,*)'LSI range for LDL = ',LSI_min_1,'-',LSI_max_1
   write(n0,*)'LSI range for HDL = ',LSI_min_2,'-',LSI_max_2
endif

NRSHMAX = Nenv                                   !maximum number of environments - 4 e.g. LDL, HDL, rest, all
!NRSHMAX = 1                                     !maximum number of environments - 1 for pure water

if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/m_time                                  !for memory allocation purposes

!orient. tcf
allocate(nOHorg(norgmax,NRSHMAX))
allocate(RXt0(norgmax,2*nmolwat,NRSHMAX),RYt0(norgmax,2*nmolwat,NRSHMAX),RZt0(norgmax,2*nmolwat,NRSHMAX))
allocate(RXt(norgmax,2*nmolwat,NRSHMAX),RYt(norgmax,2*nmolwat,NRSHMAX),RZt(norgmax,2*nmolwat,NRSHMAX))
allocate(NOHidO(norgmax,2*nmolwat,NRSHMAX),NOHidH(norgmax,2*nmolwat,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(tcf_leg_1_b(NDELS),tcf_leg_1_hs(NDELS),tcf_leg_1_th(NDELS),tcf_leg_1_nth(NDELS))
allocate(tcf_leg_2_b(NDELS),tcf_leg_2_hs(NDELS),tcf_leg_2_th(NDELS),tcf_leg_2_nth(NDELS))
allocate(tcf_leg_3_b(NDELS),tcf_leg_3_hs(NDELS),tcf_leg_3_th(NDELS),tcf_leg_3_nth(NDELS))
allocate(tcf_mean_1(NDELS,NRSHMAX) ,tcf_mean_2(NDELS,NRSHMAX)  ,tcf_mean_3(NDELS,NRSHMAX))
tcf_leg_1_b = 0.0d0;tcf_leg_1_hs = 0.0d0;tcf_leg_1_th = 0.0d0;tcf_leg_1_nth = 0.0d0
tcf_leg_2_b = 0.0d0;tcf_leg_2_hs = 0.0d0;tcf_leg_2_th = 0.0d0;tcf_leg_2_nth = 0.0d0
tcf_leg_3_b = 0.0d0;tcf_leg_3_hs = 0.0d0;tcf_leg_3_th = 0.0d0;tcf_leg_3_nth = 0.0d0
tcf_mean_1  = 0.0d0;tcf_mean_2   = 0.0d0;tcf_mean_3   = 0.0d0 
ntime = 1; nOHorg = 0

allocate(norg_0HB(NRSHMAX))
norg_0HB = 0
allocate(kv(NRSHMAX))
allocate(nOHmean(NRSHMAX))
allocate(koc(NRSHMAX))
nOHmean = 0
koc = 0
!Version 24 - List HDL/LDL
allocate(list_s(nmolwat))
allocate(list_LSI(nmolwat))
allocate(LSI(nmolwat))
allocate(NLSI(LSI_NDELS))

allocate(NRDFOO(NRDELS),NRDFOO_LDL(NRDELS),NRDFOO_HDL(NRDELS))
allocate(GROO(NRDELS),GROO_LDL(NRDELS),GROO_HDL(NRDELS))

LSI = 0
NLSI = 0
AVLSI = 0.0d0
pwnmolsol = 0
pwnatmsol = 0
pwnions   = 0  
NLSI_count = 0

!rdf
NRDFOO     = 0
NRDFOO_LDL = 0
NRDFOO_HDL = 0

write(*,*)
write(*,*)'Water reorientation tcfs calculation'
write(*,*)'Routine pure_wat_LDL_HDL_otcf'
write(*,*)'water tcfs delay time-window (ps) =',mtdelay
write(*,*)'Results are printed to Leg_LDL_HDL'
write(*,*)
write(n0,*)
write(n0,*)'Water reorientation tcfs calculation'
write(n0,*)'Routine pure_wat_LDL_HDL_otcf'
write(n0,*)'water tcfs delay time-window (ps) =',mtdelay
write(n0,*)'Results are printed to Leg_LDL_HDL'
write(n0,*)

! Start tcfs calculation
write(*,*)
write(*,*)'Starting reor. tcf calculation...'
write(*,*)'tcfs OH sampling frequency (steps) = ',m_time
write(*,*)'total number of time-steps =',nstep
write(*,*)'tcfs delay-time (steps) =',NDELS
write(*,*)
write(n0,*)
write(n0,*)'Starting reor. tcf calculation...'
write(n0,*)'tcfs OH sampling frequency (steps) = ',m_time
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'tcfs delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif      
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!End Version 9      
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   
    
    cubeh_ = cell(1)/2.0d0
    cube_  = cell(1)

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules in specific environments - pure water - single environment
   if((j==1.or.mod(j,m_time)==0).and.(j.le.nstep-NDELS))then               !time origin: sample OH groups in specific environments
      write(*,*)'time-step',j,' sampling for OH groups'
!      
!Version 24 - LDL/HDL water list           1(LDL) 2(HDL) 3(NOT LDL nor HDL)
      if(NRSHMAX > 1.and.kLDL_HDL==0)then
         do i=1,nmolwat
            read(n20,*)nw_list,list_s(i)
!check            write(n25,*)nw_list,list_s(i)
         end do
      endif       
  
      
!Version 26    ========================================================================================================
!Calculate the LSI
      if(kLDL_HDL==1)then              !calculate LDL/HDL populations on the fly
         io = 0
         do ihb = 1,natms,nwatsites
            LSI = 0
            io = io + 1
            call LSI_water(ihb,io,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,LSI,r_lsi)
!check LSI - compare with list_s(i) based on Gaia & Lars LSI cut-offs for LDL (> 0.17 Ang**2) and HDL (< 0.01 Ang**2) species            
            if(LSI(io) > LSI_min_1 .and. LSI(io) < LSI_max_1)then
               list_LSI(io) = 1               !LDL
            elseif(LSI(io) > LSI_min_2 .and. LSI(io) < LSI_max_2)then 
               list_LSI(io) = 2               !HDL
            else
               list_LSI(io) = 3               !other
            endif   
!Compute LSI distribution and average value
            nsh = LSI(io)/LSI_del+0.5D0            
!check            if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check            write(*,*)LSI(io),nsh
            if(nsh==0)then
               write(*,*)'Warning!!! LSI bin = 0; transformed to 1'
               write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
               nsh = 1
            endif   
!            nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!            write(*,*)nsh,nsh_
!            nsh = nsh_
            NLSI(nsh) = NLSI(nsh) + 1
            AVLSI = AVLSI + LSI(io)
            NLSI_count = NLSI_count + 1                         !count number of values used for distribution and average value calculation
!         end do
!calculate O-O rdf for LDL and HDL populations
            do jhb = 1,natms,nwatsites        !Loop over water oxygens
               IF(jhb.NE.ihb)THEN
!check         WRITE(*,*)ihb,jhb
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM ihb AND OXYGEN ATOMS jhb     
                  dx = x(ihb) - x(jhb)             
                  dy = y(ihb) - y(jhb)
                  dz = z(ihb) - z(jhb)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
                  R_rdf_OO = dr
!               ENDIF
!            end do       
!check            write(*,*)R_rdf_OO,cubeh_
                  if(R_rdf_OO < cubeh_)then
!bulk                  
                     NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                     NRDFOO(NSHLOO) = NRDFOO(NSHLOO)+1
                  endif   
                  if(R_rdf_OO < cubeh_.and.list_LSI(io) == 1)then
!LHDL             
                     NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                     NRDFOO_LDL(NSHLOO) = NRDFOO_LDL(NSHLOO)+1
                  endif   
                  if(R_rdf_OO < cubeh_.and.list_LSI(io) == 2)then
!HDL              
                     NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                     NRDFOO_HDL(NSHLOO) = NRDFOO_HDL(NSHLOO)+1
                  endif
               
               ENDIF
            
            end do !End jhb water oxygens loop
         
         end do    !End ihb water oxygens loop
                 
         list_s = list_LSI
      
      endif       !End calculate LDL/HDL populations on the fly
            
!End Version 26 ========================================================================================================                
            
      
      kv = 0                                                               !OH groups in each environment 
      ko = ko + 1                                                          !number of origins = number of blocks used to compute the tcfs mean and stdev
      io = 0                                              !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
      
      do i=1,natms,nwatsites                              !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...

!bulk tcf(t0) - all waters             
         do iH = 1,2
            nradsh = 1                                     !bulk environment
            kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
            kvi = kv(nradsh)                               !OH number
!check            write(*,*)j,i,iH,kvi
            nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
            NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
            NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
            dx = x(i+iH) - x(i)
            dy = y(i+iH) - y(i)
            dz = z(i+iH) - z(i)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            xoh = dx
            yoh = dy
            zoh = dz
            sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
            RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
            RYt0(ko,kvi,nradsh) = yoh/sqoh
            RZt0(ko,kvi,nradsh) = zoh/sqoh     
         end do         
!bulk tcf(t0) end

!LDL tcf(t0)
         if(NRSHMAX > 1)then
            if(list_s(io)==1)then
               do iH = 1,2
                  nradsh = 2                                     !LDL environment
                  kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                  kvi = kv(nradsh)                               !OH number
!check                  write(*,*)j,i,iH,kvi
                  nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                  NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  xoh = dx
                  yoh = dy
                  zoh = dz
                  sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                  RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = yoh/sqoh
                  RZt0(ko,kvi,nradsh) = zoh/sqoh     
               end do         
!LDL tcf(t0) end

!HDL tcf(t0)
            elseif(list_s(io)==2)then
               do iH = 1,2
                  nradsh = 3                                     !HDL environment
                  kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                  kvi = kv(nradsh)                               !OH number
!check                  write(*,*)j,i,iH,kvi
                  nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                  NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  xoh = dx
                  yoh = dy
                  zoh = dz
                  sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                  RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = yoh/sqoh
                  RZt0(ko,kvi,nradsh) = zoh/sqoh     
               end do         
!HDL tcf(t0) end

!NOT LDL/HDL tcf(t0)
            elseif(list_s(io)==3)then
               do iH = 1,2
                  nradsh = 4                                     !NOT LDL/HDL environment
                  kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                  kvi = kv(nradsh)                               !OH number
!check                  write(*,*)j,i,iH,kvi
                  nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                  NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  xoh = dx
                  yoh = dy
                  zoh = dz
                  sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                  RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = yoh/sqoh
                  RZt0(ko,kvi,nradsh) = zoh/sqoh     
               end do         
!NOT LDL/HDL tcf(t0) end
            endif           !end LDL or HDL
         endif              !end if NRSHMAX > 1
        
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nOHorg(ko,1)
   endif                                                             !end sampling new time-origin                                                                 
   
!calculate tcfs - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nOHorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then                !ntime counts the delay times for which the tcf has already been calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nOHorg(kt,kr)                                       !loop over OH groups in environment kr: nOHorg can be zero for a given time-origin
                  
               i=NOHidO(kt,iw,kr)
               ihat = NOHidH(kt,iw,kr)                                    !Hydrogen atom (ihat = 1 or 2)        
               dx = x(i+ihat) - x(i)
               dy = y(i+ihat) - y(i)
               dz = z(i+ihat) - z(i)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               xoh = dx
               yoh = dy
               zoh = dz
               sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
               RXt(kt,iw,kr) = xoh/sqoh
               RYt(kt,iw,kr) = yoh/sqoh
               RZt(kt,iw,kr) = zoh/sqoh
               X0TLEG = RXt0(kt,iw,kr)*RXt(kt,iw,kr)+RYt0(kt,iw,kr)*RYt(kt,iw,kr)+RZt0(kt,iw,kr)*RZt(kt,iw,kr)
               ARGLEG1= X0TLEG
               ARGLEG2= 0.5d0*(3.0d0*X0TLEG**2.0-1.0d0)
               ARGLEG3= 0.5d0*(5.0d0*X0TLEG**3.0-3.0d0*X0TLEG)
               nt = ntime(kt,kr)
               if(kr==1)then                                            !bulk
                  tcf_leg_1_b(nt)=tcf_leg_1_b(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_b(nt)=tcf_leg_2_b(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_b(nt)=tcf_leg_3_b(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==2)then                                        !LDL
                  tcf_leg_1_hs(nt)=tcf_leg_1_hs(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_hs(nt)=tcf_leg_2_hs(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_hs(nt)=tcf_leg_3_hs(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==3)then                                        !HDL
                  tcf_leg_1_th(nt)=tcf_leg_1_th(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_th(nt)=tcf_leg_2_th(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_th(nt)=tcf_leg_3_th(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==4)then                                        !OTHER - NOT LDL nor HDL
                  tcf_leg_1_nth(nt)=tcf_leg_1_nth(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_nth(nt)=tcf_leg_2_nth(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_nth(nt)=tcf_leg_3_nth(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean tcfs over different origins
do kr = 1,NRSHMAX         !loop over environments
   do it = 1,NDELS
!nOHorg(kt,kr) - number of OH groups found in origin (block) kt and in the environment kr               
      if(kr==1)then
         tcf_mean_1(it,kr) = tcf_leg_1_b(it)
         tcf_mean_2(it,kr) = tcf_leg_2_b(it)
         tcf_mean_3(it,kr) = tcf_leg_3_b(it)
      elseif(kr==2)then   
         tcf_mean_1(it,kr) = tcf_leg_1_hs(it)
         tcf_mean_2(it,kr) = tcf_leg_2_hs(it)
         tcf_mean_3(it,kr) = tcf_leg_3_hs(it)
      elseif(kr==3)then
         tcf_mean_1(it,kr) = tcf_leg_1_th(it)
         tcf_mean_2(it,kr) = tcf_leg_2_th(it)
         tcf_mean_3(it,kr) = tcf_leg_3_th(it)
      elseif(kr==4)then
         tcf_mean_1(it,kr) = tcf_leg_1_nth(it)
         tcf_mean_2(it,kr) = tcf_leg_2_nth(it)
         tcf_mean_3(it,kr) = tcf_leg_3_nth(it)
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of OH groups on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nOHmean(kr) = nOHmean(kr) + nOHorg(kt,kr)
   end do
end do

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a28,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko) 
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a28,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko)
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2) 
end do   
      
deallocate(tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b)
deallocate(tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs)
deallocate(tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th)
deallocate(tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth)

!================================================================================================Version 26 
!Print LSI distribution
if(kLDL_HDL==1)then
   do nkr = 1,LSI_NDELS
      if(NLSI(nkr) > 0)then
         LSI_ang2  = LSI_del*DFLOAT(nkr)
         DISTR_LSI = DFLOAT(NLSI(nkr))/(LSI_del*DFLOAT(NLSI_count))                             !Normalize
         WRITE(n15,39)LSI_ang2,DISTR_LSI
      endif
   end do  
   AVLSI = AVLSI/DFLOAT(NLSI_count)
   write(n0,*)
   write(n0,*)'LSI mean value (Ang) = ',AVLSI
   write(n0,*)

!Print rdf OO
!LOOP OVER RADIAL SHELLS - OXYGEN-OXYGEN RDF

!bulk
   WRITE(n16,'(a5)')'#bulk'
   DO kw=1,NRDELS   
      GROO(kw) = 0.0d0 
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!         GROO(kw)=DFLOAT(NRDFOO(kw))/(DSTNUM*DFLOAT(ko)*FNOM*VOLSHL)
         GROO(kw)=2.0d0*DFLOAT(NRDFOO(kw))/(DSTNUM*dfloat(nOHmean(1))*VOLSHL)
         WRITE(n16,39) RADIUS,GROO(kw)                         
      ENDIF
   END DO
   WRITE(n16,*)
   
!LDL    
   WRITE(n16,'(a4)')'#LDL'
   DO kw=1,NRDELS  
      GROO_LDL(kw) = 0.0d0 
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO_LDL(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
         GROO_LDL(kw)=2.0d0*DFLOAT(NRDFOO_LDL(kw))/(DSTNUM*dfloat(nOHmean(2))*VOLSHL)
         WRITE(n16,39) RADIUS,GROO_LDL(kw)                         
      ENDIF
   END DO
   WRITE(n16,*)

!HDL    
   WRITE(n16,'(a4)')'#HDL'
   DO kw=1,NRDELS  
      GROO_HDL(kw) = 0.0d0
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO_HDL(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
         GROO_HDL(kw)=2.0d0*DFLOAT(NRDFOO_HDL(kw))/(DSTNUM*dfloat(nOHmean(3))*VOLSHL)
         WRITE(n16,39) RADIUS,GROO_HDL(kw)                         
      ENDIF
   END DO
   WRITE(n16,*)
   

endif !Print LSI on the fly
!================================================================================================End Version 26 
      
!Print mean orient. tcfs      
write(n2,9)
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a5)')'#Bulk'
   if(kr==2)write(n1,'(a4)')'#LDL'
   if(kr==3)write(n1,'(a4)')'#HDL'
   if(kr==4)write(n1,'(a6)')'#OTHER'
   if(kr==1)write(n2,'(a5)')'#Bulk'
   if(kr==2)write(n2,'(a4)')'#LDL'
   if(kr==3)write(n2,'(a4)')'#HDL'
   if(kr==4)write(n2,'(a6)')'#OTHER'
   if(kr==1)write(n3,'(a5)')'#Bulk'
   if(kr==2)write(n3,'(a4)')'#LDL'
   if(kr==3)write(n3,'(a4)')'#HDL'
   if(kr==4)write(n3,'(a6)')'#OTHER'
   do it = 1,NDELS
      jt = it -1                                                     !print tcf starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      WRITE(n1,19)PSECS,tcf_mean_1(it,kr)/dfloat(ko)     !mean orientational tcf rank 1 over different time-origins 
      WRITE(n2,19)PSECS,tcf_mean_2(it,kr)/dfloat(ko)     !mean orientational tcf rank 2 over different time-origins   
      WRITE(n3,19)PSECS,tcf_mean_3(it,kr)/dfloat(ko)     !mean orientational tcf rank 3 over different time-origins 
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
end do

!Calculate the relaxation times - tau - trapezoid method
do kr = kr_start,NRSHMAX
   if(kr==1)then
      write(n11,'(a5)')'#Bulk'
      write(n12,'(a5)')'#Bulk'
      write(n13,'(a5)')'#Bulk'
   elseif(kr==2)then
      write(n11,'(a4)')'#LDL'
      write(n12,'(a4)')'#LDL'
      write(n13,'(a4)')'#LDL'
   elseif(kr==3)then
      write(n11,'(a4)')'#HDL'
      write(n12,'(a4)')'#HDL'
      write(n13,'(a4)')'#HDL'
   elseif(kr==4)then
      write(n11,'(a6)')'#OTHER'
      write(n12,'(a6)')'#OTHER'
      write(n13,'(a6)')'#OTHER'
   endif  
   koc(kr)=ko
!   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
         
   if(kr==1)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==2)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==3)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==4)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   endif   
end do 

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n11)
close(n12)
close(n13)
close(n20)

deallocate(nOHorg,NOHidO,NOHidH)
deallocate(ntime)
deallocate(tcf_mean_1,tcf_mean_2,tcf_mean_3)
deallocate(norg_0HB,kv,nOHmean)
deallocate(koc)

deallocate(list_s)
deallocate(list_LSI)
deallocate(LSI)
deallocate(NLSI)
deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT('#',17X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 39   FORMAT(27X,F9.4,5X,F14.5)
 
   return

END SUBROUTINE pure_wat_LDL_HDL_otcf



SUBROUTINE pure_wat_LDL_HDL_otcf_HB(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay_hb,m_time_hb,Nenv_hb,kLDL_HDL_)  
! Lars Pettersson and Gaia Camisasca collaboration - September 2018 - coded in version 25
! The program reads a list from an input file with information on LDL and HDL waters - Gaia's list
! Calculate OH reorientational tcfs - 1st, 2nd, and 3rd Legendre Polynomials - pure water LDL/HDL populations
! Every m_time_hb time-steps the code samples water molecules in the system that are either LDL, HDL or neither - output 3 orientational tcfs
! Deconvolute the otcfs for Bulk, LDL, and HDL water in terms of short and and long HBs per molecule
! The following HB definition is used: 
! O-O distance; angle OH-OO; distance H...O
! kLDL_HDL_ == 0 read LDL/HDL from input file; kLDL_HDL_ == 1 calculate LSI on the fly
!
    integer,intent(in)                               :: npurew,nbox,Nenv_hb,kLDL_HDL_
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: mtdelay_hb     ! reorientational tcf delay time (time-windows/ps)
    integer,intent(in)                               :: m_time_hb      ! time freq. to sample water OH groups in radial shells for reor. tcf calculation
   
! Local variables
!NG    real(kind=8)                               :: x(natms),y(natms),z(natms) 
!    real(kind=4)                               :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG    real(kind=4)                               :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                               :: cell(3)  
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    
    real(kind=8),dimension(:,:,:),allocatable  :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable  :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                               :: TSTEP,PSECS
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    real(kind=8)                               :: xoh,yoh,zoh,sqoh
    integer                                    :: NDELS,KINTVL,NDELS1
    integer,dimension(:,:,:),allocatable       :: NOHidO,NOHidH
    integer,dimension(:,:),allocatable         :: ntime
    integer,dimension(:,:),allocatable         :: nOHorg
    integer,dimension(:),allocatable           :: kv
    integer      :: n0, n1, n2, n3, n11, n12, n13, n15, n16, n17, n20, n25
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt, jw
    integer      :: io,jo,iox
    integer,parameter                          :: INTMAX = 2147483647
    integer                                    :: norgmax
    integer                                    :: NRSHMAX   
    integer,dimension(:),allocatable           :: norg_0HB  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!Orientational tcf
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b         !bulk-s
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs        !bulk-l
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th        !LDL-s
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth       !LDL-l
!Version 25
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_hdl_s,tcf_leg_2_hdl_s,tcf_leg_3_hdl_s
    real(kind=8),dimension(:),allocatable        :: tcf_leg_1_hdl_l,tcf_leg_2_hdl_l,tcf_leg_3_hdl_l
    
    real(kind=8),dimension(:,:),allocatable      :: tcf_mean_1   ,tcf_mean_2   ,tcf_mean_3   
    real(kind=8)                                 :: X0TLEG,ARGLEG1,ARGLEG2,ARGLEG3
    integer,dimension(:),allocatable             :: nOHmean
    integer                                      :: kr_start
    integer,dimension(:),allocatable             :: koc                  !Nb of origins corrected for origins where no waters were found in a given environment
!Version 24    
    integer                                    :: nw_list
    integer,dimension(:),allocatable           :: list_s                 !LSI index - read from input file List_W.dat n20 - Gaia
    integer,dimension(:),allocatable           :: list_LSI               !LSI index - LSI computed locally to compare with Gaia List
    integer,dimension(:,:),allocatable         :: list_HB
!    integer,dimension(:,:),allocatable         :: list_HB_bifur
    integer,dimension(:,:),allocatable         :: list_HB_HO
!Version 25
!H-Bond analysis
    integer                                    :: ihb,jhb
    real(kind=8)                               :: dxoo,dyoo,dzoo
    real(kind=8)                               :: SQDOO,DOXYG
    real(kind=8)                               :: dxoh,dyoh,dzoh
    real(kind=8)                               :: SQDOH,DOHYD
    real(kind=8)                               :: dxho,dyho,dzho
    real(kind=8)                               :: SQDHO,DHOXY
    real(kind=8)                               :: RHBOND,AHBOND,RHO_dist
    real(kind=8)                               :: DCOSNM,DCOSDN,DCOSANG,ANGDOHH
    real(kind=8),dimension(:,:),allocatable    :: d_Hi_Oj,d_Oi_Oj
    real(kind=8)                               :: RH_O_min
    logical                                    :: L_HO_decv, LHBdecv
!    logical                                    :: LHB_bifur
    logical                                    :: LHBdecv_OH
    integer                                    :: nsh_ldl,NHB_LDL, nsh_hdl,NHB_HDL,NOO_LDL,NOO_HDL,NHBA_LDL,NHBA_HDL
    integer                                    :: nkr,NRDELS,ANGDELS
    integer,dimension(:),allocatable           :: NDHBLDL,NDHBHDL,NDOOLDL,NDOOHDL,NHBALDL,NHBAHDL
    real(kind=8)                               :: AVRHBLDL,AVRHBHDL,AVROOLDL,AVROOHDL,AVRHBALDL,AVRHBAHDL
    real(kind=8)                               :: RHB_del, AHB_del
    real(kind=8)                               :: RHB_ANG,DISTR_HB,AVRHB_LDL,AVRHB_HDL  
    integer,dimension(:),allocatable           :: NOO
    integer                                    :: Nneighb
    integer                                    :: pwnmolsol,pwnatmsol,pwnions                !solution parameters - set to 0
    real(kind=8),dimension(:),allocatable      :: LSI
    real(kind=8)                               :: r_lsi,LSI_LDL_trh,LSI_HDL_trh
    

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

!Version 25 H-Bond criteria
  LHBdecv = .true.                                ! HBs per molecule deconvolution. For each molecule assume a short HB and a long HB - split OH groups equally
  LHBdecv_OH = .true.                            ! HBs deconvolution independent of molecule id - use HB threshold independent of the OH id
!  LHB_bifur = .false.                
  RHBOND = 3.5d0 
  AHBOND = 30.0d0
! HB O...H distance threshold LHBdecv_OH = .true.   
  RHO_dist = 2.0                                  !This is introduced in the HB definition only for LHBdecv_OH = .true.; deconvolution HBonded or not Hbonded OH groups
  if(LHBdecv_OH)RHBOND = 2.8d0                    !Restrict further the HB definition to match LDL for comparison
!  
  if(.not.LHBdecv)LHBdecv_OH = .false.
  if(.not.LHBdecv)RHBOND = 3.7d0                  !If HB not used for deconvolution use number of neighbors (Nneighb) up to RHBOND
  Nneighb = 5                                     !default 5: deconvolution < 5 and >= 5 up to RHBOND (use RHBOND = r_lsi = 3.7d0 )
!Version 25 HB Hi...Oj length distributions  
  RHB_del = 0.02d0
  AHB_del = 0.1d0
  NRDELS = 5.0d0/RHB_del
  ANGDELS = 180.0d0/AHB_del
  r_lsi  = 3.70d0                     !3.7d0
  LSI_LDL_trh = 0.17d0                !Gaia and Lars default values
  LSI_HDL_trh = 0.01d0                !Gaia and Lars default values

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Leg_LDL_HDL_HB/log_PW_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Leg_LDL_HDL_HB')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
!
n11 = 210
n12 = 220
n13 = 230
!
n15 = 250
n16 = 260
n17 = 270
!
n20 = 300

!check   n25 = 500

open(n0,file='Leg_LDL_HDL_HB/log_PW_tcf.dat',status='unknown',action='write')
open(n1,file='Leg_LDL_HDL_HB/wat_reor_tcf_1.dat',status='unknown',action='write')
open(n2,file='Leg_LDL_HDL_HB/wat_reor_tcf_2.dat',status='unknown',action='write')
open(n3,file='Leg_LDL_HDL_HB/wat_reor_tcf_3.dat',status='unknown',action='write')
!
open(n11,file='Leg_LDL_HDL_HB/tcf_1_rot_tau.dat',status='unknown',action='write')
open(n12,file='Leg_LDL_HDL_HB/tcf_2_rot_tau.dat',status='unknown',action='write')
open(n13,file='Leg_LDL_HDL_HB/tcf_3_rot_tau.dat',status='unknown',action='write')
!
open(n15,file='Leg_LDL_HDL_HB/P_HB_HO_dist.dat',status='unknown',action='write')
open(n16,file='Leg_LDL_HDL_HB/P_HB_OO_dist.dat',status='unknown',action='write')
open(n17,file='Leg_LDL_HDL_HB/P_HB_theta_dist.dat',status='unknown',action='write')

!Version 24 - read LDL/HDL water molecules label from input file
!Lars Pettersson - list
! 1 is LDL 
! 2 is HDL
! 3 skip

!kLDL_HDL_ == 0 read LDL/HDL from input file; kLDL_HDL_ == 1 calculate LSI on the fly
if(Nenv_hb>1.and.kLDL_HDL_==0)open(n20,file='List_W.dat',status='old',action='read')

!check   open(n25,file='Leg_LDL_HDL_HB/check_List_W.dat',status='unknown',action='write')

!time-window for calculation of each tcf
! IDINT - convert to integer - truncate

NDELS = IDINT(mtdelay_hb*1000d0/dt)                  !sets the delay time (in frames)

kr_start = 1                                      !environments loop starter - for pure water single environment 

write(*,*)'input delay-time (ps) =',mtdelay_hb
write(*,*)'delay-time (input frames units) =',NDELS
write(n0,*)'input delay-time (ps) =',mtdelay_hb
write(n0,*)'delay-time (input frames units) =',NDELS

write(*,*)
write(*,*)'Number of Environments requested = ',Nenv_hb
write(*,*)

if(Nenv_hb>1.and.LHBdecv)then
   write(*,*)'Environments: 1(Bulk short-HB), 2(Bulk long-HB), 3(LDL short-HB)'
   write(*,*)'4(LDL long-HB),5(HDL short-HB), 6(HDL long-HB)'
   write(*,*)
endif
if(kLDL_HDL_==1)write(*,*)'LDL/HDL populations calculated on the fly'
if(kLDL_HDL_==0)write(*,*)'LDL/HDL populations read from input file'
if(LHBdecv)write(*,*)'otcfs deconvolution based on HBs'
if(LHBdecv_OH)write(*,*)'WARNING!!! OH groups treated as independent'
!if(LHBdecv.and.LHB_bifur)write(*,*)'deconvolution based on bifurcated HBs'
if(.not.LHBdecv)write(*,*)'WARNING!!! otcfs deconvolution based on local neighbors Nnb = ',Nneighb
if(.not.LHBdecv)write(*,*)'Radial cut-off for local neighbours Rnb (Ang) = ', RHBOND
write(*,*)

write(n0,*)
write(n0,*)'Number of Environments requested = ',Nenv_hb
write(n0,*)

if(Nenv_hb>1.and.LHBdecv)then
   write(n0,*)'Environments: 1(Bulk short-HB), 2(Bulk long-HB), 3(LDL short-HB)'
   write(n0,*)'4(LDL long-HB),5(HDL short-HB), 6(HDL long-HB)'
   write(n0,*) 
endif
if(kLDL_HDL_==1)write(n0,*)'LDL/HDL populations calculated on the fly'
if(kLDL_HDL_==0)write(n0,*)'LDL/HDL populations read from input file'
if(LHBdecv)write(n0,*)'otcfs deconvolution based on HBs'
if(LHBdecv_OH)write(n0,*)'WARNING!!! OH groups treated as independent'
!if(LHBdecv.and.LHB_bifur)write(n0,*)'deconvolution based on bifurcated HBs'
if(.not.LHBdecv)write(n0,*)'WARNING!!! otcfs deconvolution based on local neighbors Nnb = ',Nneighb
if(.not.LHBdecv)write(n0,*)'Radial cut-off for local neighbors Rnb (Ang) = ', RHBOND
write(n0,*)

NRSHMAX = Nenv_hb                                   !maximum number of environments - 4 e.g. LDL, HDL, rest, all
!NRSHMAX = 1                                     !maximum number of environments - 1 for pure water

if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/m_time_hb                                  !for memory allocation purposes

!orient. tcf
allocate(nOHorg(norgmax,NRSHMAX))
allocate(RXt0(norgmax,2*nmolwat,NRSHMAX),RYt0(norgmax,2*nmolwat,NRSHMAX),RZt0(norgmax,2*nmolwat,NRSHMAX))
allocate(RXt(norgmax,2*nmolwat,NRSHMAX),RYt(norgmax,2*nmolwat,NRSHMAX),RZt(norgmax,2*nmolwat,NRSHMAX))
allocate(NOHidO(norgmax,2*nmolwat,NRSHMAX),NOHidH(norgmax,2*nmolwat,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(tcf_leg_1_b(NDELS),tcf_leg_1_hs(NDELS),tcf_leg_1_th(NDELS),tcf_leg_1_nth(NDELS),tcf_leg_1_hdl_s(NDELS),tcf_leg_1_hdl_l(NDELS))
allocate(tcf_leg_2_b(NDELS),tcf_leg_2_hs(NDELS),tcf_leg_2_th(NDELS),tcf_leg_2_nth(NDELS),tcf_leg_2_hdl_s(NDELS),tcf_leg_2_hdl_l(NDELS))
allocate(tcf_leg_3_b(NDELS),tcf_leg_3_hs(NDELS),tcf_leg_3_th(NDELS),tcf_leg_3_nth(NDELS),tcf_leg_3_hdl_s(NDELS),tcf_leg_3_hdl_l(NDELS))
allocate(tcf_mean_1(NDELS,NRSHMAX) ,tcf_mean_2(NDELS,NRSHMAX)  ,tcf_mean_3(NDELS,NRSHMAX))

tcf_leg_1_b = 0.0d0;tcf_leg_1_hs = 0.0d0;tcf_leg_1_th = 0.0d0;tcf_leg_1_nth = 0.0d0;tcf_leg_1_hdl_s = 0.0d0;tcf_leg_1_hdl_l = 0.0d0
tcf_leg_2_b = 0.0d0;tcf_leg_2_hs = 0.0d0;tcf_leg_2_th = 0.0d0;tcf_leg_2_nth = 0.0d0;tcf_leg_2_hdl_s = 0.0d0;tcf_leg_2_hdl_l = 0.0d0
tcf_leg_3_b = 0.0d0;tcf_leg_3_hs = 0.0d0;tcf_leg_3_th = 0.0d0;tcf_leg_3_nth = 0.0d0;tcf_leg_3_hdl_s = 0.0d0;tcf_leg_3_hdl_l = 0.0d0

tcf_mean_1  = 0.0d0;tcf_mean_2   = 0.0d0;tcf_mean_3   = 0.0d0 
ntime = 1; nOHorg = 0



allocate(norg_0HB(NRSHMAX))
norg_0HB = 0
allocate(kv(NRSHMAX))
allocate(nOHmean(NRSHMAX))
allocate(koc(NRSHMAX))
nOHmean = 0
koc = 0
!Version 24 - List HDL/LDL
allocate(list_s(nmolwat))
allocate(list_LSI(nmolwat))
allocate(LSI(nmolwat))
!Version 25 - H-BOND analysis 
allocate(list_HB(nmolwat,2))
!allocate(list_HB_bifur(nmolwat,2))
allocate(list_HB_HO(nmolwat,2))
allocate(d_Hi_Oj(nmolwat,2))
allocate(d_Oi_Oj(nmolwat,2))
allocate(NDHBLDL(NRDELS),NDHBHDL(NRDELS))
allocate(NDOOLDL(NRDELS),NDOOHDL(NRDELS))
allocate(NHBALDL(ANGDELS),NHBAHDL(ANGDELS))
allocate(NOO(nmolwat))

list_HB = 0
!list_HB_bifur = 0
list_HB_HO = 0
d_Hi_Oj = 0
d_Oi_Oj = 0

NDHBLDL = 0
NDOOLDL = 0
NHBALDL = 0

NDHBHDL = 0
NDOOHDL = 0
NHBAHDL = 0

AVRHBLDL = 0.0d0
AVRHBHDL = 0.0d0

AVROOLDL = 0.0d0
AVROOHDL = 0.0d0

AVRHBALDL = 0.0d0
AVRHBAHDL = 0.0d0

NHB_LDL = 0
NHB_HDL = 0

NOO_LDL = 0
NOO_HDL = 0

NHBA_LDL = 0
NHBA_HDL = 0

NOO = 0

LSI = 0

pwnmolsol = 0
pwnatmsol = 0
pwnions   = 0  

write(*,*)
write(*,*)'Water reorientation tcfs calculation'
write(*,*)'Routine pure_wat_LDL_HDL_otcf_HB'
write(*,*)'water tcfs delay time-window (ps) =',mtdelay_hb
write(*,*)'Results are printed to Leg_LDL_HDL_HB'
write(*,*)
write(n0,*)
write(n0,*)'Water reorientation tcfs calculation'
write(n0,*)'Routine pure_wat_LDL_HDL_otcf_HB'
write(n0,*)'water tcfs delay time-window (ps) =',mtdelay_hb
write(n0,*)'Results are printed to Leg_LDL_HDL_HB'
write(n0,*)

! Start tcfs calculation
write(*,*)
write(*,*)'Starting reor. tcf calculation...'
write(*,*)'tcfs OH sampling frequency (steps) = ',m_time_hb
write(*,*)'total number of time-steps =',nstep
write(*,*)'tcfs delay-time (steps) =',NDELS
write(*,*)
write(n0,*)
write(n0,*)'Starting reor. tcf calculation...'
write(n0,*)'tcfs OH sampling frequency (steps) = ',m_time_hb
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'tcfs delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif      
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!End Version 9      
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj

!============================================================================================ development version 25 LDL/HDL OH H-bond deconvolution

!HB analysis
!Development version 25
   
   if((j==1.or.mod(j,m_time_hb)==0).and.(j.le.nstep-NDELS))then                 !time origin: sample OH groups in specific environments
   write(*,*)'time-step',j,' sampling for OH groups'
   
!LDL/HDL list - identify LDL and HDL water molecules - Gaia Lists     
      if(NRSHMAX > 1.and.kLDL_HDL_==0)then
         do i=1,nmolwat
            read(n20,*)nw_list,list_s(i)        
!check            write(n25,*)nw_list,list_s(i)
         end do
      endif
      
!Version 26    ========================================================================================================
!Calculate the LSI
      if(kLDL_HDL_==1)then              !calculate LDL/HDL populations on the fly
         io = 0
         do ihb = 1,natms,nwatsites
            LSI = 0
            io = io + 1
            call LSI_water(ihb,io,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,LSI,r_lsi)
!check LSI - compare with list_s(i) based on Gaia & Lars LSI cut-offs for LDL (> 0.17 Ang**2) and HDL (< 0.01 Ang**2) species            
            if(LSI(io) > LSI_LDL_trh)then
               list_LSI(io) = 1               !LDL
            elseif(LSI(io) < LSI_HDL_trh)then 
               list_LSI(io) = 2               !HDL
            else
               list_LSI(io) = 3               !other
            endif   
!check            if(list_LSI(io).ne.list_s(io))then
!check               write(*,*)'Error - LSI id different from input list'
!check               write(*,*)'LSI id = ',list_LSI(io),'LSI = ',LSI(io)
!check               write(*,*)'input list LSI id = ',list_s(io)
!check               write(*,*)list_s(io),list_LSI(io),LSI(io)
!check            endif
         end do
      
         list_s = list_LSI
      endif
            
!End Version 26 ========================================================================================================               
      
      if(LHBdecv)then               !HB analysis for deconvolution of the otcfs
      
         list_HB = 0                                   !re-initialize HB list every new sampling origin
!         list_HB_bifur = 0                             !re-initialize bifurcated HB list every new sampling origin
         list_HB_HO = 0                                !re-initialize HB list every new sampling origin - OH group independent HB analysis
         io = 0                                        !oxygen atom counter 1,2,3,...
         do ihb = 1,natms,nwatsites                    !Loop over oxygen atoms 1,4,7,10,...
            jo = 0                                     !oxygen atom counter 1,2,3,...
            io = io + 1 
            do jhb = 1,natms,nwatsites                 !Loop over oxygen atoms 1,4,7,10,...
               jo = jo + 1 
               IF(jo.NE.io)THEN
!Intermolecular O-O distance: donor and acceptor calculations - DOXYG           
                  dx = x(ihb)-x(jhb)
                  dy = y(ihb)-y(jhb)
                  dz = z(ihb)-z(jhb)            
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)    
                  dxoo = dx
                  dyoo = dy
                  dzoo = dz
                  SQDOO = dxoo**2.0D0+dyoo**2.0D0+dzoo**2.0D0
                  DOXYG= DSQRT(SQDOO)
                  do iH = 1,2
!Intramolecular O(i)-H distances (iH=1 H1; iH=2 H2): donor calculations - DOHYD
                     dx = x(ihb) - x(ihb + iH)
                     dy = y(ihb) - y(ihb + iH)
                     dz = z(ihb) - z(ihb + iH)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dxoh = dx
                     dyoh = dy
                     dzoh = dz
                     SQDOH = dsqrt(dxoh*dxoh+dyoh*dyoh+dzoh*dzoh)
                     DOHYD = DSQRT(SQDOH)
!Intermolecular H(i)-O(j) distances: donor calculations
                     dx = x(ihb + iH) - x(jhb)
                     dy = y(ihb + iH) - y(jhb)
                     dz = z(ihb + iH) - z(jhb)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dxho = dx
                     dyho = dy
                     dzho = dz
                     SQDHO = dxho**2.0D0+dyho**2.0D0+dzho**2.0D0
                     DHOXY = DSQRT(SQDHO)
!Apply the COSINE rule to calculate the O-H...O angle: donor calculations
                     DCOSNM = SQDHO-SQDOH-SQDOO
                     DCOSDN = -2.0D0*DOHYD*DOXYG
                     DCOSANG = DCOSNM/DCOSDN
                     ANGDOHH = DACOSD(DCOSANG)
!Allow for an OH group to be H-bonded to more than one oxygen atom - keep shortest HB                                     
                     IF((DOXYG.LT.RHBOND).AND.(ANGDOHH.LT.AHBOND))THEN                !HB found
!check                        write(*,*)DOXYG,ANGDOHH,io,iH,jo
                        if(list_HB(io,iH)==0)then                                     !No HB was found before...
                           list_HB(io,iH) = 1
                           d_Hi_Oj(io,iH) = DHOXY                                     !HB distance
                           d_Oi_Oj(io,iH) = DOXYG                                     !oxygen-oxygen distance
                        elseif(list_HB(io,iH)==1)then                                 !HB was found before - keep the shortest HB
!                           list_HB_bifur(io,iH)=1
                           write(*,*)'Found a second HB on the same proton - checking HB length...'
                           write(*,*)'water molecule',io,' LDL/HDL id = ',list_s(io)
                           if(DHOXY < d_Hi_Oj(io,iH)) d_Hi_Oj(io,iH) = DHOXY
                           if(DOXYG < d_Oi_Oj(io,iH)) d_Oi_Oj(io,iH) = DOXYG              !if use O-O distance instead of H...O distance
                        endif
!OH group independent HB deconvolution                        
                        if(DHOXY < RHO_dist)then
                           list_HB_HO(io,iH) = 1
                        endif
!                        
                        if(list_s(io)==1)then                     !LDL
! CALCULATE Hi...Oj DISTRIBUTION AND AVERAGE VALUE LDL
                           nsh_ldl = DHOXY/RHB_del+0.5D0
                           NDHBLDL(nsh_ldl) = NDHBLDL(nsh_ldl) + 1
                           AVRHBLDL = AVRHBLDL + DHOXY
                           NHB_LDL = NHB_LDL + 1                                   !sums all HBs for every LDL water and over every origin sampled
! CALCULATE Oi...Oj DISTRIBUTION AND AVERAGE VALUE LDL
                           nsh_ldl = DOXYG/RHB_del+0.5D0
                           NDOOLDL(nsh_ldl) = NDOOLDL(nsh_ldl) + 1
                           AVROOLDL = AVROOLDL + DOXYG
                           NOO_LDL = NOO_LDL + 1                                   !sums all HBs for every LDL water and over every origin sampled  
! CALCULATE HB ANGLE DISTRIBUTION AND AVERAGE VALUE LDL
                           nsh_ldl = ANGDOHH/AHB_del+0.5D0
                           NHBALDL(nsh_ldl) = NHBALDL(nsh_ldl) + 1
                           AVRHBALDL = AVRHBALDL + ANGDOHH
                           NHBA_LDL = NHBA_LDL + 1                                 !sums all HBs for every LDL water and over every origin sampled                                 
                        elseif(list_s(io)==2)then                 !HDL
! CALCULATE Hi...Oj DISTRIBUTION AND AVERAGE VALUE HDL                     
                           nsh_hdl = DHOXY/RHB_del+0.5D0
                           NDHBHDL(nsh_hdl) = NDHBHDL(nsh_hdl) + 1
                           AVRHBHDL = AVRHBHDL + DHOXY
                           NHB_HDL = NHB_HDL + 1                                    !sums all HBs for every HDL water and over every origin sampled
! CALCULATE Oi...Oj DISTRIBUTION AND AVERAGE VALUE HDL
                           nsh_hdl = DOXYG/RHB_del+0.5D0
                           NDOOHDL(nsh_hdl) = NDOOHDL(nsh_hdl) + 1
                           AVROOHDL = AVROOHDL + DOXYG
                           NOO_HDL = NOO_HDL + 1                                    !sums all HBs for every LDL water and over every origin sampled   
! CALCULATE HB ANGLE DISTRIBUTION AND AVERAGE VALUE HDL
                           nsh_hdl = ANGDOHH/AHB_del+0.5D0
                           NHBAHDL(nsh_hdl) = NHBAHDL(nsh_hdl) + 1
                           AVRHBAHDL = AVRHBAHDL + ANGDOHH
                           NHBA_HDL = NHBA_HDL + 1                                  !sums all HBs for every LDL water and over every origin sampled                                     
                        endif
                        
                     ENDIF  !end HB found
                  end do    !end hydrogen loop iH    
               ENDIF        !end if io =/ jo
            end do          !end oxygen loop   jhb
         end do             !end oxygen loop   ihb  

!      
!HB deconvolution exceptions
!                  (i)  OH group donates 2 H-bonds
!                  (ii) None of the OH groups form H-bonds     
         io = 0
         L_HO_decv = .true.                                                                 !.true. use H...O distance; .false. use O---O distance
         do ihb = 1,natms,nwatsites             
            io = io + 1
!HB length Hi...Oj      
            if(list_HB(io,1)==1.and.list_HB(io,2)==1)then                                   !The 2 OH are H-bonded - deconvolute short and long HB
               if(L_HO_decv)then
!                  write(*,*)'2 OH groups H-bonded - deconvolute based on H...O distance'
                  if(d_Hi_Oj(io,1) < d_Hi_Oj(io,2))then                                     !if HB on OH1 is shorter than HB on OH2
                     list_HB(io,2) = 0
                  elseif(d_Hi_Oj(io,2) < d_Hi_Oj(io,1))then   
                     list_HB(io,1) = 0
                  endif
!Oxygen-Oxygen length Oi...Oj
               elseif(.not.L_HO_decv)then
!                  write(*,*)'2 OH groups H-bonded - deconvolute based on O...O distance'
                  if(d_Oi_Oj(io,1) < d_Oi_Oj(io,2))then
                     list_HB(io,2) = 0
                  elseif(d_Oi_Oj(io,2) < d_Oi_Oj(io,1))then   
                     list_HB(io,1) = 0
                  endif
               endif   
               
!               write(*,*)'Oh NO there are two HBs',io
            elseif(list_HB(io,1)==0.and.list_HB(io,2)==0)then                             !None of the OH are H-bonded - consider the shortest Hi...Oj distance
               write(*,*)'OH groups not H-Bonded - deconvolute based on H...O distance'
               write(*,*)'water molecule = ', io
               RH_O_min = 5.0d0                              !Some large HB value in Angstroms
               list_HB(io,1) = 1                             !default 
               list_HB(io,2) = 0                             !default 
               jo = 0
               do jhb = 1,natms,nwatsites                    !Loop over oxygen atoms
                  jo = jo + 1 
                  IF(jo.NE.io)THEN
!Intermolecular H(i)-O(j) distances: donor calculations               
                     dx = x(ihb + 1) - x(jhb)                  !Hi...Oj
                     dy = y(ihb + 1) - y(jhb)
                     dz = z(ihb + 1) - z(jhb)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dxho = dx
                     dyho = dy
                     dzho = dz
                     SQDHO = dxho**2.0D0+dyho**2.0D0+dzho**2.0D0
                     DHOXY = DSQRT(SQDHO) 
                     if(DHOXY < RH_O_min)then 
                        RH_O_min = DHOXY                                !Store the shortest Hi...Oj distance
                     endif                  
                  ENDIF
               end do
               write(*,*)'NO-HB r(H...O) min = ',RH_O_min
!Find Hi...Oj for the second Hi atom            
               jo = 0           
               do jhb = 1,natms,nwatsites                 !Loop over oxygen atoms
                  jo = jo + 1 
                  IF(jo.NE.io)THEN
!Intermolecular H(i)-O(j) distances: donor calculations               
                     dx = x(ihb + 2) - x(jhb)
                     dy = y(ihb + 2) - y(jhb)
                     dz = z(ihb + 2) - z(jhb)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dxho = dx
                     dyho = dy
                     dzho = dz
                     SQDHO = dxho**2.0D0+dyho**2.0D0+dzho**2.0D0
                     DHOXY = DSQRT(SQDHO) 
                     if(DHOXY < RH_O_min)then 
                        list_HB(io,1) = 0
                        list_HB(io,2) = 1 
                        write(*,*)'FOUND SHORTER NO-HB r(H...O) min = ',DHOXY
                        goto 5
                     endif                  
                  ENDIF
               end do   
               5 continue
               
            endif !end HB deconvolution conditional block        
         end do   !end loop over oxygen atoms      
   
!Check HB deconvolution

         do iox = 1,nmolwat
            if(list_HB(iox,1)==1.and.list_HB(iox,2)==1)then
               write(*,*)'Error!!! OH groups not deconvoluted 1 1'
               stop
            elseif(list_HB(iox,1)==0.and.list_HB(iox,2)==0)then
               write(*,*)'Error!!! OH groups not deconvoluted 0 0'
               stop
            endif
         end do
         
      elseif(.not.LHBdecv)then
      
!Deconvolute water according to local neighborhood - defects
        
         list_HB = 0                                   !re-initialize HB list every new sampling origin
         NOO = 0
         io = 0                                        !oxygen atom counter 1,2,3,...
         do ihb = 1,natms,nwatsites                    !Loop over oxygen atoms 1,4,7,10,...
            jo = 0                                     !oxygen atom counter 1,2,3,...
            io = io + 1 
            do jhb = 1,natms,nwatsites                 !Loop over oxygen atoms 1,4,7,10,...
               jo = jo + 1 
               IF(jo.NE.io)THEN
                  dx = x(ihb)-x(jhb)
                  dy = y(ihb)-y(jhb)
                  dz = z(ihb)-z(jhb)            
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)    
                  dxoo = dx
                  dyoo = dy
                  dzoo = dz
                  SQDOO = dxoo**2.0D0+dyoo**2.0D0+dzoo**2.0D0
                  DOXYG= DSQRT(SQDOO)
                  IF(DOXYG.LT.RHBOND)NOO(io) = NOO(io) + 1
               ENDIF
            enddo
!check            if(list_s(io)==1)write(*,*)ihb,io,NOO(io),list_s(io)
!check            if(list_s(io)==2)write(*,*)ihb,io,NOO(io),list_s(io)
            if(NOO(io) < Nneighb)then
               list_HB(io,1)=1              ! No defects
               list_HB(io,2)=1              ! No defects
            elseif(NOO(io) >= Nneighb)then
               list_HB(io,1)=0              ! defects
               list_HB(io,2)=0              ! defects
            endif   
         end do   

!====================================================================================development        

      endif                 !HB if(logical) - HB/RO5 analysis for otcfs deconvolution     
      
!End development version 25
!============================================================================================End development version 25

!check      if(LHB_bifur)list_HB=list_HB_bifur
      if(LHBdecv_OH)list_HB=list_HB_HO

!sample water molecules in specific environments - pure water - single environment
! VERSION 25 MOVED UP 

!V25   if((j==1.or.mod(j,m_time_hb)==0).and.(j.le.nstep-NDELS))then               !time origin: sample OH groups in specific environments
!V25      write(*,*)'time-step',j,' sampling for OH groups'
!      
!Version 24 - LDL/HDL water list           1(LDL) 2(HDL) 3(NOT LDL nor HDL) 
!V25      if(NRSHMAX > 1)then
!V25         do i=1,nmolwat
!V25            read(n20,*)nw_list,list_s(i)
!V25         end do
!V25      endif      


! otcf      
      kv = 0                                              !OH groups in each environment 
      ko = ko + 1                                         !number of origins = number of blocks used to compute the tcfs mean and stdev
      io = 0                                              !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
      
      do i=1,natms,nwatsites                              !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...

!bulk tcf(t0) - all waters         
!NG         if(NRSHMAX > 1)then                               !Version 25
            do iH = 1,2
               if(list_HB(io,iH)==1)then
                  nradsh = 1                                     !bulk environment - short HBs
                  kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                  kvi = kv(nradsh)                               !OH number
!check                  write(*,*)j,i,iH,kvi
                  nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                  NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  xoh = dx
                  yoh = dy
                  zoh = dz
                  sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                  RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = yoh/sqoh
                  RZt0(ko,kvi,nradsh) = zoh/sqoh
               
               elseif(list_HB(io,iH)==0)then
                  nradsh = 2                                     !bulk environment - long HBs
                  kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                  kvi = kv(nradsh)                               !OH number
!check                  write(*,*)j,i,iH,kvi
                  nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                  NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                  NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - dnint(dx/cell(1))*cell(1)
                  dy = dy - dnint(dy/cell(2))*cell(2)
                  dz = dz - dnint(dz/cell(3))*cell(3)
                  xoh = dx
                  yoh = dy
                  zoh = dz
                  sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                  RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                  RYt0(ko,kvi,nradsh) = yoh/sqoh
                  RZt0(ko,kvi,nradsh) = zoh/sqoh     
               endif
            end do               !end iH loop hydrogens loop
!bulk tcf(t0) end

!LDL tcf(t0)
            if(NRSHMAX > 2)then
               if(list_s(io)==1)then
                  do iH = 1,2
                     if(list_HB(io,iH)==1)then
                        nradsh = 3                                     !LDL environment - short HBs
                        kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                        kvi = kv(nradsh)                               !OH number
!check                        write(*,*)j,i,iH,kvi
                        nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                        NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                        NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                        dx = x(i+iH) - x(i)
                        dy = y(i+iH) - y(i)
                        dz = z(i+iH) - z(i)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        xoh = dx
                        yoh = dy
                        zoh = dz
                        sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                        RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                        RYt0(ko,kvi,nradsh) = yoh/sqoh
                        RZt0(ko,kvi,nradsh) = zoh/sqoh   
                     elseif(list_HB(io,iH)==0)then
                        nradsh = 4                                     !LDL environment - long HBs
                        kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                        kvi = kv(nradsh)                               !OH number
!check                        write(*,*)j,i,iH,kvi
                        nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                        NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                        NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                        dx = x(i+iH) - x(i)
                        dy = y(i+iH) - y(i)
                        dz = z(i+iH) - z(i)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        xoh = dx
                        yoh = dy
                        zoh = dz
                        sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                        RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                        RYt0(ko,kvi,nradsh) = yoh/sqoh
                        RZt0(ko,kvi,nradsh) = zoh/sqoh   
                     endif 
                  end do   
               endif             !End LDL
            endif
!LDL tcf(t0) end

!HDL tcf(t0)
            if(NRSHMAX > 4)then
               if(list_s(io)==2)then
                  do iH = 1,2
                     if(list_HB(io,iH)==1)then
                        nradsh = 5                                     !HDL environment - short HBs
                        kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                        kvi = kv(nradsh)                               !OH number
!check                        write(*,*)j,i,iH,kvi
                        nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                        NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                        NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                        dx = x(i+iH) - x(i)
                        dy = y(i+iH) - y(i)
                        dz = z(i+iH) - z(i)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        xoh = dx
                        yoh = dy
                        zoh = dz
                        sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                        RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                        RYt0(ko,kvi,nradsh) = yoh/sqoh
                        RZt0(ko,kvi,nradsh) = zoh/sqoh              
                     elseif(list_HB(io,iH)==0)then
                        nradsh = 6                                     !HDL environment - long HBs
                        kv(nradsh) = kv(nradsh) + 1                    !number of water OH groups (index)
                        kvi = kv(nradsh)                               !OH number
!check                        write(*,*)j,i,iH,kvi
                        nOHorg(ko,nradsh) = nOHorg(ko,nradsh) + 1      !number of OH groups in environment nradsh at origin ko
                        NOHidO(ko,kvi,nradsh) = i                      ! Oxygen id
                        NOHidH(ko,kvi,nradsh) = iH                     ! Hydrogen id 
                        dx = x(i+iH) - x(i)
                        dy = y(i+iH) - y(i)
                        dz = z(i+iH) - z(i)
                        dx = dx - dnint(dx/cell(1))*cell(1)
                        dy = dy - dnint(dy/cell(2))*cell(2)
                        dz = dz - dnint(dz/cell(3))*cell(3)
                        xoh = dx
                        yoh = dy
                        zoh = dz
                        sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
                        RXt0(ko,kvi,nradsh) = xoh/sqoh                 !origin = ko, OH number = kvi, environment type = nradsh
                        RYt0(ko,kvi,nradsh) = yoh/sqoh
                        RZt0(ko,kvi,nradsh) = zoh/sqoh  
                     endif   
                  end do         
!HDL tcf(t0) end
               endif           !end HDL water
            endif
            
            
            
            
!         endif              !end if NRSHMAX > 1        
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nOHorg(ko,1)
   endif                                                             !end sampling new time-origin                                                                 
   
!calculate tcfs - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nOHorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then                !ntime counts the delay times for which the tcf has already been calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nOHorg(kt,kr)                                       !loop over OH groups in environment kr: nOHorg can be zero for a given time-origin
                  
               i=NOHidO(kt,iw,kr)
               ihat = NOHidH(kt,iw,kr)                                    !Hydrogen atom (ihat = 1 or 2)        
               dx = x(i+ihat) - x(i)
               dy = y(i+ihat) - y(i)
               dz = z(i+ihat) - z(i)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               xoh = dx
               yoh = dy
               zoh = dz
               sqoh = dsqrt(xoh*xoh+yoh*yoh+zoh*zoh)
               RXt(kt,iw,kr) = xoh/sqoh
               RYt(kt,iw,kr) = yoh/sqoh
               RZt(kt,iw,kr) = zoh/sqoh
               X0TLEG = RXt0(kt,iw,kr)*RXt(kt,iw,kr)+RYt0(kt,iw,kr)*RYt(kt,iw,kr)+RZt0(kt,iw,kr)*RZt(kt,iw,kr)
               ARGLEG1= X0TLEG
               ARGLEG2= 0.5d0*(3.0d0*X0TLEG**2.0-1.0d0)
               ARGLEG3= 0.5d0*(5.0d0*X0TLEG**3.0-3.0d0*X0TLEG)
               nt = ntime(kt,kr)
!              
               if(kr==1)then                                            !bulk
                  tcf_leg_1_b(nt)=tcf_leg_1_b(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_b(nt)=tcf_leg_2_b(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_b(nt)=tcf_leg_3_b(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==2)then                                        !LDL
                  tcf_leg_1_hs(nt)=tcf_leg_1_hs(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_hs(nt)=tcf_leg_2_hs(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_hs(nt)=tcf_leg_3_hs(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==3)then                                        !HDL
                  tcf_leg_1_th(nt)=tcf_leg_1_th(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_th(nt)=tcf_leg_2_th(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_th(nt)=tcf_leg_3_th(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==4)then                                        !OTHER - NOT LDL nor HDL
                  tcf_leg_1_nth(nt)=tcf_leg_1_nth(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_nth(nt)=tcf_leg_2_nth(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_nth(nt)=tcf_leg_3_nth(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==5)then                                        !OTHER - NOT LDL nor HDL
                  tcf_leg_1_hdl_s(nt)=tcf_leg_1_hdl_s(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_hdl_s(nt)=tcf_leg_2_hdl_s(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_hdl_s(nt)=tcf_leg_3_hdl_s(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))
               elseif(kr==6)then                                        !OTHER - NOT LDL nor HDL
                  tcf_leg_1_hdl_l(nt)=tcf_leg_1_hdl_l(nt)+ARGLEG1/dfloat(nOHorg(kt,kr))
                  tcf_leg_2_hdl_l(nt)=tcf_leg_2_hdl_l(nt)+ARGLEG2/dfloat(nOHorg(kt,kr))
                  tcf_leg_3_hdl_l(nt)=tcf_leg_3_hdl_l(nt)+ARGLEG3/dfloat(nOHorg(kt,kr))   
                  endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean tcfs over different origins
do kr = 1,NRSHMAX         !loop over environments
   do it = 1,NDELS
!nOHorg(kt,kr) - number of OH groups found in origin (block) kt and in the environment kr               
      if(kr==1)then
         tcf_mean_1(it,kr) = tcf_leg_1_b(it)
         tcf_mean_2(it,kr) = tcf_leg_2_b(it)
         tcf_mean_3(it,kr) = tcf_leg_3_b(it)
      elseif(kr==2)then   
         tcf_mean_1(it,kr) = tcf_leg_1_hs(it)
         tcf_mean_2(it,kr) = tcf_leg_2_hs(it)
         tcf_mean_3(it,kr) = tcf_leg_3_hs(it)
      elseif(kr==3)then
         tcf_mean_1(it,kr) = tcf_leg_1_th(it)
         tcf_mean_2(it,kr) = tcf_leg_2_th(it)
         tcf_mean_3(it,kr) = tcf_leg_3_th(it)
      elseif(kr==4)then
         tcf_mean_1(it,kr) = tcf_leg_1_nth(it)
         tcf_mean_2(it,kr) = tcf_leg_2_nth(it)
         tcf_mean_3(it,kr) = tcf_leg_3_nth(it)
      elseif(kr==5)then
         tcf_mean_1(it,kr) = tcf_leg_1_hdl_s(it)
         tcf_mean_2(it,kr) = tcf_leg_2_hdl_s(it)
         tcf_mean_3(it,kr) = tcf_leg_3_hdl_s(it)
      elseif(kr==6)then
         tcf_mean_1(it,kr) = tcf_leg_1_hdl_l(it)
         tcf_mean_2(it,kr) = tcf_leg_2_hdl_l(it)
         tcf_mean_3(it,kr) = tcf_leg_3_hdl_l(it)   
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of OH groups on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nOHmean(kr) = nOHmean(kr) + nOHorg(kt,kr)
   end do
end do

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a28,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko) 
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2) 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a28,1x,F14.2)')'#Mean number of OH groups = ',dfloat(nOHmean(kr))/dfloat(ko)
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nOHmean(kr))/dfloat(ko*2) 
end do   
      
deallocate(tcf_leg_1_b  ,tcf_leg_2_b  ,tcf_leg_3_b)
deallocate(tcf_leg_1_hs ,tcf_leg_2_hs ,tcf_leg_3_hs)
deallocate(tcf_leg_1_th ,tcf_leg_2_th ,tcf_leg_3_th)
deallocate(tcf_leg_1_nth,tcf_leg_2_nth,tcf_leg_3_nth)
deallocate(tcf_leg_1_hdl_s,tcf_leg_2_hdl_s,tcf_leg_3_hdl_s)
deallocate(tcf_leg_1_hdl_l,tcf_leg_2_hdl_l,tcf_leg_3_hdl_l)

!Print HB distributions - HB length H...O; O-O distance; HB angle (< 30 deg)

if(LHBdecv)then
!LDL
   write(n15,'(a4)')'#LDL'
   do nkr = 1,NRDELS
      if(NDHBLDL(nkr) > 0)then
         RHB_ANG  = RHB_del*DFLOAT(nkr)          
         DISTR_HB = DFLOAT(NDHBLDL(nkr))/(RHB_del*DFLOAT(NHB_LDL))                    !Normalize
         WRITE(n15,39)RHB_ANG,DISTR_HB
      endif
   end do  
   AVRHB_LDL = AVRHBLDL/DFLOAT(NHB_LDL)
   write(n0,*)
   write(n0,*)'LDL mean HB length Hi...Oj (Ang) = ',AVRHB_LDL
   write(n0,*)
   
   write(n16,'(a4)')'#LDL'
   do nkr = 1,NRDELS
      if(NDOOLDL(nkr) > 0)then
         RHB_ANG  = RHB_del*DFLOAT(nkr)          
         DISTR_HB = DFLOAT(NDOOLDL(nkr))/(RHB_del*DFLOAT(NOO_LDL))                    !Normalize
         WRITE(n16,39)RHB_ANG,DISTR_HB
      endif
   end do  
   AVRHB_LDL = AVROOLDL/DFLOAT(NOO_LDL)
   write(n0,*)
   write(n0,*)'LDL mean HB OO length Oi...Oj (Ang) = ',AVRHB_LDL
   write(n0,*)
   
   write(n17,'(a4)')'#LDL'
   do nkr = 1,ANGDELS
      if(NHBALDL(nkr) > 0)then
         RHB_ANG  = AHB_del*DFLOAT(nkr)          
         DISTR_HB = DFLOAT(NHBALDL(nkr))/(AHB_del*DFLOAT(NHBA_LDL))                    !Normalize
         WRITE(n17,39)RHB_ANG,DISTR_HB
      endif
   end do  
   AVRHB_LDL = AVRHBALDL/DFLOAT(NHBA_LDL)
   write(n0,*)
   write(n0,*)'LDL mean HB angle (< 30 deg) HOi...Oj (Ang) = ',AVRHB_LDL
   write(n0,*)
   
   
   
!HDL
   WRITE(n15,*)
   write(n15,'(a4)')'#HDL'
   do nkr = 1,NRDELS
      if(NDHBHDL(nkr) > 0)then
         RHB_ANG  = RHB_del*DFLOAT(nkr)          
         DISTR_HB = DFLOAT(NDHBHDL(nkr))/(RHB_del*DFLOAT(NHB_HDL))                     !Normalize
         WRITE(n15,39)RHB_ANG,DISTR_HB
      endif
   end do  
   AVRHB_HDL = AVRHBHDL/DFLOAT(NHB_HDL)
   write(n0,*)
   write(n0,*)'HDL mean HB length Hi...Oj (Ang) = ',AVRHB_HDL
   write(n0,*)              
   
   
   WRITE(n16,*)
   write(n16,'(a4)')'#HDL'
   do nkr = 1,NRDELS
      if(NDOOHDL(nkr) > 0)then
         RHB_ANG  = RHB_del*DFLOAT(nkr)          
         DISTR_HB = DFLOAT(NDOOHDL(nkr))/(RHB_del*DFLOAT(NOO_HDL))                    !Normalize
         WRITE(n16,39)RHB_ANG,DISTR_HB
      endif
   end do  
   AVRHB_HDL = AVROOHDL/DFLOAT(NOO_HDL)
   write(n0,*)
   write(n0,*)'HDL mean HB OO length Oi...Oj (Ang) = ',AVRHB_HDL
   write(n0,*)
   
   
   WRITE(n17,*)
   write(n17,'(a4)')'#HDL'
   do nkr = 1,ANGDELS
      if(NHBAHDL(nkr) > 0)then
         RHB_ANG  = AHB_del*DFLOAT(nkr)          
         DISTR_HB = DFLOAT(NHBAHDL(nkr))/(AHB_del*DFLOAT(NHBA_HDL))                    !Normalize
         WRITE(n17,39)RHB_ANG,DISTR_HB
      endif
   end do  
   AVRHB_HDL = AVRHBAHDL/DFLOAT(NHBA_HDL)
   write(n0,*)
   write(n0,*)'HDL mean HB angle (< 30 deg) HOi...Oj (Ang) = ',AVRHB_HDL
   write(n0,*)
endif           !End LHBdecv   



!Print mean orient. tcfs      
write(n2,9)
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a9)')'#Bulk-sHB'
   if(kr==2)write(n1,'(a9)')'#Bulk-lHB'
   if(kr==3)write(n1,'(a8)')'#LDL-sHB'
   if(kr==4)write(n1,'(a8)')'#LDL-lHB'
   if(kr==5)write(n1,'(a8)')'#HDL-sHB'
   if(kr==6)write(n1,'(a8)')'#HDL-lHB'
   if(kr==1)write(n2,'(a9)')'#Bulk-sHB'
   if(kr==2)write(n2,'(a9)')'#Bulk-lHB'
   if(kr==3)write(n2,'(a8)')'#LDL-sHB'
   if(kr==4)write(n2,'(a8)')'#LDL-lHB'
   if(kr==5)write(n2,'(a8)')'#HDL-sHB'
   if(kr==6)write(n2,'(a8)')'#HDL-lHB'
   if(kr==1)write(n3,'(a9)')'#Bulk-sHB'
   if(kr==2)write(n3,'(a9)')'#Bulk-lHB'
   if(kr==3)write(n3,'(a8)')'#LDL-sHB'
   if(kr==4)write(n3,'(a8)')'#LDL-lHB'
   if(kr==5)write(n3,'(a8)')'#HDL-sHB'
   if(kr==6)write(n3,'(a8)')'#HDL-lHB'
   
   do it = 1,NDELS
      jt = it -1                                                     !print tcf starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      WRITE(n1,19)PSECS,tcf_mean_1(it,kr)/dfloat(ko)     !mean orientational tcf rank 1 over different time-origins 
      WRITE(n2,19)PSECS,tcf_mean_2(it,kr)/dfloat(ko)     !mean orientational tcf rank 2 over different time-origins   
      WRITE(n3,19)PSECS,tcf_mean_3(it,kr)/dfloat(ko)     !mean orientational tcf rank 3 over different time-origins 
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
end do

!Calculate the relaxation times - tau - trapezoid method
do kr = kr_start,NRSHMAX
   if(kr==1)then
      write(n11,'(a9)')'#Bulk-sHB'
      write(n12,'(a9)')'#Bulk-sHB'
      write(n13,'(a9)')'#Bulk-sHB'
   elseif(kr==2)then
      write(n11,'(a9)')'#LDL-sHB'
      write(n12,'(a9)')'#LDL-sHB'
      write(n13,'(a9)')'#LDL-sHB'
   elseif(kr==3)then
      write(n11,'(a8)')'#LDL-sHB'
      write(n12,'(a8)')'#LDL-sHB'
      write(n13,'(a8)')'#LDL-sHB'
   elseif(kr==4)then
      write(n11,'(a8)')'#LDL-lHB'
      write(n12,'(a8)')'#LDL-lHB'
      write(n13,'(a8)')'#LDL-lHB'
   elseif(kr==5)then
      write(n11,'(a8)')'#HDL-sHB'
      write(n12,'(a8)')'#HDL-sHB'
      write(n13,'(a8)')'#HDL-sHB'
   elseif(kr==6)then
      write(n11,'(a8)')'#HDL-lHB'
      write(n12,'(a8)')'#HDL-lHB'
      write(n13,'(a8)')'#HDL-lHB'   
   endif  
   koc(kr)=ko
!   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
!   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n11,kr,dt,tcf_mean_1,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n12,kr,dt,tcf_mean_2,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
   call tcf_solv_tau(n13,kr,dt,tcf_mean_3,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
         
   if(kr==1)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==2)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==3)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==4)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==5)then
      write(n11,*)
      write(n12,*)
      write(n13,*)
   elseif(kr==6)then
      write(n11,*)
      write(n12,*)
      write(n13,*)   
   endif   
end do 

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n11)
close(n12)
close(n13)
close(n20)

deallocate(nOHorg,NOHidO,NOHidH)
deallocate(ntime)
deallocate(tcf_mean_1,tcf_mean_2,tcf_mean_3)
deallocate(norg_0HB,kv,nOHmean)
deallocate(koc)

deallocate(list_s)
deallocate(list_HB)
deallocate(list_LSI)
deallocate(LSI)
deallocate(xyz)
deallocate(x,y,z)

  9   FORMAT('#',17X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 39   FORMAT(27X,F9.4,5X,F14.5)
 
   return

END SUBROUTINE pure_wat_LDL_HDL_otcf_HB


 
!NG HBtcf routine - version 27
SUBROUTINE pure_wat_HB_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,hbtframe,hbs_samp)
! Calculate continuous, intermittent and survival water HB tcfs - pure water
! Here each water molecule is a particle 
! The maximum number of pairs of particles HB is therefore 1/2(Nw * Nw-1)
! The tcfs are calculated for every hbs_samp(fs) time-origin for a delay-time of hbtframe(ps)
! The survive tcfs probes HBs formed at each origin - for comparison purposes

!NEW VARIABLES INTRODUCTION
    integer,intent(in)                               :: npurew,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: hbtframe                       !HB tcfs delay time
    integer,intent(in)                               :: hbs_samp
 
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG    real(kind=4)                               :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                               :: cell(3)  
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                               :: TSTEP,PSECS
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    integer                                    :: NDELS,KINTVL,NDELS1,NDELS2
    
!END NEW VARIABLES INTRODUCTION

!LOCAL VARIABLES    
    integer      :: npairw,nwatoms,nwatoms1,nwatms
    integer      :: nbwmol,kmaxstep
    integer      :: i, j, k
    integer      :: n0, n1, n2, n3, n4, n5, n6
    integer      :: it, iH, ko, jt, it1
    integer      :: kpw
    integer      :: kt, nt, jw
!debug    
!    integer      :: n1000
!    real(kind=8),dimension(:),allocatable      :: tcf_cont_m_deb,tcf_interm_m_deb,tcf_surv_m_deb  
!end debug
    
    real(kind=8)                               :: rHBOO,rHBHO,aHBOHO                    ! HB definition
    real(kind=8)                               :: sqdrOO,drOO,sqdrHO,drHO,sqdriOH,driOH
    real(kind=8)                               :: sqarHO,arHO,sqariOH,ariOH
    real(kind=8)                               :: dcosang,dcosnm,dcosdn,angdohh
    real(kind=8)                               :: acosang,acosnm,acosdn,angaohh
    logical filex
    
!NEW LOCAL VARIABLES
    real(kind=8)                               :: r_lsi
    real(kind=8)                               :: LSI_del,AVLSI,LSI_ang2,DISTR_LSI
    integer                                    :: nsh,LSI_NDELS,NLSI_count,nkr
    real(kind=8)                               :: R_rdf_OO,RDEL,cubeh_,cube_,RADIUS,PI,VOLSHL,DSTNUM,FNOM
    integer                                    :: NSHLOO,NRDELS
    integer                                    :: natcheck
!END NEW LOCAL VARIABLES

! HB tcfs    
    real(kind=8),dimension(:),allocatable      :: tcf_cont_m,tcf_interm_m,tcf_surv_m  
    real(kind=8),dimension(:),allocatable      :: tcf_react           
    real(kind=8),dimension(:),allocatable      :: Pdist_surv,pdist_cont
    
!    real(kind=8)                               :: sum_h_int,sum_h_cont
    real(kind=8)                               :: sum_h0_int,sum_h0_cont,sum_h0_surv
!
    integer                                    :: norgmax
    real(kind=8),dimension(:,:),allocatable    :: tcf_hb_int,tcf_hb_cont,tcf_hb_surv
    real(kind=8),dimension(:,:),allocatable    :: ht0_i,ht0_c,ht0_s,ht_i,ht_c,ht_s
    integer,dimension(:),allocatable           :: ntime
    integer,dimension(:,:),allocatable         :: lhist_HB,lhist_s_HB
    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!NEW CODE LSI PARAMETERS AND XTC

!LSI values
r_lsi  = 3.70d0                                 !original value of the LSI

LSI_del = 0.0005d0                              !bin for LSI distribution
LSI_NDELS= 1.0d0/LSI_del

RDEL    = 0.05d0                                !bin size in rdf/distance distributions Ang
NRDELS  = 100.0d0/RDEL

PI=3.14159265d0
FNOM = dfloat(nmolwat)

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj

!END NEW CODE LSI PARAMETERS AND XTC

nbwmol = nmolwat
kmaxstep = nstep

! time-window for calculation of each tcf
NDELS = hbtframe*1000                              !sets the delay time-window in fs
norgmax = kmaxstep/hbs_samp + 1                    !maximum number of origins
if(norgmax.gt.kmaxstep)norgmax = kmaxstep
                 
if(NDELS.ge.kmaxstep)NDELS = NDELS - 1
NDELS1 = NDELS - 1
NDELS2 = NDELS - 2

rHBOO  = 3.5
rHBHO  = 2.45
aHBOHO = 30.0

!TSTEP - TIME-STEP IN ps
TSTEP  = dt*1.0e-3
!KINTVL - INTERVAL BETWEEN TIME-STEPS
KINTVL = 1

nwatoms = nbwmol*3
nwatoms1 = (nbwmol-1)*3
npairw = nbwmol*(nbwmol-1)/2

allocate(tcf_hb_int(norgmax,NDELS),tcf_hb_cont(norgmax,NDELS),tcf_hb_surv(norgmax,NDELS))
allocate(ht0_i(norgmax,npairw),ht0_c(norgmax,npairw),ht0_s(norgmax,npairw))
allocate(ht_i(norgmax,npairw),ht_c(norgmax,npairw),ht_s(norgmax,npairw))
allocate(ntime(norgmax))
allocate(lhist_HB(norgmax,npairw),lhist_s_HB(0:norgmax,npairw))

tcf_hb_int = 0.0; tcf_hb_cont = 0.0; tcf_hb_surv = 0.0
ht0_i = 0.0; ht0_c = 0.0; ht0_s = 0.0; ht_i = 0.0; ht_c = 0.0; ht_s = 0.0
ntime = 1
lhist_HB = 1
lhist_s_HB = 1

!HB tcfs

allocate(tcf_cont_m(NDELS),tcf_interm_m(NDELS),tcf_surv_m(NDELS))
allocate(tcf_react(NDELS),Pdist_surv(NDELS),pdist_cont(NDELS))

!debug   allocate(tcf_cont_m_deb(NDELS),tcf_interm_m_deb(NDELS),tcf_surv_m_deb(NDELS))

tcf_cont_m = 0.0;tcf_interm_m = 0.0;tcf_surv_m = 0.0
tcf_react = 0.0; Pdist_surv = 0.0; pdist_cont = 0.0
!sum_h_int = 0.0;sum_h_cont = 0.0
sum_h0_int = 0.0;sum_h0_cont = 0.0;sum_h0_surv = 0.0

!debug   tcf_cont_m_deb = 0.0;tcf_interm_m_deb = 0.0;tcf_surv_m_deb = 0.0

write(*,*)
write(*,*)'Pure water HB tcf calculation'
write(*,*)'Routine pure_wat_HB_tcfs'
write(*,*)'water HB tcf delay time-window (ps) = ',hbtframe
write(*,*)'water HB sampling frequency (fs) = ',hbs_samp
write(*,*)'results are printed to pureW_HB_tcf_out'
write(*,*)

inquire(file='pureW_HB_tcf_out/log_HB_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir pureW_HB_tcf_out')

!output
n0 = 50
n1 = 100
n2 = 200
n3 = 300
n4 = 400
n5 = 500
n6 = 600

!debug
!n1000 = 1000

open(n0,file='pureW_HB_tcf_out/log_HB_tcf.dat',status='unknown',action='write')
open(n1,file='pureW_HB_tcf_out/HBtcf_cont.dat',status='unknown',action='write')
open(n2,file='pureW_HB_tcf_out/HBtcf_interm.dat',status='unknown',action='write')
open(n3,file='pureW_HB_tcf_out/HBtcf_reactive.dat',status='unknown',action='write')
open(n4,file='pureW_HB_tcf_out/HBtcf_survive.dat',status='unknown',action='write')
open(n5,file='pureW_HB_tcf_out/HB_dist_P_surv_t.dat',status='unknown',action='write')
open(n6,file='pureW_HB_tcf_out/HB_dist_p_cont_t.dat',status='unknown',action='write')

!open(n1000,file='pureW_HB_tcf_out/debug_HBtcf_cont.dat',status='unknown',action='write')

write(n0,*)
write(n0,*)'Pure water HB tcf calculation'
write(n0,*)'Routine pure_wat_HB_tcfs'
write(n0,*)'water HB tcf delay time-window (ps) = ',hbtframe
write(n0,*)'water HB sampling frequency (fs) = ',hbs_samp
write(n0,*)'results are printed to pureW_HB_tcf_out'
write(n0,*)

write(*,*) 
write(*,*)'Starting <h(0)h(t)> calculation'
write(*,*) 
write(n0,*) 
write(n0,*)'Starting <h(0)h(t)> calculation'
write(n0,*) 

!NEW READ xyz

ko = 0                                                     !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep                                                
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif      
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!End Version 9      
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   
    
    cubeh_ = cell(1)/2.0d0
    cube_  = cell(1)

   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 



!END NEW READ xyz

   if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))ko = ko + 1  
!   if(mod(j,5000)==0)write(*,*)'starting time-step ',j

   kpw = 0                                        !number of pairs of particles counter
   
! Alternative loop 
!NG        io = 0
!NG        do i = 1,natms,nwatsites                    !Loop over oxygen atoms 1,4,7,10,...   
!NG           jo = 0                                   !oxygen atom counter 1,2,3,...   
!NG           io = io + 1    
!NG           do k = 1,natms,nwatsites                 !Loop over oxygen atoms 1,4,7,10,...   
!NG              jo = jo + 1    
!NG              IF(jo.NE.io)THEN   
! End alternative loop

!sample all possible pairs - WARNING!!! this loop does not distiguish between ik and ki 
!
   do i=1,nwatoms1,nwatsites                          !nwatoms1 = (nbwmol-1)*3
      do k=i+nwatsites,nwatoms,nwatsites              !Loop over water oxygens (proton acceptors)         
         kpw = kpw + 1                                !kpw must run over the maximum number of pairs of particles
                                                      !WARNING!!! a water molecule (not OH) is considered a single particle                                                                              
! intermolecular OO distance 
         dx = x(i) - x(k)
         dy = y(i) - y(k)
         dz = z(i) - z(k)
         dx = dx - anint(dx/cell(1))*cell(1)
         dy = dy - anint(dy/cell(2))*cell(2)
         dz = dz - anint(dz/cell(3))*cell(3)
         dr2 = dx**2.0 + dy**2.0 + dz**2.0
         dr  = sqrt(dr2)
         sqdrOO = dr2
         drOO = dr
         if(drOO.gt.rHBOO)then
            if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then      !New origin
               ht0_i(ko,kpw) = 0.0                                              !at origin ko pair kpw is not HBonded
               ht0_c(ko,kpw) = 0.0
               ht0_s(ko,kpw) = 0.0
               sum_h0_int  = sum_h0_int  + ht0_i(ko,kpw)
               sum_h0_cont = sum_h0_cont + ht0_c(ko,kpw)
               sum_h0_surv = sum_h0_surv + ht0_s(ko,kpw)
            endif
            do kt = 1,ko                                                        !New instantaneous value for each origin
               ht_i(kt,kpw) = 0.0
               ht_c(kt,kpw) = 0.0
               ht_s(kt,kpw) = 0.0
               lhist_HB(kt,kpw) = 0
               lhist_s_HB(kt,kpw) = 0
            end do   
!            sum_h_int  = sum_h_int  + 0.0
!            sum_h_cont = sum_h_cont + 0.0              
                                                                                 !end if/then/else - jump to next water
         elseif(drOO.le.rHBOO)then     
            do iH = 1,2                   !check two hydrogen atoms per water
! if the first H is HBonded leave the loop over the H atoms because the second H will not be HBonded to the same O-P                
! intermolecular H...O distance  (donor)                 
               dx = x(i+iH) - x(k)
               dy = y(i+iH) - y(k)
               dz = z(i+iH) - z(k)
               dx = dx - anint(dx/cell(1))*cell(1)
               dy = dy - anint(dy/cell(2))*cell(2)
               dz = dz - anint(dz/cell(3))*cell(3)
               dr2 = dx**2.0 + dy**2.0 + dz**2.0
               dr  = sqrt(dr2)
               sqdrHO = dr2
               drHO = dr
! intermolecular O...H distance  (acceptor)                 
               dx = x(k+iH) - x(i)
               dy = y(k+iH) - y(i)
               dz = z(k+iH) - z(i)
               dx = dx - anint(dx/cell(1))*cell(1)
               dy = dy - anint(dy/cell(2))*cell(2)
               dz = dz - anint(dz/cell(3))*cell(3)
               dr2 = dx**2.0 + dy**2.0 + dz**2.0
               dr  = sqrt(dr2)
               sqarHO = dr2
               arHO = dr   
!check         if(arHO < rHBHO .and. drHO < rHBHO) write(*,*)'Warning!!! water pair "is" r(O...H) mutual donor and acceptor'  
               if((drHO.gt.rHBHO).and.(arHO.gt.rHBHO))then
!               if(drHO.gt.rHBHO)then
                  if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then     
                     if(iH==2)ht0_i(ko,kpw) = 0.0                                    !at origin ko pair kpw is not HBonded
                     if(iH==2)ht0_c(ko,kpw) = 0.0
                     if(iH==2)ht0_s(ko,kpw) = 0.0
                     if(iH==2)sum_h0_int  = sum_h0_int  + ht0_i(ko,kpw)
                     if(iH==2)sum_h0_cont = sum_h0_cont + ht0_c(ko,kpw)
                     if(iH==2)sum_h0_surv = sum_h0_surv + ht0_s(ko,kpw)
                  endif
                  if(iH==2)then              !if the first H is not HB try the second - only if neither is HB assign 0 to the water pair
                     do kt = 1,ko                             !New instantaneous value for each origin
                        ht_i(kt,kpw) = 0.0
                        ht_c(kt,kpw) = 0.0
                        ht_s(kt,kpw) = 0.0
                        lhist_HB(kt,kpw) = 0
                        lhist_s_HB(kt,kpw) = 0
                     end do
                  endif   
!                  if(iH==2)sum_h_int  = sum_h_int  + 0.0
!                  if(iH==2)sum_h_cont = sum_h_cont + 0.0                    
               elseif((drHO.le.rHBHO).or.(arHO.le.rHBHO))then
! intramolecular OH distance (donor) 
                     dx = x(i+iH) - x(i)
                     dy = y(i+iH) - y(i)
                     dz = z(i+iH) - z(i)
                     dx = dx - anint(dx/cell(1))*cell(1)
                     dy = dy - anint(dy/cell(2))*cell(2)
                     dz = dz - anint(dz/cell(3))*cell(3)
                     dr2 = dx**2.0 + dy**2.0 + dz**2.0
                     dr  = sqrt(dr2)
                     sqdriOH = dr2
                     driOH = dr 
! intramolecular OH distance (acceptor) 
                     dx = x(k+iH) - x(k)
                     dy = y(k+iH) - y(k)
                     dz = z(k+iH) - z(k)
                     dx = dx - anint(dx/cell(1))*cell(1)
                     dy = dy - anint(dy/cell(2))*cell(2)
                     dz = dz - anint(dz/cell(3))*cell(3)
                     dr2 = dx**2.0 + dy**2.0 + dz**2.0
                     dr  = sqrt(dr2)
                     sqariOH = dr2
                     ariOH = dr 
!
!Apply the COSINE rule to calculate the HB donor angle
!donnor
!                 DCOSNM = SQDHO-SQDOH-SQDOO
                  dcosnm = sqdrHO-sqdriOH-sqdrOO                   !donnor(i)
!                 ACOSNM = SQAHO-SQAOH-SQDOO
                  acosnm = sqarHO-sqariOH-sqdrOO                   !acceptor(i)
!donor                  
!                 DCOSDN = -2.0D0*DOHYD*DOXYG
                  dcosdn = -2.0*driOH*drOO
!                 DCOSANG = DCOSNM/DCOSDN
                  dcosang = dcosnm/dcosdn
!                 ANGDOHH = DACOSD(DCOSANG)
!acceptor                
!                 ACOSDN = -2.0D0*AOHYD*DOXYG
                  acosdn = -2.0*ariOH*drOO
!                 ACOSANG = ACOSNM/ACOSDN
                  acosang = acosnm/acosdn
!                 ANGAOHH = DACOSD(ACOSANG)

!acosd - returns arcosin in degrees                - donor
                  if(dcosang<-1.0)then
                     write(*,*)'Warning HB angle',dcosang
                     dcosang = -1.0
                     write(*,*)'HB angle = ',acosd(dcosang)
                     write(*,*)
                  endif
                  if(dcosang>+1.0)then
                     write(*,*)'Warning HB angle',dcosang
                     dcosang = +1.0
                     write(*,*)'HB angle = ',acosd(dcosang)
                     write(*,*)
                  endif
                  angdohh = acosd(dcosang)
!acosd - returns arcosin in degrees                  - acceptor
                  if(acosang<-1.0)then
                     write(*,*)'Warning HB angle',acosang
                     acosang = -1.0
                     write(*,*)'HB angle = ',acosd(acosang)
                     write(*,*)
                  endif
                  if(acosang>+1.0)then
                     write(*,*)'Warning HB angle',acosang
                     acosang = +1.0
                     write(*,*)'HB angle = ',acosd(acosang)
                     write(*,*)
                  endif
                  angaohh = acosd(acosang)                   
!The water pair is not HBonded                                    
                  if((angdohh.ge.aHBOHO).and.(angaohh.ge.aHBOHO))then                 
                     if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then     
                        if(iH==2)ht0_i(ko,kpw) = 0.0
                        if(iH==2)ht0_c(ko,kpw) = 0.0
                        if(iH==2)ht0_s(ko,kpw) = 0.0
                        if(iH==2)sum_h0_int  = sum_h0_int  + ht0_i(ko,kpw)
                        if(iH==2)sum_h0_cont = sum_h0_cont + ht0_c(ko,kpw)
                        if(iH==2)sum_h0_surv = sum_h0_surv + ht0_s(ko,kpw)
                     endif
                     if(iH==2)then
                        do kt = 1,ko
                           ht_i(kt,kpw) = 0.0
                           ht_c(kt,kpw) = 0.0
                           ht_s(kt,kpw) = 0.0
                           lhist_HB(kt,kpw) = 0
                           lhist_s_HB(kt,kpw) = 0
                        end do
                     endif   
!                     if(iH==2)sum_h_int  = sum_h_int  + 0.0
!                     if(iH==2)sum_h_cont = sum_h_cont + 0.0          
!The water pair is HBonded                        
                  elseif((angdohh.lt.aHBOHO).or.(angaohh.lt.aHBOHO))then
                     if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then     
                        ht0_i(ko,kpw) = 1.0
                        ht0_c(ko,kpw) = 1.0
!check                        if(lhist_HB(ko,kpw)==0)write(*,*)'Error!!!'
!if there was no HB in the previous time-step set ht0_s = 1 - newly formed HB
!if there was an HB in the previous time-step let ht0_s = 0 - count only HBs from the moment they formed
                        if(lhist_s_HB(ko-1,kpw)==0)then                 ! at j == 1 no HBs are possible because lhist_s_HB(0,kpw)==1 - subtract one origin from average
                           ht0_s(ko,kpw) = 1.0                          ! HB formed at this time origin 
                        endif   
                        sum_h0_int  = sum_h0_int  + ht0_i(ko,kpw)
                        sum_h0_cont = sum_h0_cont + ht0_c(ko,kpw)
                        sum_h0_surv = sum_h0_surv + ht0_s(ko,kpw)
                     endif                    
                     do kt = 1,ko
!HB intermittent                     
                        ht_i(kt,kpw) = 1.0
!HB continuous - HBs "present" at time t0 - average HB                           
                        if(lhist_HB(kt,kpw)==1)then                     ! lhist_HB is update every time-step - it must be continuously 1
                           ht_c(kt,kpw) = 1.0
                        else
                           ht_c(kt,kpw) = 0.0
                        endif   
!HB survival - HBs "formed" at t0 - newly generated HB                           
                        if(lhist_s_HB(kt,kpw)==1)then
                           ht_s(kt,kpw) = 1.0
                        else
                           ht_s(kt,kpw) = 0.0
                        endif
                     end do   
!                        
!                     sum_h_int  = sum_h_int  + 1.0
!                     if(lhist_HB(ko,kpw)==1)then
!                        sum_h_cont = sum_h_cont + 1.0
!                     else
!                        sum_h_cont = sum_h_cont + 0.0
!                     endif
!End Version 34                         
                     goto 5                                        !leave the loop over the H atoms because the second H will not be HB to the same O-P
                  endif  ! HB angle fulfilled                       
               endif     ! HB H...O distance fulfilled               
            end do       !end loop over H atoms                  
         endif           !HB O...O distanced fulfilled    
         5 continue
      end do                !end loop over water oxygens
   end do                   !end loop over water oxygens

!calculate tcfs - for each origin sum over different pairs
   do kt = 1,ko                                   !loop over time origins
      if(ntime(kt).le.NDELS)then                  !ntime counts the delay times for which the tcf has already been calculated
         nt = ntime(kt)                           !ntime(kt) initialized to 1
         do i = 1,npairw      
            tcf_hb_int(kt,nt)  = tcf_hb_int(kt,nt)  + ht0_i(kt,i)*ht_i(kt,i)
            tcf_hb_cont(kt,nt) = tcf_hb_cont(kt,nt) + ht0_c(kt,i)*ht_c(kt,i)
            tcf_hb_surv(kt,nt) = tcf_hb_surv(kt,nt) + ht0_s(kt,i)*ht_s(kt,i) 
         end do   
         ntime(kt) = ntime(kt) + 1                                !delay time counter
      endif
   end do                                      !end loop over time-origins
!Version 35 end

end do                      !end loop over origins   

!sum_h_int  = sum_h_int/(float(npairw)*float(kmaxstep))
!sum_h_cont = sum_h_cont/(float(npairw)*float(kmaxstep))   
sum_h0_int  = sum_h0_int/(float(npairw)*float(ko))
sum_h0_cont = sum_h0_cont/(float(npairw)*float(ko))
sum_h0_surv = sum_h0_surv/(float(npairw)*float(ko-1))                                    !subtract one origin

write(*,*)'Number of pairs possible = ',npairw
write(*,*)'Number of pairs counted = ',kpw
write(*,*)'Number of steps = ',kmaxstep
write(*,*)'Number of origins possible = ',norgmax
write(*,*)'Number of origins counted = ',ko
write(*,*)
!write(*,*)'intermittent <h> = ',sum_h_int
!write(*,*)'continuous <h> = ',sum_h_cont
write(*,*)'intermittent <h0> = ',sum_h0_int
write(*,*)'continuous <h0> = ',sum_h0_cont
write(*,*)'survival <h0> = ',sum_h0_surv
write(*,*)

write(n0,*)'Number of pairs possible = ',npairw
write(n0,*)'Number of pairs counted = ',kpw
write(n0,*)'Number of steps = ',kmaxstep
write(n0,*)'Number of origins possible = ',norgmax
write(n0,*)'Number of origins used = ',ko
write(n0,*)
!write(n0,*)'intermittent <h> = ',sum_h_int
!write(n0,*)'continuous <h> = ',sum_h_cont
write(n0,*)'intermittent <h0> = ',sum_h0_int
write(n0,*)'continuous <h0> = ',sum_h0_cont
write(n0,*)'survival <h0> = ',sum_h0_surv
write(n0,*)

!Average over number of pairs of water molecules - sum different origins tcfs
do kt = 1,ko                                               !loop over origins
   do it = 1,NDELS           
      tcf_interm_m(it) = tcf_interm_m(it) + tcf_hb_int(kt,it)/float(npairw) 
      tcf_cont_m(it)   = tcf_cont_m(it)   + tcf_hb_cont(kt,it)/float(npairw)   
      tcf_surv_m(it)   = tcf_surv_m(it)   + tcf_hb_surv(kt,it)/float(npairw) 
   end do
end do

!debug 
!tcf_cont_m_deb = tcf_cont_m
!end debug


!Divide by number of origins and normalize
do it = 1,NDELS  
   tcf_interm_m(it)=tcf_interm_m(it)/(float(ko)*sum_h0_int)
   tcf_cont_m(it)  =tcf_cont_m(it)/(float(ko)*sum_h0_cont)   
   tcf_surv_m(it)  =tcf_surv_m(it)/(float(ko-1)*sum_h0_surv) 
end do 
!End Version 34

!calculate and print mean tcfs - average over number of origins     
      write(n1,9)
      write(n2,9)
      write(n4,9)
      do it = 1,NDELS
         jt = it -1                                              !print tcf starting from time zero
         PSECS = TSTEP*float(jt)*float(KINTVL)
         WRITE(n1,19)PSECS,tcf_cont_m(it)   
         WRITE(n2,19)PSECS,tcf_interm_m(it)  
         write(n4,19)PSECS,tcf_surv_m(it)  
!debug         
!         write(n1000)PSECS,tcf_cont_m_deb(it)
      end do
      WRITE(n1,*)
      WRITE(n2,*)
      WRITE(n4,*)
!      
!Luzar reactive flux      
!      do it = 2,NDELS1
!         tcf_react(it) = -(tcf_interm_m(it+1)-tcf_interm_m(it-1))/(2.0*TSTEP)
!      end do

!Calculate the reactive flux and the distribution function every 2 time-steps
      do it = 1,NDELS2,2
         tcf_react(it)  = -(tcf_interm_m(it+2)-tcf_interm_m(it))/(2.0*TSTEP)
         Pdist_surv(it) = -(tcf_surv_m(it+2)-tcf_surv_m(it))/(2.0*TSTEP)
         pdist_cont(it) = -(tcf_cont_m(it+2)-tcf_cont_m(it))/(2.0*TSTEP)
      end do
!print reactive flux and P(t)
      write(n3,9) 
      do it = 1,NDELS2,2
         jt = it -1                                              !print tcf starting from time zero
         PSECS = TSTEP*float(jt)*float(KINTVL)
         WRITE(n3,19)PSECS,tcf_react(it)
         WRITE(n5,19)PSECS,Pdist_surv(it)
         WRITE(n6,19)PSECS,pdist_cont(it)
      end do
      write(n3,*)
      write(n5,*)
      write(n6,*)
!
!Calculate the reactive flux and the distribution function every time-step
      do it = 1,NDELS1
         tcf_react(it)  = -(tcf_interm_m(it+1)-tcf_interm_m(it))/TSTEP
         Pdist_surv(it) = -(tcf_surv_m(it+1)-tcf_surv_m(it))/TSTEP
         pdist_cont(it) = -(tcf_cont_m(it+1)-tcf_cont_m(it))/TSTEP
      end do
!print reactive flux and P(t)
      do it = 1,NDELS1
         jt = it -1                                              !print tcf starting from time zero
         PSECS = TSTEP*float(jt)*float(KINTVL)
         WRITE(n3,19)PSECS,tcf_react(it)
         WRITE(n5,19)PSECS,Pdist_surv(it)
         WRITE(n6,19)PSECS,pdist_cont(it)
      end do    
      write(n3,*)
      write(n5,*)
      write(n6,*)

deallocate(tcf_hb_int,tcf_hb_cont,tcf_hb_surv)
deallocate(ht0_i,ht0_c,ht0_s)
deallocate(ht_i,ht_c,ht_s)
deallocate(ntime)
deallocate(lhist_HB,lhist_s_HB)
deallocate(tcf_cont_m,tcf_interm_m,tcf_surv_m)
deallocate(tcf_react,Pdist_surv,pdist_cont)
      
  9   FORMAT('#',11X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
!    return

END SUBROUTINE pure_wat_HB_tcfs

!NG end 


           
SUBROUTINE pure_wat_HB_LDL_HDL_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,hbtframe,hbs_samp,Nenv_hb_HB,&
                                     LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_)
! Calculate continuous, intermittent and survival water HB tcfs - pure water
! Here each water molecule is a particle 
! The maximum number of pairs of particles HB is therefore 1/2(Nw * Nw-1)
! The tcfs are calculated for every hbs_samp(fs) time-origin for a delay-time of hbtframe(ps)
! The survive tcfs probes HBs formed at each origin - for comparison purposes
! Calculate LSI and probe the HB tcfs for three environemnts of pairs of molecules HBonded
! Both molecules forming an HB must be in the same environment to be accepted otherwise they are discarded from environments 1 and 2

!NEW VARIABLES INTRODUCTION
    integer,intent(in)                               :: npurew,nbox,Nenv_hb_HB
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: hbtframe                       !HB tcfs delay time
    integer,intent(in)                               :: hbs_samp
    real(kind=8),intent(in)                          :: LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_
 
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG    real(kind=4)                               :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                               :: cell(3)  
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                               :: TSTEP,PSECS
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    integer                                    :: NDELS,KINTVL,NDELS1,NDELS2
    
!END NEW VARIABLES INTRODUCTION

!LOCAL VARIABLES    
    integer      :: npairw,nwatoms,nwatoms1,nwatms
    integer      :: nbwmol,kmaxstep
    integer      :: i, j, k
    integer      :: n0, n1, n2, n3, n4, n5, n6
    integer      :: n10,n11
    integer      :: it, iH, ko, jt, it1
    integer      :: kpw
    integer      :: kt, nt, jw, kw
    integer      :: io,ihb,jhb,jo,kr
    integer      :: nienv, nkenv               ! molecule i and k LSI environment

    real(kind=8)                               :: rHBOO,rHBHO,aHBOHO                    ! HB definition
    real(kind=8)                               :: sqdrOO,drOO,sqdrHO,drHO,sqdriOH,driOH
    real(kind=8)                               :: sqarHO,arHO,sqariOH,ariOH
    real(kind=8)                               :: dcosang,dcosnm,dcosdn,angdohh
    real(kind=8)                               :: acosang,acosnm,acosdn,angaohh
    logical filex
    
! LSI VARIABLES

    integer,dimension(:),allocatable             :: list_s
    integer,dimension(:),allocatable             :: list_LSI               !LSI index - LSI computed locally to compare with Gaia List
    integer                                      :: pwnmolsol,pwnatmsol,pwnions                !solution parameters - set to 0
    real(kind=8),dimension(:),allocatable        :: LSI
    integer,dimension(:),allocatable             :: NLSI
    real(kind=8)                                 :: r_lsi
    real(kind=8)                                 :: LSI_del,AVLSI,LSI_ang2,DISTR_LSI
    integer                                      :: nsh,LSI_NDELS,NLSI_count,nkr
!    integer                                      :: nsh_
    real(kind=8)                                 :: R_rdf_OO,RDEL,cubeh_,cube_,RADIUS,PI,VOLSHL,DSTNUM,FNOM
    integer                                      :: NSHLOO,NRDELS
    integer                                      :: natcheck
    integer,dimension(:),allocatable             :: NRDFOO,NRDFOO_LDL,NRDFOO_HDL
    real(kind=8),dimension(:),allocatable        :: GROO,GROO_LDL,GROO_HDL
    integer                                      :: nradsh,nradsh_
    logical                                      :: LRDF
    integer,dimension(:,:),allocatable           :: k_lsi
    integer,dimension(:),allocatable             :: k_HB_org
    integer,dimension(:,:),allocatable           :: wat_id_nradsh
    integer                                      :: korg_0HB
!END LSI VARIABLES

! HB tcfs    
    real(kind=8),dimension(:,:),allocatable    :: tcf_cont_m,tcf_interm_m,tcf_surv_m  
    real(kind=8),dimension(:,:),allocatable    :: tcf_react           
    real(kind=8),dimension(:,:),allocatable    :: Pdist_surv,pdist_cont
    
!    real(kind=8),dimension(:),allocatable      :: sum_h_int,sum_h_cont
    real(kind=8),dimension(:,:),allocatable    :: sum_h0_int,sum_h0_cont,sum_h0_surv
    real(kind=8),dimension(:),allocatable      :: sum_h0_int_,sum_h0_cont_,sum_h0_surv_
!
    integer                                    :: norgmax
    integer                                    :: NRSHMAX
    real(kind=8),dimension(:,:,:),allocatable  :: tcf_hb_int,tcf_hb_cont,tcf_hb_surv
    real(kind=8),dimension(:,:,:),allocatable  :: ht0_i,ht0_c,ht0_s 
    real(kind=8),dimension(:,:),allocatable    :: ht_i,ht_c,ht_s
    integer,dimension(:),allocatable           :: ntime
    integer,dimension(:,:),allocatable         :: lhist_HB,lhist_s_HB
    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!NEW CODE LSI PARAMETERS AND XTC

!LSI values
r_lsi  = 3.70d0                                 !original value of the LSI

LSI_del = 0.0005d0                              !bin for LSI distribution
LSI_NDELS= 1.0d0/LSI_del

RDEL    = 0.05d0                                !bin size in rdf/distance distributions Ang
NRDELS  = 100.0d0/RDEL
LRDF    =.false.
PI=3.14159265d0
FNOM = dfloat(nmolwat)

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj

!END NEW CODE LSI PARAMETERS AND XTC

nbwmol = nmolwat
kmaxstep = nstep

! time-window for calculation of each tcf
NDELS = hbtframe*1000                              !sets the delay time-window in fs
norgmax = kmaxstep/hbs_samp + 1                    !maximum number of origins
if(norgmax.gt.kmaxstep)norgmax = kmaxstep
                 
if(NDELS.ge.kmaxstep)NDELS = NDELS - 1
NDELS1 = NDELS - 1
NDELS2 = NDELS - 2
                               
NRSHMAX = Nenv_hb_HB                                   !maximum number of environments = 3 [bulk environment not calculated in this routine]
NRSHMAX = 3
!check purposes NRSHMAX = 1                                            !WARNING!!! SET THIS TO 1 TO CHECK THE ROUTINE FOR A SINGLE ENVIRONMENT

!HB definition
rHBOO  = 3.5
rHBHO  = 2.45
aHBOHO = 30.0
!HB definition

!TSTEP - TIME-STEP IN ps
TSTEP  = dt*1.0e-3
!KINTVL - INTERVAL BETWEEN TIME-STEPS
KINTVL = 1

nwatoms = nbwmol*3
nwatoms1 = (nbwmol-1)*3
npairw = nbwmol*(nbwmol-1)/2

allocate(tcf_hb_int(norgmax,NDELS,NRSHMAX),tcf_hb_cont(norgmax,NDELS,NRSHMAX),tcf_hb_surv(norgmax,NDELS,NRSHMAX))
allocate(ht0_i(norgmax,npairw,NRSHMAX),ht0_c(norgmax,npairw,NRSHMAX),ht0_s(norgmax,npairw,NRSHMAX))
allocate(ht_i(norgmax,npairw),ht_c(norgmax,npairw),ht_s(norgmax,npairw))
allocate(wat_id_nradsh(norgmax,npairw))

allocate(ntime(norgmax))
allocate(lhist_HB(norgmax,npairw),lhist_s_HB(0:norgmax,npairw))
allocate(k_lsi(norgmax,NRSHMAX),k_HB_org(norgmax))

tcf_hb_int = 0.0; tcf_hb_cont = 0.0; tcf_hb_surv = 0.0
ht0_i = 0.0; ht0_c = 0.0; ht0_s = 0.0; ht_i = 0.0; ht_c = 0.0; ht_s = 0.0
ntime = 1
lhist_HB = 1
lhist_s_HB = 1
k_lsi = 0; k_HB_org = 0

!HB tcfs

allocate(tcf_cont_m(NDELS,NRSHMAX),tcf_interm_m(NDELS,NRSHMAX),tcf_surv_m(NDELS,NRSHMAX))
allocate(tcf_react(NDELS,NRSHMAX),Pdist_surv(NDELS,NRSHMAX),pdist_cont(NDELS,NRSHMAX))
allocate(sum_h0_int(norgmax,NRSHMAX),sum_h0_cont(norgmax,NRSHMAX),sum_h0_surv(norgmax,NRSHMAX))
allocate(sum_h0_int_(NRSHMAX),sum_h0_cont_(NRSHMAX),sum_h0_surv_(NRSHMAX))
!allocate(sum_h_int(NRSHMAX),sum_h_cont(NRSHMAX))

tcf_cont_m = 0.0;tcf_interm_m = 0.0;tcf_surv_m = 0.0
tcf_react = 0.0; Pdist_surv = 0.0; pdist_cont = 0.0
!sum_h_int = 0.0;sum_h_cont = 0.0
sum_h0_int = 0.0;sum_h0_cont = 0.0;sum_h0_surv = 0.0
sum_h0_int_ = 0.0;sum_h0_cont_ = 0.0;sum_h0_surv_ = 0.0

!List HDL/LDL
allocate(list_s(nmolwat))
allocate(list_LSI(nmolwat))
allocate(LSI(nmolwat))
allocate(NLSI(LSI_NDELS))

allocate(NRDFOO(NRDELS),NRDFOO_LDL(NRDELS),NRDFOO_HDL(NRDELS))
allocate(GROO(NRDELS),GROO_LDL(NRDELS),GROO_HDL(NRDELS))

LSI = 0
NLSI = 0
AVLSI = 0.0d0
pwnmolsol = 0
pwnatmsol = 0
pwnions   = 0  
NLSI_count = 0

!rdf
NRDFOO     = 0
NRDFOO_LDL = 0
NRDFOO_HDL = 0
!End List HDL/LDL

IF(NRSHMAX==1)Write(*,*)'WARNING!!! NUMBER OF ENVIRONMENTS = ', NRSHMAX

write(*,*)
write(*,*)'Pure water HB tcf calculation - LSI environments'
write(*,*)'Routine pure_wat_HB_LDL_HDL_tcfs'
write(*,*)'water HB tcf delay time-window (ps) = ',hbtframe
write(*,*)'water HB sampling frequency (fs) = ',hbs_samp
write(*,*)'results are printed to pureW_HB_tcf_out_LSI'
write(*,*)

inquire(file='pureW_HB_tcf_out_LSI/log_HB_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir pureW_HB_tcf_out_LSI')

!output
n0 = 50
n1 = 100
n2 = 200
n3 = 300
n4 = 400
n5 = 500
n6 = 600
n10 = 1000
n11 = 1100 

open(n0,file='pureW_HB_tcf_out_LSI/log_HB_tcf.dat',status='unknown',action='write')
open(n1,file='pureW_HB_tcf_out_LSI/HBtcf_cont.dat',status='unknown',action='write')
open(n2,file='pureW_HB_tcf_out_LSI/HBtcf_interm.dat',status='unknown',action='write')
open(n3,file='pureW_HB_tcf_out_LSI/HBtcf_reactive.dat',status='unknown',action='write')
open(n4,file='pureW_HB_tcf_out_LSI/HBtcf_survive.dat',status='unknown',action='write')
open(n5,file='pureW_HB_tcf_out_LSI/HB_dist_P_surv_t.dat',status='unknown',action='write')
open(n6,file='pureW_HB_tcf_out_LSI/HB_dist_p_cont_t.dat',status='unknown',action='write')

open(n10,file='pureW_HB_tcf_out_LSI/LSI_dist_HB.dat',status='unknown',action='write')
open(n11,file='pureW_HB_tcf_out_LSI/rdf_OO_HB.dat',status='unknown',action='write')

IF(NRSHMAX==1)Write(n0,*)'WARNING!!! NUMBER OF ENVIRONMENTS = ', NRSHMAX

write(n0,*)
write(n0,*)'Pure water HB tcf calculation - LSI environments'
write(n0,*)'Routine pure_wat_HB_LDL_HDL_tcfs'
write(n0,*)'water HB tcf delay time-window (ps) = ',hbtframe
write(n0,*)'water HB sampling frequency (fs) = ',hbs_samp
write(n0,*)'results are printed to pureW_HB_tcf_out_LSI'
write(n0,*)

write(*,*) 
write(*,*)'Starting <h(0)h(t)> calculation'
write(*,*) 
write(n0,*) 
write(n0,*)'Starting <h(0)h(t)> calculation'
write(n0,*) 

!NEW READ xyz

ko = 0                                                               !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep                                                
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif      
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!End Version 9      
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   
    
    cubeh_ = cell(1)/2.0d0
    cube_  = cell(1)

   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!END NEW READ xyz

!   if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))ko = ko + 1                       !origins counter
!   if(mod(j,5000)==0)write(*,*)'starting time-step ',j

   kpw = 0                    !number of pairs of particles counter                                                                  
   
   if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then  
      ko = ko + 1                                                                         !origins counter      
!Calculate LSI ========================================================================================================
      if(NRSHMAX > 1)then                                                                 !to check this routine for a single environment
         io = 0                                                                           !water molecule number 1,2,3,4   ihb = 1,4,7,10
         do ihb = 1,natms,nwatsites
            LSI = 0
            io = io + 1                                                                    
!            call LSI_water(ihb,io,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,LSI,r_lsi)
             call LSI_water_fast(ihb,io,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,LSI,r_lsi)
!check LSI - compare with list_s(i) based on Gaia & Lars LSI cut-offs for LDL (> 0.17 Ang**2) and HDL (< 0.01 Ang**2) species            
            if(LSI(io) > LSI_min_1_ .and. LSI(io) < LSI_max_1_)then
               list_LSI(io) = 1               !LDL
            elseif(LSI(io) > LSI_min_2_ .and. LSI(io) < LSI_max_2_)then 
               list_LSI(io) = 2               !HDL
            else
               list_LSI(io) = 3               !other
            endif   
!Compute LSI distribution and average value
            nsh = LSI(io)/LSI_del+0.5D0            
!check         if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check         write(*,*)LSI(io),nsh
            if(nsh==0)then
               write(*,*)'Warning!!! LSI bin = 0; transformed to 1'
               write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
               nsh = 1
            endif   
!            nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!            write(*,*)nsh,nsh_
!            nsh = nsh_
            NLSI(nsh) = NLSI(nsh) + 1
            AVLSI = AVLSI + LSI(io)
            NLSI_count = NLSI_count + 1                         !count number of values used for distribution and average value calculation
!         end do

!calculate O-O rdf for LDL and HDL populations
            if(LRDF)then
               do jhb = 1,natms,nwatsites        !Loop over water oxygens
                  IF(jhb.NE.ihb)THEN
!check            WRITE(*,*)ihb,jhb
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM ihb AND OXYGEN ATOMS jhb     
                     dx = x(ihb) - x(jhb)             
                     dy = y(ihb) - y(jhb)
                     dz = z(ihb) - z(jhb)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dr2 = dx**2 + dy**2 + dz**2
                     dr  = dsqrt(dr2)
                     R_rdf_OO = dr
!                  ENDIF
!               end do       
!check               write(*,*)R_rdf_OO,cubeh_
                     if(R_rdf_OO < cubeh_)then
!bulk                     
                        NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                        NRDFOO(NSHLOO) = NRDFOO(NSHLOO)+1
                     endif   
                     if(R_rdf_OO < cubeh_.and.list_LSI(io) == 1)then
!LHDL                
                        NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                        NRDFOO_LDL(NSHLOO) = NRDFOO_LDL(NSHLOO)+1
                     endif   
                     if(R_rdf_OO < cubeh_.and.list_LSI(io) == 2)then
!HDL                 
                        NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                        NRDFOO_HDL(NSHLOO) = NRDFOO_HDL(NSHLOO)+1
                     endif               
                  ENDIF            
               end do !End jhb water oxygens loop
            endif                                     !END RDF    
!End rdf          
          
         end do    !End ihb water oxygens loop
             
         list_s = list_LSI
      else
         list_s = 1                          !water - single environment
      endif                                  !End if NRSHMAX > 1  
!END LSI CALCULATION ======================================================================================================== 
   endif                    !End sampling ---> The above block stored each water molecule LSI environment for the last sampling origin

   

!HB TCF CALCULATION 
   io = 0                 !water molecule number 1,2,3,4   ihb = 1,4,7,10
   jo = 0
!sample all possible pairs - WARNING!!! this loop does not distiguish between ik and ki - check for donor and acceptor HBs    
   do i=1,nwatoms1,nwatsites                          !nwatoms1 = (nbwmol-1)*3
      io = io + 1
      do k=i+nwatsites,nwatoms,nwatsites              !Loop over water oxygens (proton acceptors)  
         jo = io + 1
         kpw = kpw + 1                                !kpw must run over the maximum number of pairs of particles
                                                      !WARNING!!! a water molecule (not OH) is considered a single particle
                                                      
!Check environment at t0 - molecules may switch environment throught the tcf calculation for each oeigin        
         if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then 
            nienv = list_s(io)                           !LSI environemnt of water i [not pairs]
            nkenv = list_s(jo)                           !LSI environemnt of water k [not pairs]
!check            write(*,*)'origin',ko,kpw,nienv,nkenv
            if(nienv/=nkenv)nienv = 3
            k_lsi(ko,nienv)  = k_lsi(ko,nienv) + 1     !number of pairs in a given environment at origin ko
            wat_id_nradsh(ko,kpw)= nienv
            
         endif
         nradsh = wat_id_nradsh(ko,kpw)
!check         write(*,*)'step',j,'origin',ko,'env',nradsh,'pairs',k_lsi(ko,nradsh)
! intermolecular OO distance 
         dx = x(i) - x(k)
         dy = y(i) - y(k)
         dz = z(i) - z(k)
         dx = dx - anint(dx/cell(1))*cell(1)
         dy = dy - anint(dy/cell(2))*cell(2)
         dz = dz - anint(dz/cell(3))*cell(3)
         dr2 = dx**2.0 + dy**2.0 + dz**2.0
         dr  = sqrt(dr2)
         sqdrOO = dr2
         drOO = dr
         if(drOO.gt.rHBOO)then
            if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then             !New origin
               ht0_i(ko,kpw,nradsh) = 0.0                                              !at origin ko pair kpw in environment nradsh is not HBonded
               ht0_c(ko,kpw,nradsh) = 0.0
               ht0_s(ko,kpw,nradsh) = 0.0
               sum_h0_int(ko,nradsh)  = sum_h0_int(ko,nradsh)  + ht0_i(ko,kpw,nradsh)
               sum_h0_cont(ko,nradsh) = sum_h0_cont(ko,nradsh) + ht0_c(ko,kpw,nradsh)
               sum_h0_surv(ko,nradsh) = sum_h0_surv(ko,nradsh) + ht0_s(ko,kpw,nradsh)
            endif
            do kt = 1,ko                                                                !New instantaneous value for each origin
!               nradsh_= wat_id_nradsh(kt,kpw)
               ht_i(kt,kpw) = 0.0
               ht_c(kt,kpw) = 0.0
               ht_s(kt,kpw) = 0.0
               lhist_HB(kt,kpw) = 0                                                     !Warning!!! initialized to 1
               lhist_s_HB(kt,kpw) = 0
            end do   
!            sum_h_int(nradsh)= sum_h_int (nradsh)+ 0.0
!            sum_h_cont(nradsh)= sum_h_cont(nradsh)+ 0.0              
                                                                                 !end if/then/else - jump to next water
         elseif(drOO.le.rHBOO)then     
            do iH = 1,2                   !check two hydrogen atoms per water
! if the first H is HBonded leave the loop over the H atoms                 
! intermolecular H...O distance  (donor)                 
               dx = x(i+iH) - x(k)
               dy = y(i+iH) - y(k)
               dz = z(i+iH) - z(k)
               dx = dx - anint(dx/cell(1))*cell(1)
               dy = dy - anint(dy/cell(2))*cell(2)
               dz = dz - anint(dz/cell(3))*cell(3)
               dr2 = dx**2.0 + dy**2.0 + dz**2.0
               dr  = sqrt(dr2)
               sqdrHO = dr2
               drHO = dr
! intermolecular O...H distance  (acceptor)                 
               dx = x(k+iH) - x(i)
               dy = y(k+iH) - y(i)
               dz = z(k+iH) - z(i)
               dx = dx - anint(dx/cell(1))*cell(1)
               dy = dy - anint(dy/cell(2))*cell(2)
               dz = dz - anint(dz/cell(3))*cell(3)
               dr2 = dx**2.0 + dy**2.0 + dz**2.0
               dr  = sqrt(dr2)
               sqarHO = dr2
               arHO = dr   
!check         if(arHO < rHBHO .and. drHO < rHBHO) write(*,*)'Warning!!! water pair "is" r(O...H) mutual donor and acceptor'  
               if((drHO.gt.rHBHO).and.(arHO.gt.rHBHO))then
!               if(drHO.gt.rHBHO)then
                  if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then     
                     if(iH==2)ht0_i(ko,kpw,nradsh) = 0.0                                    !at origin ko pair kpw is not HBonded
                     if(iH==2)ht0_c(ko,kpw,nradsh) = 0.0
                     if(iH==2)ht0_s(ko,kpw,nradsh) = 0.0
                     if(iH==2)sum_h0_int(ko,nradsh)  = sum_h0_int(ko,nradsh) + ht0_i(ko,kpw,nradsh)
                     if(iH==2)sum_h0_cont(ko,nradsh) = sum_h0_cont(ko,nradsh) + ht0_c(ko,kpw,nradsh)
                     if(iH==2)sum_h0_surv(ko,nradsh) = sum_h0_surv(ko,nradsh) + ht0_s(ko,kpw,nradsh)
                  endif
                  if(iH==2)then              !if the first H is not HB try the second - only if neither is HB assign 0 to the water pair
                     do kt = 1,ko                             !New instantaneous value for each origin
!                        nradsh_= wat_id_nradsh(kt,kpw)
                        ht_i(kt,kpw) = 0.0
                        ht_c(kt,kpw) = 0.0
                        ht_s(kt,kpw) = 0.0
                        lhist_HB(kt,kpw) = 0
                        lhist_s_HB(kt,kpw) = 0
                     end do
                  endif   
!                  if(iH==2)sum_h_int(nradsh) = sum_h_int(nradsh)  + 0.0
!                  if(iH==2)sum_h_cont(nradsh) = sum_h_cont(nradsh) + 0.0     
! at least either a donor or acceptor HB exists - calculate intramolecular OH and then HB angle                  
               elseif((drHO.le.rHBHO).or.(arHO.le.rHBHO))then
! intramolecular OH distance (donor) 
                  dx = x(i+iH) - x(i)
                  dy = y(i+iH) - y(i)
                  dz = z(i+iH) - z(i)
                  dx = dx - anint(dx/cell(1))*cell(1)
                  dy = dy - anint(dy/cell(2))*cell(2)
                  dz = dz - anint(dz/cell(3))*cell(3)
                  dr2 = dx**2.0 + dy**2.0 + dz**2.0
                  dr  = sqrt(dr2)
                  sqdriOH = dr2
                  driOH = dr 
! intramolecular OH distance (acceptor) 
                  dx = x(k+iH) - x(k)
                  dy = y(k+iH) - y(k)
                  dz = z(k+iH) - z(k)
                  dx = dx - anint(dx/cell(1))*cell(1)
                  dy = dy - anint(dy/cell(2))*cell(2)
                  dz = dz - anint(dz/cell(3))*cell(3)
                  dr2 = dx**2.0 + dy**2.0 + dz**2.0
                  dr  = sqrt(dr2)
                  sqariOH = dr2
                  ariOH = dr 
!
!Apply the COSINE rule to calculate the HB donor angle
!                 DCOSNM = SQDHO-SQDOH-SQDOO
                  dcosnm = sqdrHO-sqdriOH-sqdrOO                   !donnor(i)
!                 ACOSNM = SQAHO-SQAOH-SQDOO
                  acosnm = sqarHO-sqariOH-sqdrOO                   !acceptor(i)
!donor                  
!                 DCOSDN = -2.0D0*DOHYD*DOXYG
                  dcosdn = -2.0*driOH*drOO
!                 DCOSANG = DCOSNM/DCOSDN
                  dcosang = dcosnm/dcosdn
!                 ANGDOHH = DACOSD(DCOSANG)
!acceptor                
!                 ACOSDN = -2.0D0*AOHYD*DOXYG
                  acosdn = -2.0*ariOH*drOO
!                 ACOSANG = ACOSNM/ACOSDN
                  acosang = acosnm/acosdn
!                 ANGAOHH = DACOSD(ACOSANG)

!acosd - returns arcosin in degrees                - donor
                  if(dcosang<-1.0)then
                     write(*,*)'Warning HB angle',dcosang
                     dcosang = -1.0
                     write(*,*)'HB angle = ',acosd(dcosang)
                     write(*,*)
                  endif
                  if(dcosang>+1.0)then
                     write(*,*)'Warning HB angle',dcosang
                     dcosang = +1.0
                     write(*,*)'HB angle = ',acosd(dcosang)
                     write(*,*)
                  endif
                  angdohh = acosd(dcosang)
!acosd - returns arcosin in degrees                  - acceptor
                  if(acosang<-1.0)then
                     write(*,*)'Warning HB angle',acosang
                     acosang = -1.0
                     write(*,*)'HB angle = ',acosd(acosang)
                     write(*,*)
                  endif
                  if(acosang>+1.0)then
                     write(*,*)'Warning HB angle',acosang
                     acosang = +1.0
                     write(*,*)'HB angle = ',acosd(acosang)
                     write(*,*)
                  endif
                  angaohh = acosd(acosang)                   
!The water pair is not HBonded                                    
                  if((angdohh.ge.aHBOHO).and.(angaohh.ge.aHBOHO))then                 
                     if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then     
                        if(iH==2)ht0_i(ko,kpw,nradsh) = 0.0
                        if(iH==2)ht0_c(ko,kpw,nradsh) = 0.0
                        if(iH==2)ht0_s(ko,kpw,nradsh) = 0.0
                        if(iH==2)sum_h0_int(ko,nradsh)  = sum_h0_int(ko,nradsh)  + ht0_i(ko,kpw,nradsh)
                        if(iH==2)sum_h0_cont(ko,nradsh) = sum_h0_cont(ko,nradsh) + ht0_c(ko,kpw,nradsh)
                        if(iH==2)sum_h0_surv(ko,nradsh) = sum_h0_surv(ko,nradsh) + ht0_s(ko,kpw,nradsh)
                     endif
                     if(iH==2)then
                        do kt = 1,ko
!                           nradsh_= wat_id_nradsh(kt,kpw)
                           ht_i(kt,kpw) = 0.0
                           ht_c(kt,kpw) = 0.0
                           ht_s(kt,kpw) = 0.0
                           lhist_HB(kt,kpw) = 0
                           lhist_s_HB(kt,kpw) = 0
                        end do
                     endif   
!                     if(iH==2)sum_h_int(nradsh)  = sum_h_int(nradsh)  + 0.0
!                     if(iH==2)sum_h_cont(nradsh) = sum_h_cont(nradsh) + 0.0 
!check                     if(nradsh == 1) write(*,*)j,kt,kpw,nradsh,ht_i(kt,kpw)
!The water pair is HBonded                        
                  elseif((angdohh.lt.aHBOHO).or.(angaohh.lt.aHBOHO))then
                     if((j==1.or.mod(j,hbs_samp)==0).and.(j.le.kmaxstep-NDELS))then     
                        ht0_i(ko,kpw,nradsh) = 1.0
                        ht0_c(ko,kpw,nradsh) = 1.0
!check                        if(lhist_HB(ko,kpw)==0)write(*,*)'Error!!!'
!if there was no HB in the previous time-step set ht0_s = 1 - newly formed HB
!if there was an HB in the previous time-step let ht0_s = 0 - count only HBs from the moment they formed
                        if(lhist_s_HB(ko-1,kpw)==0)then     ! at j == 1 no HBs are possible because lhist_s_HB(0,kpw)==1 - subtract one origin from average
                           ht0_s(ko,kpw,nradsh) = 1.0                          ! HB formed at this time origin 
                        endif   
                        sum_h0_int(ko,nradsh)  = sum_h0_int(ko,nradsh)  + ht0_i(ko,kpw,nradsh)
                        sum_h0_cont(ko,nradsh) = sum_h0_cont(ko,nradsh) + ht0_c(ko,kpw,nradsh)
                        sum_h0_surv(ko,nradsh) = sum_h0_surv(ko,nradsh) + ht0_s(ko,kpw,nradsh)
                     endif                    
                     do kt = 1,ko
!HB intermittent         
!                        nradsh_= wat_id_nradsh(kt,kpw)
                        ht_i(kt,kpw) = 1.0
!check                        if(nradsh_ == 1) write(*,*)j,kt,kpw,nradsh_,ht_i(kt,kpw)
!HB continuous - HBs "present" at time t0 - average HB                           
                        if(lhist_HB(kt,kpw)==1)then                     ! lhist_HB is update every time-step - it must be continuously 1
                           ht_c(kt,kpw) = 1.0
                        else
                           ht_c(kt,kpw) = 0.0
                        endif   
!HB survival - HBs "formed" at t0 - newly generated HB                           
                        if(lhist_s_HB(kt,kpw)==1)then
                           ht_s(kt,kpw) = 1.0
                        else
                           ht_s(kt,kpw) = 0.0
                        endif
                     end do   
!                        
!                     sum_h_int(nradsh)  = sum_h_int(nradsh)  + 1.0
!                     if(lhist_HB(ko,kpw)==1)then
!                        sum_h_cont(nradsh) = sum_h_cont(nradsh) + 1.0
!                     else
!                        sum_h_cont(nradsh) = sum_h_cont(nradsh) + 0.0
!                     endif
!End Version 34                         
                     goto 5                                        !leave the loop over the H atoms because the second H will not be HB to the same O-P
                  endif  ! HB angle fulfilled                       
               endif     ! HB H...O distance fulfilled               
            end do       !end loop over H atoms                  
         endif           !HB O...O distanced fulfilled    
         5 continue
      end do                !end loop over water oxygens
   end do                   !end loop over water oxygens

   
!END MAIN HB LOOP   
     
     
!calculate tcfs - for each origin sum over different pairs
   do kt = 1,ko                                   !loop over time origins   
      do kr = 1,NRSHMAX                           !loop over environment
         nt = ntime(kt)                           !ntime(kt) initialized to 1
         if((k_lsi(kt,kr).gt.0).and.(ntime(kt).le.NDELS))then                     !ntime counts the delay times for which the tcf has already been calculated
            do i = 1,npairw                                                       !loop over every pair possible - if a pair in not on a given environment the value is zero
               tcf_hb_int(kt,nt,kr)  = tcf_hb_int(kt,nt,kr)  + ht0_i(kt,i,kr)*ht_i(kt,i)
               tcf_hb_cont(kt,nt,kr) = tcf_hb_cont(kt,nt,kr) + ht0_c(kt,i,kr)*ht_c(kt,i)
               tcf_hb_surv(kt,nt,kr) = tcf_hb_surv(kt,nt,kr) + ht0_s(kt,i,kr)*ht_s(kt,i)
            end do
         endif
         if(k_lsi(kt,kr)==0)k_HB_org(kr) = k_HB_org(kr) + 1
      end do
      ntime(kt) = ntime(kt) + 1                                !delay time counter      
   end do                                      !end loop over time-origins

end do                                         !end loop over origins   



!Print LSI distribution
!if(kLDL_HDL==1)then
   do nkr = 1,LSI_NDELS
      if(NLSI(nkr) > 0)then
         LSI_ang2  = LSI_del*DFLOAT(nkr)
         DISTR_LSI = DFLOAT(NLSI(nkr))/(LSI_del*DFLOAT(NLSI_count))                             !Normalize
         WRITE(n10,39)LSI_ang2,DISTR_LSI
      endif
   end do  
   AVLSI = AVLSI/DFLOAT(NLSI_count)
   write(n0,*)
   write(n0,*)'LSI mean value (Ang) = ',AVLSI
   write(n0,*)

!Print rdf OO
!LOOP OVER RADIAL SHELLS - OXYGEN-OXYGEN RDF
if(LRDF)then
!bulk
   WRITE(n11,'(a5)')'#bulk'
   DO kw=1,NRDELS   
      GROO(kw) = 0.0d0 
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!         GROO(kw)=DFLOAT(NRDFOO(kw))/(DSTNUM*DFLOAT(ko)*FNOM*VOLSHL)
!ok         GROO(kw)=2.0d0*DFLOAT(NRDFOO(kw))/(DSTNUM*dfloat(nOHmean(1))*VOLSHL)
!ok         WRITE(n11,39) RADIUS,GROO(kw)                         
      ENDIF
   END DO
   WRITE(n11,*)
   
!LDL    
   WRITE(n11,'(a4)')'#LDL'
   DO kw=1,NRDELS  
      GROO_LDL(kw) = 0.0d0 
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO_LDL(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!ok         GROO_LDL(kw)=2.0d0*DFLOAT(NRDFOO_LDL(kw))/(DSTNUM*dfloat(nOHmean(2))*VOLSHL)
!ok         WRITE(n11,39) RADIUS,GROO_LDL(kw)                         
      ENDIF
   END DO
   WRITE(n11,*)

!HDL    
   WRITE(n11,'(a4)')'#HDL'
   DO kw=1,NRDELS  
      GROO_HDL(kw) = 0.0d0
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO_HDL(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!ok         GROO_HDL(kw)=2.0d0*DFLOAT(NRDFOO_HDL(kw))/(DSTNUM*dfloat(nOHmean(3))*VOLSHL)
!ok         WRITE(n11,39) RADIUS,GROO_HDL(kw)                         
      ENDIF
   END DO
   WRITE(n11,*)
   
endif
!endif !Print LSI on the fly



!loop over different origins and environments - average over number of particles per origin and number of origins
do kt = 1,ko                                   !loop over time origins  
   do kr = 1,NRSHMAX   
      if(k_lsi(kt,kr) > 0)then
         sum_h0_int_(kr)  = sum_h0_int_(kr)  + sum_h0_int(kt,kr)/float(k_lsi(kt,kr))
         sum_h0_cont_(kr) = sum_h0_cont_(kr) + sum_h0_cont(kt,kr)/float(k_lsi(kt,kr))
         sum_h0_surv_(kr) = sum_h0_surv_(kr) + sum_h0_surv(kt,kr)/float(k_lsi(kt,kr)) 
      endif 
   end do
end do

do kr = 1,NRSHMAX   
   korg_0HB = k_HB_org(kr)
   sum_h0_int_(kr)  = sum_h0_int_(kr)/float(ko-korg_0HB)
   sum_h0_cont_(kr) = sum_h0_cont_(kr)/float(ko-korg_0HB)
   sum_h0_surv_(kr) = sum_h0_surv_(kr)/float(ko-1-korg_0HB)                                   !subtract one origin
end do

!LSI
write(*,*) 'LSI range env 1 = ',LSI_min_1_,'-',LSI_max_1_
write(*,*) 'LSI range env 2 = ',LSI_min_2_,'-',LSI_max_2_
write(n0,*)'LSI range env 1 = ',LSI_min_1_,'-',LSI_max_1_
write(n0,*)'LSI range env 2 = ',LSI_min_2_,'-',LSI_max_2_
!END LSI 


write(*,*)'Number of pairs possible = ',npairw
write(*,*)'Number of pairs counted = ',kpw
write(*,*)'Number of steps = ',kmaxstep
write(*,*)'Number of origins possible = ',norgmax
write(*,*)'Number of origins counted = ',ko
write(*,*)
do kr = 1,NRSHMAX
   write(*,*)'Environment ',kr
   write(*,*)'==============='
   write(*,*)
   do kt = 1,ko 
      write(*,*)'origin nb = ',kt,'Number of pairs = ',k_lsi(kt,kr)
   end do   
!   write(*,*)'intermittent <h> = ',sum_h_int(kr)
!   write(*,*)'continuous <h> = ',sum_h_cont(kr)
   write(*,*)'intermittent <h0> = ',sum_h0_int_(kr)
   write(*,*)'continuous <h0> = ',sum_h0_cont_(kr)
   write(*,*)'survival <h0> = ',sum_h0_surv_(kr)
   write(*,*)
end do
write(n0,*)'Number of pairs possible = ',npairw
write(n0,*)'Number of pairs counted = ',kpw
write(n0,*)'Number of steps = ',kmaxstep
write(n0,*)'Number of origins possible = ',norgmax
write(n0,*)'Number of origins used = ',ko
write(n0,*)
do kr = 1,NRSHMAX
   write(n0,*)'Environment ',kr
   write(n0,*)'==============='
   write(n0,*)
   do kt = 1,ko 
      write(n0,*)'origin nb = ',kt,'Number of pairs = ',k_lsi(kt,kr)
   end do
!   write(n0,*)'intermittent <h> = ',sum_h_int(kr)
!   write(n0,*)'continuous <h> = ',sum_h_cont(kr)
   write(n0,*)'intermittent <h0> = ',sum_h0_int_(kr)
   write(n0,*)'continuous <h0> = ',sum_h0_cont_(kr)
   write(n0,*)'survival <h0> = ',sum_h0_surv_(kr)
   write(n0,*)
end do

!Average over number of pairs
do kr = 1,NRSHMAX     !loop over environments
   do kt = 1,ko       !loop over origins
      if(k_lsi(kt,kr) > 0)then
         do it = 1,NDELS           
            tcf_interm_m(it,kr) = tcf_interm_m(it,kr) + tcf_hb_int(kt,it,kr)/float(k_lsi(kt,kr)) 
            tcf_cont_m(it,kr)   = tcf_cont_m(it,kr)   + tcf_hb_cont(kt,it,kr)/float(k_lsi(kt,kr))  
            tcf_surv_m(it,kr)   = tcf_surv_m(it,kr)   + tcf_hb_surv(kt,it,kr)/float(k_lsi(kt,kr))
         end do
      endif 
   end do
end do

!Divide by number of origins and normalize
do kr = 1,NRSHMAX
   korg_0HB = k_HB_org(kr)
   do it = 1,NDELS  
      tcf_interm_m(it,kr)=tcf_interm_m(it,kr)/(float(ko-korg_0HB)*sum_h0_int_(kr))
      tcf_cont_m(it,kr)=tcf_cont_m(it,kr)/(float(ko-korg_0HB)*sum_h0_cont_(kr))   
      tcf_surv_m(it,kr)=tcf_surv_m(it,kr)/(float(ko-1-korg_0HB)*sum_h0_surv_(kr)) 
   end do 
end do

!calculate and print mean tcfs - average over number of origins     
      write(n1,9)
      write(n2,9)
      write(n4,9)
      do kr = 1,NRSHMAX
         do it = 1,NDELS
            jt = it -1                                              !print tcf starting from time zero
            PSECS = TSTEP*float(jt)*float(KINTVL)
            WRITE(n1,19)PSECS,tcf_cont_m(it,kr)   
            WRITE(n2,19)PSECS,tcf_interm_m(it,kr)  
            write(n4,19)PSECS,tcf_surv_m(it,kr)  
         end do
      WRITE(n1,*)
      WRITE(n2,*)
      WRITE(n4,*)   
      end do
      
!      
!Luzar reactive flux      
!      do it = 2,NDELS1
!         tcf_react(it) = -(tcf_interm_m(it+1)-tcf_interm_m(it-1))/(2.0*TSTEP)
!      end do

!Calculate the reactive flux and the distribution function every 2 time-steps
      do kr = 1,NRSHMAX
         do it = 1,NDELS2,2                 !numerical derivative - intermediate value
            tcf_react(it,kr) = -(tcf_interm_m(it+2,kr)- tcf_interm_m(it,kr))/(2.0*TSTEP)
            Pdist_surv(it,kr) = -(tcf_surv_m(it+2,kr) - tcf_surv_m(it,kr))/(2.0*TSTEP)
            pdist_cont(it,kr) = -(tcf_cont_m(it+2,kr) - tcf_cont_m(it,kr))/(2.0*TSTEP)
         end do
      end do
      
!print reactive flux and P(t)      
      do kr = 1,NRSHMAX
         write(n3,9)
         do it = 1,NDELS2,2
            jt = it -1                                              !print tcf starting from time zero
            PSECS = TSTEP*float(jt)*float(KINTVL)
            WRITE(n3,19)PSECS,tcf_react(it,kr)
            WRITE(n5,19)PSECS,Pdist_surv(it,kr)
            WRITE(n6,19)PSECS,pdist_cont(it,kr)
         end do
         write(n3,*)
         write(n5,*)
         write(n6,*)
      end do
      write(n3,*)
      write(n5,*)
      write(n6,*)
!
!Calculate the reactive flux and the distribution function every time-step
      do kr = 1,NRSHMAX
         do it = 1,NDELS1
            tcf_react(it,kr) = -(tcf_interm_m(it+1,kr)-tcf_interm_m(it,kr))/TSTEP
            Pdist_surv(it,kr) = -(tcf_surv_m(it+1,kr)-tcf_surv_m(it,kr))/TSTEP
            pdist_cont(it,kr) = -(tcf_cont_m(it+1,kr)-tcf_cont_m(it,kr))/TSTEP
         end do
      end do
!print reactive flux and P(t)
      do kr = 1,NRSHMAX
         do it = 1,NDELS1
            jt = it -1                                              !print tcf starting from time zero
            PSECS = TSTEP*float(jt)*float(KINTVL)
            WRITE(n3,19)PSECS,tcf_react(it,kr)
            WRITE(n5,19)PSECS,Pdist_surv(it,kr)
            WRITE(n6,19)PSECS,pdist_cont(it,kr)
         end do
         write(n3,*)
         write(n5,*)
         write(n6,*)
      end do
      write(n3,*)
      write(n5,*)
      write(n6,*)

      
deallocate(tcf_hb_int,tcf_hb_cont,tcf_hb_surv)
deallocate(ht0_i,ht0_c,ht0_s)
deallocate(ht_i,ht_c,ht_s)
deallocate(wat_id_nradsh)

deallocate(ntime)
deallocate(lhist_HB,lhist_s_HB)
deallocate(k_lsi,k_HB_org)

deallocate(tcf_cont_m,tcf_interm_m,tcf_surv_m)
deallocate(tcf_react,Pdist_surv,pdist_cont)
deallocate(sum_h0_int,sum_h0_cont,sum_h0_surv)
deallocate(sum_h0_int_,sum_h0_cont_,sum_h0_surv_)

deallocate(list_s)
deallocate(list_LSI)
deallocate(LSI)
deallocate(NLSI)

deallocate(NRDFOO,NRDFOO_LDL,NRDFOO_HDL)
deallocate(GROO,GROO_LDL,GROO_HDL)     
      
  9   FORMAT('#',11X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 39   FORMAT(27X,F9.4,5X,F14.5)
!    return

END SUBROUTINE pure_wat_HB_LDL_HDL_tcfs

!NG end 


!start development **************************************************************************************


SUBROUTINE pure_wat_LSI_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,lsitframe,lsi_samp,Nenv_lsi,&
                                     LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II)
! Calculate continuous and intermittent Local Structure Index (LSI) tcfs - pure water
! Here each water molecule is a particle 
! The tcfs are calculated for every lsi_samp(fs) time-origin for a delay-time of lsitframe(ps)

!NEW VARIABLES INTRODUCTION
    integer,intent(in)                               :: npurew,nbox,Nenv_lsi
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: lsitframe                       !HB tcfs delay time
    integer,intent(in)                               :: lsi_samp
    real(kind=8),intent(in)                          :: LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II
 
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG    real(kind=4)                               :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                               :: cell(3)  
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                               :: TSTEP,PSECS
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    integer                                    :: NDELS,KINTVL
    
!END NEW VARIABLES INTRODUCTION

!LOCAL VARIABLES    
    integer      :: kmaxstep
    integer      :: i, j, k
    integer      :: n0, n1, n2
    integer      :: n10,n11
    integer      :: it, iH, ko, jt, it1
    integer      :: kpw
    integer      :: kt, nt, jw, kw
    integer      :: io,ihb,jhb,kr

    logical filex
    
! LSI VARIABLES
    integer,dimension(:),allocatable             :: list_LSI                                   !LSI index - LSI computed locally to compare with Gaia List
    integer                                      :: pwnmolsol,pwnatmsol,pwnions                !solution parameters - set to 0
    real(kind=8),dimension(:),allocatable        :: LSI
    integer,dimension(:),allocatable             :: NLSI
    real(kind=8)                                 :: r_lsi
    real(kind=8)                                 :: LSI_del,AVLSI,LSI_ang2,DISTR_LSI
    integer                                      :: nsh,LSI_NDELS,NLSI_count,nkr
!    integer                                      :: nsh_
    real(kind=8)                                 :: R_rdf_OO,RDEL,cubeh_,cube_,RADIUS,PI,VOLSHL,DSTNUM,FNOM
    integer                                      :: NSHLOO,NRDELS
    integer                                      :: natcheck
    integer,dimension(:),allocatable             :: NRDFOO,NRDFOO_LDL,NRDFOO_HDL
    real(kind=8),dimension(:),allocatable        :: GROO,GROO_LDL,GROO_HDL
    integer                                      :: nradsh,nradsh_
    logical                                      :: LRDF
    integer,dimension(:,:),allocatable           :: k_lsi
    integer,dimension(:),allocatable             :: k_HB_org
    integer                                      :: korg_0HB
!END LSI VARIABLES

! LSI tcfs    
    real(kind=8),dimension(:,:),allocatable    :: tcf_cont_m,tcf_interm_m    
    
!    real(kind=8),dimension(:),allocatable      :: sum_h_int,sum_h_cont
    real(kind=8),dimension(:,:),allocatable    :: sum_h0_int,sum_h0_cont
    real(kind=8),dimension(:),allocatable      :: sum_h0_int_,sum_h0_cont_
!
    integer                                    :: norgmax
    integer                                    :: NRSHMAX
    real(kind=8),dimension(:,:,:),allocatable  :: tcf_hb_int,tcf_hb_cont
    real(kind=8),dimension(:,:,:),allocatable  :: ht0_i,ht0_c
    real(kind=8),dimension(:,:,:),allocatable  :: ht_i,ht_c
    integer,dimension(:),allocatable           :: ntime
    integer,dimension(:,:,:),allocatable       :: lhist_HB
    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!NEW CODE LSI PARAMETERS AND XTC

!LSI values
r_lsi  = 3.70d0                                 !original value of the LSI
LSI_del = 0.0005d0                              !bin for LSI distribution
LSI_NDELS= 1.0d0/LSI_del

RDEL    = 0.05d0                                !bin size in rdf/distance distributions Ang
NRDELS  = 100.0d0/RDEL
LRDF    =.false.
PI      = 3.14159265d0
FNOM = dfloat(nmolwat)

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj

!END NEW CODE LSI PARAMETERS AND XTC

kmaxstep = nstep

! time-window for calculation of each tcf
NDELS = lsitframe*1000                              !sets the delay time-window in fs
norgmax = kmaxstep/lsi_samp + 1                    !maximum number of origins
if(norgmax.gt.kmaxstep)norgmax = kmaxstep
                 
if(NDELS.ge.kmaxstep)NDELS = NDELS - 1
                               
NRSHMAX = Nenv_lsi                                   !maximum number of environments = 3 [bulk environment not calculated in this routine]
NRSHMAX = 3
!NRSHMAX = 2

!TSTEP - TIME-STEP IN ps
TSTEP  = dt*1.0e-3
!KINTVL - INTERVAL BETWEEN TIME-STEPS
KINTVL = 1

allocate(tcf_hb_int(norgmax,NDELS,NRSHMAX),tcf_hb_cont(norgmax,NDELS,NRSHMAX))
allocate(ht0_i(norgmax,nmolwat,NRSHMAX),ht0_c(norgmax,nmolwat,NRSHMAX))
allocate(ht_i(norgmax,nmolwat,NRSHMAX),ht_c(norgmax,nmolwat,NRSHMAX))

allocate(ntime(norgmax))
allocate(lhist_HB(norgmax,nmolwat,NRSHMAX))
allocate(k_lsi(norgmax,NRSHMAX),k_HB_org(norgmax))

tcf_hb_int = 0.0; tcf_hb_cont = 0.0
ht0_i = 0.0; ht0_c = 0.0; ht_i = 0.0; ht_c = 0.0
ntime = 1
lhist_HB = 1
k_lsi = 0; k_HB_org = 0

!HB tcfs

allocate(tcf_cont_m(NDELS,NRSHMAX),tcf_interm_m(NDELS,NRSHMAX))
allocate(sum_h0_int(norgmax,NRSHMAX),sum_h0_cont(norgmax,NRSHMAX))
allocate(sum_h0_int_(NRSHMAX),sum_h0_cont_(NRSHMAX))
!allocate(sum_h_int(NRSHMAX),sum_h_cont(NRSHMAX))

tcf_cont_m = 0.0;tcf_interm_m = 0.0
sum_h0_int = 0.0;sum_h0_cont = 0.0
sum_h0_int_ = 0.0;sum_h0_cont_ = 0.0

!List HDL/LDL
allocate(list_LSI(nmolwat))
allocate(LSI(nmolwat))
allocate(NLSI(LSI_NDELS))

allocate(NRDFOO(NRDELS),NRDFOO_LDL(NRDELS),NRDFOO_HDL(NRDELS))
allocate(GROO(NRDELS),GROO_LDL(NRDELS),GROO_HDL(NRDELS))

LSI = 0
NLSI = 0
AVLSI = 0.0d0
pwnmolsol = 0
pwnatmsol = 0
pwnions   = 0  
NLSI_count = 0

!rdf
NRDFOO     = 0
NRDFOO_LDL = 0
NRDFOO_HDL = 0
!End List HDL/LDL

write(*,*)
write(*,*)'Pure water LSI tcf calculation - LSI environments'
write(*,*)'Routine pure_wat_LSI_tcfs'
write(*,*)'Maximum number of origins = ',norgmax
write(*,*)'water LSI tcf delay time-window (ps) = ',lsitframe
write(*,*)'water LSI sampling frequency (fs) = ',lsi_samp
write(*,*)'results are printed to pureW_LSI_tcf_out'
write(*,*)'WARNING!!! NUMBER OF ENVIRONMENTS = ',NRSHMAX
write(*,*)

inquire(file='pureW_LSI_tcf_out/log_LSI_tcf.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir pureW_LSI_tcf_out')

!output
n0 = 50
n1 = 100
n2 = 200
n10 = 1000
n11 = 1100 

open(n0,file='pureW_LSI_tcf_out/log_LSI_tcf.dat',status='unknown',action='write')
open(n1,file='pureW_LSI_tcf_out/LSI_tcf_cont.dat',status='unknown',action='write')
open(n2,file='pureW_LSI_tcf_out/LSI_tcf_interm.dat',status='unknown',action='write')

open(n10,file='pureW_LSI_tcf_out/LSI_dist.dat',status='unknown',action='write')
open(n11,file='pureW_LSI_tcf_out/rdf_LSI_OO.dat',status='unknown',action='write')

write(n0,*)
write(n0,*)'Pure water LSI tcf calculation - LSI environments'
write(n0,*)'Routine pure_wat_LSI_tcfs'
write(n0,*)'water LSI tcf delay time-window (ps) = ',lsitframe
write(n0,*)'water LSI sampling frequency (fs) = ',lsi_samp
write(n0,*)'results are printed to pureW_LSI_tcf_out'
write(n0,*)'WARNING!!! NUMBER OF ENVIRONMENTS = ',NRSHMAX
write(n0,*)

write(*,*) 
write(*,*)'Starting <h_LSI(0)h_LSI(t)> calculation'
write(*,*) 
write(n0,*) 
write(n0,*)'Starting <h_LSI(0)h_LSI(t)> calculation'
write(n0,*) 

!NEW READ xyz

ko = 0                                                               !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep                                                
   if(mod(j,500)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms in solv_inp is wrong'
         stop
      endif      
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!End Version 9      
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   
    
    cubeh_ = cell(1)/2.0d0
    cube_  = cell(1)

   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 
!END NEW READ xyz


!Calculate LSI 
   io = 0                                                                           !water molecule number 1,2,3,4   ihb = 1,4,7,10
   ht_i = 0                                                                         !re-initialize ht_i
   ht_c = 0
   if((j==1.or.mod(j,lsi_samp)==0).and.(j.le.kmaxstep-NDELS))ko = ko + 1
!   
   do ihb = 1,natms,nwatsites
      LSI = 0
      io = io + 1 
!      write(*,*)'Determining LSI - step = ',j,'water = ',ihb 
!NG         call LSI_water(ihb,io,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,LSI,r_lsi)
      call LSI_water_fast(ihb,io,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,LSI,r_lsi)
!      write(*,*)'Ended determining LSI - step = ',j,'water = ',ihb  
      if(LSI(io) > LSI_min_I .and. LSI(io) < LSI_max_I)then
         list_LSI(io) = 1               !LDL
         nradsh = 1
      elseif(LSI(io) > LSI_min_II .and. LSI(io) < LSI_max_II)then 
         list_LSI(io) = 2               !HDL
         nradsh = 2
      else
         list_LSI(io) = 3               !other
         nradsh = 3
      endif
      
      if((j==1.or.mod(j,lsi_samp)==0).and.(j.le.kmaxstep-NDELS))then               !New origin
!         ko = ko + 1                                                            !origins counter
!check         write(*,*)j,ko
         ht0_i(ko,io,nradsh) = 1.0                                              !at origin ko water io is in environment nradsh
         ht0_c(ko,io,nradsh) = 1.0
         sum_h0_int(ko,nradsh)  = sum_h0_int(ko,nradsh)  + ht0_i(ko,io,nradsh)
         sum_h0_cont(ko,nradsh) = sum_h0_cont(ko,nradsh) + ht0_c(ko,io,nradsh)
         k_lsi(ko,nradsh)  = k_lsi(ko,nradsh) + 1                               !number of waters in a given environment at origin k
      endif
!      
      do kt = 1,ko                                                                !New instantaneous value for each origin
         ht_i(kt,io,nradsh) = 1.0                   !Re-initialized to 0 every time-step

!Update list for continuous tcf    
!!Warning!!! lhist_HB initialized to 1
         if(nradsh == 1)then
            lhist_HB(kt,io,2) = 0                    
            lhist_HB(kt,io,3) = 0
         elseif(nradsh == 2)then   
            lhist_HB(kt,io,1) = 0                    
            lhist_HB(kt,io,3) = 0
         elseif(nradsh == 3)then   
            lhist_HB(kt,io,1) = 0                    
            lhist_HB(kt,io,2) = 0
         else
            write(*,*)'Error - LSI'
            stop
         endif         
         
!tcf continuous    
!         do kr=1,NRSHMAX
         if(lhist_HB(kt,io,nradsh)==1)then                     !lhist_HB is updated every time-step - it must be continuously 1
            ht_c(kt,io,nradsh) = 1.0                           !Re-initialized to 0 every time-step
         else
            ht_c(kt,io,nradsh) = 0.0
         endif               
         
!check         if(ht0_c(1,io,3) == 1.0)write(*,*)j,ht_c(1,io,3),ht_i(1,io,3)
         
         if((ht_c(kt,io,nradsh) == 1.0).and.(ht_i(kt,io,nradsh) == 0.0))then
            write(*,*)'Error LSI tcf-continuous'
            stop
!check         elseif((ht_c(kt,io,nradsh) == 0.0).and.(ht_i(kt,io,nradsh) == 1.0))then
!check         write(*,*)'LSI environment re-crossing','step = ',j,'Environment = ',nradsh
         endif   
!         end do
      end do  
      
!Compute LSI distribution and average value
      nsh = LSI(io)/LSI_del+0.5D0            
!check         if(nsh==0)write(*,*)'ERROR!!! LSI =',LSI(io),'LSI_del =',LSI_del
!check         write(*,*)LSI(io),nsh
      if(nsh==0)then
         write(*,*) 'Warning!!! LSI bin = 0; transformed to 1'
         write(n0,*)'Warning!!! LSI bin = 0; transformed to 1'
         nsh = 1
      endif   
!      nsh_ = IDNINT(LSI(io)/LSI_del)               !ok this is the same as: nsh = LSI(io)/LSI_del+0.5D0
!      write(*,*)nsh,nsh_
!      nsh = nsh_
      NLSI(nsh) = NLSI(nsh) + 1
      AVLSI = AVLSI + LSI(io)
      NLSI_count = NLSI_count + 1                         !count number of values used for distribution and average value calculation
!      end do

!calculate O-O rdf for LDL and HDL populations
      if(LRDF)then
         do jhb = 1,natms,nwatsites        !Loop over water oxygens
            IF(jhb.NE.ihb)THEN
!check            WRITE(*,*)ihb,jhb
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM ihb AND OXYGEN ATOMS jhb     
               dx = x(ihb) - x(jhb)             
               dy = y(ihb) - y(jhb)
               dz = z(ihb) - z(jhb)
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2
               dr  = dsqrt(dr2)
               R_rdf_OO = dr
!            ENDIF
!         end do       
!check         write(*,*)R_rdf_OO,cubeh_
               if(R_rdf_OO < cubeh_)then
!bulk               
                  NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                  NRDFOO(NSHLOO) = NRDFOO(NSHLOO)+1
               endif   
               if(R_rdf_OO < cubeh_.and.list_LSI(io) == 1)then
!LHDL          
                  NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                  NRDFOO_LDL(NSHLOO) = NRDFOO_LDL(NSHLOO)+1
               endif   
               if(R_rdf_OO < cubeh_.and.list_LSI(io) == 2)then
!HDL           
                  NSHLOO = R_rdf_OO/RDEL+0.5D0                                        
                  NRDFOO_HDL(NSHLOO) = NRDFOO_HDL(NSHLOO)+1
               endif               
            ENDIF            
         end do !End jhb water oxygens loop
      endif                                     !END RDF    
!End rdf              
   end do    !End ihb water oxygens loop
          
!END LSI CALCULATION 
     
!calculate tcfs - for each origin sum over different pairs
!   write(*,*)'Accumulate tcf - step = ',j
   do kt = 1,ko                                   !loop over time origins   
      do kr = 1,NRSHMAX                           !loop over environment
         nt = ntime(kt)                           !ntime(kt) initialized to 1
         if((k_lsi(kt,kr).gt.0).and.(ntime(kt).le.NDELS))then                     !ntime counts the delay times for which the tcf has already been calculated
            do i = 1,nmolwat                                                       !loop over every pair possible - if a pair in not on a given environment the value is zero
               tcf_hb_int(kt,nt,kr)  = tcf_hb_int(kt,nt,kr)  + ht0_i(kt,i,kr)*ht_i(kt,i,kr)
               tcf_hb_cont(kt,nt,kr) = tcf_hb_cont(kt,nt,kr) + ht0_c(kt,i,kr)*ht_c(kt,i,kr)
            end do
         endif
         if(k_lsi(kt,kr)==0)k_HB_org(kr) = k_HB_org(kr) + 1
      end do
      ntime(kt) = ntime(kt) + 1                                !delay time counter      
   end do                                      !end loop over time-origins
!   write(*,*)'End accumulating tcf - step = ',j

end do                                         !end loop over origins   

!Print LSI distribution
!if(kLDL_HDL==1)then
   do nkr = 1,LSI_NDELS
      if(NLSI(nkr) > 0)then
         LSI_ang2  = LSI_del*DFLOAT(nkr)
         DISTR_LSI = DFLOAT(NLSI(nkr))/(LSI_del*DFLOAT(NLSI_count))                             !Normalize
         WRITE(n10,39)LSI_ang2,DISTR_LSI
      endif
   end do  
   AVLSI = AVLSI/DFLOAT(NLSI_count)
   write(n0,*)
   write(n0,*)'LSI mean value (Ang) = ',AVLSI
   write(n0,*)

!Print rdf OO
!LOOP OVER RADIAL SHELLS - OXYGEN-OXYGEN RDF
if(LRDF)then
!bulk
   WRITE(n11,'(a5)')'#bulk'
   DO kw=1,NRDELS   
      GROO(kw) = 0.0d0 
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!         GROO(kw)=DFLOAT(NRDFOO(kw))/(DSTNUM*DFLOAT(ko)*FNOM*VOLSHL)
!ok         GROO(kw)=2.0d0*DFLOAT(NRDFOO(kw))/(DSTNUM*dfloat(nOHmean(1))*VOLSHL)
!ok         WRITE(n11,39) RADIUS,GROO(kw)                         
      ENDIF
   END DO
   WRITE(n11,*)
   
!LDL    
   WRITE(n11,'(a4)')'#LDL'
   DO kw=1,NRDELS  
      GROO_LDL(kw) = 0.0d0 
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO_LDL(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!ok         GROO_LDL(kw)=2.0d0*DFLOAT(NRDFOO_LDL(kw))/(DSTNUM*dfloat(nOHmean(2))*VOLSHL)
!ok         WRITE(n11,39) RADIUS,GROO_LDL(kw)                         
      ENDIF
   END DO
   WRITE(n11,*)

!HDL    
   WRITE(n11,'(a4)')'#HDL'
   DO kw=1,NRDELS  
      GROO_HDL(kw) = 0.0d0
      DSTNUM = FNOM/(cube_**3.0d0)
      IF (NRDFOO_HDL(kw) > 0) THEN                                   
         RADIUS=RDEL*DFLOAT(kw)                               
         VOLSHL=4.0D0*PI*RDEL*RADIUS*RADIUS+(PI*RDEL**3.0)/3.0D0         
!ok         GROO_HDL(kw)=2.0d0*DFLOAT(NRDFOO_HDL(kw))/(DSTNUM*dfloat(nOHmean(3))*VOLSHL)
!ok         WRITE(n11,39) RADIUS,GROO_HDL(kw)                         
      ENDIF
   END DO
   WRITE(n11,*)
   
endif
!endif !Print LSI on the fly



!loop over different origins and environments - average over number of particles per origin and number of origins
do kt = 1,ko                                   !loop over time origins  
   do kr = 1,NRSHMAX   
      if(k_lsi(kt,kr) > 0)then
         sum_h0_int_(kr)  = sum_h0_int_(kr)  + sum_h0_int(kt,kr)/float(k_lsi(kt,kr))
         sum_h0_cont_(kr) = sum_h0_cont_(kr) + sum_h0_cont(kt,kr)/float(k_lsi(kt,kr))
      endif 
   end do
end do

do kr = 1,NRSHMAX   
   korg_0HB = k_HB_org(kr)
   sum_h0_int_(kr)  = sum_h0_int_(kr)/float(ko-korg_0HB)
   sum_h0_cont_(kr) = sum_h0_cont_(kr)/float(ko-korg_0HB)
end do

!LSI
write(*,*) 'LSI range env 1 = ',LSI_min_I,'-',LSI_max_I
write(*,*) 'LSI range env 2 = ',LSI_min_II,'-',LSI_max_II
write(n0,*)'LSI range env 1 = ',LSI_min_I,'-',LSI_max_I
write(n0,*)'LSI range env 2 = ',LSI_min_II,'-',LSI_max_II
!END LSI 


write(*,*)'Number of steps = ',kmaxstep
write(*,*)'Number of origins possible = ',norgmax
write(*,*)'Number of origins counted = ',ko
write(*,*)

do kr = 1,NRSHMAX
   write(*,*)'Environment ',kr
   write(*,*)'==============='
   write(*,*)
   do kt = 1,ko 
      write(*,*)'origin nb = ',kt,'Number of waters = ',k_lsi(kt,kr)
   end do   
!   write(*,*)'intermittent <h> = ',sum_h_int(kr)
!   write(*,*)'continuous <h> = ',sum_h_cont(kr)
   write(*,*)'intermittent <h0> = ',sum_h0_int_(kr)
   write(*,*)'continuous <h0> = ',sum_h0_cont_(kr)
   write(*,*)
end do

write(n0,*)'Number of steps = ',kmaxstep
write(n0,*)'Number of origins possible = ',norgmax
write(n0,*)'Number of origins used = ',ko
write(n0,*)
do kr = 1,NRSHMAX
   write(n0,*)'Environment ',kr
   write(n0,*)'==============='
   write(n0,*)
   do kt = 1,ko 
      write(n0,*)'origin nb = ',kt,'Number of waters = ',k_lsi(kt,kr)
   end do
!   write(n0,*)'intermittent <h> = ',sum_h_int(kr)
!   write(n0,*)'continuous <h> = ',sum_h_cont(kr)
   write(n0,*)'intermittent <h0> = ',sum_h0_int_(kr)
   write(n0,*)'continuous <h0> = ',sum_h0_cont_(kr)
   write(n0,*)
end do

!Average over number of pairs
do kr = 1,NRSHMAX     !loop over environments
   do kt = 1,ko       !loop over origins
      if(k_lsi(kt,kr) > 0)then
         do it = 1,NDELS           
            tcf_interm_m(it,kr) = tcf_interm_m(it,kr) + tcf_hb_int(kt,it,kr)/float(k_lsi(kt,kr)) 
            tcf_cont_m(it,kr)   = tcf_cont_m(it,kr)   + tcf_hb_cont(kt,it,kr)/float(k_lsi(kt,kr))  
         end do
      endif 
   end do
end do

!Divide by number of origins and normalize
do kr = 1,NRSHMAX
   korg_0HB = k_HB_org(kr)
   do it = 1,NDELS  
      tcf_interm_m(it,kr)=tcf_interm_m(it,kr)/(float(ko-korg_0HB)*sum_h0_int_(kr))
      tcf_cont_m(it,kr)=tcf_cont_m(it,kr)/(float(ko-korg_0HB)*sum_h0_cont_(kr))   
   end do 
end do

!calculate and print mean tcfs - average over number of origins     
      write(n1,9)
      write(n2,9)
      do kr = 1,NRSHMAX
         do it = 1,NDELS
            jt = it -1                                              !print tcf starting from time zero
            PSECS = TSTEP*float(jt)*float(KINTVL)
            WRITE(n1,19)PSECS,tcf_cont_m(it,kr)   
            WRITE(n2,19)PSECS,tcf_interm_m(it,kr)  
         end do
      WRITE(n1,*)
      WRITE(n2,*)   
      end do            

      
deallocate(tcf_hb_int,tcf_hb_cont)
deallocate(ht0_i,ht0_c)
deallocate(ht_i,ht_c)
deallocate(ntime)
deallocate(lhist_HB)
deallocate(k_lsi,k_HB_org)
deallocate(tcf_cont_m,tcf_interm_m)
deallocate(sum_h0_int,sum_h0_cont)
deallocate(sum_h0_int_,sum_h0_cont_)
deallocate(list_LSI)
deallocate(LSI)
deallocate(NLSI)
deallocate(NRDFOO,NRDFOO_LDL,NRDFOO_HDL)
deallocate(GROO,GROO_LDL,GROO_HDL)      
      
      
  9   FORMAT('#',11X,'time(ps)',9X,'tcf'/)  
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 39   FORMAT(27X,F9.4,5X,F14.5)
!    return

END SUBROUTINE pure_wat_LSI_tcfs



SUBROUTINE pure_wat_tetra(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,tetsampw,QMMM_pureW) 
           
! Calculate water's tetrahedrality in different environments - single enviroment for pure water

    integer,intent(in)                               :: npurew,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    integer,intent(in)                               :: tetsampw
    integer,intent(in)                               :: QMMM_pureW
    character(len=4),dimension(nwatsites)                :: atom_name_xyz 
    character(len=4),dimension(natms)                :: atomname

! Local variables
! Tetrahedrality and distance distributions
    real(kind=8)                                     :: SG,CHHASG,RO4
    real(kind=8)                                     :: RO_nb(5)
    real(kind=8)                                     :: RH2,RHO1
    integer                                          :: nsh
    real(kind=8)                                     :: SGDEL,RDEL
    integer                                          :: NWTDNB(4)
    REAL(KIND=4)                                     :: RNTDWH(2)
    integer,dimension(:,:),allocatable               :: NPDSG,NDPRO4,NDPRH2,NDPRHO1
    integer,dimension(:,:,:),allocatable             :: NDPRO
    real(kind=8),dimension(:),allocatable            :: AVWSG,AVCHHASG,AVRO4,AVRH2,AVRHO1
    real(kind=8),dimension(:,:),allocatable          :: AVRO
    integer,parameter                                :: nshdist = 50000                      !max. number of bins in tetrahedrality and distance dist.
    integer                                          :: NTDELS,NRDELS
    real(kind=8)                                     :: TETHQ,RD_ANG
    real(kind=8)                                     :: TETHDIST,DISTR_R 
! 
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=4),dimension(:),allocatable            :: xtr,ytr,ztr
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index
    integer,dimension(:),allocatable                 :: kv
    integer      :: n0, n1, n2, n3,n4, n5, n6, n7, n8, n9, n10, n11
    integer      :: i, j, k
    integer      :: kw, iw, kt, kr, kvi, ki, jw, ik
    integer      :: iH, ihat, nt, ko, jt
    integer      :: kt1, kt2, ni, nc
    integer      :: io
    integer      :: L, KN, IT, NKR
    integer,parameter                                 :: INTMAX = 2147483647
    integer,dimension(:),allocatable                  :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable                  :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                                      :: dr_bulk,dr4nb
    integer                                           :: NRSHMAX   
    integer                                           :: nradsh    
    integer                                           :: natcheck
    integer                                           :: kr_start
    integer                                           :: kstps
    logical                                           :: filex 
    integer                                           :: pwnmolsol,pwnatmsol,pwnions
!Version 9    
    real(kind=8)                                      :: ROMAX,ROMIN
    integer,dimension(:),allocatable                  :: LSTORD 
    integer                                           :: LO
    real(kind=4)                                      :: dx1,dy1,dz1,wdsq1
    real(kind=4),dimension(:),allocatable             :: wd1
    logical                                           :: lgas
!end Version 9
    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!Version 9
atom_name_xyz(1) ='O'
atom_name_xyz(2) ='H'
atom_name_xyz(3) ='H'
if(nwatsites ==4)atom_name_xyz(4) = 'MW'
!End Version 9

lgas=.false.
if(nmolwat==1)then
   lgas=.true.
   write(*,*)'The system is a gas - tetrahedrality calculation will be skipped'
endif   

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
allocate(xtr(natms),ytr(natms),ztr(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Tet_PW_out/log_PW_tet.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Tet_PW_out')   

pwnmolsol = 0
pwnatmsol = 0
pwnions = 0
    
n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140
n5 = 150
n6 = 160
n7 = 170
n8 = 180
n9 = 190
n10 = 200
n11 = 210

open(n0,file='Tet_PW_out/log_env_tetra.dat',status='unknown',action='write')
open(n1,file='Tet_PW_out/wat_tetra_dist.dat',status='unknown',action='write')
open(n2,file='Tet_PW_out/wat_rO4_dist.dat',status='unknown',action='write')
open(n3,file='Tet_PW_out/wat_rOH2_dist.dat',status='unknown',action='write')
open(n4,file='Tet_PW_out/wat_rHO1_dist.dat',status='unknown',action='write')
open(n5,file='Tet_PW_out/wat_rO1_nb_dist.dat',status='unknown',action='write')
open(n6,file='Tet_PW_out/wat_rO2_nb_dist.dat',status='unknown',action='write')
open(n7,file='Tet_PW_out/wat_rO3_nb_dist.dat',status='unknown',action='write')
open(n8,file='Tet_PW_out/wat_rO4_nb_dist.dat',status='unknown',action='write')
open(n9,file='Tet_PW_out/wat_rO5_nb_dist.dat',status='unknown',action='write')

open(n10,file='Tet_PW_out/QMMM_conf.xyz',status='unknown',action='write')

open(n11,file='Tet_PW_out/QMMM_check.xyz',status='unknown',action='write')
    
kr_start = 1                                     !environments loop started

NRSHMAX = 1                                      !maximum number of environments - single environment for pure water
SGDEL   = 0.01D0                                 !tetrahedrality bin size 
RDEL    = 0.01D0                                 !bin size in distance distributions Ang
NRDELS  = 100.0d0/RDEL                            !maximum number of bins in distance distributions - set to 10 Ang

if(NRDELS>nshdist)then
   write(*,*)'error - routine pure_wat_tetra - too many bins'
   stop
endif   

!Tetrahedrality and distances
allocate(NPDSG(NRSHMAX,nshdist),NDPRO4(NRSHMAX,nshdist),NDPRH2(NRSHMAX,nshdist),NDPRHO1(NRSHMAX,nshdist))
allocate(NDPRO(NRSHMAX,5,nshdist))
allocate(AVWSG(NRSHMAX),AVCHHASG(NRSHMAX),AVRO4(NRSHMAX),AVRH2(NRSHMAX),AVRHO1(NRSHMAX))
allocate(AVRO(NRSHMAX,5))
NPDSG=0;NDPRO4=0;NDPRH2=0;NDPRHO1=0
NDPRO=0
AVWSG=0.0d0;AVCHHASG=0.0d0;AVRO4=0.0d0;AVRH2=0.0d0;AVRHO1=0.0d0
AVRO=0.0d0

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))   
allocate(kv(NRSHMAX))
allocate(LSTORD(nmolwat)) 
allocate(wd1(nmolwat))
NTDW4O = 0; NONTDW = 0; LSTORD = 0; wd1 = 0.0
kv = 0
kstps = 0
nc = 1                     !counter of the QM/MM (xyz) configurations

write(*,*)
write(*,*)'Water tetrahedrality'
write(*,*)'Routine pure_wat_tetra'
write(*,*)'Results are printed to Tet_PW_out'
write(*,*)

! Start tetrahedrality calculation
write(*,*)
write(*,*)'Starting tetrahedrality calculation...'
write(*,*)
write(n0,*)
write(n0,*)'Starting tetrahedrality calculation...'
write(n0,*)

do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
   NTDW4O = 0; NONTDW = 0
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!check      
!      write(n11,*)natms
!      write(n11,*)
!      do i = 1,natms,nwatsites
!         write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),x(i),y(i),z(i)                                !print first QM/MM molecule
!         do jw = 1,nwatsites-1                
!            write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),x(i+jw),y(i+jw),z(i+jw) 
!         end do
!      end do   
      
!End Version 9

   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
   
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
!check       write(*,*)x(kp),y(kp),z(kp)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules
!sample frequency - calculate the tetrahedrality every tetsampw time-steps  

   if(.not.lgas.and.(j==1.or.mod(j,tetsampw)==0))then
      write(*,*)'time-step',j,' sampling water molecules'
      kstps = kstps + 1                                                                            !count number of origins sampled

!Version 9 check
!      write(n11,*)natms
!      write(n11,*)'# comment line',j
!End version 9 check
      
      io = 0
      do i=1,natms,nwatsites        !Loop over water oxygens
         io = io + 1                                         !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water                                                   
            nradsh = 1                                                                              !hydration shell environment
            kv(nradsh) = kv(nradsh) + 1                                                             !number of waters
!NG            call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
!NG                             SG,CHHASG,RO4,RO_nb,RH2,RHO1)
            
            call tetrahedron(i,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,&
                             SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)  
                             
!Version 9
!             write(*,*)'Here we are now',j,i,SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH
!check    
!            write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),x(i),y(i),z(i)                                
!            do jw = 1,nwatsites-1                
!               write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),x(i+jw),y(i+jw),z(i+jw) 
!            end do
!end check            
!write xyz file with each water i at the origin (0,0,0) and the remaining waters ordered by their distance to water i                             
            if(QMMM_pureW == 1.and.mod(j,10*tetsampw)==0.and.nc.le.500)then
               write(*,*)
               write(*,*)'writing 500 configuration for QM/MM analysis - step ',j
               write(*,*)
               write(n10,*)natms
               write(n10,'(a19,1x,f7.3,2(1x,a22,f7.3))')'# Tetrah water 1 = ', SG,' R_inter_O1-Ha (A) = ',RNTDWH(1),' R_inter_O1-Hb (A) = ',RNTDWH(2)       
               do jw=1,natms
!translate water molecules to put water i at the origin
                  xtr(jw) = x(jw) - x(i) 
                  ytr(jw) = y(jw) - y(i) 
                  ztr(jw) = z(jw) - z(i)
               end do   
               
!apply pbc if oxygen is out of the box
               do jw=1,natms,nwatsites
                  if(xtr(jw).lt.-cell(1)/2.0)then
                     do jt=0,nwatsites-1
                        xtr(jw+jt)=xtr(jw+jt)+cell(1)
                     end do
                  endif
                  if(xtr(jw).gt.cell(1)/2.0)then
                     do jt=0,nwatsites-1
                        xtr(jw+jt)=xtr(jw+jt)-cell(1)
                     end do
                  endif
                  
                  if(ytr(jw).lt.-cell(2)/2.0)then
                     do jt=0,nwatsites-1
                        ytr(jw+jt)=ytr(jw+jt)+cell(2)
                     end do
                  endif
                  if(ytr(jw).gt.cell(2)/2.0)then
                     do jt=0,nwatsites-1
                        ytr(jw+jt)=ytr(jw+jt)-cell(2)
                     end do
                  endif
                  
                  if(ztr(jw).lt.-cell(3)/2.0)then
                     do jt=0,nwatsites-1
                        ztr(jw+jt)=ztr(jw+jt)+cell(3)
                     end do
                  endif
                  if(ztr(jw).gt.cell(3)/2.0)then
                     do jt=0,nwatsites-1
                        ztr(jw+jt)=ztr(jw+jt)-cell(3)
                     end do
                  endif                 
               end do

!order water molecules (oxygen atoms) around water (oxygen atom) i - calculate distances 
               LO=1
               do jw=1,natms,nwatsites
                  dx1 = xtr(jw) - xtr(i)
                  dy1 = ytr(jw) - ytr(i) 
                  dz1 = ztr(jw) - ztr(i)
                  wdsq1 = dx1**2.0+dy1**2.0+dz1**2.0
                  wd1(LO) = sqrt(wdsq1)
                  LO = LO + 1
               end do
               
!ORDER WATER MOLECULES - LSTORD(LO) LISTS THE MOLECULES BY THEIR DISTANCE TO water i    
!THE LAST POSITION OF THE ARRAY WILL BE ZERO SINCE WATER i IS THE ORIGIN AND THEREFORE IS NOT ORDERED
               ROMAX = 50.0D0
               ROMIN = 0.0D0
               LO=1
               DO WHILE(LO.LE.nmolwat)
                  ni = 1
                  DO ko=1,natms,nwatsites
                     IF((wd1(ni).GT.ROMIN).AND.(wd1(ni).LT.ROMAX))THEN
                        ROMAX = wd1(ni)
                        LSTORD(LO)=ko
                     ENDIF
                     ni = ni + 1
                  END DO
                  ROMIN = ROMAX
                  ROMAX = 50.0D0
                  LO = LO+1
               END DO
!check               WRITE(*,*)NWTDNB
!check               write(*,*)LSTORD
!
!end order waters
               
!central QM water               
               write(n10,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),xtr(i),ytr(i),ztr(i)                                !print first QM/MM molecule
               do jw = 1,nwatsites-1                
                  write(n10,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),xtr(i+jw),ytr(i+jw),ztr(i+jw) 
               end do
!four nearest waters - QM neighbors               
!               do ik = 1,4
!                  kt = NWTDNB(ik)
!                  write(n10,*)atom_name_xyz(1),xtr(kt),ytr(kt),ztr(kt)
!                  do jw = 1,nwatsites-1   
!                     write(n10,*)atom_name_xyz(1+jw),xtr(kt+jw),ytr(kt+jw),ztr(kt+jw)
!                  end do
!               end do   
!Nw - 1 ordered water molecules
               do ik = 1,nmolwat-1
                  kt = LSTORD(ik)
                  write(n10,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),xtr(kt),ytr(kt),ztr(kt)
                  do jw = 1,nwatsites-1   
                     write(n10,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),xtr(kt+jw),ytr(kt+jw),ztr(kt+jw)
                  end do
               end do                  
!MM waters               
!               do iw=1,natms,nwatsites                                                    !print the remaining waters
!                  if((iw.ne.i).and.(iw.ne.NWTDNB(1)).and.(iw.ne.NWTDNB(2)).and.(iw.ne.NWTDNB(3)).and.(iw.ne.NWTDNB(4)))then
!                     write(n10,*)atom_name_xyz(1),xtr(iw),ytr(iw),ztr(iw)                       !oxygen
!                     do jw = 1,nwatsites-1                
!                        write(n10,*)atom_name_xyz(1+jw),xtr(iw+jw),ytr(iw+jw),ztr(iw+jw)         ! Hs and M
!                     end do
!                  endif   
!               end do

!check ordering
                do ik = 1,4
                   kt1 = NWTDNB(ik)
                   kt2 = LSTORD(ik) 
                   if(kt1.ne.kt2)then
                      write(*,*)'Error - tetrahedrality - QM/MM water ordering wrong'
                      write(*,*)ik,kt1,kt2
                      stop
                   endif   
                end do   
                nc = nc + 1                   !count number of QM/MM configurations already written nc_max = 500 (define in input if required)
            endif                  
!End Version 9                             
                             
! distributions amd means
!            nsh = IDNINT(3.0D0/SGDEL+SG/SGDEL+0.5D0)
            nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0
            NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh)+1
            AVWSG(nradsh) = AVWSG(nradsh) + SG
            AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
            SG = 0.0D0
! CALCULATE THE RO4 DISTRIBUTION AND AVERAGE VALUE
            nsh = RO4/RDEL+0.5D0
            NDPRO4(nradsh,nsh) = NDPRO4(nradsh,nsh) + 1
            AVRO4(nradsh) = AVRO4(nradsh) + RO4
            RO4 = 0.0D0
! CALCULATE THE RH2 DISTRIBUTION AND AVERAGE VALUE
            nsh = RH2/RDEL+0.5D0
            NDPRH2(nradsh,nsh) = NDPRH2(nradsh,nsh) + 1
            AVRH2(nradsh) = AVRH2(nradsh) + RH2
            RH2 = 0.0D0
! CALCULATE THE RHO1 DISTRIBUTION AND AVERAGE VALUE
            nsh = RHO1/RDEL+0.5D0
            NDPRHO1(nradsh,nsh) = NDPRHO1(nradsh,nsh) + 1
            AVRHO1(nradsh) = AVRHO1(nradsh) + RHO1
            RHO1 = 0.0D0
! CALCULATE THE RO DISTRIBUTIONS AND AVERAGE VALUES
            do KN = 1,5                                                          !run over five distinct oxygen-oxygen distances
               nsh = RO_nb(KN)/RDEL+0.5D0
               NDPRO(nradsh,KN,nsh) = NDPRO(nradsh,KN,nsh) + 1
               AVRO(nradsh,KN) = AVRO(nradsh,KN) + RO_nb(KN)
               RO_nb(KN) = 0.0D0       
            end do       
! End distributions and means                     

      end do                                                                        !end loop over water oxygens
   endif                                                                            !end sampling of waters
   
   if(lgas.and.(j==1.or.mod(j,tetsampw)==0))then
      kstps = kstps + 1
      write(n10,*)natms
      write(n10,'(a22,1x,2(1x,a23))')'# Tetrah water 1 = N/A','R_inter_O1-Ha (A) = N/A','R_inter_O1-Hb (A) = N/A'
!Single QM water 
      do i=1,natms,nwatsites
         write(n10,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),x(i),y(i),z(i)                                !print first QM/MM molecule
         do jw = 1,nwatsites-1                
            write(n10,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),x(i+jw),y(i+jw),z(i+jw) 
         end do
      end do   
   endif
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!Average number of waters on each environment

write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps

if(.not.lgas)then

   do kr = 1,NRSHMAX 
      write(*,*)      
      write(*,'(a15,1x,i4)')'#Environment = ',kr 
!   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep) 
      write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps) 
      write(n0,*)
      write(n0,'(a15,1x,i4)')'#Environment = ',kr 
!   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep)
      write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps)
   end do   

!Print distributions and mean values 

   write(n0,*)
   write(n0,*)'TETRAHEDRALITY AND DISTANCE DISTRIBUTIONS'
   write(n0,*)'========================================='
   write(n0,*)
        
   do kr = kr_start,NRSHMAX                                                      !loop over environments
      IF(kr==1)WRITE(n0,*) 
      IF(kr==1)WRITE(n0,*)'Bulk'
      IF(kr==1)WRITE(n0,*)'****'
      IF(kr==2)WRITE(n0,*)
      IF(kr==2)WRITE(n0,*)'Hydration Shell'
      IF(kr==2)WRITE(n0,*)'***************'
      IF(kr==3)WRITE(n0,*)
      IF(kr==3)WRITE(n0,*)'HSh-tetrahedral'
      IF(kr==3)WRITE(n0,*)'***************'
      IF(kr==4)WRITE(n0,*)
      IF(kr==4)WRITE(n0,*)'HSh-non-tetrahedral' 
      IF(kr==4)WRITE(n0,*)'*******************' 
   
! q tetrahedrality
      IF(kr==1)write(n1,'(a5)')'#Bulk'
      IF(kr==2)write(n1,'(a4)')'#HSh'
      IF(kr==3)write(n1,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n1,'(a15)')'#HSh-non-tetrah'   
! q VARIES BETWEEN -3 AND 1 FOR THE DEBENEDETTI NORMALIZATION IMPLEMENTED HERE
      NTDELS = 4.0D0/SGDEL
      DO IT = 1,NTDELS
         IF(NPDSG(kr,IT).GT.0)THEN
            TETHQ    = -3.01D0+SGDEL*DFLOAT(IT)
            TETHDIST = DFLOAT(NPDSG(kr,IT))/DFLOAT(kv(kr))                    !kv is the sum of the number of waters over the number of origins 
            WRITE(n1,9)TETHQ,TETHDIST
         ENDIF
      END DO
      WRITE(n1,*)
!
      AVWSG(kr)    = AVWSG(kr)/DFLOAT(kv(kr))
      AVCHHASG(kr) = AVCHHASG(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'Errington & Debenedetti perfect tetrahedron <q> = 1'
      WRITE(n0,*)'Average Errington & Debenedetti <q> = ',AVWSG(kr)
      WRITE(n0,*)
      WRITE(n0,*)'Chau & Hardwick perfect tetrahedron <q> = 0 '
      WRITE(n0,*)'Average Chau & Hardwick <Sg> =  <q> = ',AVCHHASG(kr)
      WRITE(n0,*)
   
! rO4 - mean over the 4 nearest oxygens 
      IF(kr==1)write(n2,'(a5)')'#Bulk'
      IF(kr==2)write(n2,'(a4)')'#HSh'
      IF(kr==3)write(n2,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n2,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRO4(kr,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRO4(kr,NKR))/DFLOAT(kv(kr))
            WRITE(n2,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n2,*)
!
      AVRO4(kr) = AVRO4(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RO_O4> Ang = ',AVRO4(kr)
      WRITE(n0,*)
   
! rOH2 = rH2 - mean over the 2 nearest hydrogens 
      IF(kr==1)write(n3,'(a5)')'#Bulk'
      IF(kr==2)write(n3,'(a4)')'#HSh'
      IF(kr==3)write(n3,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n3,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRH2(kr,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRH2(kr,NKR))/DFLOAT(kv(kr))
            WRITE(n3,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n3,*)
!
      AVRH2(kr) = AVRH2(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RO_H2> Ang = ',AVRH2(kr)
      WRITE(n0,*)

! rHO1 - mean over the nearest oxygen of each hydrogen 
      IF(kr==1)write(n4,'(a5)')'#Bulk'
      IF(kr==2)write(n4,'(a4)')'#HSh'
      IF(kr==3)write(n4,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n4,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRHO1(kr,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRHO1(kr,NKR))/DFLOAT(kv(kr))
            WRITE(n4,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n4,*)
!
      AVRHO1(kr) = AVRHO1(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RH_O1> Ang = ',AVRHO1(kr)
      WRITE(n0,*)

!rOi [i = 1,2,3,4,5] - mean over the nearest, the second nearest, the third nearest...etc oxygen
      do KN = 1,5
         IF(kr==1)write(n4+10*KN,'(a5)')'#Bulk'
         IF(kr==2)write(n4+10*KN,'(a4)')'#HSh'
         IF(kr==3)write(n4+10*KN,'(a11)')'#HSh-tetrah'
         IF(kr==4)write(n4+10*KN,'(a15)')'#HSh-non-tetrah'    
         DO NKR =1,NRDELS
            IF(NDPRO(kr,KN,NKR).GT.0)THEN
               RD_ANG = RDEL*DFLOAT(NKR)          
               DISTR_R=DFLOAT(NDPRO(kr,KN,NKR))/DFLOAT(kv(kr))
               WRITE(n4+10*KN,9)RD_ANG,DISTR_R
            ENDIF
         END DO
         WRITE(n4+10*KN,*)
!
         AVRO(kr,KN) = AVRO(kr,KN)/DFLOAT(kv(kr))
         WRITE(n0,*)
         WRITE(n0,*)
         WRITE(n0,*)'AVERAGE <RO_O',KN,'> Ang = ',AVRO(kr,KN)
         WRITE(n0,*)                    
      end do   

   end do            !end loop over environments

endif                !end lgas

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n4) 
close(n5) 
close(n6) 
close(n7) 
close(n8) 
close(n9) 
close(n10)
close(n11)

deallocate(NPDSG,NDPRO4,NDPRH2,NDPRHO1,NDPRO)
deallocate(AVWSG,AVCHHASG,AVRO4,AVRH2,AVRHO1)
deallocate(AVRO)
deallocate(NTDW4O)
deallocate(NONTDW)
deallocate(LSTORD)
deallocate(wd1)
deallocate(kv)

deallocate(xyz)
deallocate(x,y,z)


  9   FORMAT(27X,F9.4,5X,F14.5)
    return

END SUBROUTINE pure_wat_tetra



SUBROUTINE pure_wat_tetra_Env(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,tetsampw,R5Tet) 

           
! Calculate water's tetrahedrality in different environments - single enviroment for pure water

    integer,intent(in)                               :: npurew,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    integer,intent(in)                               :: tetsampw
    character(len=4),dimension(nwatsites)            :: atom_name_xyz 
    character(len=4),dimension(natms)                :: atomname
    real(kind=8),intent(in)                          :: R5Tet
!    integer,intent(in)                               :: nbulk

! Local variables
! Tetrahedrality and distance distributions
    real(kind=8)                                     :: SG,CHHASG,RO4
    real(kind=8)                                     :: RO_nb(5)
    real(kind=8)                                     :: RH2,RHO1
    integer                                          :: nsh
    real(kind=8)                                     :: SGDEL,RDEL
    integer                                          :: NWTDNB(4)
    REAL(KIND=4)                                     :: RNTDWH(2)
    integer,dimension(:,:),allocatable               :: NPDSG,NDPRO4,NDPRH2,NDPRHO1
    integer,dimension(:,:,:),allocatable             :: NDPRO
    real(kind=8),dimension(:),allocatable            :: AVWSG,AVCHHASG,AVRO4,AVRH2,AVRHO1
    real(kind=8),dimension(:,:),allocatable          :: AVRO
    integer,parameter                                :: nshdist = 50000                      !max. number of bins in tetrahedrality and distance dist.
    integer                                          :: NTDELS,NRDELS
    real(kind=8)                                     :: TETHQ,RD_ANG
    real(kind=8)                                     :: TETHDIST,DISTR_R 
! 
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=4),dimension(:),allocatable            :: xtr,ytr,ztr
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index
    integer,dimension(:),allocatable                 :: kv
    integer      :: n0, n1, n2, n3,n4, n5, n6, n7, n8, n9, n10, n11, n12
    integer      :: i, j, k
    integer      :: kw, iw, kt, kr, kvi, ki, jw, ik
    integer      :: iH, ihat, nt, ko, jt
    integer      :: kt1, kt2, ni, nc
    integer      :: io
    integer      :: L, KN, IT, NKR
    integer,parameter                                 :: INTMAX = 2147483647
    integer,dimension(:),allocatable                  :: NTDW4O                     !store the four nearest O atoms of each water molecule
    integer,dimension(:),allocatable                  :: NONTDW                     !store the waters next to the solute that are non-tetrahedral
    real(kind=8)                                      :: dr_bulk,dr4nb
    integer                                           :: NRSHMAX   
    integer                                           :: nradsh    
    integer                                           :: natcheck
    integer                                           :: kr_start
    integer                                           :: kstps
    logical                                           :: filex 
    integer                                           :: pwnmolsol,pwnatmsol,pwnions
!Version 9    
    real(kind=8)                                      :: ROMAX,ROMIN
    integer,dimension(:),allocatable                  :: LSTORD 
    integer                                           :: LO
    real(kind=4)                                      :: dx1,dy1,dz1,wdsq1
    real(kind=4),dimension(:),allocatable             :: wd1
    logical                                           :: lgas
!end Version 9

!v34 VARIABLES
    real(kind=8)                                 :: RADMAX,RADMIN
    real(kind=8)                                 :: dx,dy,dz,dr,dr2,dr5
    integer                                      :: KO5
    INTEGER                                      :: NTD5O(5)
!v34    
    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!Version 9
atom_name_xyz(1) ='O'
atom_name_xyz(2) ='H'
atom_name_xyz(3) ='H'
if(nwatsites ==4)atom_name_xyz(4) = 'MW'
!End Version 9

lgas=.false.
if(nmolwat==1)then
   lgas=.true.
   write(*,*)'The system is a gas - tetrahedrality calculation will be skipped'
endif   

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
allocate(xtr(natms),ytr(natms),ztr(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Tet_PW_out/log_PW_tet.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Tet_PW_out')   

pwnmolsol = 0
pwnatmsol = 0
pwnions = 0
    
n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140
n5 = 150
n6 = 160
n7 = 170
n8 = 180
n9 = 190
n10 = 200
n11 = 210
!n12 = 220

open(n0,file='Tet_PW_out/log_env_tetra.dat',status='unknown',action='write')
open(n1,file='Tet_PW_out/wat_tetra_dist.dat',status='unknown',action='write')
open(n2,file='Tet_PW_out/wat_rO4_dist.dat',status='unknown',action='write')
open(n3,file='Tet_PW_out/wat_rOH2_dist.dat',status='unknown',action='write')
open(n4,file='Tet_PW_out/wat_rHO1_dist.dat',status='unknown',action='write')
open(n5,file='Tet_PW_out/wat_rO1_nb_dist.dat',status='unknown',action='write')
open(n6,file='Tet_PW_out/wat_rO2_nb_dist.dat',status='unknown',action='write')
open(n7,file='Tet_PW_out/wat_rO3_nb_dist.dat',status='unknown',action='write')
open(n8,file='Tet_PW_out/wat_rO4_nb_dist.dat',status='unknown',action='write')
open(n9,file='Tet_PW_out/wat_rO5_nb_dist.dat',status='unknown',action='write')

!open(n10,file='Tet_PW_out/QMMM_conf.xyz',status='unknown',action='write')

!open(n11,file='Tet_PW_out/QMMM_check.xyz',status='unknown',action='write')
!v34
!open(n12,file='Tet_PW_out/wat_tetra_dist_env.dat',status='unknown',action='write')
!End v34
    
kr_start = 1                                     !environments loop started
!if(nbulk == 0) kr_start = 2                      !skip bulk to speed up calculations 

!NRSHMAX = 1                                      !maximum number of environments - single environment for pure water
NRSHMAX = 3                                      !Environments version
SGDEL   = 0.01D0                                 !tetrahedrality bin size 
RDEL    = 0.01D0                                 !bin size in distance distributions Ang
NRDELS  = 100.0d0/RDEL                            !maximum number of bins in distance distributions - set to 10 Ang

if(NRDELS>nshdist)then
   write(*,*)'error - routine pure_wat_tetra_Env - too many bins'
   stop
endif   

!Tetrahedrality and distances
allocate(NPDSG(NRSHMAX,nshdist),NDPRO4(NRSHMAX,nshdist),NDPRH2(NRSHMAX,nshdist),NDPRHO1(NRSHMAX,nshdist))
allocate(NDPRO(NRSHMAX,5,nshdist))
allocate(AVWSG(NRSHMAX),AVCHHASG(NRSHMAX),AVRO4(NRSHMAX),AVRH2(NRSHMAX),AVRHO1(NRSHMAX))
allocate(AVRO(NRSHMAX,5))
NPDSG=0;NDPRO4=0;NDPRH2=0;NDPRHO1=0
NDPRO=0
AVWSG=0.0d0;AVCHHASG=0.0d0;AVRO4=0.0d0;AVRH2=0.0d0;AVRHO1=0.0d0
AVRO=0.0d0

!Water identifiers
allocate(NTDW4O(4))  
allocate(NONTDW(nmolwat))   
allocate(kv(NRSHMAX))
allocate(LSTORD(nmolwat)) 
allocate(wd1(nmolwat))
NTDW4O = 0; NONTDW = 0; LSTORD = 0; wd1 = 0.0
kv = 0
kstps = 0
nc = 1                     !counter of the QM/MM (xyz) configurations

write(*,*)
write(*,*)'Water tetrahedrality'
write(*,*)'Routine pure_wat_tetra_Env'
write(*,*)'Results are printed to Tet_PW_out'
write(*,*)

! Start tetrahedrality calculation
write(*,*)
write(*,*)'Starting tetrahedrality calculation...'
write(*,*)
write(n0,*)
write(n0,*)'Starting tetrahedrality calculation...'
write(n0,*)

do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
   NTDW4O = 0; NONTDW = 0
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!check      
!      write(n11,*)natms
!      write(n11,*)
!      do i = 1,natms,nwatsites
!         write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),x(i),y(i),z(i)                                !print first QM/MM molecule
!         do jw = 1,nwatsites-1                
!            write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),x(i+jw),y(i+jw),z(i+jw) 
!         end do
!      end do   
      
!End Version 9

   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
   
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
!check       write(*,*)x(kp),y(kp),z(kp)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules
!sample frequency - calculate the tetrahedrality every tetsampw time-steps  

   if(.not.lgas.and.(j==1.or.mod(j,tetsampw)==0))then
      write(*,*)'time-step',j,' sampling water molecules'
      kstps = kstps + 1                                                                            !count number of origins sampled

!Version 9 check
!      write(n11,*)natms
!      write(n11,*)'# comment line',j
!End version 9 check
      
      io = 0
      do i=1,natms,nwatsites                                 !Loop over water oxygens
         io = io + 1                                         !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!bulk water                                                   
            nradsh = 1                                                                              !pure water environment
            kv(nradsh) = kv(nradsh) + 1                                                             !number of waters
!NG            call tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
!NG                             SG,CHHASG,RO4,RO_nb,RH2,RHO1)
            
            call tetrahedron(i,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,cell,&
                             SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)  
                             
                             
!v34
! Find the 5 closest oxygen atoms to the oxygen i
            RADMAX  = 25.0D0
            RADMIN  =  0.0D0
            JT      = 1
            DO WHILE(JT.LE.5)
               do jw=1,natms,nwatsites                                       !Loop over water oxygens
                  IF(jw.NE.i)THEN
            ! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND OXYGEN ATOMS jw     
                     dx = x(i) - x(jw)             
                     dy = y(i) - y(jw)
                     dz = z(i) - z(jw)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dr2 = dx**2 + dy**2 + dz**2
                     dr  = dsqrt(dr2)        
            !FIND THE NEAREST 5 OXYGEN ATOMS jw TO THE OXYGEN ATOM i
                     IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
                        RADMAX = dr
            !NTD5O STORES THE O ATOMS NEIGHBORS OF I
            !JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR ; JT = 3 THIRD NEIGHBOR ; JT = 4 FOURTH NEIGHBOR ; JT = 5 FIFTH NEIGHBOR
                        NTD5O(JT) = jw
                     ENDIF
                  ENDIF
               end do
               JT = JT + 1
               RADMIN = RADMAX
               RADMAX = 25.0D0
            END DO    
            
!distance to the fifth neighbor         
            KO5 = NTD5O(5)
            dx  = x(i) - x(KO5)             
            dy  = y(i) - y(KO5)
            dz  = z(i) - z(KO5)
            dx  = dx - dnint(dx/cell(1))*cell(1)
            dy  = dy - dnint(dy/cell(2))*cell(2)
            dz  = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2) 
            dr5 = dr
!End v34
                             
! distributions amd means
! Tetrahedrality
!            nsh = IDNINT(3.0D0/SGDEL+SG/SGDEL+0.5D0)
            nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0
            NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh)+1
            AVWSG(nradsh) = AVWSG(nradsh) + SG
            AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
!v34            SG = 0.0D0
! CALCULATE THE RO4 DISTRIBUTION AND AVERAGE VALUE
            nsh = RO4/RDEL+0.5D0
            NDPRO4(nradsh,nsh) = NDPRO4(nradsh,nsh) + 1
            AVRO4(nradsh) = AVRO4(nradsh) + RO4
            RO4 = 0.0D0
! CALCULATE THE RH2 DISTRIBUTION AND AVERAGE VALUE
            nsh = RH2/RDEL+0.5D0
            NDPRH2(nradsh,nsh) = NDPRH2(nradsh,nsh) + 1
            AVRH2(nradsh) = AVRH2(nradsh) + RH2
            RH2 = 0.0D0
! CALCULATE THE RHO1 DISTRIBUTION AND AVERAGE VALUE
            nsh = RHO1/RDEL+0.5D0
            NDPRHO1(nradsh,nsh) = NDPRHO1(nradsh,nsh) + 1
            AVRHO1(nradsh) = AVRHO1(nradsh) + RHO1
            RHO1 = 0.0D0
! CALCULATE THE RO DISTRIBUTIONS AND AVERAGE VALUES
            do KN = 1,5                                                          !run over five distinct oxygen-oxygen distances
               nsh = RO_nb(KN)/RDEL+0.5D0
               NDPRO(nradsh,KN,nsh) = NDPRO(nradsh,KN,nsh) + 1
               AVRO(nradsh,KN) = AVRO(nradsh,KN) + RO_nb(KN)
               RO_nb(KN) = 0.0D0       
            end do     
            
!v34
!Tetrahedrality environments
            nsh = 3.0D0/SGDEL+SG/SGDEL+0.5D0
            if(dr5.gt.R5Tet)then
               nradsh = 2                                                                              !environment 2
               kv(nradsh) = kv(nradsh) + 1                                                             !number of waters
               NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh)+1
               AVWSG(nradsh) = AVWSG(nradsh) + SG
               AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG
            elseif(dr5.lt.R5Tet)then
               nradsh = 3                                                                              !environment 3
               kv(nradsh) = kv(nradsh) + 1                                                             !number of waters
               NPDSG(nradsh,nsh) = NPDSG(nradsh,nsh)+1
               AVWSG(nradsh) = AVWSG(nradsh) + SG
               AVCHHASG(nradsh) = AVCHHASG(nradsh) + CHHASG               
            endif
            SG = 0.0D0
!End v34            
            
! End distributions and means                     

      end do                                                                        !end loop over water oxygens
   endif                                                                            !end sampling of waters
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!Average number of waters on each environment

write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps

if(.not.lgas)then

   do kr = 1,NRSHMAX 
      write(*,*)      
      write(*,'(a15,1x,i4)')'#Environment = ',kr 
!   write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep) 
      write(*,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps) 
      write(n0,*)
      write(n0,'(a15,1x,i4)')'#Environment = ',kr 
!   write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(nstep)
      write(n0,'(a25,1x,F9.2)')'#Mean number of waters = ',dfloat(kv(kr))/dfloat(kstps)
   end do   

!Print distributions and mean values 

   write(n0,*)
   write(n0,*)'TETRAHEDRALITY AND DISTANCE DISTRIBUTIONS'
   write(n0,*)'========================================='
   write(n0,*)
        
   do kr = kr_start,NRSHMAX                                                      !loop over environments
      IF(kr==1)WRITE(n0,*) 
      IF(kr==1)WRITE(n0,*)'Bulk'
      IF(kr==1)WRITE(n0,*)'****'
      IF(kr==2)WRITE(n0,*)
      IF(kr==2)WRITE(n0,*)'No Interstitial Waters'
      IF(kr==2)WRITE(n0,*)'**********************'
      IF(kr==3)WRITE(n0,*)
      IF(kr==3)WRITE(n0,*)'Interstitial Waters'
      IF(kr==3)WRITE(n0,*)'*******************'
      IF(kr==4)WRITE(n0,*)
      IF(kr==4)WRITE(n0,*)'Undefined Water' 
      IF(kr==4)WRITE(n0,*)'***************' 
   
! q tetrahedrality
      IF(kr==1)write(n1,'(a5)')'#Bulk'
      IF(kr==2)write(n1,'(a4)')'#RO4'
      IF(kr==3)write(n1,'(a4)')'#RO5'
      IF(kr==4)write(n1,'(a9)')'#Undefined'   
! q VARIES BETWEEN -3 AND 1 FOR THE DEBENEDETTI NORMALIZATION IMPLEMENTED HERE
      NTDELS = 4.0D0/SGDEL
      DO IT = 1,NTDELS
         IF(NPDSG(kr,IT).GT.0)THEN
            TETHQ    = -3.01D0+SGDEL*DFLOAT(IT)
            TETHDIST = DFLOAT(NPDSG(kr,IT))/DFLOAT(kv(kr))                    !kv is the sum of the number of waters over the number of origins 
            WRITE(n1,9)TETHQ,TETHDIST
!            WRITE(n12,9)TETHQ,TETHDIST
         ENDIF
      END DO
      WRITE(n1,*)
!      WRITE(n12,*)
      AVWSG(kr)    = AVWSG(kr)/DFLOAT(kv(kr))
      AVCHHASG(kr) = AVCHHASG(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'Errington & Debenedetti perfect tetrahedron <q> = 1'
      WRITE(n0,*)'Average Errington & Debenedetti <q> = ',AVWSG(kr)
      WRITE(n0,*)
      WRITE(n0,*)'Chau & Hardwick perfect tetrahedron <q> = 0 '
      WRITE(n0,*)'Average Chau & Hardwick <Sg> =  <q> = ',AVCHHASG(kr)
      WRITE(n0,*)
   end do   
!
!   kr_start = 1
   NRSHMAX = 1
   do kr = kr_start,NRSHMAX      
! rO4 - mean over the 4 nearest oxygens 
      IF(kr==1)write(n2,'(a5)')'#Bulk'
      IF(kr==2)write(n2,'(a4)')'#HSh'
      IF(kr==3)write(n2,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n2,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRO4(kr,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRO4(kr,NKR))/DFLOAT(kv(kr))
            WRITE(n2,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n2,*)
!
      AVRO4(kr) = AVRO4(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RO_O4> Ang = ',AVRO4(kr)
      WRITE(n0,*)
   
! rOH2 = rH2 - mean over the 2 nearest hydrogens 
      IF(kr==1)write(n3,'(a5)')'#Bulk'
      IF(kr==2)write(n3,'(a4)')'#HSh'
      IF(kr==3)write(n3,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n3,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRH2(kr,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRH2(kr,NKR))/DFLOAT(kv(kr))
            WRITE(n3,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n3,*)
!
      AVRH2(kr) = AVRH2(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RO_H2> Ang = ',AVRH2(kr)
      WRITE(n0,*)

! rHO1 - mean over the nearest oxygen of each hydrogen 
      IF(kr==1)write(n4,'(a5)')'#Bulk'
      IF(kr==2)write(n4,'(a4)')'#HSh'
      IF(kr==3)write(n4,'(a11)')'#HSh-tetrah'
      IF(kr==4)write(n4,'(a15)')'#HSh-non-tetrah'    
      DO NKR =1,NRDELS
         IF(NDPRHO1(kr,NKR).GT.0)THEN
            RD_ANG = RDEL*DFLOAT(NKR)          
            DISTR_R=DFLOAT(NDPRHO1(kr,NKR))/DFLOAT(kv(kr))
            WRITE(n4,9)RD_ANG,DISTR_R
         ENDIF
      END DO
      WRITE(n4,*)
!
      AVRHO1(kr) = AVRHO1(kr)/DFLOAT(kv(kr))
      WRITE(n0,*)
      WRITE(n0,*)
      WRITE(n0,*)'AVERAGE <RH_O1> Ang = ',AVRHO1(kr)
      WRITE(n0,*)

!rOi [i = 1,2,3,4,5] - mean over the nearest, the second nearest, the third nearest...etc oxygen
      do KN = 1,5
         IF(kr==1)write(n4+10*KN,'(a5)')'#Bulk'
         IF(kr==2)write(n4+10*KN,'(a4)')'#HSh'
         IF(kr==3)write(n4+10*KN,'(a11)')'#HSh-tetrah'
         IF(kr==4)write(n4+10*KN,'(a15)')'#HSh-non-tetrah'    
         DO NKR =1,NRDELS
            IF(NDPRO(kr,KN,NKR).GT.0)THEN
               RD_ANG = RDEL*DFLOAT(NKR)          
               DISTR_R=DFLOAT(NDPRO(kr,KN,NKR))/DFLOAT(kv(kr))
               WRITE(n4+10*KN,9)RD_ANG,DISTR_R
            ENDIF
         END DO
         WRITE(n4+10*KN,*)
!
         AVRO(kr,KN) = AVRO(kr,KN)/DFLOAT(kv(kr))
         WRITE(n0,*)
         WRITE(n0,*)
         WRITE(n0,*)'AVERAGE <RO_O',KN,'> Ang = ',AVRO(kr,KN)
         WRITE(n0,*)                    
      end do   

   end do            !end loop over environments

endif                !end lgas

close(n0) 
close(n1) 
close(n2) 
close(n3) 
close(n4) 
close(n5) 
close(n6) 
close(n7) 
close(n8) 
close(n9) 
close(n10)
close(n11)
!close(n12)

deallocate(NPDSG,NDPRO4,NDPRH2,NDPRHO1,NDPRO)
deallocate(AVWSG,AVCHHASG,AVRO4,AVRH2,AVRHO1)
deallocate(AVRO)
deallocate(NTDW4O)
deallocate(NONTDW)
deallocate(LSTORD)
deallocate(wd1)
deallocate(kv)

deallocate(xyz)
deallocate(x,y,z)


  9   FORMAT(27X,F9.4,5X,F14.5)
    return

END SUBROUTINE pure_wat_tetra_Env



!Start Voronoi Laboratory
SUBROUTINE pure_wat_Voronoi(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,Voronsampw,rcut)
           
! Calculate water Voronoi polyhedra in different environments - single enviroment for pure water
! This is based on the algorithm available at ccp5 library f35 - http://www.ccp5.ac.uk/software/allen_tildersley
! A second method was implemented based on the description provided in the paper: Kazmierczak and Swiatla-Wojcik RSC Adv., 2014, 4, 41812
! This second method has some unknow problem - check

    integer,intent(in)                               :: npurew,nbox
    real(kind=4),intent(in)                          :: rcut
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    integer,intent(in)                               :: Voronsampw
    character(len=4),dimension(nwatsites)            :: atom_name_xyz 
    character(len=4),dimension(natms)                :: atomname
    
    real(kind=4)                                     :: Area_Vor, Vol_Vor

! Local variables

!NG    real(kind=8)                                   :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index
    integer                                          :: kv
    integer      :: n0, n1, n2, n3, n4, n5, n6
    integer      :: i, j, k
    integer      :: kw, iw, kt, kr, kvi, ki, jw, ik
    integer      :: iH, ihat, nt, ko, jt
    integer      :: kt1, kt2, ni, nc
    integer      :: io,jo
    integer      :: isort
    integer      :: L, KN, IT, NKR
    integer,parameter                                 :: INTMAX = 2147483647
    real(kind=4),parameter                            :: PI = 3.141519
    integer                                           :: natcheck
    integer                                           :: kr_start
    integer                                           :: kstps, kskip
    logical                                           :: filex 
    integer                                           :: pwnmolsol,pwnatmsol,pwnions

!Local Variables
integer                                            :: MAXVER,MAXCAN,MAXN
real(kind=8)                                       :: dx,dy,dz,dr,dr2   
!NG real(kind=4)                                       :: BOX,RCUTSQ, BOXINV
real(kind=4)                                       :: RCUTSQ,COORD
integer                                            :: INAB, JNAB
logical                                            :: print_voron
!
integer,dimension(:,:),allocatable                 :: NABLST
integer,dimension(:),allocatable                   :: NNAB
real(kind=4),dimension(:),allocatable              :: PX, PY, PZ, PS
integer,dimension(:),allocatable                   :: TAG, EDGES
real(kind=4),dimension(:),allocatable              :: RXVER, RYVER, RZVER
integer,dimension(:),allocatable                   :: IVER, JVER, KVER
integer                                            :: CAN, VER                                                !candidates
INTEGER                                            :: NCAN, NVER, NEDGE, NFACE, NCOORD
LOGICAL                                            :: OK

!Vertices distribution
integer,parameter                                  :: nshdist = 50000                      !max. number of bins
integer                                            :: KRT
real(kind=4)                                       :: VDEL, ADEL
integer                                            :: nsh, NVDELS, NADELS
real(kind=4)                                       :: AVVRT, AVVOL, AVARA, AVASP
integer,dimension(:),allocatable                   :: NDPVRT, NDPARA, NDPVOL, NDPASP
real(kind=4)                                       :: AVVRT_, AVVOL_, AVARA_, AVASP_
integer,dimension(:),allocatable                   :: NDPVRT_, NDPARA_, NDPVOL_, NDPASP_
real(kind=4)                                       :: DISTR_V, VD_NUMB, ASPHER
integer                                            :: NDPASP_SUM
real(kind=4),dimension(:),allocatable              :: xn, yn, zn
real(kind=4)                                       :: xi, yi, zi
real(kind=4)                                       :: norm_i, norm_j
real(kind=4),dimension(:),allocatable              :: PD 
logical                                            :: LEULER
            
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!rcut   = 6.5                      !Neighbors cut-off in Ang  [optimized >= 8.0]  - moved to input                
MAXVER = 100                    !Maximum vertices
MAXCAN = 1000                    !maximum candidates
MAXN   = nmolwat*nwatsites
RCUTSQ = RCUT**2.0
VDEL   = 0.5                      !VERTICE DISTRIBUTION BIN SIZE IN NUMBER OF VERTICES
NVDELS = 200.0d0/VDEL            !maximum number of bins

ADEL   = 0.01 
NADELS = 200.0d0/ADEL
print_voron=.false.

if( NVDELS > nshdist .or. NADELS > nshdist )then
   write(*,*)'Error - routine pure_wat_Voronoi - too many bins'
   stop
endif  

pwnmolsol = 0
pwnatmsol = 0
pwnions   = 0

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

!Local variables
allocate(NABLST(MAXVER,MAXN))
allocate(NNAB(MAXN))
allocate(PX(MAXCAN), PY(MAXCAN), PZ(MAXCAN), PS(MAXCAN), PD(MAXCAN))
allocate(TAG(MAXCAN), EDGES(MAXCAN))
allocate(RXVER(MAXVER), RYVER(MAXVER), RZVER(MAXVER))
allocate(IVER(MAXVER), JVER(MAXVER), KVER(MAXVER))
allocate(xn(maxcan), yn(maxcan), zn(maxcan))
allocate( NDPVRT(nshdist), NDPVOL(nshdist), NDPARA(nshdist), NDPASP(nshdist) )
allocate( NDPVRT_(nshdist), NDPVOL_(nshdist), NDPARA_(nshdist), NDPASP_(nshdist) )

AVVRT = 0; AVVOL = 0; AVARA = 0; NDPVRT = 0; NDPVOL = 0; NDPARA = 0; NDPASP = 0; ASPHER = 0; AVASP = 0
AVVRT_= 0; AVVOL_= 0; AVARA_= 0; NDPVRT_= 0; NDPVOL_= 0; NDPARA_= 0; NDPASP_=0;              AVASP_= 0

PX = 0.0; PY = 0.0; PZ = 0.0
NDPASP_SUM = 0

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
inquire(file='Voronoi_out/log_Voronoi.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p Voronoi_out')   

n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140
n5 = 150

open(n0,file='Voronoi_out/log_Voronoi.dat',status='unknown',action='write')
open(n1,file='Voronoi_out/Voronoi_Asphericity.dat',status='unknown',action='write')
open(n2,file='Voronoi_out/Voronoi_Vertex_dist.dat',status='unknown',action='write')
open(n3,file='Voronoi_out/Voronoi_Area_dist.dat',status='unknown',action='write')
open(n4,file='Voronoi_out/Voronoi_Volume_dist.dat',status='unknown',action='write')

   
!kr_start = 1                                     !environments loop started

kstps = 0
kskip = 0

write(*,*)
write(*,*)'Water Voronoi polyhedra'
write(*,*)'Routine pure_wat_Voronoi'
write(*,*)'Results are printed to Wat_Voronoi.dat'
write(*,*)'Voronoi cut-off set to',rcut,' Ang'
write(*,*)

! Start Voronoi calculation
write(*,*)
write(*,*)'Starting Voronoi calculation...'
write(*,*)
write(n0,*)
write(n0,*)'Starting Voronoi calculation...'
write(n0,*)'Voronoi cut-off set to',rcut,' Ang'
write(n0,*)
write(n2,'(A30)')'#VORONOI vertices distribution'
write(n3,'(A30)')'#VORONOI surface  distribution'
write(n4,'(A30)')'#VORONOI volumes  distribution'



do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
   LEULER = .true.
   NDPVRT_ = 0
   AVVRT_  = 0.0
   NDPASP_ = 0
   AVASP_  = 0.0
   NDPARA_ = 0
   AVARA_  = 0.0
   NDPVOL_ = 0
   AVVOL_  = 0.0

!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!read solute and water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!Version 9 - CP2K x,y,z - waters outside the box; apply pbc
      if(inputformat.eq.'BOMD   ')then
         do jw=1,natms,nwatsites
            if(x(jw).lt.-cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)+cell(1)
               end do
            endif
            if(x(jw).gt.cell(1)/2.0)then
               do jt=0,nwatsites-1
                  x(jw+jt)=x(jw+jt)-cell(1)
               end do
            endif
            
            if(y(jw).lt.-cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)+cell(2)
               end do
            endif
            if(y(jw).gt.cell(2)/2.0)then
               do jt=0,nwatsites-1
                  y(jw+jt)=y(jw+jt)-cell(2)
               end do
            endif
            
            if(z(jw).lt.-cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)+cell(3)
               end do
            endif
            if(z(jw).gt.cell(3)/2.0)then
               do jt=0,nwatsites-1
                  z(jw+jt)=z(jw+jt)-cell(3)
               end do
            endif                 
         end do
      endif  !end cp2k BOMD PBC
!check      
!      write(n11,*)natms
!      write(n11,*)
!      do i = 1,natms,nwatsites
!         write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1),x(i),y(i),z(i)                                !print first QM/MM molecule
!         do jw = 1,nwatsites-1                
!            write(n11,'(a3,2x,3(f7.3,2x))')atom_name_xyz(1+jw),x(i+jw),y(i+jw),z(i+jw) 
!         end do
!      end do   
      
!End Version 9

   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
   
    cell(1)=dble(box(1))
    cell(2)=dble(box(2)) 
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
!check       write(*,*)x(kp),y(kp),z(kp)
       kp = kp + 1
    end do   

   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 

!sample water molecules
!sample frequency - calculate Voronoi every Voronsampw time-steps  

!   if(j==1.or.mod(j,Voronsampw)==0.and.j /= 40)then
   if(j==1.or.mod(j,Voronsampw)==0)then
      write(*,*)'time-step',j,' sampling water molecules'
      kstps = kstps + 1                                                                            !count number of origins sampled
      
      NABLST = 0
      NNAB   = 0
      
      io = 0
      do i=1,natms,nwatsites                                 !Loop over water oxygens
         io = io + 1                                         !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
!NG Keep coordinates of central particle         
         xi = x(i)
         yi = y(i)
         zi = z(i)
!         
!         call VORON3(i,natms,pwnmolsol,pwnatmsol,nmolwat,nwatsites,pwnions,x,y,z,MAXVER,MAXCAN,MAXN,cell,rcut)
!
! Put subroutine VORON3 below

!NG   *************************   Routine VORON3   ************************* 

         CAN    = 0

!                             ** SELECT CANDIDATES **

!do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         jo = 0
         do jw = 1,natms,nwatsites        !Loop over water oxygens
            jo = jo + 1                   !Oxygen atoms number - jo = 1,2,3,4,5,...
            IF(jw.NE.i)THEN  
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND OXYGEN ATOMS jw     
!               dx = x(i) - x(jw)             
!               dy = y(i) - y(jw)
!               dz = z(i) - z(jw)
               dx = x(jw) - x(i)              
               dy = y(jw) - y(i)
               dz = z(jw) - z(i) 
               dx = dx - dnint(dx/cell(1))*cell(1)
               dy = dy - dnint(dy/cell(2))*cell(2)
               dz = dz - dnint(dz/cell(3))*cell(3)
               dr2 = dx**2 + dy**2 + dz**2               
!      dr  = dsqrt(dr2)
               IF ( dr2 .LT. RCUTSQ ) THEN      
                  CAN = CAN + 1
                  IF ( CAN .GT. MAXCAN ) THEN
                     WRITE(*,'('' TOO MANY CANDIDATES '')')
                     STOP
                  ENDIF
!check ok         write(*,*) i,jw,sqrt(dr2)
!WARNING!!! CAN takes values 1, 2, 3,... etc not the oxygen index (i)
                  PX(CAN)  = dx
                  PY(CAN)  = dy
                  PZ(CAN)  = dz
                  PS(CAN)  = dr2
!                  TAG(CAN) = jw                   !WARNING TAG(CAN) equal to the real oxygen atom number (1, 4, 7, ...) not 1,2,3,4, ...
                  TAG(CAN) = jo
!check                  write(*,*)j,i,CAN,sqrt(PS(CAN))

!NG   Lab - THE METHOD IS CHOSEN IN THE ROUTINE WORK
!METHOD 2 - THERE IS SOME PROBLEM - EQUATIONS SHOULD BE INCORRECT
! Voronoi - Method 2 - bisector plane 
                  norm_i  =  x(i)**2.0  + y(i)**2.0  + z(i)**2.0 
                  norm_j  = x(jw)**2.0 + y(jw)**2.0 + z(jw)**2.0 
                  PD(CAN) = 0.5*(norm_i - norm_j)
! Keep coordinates of nearest neighbors                  
                  xn(CAN)=x(jw)
                  yn(CAN)=y(jw)
                  zn(CAN)=z(jw)
!NG   End Lab         
               ENDIF
            ENDIF   
         end do
         
!check         write(*,*)j,i,CAN
         
!                             ** CANDIDATES HAVE BEEN SELECTED **

         NCAN = CAN

!                             ** SORT INTO ASCENDING ORDER OF DISTANCE **
!                             ** THIS SHOULD IMPROVE EFFICIENCY        **

         CALL SORT ( MAXCAN, PX, PY, PZ, PS, TAG, NCAN, xn, yn, zn, PD )                          !called for each i         
         
!check sort ok         do isort= 1,NCAN
!check sort ok            write(*,*)j,i,sqrt(PS(isort))
!check sort ok         end do   
         

!                             ** PERFORM VORONOI ANALYSIS **

!!called for each i
         if( NCAN < 5 ) then
            write(*,'(a19,1x,i5,1x,a24,f5.1,1x,a3)')'Warning!!! Molecule', i, 'has < 4 neighbors within', rcut, 'Ang'
            stop
         endif
         
CALL WORK ( i,MAXCAN,MAXVER,NCAN,NVER,NEDGE,NFACE,PX,PY,PZ,PS,EDGES,RXVER,RYVER,RZVER,IVER,JVER,KVER,Area_Vor,Vol_Vor,xi,yi,zi,xn,yn,zn,PD,LEULER )   
         
!NG   *************************   End Routine VORON3   *************************

!NG START DISTRIBUTION FUNCTIONS CALCULATION
! RXVER, RYVER, RZVER - vertices coordinates
! NVER - Number of vertices (total) of the polyhedron of molecule jw

!PRELIMINARY DISTRIBUTIONS
         if(LEULER)then                                    !IF A VALID POLYHEDRON WAS FOUND
! CALCULATE THE VERTICES DISTRIBUTION AND AVERAGE VALUE
            nsh = dfloat(NVER)/ADEL + 0.5D0
            NDPVRT_(nsh) = NDPVRT_(nsh) + 1
            AVVRT_ = AVVRT_ + NVER
! CALCULATE THE TOTAL AREA DISTRIBUTION AND AVERAGE VALUE         
            nsh = Area_Vor/VDEL + 0.5D0
            NDPARA_(nsh) = NDPARA_(nsh) + 1
            AVARA_ = AVARA_ + Area_Vor         
!            write(*,*)i,Area_Vor
! CALCULATE THE VOLUME DISTRIBUTION AND AVERAGE VALUE         
            nsh = Vol_Vor/VDEL + 0.5D0
            NDPVOL_(nsh) = NDPVOL_(nsh) + 1
            AVVOL_ = AVVOL_ + Vol_Vor
! CALCULATE THE ASPHERICITY DISTRIBUTION AND AVERAGE VALUE    
            ASPHER = (Area_Vor**3.0)/(36.0*pi*Vol_Vor**2.0)
            nsh = ASPHER/ADEL + 0.5D0
            NDPASP_(nsh) = NDPASP_(nsh) + 1
            AVASP_ = AVASP_ + ASPHER
         else
            kstps = kstps - 1
            kskip = kskip + 1
            write(*,*)'Step ',j,' skipped' 
            goto 1150
         endif
         
!End preliminary dsitributions         
!DISTRIBUTIONS - check if a polyhedron was found for each molecule        
         if( io == nmolwat .and. LEULER )then
! CALCULATE THE VERTICES AND THE ASPHERICITY DISTRIBUTION AND AVERAGE VALUE
            do NKR =1,NADELS                   
               NDPVRT(NKR) = NDPVRT(NKR) + NDPVRT_(NKR)               
               NDPASP(NKR) = NDPASP(NKR) + NDPASP_(NKR)               
!              NDPASP_SUM = NDPASP_SUM + NDPASP(NKR)               
            end do
            AVVRT = AVVRT + AVVRT_
            AVASP = AVASP + AVASP_
! CALCULATE THE TOTAL AREA AND THE VOLUME DISTRIBUTION AND AVERAGE VALUE
            do NKR =1,NVDELS            
               NDPARA(NKR) = NDPARA(NKR) + NDPARA_(NKR)                
               NDPVOL(NKR) = NDPVOL(NKR) + NDPVOL_(NKR)
            end do
            AVARA = AVARA + AVARA_
            AVVOL = AVVOL + AVVOL_
         endif 
         
!NG END DISTRIBUTION FUNCTIONS CALCULATION

!                             ** WRITE OUT RESULTS **

!DEFAULT is to set print_voron F - this is a gigantic output

         if(print_voron)then
         
            WRITE(n0,'(/1X,''CENTRAL MOLECULE '',2(I5))') io, i
            WRITE(n0,'(/1X,''NUMBER OF NEIGHBOURS = NUMBER OF FACES '',I5)') NFACE
            WRITE(n0,'(/1X,''NEIGHBOUR LIST '')')
            WRITE(n0,10001)
!10001   FORMAT(/1X,'ATOM ',3X,'FACE ',/1X,'INDEX',3X,'EDGES',3X,'            RELATIVE POSITION         ',3X,'  DISTANCE')         

!NG   THE FOLLOWING "RELATIVE POSITION" IS THE POSITION OF THE NEIGHBOR PX(CAN), PY(CAN), PZ(CAN) AND THE DISTANCE IN Ang PS(CAN)
            DO CAN = 1, NCAN
               IF (EDGES(CAN) .NE. 0) THEN
                  PS(CAN) = SQRT ( PS(CAN) )
                  WRITE(n0,'(1X,I5,3X,I5,3X,3F12.5,3X,F12.5)')TAG(CAN), EDGES(CAN), PX(CAN), PY(CAN), PZ(CAN), PS(CAN)
!                  NNAB(J) = NNAB(J) + 1          !original routine J = i
!                  NABLST(NNAB(J),J) = TAG(CAN)
            
                  NNAB(io) = NNAB(io) + 1                        !NG   io is the central particle of the VORONOI polyhedral
                  NABLST(NNAB(io),io) = TAG(CAN)    
!write(*,*)      i,io,NNAB(io)               
               ENDIF
            END DO
      
            WRITE(n0,'(/1X,''NUMBER OF EDGES '',I5)') NEDGE
            WRITE(n0,'(/1X,''NUMBER OF VERTICES '',I5)') NVER
            WRITE(n0,'(/1X,''VERTEX LIST '')')
            WRITE(n0,10002)
!10002   FORMAT(/1X,'      INDICES           RELATIVE POSITION')         
            DO VER = 1, NVER
            
               WRITE(n0,'(1X,3I5,3X,3F12.5)') TAG(IVER(VER)), TAG(JVER(VER)), TAG(KVER(VER)), RXVER(VER), RYVER(VER), RZVER(VER)
            
            END DO
            
         endif      !end if print_voron
      end do !End Loop over water oxygens i
     
!    *******************************************************************
!    ** MAIN LOOP ENDS                                                **
!    *******************************************************************

      if(print_voron) WRITE(n0,'(1H1,''FINAL SUMMARY'')')
      if(print_voron) WRITE(n0,10003)
!10003   FORMAT(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ') 

      NCOORD = 0

!      DO J = 1, N          !original
      ko = 0
      DO K = 1,natms,nwatsites

         ko = ko + 1
         NCOORD = NCOORD + NNAB(ko)
         
         if(print_voron) WRITE(n0,'(1X,I5,3X,I5,3X,30I4)') ko, NNAB(ko),( NABLST(INAB,ko), INAB = 1, NNAB(ko) )

!        ** CHECK THAT IF I IS A NEIGHBOUR OF K **
!        ** THEN K IS ALSO A NEIGHBOUR OF I     **

         DO INAB = 1, NNAB(ko)
            
            I = NABLST(INAB,ko)
!            write(*,*)ko, NNAB(ko),I
            OK = .FALSE.
            JNAB = 1
!            JNAB = 0       !NG   THIS IS IN THE ORIGINAL CODE - IT GIVES AN ERROR - THE NABLST INDEX DOES NOT START AT 0 

1200        IF ( ( .NOT. OK ) .AND. ( JNAB .LE. NNAB(I) ) ) THEN
               OK = ( ko .EQ. NABLST(JNAB,I) )
               JNAB = JNAB + 1
               GOTO 1200
            ENDIF

            IF ( .NOT. OK ) THEN
               WRITE(n0,'(1X,I3,'' IS NOT A NEIGHBOUR OF '',I3)') ko, I
               WRITE(* ,'(1X,I3,'' IS NOT A NEIGHBOUR OF '',I3)') ko, I
            ENDIF
    
         END DO
    
      END DO
    
      COORD = REAL ( NCOORD ) / REAL ( nmolwat )
    
      if(print_voron) WRITE(n0,'(/1X,'' AVERAGE COORDINATION NUMBER = '',F10.5)') COORD
!NG     
10001   FORMAT(/1X,'ATOM ',3X,'FACE ',/1X,'INDEX',3X,'EDGES',3X,'            RELATIVE POSITION         ',3X,'  DISTANCE')
10002   FORMAT(/1X,'      INDICES           RELATIVE POSITION')
10003   FORMAT(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ')        
         
!          end do                                                                        !end loop over water oxygens
      endif                                                                            !end sampling of waters
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
1150   continue 
end do                                                                              !end time-step main loop


!NG Averages
!Vertices
do NKR =1,NADELS
   if(NDPVRT(NKR) > 0)then
      VD_NUMB = ADEL*dfloat(NKR)          
      DISTR_V = dfloat(NDPVRT(NKR))/(dfloat(kstps)*dfloat(nmolwat))
      write(n2,9)VD_NUMB,DISTR_V
   endif
end do
AVVRT = AVVRT/(dfloat(kstps)*dfloat(nmolwat))
write(n0,*)'Average number of vertices = ', AVVRT
write(n2,*)

!Area
do NKR =1,NVDELS
   if(NDPARA(NKR) > 0)then
      VD_NUMB = VDEL*dfloat(NKR)          
      DISTR_V = dfloat(NDPARA(NKR))/(dfloat(kstps)*dfloat(nmolwat))
      write(n3,9)VD_NUMB,DISTR_V
   endif
end do
AVARA = AVARA/(dfloat(kstps)*dfloat(nmolwat))
write(n0,*)'Average area = ', AVARA
write(n3,*)

!Volume
do NKR =1,NVDELS
   if(NDPVOL(NKR) > 0)then
      VD_NUMB = VDEL*dfloat(NKR)          
      DISTR_V = dfloat(NDPVOL(NKR))/(dfloat(kstps)*dfloat(nmolwat))
      write(n4,9)VD_NUMB,DISTR_V
   endif
end do
AVVOL = AVVOL/(dfloat(kstps)*dfloat(nmolwat))
write(n0,*)'Average volume = ', AVVOL
write(n4,*)

!Asphericity
do NKR =1,NADELS
   if(NDPASP(NKR) > 0)then
      VD_NUMB = ADEL*dfloat(NKR)          
      DISTR_V = dfloat(NDPASP(NKR))/(dfloat(kstps)*dfloat(nmolwat))
!      DISTR_V = dfloat(NDPASP(NKR))/(dfloat(NDPASP_SUM)*ADEL)
      write(n1,9)VD_NUMB,DISTR_V
   endif
end do
AVASP = AVASP/(dfloat(kstps)*dfloat(nmolwat))
write(n0,*)'Average asphericity = ', AVASP
write(n1,*)

!NG End Averages




!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!Average number of waters on each environment

write(n0,*)'Number of molecules = ',nmolwat 
write(*,*)'Number of origins sampled =',kstps
write(n0,*)'Number of origins sampled =',kstps
write(*,*)'Number of origins skipped =',kskip
write(n0,*)'Number of origins skipped =',kskip


close(n0) 
close(n1)
close(n2)
close(n3)
close(n4)


deallocate(xyz)
deallocate(x,y,z)

deallocate(NABLST)
deallocate(NNAB)
deallocate(PX, PY, PZ, PS)
deallocate(TAG, EDGES)
deallocate(RXVER, RYVER, RZVER)
deallocate(IVER, JVER, KVER)


  9   FORMAT(27X,F9.4,5X,F14.5)

    return

END SUBROUTINE pure_wat_Voronoi



!End Voronoi Laboratory



SUBROUTINE pure_wat_diffus_msd(trjfile,npurew,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,msd_delay,difsampw) 
           
! Calculate water's mean square displacement - Pure Water
! Every difsampw time-steps the code samples every water in the system 
    integer,intent(in)                               :: npurew,nbox,nsyst
    character(len=7),intent(in)                      :: inputformat           ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                  ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)               ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: msd_delay             ! msd delay time (time-windows/ps)
    integer,intent(in)                               :: difsampw              ! time freq. to sample waters for msd calculation
   
! Local variables
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
                                                     
    real(kind=8),dimension(:,:,:),allocatable        :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable        :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                                     :: TSTEP,PSECS    
    integer                                    :: NDELS,KINTVL,NDELS1    
    integer,dimension(:,:,:),allocatable       :: NH2OidO    
    integer,dimension(:,:),allocatable         :: ntime
    integer,dimension(:,:),allocatable         :: nH2Oorg                    ! Number of waters used for msd calculation at each origin
    integer,dimension(:),allocatable           :: kv
    integer      :: n0, n1, n2, n3, n4, nf
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer                                    :: norgmax
    integer                                    :: NRSHMAX  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!MSD   
    real(kind=8),dimension(:),allocatable        :: msdx_b,  msdy_b,  msdz_b                         !bulk
    real(kind=8),dimension(:),allocatable        :: msdx_hs, msdy_hs, msdz_hs                        !hsh
    real(kind=8),dimension(:),allocatable        :: msdx_th, msdy_th, msdz_th                        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: msdx_nth,msdy_nth,msdz_nth                       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: msd_mean,msd_mean_x,msd_mean_y,msd_mean_z 
    integer,dimension(:),allocatable             :: nH2Omean
    integer                                      :: kr_start
    integer                                      :: LSF_begin,LSF_end
    integer                                      :: n_LSF_begin,n_LSF_end
    real(kind=8)                                 :: a,b,d,r                    !Y = a + b*X; d = standard deviation of the fit; r = correlation coeficient
    integer                                      :: n                          !number of points used in the fit: range ]n_LSF_begin,n_LSF_end[
    integer                                      :: nroute               !integer to decide which files to open to compute the Einstein diffusion from LSF

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
!define internally range for Linear Squares Fit (ps)

LSF_begin = msd_delay/2                               !Start fit after 10 ps         (move to input if needed)
LSF_end = msd_delay - 1                     !End fit 10 ps before msd end  (move to input if needed)
nroute = 2                                     !Keep always two (2)

n_LSF_begin = IDINT(LSF_begin*1000d0/dt) 
n_LSF_end = IDINT(LSF_end*1000d0/dt) 
    
inquire(file='MSD_PW_out/log_PW_msd.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p MSD_PW_out')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140

open(n0,file='MSD_PW_out/log_PW_msd.dat',status='unknown',action='write')
open(n1,file='MSD_PW_out/msdx.dat',status='unknown',action='write')
open(n2,file='MSD_PW_out/msdy.dat',status='unknown',action='write')
open(n3,file='MSD_PW_out/msdz.dat',status='unknown',action='write')
open(n4,file='MSD_PW_out/msd_mean.dat',status='unknown',action='write')

!time-window for calculation of the mean square displacement
! IDINT - convert to integer - truncate

NDELS = IDINT(msd_delay*1000d0/dt)                  !sets the delay time (in frames)

kr_start = 1                                        !environments loop starter - for pure water single environment 

write(*,*)'input delay-time (ps) =',msd_delay
write(*,*)'delay-time (input frames units) =',NDELS
write(n0,*)'input delay-time (ps) =',msd_delay
write(n0,*)'delay-time (input frames units) =',NDELS

!NRSHMAX = 4                                      !maximum number of environments
NRSHMAX = 1                                       !maximum number of environments - 1 for pure water
if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/difsampw                                  !for memory allocation purposes

!mean square displacement
allocate(nH2Oorg(norgmax,NRSHMAX))
allocate(msdx_b(NDELS),  msdy_b(NDELS),  msdz_b(NDELS))
allocate(msdx_hs(NDELS), msdy_hs(NDELS), msdz_hs(NDELS))
allocate(msdx_th(NDELS) ,msdy_th(NDELS), msdz_th(NDELS))
allocate(msdx_nth(NDELS),msdy_nth(NDELS),msdz_nth(NDELS))
allocate(msd_mean(NDELS,NRSHMAX),msd_mean_x(NDELS,NRSHMAX),msd_mean_y(NDELS,NRSHMAX),msd_mean_z(NDELS,NRSHMAX))
msdx_b=0.0d0;   msdy_b=0.0d0;   msdz_b=0.0d0
msdx_hs=0.0d0;  msdy_hs=0.0d0;  msdz_hs=0.0d0
msdx_th=0.0d0;  msdy_th=0.0d0;  msdz_th=0.0d0
msdx_nth=0.0d0; msdy_nth=0.0d0; msdz_nth=0.0d0
msd_mean=0.0d0;msd_mean_x=0.0d0;msd_mean_y=0.0d0;msd_mean_z=0.0d0


allocate(RXt0(norgmax,nmolwat,NRSHMAX),RYt0(norgmax,nmolwat,NRSHMAX),RZt0(norgmax,nmolwat,NRSHMAX))
allocate(RXt(norgmax,nmolwat,NRSHMAX),RYt(norgmax,nmolwat,NRSHMAX),RZt(norgmax,nmolwat,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(NH2OidO(norgmax,nmolwat,NRSHMAX))
ntime = 1; nH2Oorg = 0
allocate(kv(NRSHMAX))
allocate(nH2Omean(NRSHMAX))
nH2Omean = 0

write(*,*)
write(*,*)'Pure Water mean square displacement/Einstein diffusion calculation'
write(*,*)'Routine pure_wat_diffus_msd'
write(*,*)'water msd delay time-window (ps) =',msd_delay
write(*,*)'Results are printed to MSD_PW_out'
write(*,*)
write(n0,*)
write(n0,*)'Pure Water  mean square displacement/Einstein diffusion calculation'
write(n0,*)'Routine pure_wat_diffus_msd'
write(n0,*)'water msd delay time-window (ps) =',msd_delay
write(n0,*)'Results are printed to MSD_PW_out'
write(n0,*)

! Start msd calculation
write(*,*)
write(*,*)'Starting msd calculation...'
write(*,*)'msd water sampling frequency (steps) = ',difsampw
write(*,*)'total number of time-steps =',nstep
write(*,*)'msd delay-time (steps) =',NDELS
write(*,*)
write(n0,*)
write(n0,*)'Starting msd calculation...'
write(n0,*)'msd water sampling frequency (steps) = ',difsampw
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'msd delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!End Version 9.0     
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 
                                             
!sample water molecules in specific environments - pure water - single environment
   if((j==1.or.mod(j,difsampw)==0).and.(j.le.nstep-NDELS))then               !time origin: sample waters in specific environments
      write(*,*)'time-step',j,' sampling water molecules'
      kv = 0                                                               !waters in each environment 
      ko = ko + 1                                                          !number of origins = number of blocks used to compute the msd
      io = 0
      do i=1,natms,nwatsites                              !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...

!bulk mean square displacement             
         nradsh = 1                                       !bulk environment
         kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
         kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check         write(*,*)j,i,kvi
         nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
         NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
         RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
         RYt0(ko,kvi,nradsh) = y(i)
         RZt0(ko,kvi,nradsh) = z(i)            
!bulk msd(t0) end                     
         
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nH2Oorg(ko,1)
   endif                                                             !end sampling new time-origin      
   
!calculate msd - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nH2Oorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then               !ntime counts the delay times for which the msd was already calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nH2Oorg(kt,kr)                                       !loop over O in environment kr: nH2Oorg can be zero for a given time-origin
                  
               i=NH2OidO(kt,iw,kr)
               RXt(kt,iw,kr) = x(i)
               RYt(kt,iw,kr) = y(i)
               RZt(kt,iw,kr) = z(i)
               nt = ntime(kt,kr)
               if(kr==1)then                                                                            !bulk
                  msdx_b(nt) = msdx_b(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_b(nt) = msdy_b(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_b(nt) = msdz_b(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==2)then                                                                        !hshell
                  msdx_hs(nt) = msdx_hs(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_hs(nt) = msdy_hs(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_hs(nt) = msdz_hs(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==3)then                                                                         !hshell-tetrahedral
                  msdx_th(nt) = msdx_th(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_th(nt) = msdy_th(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_th(nt) = msdz_th(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==4)then                                                                         !hShell-non-tetrahedral
                  msdx_nth(nt) = msdx_nth(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_nth(nt) = msdy_nth(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_nth(nt) = msdz_nth(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean msd over different origins - define new variables "msd_mean"
do kr = 1,NRSHMAX                                                             !loop over environments
   do it = 1,NDELS
!nH2Oorg(kt,kr) - number of O atoms found in origin (block) kt and in the environment kr               
      if(kr==1)then      
         msd_mean_x(it,kr) = msdx_b(it)
         msd_mean_y(it,kr) = msdy_b(it)
         msd_mean_z(it,kr) = msdz_b(it)
         msd_mean(it,kr)   = (msdx_b(it)+msdy_b(it)+msdz_b(it))
      elseif(kr==2)then   
         msd_mean_x(it,kr) = msdx_hs(it)
         msd_mean_y(it,kr) = msdy_hs(it)
         msd_mean_z(it,kr) = msdz_hs(it)
         msd_mean(it,kr)   = (msdx_hs(it)+msdy_hs(it)+msdz_hs(it))
      elseif(kr==3)then
         msd_mean_x(it,kr) = msdx_th(it)
         msd_mean_y(it,kr) = msdy_th(it)
         msd_mean_z(it,kr) = msdz_th(it)
         msd_mean(it,kr)   = (msdx_th(it)+msdy_th(it)+msdz_th(it))
      elseif(kr==4)then
         msd_mean_x(it,kr) = msdx_nth(it)
         msd_mean_y(it,kr) = msdy_nth(it)
         msd_mean_z(it,kr) = msdz_nth(it)
         msd_mean(it,kr)   = (msdx_nth(it)+msdy_nth(it)+msdz_nth(it))
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of O atoms on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nH2Omean(kr) = nH2Omean(kr) + nH2Oorg(kt,kr)
   end do
end do

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko)                      !remove /2 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko)                      !remove /2 
end do   
      
deallocate(msdx_b,  msdy_b,  msdz_b)
deallocate(msdx_hs, msdy_hs, msdz_hs)
deallocate(msdx_th ,msdy_th, msdz_th)

write(n1,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n2,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n3,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n4,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'

!Print mean msd      
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a7)')'#Bulk-x'
   if(kr==1)write(n2,'(a7)')'#Bulk-y'
   if(kr==1)write(n3,'(a7)')'#Bulk-z'
   if(kr==1)write(n4,'(a5)')'#Bulk'
   do it = 1,NDELS
      jt = it -1                                                     !print msd starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      WRITE(n1,19)PSECS,msd_mean_x(it,kr)/dfloat(ko)
      WRITE(n2,19)PSECS,msd_mean_y(it,kr)/dfloat(ko)
      WRITE(n3,19)PSECS,msd_mean_z(it,kr)/dfloat(ko)
      WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(ko)
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
   write(n4,*)
end do

close(n1) 
close(n2) 
close(n3) 
close(n4)

!Calculate the Einstein self-diffusion coefficient - MSD Linear Squares fit
write(n0,*)
write(n0,*)'Einstein self-diffusion coefficient'
write(n0,*)
do kr = kr_start,NRSHMAX
   nf = n4
   if(kr==1)then
      write(n0,'(a5)')'#Bulk'
   elseif(kr==2)then
      write(n0,'(a4)')'#HSh'
   elseif(kr==3)then
      write(n0,'(a11)')'#HSh-tetrah'
   elseif(kr==4)then
      write(n0,'(a15)')'#HSh-non-tetrah'
   endif  

   call msd_LSF_Einstein(nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end,a,b,d,r,n)
   
   if(kr==1)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==2)then
      write(n0,*)
   elseif(kr==3)then
      write(n0,*)
   elseif(kr==4)then
      write(n0,*)
   endif   
end do 

close(n0)

deallocate(nH2Oorg,NH2OidO)
deallocate(ntime)
deallocate(msd_mean,msd_mean_x,msd_mean_y,msd_mean_z)
deallocate(kv,nH2Omean)

deallocate(xyz)
deallocate(x,y,z)

!  9   FORMAT('#',15X,'time(ps)',7X,'msd(Ang**2/ps)'/) 
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
   return

END SUBROUTINE pure_wat_diffus_msd



SUBROUTINE pure_wat_diffus_msd_Env(trjfile,npurew,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,msd_delay,difsampw,R5Dif) 
           
! Calculate water's mean square displacement - Pure Water
! Every difsampw time-steps the code samples every water in the system
!v33 
! Introduced the possibility of calculating the diffusion of water molecules with 3 (or less) 4 (or less) and 5 (or more) neighbors within some cut-off
! Environment 1 - Bulk; Environment 2 - 4LNbE [4 or less neighbors environment]; Environment 3 - 3LNbE; Environment 4 - 5MNbE [5 or more neighbors environment] 

    integer,intent(in)                               :: npurew,nbox,nsyst
    character(len=7),intent(in)                      :: inputformat           ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                  ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)               ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: msd_delay             ! msd delay time (time-windows/ps)
    integer,intent(in)                               :: difsampw              ! time freq. to sample waters for msd calculation
    real(kind=8),intent(in)                          :: R5Dif              ! water neighbors cut-off
   
! Local variables
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
                                                     
    real(kind=8),dimension(:,:,:),allocatable        :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable        :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                                     :: TSTEP,PSECS    
    integer                                    :: NDELS,KINTVL,NDELS1    
    integer,dimension(:,:,:),allocatable       :: NH2OidO    
    integer,dimension(:,:),allocatable         :: ntime
    integer,dimension(:,:),allocatable         :: nH2Oorg                    ! Number of waters used for msd calculation at each origin
    integer,dimension(:),allocatable           :: kv
    integer      :: n0, n1, n2, n3, n4, n5, n6, n7, n8, nf
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt, jw
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer                                    :: norgmax
    integer                                    :: NRSHMAX  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!MSD   
    real(kind=8),dimension(:),allocatable        :: msdx_b,  msdy_b,  msdz_b                         !bulk
    real(kind=8),dimension(:),allocatable        :: msdx_hs, msdy_hs, msdz_hs                        !hsh
    real(kind=8),dimension(:),allocatable        :: msdx_th, msdy_th, msdz_th                        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: msdx_nth,msdy_nth,msdz_nth                       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: msd_mean,msd_mean_x,msd_mean_y,msd_mean_z 
    integer,dimension(:),allocatable             :: nH2Omean
    integer                                      :: kr_start
    integer                                      :: LSF_begin,LSF_end
    integer                                      :: n_LSF_begin,n_LSF_end
    real(kind=8)                                 :: a,b,d,r                    !Y = a + b*X; d = standard deviation of the fit; r = correlation coeficient
    integer                                      :: n                          !number of points used in the fit: range ]n_LSF_begin,n_LSF_end[
    integer                                      :: nroute               !integer to decide which files to open to compute the Einstein diffusion from LSF

!v33 VARIABLES

    logical                                      :: LARRHEN              ! .TRUE. calculate self-diffusion of water ensembles with 4 or 5 water neighbors within cut-off (R5Dif)
    real(kind=8)                                 :: RADMAX,RADMIN
    integer                                      :: KO5,KO4
    real(kind=8)                                 :: dx,dy,dz,dr,dr2,dr4,dr5
    INTEGER                                      :: NTD5O(5)
    
!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj

!v33 VARIABLES
LARRHEN = .true.
!R5Dif = 3.7        !cut-off in Angstroms  [defined in the input file]

!End v33



!define internally range for Linear Squares Fit (ps)

LSF_begin = msd_delay/2                               !Start fit after 10 ps         (move to input if needed)
LSF_end = msd_delay - 1                               !End fit 10 ps before msd end  (move to input if needed)
nroute = 2                                            !Keep always two (2)
if(LARRHEN)nroute = 4

n_LSF_begin = IDINT(LSF_begin*1000d0/dt) 
n_LSF_end = IDINT(LSF_end*1000d0/dt) 
    
inquire(file='MSD_PW_out/log_PW_msd.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p MSD_PW_out')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140

n5 = 150
n6 = 160
n7 = 170
n8 = 180

open(n0,file='MSD_PW_out/log_PW_msd.dat',status='unknown',action='write')
open(n1,file='MSD_PW_out/msdx.dat',status='unknown',action='write')
open(n2,file='MSD_PW_out/msdy.dat',status='unknown',action='write')
open(n3,file='MSD_PW_out/msdz.dat',status='unknown',action='write')
open(n4,file='MSD_PW_out/msd_mean.dat',status='unknown',action='write')

open(n5,file='MSD_PW_out/msd_bulk.dat',status='unknown',action='write')
open(n6,file='MSD_PW_out/msd_4LNbE.dat',status='unknown',action='write')
open(n7,file='MSD_PW_out/msd_3LNbE.dat',status='unknown',action='write')
open(n8,file='MSD_PW_out/msd_5MNbE.dat',status='unknown',action='write')


!time-window for calculation of the mean square displacement
! IDINT - convert to integer - truncate

NDELS = IDINT(msd_delay*1000d0/dt)                  !sets the delay time (in frames)

kr_start = 1                                        !environments loop starter - for pure water single environment 

write(*,*)'input delay-time (ps) =',msd_delay
write(*,*)'delay-time (input frames units) =',NDELS
write(n0,*)'input delay-time (ps) =',msd_delay
write(n0,*)'delay-time (input frames units) =',NDELS

if(LARRHEN) NRSHMAX = 4                                      !maximum number of environments - 4 for pure water - analyse environments
if(.not.LARRHEN)NRSHMAX = 1                                  !maximum number of environments - 1 for pure water - bulk


if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/difsampw                                  !for memory allocation purposes

!mean square displacement
allocate(nH2Oorg(norgmax,NRSHMAX))
allocate(msdx_b(NDELS),  msdy_b(NDELS),  msdz_b(NDELS))
allocate(msdx_hs(NDELS), msdy_hs(NDELS), msdz_hs(NDELS))
allocate(msdx_th(NDELS) ,msdy_th(NDELS), msdz_th(NDELS))
allocate(msdx_nth(NDELS),msdy_nth(NDELS),msdz_nth(NDELS))
allocate(msd_mean(NDELS,NRSHMAX),msd_mean_x(NDELS,NRSHMAX),msd_mean_y(NDELS,NRSHMAX),msd_mean_z(NDELS,NRSHMAX))
msdx_b=0.0d0;   msdy_b=0.0d0;   msdz_b=0.0d0
msdx_hs=0.0d0;  msdy_hs=0.0d0;  msdz_hs=0.0d0
msdx_th=0.0d0;  msdy_th=0.0d0;  msdz_th=0.0d0
msdx_nth=0.0d0; msdy_nth=0.0d0; msdz_nth=0.0d0
msd_mean=0.0d0;msd_mean_x=0.0d0;msd_mean_y=0.0d0;msd_mean_z=0.0d0


allocate(RXt0(norgmax,nmolwat,NRSHMAX),RYt0(norgmax,nmolwat,NRSHMAX),RZt0(norgmax,nmolwat,NRSHMAX))
allocate(RXt(norgmax,nmolwat,NRSHMAX),RYt(norgmax,nmolwat,NRSHMAX),RZt(norgmax,nmolwat,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(NH2OidO(norgmax,nmolwat,NRSHMAX))
ntime = 1; nH2Oorg = 0
allocate(kv(NRSHMAX))
allocate(nH2Omean(NRSHMAX))
nH2Omean = 0

write(*,*)
write(*,*)'Pure Water mean square displacement/Einstein diffusion - Environments version'
write(*,*)'Routine pure_wat_diffus_msd_Env'
if(LARRHEN)then
   write(*,*)'Self-diffusion will be calculated for 4 environments'
   write(*,*)'interstitial cut-off = ',R5Dif,' Ang'
else
   write(*,*)'Self-diffusion of water - no environment deconvolution'
endif
write(*,*)'water msd delay time-window (ps) =',msd_delay
write(*,*)'Results are printed to MSD_PW_out'
write(*,*)
write(n0,*)
write(n0,*)'Pure Water mean square displacement/Einstein diffusion - Environments version'
write(n0,*)'Routine pure_wat_diffus_msd_Env'
if(LARRHEN)then
   write(n0,*)'Self-diffusion will be calculated for 4 environments'
   write(n0,*)'interstitial cut-off = ',R5Dif,' Ang'
else
   write(n0,*)'Self-diffusion of water - no environment deconvolution'
endif
write(n0,*)'water msd delay time-window (ps) =',msd_delay
write(n0,*)'Results are printed to MSD_PW_out'
write(n0,*)

! Start msd calculation
write(*,*)
write(*,*)'Starting msd calculation...'
write(*,*)'msd water sampling frequency (steps) = ',difsampw
write(*,*)'total number of time-steps =',nstep
write(*,*)'msd delay-time (steps) =',NDELS
write(*,*)
write(n0,*)
write(n0,*)'Starting msd calculation...'
write(n0,*)'msd water sampling frequency (steps) = ',difsampw
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'msd delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!End Version 9.0     
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 
                                             
!sample water molecules in specific environments - pure water - single environment
   if((j==1.or.mod(j,difsampw)==0).and.(j.le.nstep-NDELS))then               !time origin: sample waters in specific environments
      write(*,*)'time-step',j,' sampling water molecules'
      kv = 0                                                               !waters in each environment 
      ko = ko + 1                                                          !number of origins = number of blocks used to compute the msd
      io = 0
      do i=1,natms,nwatsites                              !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...

!bulk mean square displacement             
         nradsh = 1                                       !bulk environment
         kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
         kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check         write(*,*)j,i,kvi
         nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
         NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
         RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
         RYt0(ko,kvi,nradsh) = y(i)
         RZt0(ko,kvi,nradsh) = z(i)            
!bulk msd(t0) end

!version 33 - Move the next part to a a separate routine in water_solv.f90 when validated

         if(LARRHEN)then                   !multiple environments
! Find the 5 closest oxygen atoms to the oxygen i
            RADMAX  = 25.0D0
            RADMIN  =  0.0D0
            JT      = 1
            DO WHILE(JT.LE.5)
               do jw=1,natms,nwatsites                                       !Loop over water oxygens
                  IF(jw.NE.i)THEN
            ! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND OXYGEN ATOMS jw     
                     dx = x(i) - x(jw)             
                     dy = y(i) - y(jw)
                     dz = z(i) - z(jw)
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                     dr2 = dx**2 + dy**2 + dz**2
                     dr  = dsqrt(dr2)        
            !FIND THE NEAREST 5 OXYGEN ATOMS jw TO THE OXYGEN ATOM i
                     IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
                        RADMAX = dr
            !NTD5O STORES THE O ATOMS NEIGHBORS OF I
            !JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR ; JT = 3 THIRD NEIGHBOR ; JT = 4 FOURTH NEIGHBOR ; JT = 5 FIFTH NEIGHBOR
                        NTD5O(JT) = jw
                     ENDIF
                  ENDIF
               end do
               JT = JT + 1
               RADMIN = RADMAX
               RADMAX = 25.0D0
            END DO    
            
!distance to the fifth neighbor         
            KO5 = NTD5O(5)
            dx = x(i) - x(KO5)             
            dy = y(i) - y(KO5)
            dz = z(i) - z(KO5)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2) 
            dr5 = dr
            
!distance to the fourth neighbor         
            KO4 = NTD5O(4)
            dx = x(i) - x(KO4)             
            dy = y(i) - y(KO4)
            dz = z(i) - z(KO4)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2) 
            dr4 = dr            
            
            if(dr5 > R5Dif)then                               !fifth neighbor beyond cut-off
!four or less neighbors mean square displacement             
               nradsh = 2                                       !four or less neighbors environment
               kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
               kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check               write(*,*)j,i,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)
               
            elseif(dr5 > R5Dif .and. dr4 > R5Dif)then       !fifth and fourth neighbors beyond cut-off
!three or less neighbors mean square displacement             
               nradsh = 3                                       !three or less neighbors environment
               kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
               kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check               write(*,*)j,i,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)
            elseif(dr5 < R5Dif)then                           !fifth neighbor within cut-off
!five or more neighbors mean square displacement             
               nradsh = 4                                       !five or more neighbors environment
               kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
               kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check               write(*,*)j,i,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)  
            endif
!environments msd(t0) end            
                        
         endif
         
!end version 33


      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nH2Oorg(ko,1)
   endif                                                             !end sampling new time-origin      
   
!calculate msd - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nH2Oorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then               !ntime counts the delay times for which the msd was already calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nH2Oorg(kt,kr)                                       !loop over O in environment kr: nH2Oorg can be zero for a given time-origin
                  
               i=NH2OidO(kt,iw,kr)
               RXt(kt,iw,kr) = x(i)
               RYt(kt,iw,kr) = y(i)
               RZt(kt,iw,kr) = z(i)
               nt = ntime(kt,kr)
               if(kr==1)then                                                                            !bulk
                  msdx_b(nt) = msdx_b(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_b(nt) = msdy_b(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_b(nt) = msdz_b(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==2)then                                                                        !self-HShell - 4 or less neighbors [4LNbE]
                  msdx_hs(nt) = msdx_hs(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_hs(nt) = msdy_hs(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_hs(nt) = msdz_hs(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==3)then                                                                         !self-HShell - 3 or less neighbors [3LNbE]
                  msdx_th(nt) = msdx_th(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_th(nt) = msdy_th(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_th(nt) = msdz_th(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==4)then                                                                         !self-HShell - 5 or more neighbors [5MNbE]
                  msdx_nth(nt) = msdx_nth(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_nth(nt) = msdy_nth(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_nth(nt) = msdz_nth(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean msd over different origins - define new variables "msd_mean"
do kr = 1,NRSHMAX                                                             !loop over environments
   do it = 1,NDELS
!nH2Oorg(kt,kr) - number of O atoms found in origin (block) kt and in the environment kr               
      if(kr==1)then      
         msd_mean_x(it,kr) = msdx_b(it)
         msd_mean_y(it,kr) = msdy_b(it)
         msd_mean_z(it,kr) = msdz_b(it)
         msd_mean(it,kr)   = (msdx_b(it)+msdy_b(it)+msdz_b(it))
      elseif(kr==2)then   
         msd_mean_x(it,kr) = msdx_hs(it)
         msd_mean_y(it,kr) = msdy_hs(it)
         msd_mean_z(it,kr) = msdz_hs(it)
         msd_mean(it,kr)   = (msdx_hs(it)+msdy_hs(it)+msdz_hs(it))
      elseif(kr==3)then
         msd_mean_x(it,kr) = msdx_th(it)
         msd_mean_y(it,kr) = msdy_th(it)
         msd_mean_z(it,kr) = msdz_th(it)
         msd_mean(it,kr)   = (msdx_th(it)+msdy_th(it)+msdz_th(it))
      elseif(kr==4)then
         msd_mean_x(it,kr) = msdx_nth(it)
         msd_mean_y(it,kr) = msdy_nth(it)
         msd_mean_z(it,kr) = msdz_nth(it)
         msd_mean(it,kr)   = (msdx_nth(it)+msdy_nth(it)+msdz_nth(it))
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of O atoms on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nH2Omean(kr) = nH2Omean(kr) + nH2Oorg(kt,kr)
   end do
end do

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko)                      !remove /2 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko)                      !remove /2 
end do   
      
deallocate(msdx_b,  msdy_b,  msdz_b)
deallocate(msdx_hs, msdy_hs, msdz_hs)
deallocate(msdx_th ,msdy_th, msdz_th)

write(n1,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n2,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n3,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n4,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n6,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n7,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n8,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'

!Print mean msd      
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a7)')'#Bulk-x'
   if(kr==1)write(n2,'(a7)')'#Bulk-y'
   if(kr==1)write(n3,'(a7)')'#Bulk-z'
   if(kr==2)write(n1,'(a7)')'#4LNbE-x'
   if(kr==2)write(n2,'(a7)')'#4LNbE-y'
   if(kr==2)write(n3,'(a7)')'#4LNbE-z'
   if(kr==3)write(n1,'(a7)')'#3LNbE-x'
   if(kr==3)write(n2,'(a7)')'#3LNbE-y'
   if(kr==3)write(n3,'(a7)')'#3LNbE-z'
   if(kr==4)write(n1,'(a7)')'#5MNbE-x'
   if(kr==4)write(n2,'(a7)')'#5MNbE-y'
   if(kr==4)write(n3,'(a7)')'#5MNbE-z'
   
   if(kr==1)write(n4,'(a5)')'#Bulk'
   if(kr==2)write(n4,'(a6)')'#4LNbE'
   if(kr==3)write(n4,'(a6)')'#3LNbE'
   if(kr==4)write(n4,'(a6)')'#5MNbE'
   
   if(kr==1)write(n5,'(a5)')'#Bulk'
   if(kr==2)write(n6,'(a6)')'#4LNbE'
   if(kr==3)write(n7,'(a6)')'#3LNbE'
   if(kr==4)write(n8,'(a6)')'#5MNbE'
!   
   do it = 1,NDELS
      jt = it -1                                                     !print msd starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      WRITE(n1,19)PSECS,msd_mean_x(it,kr)/dfloat(ko)
      WRITE(n2,19)PSECS,msd_mean_y(it,kr)/dfloat(ko)
      WRITE(n3,19)PSECS,msd_mean_z(it,kr)/dfloat(ko)
      WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      
      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)          !bulk
      if(kr==2)WRITE(n6,19)PSECS,msd_mean(it,kr)/dfloat(ko)          !4LNbE
      if(kr==3)WRITE(n7,19)PSECS,msd_mean(it,kr)/dfloat(ko)          !3LNbE
      if(kr==4)WRITE(n8,19)PSECS,msd_mean(it,kr)/dfloat(ko)          !5MNbE
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
   write(n4,*)
end do

close(n1) 
close(n2) 
close(n3) 
close(n4)
close(n5)
close(n6)
close(n7)
close(n8)

!Calculate the Einstein self-diffusion coefficient - MSD Linear Squares fit
write(n0,*)
write(n0,*)'Einstein self-diffusion coefficient'
write(n0,*)
do kr = kr_start,NRSHMAX
   if(kr==1)then
      nf = n5
      write(n0,'(a5)')'#Bulk'
      write(n0,*)
   elseif(kr==2)then
      nf = n6
      write(n0,'(a6)')'#4LNbE'
      write(n0,*)
   elseif(kr==3)then
      nf = n7
      write(n0,'(a6)')'#3LNbE'
      write(n0,*)
   elseif(kr==4)then
      nf = n8
      write(n0,'(a6)')'#5MNbE'
      write(n0,*)
   endif  

   call msd_LSF_Einstein(nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end,a,b,d,r,n)
   
   if(kr==1)then
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==2)then
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==3)then
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==4)then
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(n0,*)
      write(*,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   endif   
end do 

close(n0)

deallocate(nH2Oorg,NH2OidO)
deallocate(ntime)
deallocate(msd_mean,msd_mean_x,msd_mean_y,msd_mean_z)
deallocate(kv,nH2Omean)

deallocate(xyz)
deallocate(x,y,z)

!  9   FORMAT('#',15X,'time(ps)',7X,'msd(Ang**2/ps)'/) 
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
   return

END SUBROUTINE pure_wat_diffus_msd_Env



SUBROUTINE pure_wat_LDL_HDL_msd(trjfile,npurew,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,msd_delay,difsampw,Nenvrm) 
! Lars Pettersson and Gaia Camisasca collaboration - June 2018 - coded in version 24
! The program reads a list from an input file with information on LDL and HDL waters           
! Calculate water's mean square displacement - pure Water
! Every difsampw time-steps the code samples water molecules in the system that are either LDL, HDL or neither
    integer,intent(in)                               :: npurew,nbox,nsyst,Nenvrm
    character(len=7),intent(in)                      :: inputformat           ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt
    integer,intent(in)                               :: nens                  ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)               ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname
    integer,intent(in)                               :: msd_delay             ! msd delay time (time-windows/ps)
    integer,intent(in)                               :: difsampw              ! time freq. to sample waters for msd calculation
   
! Local variables
!NG    real(kind=8)                                     :: x(natms),y(natms),z(natms) 
!    real(kind=4)                                     :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable            :: x,y,z
!NG    real(kind=4)                                     :: xx(natms),yy(natms),zz(natms)
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                    !atomic index (solute and solvent index)
                                                     
    real(kind=8),dimension(:,:,:),allocatable        :: RXt0,RYt0,RZt0            !OH vectors at t0
    real(kind=8),dimension(:,:,:),allocatable        :: RXt,RYt,RZt               !OH vectors at t
    real(kind=8)                                     :: TSTEP,PSECS    
    integer                                    :: NDELS,KINTVL,NDELS1    
    integer,dimension(:,:,:),allocatable       :: NH2OidO    
    integer,dimension(:,:),allocatable         :: ntime
    integer,dimension(:,:),allocatable         :: nH2Oorg                    ! Number of waters used for msd calculation at each origin
    integer,dimension(:),allocatable           :: kv
    integer      :: n0, n1, n2, n3, n4, nf
    integer      :: n5, n6, n7, n8
    integer      :: n20
    integer      :: i, j, k, L
    integer      :: kw, iw, kt, kr, kvi, ki
    integer      :: it, iH, ihat, nt, ko, jt
    integer      :: io
    integer,parameter                          :: INTMAX = 2147483647
    integer                                    :: norgmax
    integer                                    :: NRSHMAX  
    integer                                    :: nradsh    
    integer                                    :: natcheck

    logical filex
    
!MSD   
    real(kind=8),dimension(:),allocatable        :: msdx_b,  msdy_b,  msdz_b                         !bulk
    real(kind=8),dimension(:),allocatable        :: msdx_hs, msdy_hs, msdz_hs                        !hsh
    real(kind=8),dimension(:),allocatable        :: msdx_th, msdy_th, msdz_th                        !hsh-tetrahedral
    real(kind=8),dimension(:),allocatable        :: msdx_nth,msdy_nth,msdz_nth                       !hsh-non-tetrahedral
    real(kind=8),dimension(:,:),allocatable      :: msd_mean,msd_mean_x,msd_mean_y,msd_mean_z 
    integer,dimension(:),allocatable             :: nH2Omean
    integer                                      :: kr_start
    integer                                      :: LSF_begin,LSF_end
    integer                                      :: n_LSF_begin,n_LSF_end
    real(kind=8)                                 :: a,b,d,r                    !Y = a + b*X; d = standard deviation of the fit; r = correlation coeficient
    integer                                      :: n                          !number of points used in the fit: range ]n_LSF_begin,n_LSF_end[
    integer                                      :: nroute               !integer to decide which files to open to compute the Einstein diffusion from LSF
!Version 24    
    integer                                      :: nw_list
    integer,dimension(:),allocatable             :: list_s    
    

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj    
    
!define internally range for Linear Squares Fit (ps)

LSF_begin = msd_delay/2                               !Start fit after 10 ps         (move to input if needed)
LSF_end = msd_delay - 1                     !End fit 10 ps before msd end  (move to input if needed)
nroute = 3                                     !Keep always three (3)

n_LSF_begin = IDINT(LSF_begin*1000d0/dt) 
n_LSF_end = IDINT(LSF_end*1000d0/dt) 
    
inquire(file='MSD_LDL_HDL/log_PW_msd.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p MSD_LDL_HDL')

n0 = 100
n1 = 110
n2 = 120
n3 = 130
n4 = 140

n5 = 150
n6 = 160
n7 = 170
n8 = 180

n20 = 300

open(n0,file='MSD_LDL_HDL/log_PW_msd.dat',status='unknown',action='write')
open(n1,file='MSD_LDL_HDL/msdx.dat',status='unknown',action='write')
open(n2,file='MSD_LDL_HDL/msdy.dat',status='unknown',action='write')
open(n3,file='MSD_LDL_HDL/msdz.dat',status='unknown',action='write')
open(n4,file='MSD_LDL_HDL/msd_mean.dat',status='unknown',action='write')

open(n5,file='MSD_LDL_HDL/Bulk.dat',status='unknown',action='write')
open(n6,file='MSD_LDL_HDL/LDL.dat',status='unknown',action='write')
open(n7,file='MSD_LDL_HDL/HDL.dat',status='unknown',action='write')
open(n8,file='MSD_LDL_HDL/other.dat',status='unknown',action='write')

!Version 24 - read LDL/HDL water molecules label from input file
!Lars Pettersson - list
! 1 is LDL 
! 2 is HDL
! 3 skip
if(Nenvrm>1)open(n20,file='List_W.dat',status='old',action='read')

!time-window for calculation of the mean square displacement
! IDINT - convert to integer - truncate

NDELS = IDINT(msd_delay*1000d0/dt)                  !sets the delay time (in frames)

kr_start = 1                                        !environments loop starter - for pure water single environment 

write(*,*)'input delay-time (ps) =',msd_delay
write(*,*)'delay-time (input frames units) =',NDELS
write(n0,*)'input delay-time (ps) =',msd_delay
write(n0,*)'delay-time (input frames units) =',NDELS

if(Nenvrm>1)then
   write(*,*)
   write(*,*)'Environments: 1 (Bulk), 2(LDL), 3(HDL), 4(other)'
   write(*,*)
   write(n0,*)
   write(n0,*)'Environments: 1 (Bulk), 2(LDL), 3(HDL), 4(other)'
   write(n0,*)
endif 

NRSHMAX = Nenvrm                                       !maximum number of environments - 4 e.g. LDL, HDL, rest, all
!NRSHMAX = 1                                           !maximum number of environments - 1 for pure water

if(NDELS.ge.nstep)NDELS = nstep - 1
NDELS1 = NDELS - 1
!TSTEP - time-step in ps
TSTEP  = dt*1.0d-3
!KINTVL - interval between frames
KINTVL =  1
norgmax = 1 + nstep/difsampw                                  !for memory allocation purposes

!mean square displacement
allocate(nH2Oorg(norgmax,NRSHMAX))
allocate(msdx_b(NDELS),  msdy_b(NDELS),  msdz_b(NDELS))
allocate(msdx_hs(NDELS), msdy_hs(NDELS), msdz_hs(NDELS))
allocate(msdx_th(NDELS) ,msdy_th(NDELS), msdz_th(NDELS))
allocate(msdx_nth(NDELS),msdy_nth(NDELS),msdz_nth(NDELS))
allocate(msd_mean(NDELS,NRSHMAX),msd_mean_x(NDELS,NRSHMAX),msd_mean_y(NDELS,NRSHMAX),msd_mean_z(NDELS,NRSHMAX))
msdx_b=0.0d0;   msdy_b=0.0d0;   msdz_b=0.0d0
msdx_hs=0.0d0;  msdy_hs=0.0d0;  msdz_hs=0.0d0
msdx_th=0.0d0;  msdy_th=0.0d0;  msdz_th=0.0d0
msdx_nth=0.0d0; msdy_nth=0.0d0; msdz_nth=0.0d0
msd_mean=0.0d0;msd_mean_x=0.0d0;msd_mean_y=0.0d0;msd_mean_z=0.0d0


allocate(RXt0(norgmax,nmolwat,NRSHMAX),RYt0(norgmax,nmolwat,NRSHMAX),RZt0(norgmax,nmolwat,NRSHMAX))
allocate(RXt(norgmax,nmolwat,NRSHMAX),RYt(norgmax,nmolwat,NRSHMAX),RZt(norgmax,nmolwat,NRSHMAX))
allocate(ntime(norgmax,NRSHMAX))
allocate(NH2OidO(norgmax,nmolwat,NRSHMAX))
ntime = 1; nH2Oorg = 0
allocate(kv(NRSHMAX))
allocate(nH2Omean(NRSHMAX))
nH2Omean = 0

!Version 24 - List HDL/LDL
allocate(list_s(nmolwat))

write(*,*)
write(*,*)'Pure Water mean square displacement/Einstein diffusion calculation'
write(*,*)'Routine pure_wat_LDL_HDL_msd'
write(*,*)'water msd delay time-window (ps) =',msd_delay
write(*,*)'Results are printed to MSD_LDL_HDL'
write(*,*)
write(n0,*)
write(n0,*)'Pure Water  mean square displacement/Einstein diffusion calculation'
write(n0,*)'Routine pure_wat_LDL_HDL_msd'
write(n0,*)'water msd delay time-window (ps) =',msd_delay
write(n0,*)'Results are printed to MSD_LDL_HDL'
write(n0,*)

! Start msd calculation
write(*,*)
write(*,*)'Starting msd calculation...'
write(*,*)'msd water sampling frequency (steps) = ',difsampw
write(*,*)'total number of time-steps =',nstep
write(*,*)'msd delay-time (steps) =',NDELS
write(*,*)
write(n0,*)
write(n0,*)'Starting msd calculation...'
write(n0,*)'msd water sampling frequency (steps) = ',difsampw
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'msd delay-time (steps) =',NDELS
write(n0,*)

ko = 0                   !sampling origins counter - counts the number of origins used to sample waters
do j = 1,nstep
   if(mod(j,5000)==0)write(*,*)'starting time-step ',j
!/Trajectory Reading Formats/ 
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2) = cell(1)
!         cell(3) = cell(1)
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
      endif   
!check   write(*,*)cell(1)
!read water coordinates
      read(npurew,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif
!Version 9.0      
      if(inputformat.eq.'BOMD   ')read(npurew,*)                                                    !comment line .xyz version 9.0
      do i = 1,natms
          if(inputformat.eq.'TINKER ')read(npurew,*)natmid(i),atomname(i),x(i),y(i),z(i)
          if(inputformat.eq.'BOMD   ')read(npurew,*)atomname(i),x(i),y(i),z(i)
      end do
!End Version 9.0     
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
    status = 1   ! status=1 on entry means read
    call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
    
    cell(1)=dble(box(1))
    cell(2)=dble(box(2))
    cell(3)=dble(box(3))
    kp=1
    do i = 1,npart*3,3
       x(kp)=xyz(i)
       y(kp)=xyz(i+1)
       z(kp)=xyz(i+2)
       kp = kp + 1
    end do   

!NG-version8.0      read(npurew)cell(1),cell(2),cell(3)
!read solute and water coordinates     
!NG-version8.0       do i = 1,natms
!NG-version8.0         read(npurew)x(i),y(i),z(i)
!NG         read(npurew)xx(i),yy(i),zz(i)
!NG         x(i) = dble(xx(i))
!NG         y(i) = dble(yy(i))
!NG         z(i) = dble(zz(i))
!NG-version8.0       end do      
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj 
                                             
!sample water molecules in specific environments - pure water - single environment
   if((j==1.or.mod(j,difsampw)==0).and.(j.le.nstep-NDELS))then               !time origin: sample waters in specific environments
      write(*,*)'time-step',j,' sampling water molecules'
      
!Version 24 - LDL/HDL water list           1(LDL) 2(HDL) 3(NOT LDL nor HDL)      
      do i=1,nmolwat
         read(n20,*)nw_list,list_s(i)
      end do           
      
      kv = 0                                                               !waters in each environment 
      ko = ko + 1                                                          !number of origins = number of blocks used to compute the msd
      io = 0                                              !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...
      do i=1,natms,nwatsites                              !Loop over water oxygens
         io = io + 1                                      !Oxygen atoms number - io = 1,2,3,4,5,...; i = 1,4,7,10,13,...

!bulk mean square displacement - all waters            
         nradsh = 1                                       !bulk environment
         kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
         kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check         write(*,*)j,i,kvi
         nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
         NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
         RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
         RYt0(ko,kvi,nradsh) = y(i)
         RZt0(ko,kvi,nradsh) = z(i)            
!bulk msd(t0) end    


!LDL/HDL and NOT LDL nor HDL msd
         if(NRSHMAX > 1)then
            if(list_s(io)==1)then
!LDL msd(t0)           
               nradsh = 2                                       !LDL environment
               kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
               kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check         write(*,*)j,i,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)            
!LDL msd(t0) end                  
!HDL msd(t0)
            elseif(list_s(io)==2)then
               nradsh = 3                                       !HDL environment
               kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
               kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check         write(*,*)j,i,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)            
!HDL msd(t0) end                  
!NOT LDL/HDL msd(t0)
            elseif(list_s(io)==3)then
               nradsh = 4                                       !HDL environment
               kv(nradsh) = kv(nradsh) + 1                      !number of oxygens (index)
               kvi = kv(nradsh)                                 !O number 1,2,3,4,5,...
!check         write(*,*)j,i,kvi
               nH2Oorg(ko,nradsh) = nH2Oorg(ko,nradsh) + 1      !number of waters in environment nradsh at origin ko
               NH2OidO(ko,kvi,nradsh) = i                       !Oxygen id
               RXt0(ko,kvi,nradsh) = x(i)                       !origin = ko, O number = kvi, environment type = nradsh
               RYt0(ko,kvi,nradsh) = y(i)
               RZt0(ko,kvi,nradsh) = z(i)            
!NOT LDL/HDL msd(t0) end                  
            endif
         endif   
         
      end do                                                         !end loop over water oxygens
!check      write(*,*)ko,nH2Oorg(ko,1)

   endif                                                             !end sampling new time-origin      
   
!calculate msd - for each origin average over different waters found on each environment
!check   write(*,*)'time-step = ',j,'number of origins = ',ko
   do kt = 1,ko                                                           !loop over time origins
      do kr = 1,NRSHMAX                                                   !loop over environments
         if(nH2Oorg(kt,kr)>0.and.ntime(kt,kr).le.NDELS)then               !ntime counts the delay times for which the msd was already calculated
!check         write(*,*)ntime(kt,kr)
            do iw = 1,nH2Oorg(kt,kr)                                       !loop over O in environment kr: nH2Oorg can be zero for a given time-origin
                  
               i=NH2OidO(kt,iw,kr)
               RXt(kt,iw,kr) = x(i)
               RYt(kt,iw,kr) = y(i)
               RZt(kt,iw,kr) = z(i)
               nt = ntime(kt,kr)
               if(kr==1)then                                                                            !bulk
                  msdx_b(nt) = msdx_b(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_b(nt) = msdy_b(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_b(nt) = msdz_b(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==2)then                                                                        !LDL
                  msdx_hs(nt) = msdx_hs(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_hs(nt) = msdy_hs(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_hs(nt) = msdz_hs(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==3)then                                                                         !HDL
                  msdx_th(nt) = msdx_th(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_th(nt) = msdy_th(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_th(nt) = msdz_th(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               elseif(kr==4)then                                                                         !NOT LDL nor HDL
                  msdx_nth(nt) = msdx_nth(nt) + (RXt0(kt,iw,kr)-RXt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdy_nth(nt) = msdy_nth(nt) + (RYt0(kt,iw,kr)-RYt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
                  msdz_nth(nt) = msdz_nth(nt) + (RZt0(kt,iw,kr)-RZt(kt,iw,kr))**2/dfloat(nH2Oorg(kt,kr))
               endif
            end do   
            ntime(kt,kr) = ntime(kt,kr) + 1                               !delay time counter
         endif
      end do                                                              !end loop over environments
   end do                                                                 !end loop over time-origins
   
!xtc2fortran_trj
   if(inputformat.eq.'GROMACS'.and.status.eq.0) then
      write(*,*)'error on reading trajectory file - step',i
      stop
   endif
!end xtc2fortran_trj    
 
end do                                                                              !end time-step main loop

!xtc2fortran_trj
IF(inputformat.eq.'GROMACS')THEN
   call f77_molfile_finish
   write(*,*)'finished reading xtc trj...'
ENDIF   
!end xtc2fortran_trj                

!calculate mean msd over different origins - define new variables "msd_mean"
do kr = 1,NRSHMAX                                                             !loop over environments
   do it = 1,NDELS
!nH2Oorg(kt,kr) - number of O atoms found in origin (block) kt and in the environment kr               
      if(kr==1)then      
         msd_mean_x(it,kr) = msdx_b(it)
         msd_mean_y(it,kr) = msdy_b(it)
         msd_mean_z(it,kr) = msdz_b(it)
         msd_mean(it,kr)   = (msdx_b(it)+msdy_b(it)+msdz_b(it))
      elseif(kr==2)then   
         msd_mean_x(it,kr) = msdx_hs(it)
         msd_mean_y(it,kr) = msdy_hs(it)
         msd_mean_z(it,kr) = msdz_hs(it)
         msd_mean(it,kr)   = (msdx_hs(it)+msdy_hs(it)+msdz_hs(it))
      elseif(kr==3)then
         msd_mean_x(it,kr) = msdx_th(it)
         msd_mean_y(it,kr) = msdy_th(it)
         msd_mean_z(it,kr) = msdz_th(it)
         msd_mean(it,kr)   = (msdx_th(it)+msdy_th(it)+msdz_th(it))
      elseif(kr==4)then
         msd_mean_x(it,kr) = msdx_nth(it)
         msd_mean_y(it,kr) = msdy_nth(it)
         msd_mean_z(it,kr) = msdz_nth(it)
         msd_mean(it,kr)   = (msdx_nth(it)+msdy_nth(it)+msdz_nth(it))
      endif   
   end do
end do

deallocate(RXt0,RYt0,RZt0,RXt,RYt,RZt)

!Average number of O atoms on each environment, kr, over the number of origins, ko 
do kr = 1,NRSHMAX
   do kt = 1, ko
      nH2Omean(kr) = nH2Omean(kr) + nH2Oorg(kt,kr)
   end do
end do

do kr = kr_start,NRSHMAX 
   write(*,*)      
   write(*,'(a15,1x,i4)')'#Environment = ',kr 
   write(*,'(a21,1x,i7)')'#Number of origins = ',ko
   write(*,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko)                      !remove /2 
   write(n0,*)
   write(n0,'(a15,1x,i4)')'#Environment = ',kr 
   write(n0,'(a21,1x,i7)')'#Number of origins = ',ko
   write(n0,'(a25,1x,F14.2)')'#Mean number of waters = ',dfloat(nH2Omean(kr))/dfloat(ko)                      !remove /2 
end do   
      
deallocate(msdx_b,  msdy_b,  msdz_b)
deallocate(msdx_hs, msdy_hs, msdz_hs)
deallocate(msdx_th ,msdy_th, msdz_th)

write(n1,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n2,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n3,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n4,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'

write(n5,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n6,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n7,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'
write(n8,'(a1,15x,a8,7x,a14)')'#','time(ps)','msd(Ang**2/ps)'

!Print mean msd      
do kr = kr_start,NRSHMAX                                                      !loop over environments
   if(kr==1)write(n1,'(a7)')'#Bulk-x'
   if(kr==1)write(n2,'(a7)')'#Bulk-y'
   if(kr==1)write(n3,'(a7)')'#Bulk-z'
   if(kr==1)write(n4,'(a5)')'#Bulk'
   if(kr==1)write(n5,'(a5)')'#Bulk'
   if(kr==2)write(n1,'(a6)')'#LDL-x'
   if(kr==2)write(n2,'(a6)')'#LDL-y'
   if(kr==2)write(n3,'(a6)')'#LDL-z'
   if(kr==2)write(n4,'(a4)')'#LDL'
   if(kr==2)write(n6,'(a4)')'#LDL'
   if(kr==3)write(n1,'(a6)')'#HDL-x'
   if(kr==3)write(n2,'(a6)')'#HDL-y'
   if(kr==3)write(n3,'(a6)')'#HDL-z'
   if(kr==3)write(n4,'(a4)')'#HDL'
   if(kr==3)write(n7,'(a4)')'#HDL'
   if(kr==4)write(n1,'(a8)')'#OTHER-x'
   if(kr==4)write(n2,'(a8)')'#OTHER-y'
   if(kr==4)write(n3,'(a8)')'#OTHER-z'
   if(kr==4)write(n4,'(a6)')'#OTHER'
   if(kr==4)write(n8,'(a6)')'#OTHER'
   do it = 1,NDELS
      jt = it -1                                                     !print msd starting from time zero
      PSECS = TSTEP*dfloat(jt)*dfloat(KINTVL)
      WRITE(n1,19)PSECS,msd_mean_x(it,kr)/dfloat(ko)
      WRITE(n2,19)PSECS,msd_mean_y(it,kr)/dfloat(ko)
      WRITE(n3,19)PSECS,msd_mean_z(it,kr)/dfloat(ko)
      WRITE(n4,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      
      if(kr==1)WRITE(n5,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      if(kr==2)WRITE(n6,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      if(kr==3)WRITE(n7,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      if(kr==4)WRITE(n8,19)PSECS,msd_mean(it,kr)/dfloat(ko)
      
   end do
   write(n1,*)
   write(n2,*)
   write(n3,*)
   write(n4,*)
   
end do

write(n5,*)
write(n6,*)
write(n7,*)
write(n8,*)

close(n1) 
close(n2) 
close(n3) 
close(n4)

close(n5)
close(n6)
close(n7)
close(n8)

!Calculate the Einstein self-diffusion coefficient - MSD Linear Squares fit
write(n0,*)
write(n0,*)'Einstein self-diffusion coefficient'
write(n0,*)
do kr = kr_start,NRSHMAX
   nf = n4
   if(kr==1)then
      write(*,'(a5)')'#Bulk'
      write(n0,'(a5)')'#Bulk'
      nf = n5
   elseif(kr==2)then
      write(*,'(a4)')'#LDL'
      write(n0,'(a4)')'#LDL'
      nf = n6
   elseif(kr==3)then
      write(*,'(a4)')'#HDL'
      write(n0,'(a4)')'#HDL'
      nf = n7
   elseif(kr==4)then
      write(*,'(a6)')'#OTHER'
      write(n0,'(a6)')'#OTHER'
      nf = n8
   endif  

   call msd_LSF_Einstein(nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end,a,b,d,r,n)
   
   if(kr==1)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(*,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==2)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(*,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==3)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(*,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   elseif(kr==4)then
      write(n0,*)
      write(n0,*)'Number of points used in the fit = ',n
      write(n0,*)'Fit range (ps) = [',LSF_begin,'-',LSF_end,'['
      write(n0,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(n0,*)'Slope (10-5 cm2/s) = ',b*10.0d0
      write(n0,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(n0,*)'correlation coeficient = ',r 
      write(n0,*)'Standard deviation of the fit (10-5 cm2/s) = ',d*10.0d0
      write(*,*)
      write(*,*)'Ordinate intercept (10-5 cm2/s) = ',a*10.0d0
      write(*,*)'Self-Diffusion (slope/6) (10-5 cm2/s) = ',b*10.0d0/6.0d0
      write(*,*)'correlation coeficient = ',r 
   endif   
end do 

close(n0)

deallocate(nH2Oorg,NH2OidO)
deallocate(ntime)
deallocate(msd_mean,msd_mean_x,msd_mean_y,msd_mean_z)
deallocate(kv,nH2Omean)

deallocate(list_s)
deallocate(xyz)
deallocate(x,y,z)

!  9   FORMAT('#',15X,'time(ps)',7X,'msd(Ang**2/ps)'/) 
 19   FORMAT(9X,F14.5,4X,1PE14.4)
 29   FORMAT(9X,F14.5,4X,1PE14.4,1PE14.4)
 
   return

END SUBROUTINE pure_wat_LDL_HDL_msd
 



!Last subroutine - vmd_trj

SUBROUTINE vmd_trj(trjfile,ninput,inputformat,nbox,nens,cube,natms,nstep,dt,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,vmd_samp,nsoltypes,&
                   Zattype,mattype,ZMBIO,natHSh,natsol_id,oatsol_hs,ratsol_hs,resindex,resname,atomtype,chrg,nchains,nchain_nat,nchain_lab)
! Print system xyz confgs for vmd visualization 

    integer,intent(in)                              :: ninput,nbox
    character(len=7),intent(in)                     :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    integer,intent(in)                              :: nens        ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                         :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                              :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions
    real(kind=8),intent(in)                         :: dt
!    character(len=4),dimension(natms),intent(in)    :: atomname
     character(len=4),dimension(natms)              :: atomname
    integer,intent(in)                              :: vmd_samp
!    integer                                         :: vmd_samp
    
    integer,intent(in)                              :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in):: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)    :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    real(kind=8),dimension(natmsol),intent(in)      :: ZMBIO
    integer,intent(in)                              :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
    integer,dimension(natHsh),intent(in)            :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
    real(kind=8),dimension(natHsh),intent(in)       :: ratsol_hs   ! Hydration shell radius the solute atoms to analyze
    real(kind=8),dimension(natHsh),intent(in)       :: oatsol_hs   ! Hydration shell onset the solute atoms to analyze
    integer,dimension(natms),intent(in)             :: resindex 
    character(len=4),dimension(natms),intent(in)    :: resname
    character(len=4),dimension(natms),intent(in)    :: atomtype 
    real(kind=4),dimension(natms),intent(in)        :: chrg
    integer,intent(in)                              :: nchains     ! Number of protein chains (e.g., HgbS = 8)
    integer,dimension(nchains),intent(in)           :: nchain_nat  ! Number of atoms in each chain
    character(len=1),dimension(nchains),intent(in)  :: nchain_lab  ! Chain label (e.g., A, B, C etc)

! Local variables
!NG    real(kind=8)                             :: x(natms),y(natms),z(natms) 
!    real(kind=4)                               :: x(natms),y(natms),z(natms)
    real(kind=4),dimension(:),allocatable      :: x,y,z
!NG     real(kind=4)                            :: xx(natms),yy(natms),zz(natms) 
    real(kind=4),dimension(:),allocatable      :: xnt,ynt,znt
    real(kind=8)                               :: cell(3) 
    real(kind=8)                               :: cmxsol,cmysol,cmzsol
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    real(kind=8)                               :: cubeh(3) 
    integer,dimension(natms)                   :: natmid                    !atomic index (solute and solvent index)
    real(kind=8)                               :: timeprt  
    integer                                    :: max_ncfg
    integer                                    :: natcheck
    
    integer,dimension(:),allocatable         :: kwlab      !index of waters in the HSh
    integer                                  :: io
    
    integer      :: i, j, k, is
    integer      :: iw,jw,ki
    integer      :: n1,n2,n3,n4,n5,n6,n7
    integer      :: niw
    real(kind=4) :: occ
    
    integer                       :: at_id       ! solute atomic index to analyse the respective solvation environment
    integer                       :: nwat_hsh
    integer,dimension(5000)       :: wat_id_hsh  ! allcate a maximum of 500 waters in HSh
    real(kind=4)                  :: vx,vy,vz,rnm
    
!    logical filex

!xtc2fortran_trj
    character(len=15),intent(in)            :: trjfile    
!   vmd plugin variables
!   'xyx' is an array of type real*4 and will contain the coordinates in angstrom after the successful read.  
!   the coordinates are stored in the order x(1),y(1),z(1),x(2),y(2),z(2),...
!   so xyz(:) has to be at least of size 3*npart.

    integer*4                               :: npart, maxatom, handle(4), status
    real(kind=4)                            :: box(6)
    real(kind=4),dimension(:),allocatable   :: xyz
    character                               :: infile*200, intype*10
    integer                                 :: kp
!end xtc2fortran_trj

!xtc2fortran_trj

allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
allocate(xnt(natms),ynt(natms),znt(natms))

allocate(kwlab(nmolwat))        !array of waters in each environment; takes values of 0 (water not sampled) and 1 (water already sampled) 

IF(inputformat.eq.'GROMACS')THEN
   infile = TRIM(trjfile)
   intype = 'auto'
   npart  = -1
   handle(1) = -1
   handle(2) = -1
   handle(3) = -1
   handle(4) = -1
!set up everything and register all static plugins
!NG   call f77_molfile_init
   call f77_molfile_open_read(handle(1),npart,infile,intype)
   if (handle(1).lt.0) then
      print*,'file type unknown or not registered'
      stop
   else
      print*,'file successfully opened:'
      print*,'handle:',handle(1)
      print*,'npart: ',npart
      print*,'nsteps:',nstep
   end if
     
   if(npart /= natms)then
     write(*,*)'error - number of atoms is wrong in solv_inp.dat'
     stop
   endif  
ENDIF    

!end xtc2fortran_trj

occ = 1.00
wat_id_hsh = 0
nwat_hsh = 0
 cmxsol = 0.0d0; cmysol = 0.0d0; cmzsol = 0.0d0
vx= 0.0; vy=0.0; vz=0.0
rnm = 0.1                           !convert Ang to nm
is=1

!inquire(file='solv_out/sol_transl.xyz',exist=filex)
!if(.not.filex)call SYSTEM('mkdir -p solv_out')    

timeprt = 1.0d0                         !Maximum time lenght to sample configurations to print - normally between 1.0 to 10.0 ps for LARGE systems
if(natms>=5000)timeprt = 0.5d0          !maximum trj time lenght in ps to prt - set to 0.5 ps for LARGE systems
max_ncfg = IDINT(timeprt/(dt*1.0d-3))      !convert time in ps to maximum trj lenght to print - beyond max_ncfg configuration it will stop printing every vmd_samp                 

n2 = 200
open(n2,file='solv_out/log_visual.log',status='unknown',action='write')

write(*,*)
write(*,*)'Routine vmd_trj'
write(*,*)'Results will be printed to solv_out'
write(*,*)'Print xyz configurations every',vmd_samp,'time-steps'
write(*,*)
write(*,*)'the first',timeprt,' ps will be printed out...'
write(*,*)'number of configurations to print = ',max_ncfg/vmd_samp 
write(*,*)
write(*,*)'sol_transl_first_cfg.xyz: solute translated to the origin (0,0,0)'
write(*,*)'HSh_trj_transl.xyz: HSh printed out for visualization - solute translated to origin'
write(*,*)'HSh_trj_no_transl.xyz: HSh printed out for visualization - original xyz'
write(*,*)

write(n2,*)
write(n2,*)'Routine vmd_trj'
write(n2,*)'Results will be printed to solv_out'
write(n2,*)'Print xyz configurations every',vmd_samp,'time-steps'
write(n2,*)'the first',timeprt,' ps will be printed out...'
write(n2,*)'number of configurations to print = ',max_ncfg/vmd_samp 
write(n2,*)'sol_transl_first_cfg.xyz: solute translated to the origin (0,0,0)'
write(n2,*)'HSh_trj_transl.xyz: HSh printed out for visualization - solute translated to origin'
write(n2,*)'HSh_trj_no_transl.xyz: HSh printed out for visualization - original xyz'
write(n2,*)

!Solute translation
n1 = 150              ! solute + water full + ions .xyz format - print first configuration before and after translation
n3 = 300              ! Solute + HSh + ions: .xyz format - translate system - waters in the HSh - print every vmd_samp configuration
!No solute translation
n4 = 350              ! Solute + HSh + ions: .xyz format - do not translate system - waters in the HSh - print every vmd_samp configuration
!pdb
n5 = 450              ! Solute .pdb with force field charges - do not translate system - first configuration
n6 = 550              ! Solute (divide solute into two chains) - chain I .pdb with force field charges - do not translate system - first configuration
n7 = 650              ! Solute (divide solute into two chains) - chain II .pdb with force field charges - do not translate system - first configuration

open(n1,file='solv_out/sol_transl_first_cfg.xyz',status='unknown',action='write')
open(n3,file='solv_out/HSh_trj_transl.xyz',status='unknown',action='write')
open(n4,file='solv_out/HSh_trj_no_transl.xyz',status='unknown',action='write')

open(n5,file='solv_out/sol_chrg.pdb',status='unknown',action='write')
open(n6,file='solv_out/sol_chrg_prot_I.pdb',status='unknown',action='write')
open(n7,file='solv_out/sol_chrg_prot_II.pdb',status='unknown',action='write')

write(n1,'(i6)') natms
write(n1,'(a45)') 'First configuration before solute translation'

!Start vmd trajectory
do j=1,nstep
!   nwat_hsh = 0                                                    !Number of waters in the HSh
   if(mod(j,5000)==0)write(*,*)'starting step',j,' out of',nstep
   
!/Trajectory Reading Formats/    
   IF(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')THEN
      if(nens==1)then
         read(nbox,*)cell(1),cell(2),cell(3)
!         cell(2)   = cell(1)
!         cell(3)   = cell(1)
         cubeh     = cell/2.0d0
      elseif(nens==0)then
         cell(1) = cube(1)
         cell(2) = cube(2)
         cell(3) = cube(3)
         cubeh = cell/2.0d0
      endif   
!read solute and water coordinates
      read(ninput,*)natcheck
      if(natcheck.ne.natms)then
         write(*,*)'error - number of atoms from solv_inp is wrong'
         stop
      endif   
      do i = 1,natms
         read(ninput,*)natmid(i),atomname(i),x(i),y(i),z(i)
         if(j==1)write(n1,'(a4,3(3x,f12.5))')atomname(i),x(i),y(i),z(i)            !confg 1: xyz before solute translation to the origin 
      end do
   ELSEIF(inputformat.eq.'GROMACS')THEN
!xtc2fortran_trj
 
!   'status' is of type integer*4 and if set to zero on enter will instruct the reader to skip this frame.
!    on exit it will be set to zero if there was a problem or to 1 on success.         
      status = 1   ! status=1 on entry means read
      call f77_molfile_read_next(handle(1),npart,xyz(1),box,status);
      
      cell(1)=dble(box(1))
      cell(2)=dble(box(2))
      cell(3)=dble(box(3))
      cubeh = cell/2.0d0
      kp=1
      do i = 1,npart*3,3
         x(kp)=xyz(i)
         y(kp)=xyz(i+1)
         z(kp)=xyz(i+2)
         kp = kp + 1 
      end do
    
! write first configuration 
      if(j==1)then
         do i = 1,natms
            write(n1,'(a4,3(3x,f12.5))')atomname(i),x(i),y(i),z(i)
         end do
! write first configuration only - force field charge visualization - charges are written as beta temperatures in .pdb format  
!NG    do i=1,nmolsol*natmsol
!NG       write(n5,'(a6,i5,1x,a5,a4,i5,4x,3(f8.3),2f6.2)')'ATOM  ',i,atomtype(i),resname(i),resindex(i),x(i),y(i),z(i),occ,chrg(i)
!NG    end do
         is = 1
         do k = 1,nchains
            do i = 1,nchain_nat(k)
               write(n5,'(a6,i5,1x,a5,a4,a1,i5,2x,3(f8.3),2f6.2)')'ATOM  ',is,atomtype(is),resname(is),nchain_lab(k),resindex(is),x(is),y(is),z(is),occ,chrg(is)
               is=is+1
            end do
         end do
         is = 1
         do k = 1,nchains/2
            do i = 1,nchain_nat(k)
               write(n6,'(a6,i5,1x,a5,a4,a1,i5,2x,3(f8.3),2f6.2)')'ATOM  ',is,atomtype(is),resname(is),nchain_lab(k),resindex(is),x(is),y(is),z(is),occ,chrg(is)
               is=is+1
            end do
         end do
         is = 1
         do k = nchains/2+1,nchains
            do i = 1,nchain_nat(k)
               write(n7,'(a6,i5,1x,a5,a4,a1,i5,2x,3(f8.3),2f6.2)')'ATOM  ',is,atomtype(is),resname(is),nchain_lab(k),resindex(is),x(is),y(is),z(is),occ,chrg(is)
               is=is+1
            end do
         end do
      endif   
!
   ENDIF
!End - /Trajectory Reading Formats/ 
!end xtc2fortran_trj
!
!Find water molecules in the HSh at time zero - first MD configuration
!
   kwlab=0
   io=0
!
   if(j==1)then
      do i=nmolsol*natmsol+1,natms-nions,nwatsites              !Loop over water oxygens  
         io=io+1
!Hydration shell        
         do at_id=1,natHSh                                      !Loop over atomic species to analyze the environment         
            k = natsol_id(at_id)                                !Solute atoms to study hydration
            dx = x(i) - x(k)
            dy = y(i) - y(k)
            dz = z(i) - z(k)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2.0d0 + dy**2.0d0 + dz**2.0d0
            dr  = dsqrt(dr2)
            if(dr>=oatsol_hs(at_id).and.dr<=ratsol_hs(at_id))then     !water within hydration shell radius
            
               if(kwlab(io)==0)then
                  kwlab(io)=1                                  !oxygen id for correction of the nb of waters
               else
                  goto 10                                      !sample another water molecule            
               endif
            
!Hydration shell             
               nwat_hsh = nwat_hsh + 1 
               wat_id_hsh(nwat_hsh) = i                        !keep waters id to write HSh waters
            endif
         end do
         10 continue
      end do                                                    !End loop over oxygens
      write(*,*)
      write(*,*)'configuration',j,'number of HSh waters =',nwat_hsh
      write(*,*)
      write(n2,*)
      write(n2,*)'configuration',j,'number of HSh waters =',nwat_hsh
      write(n2,*)
   endif
   
! 
!Find centre of mass   
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
      call sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
   endif
!   
!System translation - put the solute centre of mass at the origin

!   if(j<=max_ncfg)then                                                   !print continuous trj - first max_ncfg configurations (file too big )
!Version 10
!   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
   if(j==1)then
      write(n1,'(i6)') natms
!      write(n1,'(a16,i7,a25)')'MD configuration',j,' after solute translation'
      write(n1,'(a44)') 'First configuration after solute translation'
   endif
!   
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
      write(n3,'(i5)') nmolsol*natmsol + nwat_hsh*nwatsites + nions                           ! print only solute, HSh waters, and ions
      write(n3,'(a8,i9)') 'MD confg',j
      write(n4,'(i5)') nmolsol*natmsol + nwat_hsh*nwatsites + nions                           ! print only solute, HSh waters, and ions
      write(n4,'(a8,i9)') 'MD confg',j
   endif 

!translate solute atoms
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
      i=1
      do while (i.le.nmolsol*natmsol)
!         x(i) = x(i) - cell(1)*dnint( (x(i)-cmxsol)/cell(1) )    
!         y(i) = y(i) - cell(2)*dnint( (y(i)-cmysol)/cell(2) )
!         z(i) = z(i) - cell(3)*dnint( (z(i)-cmzsol)/cell(3) )
         xnt(i) = x(i)     
         ynt(i) = y(i)
         znt(i) = z(i)
         x(i) = x(i) - cmxsol    
         y(i) = y(i) - cmysol
         z(i) = z(i) - cmzsol
         if(j==1)write(n1,'(a4,3(3x,f12.5))')atomname(i),x(i),y(i),z(i) 
         write(n3,'(a4,3(3x,f12.5))')atomname(i),x(i),y(i),z(i) 
         write(n4,'(a4,3(3x,f12.5))')atomname(i),xnt(i),ynt(i),znt(i) 
!         if(j==1)write(n4,'(i5,a3,2x,a5,i5,3f8.3,3f8.4)')&                                                              ! .gro format
!         resindex(i),resname(i),trim(adjustl(atomname(i))),i,rnm*x(i),rnm*y(i),rnm*z(i),vx,vy,vz  
         i=i+1
      enddo
   endif
    
!water - PBC  
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
      do i = nmolsol*natmsol + 1,natms-nions,nwatsites                          !Loop over oxygens  
!water oxygen atoms translation
         do k = 1,nwatsites
            iw = k - 1
            xnt(i+iw) = x(i+iw)    
            ynt(i+iw) = y(i+iw)
            znt(i+iw) = z(i+iw)
!            
            x(i+iw) = x(i+iw) - cmxsol    
            y(i+iw) = y(i+iw) - cmysol
            z(i+iw) = z(i+iw) - cmzsol
         end do   
!The centre of the MD box was originally (cubeh(1),cubeh(2),cubeh(3)) - positive coordinates
!The centre of the MD box after translation is (0,0,0) - positive and negative coordinates
!Periodic boundary conditions - water
         if (x(i).LE.-cubeh(1))then
            do k = 1,nwatsites
               iw = k - 1
               x(i+iw) = x(i+iw) + cell(1) 
            end do   
         endif   
         if (x(i).GT.cubeh(1)) then
            do k = 1,nwatsites
               iw = k - 1
               x(i+iw) = x(i+iw) - cell(1)
            end do  
         endif   
        if (y(i).LE.-cubeh(2))then
            do k = 1,nwatsites
               iw = k - 1
               y(i+iw) = y(i+iw) + cell(2) 
            end do   
         endif   
         if (y(i).GT.cubeh(2)) then
            do k = 1,nwatsites
               iw = k - 1
               y(i+iw) = y(i+iw) - cell(2)
            end do  
         endif
         if (z(i).LE.-cubeh(3))then
            do k = 1,nwatsites
               iw = k - 1
               z(i+iw) = z(i+iw) + cell(3) 
            end do   
         endif   
         if (z(i).GT.cubeh(3)) then
            do k = 1,nwatsites
               iw = k - 1
               z(i+iw) = z(i+iw) - cell(3)
            end do  
         endif
      end do       !end loop over oxygen atoms
   endif   

!print water
   if(j==1)then
      do i = nmolsol*natmsol + 1,natms-nions,nwatsites
         do k = 1,nwatsites
            iw = k - 1
            write(n1,'(a4,3(3x,f12.5))')atomname(i+iw),x(i+iw),y(i+iw),z(i+iw)
         end do
      end do   
   endif   
   
!HSh - xyz   
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then           !print HSh waters
      do i = 1,nwat_hsh
         jw = wat_id_hsh(i)
         write(n3,'(a4,3(3x,f12.5))')atomname(jw),x(jw),y(jw),z(jw)
         write(n3,'(a4,3(3x,f12.5))')atomname(jw+1),x(jw+1),y(jw+1),z(jw+1)
         write(n3,'(a4,3(3x,f12.5))')atomname(jw+2),x(jw+2),y(jw+2),z(jw+2)
!         
         write(n4,'(a4,3(3x,f12.5))')atomname(jw),xnt(jw),ynt(jw),znt(jw)
         write(n4,'(a4,3(3x,f12.5))')atomname(jw+1),xnt(jw+1),ynt(jw+1),znt(jw+1)
         write(n4,'(a4,3(3x,f12.5))')atomname(jw+2),xnt(jw+2),ynt(jw+2),znt(jw+2)
      end do                                                 
   endif
!HSh - gro
!   if(j==1)then
!      niw = nmolsol*natmsol
!      do i = 1,nwat_hsh
!         niw = niw + 1
!         jw = wat_id_hsh(i)
!         write(n4,'(i5,a3,2x,a5,i5,3f8.3,3f8.4)')&
!         resindex(jw),resname(jw),trim(adjustl(atomname(jw))),niw,rnm*x(jw),rnm*y(jw),rnm*z(jw),vx,vy,vz
!         niw = niw + 1
!         write(n4,'(i5,a3,2x,a5,i5,3f8.3,3f8.4)')&
!         resindex(jw+1),resname(jw+1),trim(adjustl(atomname(jw+1))),niw,rnm*x(jw+1),rnm*y(jw+1),rnm*z(jw+1),vx,vy,vz
!         niw = niw + 1
!         write(n4,'(i5,a3,2x,a5,i5,3f8.3,3f8.4)')&
!         resindex(jw+2),resname(jw+2),trim(adjustl(atomname(jw+2))),niw,rnm*x(jw+2),rnm*y(jw+2),rnm*z(jw+2),vx,vy,vz       
!      end do            
!   endif   
   
!ions - PBC
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
      do i = natms-nions+1,natms                          !Loop over oxygens  
!ions translation
            xnt(i) = x(i) 
            ynt(i) = y(i)
            znt(i) = z(i)
!            
            x(i) = x(i) - cmxsol    
            y(i) = y(i) - cmysol
            z(i) = z(i) - cmzsol 
!The centre of the MD box was originally (cubeh(1),cubeh(2),cubeh(3)) - positive coordinates
!The centre of the MD box after translation is (0,0,0) - positive and negative coordinates
!Periodic boundary conditions - ions
                  
         if (x(i).LE.-cubeh(1))then
               x(i) = x(i) + cell(1) 
         endif   
         if (x(i).GT.cubeh(1)) then
               x(i) = x(i) - cell(1)
         endif   
         if (y(i).LE.-cubeh(2))then
               y(i) = y(i) + cell(2) 
         endif   
         if (y(i).GT.cubeh(2)) then
               y(i) = y(i) - cell(2)  
         endif
         if (z(i).LE.-cubeh(3))then
               z(i) = z(i) + cell(3)    
         endif   
         if (z(i).GT.cubeh(3)) then
               z(i) = z(i) - cell(3)
         endif
      end do 
   endif   
         

!print ions
!xyz
   if((j==1.or.mod(j,vmd_samp)==0).and.(j<=max_ncfg))then
      ki = nmolsol*natmsol + nmolwat*nwatsites
      do i = 1,nions
         ki = ki + 1
         if(j==1)write(n1,'(a4,3(3x,f12.5))')atomname(ki),x(ki),y(ki),z(ki)
         write(n3,'(a4,3(3x,f12.5))')atomname(ki),x(ki),y(ki),z(ki)
         write(n4,'(a4,3(3x,f12.5))')atomname(ki),xnt(ki),ynt(ki),znt(ki)
      end do
   endif 
!gro 
!   if(j==1)then
!      ki = nmolsol*natmsol + nmolwat*nwatsites
!      do i = 1,nions
!         niw = niw + 1
!         ki = ki + 1
!         write(n4,'(i5,a3,2x,a5,i5,3f8.3,3f8.4)')&
!         resindex(ki),resname(ki),trim(adjustl(atomname(ki))),niw,rnm*x(ki),rnm*y(ki),rnm*z(ki),vx,vy,vz 
!      end do
!      write(n4,'(1x,4f10.5)')cell
!      write(n4,*)
!   endif
!   
!!xtc2fortran_trj
!   if(status.eq.0) then
!      write(*,*)'error on reading trajectory file - step',i
!      stop
!   endif
!end xtc2fortran_trj      
      
end do          !end loop over time-steps

!xtc2fortran_trj
call f77_molfile_finish
write(*,*)'finished reading xtc trj...'
!end xtc2fortran_trj

deallocate(xyz)
deallocate(x,y,z)
deallocate(xnt,ynt,znt)

   return

END SUBROUTINE vmd_trj



END MODULE inp_out_solv
