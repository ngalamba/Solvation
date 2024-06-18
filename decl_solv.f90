MODULE decl_solv

implicit none

!general
integer    :: i,j,k
!files
integer :: nrep=7,ngrmx=8,ntrjunf=9,nbox=10,npurew=11,nsolprm=12,ngroitp=13,ncell=14,nLJitp=15
!input file

!input MD parameters
character(len=7) :: inputformat                         ! input format type: TINKER, GROMACS, BOMD
character(len=9) :: ffname                 ! ffname: GROMOS, AMBER, CHARMM, AMOEBA, OPLS, etc
character(len=15):: trjfile,topfile,LJfile,boxfile      ! GROMACS trj file (.xtc); GROMACS topology file (.top or .itp)
character(len=14):: unfoldtrjfile                         ! internal - unfolded trajectory file name for Einstein diffusion trj_unfold.xtc
integer          :: ksLJ                                ! 0: do not read LJ 12-6 parmts; 1: read
!integer          :: mxtctrj                            ! 0: build fortran unformatted input trj file from .xtc; 1: fortran unformatted trj file available
integer          :: nstep                        ! Number of steps (frames)
real(kind=8)     :: dt                           ! time-step (fs)
integer          :: dfr                          ! Number of timesteps between two frames
integer          :: natms                        ! Total number of atoms
integer          :: nsyst                        ! 0 Pure Water; 1 Solute and Water System
integer          :: nens                         ! ensemble 0: NVE, NVT; 1: NPT
integer          :: nmolsol                      ! Number of molecules of solute
integer          :: natmsol                      ! Number of atomos per solute molecule
character(len=7) :: WATMODEL                     ! Water model used in aqueous solution systems
integer          :: nmolwat                      ! Number of water molecules
integer          :: nions                        ! Number of ions
integer          :: nwatsites                    ! Number of water sites (SPC/E = 3; TIP4P-Ew = 4)
integer          :: kbiowater                    ! Flag for bio-bulk water histograms
integer          :: biowsamp                     ! time freq. to sample
real(kind=8)     :: biowRadius                   ! bio water cut-off (Ang)
integer          :: kwtcfcrsh                    ! Flag for radial reor. tcfs 
integer          :: m_time,mtdelay               ! time freq. to sample water OH groups in radial shells for reor. tcf calculation; delay time
integer          :: m_time_hb,mtdelay_hb         ! time freq. to sample water OH groups in radial shells for reor. tcf calculation; delay time
real(kind=8)     :: rwtcfbulk                    ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk reorientational tcf
integer          :: khbtcf                       ! Flag for HB time-correlation functions  
integer          :: hbs_time,hbtdelay            ! time freq. to sample HBs and HB tcfs delay time
integer          :: kwtetra                      ! Flag for tetrahedrality calculation - solute heavy atoms excluded
real(kind=8)     :: rwtetbulk                    ! onset (Ang) for bulk water relative to the solute centre of mass
integer          :: kwsoltetra                   ! Flag for tetrahedrality calculation - solute heavy atoms included
real(kind=8)     :: rwsoltetbulk                 ! onset (Ang) for bulk water relative to the solute centre of mass
integer          :: tetsampw,tetsampws           ! interval for computing the tetrahedrality
integer          :: Voronsampw                   ! interval for computing Voronoi for pure water
integer          :: KWLSI                        ! 0: do not calculate; 1: calculate - LSI for solutions    
real(kind=8)     :: rwlsibulk                    ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk LSI
integer          :: LSI_sample                   ! LSI sampling interval - time-steps
integer          :: kscontact,intsf              ! Flag for contact map between proteins; trajectory sampling frequency
real(kind=8)     :: rcont                        ! distance for contact map - residues at a r < rcont are accounted
integer          :: kMIC,kMIC_                   ! apply Minimum Image Convention or not to proteins
integer          :: intsf_                       ! contact map trajectory sample frequency (steps) [routine potent_map]
real(kind=8)     :: rcont_                       ! distance cut-off for contact map - residues at a r < rcont are accounted [routine potent_map]
integer          :: kspotential,ncontsamp        ! Flag for contact map electrostatics and van der Waals; trajectory sampling frequency
integer          :: mcontsamp                    ! trajectory sampling frequency - protein-water and water-water potential energy
integer          :: ksprotwat                    ! Flag for calling protein-water routine - calculated the protein-water and water-water PE
integer          :: kprot_wat,kwat_wat           ! Flags for calculation protein-water and/or wate-water PE
integer          :: kumbrella_trj                ! 0: Trjs are not from umbrella sampling 1: Trjs are from umbrella sampling [no averages over trjs will be calculated]
integer          :: kntrj                        ! Number of MD trajectories to build a mean electrostatic contact map 
!real(kind=8)     :: Coul_th,LJ126_th             ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!real(kind=8)     :: Coul_th_min,LJ126_th_min        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!real(kind=8)     :: Coul_th_max,LJ126_th_max        ! Coulombic and Lennard-Jones 12-6 thresholds for printing res-res interactions
!Version 15 - High energy RES-RES contacts will be chosen based on a single criterion so that Coulombic and van der Waals pairs are the same
real(kind=8)     :: Coul_th_min,Coul_th_max      ! Potential energy threshold for printing res-res interactions [Min,Max]
real(kind=8)     :: LJ126_th_min,LJ126_th_max 
integer          :: nattp                        ! Number of distinct atomtypes in the force field file (ffnonbonded.itp - Lennard-Jones parm)
integer          :: ncrossLJ                     !number of LJ126 crossed IJ terms - for reading non-geometric C12
integer          :: kdifwsol                     ! Flag for diffusion calculation - msd Einstein
integer          :: difsampw,msd_delay           ! time freq. to sample waters for diffusion calculation; delay time
real(kind=8)     :: difbulk                      ! onset (Ang) for bulk water relative to the solute centre of mass - for bulk diffusion
integer          :: krmvmd,vmd_samp              ! Flag for printing vmd-xyz snapshots and print freq. (time-steps) 
!Pure water
!integer          :: nbwmol                      ! Number of wat. mol.; number of MD time-steps; delay time to calculate pure water's freq. tcf
!integer          :: nwatpoints                  ! Number of water sites (SPC/E = 3; TIP4P-Ew = 4) - pure water
!real(kind=8)     :: wside                       ! Pure Water MD box side length/Ang
integer          :: ktcfpureW                    ! Flag for reorientational tcfs - pure water
integer          :: ktcfhbW                      ! Flag for reorientational tcfs - pure water - use HB definition for OH tcf deconvolution
!integer          :: tcfframe,tcf_samp           ! reorentional tcf delay time; reorientational tcf sampling time - pure water
integer          :: khbpureW                     ! Flag for HB time-correlation functions - pure water
integer          :: hbtframe,hbs_samp            ! HB tcfs delay time; time freq. to sample HBs - pure water
integer          :: ktetpureW                    ! Flag for tetrahedrality calculation of pure water
real(kind=4)     :: rcut                         ! cut-off (Ang) for neighbor search in Voronoi
integer          :: kVoronoipureW                ! Flago for Voronoi calculation of pure water
integer          :: kdifpureW                    ! 0: do not calculate; 1: calculate - diffusion of pure water
integer          :: unfold_trj                   ! 0: unfold trj for msd calculation; 1: do not unfold trajectory
logical          :: ldiff                        ! diffusion calculation - will re-name the (unfolded) unformatted trj file
integer          :: QMMM_pureW                   ! 0: do not print; 1: print QMMM conf. for pure water
integer          :: Nenv,Nenvrm,Nenvt            ! Number of environments to calculate the orientational tcfs/diffusion for pure water (LDL/HDL)
real(kind=8)     :: R5Tet                        ! cut-off to study tetrahedrality in different environments in pure water
real(kind=8)     :: R5Dif                        ! cut-off to study diffusion in different environments in pure water
!integer          :: nbulk                        ! 0 skip tetrahedrality calculation for bulk water; 1: include bulk tetrahedrality
integer          :: Nenv_hb,Nenv_hb_HB           ! Number of environments to calculate the orientational tcfs [HB deconvolution version 24/25]
integer          :: kLDL_HDL,kLDL_HDL_           ! 0:read LDL/HDL input list; 1:calculate LDL/HDL on the fly
real(kind=8)     :: LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2         !
real(kind=8)     :: LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_
integer          :: klsipureW                ! 0: do not calculate; 1: calculate - LSI tcfs (continuous and intermittent) - pure water
integer          :: lsitframe,lsi_samp,Nenv_lsi
real(kind=8)     :: LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II

!input solute
integer                                       :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
character(len=2),dimension(:),allocatable     :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
real(kind=8),dimension(:),allocatable         :: mattype     ! Mass of atoms of type (C, H, O, N etc)
real(kind=8),dimension(:),allocatable         :: ZMBIO       ! Solute atomic masses
integer,dimension(:),allocatable              :: Natresid    ! Number of atoms of each residue
integer                                       :: niontypes   ! Number of ionic distinct species (Na+, Cl- etc)
integer,dimension(:),allocatable              :: niontype    ! Number of ions of type (Na+, Cl- etc)
real(kind=8),dimension(:),allocatable         :: miontype    ! Mass of ions of type (Na+, Cl- etc)
character(len=2),dimension(:),allocatable     :: Ziontype    ! Atomic symbol of ion types (Na+, Cl- etc)
integer                                       :: nchains     ! Number of protein chains (e.g., HgbS = 8)
integer,dimension(:),allocatable              :: nchain_nat  ! Number of atoms in each chain
character(len=1),dimension(:),allocatable     :: nchain_lab 
character(len=1),dimension(:),allocatable     :: chain_label_res ! Chain label (e.g., A, B, C etc) internal code use
integer,dimension(:),allocatable              :: nb_res_ch   ! number of residues of each chain
integer                                       :: natHSh      ! Number of solute atoms to analyse the respective solvent environment
integer,dimension(:),allocatable              :: natsol_id   ! Indexes of atoms to analyse the respective solvent environment
real(kind=8),dimension(:),allocatable         :: ratsol_hs   ! Hydration shell radius of the solute atoms to analyze
real(kind=8),dimension(:),allocatable         :: oatsol_hs   ! Hydration shell onset of the solute atoms to analyze
integer                                       :: at_id       ! solute atomic index to analyse the respective solvation environment
real(kind=8),dimension(:),allocatable         :: C6kjmol,C12kjmol  !Solute Lennard-Jones parameters kJ/mol
integer,dimension(:,:),allocatable            :: res_cont_pair_t0  !List of residue-residue in the contact map
integer                                       :: nres_cont_pair_t0 !Number of contacts in the contact map
real(kind=8),dimension(:,:),allocatable       :: C6_IJ,C12_IJ      ! Lennard-Jones C6 and C12 crossed IJ terms solute-solute
real(kind=8),dimension(:),allocatable       :: C6_IW,C12_IW      ! Lennard-Jones C6 and C12 crossed IJ terms solute-water
!trajectory
real(kind=8) :: cell(3)
real(kind=8) :: cube(3)
!NG real(kind=8),dimension(:),allocatable :: x,y,z
real(kind=4),dimension(:),allocatable :: x,y,z
!calculations
real(kind=8) :: dx,dy,dz,dr2           ! interatomic distances
!atom type
character(len=4),dimension(:),allocatable     :: resname
real(kind=4),dimension(:),allocatable         :: chrg
character(len=4),dimension(:),allocatable     :: atomname
character(len=4),dimension(:),allocatable     :: atomname_d
character(len=4),dimension(:),allocatable     :: atomtype
character(len=3),dimension(:),allocatable     :: segname
integer,dimension(:),allocatable              :: resindex    
integer                                       :: nb_res
character(len=4),dimension(:),allocatable     :: name_res
integer,dimension(:),allocatable              :: numb_res               ! Number of the residue: from 1 to the total number of residues of each chain
!parameters
real(kind=8),parameter :: pi=3.14159265d0

END MODULE decl_solv
