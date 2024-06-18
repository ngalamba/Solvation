!
! *************************************************************************************************
! Program Solvation
! 
! Dynamic and Structural Solvation Analysis
!
! Nuno Galamba 
! 
! January 2019                                       Version 36.0
!
! *************************************************************************************************
!
!               Future modifications of the code - version 36.0:
!               ===============================================
!
!               . add bulk definition based on the distance to any solute atom instead of the center of mass of the solute
!               . print mean number of water molecules within the HSh of each solute atom
!
! Version 36.0  . implemented routine to distinguish biological and bulk water in protein/co-soute/Water (e.g. Ubiquitin/DES/Water)
!                 biological water: < cut-off distance from the protein heavy atoms; bulk water > cut-off distance from the protein heavy atoms 
!
! Version 35.0  . implemented Voronoi polyhedra calculation in pure water - not validated
!               . allowed printing of tetrahedrons for visualization in the routine wat_tetra/tetrahedron
!
! Version 34.0  . calculating the tetrahedrality of water molecules with 4 and 5 neighbors within some cut-off - pure water only
!
! Version 33.0  . calculating the diffusion of water molecules with 3 (or less) 4 (or less) and 5 (or more) neighbors within some cut-off - pure water only
!
! Version 32.0  . print residue-residue distances to build a countour map - interaction - distance - potential energy             
!
! Version 31.0  . implemented the calculation of water-water pair interaction energies P(O...On) with n = 1 to 5, in solute-water systems
!
! Version 30.0  . implemented the calculation of the LSI for solute-water systems
!
! Version 29.0  . implemented the routine top_gro_sol to handle non-protein solutes topology reading
!
! Version 28.0  . implemented LSI tcfs for pure water - continuous and intermittent
!
! Version 27.0  . implemented HB tcfs for pure water - Luzar and Chandler
!                 HB tcfs for different environments based on the LSI index 
!
! Version 26.0  . implemented on the fly calculation of the LSI index for pure water in the routine pure_wat_LDL_HDL_otcf_HB
!                 implemented on the fly calculation of the LSI index for pure water in the routine pure_wat_LDL_HDL_otcf
!                 routine pure_wat_LDL_HDL_otcf - calculate LSI distribution
!
! Version 25.0  . corrected LDL/HDL labels in Leg tcfs output
!                 Introduce a new routine pure_wat_LDL_HDL_otcf_HB that allows deconvolution of the otcfs based on:
!                       (i) OH dependent short/long HBs (ii) OH independent HBonded/non-HBonded OH group (iii) Numb. of neighbors < 5 or > 5 up to r Ang
!                 (i)   For each molecule assess the OH forming a short and the OH forming a long HB - calcultate the respective otcfs
!                 (ii)  For each OH assess whether it is HBonded or not and calculate the the respective otcfs - use stricter HB definition
!                 (iii) For each water assess whether it has < 5 or >=5 water neighbors up to r and calculate respective otcfs
!                 
! Version 24.0  . coded the routines pure_wat_LDL_HDL_otcf and pure_wat_LDL_HDL_msd to calculate the orientational tcfs and diffusion for: 
!                 bulk (all), LDL, HDL, and other water molecules in pure water - Lars Pettersson and Gaia Camisasca collaboration - June 2018
!
! Version 23.0  . re-named the input varibale ncontsamp to mcontsamp on the protein-water and water-water rotine: prot_water_umbrella
!                 the same variable was used in potent_map and pot_map_umbrella
!
! Version 22.0  . calculate protein-water and/or water-water potential energy for single trajectories - NO EWALD SUM
!
! Version 21.0  . print standard deviations of the number of waters found in the different regions used to calculate the orientational dyn and the msd
!
! Version 20.0  . Defined a vdW threshold for studying the effect of vdW interactions in the potential energy
!
! Version 19.0  . print umbrella sampling intermolecular energy analysis as a function of the reaction coordinate in addition to the COM-COM
!
! Version 18.0  . corrected a minor bug regarding the time corresponding to each origin when ncontsamp =/ 1
!                 increased the number of decimal places for time/ns in routine potent_map
!
! Version 17.0  . corrected a minor bug on the stdevs in routines potent_map and pot_map_umbrella when sampling every origin of the trj [ncontsamp = 1]
!
! Version 16.0  . Exended the analysis to umbrella sampling trajectories
!                 Fixed a number of bugs in the energy calculations for SMD trajectories
!
! Version 15.0  . Energy decomposition based on a single threshold interval and only for the pot en
!                 Coulomb and van der Waals en are printed if the pot en (Coulomb + van der Waals) is above or below threshold
!                 Thus the highest or lowest van der Waals contact pairs will not be necessarily printed
!
! Version 14.0  . Fixed a bug in the Lennard-Jones potential map when more than one trajectory (SMD) is analysed
!               . Allowed to include threshold res-res interactions based on the largest energy observed through the full trj 
!                 . instead of the time zero t0
!
! Version 13.0  . Fixed a bug in the calculation of the msd and diffusion of solvation waters
!                 A PBC trj is now used to find waters in the solvation shell while the unfolded trj is only used to compute the msd
!               . Introduced a maximum number of waters in the bulk environment to be studied in all solvation routines  
!
! Version 12.0  . Protein residue-residue electrostatic and van der Waals potential map - allow averaging over different trjs
!
! Version 11.0  . Protein residue-residue electrostatic and van der Waals potential map
!
! Version 10.0  . Protein residue-residue contact map                  
!
! Version 9.0   . identify reference solute atoms neigbors (solute atoms) for non-tetrahedral waters 
!                 reorientational tcf calculation
!
! Version 8.0   . read trj directly from the .xtc
!               . define different input shell radius for different solute atoms 
!               . account for waters in the hydration shells of several atoms of solute
!               . account for time origins where no water molecules are found in a given environment
!               . read topology for any system based on labels/rewind
!
! Version 7.0   . extended the code to study the solvation environment of multiple groups  
!                 (e.g., amino acids) of the solute (e.g., proteins) in ionic aqueous solutions
!               . the order of the species in the input xyz file is: protein,water,ions
!               . converted the internal x,y,z to single precision to reduce the allocated  
!                 memory; xtc2fortran_trj writes the xyz in single precision (binary)
!
!               . deleted the following routines from previous versions: 
!                 wat_env_tcfs;  wat_env_tcfs_mem; wat_env_tcfs_t
!               . the former wat_env_tcfs_ram routine was renamed wat_env_orient_tcfs
!            
!                 Orientational tcf:   (i)  Tetrahedral waters have 4 waters neighbors closer than
!                                           any solute atom - ions are neglected (i.e., not analyzed)
!                                      (ii) Non-tetrahedral waters have a solute atom closer than 
!                                           the 4th nearest water neighbor - ions are neglected
!                 
!                 Tetrahedrality:      (i)  Tetrahedral waters have 4 waters neighbors closer than
!                                           any solute atom or ion 
!                                      (ii) Non-tetrahedral waters have a solute atom or ion closer  
!                                           than the 4th nearest water neighbor - tetrahedrality is
!                                           computed only for waters - solute and ions neglected
!                 
!                 Tetrah_water-solute: (i)  Tetrahedral waters have 4 waters neighbors closer than
!                                           any solute atom or ion 
!                                      (ii) Non-tetrahedral waters have a solute atom or ion closer  
!                                           than the 4th nearest water neighbor
!                                           tetrahedrality is computed considering the solute heavy
!                                           atoms and ions as possible tetrahedral vertices
!
!               . implemented orientational dynamics, tetrahedrality, and Einstein diffusion for 
!                 pure water
!               . implemented the mean square displacement and Einstein water diffusion in solution
!
! Version 6.0   . introduced the routine wat_env_tcfs_t which calculates the tcfs for waters in 
!               . a given environment both only at t0 and at every t
!
! Version 5.0   . introduced the routine wat_env_tcfs_ram with further lower memory requirements 
!                 than wat_env_tcfs and wat_env_tcfs_mem; it no longer calculates the tcfs stdev
!                 if enough RAM available use wat_env_tcfs_mem
!
! Version 4.0   . introduced the routine wat_env_tcfs_mem with lower memory requirements than
!                 wat_env_tcfs and added the calculation of the Legendre polynomials tcfs 1 and 3
!
! Version 3.0   . introduced the routines wat_tetra and wat_tetra_sol
!
! *************************************************************************************************
!
PROGRAM Solvation

![trajectory inputformat = TINKER  xyz: ********.arc ]
![trajectory inputformat = GROMACS xyz: ********.xtc ]
![trajectory inputformat = BOMD    xyz: ********.xyz ]

![trajectory atomic order: solute atoms, water atoms [O,H1,H2,M], ions]

!===========
!input files
!===========
! solv_inp/solv_inp.dat              !MD parameters
! solv_inp/sol_sites.dat             !solute parameters
! solv_inp/*******.arc               !xyz (Ang)                   [TINKER/BOMD]
! solv_inp/*******.ang               !box side lengths (Ang)      [TINKER/BOMD]
! solv_inp/*******.xtc               !xyz file (Ang)              [GROMACS]
! solv_inp/*******.top[itp]          !topology file               [GROMACS]

USE decl_solv
USE inp_out_solv

implicit none

integer :: nversion
logical :: filex

nversion = 36

write(*,*)
write(*,'(5x,a33)')'*********************************'
write(*,'(5x,a28,1x,i3)')'Program Solvation - version', nversion
write(*,'(5x,a33)')'*********************************'
write(*,*)

open(nrep,file='calc_details.out',status='unknown',action='write')

write(nrep,*)
write(nrep,'(5x,a31)')'******************************'
write(nrep,'(5x,a26,1x,i3)')'Program Solvation version', nversion
write(nrep,'(5x,a31)')'******************************'
write(nrep,*)

!general output data

!NG inquire(file='solv_out/sol_transl.xyz',exist=filex)
inquire(file='solv_out/box_out.dat',exist=filex)
if(.not.filex)call SYSTEM('mkdir -p solv_out')

!input data file
open(ngrmx,file='solv_inp/solv_inp.dat',status='old',action ='read')

call read_parm(ngrmx,ffname,inputformat,trjfile,unfoldtrjfile,unfold_trj,topfile,LJfile,ksLJ,boxfile,nsyst,nens,cube,nstep,dt,dfr,nmolsol,natmsol,&
               nmolwat,nions,nwatsites,kbiowater,biowsamp,biowRadius,kwtcfcrsh,mtdelay,m_time,mtdelay_hb,m_time_hb,rwtcfbulk,khbtcf,hbtdelay,hbs_time,kwtetra,rwtetbulk,&
               kwsoltetra,rwsoltetbulk,tetsampw,tetsampws,Voronsampw,khbpureW,hbtframe,hbs_samp,kscontact,intsf,rcont,kspotential,kumbrella_trj,kntrj,&
               ncontsamp,kMIC,kMIC_,intsf_,rcont_,Coul_th_min,Coul_th_max,LJ126_th_min,LJ126_th_max,nattp,ncrossLJ,mcontsamp,ksprotwat,kprot_wat,kwat_wat,&
               kdifwsol,msd_delay,difsampw,difbulk,krmvmd,vmd_samp,ktcfpureW,ktcfhbW,ktetpureW,kVoronoipureW,rcut,Nenvt,R5Tet,kdifpureW,QMMM_pureW,Nenv,Nenvrm,R5Dif,Nenv_hb,&
               kLDL_HDL,kLDL_HDL_,LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2,Nenv_hb_HB,LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_,&
               klsipureW,lsitframe,lsi_samp,Nenv_lsi,LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II,KWLSI,LSI_sample,rwlsibulk,WATMODEL)

write(nrep,*)'==========================='
write(nrep,*)'TRAJECTORY INPUT PARAMETERS'
write(nrep,*)'==========================='
write(nrep,*)
write(nrep,*)'Format input = ',inputformat
write(nrep,*)'Force Field = ',ffname
write(nrep,*)'trajectory file = ',trjfile
if(inputformat.eq.'GROMACS')write(nrep,*)'topology file = ',topfile
if(inputformat.eq.'TINKER ')write(nrep,*)'box side lengths file = ',boxfile
if(nsyst==0)write(nrep,*)'System type = PURE WATER'
if(nsyst==1)write(nrep,*)'System type = SOLUTION'
if(nens==0)write(nrep,*)'Ensemble type = NVE/NVT'
if(nens==1)write(nrep,*)'Ensemble type = NPT'
write(nrep,*)'Number of frames = ',nstep
write(nrep,*)'Time-step (fs) = ',dt
write(nrep,*)'Interval between frames = ',dfr
!
write(*,*)'Format input = ',inputformat
write(*,*)'Force Field = ',ffname
write(*,*)'Number of frames to analyze = ',nstep
write(*,*)'Time-step (fs) = ',dt
write(*,*)'Interval between frames = ',dfr
!
if(nsyst==1)then                                         !Solution
   write(nrep,*)
   write(nrep,*)'================='
   write(nrep,*)'SOLUTION ANALYSIS'
   write(nrep,*)'================='
   write(nrep,*)
   write(nrep,*)'Number of solute molecules = ',nmolsol
   write(nrep,*)'Number of atoms/solute molecule = ',natmsol
   write(nrep,*)'Number of water molecules = ',nmolwat
   write(*,*)   'Number of water molecules = ',nmolwat
   write(nrep,*)'Number of water sites = ',nwatsites
   write(*,*)   'Number of water sites = ',nwatsites
   write(nrep,*)'Number of ions = ',nions
   if(kwtcfcrsh==1)then
      write(nrep,*)'Calculate water orient. tcfs = YES'
      write(nrep,*)'   delay-time (ps) = ',mtdelay
      write(nrep,*)'   sample freq. (steps) = ',m_time
      write(nrep,*)'   reorientational tcf bulk water onset (Ang)', rwtcfbulk
   endif
   if(khbtcf==1)then
      write(nrep,*)'Calculate Solvent-Water HB dynamics = YES'
      write(nrep,*)'   Solvent-Water HB dynamics delay-time (ps)',hbtdelay 
      write(nrep,*)'   sample freq. (steps)',hbs_time
   endif
   if(kwtetra==1)then
      write(nrep,*)'Calculate Solvent-Water tetrahedrality = YES'
      write(nrep,*)'   Tetrahedrality bulk water onset (Ang)',rwtetbulk
   endif
   if(kwsoltetra==1)then
      write(nrep,*)'Calculate Solute-Water tetrahedrality = YES'
      write(nrep,*)'   Tetrahedrality bulk water onset (Ang)',rwsoltetbulk
   endif  
   if(KWLSI==1)then
      write(nrep,*)'Calculate Solute-Water LSI = YES'
      write(nrep,*)'   LSI bulk water onset (Ang)',rwlsibulk
   endif     
   if(kscontact==1)then
      write(nrep,*)'Calculate protein contact map = YES'
      write(nrep,*)'   residue-residue sample freq. (steps)',intsf
      write(nrep,*)'   residue-residue distance cut-off (Ang)',rcont
   endif
   if(kspotential==1)then
      write(nrep,*)'Calculate protein Coul. and van der Waals contact map = YES'
      if(kMIC==0)write(nrep,*)'MIC will not be applied - no jump protein'
      if(kMIC==1)write(nrep,*)'MIC will be applied'      
      write(nrep,*)'   number of trajectories to average',kntrj
      write(nrep,*)'   energy calculation frequency (steps)',ncontsamp
      write(nrep,*)'   residue-residue sample freq. (steps)',intsf_
      write(nrep,*)'   residue-residue distance cut-off (Ang)',rcont_
      write(nrep,*)'   Coulomb energy threshold (kJ/mol)',Coul_th_min,'-',Coul_th_max 
!      write(nrep,*)'   van der Waals energy threshold (kJ/mol)',LJ126_th_min,'-',LJ126_th_max
      write(nrep,*)'   potential energy threshold (kJ/mol)',Coul_th_min,'-',Coul_th_max 
   endif
   if(ksprotwat==1)then
      write(nrep,*)'Calculate protein-water potential energy = YES'
      if(kMIC_==0)write(nrep,*)'MIC will not be applied - no jump protein'
      if(kMIC_==1)write(nrep,*)'MIC will be applied'  
   endif 
   if(kdifwsol==1)then
      write(nrep,*)'Calculate Solvent-Water diffusion = YES'
      write(nrep,*)'   delay-time (ps)',msd_delay
      write(nrep,*)'   sample freq. (steps) = ',difsampw
      write(nrep,*)'   diffusion bulk water onset (Ang)',difbulk
   endif
   if(krmvmd==1)then
      write(nrep,*)'Print VMD trj (xyz) = YES'
      write(nrep,*)'   sample freq. (steps)',vmd_samp
   endif
   write(nrep,*)
elseif(nsyst==0)then
   write(nrep,*)
   write(nrep,*)'==================='
   write(nrep,*)'PURE WATER ANALYSIS'
   write(nrep,*)'==================='
   write(nrep,*)
   write(nrep,*)'Number of waters, water sites',nmolwat,nwatsites
   write(nrep,*)'Box side length (NVE/T)',cube
   if(ktcfpureW==1)then
      write(nrep,*)'Calculate pure water orientational tcf = YES' 
      write(nrep,*)'   Water orientational tcf delay-time (ps)',mtdelay
      write(nrep,*)'   sample freq. (steps)',m_time
   endif
    if(ktcfhbW==1)then
      write(nrep,*)'Calculate pure water orientational tcf = YES' 
      write(nrep,*)'   Water orientational tcf delay-time (ps)',mtdelay
      write(nrep,*)'   sample freq. (steps)',m_time
      write(nrep,*)'   otcf are deconvoluted using HB definition'
   endif
   if(khbpureW==1)then
      write(nrep,*)'Calculate Pure Water HB dynamics = YES' 
      write(nrep,*)'   Water HB dynamics delay-time (ps)',hbtframe
      write(nrep,*)'   sample freq. (steps)',hbs_samp
   endif   
   if(ktetpureW==1)then
      write(nrep,*)'Calculate pure water tetrahedrality = YES' 
      write(nrep,*)'   sample freq. (steps)',tetsampw
   endif
   if(kVoronoipureW==1)then
      write(nrep,*)'Calculate pure water Voronoi = YES' 
      write(nrep,*)'   sample freq. (steps)',Voronsampw
   endif
   if(kdifpureW==1)then
      write(nrep,*)'Calculate pure water diffusion = YES'
      write(nrep,*)'   delay-time (ps)',msd_delay
      write(nrep,*)'   sample freq. (steps)',difsampw
   endif
   if(klsipureW==1)then
      write(nrep,*)'Calculate pure water LSI tcfs = YES'
      write(nrep,*)'   delay-time (ps)', lsitframe
      write(nrep,*)'   sample freq. (steps)',lsi_samp
      write(nrep,*)'   number of environments',Nenv_lsi
   endif
   
endif

write(nrep,*)
write(nrep,*)
write(nrep,*)
!
write(*,*)
write(*,*)'Finished reading input solv_inp/solv_inp.dat'
write(*,*)

!Solute and Solvent System
if(nsyst==1)then
   open(nsolprm,file='solv_inp/sol_sites.dat',status='old',action ='read')
   call read_sol_parm(nsolprm,nsoltypes,Zattype,mattype,niontypes,niontype,miontype,Ziontype,nchains,&
                      nchain_nat,nchain_lab,nb_res_ch,chain_label_res,natHSh,natsol_id,oatsol_hs,ratsol_hs,nb_res)
   write(nrep,*)
   write(nrep,*)'======================='
   write(nrep,*)'SOLUTE INPUT PARAMETERS'
   write(nrep,*)'======================='
   write(nrep,*)   
   write(nrep,*)'Number of solute atom types & mass'
   write(nrep,*)nsoltypes
   do i = 1,nsoltypes
      write(nrep,*)Zattype(i),mattype(i)
   end do
   write(nrep,*)
   write(nrep,*)'Number of ion types & mass'
   write(nrep,*)niontypes
   do i = 1,niontypes
      write(nrep,*)Ziontype(i),niontype(i),miontype(i)
   end do
   write(nrep,*)
!protein residues defined in the input   
   write(*,*)'Number of residues:',nb_res
   write(*,'(a72,i6)')'Number of solute atoms to analyze the respective solvation environment:',natHSh
!   write(*,*)natHSh
   write(*,*)
   write(nrep,*)'Number of residues:',nb_res
   write(nrep,'(a72,i6)')'Number of solute atoms to analyze the respective solvation environment:',natHSh
!   write(nrep,*)natHSh
   write(nrep,*)
endif

!Redefine the MD time-step to the step between frames in input
dt = dt*dble(dfr)                                            !Attention!!! dt becomes the time between frames and dfr will no longer be used 

if(nsyst==0)then
   write(nrep,*)
   write(nrep,*) '================='
   write(nrep,*) 'PURE WATER SYSTEM'
   write(nrep,*) '================='
   write(nrep,*)
   write(nrep,*) 'No. of frames in trj file = ',nstep
   write(nrep,*) 'No. of time-steps between frames = ',dfr
   write(nrep,*) 'time between frames in fs = ',dt
   write(nrep,*) 'No. of molecules of water = ',nmolwat
   natms = nmolwat*nwatsites 
   write(nrep,*) 'No. of atoms = ',natms
elseif(nsyst==1)then  
   write(nrep,*)
   write(nrep,*) '========================='
   write(nrep,*) 'SOLUTE AND SOLVENT SYSTEM'
   write(nrep,*) '========================='
   write(nrep,*)
   write(nrep,*) 'No. of frames in trj file = ',nstep
   write(nrep,*) 'No. of time-steps between frames = ',dfr
   write(nrep,*) 'time between frames in fs = ',dt
   write(nrep,*) 'No. of molecules of solute = ',nmolsol
   write(nrep,*) 'No. of molecules of water = ',nmolwat
   write(nrep,*) 'No. of ions = ',nions
   natms = nmolsol*natmsol + nmolwat*nwatsites + nions
   write(nrep,*) 'No. of atoms = ',natms
endif

!NG allocate(x(natms),y(natms),z(natms))
allocate(resindex(natms),resname(natms),segname(natms),atomname(natms),atomname_d(natms),atomtype(natms),chrg(natms))

! MD trajectory and system topology

!trajectory
if(inputformat.eq.'GROMACS')then
!  open (ntrjunf,file='solv_inp/unftrj',access = 'sequential',status='unknown',form = 'unformatted') 
!  open (ntrjunf,file='solv_inp/trj_unfolded',access = 'sequential',status='unknown',form = 'unformatted') 
!  open (ncell,file='solv_out/box_out.dat',status='unknown',action='write') 
!  call xtc2fortran_trj(ntrjunf,ncell,natms,nstep,trjfile)
!NG   call f77_molfile_init  
!topology   
   if(nsyst==1)then
      open(ngroitp,file= 'solv_inp/'//topfile,status='old',action='read')
      if(ksLJ==1)then
         open(nLJitp,file= 'solv_inp/'//LJfile,status='old',action='read')
      endif
      if(nb_res > 0)then
      call top_gro(ngroitp,ffname,nLJitp,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,niontypes,niontype,nb_res,nb_res_ch,nchains,segname,resindex,resname,&
                   chrg,atomname,atomname_d,atomtype,nsoltypes,Zattype,mattype,ZMBIO,Ziontype,chain_label_res,Natresid,name_res,numb_res,ksLJ,nattp,ncrossLJ,&
                   C6kjmol,C12kjmol,C6_IJ,C12_IJ,C6_IW,C12_IW)
      else
      call top_gro_sol(ngroitp,ffname,nLJitp,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,niontypes,niontype,nb_res,nb_res_ch,nchains,segname,resindex,resname,&
                       chrg,atomname,atomtype,nsoltypes,Zattype,mattype,ZMBIO,Ziontype,chain_label_res,Natresid,name_res,numb_res,ksLJ,nattp,ncrossLJ,&
                       C6kjmol,C12kjmol,C6_IJ,C12_IJ,C6_IW,C12_IW)             
      endif                   
      write(nrep,*)'List of solute atoms to study solvation and cut-offs (Ang):'
      write(*,*)'List of solute atoms to study solvation and cut-offs (Ang):'
      do i = 1,natHSh
         k = natsol_id(i)
         write(nrep,'(i6,1x,a5,a5,a1,F7.2,a1,F7.2,a1)')natsol_id(i),atomname(k),resname(k),'[',oatsol_hs(i),'-',ratsol_hs(i),']'
         write(*,'(i6,1x,a5,a5,a1,F7.2,a1,F7.2,a1)')natsol_id(i),atomname(k),resname(k),'[',oatsol_hs(i),'-',ratsol_hs(i),']'
      end do
      write(*,*)
      write(nrep,*)             
   endif             
elseif(inputformat.eq.'TINKER '.or.inputformat.eq.'BOMD   ')then
!formatted trajectory
   open(ntrjunf,file='solv_inp/'//trjfile,access ='sequential',status='old',action='read')
!if(nens==1) open(nbox,file='solv_inp/tinkbox.ang',access ='sequential',status='old',action='read')          !NPT ensemble     
if(nens==1) open(nbox,file='solv_inp/'//boxfile,access ='sequential',status='old',action='read')     
else
   write(*,*)'Input files should be in TINKER, GROMACS or BOMD format.'
   stop
endif

!Input warnings
!system 0: pure water; 1 solution
if(nsyst==1.and.khbpureW==1)then
   write(*,*)'input error - system solution/pure water HB dyn.'
   stop
elseif(khbpureW==1.and.kwtcfcrsh==1)then
   write(*,*)'input error - pure water HB dyn./solvent-water orient. dyn.'
   stop
elseif(khbpureW==0.and.khbtcf==1)then
   write(*,*)'input error - pure water HB dyn./solvent-water HB dyn.'
   stop   
endif   
   
   
!START CALCULATIONS   
! Solutions Routines

! biological-bulk water histograms
if(nsyst==1.and.kbiowater==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting biological water analysis routine ...'
   write(*,*)'=============================================='
   write(*,*)
   call  wat_bio_histo(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                  biowsamp,biowRadius,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO,WATMODEL)  
   write(*,*)
   write(*,*)'Normal termination of biological water analysis routine ...'
endif

!ORIENTATIONAL TIME-CORRELATION FUNCTION of SOLVENT-WATER

if(nsyst==1.and.kwtcfcrsh==1)then 
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting reor. tcf radial analysis routine ...'
   write(*,*)'=============================================='
   write(*,*)
   call wat_env_orient_tcfs(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nions,nwatsites,atomname,&
                            mtdelay,m_time,rwtcfbulk,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO,&
                            resindex,resname,atomtype,chrg,nb_res,chain_label_res)
   write(*,*)
   write(*,*)'Normal termination of reor. tcf radial analysis routine ...'
endif

!HB TIME-CORRELATION FUNCTIONS - CONTINUOUS AND INTERMITTENT

if(nsyst==1.and.khbtcf==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting HB tcf radial analysis routine ...'
   write(*,*)'==========================================='
   write(*,*)
!   call wat_HB_tcfs(trjfile,ntrjunf,nbox,dt,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,hbtdelay,hbs_time,&
!                    nsoltypes,Zattype,mattype,natsol_id,oatsol_hs,ratsol_hs)   
   write(*,*)
   write(*,*)'Normal termination of HB tcfs analysis routine ...'
endif

!TETRAHEDRALITY of the SOLVENT-WATER - Tetrahedrality calculated excluding solute heavy atoms

if(nsyst==1.and.kwtetra==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting tetrahedrality analysis routine ...'
   write(*,*)'============================================'
   write(*,*)
   call wat_tetra(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                  nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,rwtetbulk,tetsampw,ZMBIO,WATMODEL)  
   write(*,*)
   write(*,*)'Normal termination of tetrahedrality analysis routine ...'
endif

!TETRAHEDRALITY of the SOLVENT-WATER - Tetrahedrality calculated including solute heavy atoms

if(nsyst==1.and.kwsoltetra==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting water-solute tetrahedrality analysis routine ...'
   write(*,*)'========================================================='
   write(*,*)
   call wat_tetra_sol(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                      nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,rwsoltetbulk,tetsampws,ZMBIO)  
   write(*,*)
   write(*,*)'Normal termination of water-solute tetrahedrality analysis routine ...'
endif

!LSI of the SOLVENT-WATER - LOCAL STRUCTURE INDEX

if(nsyst==1.and.KWLSI==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting water-solute LSI analysis routine ...'
   write(*,*)'=============================================='
   write(*,*)
   call wat_LSI(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,&
                nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,LSI_sample,rwlsibulk,ZMBIO)
   write(*,*)
   write(*,*)'Normal termination of water-solute LSI analysis routine ...'
endif

!PROTEINS CONTACT MAP calculation

if(nsyst==1.and.kscontact==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init                              !This was passed inside the routine to allow opening multiple trajectories for analysis
   write(*,*)
   write(*,*)'Starting protein contact map analysis routine ...'
   write(*,*)'================================================='
   write(*,*)
   call prot_contact(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,&
                      nsoltypes,intsf,rcont,atomtype,nb_res,chain_label_res,Natresid,name_res,numb_res,nres_cont_pair_t0,res_cont_pair_t0)                   
   write(*,*)
   write(*,*)'Normal termination of protein contact map analysis routine ...'
endif

!PROTEINS ELECTROSTATIC AND VAN DER WAALS POTENTIAL MAP calculation

!if(nsyst==1.and.kspotential==1.and.Coul_th/=0.and.LJ126_th/=0)then
if(nsyst==1.and.kspotential==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init                              !This was also passed inside the routine to allow opening multiple trajectories for analysis - do not commet here
   write(*,*)
   write(*,*)'Starting potential map analysis routine ...'
   write(*,*)'==========================================='
   write(*,*)
   if(kumbrella_trj==0)then
      write(*,*)
      write(*,*)'WARNING!!! No umbrella trjs  - averages over trjs will be calculated'
      write(*,*)
!VERSION_12      call potent_map(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,Zattype,mattype,ZMBIO,&
!VERSION_12                     nsoltypes,ncontsamp,kumbrella_trj,kntrj,Coul_th,LJ126_th,atomtype,chrg,nb_res,chain_label_res,Natresid,name_res,numb_res,&
!VERSION_12                     C6kjmol,C12kjmol,nres_cont_pair_t0,res_cont_pair_t0,C6_IJ,C12_IJ)  
      
      call potent_map(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,atomname_d,Zattype,mattype,ZMBIO,&
                     nsoltypes,kumbrella_trj,ncontsamp,kntrj,Coul_th_min,Coul_th_max,LJ126_th_min,LJ126_th_max,atomtype,chrg,nb_res,chain_label_res,&
                     Natresid,name_res,numb_res,kMIC,intsf_,rcont_,C6kjmol,C12kjmol,C6_IJ,C12_IJ)  
                     
   elseif(kumbrella_trj==1)then
!Analyse different trjs without averaging as each trj corresponds to a different distance between proteins - umbrella windows
      write(*,*)
      write(*,*)'WARNING!!! umbrella trjs  - no averages over trjs will be calculated'
      write(*,*)
      call pot_map_umbrella(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,Zattype,mattype,ZMBIO,&
                            nsoltypes,kumbrella_trj,ncontsamp,kntrj,Coul_th_min,Coul_th_max,atomtype,chrg,nb_res,chain_label_res,&
                            Natresid,name_res,numb_res,kMIC,intsf_,rcont_,C6kjmol,C12kjmol,C6_IJ,C12_IJ) 
   endif
   write(*,*)
   write(*,*)'Normal termination of potential map analysis routine ...'
endif

!PROTEIN-WATER and WATER-WATER POTENTIAL ENERGY - NO EWALD SUM

if(nsyst==1.and.ksprotwat==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init                              !This was also passed inside the routine to allow opening multiple trajectories for analysis - do not commet here
   write(*,*)
   write(*,*)'Starting protein-water analysis routine ...'
   write(*,*)'==========================================='
   write(*,*)
   call prot_water_umbrella(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,Zattype,mattype,ZMBIO,&
                            nsoltypes,nions,mcontsamp,atomtype,chrg,nb_res,chain_label_res,Natresid,name_res,numb_res,kMIC_,C6kjmol,C12kjmol,kprot_wat,kwat_wat,&
                            C6_IW,C12_IW)
   write(*,*)
   write(*,*)'Normal termination of protein-water analysis routine ...'                         
endif

!VMD trajectory

if(nsyst==1.and.krmvmd==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting VMD-xyz trajectory routine ...'
   write(*,*)'======================================='
   write(*,*)
   call vmd_trj(trjfile,ntrjunf,inputformat,nbox,nens,cube,natms,nstep,dt,nmolsol,natmsol,nmolwat,nwatsites,nions,atomname,vmd_samp,nsoltypes,&
                Zattype,mattype,ZMBIO,natHSh,natsol_id,oatsol_hs,ratsol_hs,resindex,resname,atomtype,chrg,nchains,nchain_nat,nchain_lab)
   write(*,*)
   write(*,*)'Finished writing VMD-xyz trajectory ...'
endif


!DIFFUSION of the SOLVENT-WATER - routine call wat_env_diffusion computes the msd using the same number of points at each delay time
!Warning!!! the trajectory (.xtc) input file is re-named - Last routine called [fix] - THIS HAS BEEN FIXED in VERSION 13
if(nsyst==1.and.kdifwsol==1)then 
   rewind(ntrjunf)
   rewind(nbox)
!Unfolf trajectory for mean square displacement calculation - use gromacs   
   if(unfold_trj==0)then
      write(*,*)'Start unfolding trj for MSD calculation'
      write(*,*)'---> required input files: md_trj.xtc and md.gro'
      write(*,*)
      write(nrep,*)'Start unfolding trj for MSD calculation'
      write(nrep,*)
      call SYSTEM('echo -e "0\n" | gmx_d trjconv -f md_trj.xtc -s md.gro -o trj_unfold.xtc -pbc nojump')
   endif
!   write(trjfile,'(a)')'trj_unfold.xtc'
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting diffusion analysis routine ...'
   write(*,*)'======================================='
   write(*,*)
   call wat_env_diffusion(trjfile,ntrjunf,unfoldtrjfile,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nions,nwatsites,atomname,&
                          msd_delay,difsampw,difbulk,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO)
   write(*,*) 
   write(*,*)'Normal termination of diffusion analysis routine ...'
endif


!DIFFUSION of the SOLVENT-WATER - routine call wat_env_diff_stat computes the msd using the maximum number of origins - higher statistics at lower delay times
!THIS ROUTINE is NOW OUTDATED and HAS A BUG - IT USES THE UNFOLDED TRJ TO DETERMINE SOLVATION WATERS USING THE MIC
!Warning!!! the trajectory (.xtc) input file is re-named - Last routine called [fix]
!if(nsyst==1.and.kdifwsol==1)then 
!   rewind(ntrjunf)
!   rewind(nbox)
!Unfolf trajectory for mean square displacement calculation - use gromacs   
!NG   if(unfold_trj==0)then
!NG      write(*,*)'Start unfolding trj for MSD calculation'
!NG      write(*,*)'---> required input files: md_trj.xtc and md.gro'
!NG      write(*,*)
!NG      write(nrep,*)'Start unfolding trj for MSD calculation'
!NG      write(nrep,*)
!NG      call SYSTEM('echo -e "0\n" | gmx_d trjconv -f md_trj.xtc -s md.gro -o trj_unfold.xtc -pbc nojump')
!NG   endif
!NG   write(trjfile,'(a)')'trj_unfold.xtc'
!   call f77_molfile_init
!   write(*,*)
!   write(*,*)'Starting diffusion analysis routine ...'
!   write(*,*)'======================================='
!   write(*,*)
!   call wat_env_diff_stat(trjfile,ntrjunf,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nions,nwatsites,atomname,&
!                          msd_delay,difsampw,difbulk,nsoltypes,Zattype,mattype,natHSh,natsol_id,oatsol_hs,ratsol_hs,ZMBIO)                          
!   write(*,*) 
!   write(*,*)'Normal termination of diffusion analysis routine ...'
!endif

!*******************************************************PURE WATER*****************************************************************
   
! Pure Water Routines 

!Calculate orientational tcfs for pure water

if(nsyst==0.and.ktcfpureW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting reor. tcf radial analysis routine ...'
   write(*,*)'=============================================='
   write(*,*)
   if(Nenv==1)call pure_wat_orient_tcfs(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay,m_time)
! pure_wat_LDL_HDL_otcf can also be used for Nenv==1 - pure_wat_orient_tcfs can be commented and always call pure_wat_LDL_HDL_otcf 
!ok if(Nenv>=1)call pure_wat_LDL_HDL_otcf(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay,m_time,Nenv)
   if(Nenv>1)call pure_wat_LDL_HDL_otcf(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay,m_time,Nenv,kLDL_HDL,&
                                        LSI_min_1,LSI_max_1,LSI_min_2,LSI_max_2) 
   write(*,*)
   write(*,*)'Normal termination of orientational tcfs analysis routine ...'
endif

!Calculate orientational tcfs for pure water - deconvolute tcfs based on HB length
if(nsyst==0.and.ktcfhbW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting reor. tcf radial analysis routine [HB deconvolution] ...'
   write(*,*)'================================================================='
   write(*,*)
   call pure_wat_LDL_HDL_otcf_HB(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,mtdelay,m_time,Nenv_hb,kLDL_HDL_) 
   write(*,*)
   write(*,*)'Normal termination of orientational tcfs analysis routine ...'
endif


!HB TIME-CORRELATION FUNCTIONS OF PURE WATER - CONTINUOUS AND INTERMITTENT

if(nsyst==0.and.khbpureW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting HB tcf analysis routine ...'
   write(*,*)'===================================='
   write(*,*)
   if(Nenv_hb_HB<=1)call pure_wat_HB_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,hbtframe,hbs_samp)
   if(Nenv_hb_HB>1) call pure_wat_HB_LDL_HDL_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,hbtframe,hbs_samp,Nenv_hb_HB,&
                                                  LSI_min_1_,LSI_max_1_,LSI_min_2_,LSI_max_2_)
   write(*,*)
   write(*,*)'Normal termination of HB tcfs analysis routine ...'
endif



!LSI TIME-CORRELATION FUNCTIONS OF PURE WATER - CONTINUOUS AND INTERMITTENT
if(nsyst==0.and.klsipureW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting LSI tcf analysis routine ...'
   write(*,*)'====================================='
   write(*,*)
   call pure_wat_LSI_tcfs(trjfile,npurew,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,lsitframe,lsi_samp,Nenv_lsi,&
                          LSI_min_I,LSI_max_I,LSI_min_II,LSI_max_II)
   write(*,*)
   write(*,*)'Normal termination of LSI tcfs analysis routine ...'
endif



!TETRAHEDRALITY of PURE WATER 

if(nsyst==0.and.ktetpureW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting tetrahedrality analysis routine ...'
   write(*,*)'============================================'
   write(*,*)
   if(Nenvt == 1)call pure_wat_tetra(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,tetsampw,QMMM_pureW) 
   if(Nenvt > 1) call pure_wat_tetra_Env(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,tetsampw,R5Tet) 
   write(*,*)
   write(*,*)'Normal termination of tetrahedrality analysis routine ...'
endif


!VORONOI of PURE WATER 

if(nsyst==0.and.kVoronoipureW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting Voronoi analysis routine ...'
   write(*,*)'====================================='
   write(*,*)
   call pure_wat_Voronoi(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,Voronsampw,rcut) 
!not implemented   if(Nenvt > 1) call pure_wat_Voronoi_Env(trjfile,ntrjunf,inputformat,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,Voronsampw,R5Voron) 
   write(*,*)
   write(*,*)'Normal termination of Voronoi analysis routine ...'
endif



!MSD and Einstein diffusion of PURE WATER 
!Warning!!! the trajectory (.xtc) input file is re-named - Last routine called [fix]
if(nsyst==0.and.kdifpureW==1)then
   rewind(ntrjunf)
   rewind(nbox)
   if(unfold_trj==0)then
      write(*,*)'Start unfolding trj for MSD calculation'
      write(*,*)'---> required input files: md_trj.xtc and md.gro'
      write(*,*)
      write(nrep,*)'Start unfolding trj for MSD calculation'
      write(nrep,*)
      call SYSTEM('echo -e "0\n" | gmx_d trjconv -f md_trj.xtc -s md.gro -o trj_unfold.xtc -pbc nojump')
   endif   
   write(trjfile,'(a)')'trj_unfold.xtc'
!Unfolf trajectory for mean square displacement calculation - use gromacs
   call f77_molfile_init
   write(*,*)
   write(*,*)'Starting diffusion analysis routine ...'
   write(*,*)'======================================='
   write(*,*)
   if(Nenvrm==1)call pure_wat_diffus_msd(trjfile,ntrjunf,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,msd_delay,difsampw)
!Gaia and Lars   if(Nenvrm>1)call pure_wat_LDL_HDL_msd(trjfile,npurew,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,msd_delay,difsampw,Nenvrm)
   if(Nenvrm>1)call pure_wat_diffus_msd_Env(trjfile,ntrjunf,inputformat,nsyst,nbox,dt,nens,cube,natms,nstep,nmolwat,nwatsites,msd_delay,difsampw,R5Dif)
   write(*,*)
   write(*,*)'Normal termination of diffusion analysis routine ...'
endif





END PROGRAM Solvation
