MODULE water_solv

! contains subroutines:
! sol_com              - calculate the centre of mass of the solute
! sol_com_dist         - calculate the solute centre of mass
! prot_res_cont        - calculate protein-protein contact map
! tetrahedron          - calculate the tetrahedrality of water 
! sol_tetrahedron      - calculate the tetrahedrality of water - include the solute heavy atoms 
! tcf_solv_tau         - calculate the integral of any tcf
! msd_LSF_Einstein     - Fit the mean square displacement of water to compute the Einstein diffusion
! LSI_water            - calculate the LSI index for pure water
! LSI_water_fast       - calculate the LSI index for pure water - only sort the 10 nearest neighbors to speed up calculations
! water_water_En       - calculate the interaction energy between a water molecule and its five nearest neighbors
! VORON3               - calculate voronoi polyhedra
! SORT                 - voronoi polyhedra
! WORK                 - voronoi polyhedra

implicit none

 contains

SUBROUTINE sol_com(natms,nmolsol,natmsol,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
!Calculate solute centre of mass for bulk water tcfs and tetrahedrality calculation
!The centre of mass of the associated solute molecules is calculated - not the com of the individual solute molecules
    integer,intent(in)                                        :: natms,nmolsol,natmsol
    integer,intent(in)                                        :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in)          :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)              :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    real(kind=8),dimension(natmsol),intent(in)                :: ZMBIO       ! Solute aomic masses
!NG    real(kind=8),intent(in)                                   :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                                   :: x(natms),y(natms),z(natms)
    real(kind=8),intent(out)                                  :: cmxsol,cmysol,cmzsol
!Local    
    real(kind=8)                                              :: ZSOLMASS,CMXS,CMYS,CMZS
    integer                                                   :: i

 ZSOLMASS = 0.0D0
 CMXS=0.0D0
 CMYS=0.0D0
 CMZS=0.0D0
do i=1,nmolsol*natmsol                !assumes solute agregation
   CMXS     = CMXS+x(i)*ZMBIO(i)
   CMYS     = CMYS+y(i)*ZMBIO(i)
   CMZS     = CMZS+z(i)*ZMBIO(i)
   ZSOLMASS = ZSOLMASS + ZMBIO(i)
end do
 cmxsol = CMXS/ZSOLMASS
 cmysol = CMYS/ZSOLMASS
 cmzsol = CMZS/ZSOLMASS

return

END SUBROUTINE sol_com

SUBROUTINE sol_com_dist(natms,natmsol,nstart,nend,nsoltypes,Zattype,mattype,ZMBIO,x,y,z,cmxsol,cmysol,cmzsol)
!Calculate solute centre of mass 
    integer,intent(in)                                        :: natms,natmsol
    integer,intent(in)                                        :: nstart,nend
    integer,intent(in)                                        :: nsoltypes   ! Number of solute distinct atomic species (C, H, O, N etc)
    character(len=2),dimension(nsoltypes),intent(in)          :: Zattype     ! Atomic symbol of the solute atom types (C, H, O, N etc)
    real(kind=8),dimension(nsoltypes),intent(in)              :: mattype     ! Mass of atoms of type (C, H, O, N etc)
    real(kind=8),dimension(natmsol),intent(in)                :: ZMBIO       ! Solute aomic masses
!NG    real(kind=8),intent(in)                                   :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                                   :: x(natms),y(natms),z(natms)
    real(kind=8),intent(out)                                  :: cmxsol,cmysol,cmzsol
!Local    
    real(kind=8)                                              :: ZSOLMASS,CMXS,CMYS,CMZS
    integer                                                   :: i

 ZSOLMASS = 0.0D0
 CMXS=0.0D0
 CMYS=0.0D0
 CMZS=0.0D0
 
do i=nstart,nend            
   CMXS     = CMXS+x(i)*ZMBIO(i)
   CMYS     = CMYS+y(i)*ZMBIO(i)
   CMZS     = CMZS+z(i)*ZMBIO(i)
   ZSOLMASS = ZSOLMASS + ZMBIO(i)
end do
 cmxsol = CMXS/ZSOLMASS
 cmysol = CMYS/ZSOLMASS
 cmzsol = CMZS/ZSOLMASS

return

END SUBROUTINE sol_com_dist



SUBROUTINE prot_res_cont(trjfile,ninput,inputformat,nbox,dt,nens,cube,natms,nstep,nmolsol,natmsol,nmolwat,nwatsites,atomname,&
                         nsoltypes,kumbrella_trj,kMIC,intsf_,rcont_,atomtype,nb_res,chain_label_res,Natresid,name_res,numb_res,nres_cont_pair_t0,res_cont_pair_t0)  

! Calculate protein-protein contact map 
! Finds all residue-residue contacts during a trajectory based on a treshold distance "rcont_"; sample freq. is "intsf_" steps
! If rcont_ is small (i.e., not all pairs shall be considered) then intsf_ should be be set to a small number (e.g., 1) to increase the number of
! contacts that existed at least for a single time-step within rcont_ 

    integer,intent(in)                               :: ninput,nbox
    character(len=7),intent(in)                      :: inputformat          ! type of MD input: TINKER, GROMACS, BOMD
    real(kind=8),intent(in)                          :: dt                   ! time between frames (fs)
    integer,intent(in)                               :: nens                 ! ensemble: 0: NVE, NVT; 1: NPT
    real(kind=8),intent(in)                          :: cube(3)              ! box side length NVE and NVT
    integer,intent(in)                               :: natms,nstep,nmolsol,natmsol,nmolwat,nwatsites
    character(len=4),dimension(natms)                :: atomname    
    integer,intent(in)                               :: nsoltypes               ! Number of solute distinct atomic species (C, H, O, N etc)
    integer,intent(in)                               :: kumbrella_trj           ! 0: Trjs are not from umbrella sampling 1: Trjs are from umbrella sampling
    integer,intent(in)                               :: kMIC                    ! apply Minimum Image Convention or not to proteins
    integer,intent(in)                               :: intsf_                   ! Protein contact map sampling frequency
    real(kind=8),intent(in)                          :: rcont_                   ! distance for contact map - residues at a r < rcont_ are accounted
    character(len=4),dimension(natms),intent(in)     :: atomtype 
    integer,intent(in)                               :: nb_res                 ! Number of residues [GROMCAS]
    character(len=1),dimension(nb_res),intent(in)    :: chain_label_res        ! Chain label (e.g., A, B, C etc) internal code use
    integer,dimension(nb_res),intent(in)             :: Natresid               ! Number of atoms of each residue
    character(len=4),dimension(nb_res),intent(in)    :: name_res
    integer,dimension(nb_res),intent(in)             :: numb_res               ! Number of the residue: from 1 to the total number of residues of each chain

!Version_12    integer,dimension(:,:),allocatable,intent(out)   :: res_cont_pair_t0
    integer,dimension(nb_res,nb_res),intent(out)     :: res_cont_pair_t0
    integer,intent(out)                              :: nres_cont_pair_t0
    
!Local    
    real(kind=4),dimension(:),allocatable            :: x,y,z
    real(kind=8)                                     :: cell(3)  
    integer,dimension(natms)                         :: natmid                 !atomic index (solute and solvent index)
    integer      :: n0, n1
    integer      :: i, j, k
    integer      :: ic1,ic2,iats1,iats2,k1,k2,it,jt,jats1,jats2,kj1,kj2                                         !chain indexes for contact map
    real(kind=8)                                     :: dx,dy,dz,dr,dr2 
    integer                                          :: natcheck
    logical                                          :: filex
    
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

!xtc2fortran
allocate(xyz(natms*3))
allocate(x(natms),y(natms),z(natms))
!end xtc2fortran

n0 = 80
n1 = 90

if(kumbrella_trj==0)then
   open(n0,file='Prot_potent_out/log_res_contact.dat',status='unknown',action='write')
   open(n1,file='Prot_potent_out/prot_res_map.dat',status='unknown',action='write')
elseif(kumbrella_trj==1)then   
   open(n0,file='Prot_potent_umbrella_out/log_res_contact.dat',status='unknown',action='write')
   open(n1,file='Prot_potent_umbrella_out/prot_res_map.dat',status='unknown',action='write')
endif

!Print headers

write(n1,*)'contact pair (FIRST) : RES_in_system; atom; RES_name; RES_in_chain; chain'
write(n1,*)'contact pair (SECOND): RES_in_system; atom; RES_name; RES_in_chain; chain'
write(n1,*)'total nb of residues = ',nb_res,' total nb of protein atoms = ',natmsol
write(n1,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(n1,*)

write(*,*)
write(*,*)'Starting protein contact map calculation...'
write(*,*)'Routine prot_res_cont'
if(kMIC==0)write(*,*)'MIC will not be applied - no jump protein'
if(kMIC==1)write(*,*)'MIC will be applied'
write(*,*)'Results are printed to Prot_potent_out'
write(*,*)'traj sampling freq for contacts (steps) = ',intsf_
write(*,*)'max. residue-residue contact dist. (Ang) = ',rcont_
write(*,*)'total number of time-steps = ',nstep
write(*,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(*,*)
write(n0,*)
write(n0,*)'Starting protein contact map calculation...'
write(n0,*)'Routine prot_res_cont'
if(kMIC==0)write(n0,*)'MIC will not be applied - no jump protein'
if(kMIC==1)write(n0,*)'MIC will be applied'
write(n0,*)'Results are printed to Prot_potent_out'
write(n0,*)'traj sampling freq for contact map (steps) = ',intsf_
write(n0,*)'max. residue-residue contact dist. (Ang) = ',rcont_
write(n0,*)'total number of time-steps =',nstep
write(n0,*)'Maximum nb of contact pairs possible = ',(nb_res/2)**2
write(n0,*)
             
kj1 = 0
kj2 = 0

res_cont_pair_t0 = 0
nres_cont_pair_t0 = 0

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

!Calculate protein contacts between chains ABCD and EFGH EVERY "intsf_" steps
!The contacts will change along time and these will be accumulated in the array res_cont_pair_t0(ic1,ic2)...

   if(j==1.or.mod(j,intsf_)==0)then                     !sample new contacts every intsf_ steps 
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
                  if(kMIC==1)then
                     dx = dx - dnint(dx/cell(1))*cell(1)
                     dy = dy - dnint(dy/cell(2))*cell(2)
                     dz = dz - dnint(dz/cell(3))*cell(3)
                  endif
                  dr2 = dx**2 + dy**2 + dz**2
                  dr  = dsqrt(dr2)
!check                  write(*,*)dr
                  IF((dr.LE.rcont_).and.(res_cont_pair_t0(ic1,ic2)==0))THEN             
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

deallocate(xyz)
deallocate(x,y,z)
!deallocate(res_cont_pair_t0)

return

END SUBROUTINE prot_res_cont



SUBROUTINE tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                       SG,CHHASG,RO4,RO_nb,RH2,RHO1,NWTDNB,RNTDWH)
!TETRAHEDRALITY and DISTANCE PARAMETERS CALCULATION
    integer,intent(in)                         :: i
    integer,intent(in)                         :: natms,nmolsol,natmsol,nmolwat,nwatsites,nions
!NG    real(kind=8),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=8),intent(in)                    :: cell(3) 
    
    real(kind=8),intent(out)                   :: SG,CHHASG,RO4
    real(kind=8),intent(out)                   :: RO_nb(5)
    real(kind=8),intent(out)                   :: RH2,RHO1
    INTEGER                                    :: NWTDNB(4)
    REAL(KIND=4)                               :: RNTDWH(2)
    INTEGER                                    :: NTDWH(2)
!Local
    integer                                    :: jw,kw
    integer                                    :: JT   
    integer                                    :: IE,JE
    integer                                    :: KT,KI,KO,KH,KN,KHT
    real(kind=8)                               :: RADMAX,RADMIN
    real(kind=8)                               :: dx,dy,dz,dr,dr2   
!    real(kind=8)                               :: NTDWO(4),NTD5O(5),NTDWOO(2),NTDWOH(2),NTDWHO(2)                  !Version 9 correction
    INTEGER                                    :: NTDWO(4),NTD5O(5),NTDWOO(2),NTDWOH(2),NTDWHO(2)
    real(kind=8)                               :: THDP
    real(kind=8)                               :: XT(4),YT(4),ZT(4)
    
RO_nb = 0.0d0;SG=0.0d0;CHHASG=0.0d0;RO4=0.0d0;RH2=0.0d0;RHO1=0.0d0

! I FIND THE FOUR NEAREST O ATOMS jw TO EACH OXYGEN i and CALCULATE the TETRAHEDRALITY

RADMAX  = 25.0D0
RADMIN  =  0.0D0
JT      = 1
DO WHILE(JT.LE.4) 
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
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
 
! FIND THE 4 NEAREST OXYGEN ATOMS jw TO THE OXYGEN ATOM i
     
         IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
            RADMAX = dr
! NTDWO STORES THE O ATOMS NEIGHBORS OF I
! JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR ; JT = 3 THIRD NEIGHBOR ; JT = 4 FOURTH NEIGHBOR
            NTDWO(JT) = jw
         ENDIF
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 25.0D0
END DO

! Version 9
NWTDNB = NTDWO

! I FIND THE TWO NEAREST H ATOMS TO EACH OXYGEN i

RADMAX  = 25.0D0
RADMIN  =  0.0D0
JT      = 1
DO WHILE(JT.LE.2) 
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
      IF(jw.NE.i)THEN
         do kw = 1,2
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND HYDROGEN ATOMS jw+1 and jw+2     
            dx = x(i) - x(jw+kw)             
            dy = y(i) - y(jw+kw)
            dz = z(i) - z(jw+kw)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)
 
! FIND THE 2 NEAREST HYDROGEN ATOMS jw TO THE OXYGEN ATOM i
     
            IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
               RADMAX = dr
! NTDWO STORES THE O ATOMS NEIGHBORS OF I
! JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR 
               NTDWH(JT) = jw
               RNTDWH(JT) = dr
            ENDIF
         end do   
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 25.0D0
END DO

! End Version 9 

! LOOK FOR ERRORS
DO IE=1,3
   DO JE = IE+1,4
      if(NTDWO(IE).EQ.NTDWO(JE))then
         write(*,*)'Error - water 4 nearest neighbors ',i
         stop
      endif   
   END DO
END DO

! CALCULATE FOR EACH WATER MOLECULE THE TETRAHEDRAL PARAMETER
! CHECK FOR SIX ANGLES FOR EACH WATER MOLECULE BASED UPON DOT PRODUCTS OF THE FOUR UNIT VECTORS
! DEFINING THE VERTICES OF A TETRAHEDRON WHERE EACH OXYGEN I IS AT THE TETRAHEDRON CENTRE
! CONSIDER FOUR OXYGEN ATOMS IN THE VERTICES OF THE TETRAHEDRON

 RO4    = 0.0D0                       !distance of the four nearest neighbors
 SG     = 0.0D0                       !tetrahedrality
 CHHASG = 0.0D0                       !tetrahedrality

DO KT = 1,4
! oxygen ID            
   KI=NTDWO(KT)
   dx = x(i)-x(KI)
   dy = y(i)-y(KI)
   dz = z(i)-z(KI)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2)
   RO4 = RO4 + dr
! CALCULATE THE FOUR UNIT VECTORS
   XT(KT)=dx/dr
   YT(KT)=dy/dr
   ZT(KT)=dz/dr
END DO

RO4 = RO4/4.0D0

DO JT = 1,3
   DO KT = JT+1,4
! dot product
      THDP=XT(JT)*XT(KT)+YT(JT)*YT(KT)+ZT(JT)*ZT(KT)
      SG = SG + (THDP+1.0D0/3.0D0)**2.0
   END DO
END DO
! Chau and Hardwick normalization
 CHHASG=SG*3.0D0/32.0D0
! Errington and Debenedetti normalization
 SG=1.0D0-SG*3.0D0/8.0D0
 
! II FIND THE 5 CLOSEST OXYGEN ATOMS TO COMPUTE THE RO_nb DISTRIBUTIONS - DISTRIBUTION OF THE KNth NEAREST NEIGHBOR DISTANCE (KN = 1,2,3,4,5)
RADMAX  = 25.0D0
RADMIN  =  0.0D0
JT      = 1
DO WHILE(JT.LE.5)
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
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
!  NTD5O STORES THE O ATOMS NEIGHBORS OF I
!  JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR ; JT = 3 THIRD NEIGHBOR ; JT = 4 FOURTH NEIGHBOR ; JT = 5 FIFTH NEIGHBOR
            NTD5O(JT) = jw
         ENDIF
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 25.0D0
END DO

! COMPUTE THE RO1, RO2, RO3, RO4 AND RO5 DISTANCES
! DISTANCES OF THE FIRST, SECOND, THIRD, FOURTH AND FIFTH NEIGHBORS 

RO_nb = 0.0D0
do KN = 1,5
   KO = NTD5O(KN)
   dx = x(i)-x(KO)
   dy = y(i)-y(KO)
   dz = z(i)-z(KO)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2) 
   RO_nb(KN) = dr
end do   



! III FIND THE TWO NEAREST H ATOMS KH TO EACH OXYGEN i - THIS IS FOR THE CALCULATION OF RH2

RADMAX  = 25.0D0
RADMIN  =  0.0D0
JT      = 1
DO WHILE(JT.LE.2) 
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
      IF(jw.NE.i)THEN
         DO KH = 1,2
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND HYDROGEN ATOMS KH OF THE jw WATER MOLECULES     
            dx = x(i) - x(jw+KH)             
            dy = y(i) - y(jw+KH)
            dz = z(i) - z(jw+KH)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)        
! FIND THE NEAREST TWO HYDROGEN ATOMS KH TO THE OXYGEN ATOM i              
            IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
               RADMAX = dr
!  NTDWOO STORES THE O ATOMS jw NEAREST NEIGHBORS OF i                     
!  NTDWOH STORES THE RESPECTIVE H ATOM - KH = 1 or 2
!  JT = 1 FIRST H NEIGHBOR ; JT = 2 SECOND H NEIGHBOR
               NTDWOO(JT) = jw
               NTDWOH(JT) = KH
            ENDIF          
         END DO
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 25.0D0
END DO
         
RH2 = 0.0D0
DO KHT = 1,2
! OXYGEN ID
   KO = NTDWOO(KHT)
! HYDROGEN ID            
   KH = NTDWOH(KHT)
   dx = x(i) - x(KO+KH)             
   dy = y(i) - y(KO+KH)
   dz = z(i) - z(KO+KH)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2) 
   RH2   = RH2 + dr
END DO
RH2   = RH2/2.0D0



!IV FIND THE NEAREST O ATOM TO EACH HYDROGEN ATOM - THIS IS FOR THE CALCULATION OF RHO1

DO KH = 1,2
   RADMAX = 25.0D0
   RADMIN = 0.0D0
   JT     = 1
   DO WHILE(JT.LE.1) 
      do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
         IF(jw.NE.i)THEN
! COMPONENTS OF VECTOR DISTANCE BETWEEN THE H ATOMS OF i AND THE O ATOMS OF THE jw WATER MOLECULES
            dx = x(i+KH) - x(jw)             
            dy = y(i+KH) - y(jw)
            dz = z(i+KH) - z(jw)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)                          
            IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
               RADMAX = dr
! NTDWHO STORES THE CLOSEST O ATOM jw TO THE H ATOMS OF THE WATER i                     
               NTDWHO(KH) = jw
            ENDIF
         ENDIF
      end do
      JT = JT + 1     
!      RADMIN = RADMAX
!      RADMAX = 25.0D0
   END DO
END DO   
         
RHO1 = 0.0D0
do KH = 1,2
   KO = NTDWHO(KH)
   dx = x(i+KH) - x(KO)             
   dy = y(i+KH) - y(KO)
   dz = z(i+KH) - z(KO)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2)     
   RHO1 = RHO1 + dr
END DO         
RHO1 = RHO1/2.0D0

return
END SUBROUTINE tetrahedron


!version 31

SUBROUTINE water_water_En(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,WATMODEL,ENvdw_nb,ENcoul_nb,ENpot_nb)

!Calculate P_energy(water-water_n) for n = 1 to 5 for water molecules in the bulk and in the HSh - deconvolute molecules with 4MWN and L4WN      

    integer,intent(in)                         :: i
    integer,intent(in)                         :: natms,nmolsol,natmsol,nmolwat,nwatsites,nions
    character(len=7),intent(in)                :: WATMODEL
!NG    real(kind=8),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=8),intent(in)                    :: cell(3) 
    real(kind=8),intent(out)                   :: ENvdw_nb(5),ENcoul_nb(5),ENpot_nb(5)
    real(kind=8)                               :: AVSNO,PMTTV,ELCH
!version 31
!    PARAMETER(BOHR  = 0.5291772108D0)
    PARAMETER(AVSNO = 6.02214199D+23)
!c 1/(4*PI*EPS0) C-2JAng 8.9875D19      
!    PARAMETER(PMTTV=8.9875D19/BOHR)
    PARAMETER(PMTTV=8.9875D19)
!c electron (C)      
    PARAMETER(ELCH=1.602177D-019)
!end version 31    
    INTEGER                                    :: NWTDNB(4)
    REAL(KIND=4)                               :: RNTDWH(2)
    INTEGER                                    :: NTDWH(2)
!Local
    integer                                    :: jw,kw
    integer                                    :: JT   
    integer                                    :: IE,JE
    integer                                    :: KT,KI,KO,KH,KN,KHT
    real(kind=8)                               :: RADMAX,RADMIN
    real(kind=8)                               :: dx,dy,dz,dr,dr2   
    integer                                    :: NTD5O(5)
    real(kind=8)                               :: THDP
    real(kind=8)                               :: XT(4),YT(4),ZT(4)
!version 31    
    real(kind=8)                               :: r_oxy
    real(kind=8)                               :: CHOXM,CHHYD,ESPLJ,SIGMLJ,FCONVJ
    real(kind=8)                               :: CH_MH(3)
!end version 31
    
if(WATMODEL=='TIP4PEW')then
!check   write(*,*)'WATER MODEL = ',WATMODEL 
!TIP4P-EW PARAMETERS
   CHOXM = -1.04844D0*ELCH
   CHHYD = +0.52422D0*ELCH
! Value in joules [value in kJ/mol = 0.680946]
   ESPLJ = 1.13073720468687D-021
! Value in Ang 3.16435       
   SIGMLJ= 3.16435D0
!TIP4P-Ew charges H, H, M
   CH_MH(1) = CHHYD
   CH_MH(2) = CHHYD
   CH_MH(3) = CHOXM   
elseif(WATMODEL=='TIP4P05')then 
!TIP4P-2005 PARAMETERS
   CHOXM = -1.1128D0*ELCH
   CHHYD = +0.5564D0*ELCH
! Value in joules [value in kJ/mol = 0.7749]   
   ESPLJ = 1.28675169323919E-021
! Value in Ang 3.1589       
   SIGMLJ= 3.1589D0
!TIP4P-2005 charges H, H, M
   CH_MH(1) = CHHYD
   CH_MH(2) = CHHYD
   CH_MH(3) = CHOXM      
elseif(WATMODEL=='SPCE  ')then    
!SPCE PARAMETERS
   CHOXM = -0.8476*ELCH
   CHHYD = +0.4238*ELCH
! Value in joules [value in kJ/mol = 0.650]   
   ESPLJ = 1.07935036857075E-021
! Value in Ang 3.1589       
   SIGMLJ= 3.166
!TIP4P-2005 charges H, H, M
   CH_MH(1) = CHOXM
   CH_MH(2) = CHHYD
   CH_MH(3) = CHHYD        
endif  
!
 FCONVJ = 1.0D-3*AVSNO    
    
ENvdw_nb = 0.0d0;ENcoul_nb = 0.0d0;ENpot_nb = 0.0d0
 
! FIND THE 5 CLOSEST OXYGEN ATOMS TO COMPUTE THE Energy_nb DISTRIBUTIONS - DISTRIBUTION OF THE KNth NEAREST NEIGHBOR interaction energy (KN = 1,2,3,4,5)
! FIND THE 5 NEAREST O ATOMS jw TO EACH OXYGEN i
RADMAX  = 25.0D0
RADMIN  =  0.0D0
JT      = 1
DO WHILE(JT.LE.5)
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
      IF(jw.NE.i)THEN
!COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND OXYGEN ATOMS jw     
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
!  NTD5O STORES THE O ATOMS NEIGHBORS OF I
!  JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR ; JT = 3 THIRD NEIGHBOR ; JT = 4 FOURTH NEIGHBOR ; JT = 5 FIFTH NEIGHBOR
            NTD5O(JT) = jw
         ENDIF
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 25.0D0
END DO

! COMPUTE THE RO1, RO2, RO3, RO4 AND RO5 DISTANCES
! DISTANCES OF THE FIRST, SECOND, THIRD, FOURTH AND FIFTH NEIGHBORS 
if(nwatsites==4)then                 !TIP4P potentials
   do KN = 1,5                       !5 neighbors
      KO = NTD5O(KN)
!oxygen(i)-oxgen distance(jw) - van der Waals interactions
!compute the van der Waals interaction between i and jw = 1, 2, 3, 4, and 5 individually
      dx = x(i)-x(KO)
      dy = y(i)-y(KO)
      dz = z(i)-z(KO)
      dx = dx - dnint(dx/cell(1))*cell(1)
      dy = dy - dnint(dy/cell(2))*cell(2)
      dz = dz - dnint(dz/cell(3))*cell(3)
      dr2 = dx**2 + dy**2 + dz**2
      dr  = dsqrt(dr2) 
      r_oxy = dr
      ENvdw_nb(KN) = FCONVJ*4.0D0*ESPLJ*((SIGMLJ/dr)**12.0D0-(SIGMLJ/dr)**6.0D0)
!Hatom1(i), Hatom2(i) & Msite(i) interactions with Hatom1(jw) Hatom2(jw)and MS-site(jw) 
      do jw = 1,nwatsites-1 
         do kw = 1,nwatsites-1
            dx = x(i+jw)-x(KO+kw)
            dy = y(i+jw)-y(KO+kw)
            dz = z(i+jw)-z(KO+kw)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)
!calculate electrostatic interactions - sum 9 interactions for each water molecule i with the KNth neighbor         
            ENcoul_nb(KN) = ENcoul_nb(KN) + PMTTV*CH_MH(jw)*CH_MH(kw)/dr
!check         write(*,*)jw,kw,dr,PMTTV*CH_MH(jw)*CH_MH(kw)/dr
         end do 
      end do
   !calculate potential energy between water molecule i and each of the 5 nearest neighbors - energies are not summed   
      ENpot_nb(KN) = ENvdw_nb(KN) + FCONVJ*ENcoul_nb(KN)
   !check   write(*,*)r_oxy,ENvdw_nb(KN),ENcoul_nb(KN),ENpot_nb(KN)
   end do
elseif(nwatsites==3)then                               !SPC/E
   do KN = 1,5                       !5 neighbors
      KO = NTD5O(KN)
   !oxygen(i)-oxgen distance(jw) - van der Waals interactions
   !compute the van der Waals interaction between i and jw = 1, 2, 3, 4, and 5 individually
      dx = x(i)-x(KO)
      dy = y(i)-y(KO)
      dz = z(i)-z(KO)
      dx = dx - dnint(dx/cell(1))*cell(1)
      dy = dy - dnint(dy/cell(2))*cell(2)
      dz = dz - dnint(dz/cell(3))*cell(3)
      dr2 = dx**2 + dy**2 + dz**2
      dr  = dsqrt(dr2) 
      r_oxy = dr
      ENvdw_nb(KN) = FCONVJ*4.0D0*ESPLJ*((SIGMLJ/dr)**12.0D0-(SIGMLJ/dr)**6.0D0)
!Hatom1(i), Hatom2(i) & Msite(i) interactions with Hatom1(jw) Hatom2(jw)and MS-site(jw) 
      do jw = 0,nwatsites-1 
         do kw = 0,nwatsites-1
            dx = x(i+jw)-x(KO+kw)
            dy = y(i+jw)-y(KO+kw)
            dz = z(i+jw)-z(KO+kw)
            dx = dx - dnint(dx/cell(1))*cell(1)
            dy = dy - dnint(dy/cell(2))*cell(2)
            dz = dz - dnint(dz/cell(3))*cell(3)
            dr2 = dx**2 + dy**2 + dz**2
            dr  = dsqrt(dr2)
!calculate electrostatic interactions - sum 9 interactions for each water molecule i with the KNth neighbor         
            ENcoul_nb(KN) = ENcoul_nb(KN) + PMTTV*CH_MH(jw)*CH_MH(kw)/dr
!check         write(*,*)jw,kw,dr,PMTTV*CH_MH(jw)*CH_MH(kw)/dr
         end do 
      end do
!calculate potential energy between water molecule i and each of the 5 nearest neighbors - energies are not summed   
      ENpot_nb(KN) = ENvdw_nb(KN) + FCONVJ*ENcoul_nb(KN)
!check   write(*,*)r_oxy,ENvdw_nb(KN),ENcoul_nb(KN),ENpot_nb(KN)
   end do

endif

return

END SUBROUTINE water_water_En
      
!End version 31



SUBROUTINE LSI_water(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
    integer,intent(in)                         :: iw,io
    integer,intent(in)                         :: natms,nmolsol,natmsol,nmolwat,nwatsites,nions
!NG    real(kind=8),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=8),intent(in)                    :: cell(3)
    real(kind=8),intent(in)                    :: r_lsi
    real(kind=8),intent(out)                   :: LSI(nmolwat)
!Local    
    integer                                    :: jw, JT, kw
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    real(kind=8)                               :: RADMAX,RADMIN
    integer                                    :: NTDWO(nmolwat-1)
    real(kind=8)                               :: LSI_delta, LSI_delta_mean
    integer                                    :: n_ut                       !number of oxygens included in the calculation of the LSI
    real(kind=8)                               :: dr_u(nmolwat-1)            !ordered distances from iw (u is used in the paper by Shiratani and Sasai, 1996)
    

!check WRITE(*,*)'Entering LSI_water'
   
!ORDER OXYGEN ATOMS jw AROUND EACH OXYGEN iw and CALCULATE the LSI INDEX

RADMAX  = 500.0D0                     ! Angstroms
RADMIN  =  0.0D0                      !
JT      = 1
n_ut    = 0                           !Number of waters included in the LSI calculation
LSI_delta_mean = 0.0d0
LSI_delta      = 0.0d0
LSI            = 0.0d0

DO WHILE(JT.LE.nmolwat-1) 
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
      IF(jw.NE.iw)THEN
!check         WRITE(*,*)iw,jw
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM iw AND OXYGEN ATOMS jw     
         dx = x(iw) - x(jw)             
         dy = y(iw) - y(jw)
         dz = z(iw) - z(jw)
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)     
         IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
            RADMAX = dr
! NTDWO STORES THE O ATOMS NEIGHBORS OF I
! JT = JT NEIGHBOR with JT = 1,2,3,...,nmolwat-1
            NTDWO(JT) = jw                               !stores nmolwat-1 oxygen atoms id
         ENDIF
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 500.0D0
END DO

do kw = 1,nmolwat-1
   jw = NTDWO(kw)                   !jw is not 1, 2, 3, 4,... it is equivalent to ihb or jhb
!check   write(*,*)iw,jw
   dx = x(iw) - x(jw)             
   dy = y(iw) - y(jw)
   dz = z(iw) - z(jw)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2)
   dr_u(kw) = dr                     !oxygen-oxygen distances ordered from the shortest to the longest from water oxygen iw
!check   write(*,*)iw,jw,dr
end do  

!calculate delta_mean of LSI
do kw = 1,nmolwat-2                 !Each oxygen has a maximum of nmolwat - 1 neighbours; the -2 is because of the difference r_i+1 - ri
   if(dr_u(kw) < r_lsi)then
      n_ut = n_ut + 1
      LSI_delta_mean = LSI_delta_mean + (dr_u(kw+1) - dr_u(kw))
   endif
!check   if(iw==xxx)write(*,*)iw,kw,n_ut,dr_u(kw)
!  if (n_ut==0)write(*,*)'No coordination sphere' 
end do
LSI_delta_mean = LSI_delta_mean/dfloat(n_ut)
!check   write(*,*)io,n_ut,LSI_delta_mean

if(LSI_delta_mean < 0.0d0)then
   write(*,*)'Error!!! LSI_delta_mean < 0'
   stop
endif

!calculate LSI
n_ut = 0                            !Re-count oxygens
do kw = 1,nmolwat-2                 !Each oxygen has a maximum of nmolwat - 1 neighbours; the -2 is because of the difference r_i+1 - ri
   if(dr_u(kw) < r_lsi)then
      n_ut = n_ut + 1
      LSI_delta = LSI_delta + ((dr_u(kw+1) - dr_u(kw)) - LSI_delta_mean)**2.0d0
   endif
end do
LSI_delta = LSI_delta/dfloat(n_ut)
!check   write(*,*)io,n_ut,LSI_delta
LSI(io) = LSI_delta
    
return
END SUBROUTINE LSI_water


SUBROUTINE LSI_water_fast(iw,io,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,LSI,r_lsi)
!TO SPEED UP CALCULATION WE ONLY SORT THE 9 NEAREST WATER NEIGHBORS
!REPLACED nmolwat by NW_sort
    integer,intent(in)                         :: iw,io
    integer,intent(in)                         :: natms,nmolsol,natmsol,nmolwat,nwatsites,nions
!NG    real(kind=8),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                    :: x(natms),y(natms),z(natms)
    real(kind=8),intent(in)                    :: cell(3)
    real(kind=8),intent(in)                    :: r_lsi
    real(kind=8),intent(out)                   :: LSI(nmolwat)
!Local    
    integer                                    :: jw, JT, kw
    integer                                    :: NW_sort
    real(kind=8)                               :: dx,dy,dz,dr,dr2
    real(kind=8)                               :: RADMAX,RADMIN
    integer                                    :: NTDWO(nmolwat-1)
    real(kind=8)                               :: LSI_delta, LSI_delta_mean
    integer                                    :: n_ut                       !number of oxygens included in the calculation of the LSI
    real(kind=8)                               :: dr_u(nmolwat-1)            !ordered distances from iw (u is used in the paper by Shiratani and Sasai, 1996)
    

!check WRITE(*,*)'Entering LSI_water'
   
!ORDER OXYGEN ATOMS jw AROUND EACH OXYGEN iw and CALCULATE the LSI INDEX
NW_sort = 20
RADMAX  = 500.0D0                     ! Angstroms
RADMIN  =  0.0D0                      !
JT      = 1
n_ut    = 0                           !Number of waters included in the LSI calculation
LSI_delta_mean = 0.0d0
LSI_delta      = 0.0d0
LSI            = 0.0d0

!NG DO WHILE(JT.LE.nmolwat-1) 
DO WHILE(JT.LE.NW_sort-1)
   do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
      IF(jw.NE.iw)THEN
!check         WRITE(*,*)iw,jw
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM iw AND OXYGEN ATOMS jw     
         dx = x(iw) - x(jw)             
         dy = y(iw) - y(jw)
         dz = z(iw) - z(jw)
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)     
         IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
            RADMAX = dr
! NTDWO STORES THE O ATOMS NEIGHBORS OF I
! JT = JT NEIGHBOR with JT = 1,2,3,...,nmolwat-1
            NTDWO(JT) = jw                               !stores nmolwat-1 oxygen atoms id
         ENDIF
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 500.0D0
END DO

!NG   do kw = 1,nmolwat-1
do kw = 1,NW_sort-1
   jw = NTDWO(kw)                   !jw is not 1, 2, 3, 4,... it is equivalent to ihb or jhb
!check   write(*,*)iw,jw
   dx = x(iw) - x(jw)             
   dy = y(iw) - y(jw)
   dz = z(iw) - z(jw)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2)
   dr_u(kw) = dr                     !oxygen-oxygen distances ordered from the shortest to the longest from water oxygen iw
!check   write(*,*)iw,jw,dr
end do  

!calculate delta_mean of LSI
!NG   do kw = 1,nmolwat-2                 !Each oxygen has a maximum of nmolwat - 1 neighbours; the -2 is because of the difference r_i+1 - ri
do kw = 1,NW_sort-2
   if(dr_u(kw) < r_lsi)then
      n_ut = n_ut + 1
      LSI_delta_mean = LSI_delta_mean + (dr_u(kw+1) - dr_u(kw))
   endif
!check   if(iw==xxx)write(*,*)iw,kw,n_ut,dr_u(kw)
!  if (n_ut==0)write(*,*)'No coordination sphere' 
end do
LSI_delta_mean = LSI_delta_mean/dfloat(n_ut)
!check   write(*,*)io,n_ut,LSI_delta_mean

if(LSI_delta_mean < 0.0d0)then
   write(*,*)'Error!!! LSI_delta_mean < 0'
   stop
endif

!calculate LSI
n_ut = 0                            !Re-count oxygens
!NG   do kw = 1,nmolwat-2                 !Each oxygen has a maximum of nmolwat - 1 neighbours; the -2 is because of the difference r_i+1 - ri
do kw = 1,NW_sort-2
   if(dr_u(kw) < r_lsi)then
      n_ut = n_ut + 1
      LSI_delta = LSI_delta + ((dr_u(kw+1) - dr_u(kw)) - LSI_delta_mean)**2.0d0
   endif
end do
LSI_delta = LSI_delta/dfloat(n_ut)
!check   write(*,*)io,n_ut,LSI_delta
LSI(io) = LSI_delta
    
return
END SUBROUTINE LSI_water_fast




SUBROUTINE sol_tetrahedron(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,cell,&
                           SG,CHHASG,atomname)                           
!TETRAHEDRALITY CALCULATION - include solute heavy atoms as tetrahedron vertices
    integer,intent(in)                                 :: i
    integer,intent(in)                                 :: natms,nmolsol,natmsol,nmolwat,nwatsites,nions
!NG    real(kind=8),intent(in)                            :: x(natms),y(natms),z(natms)
    real(kind=4),intent(in)                            :: x(natms),y(natms),z(natms)
    real(kind=8),intent(in)                            :: cell(3) 
    character(len=4),dimension(natms),intent(in)       :: atomname    
!    character(len=1),dimension(natms),intent(in)       :: atom_symb    
    real(kind=8),intent(out)                           :: SG,CHHASG
!Local
    integer                                            :: jw
    integer                                            :: JT   
    integer                                            :: IE,JE
    integer                                            :: KT,KI,KO,KH,KN,KHT
    real(kind=8)                                       :: RADMAX,RADMIN
    real(kind=8)                                       :: dx,dy,dz,dr,dr2   
!    real(kind=8)                                       :: NTDWO(4)
    INTEGER                                            :: NTDWO(4)
    real(kind=8)                                       :: THDP
    real(kind=8)                                       :: XT(4),YT(4),ZT(4)
    character(len=1)                                   :: atom_symb
    
SG=0.0d0;CHHASG=0.0d0

! FIND THE FOUR NEAREST WATER OXYGENS, SOLUTE HEAVY ATOMS or IONS jw TO EACH OXYGEN i and CALCULATE the TETRAHEDRALITY

RADMAX  = 25.0D0
RADMIN  =  0.0D0
JT      = 1
DO WHILE(JT.LE.4) 
!   do jw = nmolsol*natmsol+1,natms,nwatsites                                      !Loop over water oxygens
   do jw = 1,natms                                                                 !Loop over solute heavy atoms, water oxygens, and ions
      atom_symb=TRIM(atomname(jw))
!check      write(*,*)atom_symb
!      IF(jw/=i.and.atom_symb(jw)/='H'.and.atom_symb(jw)/='M')THEN                  !Exclude water and solute H atoms and M sites of TIP4P-ew water 
      IF(jw/=i.and.atom_symb/='H'.and.atom_symb/='M')THEN                           !Exclude water and solute H atoms and M sites of TIP4P-ew water 
!check         write(*,*)atom_symb
! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND ATOMS jw     
         dx = x(i) - x(jw)             
         dy = y(i) - y(jw)
         dz = z(i) - z(jw)
         dx = dx - dnint(dx/cell(1))*cell(1)
         dy = dy - dnint(dy/cell(2))*cell(2)
         dz = dz - dnint(dz/cell(3))*cell(3)
         dr2 = dx**2 + dy**2 + dz**2
         dr  = dsqrt(dr2)
 
! FIND THE 4 NEAREST OXYGEN ATOMS jw TO THE OXYGEN ATOM i
     
         IF((dr.GT.RADMIN).AND.(dr.LT.RADMAX))THEN
            RADMAX = dr
! NTDWO STORES THE O ATOMS NEIGHBORS OF I
! JT = 1 FIRST NEIGHBOR ; JT = 2 SECOND NEIGHBOR ; JT = 3 THIRD NEIGHBOR ; JT = 4 FOURTH NEIGHBOR
            NTDWO(JT) = jw
         ENDIF
      ENDIF
   end do
   JT = JT + 1
   RADMIN = RADMAX
   RADMAX = 25.0D0
END DO


! LOOK FOR ERRORS
DO IE=1,3
   DO JE = IE+1,4
      if(NTDWO(IE).EQ.NTDWO(JE))then
         write(*,*)'Error - water 4 nearest neighbors ',i
         stop
      endif   
   END DO
END DO

! CALCULATE FOR EACH WATER MOLECULE THE TETRAHEDRAL PARAMETER
! CHECK FOR SIX ANGLES FOR EACH WATER MOLECULE BASED UPON DOT PRODUCTS OF THE FOUR UNIT VECTORS
! DEFINING THE VERTICES OF A TETRAHEDRON WHERE EACH OXYGEN I IS AT THE TETRAHEDRON CENTRE
! CONSIDER FOUR OXYGEN ATOMS IN THE VERTICES OF THE TETRAHEDRON

 SG     = 0.0D0                       !tetrahedrality
 CHHASG = 0.0D0                       !tetrahedrality

DO KT = 1,4
! oxygen or solute heavy atom ID            
   KI=NTDWO(KT)
   dx = x(i)-x(KI)
   dy = y(i)-y(KI)
   dz = z(i)-z(KI)
   dx = dx - dnint(dx/cell(1))*cell(1)
   dy = dy - dnint(dy/cell(2))*cell(2)
   dz = dz - dnint(dz/cell(3))*cell(3)
   dr2 = dx**2 + dy**2 + dz**2
   dr  = dsqrt(dr2)
! CALCULATE THE FOUR UNIT VECTORS
   XT(KT)=dx/dr
   YT(KT)=dy/dr
   ZT(KT)=dz/dr
END DO

DO JT = 1,3
   DO KT = JT+1,4
! dot product
      THDP=XT(JT)*XT(KT)+YT(JT)*YT(KT)+ZT(JT)*ZT(KT)
      SG = SG + (THDP+1.0D0/3.0D0)**2.0
   END DO
END DO
! Chau and Hardwick normalization
 CHHASG=SG*3.0D0/32.0D0
! Errington and Debenedetti normalization
 SG=1.0D0-SG*3.0D0/8.0D0

return
END SUBROUTINE sol_tetrahedron


!NG SUBROUTINE tcf_solv_tau(nprt,kr,dt,tcf_mean,ko,norg_0HB,TSTEP,NDELS,NRSHMAX)
SUBROUTINE tcf_solv_tau(nprt,kr,dt,tcf_mean,koc,norg_0HB,TSTEP,NDELS,NRSHMAX)
!Calculate the integral of any tcf
!integer,intent(in)                                    :: nprt,kr,ko
integer,intent(in)                                    :: nprt,kr
integer,dimension(NRSHMAX),intent(in)                 :: koc
real(kind=8),intent(in)                               :: dt
real(kind=8),dimension(NDELS,NRSHMAX),intent(in)      :: tcf_mean 
integer,dimension(NRSHMAX),intent(in)                 :: norg_0HB
real(kind=8),intent(in)                               :: TSTEP
integer,intent(in)                                    :: NDELS,NRSHMAX

!Local
integer                   :: J, K
real(kind=8)              :: PSECS,STEP
integer                   :: NNDELS,KINTVL
real(kind=8)              :: AUNDER1,AUNDER2,AUNDER3,AUNDER4
real(kind=8)              :: AABOVE1,AABOVE2,AABOVE3,AABOVE4
real(kind=8)              :: AABOVE,AUNDER,AABOVEERR,AUNDERERR
real(kind=8)              :: ROT_tau,ERROR_tau

KINTVL = 1                         !no longer used
STEP = dt*1.0d-3                   !time-step in ps
!
!DEFINITION OF THE NUMBER OF DELAY TIMES FOR THE FIRST LOOP
!
NNDELS = 10  
70 IF(NNDELS>NDELS)NNDELS=NDELS

AUNDER1 = 0.0D0
! AABOVE1 = (KINTVL*STEP)*XSQAV            !if j = 1 <=> t = dt ps               XSQAV = 1; KINTVL = 1
AABOVE1 = 0.0D0
AUNDER2 = 0.0D0
AABOVE2 = 0.0D0
AUNDER3 = 0.0D0
AABOVE3 = 0.0D0
AUNDER4 = 0.0D0
AABOVE4 = 0.0D0

DO 5 J=1,NNDELS-1                          !j = 1 <=> t = 0 ps

!NG   IF((tcf_mean(J,kr).GT.(0.0)).AND.(tcf_mean(J,kr).GT.tcf_mean(J+1,kr)))THEN
!NG      AABOVE1=AABOVE1+(KINTVL*STEP)*tcf_mean(J,kr)/dfloat(ko-norg_0HB(kr))
!NG   ELSEIF((tcf_mean(J,kr).LT.(0.0)).AND.(tcf_mean(J,kr).GT.tcf_mean(J+1,kr)))THEN
!NG      AABOVE2=AABOVE2+(KINTVL*STEP)*DABS(tcf_mean(J,kr))/dfloat(ko-norg_0HB(kr))
!NG   ELSEIF((tcf_mean(J,kr).LT.(0.0)).AND.(tcf_mean(J,kr).LT.tcf_mean(J+1,kr)))THEN
!NG      AABOVE3=AABOVE3+(KINTVL*STEP)*DABS(tcf_mean(J,kr))/dfloat(ko-norg_0HB(kr)) 
!NG   ELSEIF((tcf_mean(J,kr).GT.(0.0)).AND.(tcf_mean(J,kr).LT.tcf_mean(J+1,kr)))THEN
!NG      AABOVE4=AABOVE4+(KINTVL*STEP)*tcf_mean(J,kr)/dfloat(ko-norg_0HB(kr))
!NG   ENDIF
   IF((tcf_mean(J,kr).GT.(0.0)).AND.(tcf_mean(J,kr).GT.tcf_mean(J+1,kr)))THEN
      AABOVE1=AABOVE1+(KINTVL*STEP)*tcf_mean(J,kr)/dfloat(koc(kr)-norg_0HB(kr))
   ELSEIF((tcf_mean(J,kr).LT.(0.0)).AND.(tcf_mean(J,kr).GT.tcf_mean(J+1,kr)))THEN
      AABOVE2=AABOVE2+(KINTVL*STEP)*DABS(tcf_mean(J,kr))/dfloat(koc(kr)-norg_0HB(kr))
   ELSEIF((tcf_mean(J,kr).LT.(0.0)).AND.(tcf_mean(J,kr).LT.tcf_mean(J+1,kr)))THEN
      AABOVE3=AABOVE3+(KINTVL*STEP)*DABS(tcf_mean(J,kr))/dfloat(koc(kr)-norg_0HB(kr)) 
   ELSEIF((tcf_mean(J,kr).GT.(0.0)).AND.(tcf_mean(J,kr).LT.tcf_mean(J+1,kr)))THEN
      AABOVE4=AABOVE4+(KINTVL*STEP)*tcf_mean(J,kr)/dfloat(koc(kr)-norg_0HB(kr))
   ENDIF
   

5 CONTINUE

DO 15 K=2,NNDELS                           !j = 2 <=> t = dt ps

   if(K < NNDELS)then

!NG      IF((tcf_mean(K,kr).GT.(0.0)).AND.(tcf_mean(K,kr).GT.tcf_mean(K+1,kr)))THEN
!NG         AUNDER1=AUNDER1+(KINTVL*STEP)*tcf_mean(K,kr)/dfloat(ko-norg_0HB(kr))
!NG      ELSEIF((tcf_mean(K,kr).LT.(0.0)).AND.(tcf_mean(K,kr).GT.tcf_mean(K+1,kr)))THEN
!NG         AUNDER2=AUNDER2+(KINTVL*STEP)*DABS(tcf_mean(K,kr))/dfloat(ko-norg_0HB(kr))
!NG      ELSEIF((tcf_mean(K,kr).LT.(0.0)).AND.(tcf_mean(K,kr).LT.tcf_mean(K+1,kr)))THEN
!NG         AUNDER3=AUNDER3+(KINTVL*STEP)*DABS(tcf_mean(K,kr))/dfloat(ko-norg_0HB(kr))
!NG      ELSEIF((tcf_mean(K,kr).GT.(0.0)).AND.(tcf_mean(K,kr).LT.tcf_mean(K+1,kr)))THEN
!NG         AUNDER4=AUNDER4+(KINTVL*STEP)*tcf_mean(K,kr)/dfloat(ko-norg_0HB(kr))
!NG      ENDIF
      IF((tcf_mean(K,kr).GT.(0.0)).AND.(tcf_mean(K,kr).GT.tcf_mean(K+1,kr)))THEN
         AUNDER1=AUNDER1+(KINTVL*STEP)*tcf_mean(K,kr)/dfloat(koc(kr)-norg_0HB(kr))
      ELSEIF((tcf_mean(K,kr).LT.(0.0)).AND.(tcf_mean(K,kr).GT.tcf_mean(K+1,kr)))THEN
         AUNDER2=AUNDER2+(KINTVL*STEP)*DABS(tcf_mean(K,kr))/dfloat(koc(kr)-norg_0HB(kr))
      ELSEIF((tcf_mean(K,kr).LT.(0.0)).AND.(tcf_mean(K,kr).LT.tcf_mean(K+1,kr)))THEN
         AUNDER3=AUNDER3+(KINTVL*STEP)*DABS(tcf_mean(K,kr))/dfloat(koc(kr)-norg_0HB(kr))
      ELSEIF((tcf_mean(K,kr).GT.(0.0)).AND.(tcf_mean(K,kr).LT.tcf_mean(K+1,kr)))THEN
         AUNDER4=AUNDER4+(KINTVL*STEP)*tcf_mean(K,kr)/dfloat(koc(kr)-norg_0HB(kr))
      ENDIF      
      
   else
!      AUNDER1=AUNDER1+(KINTVL*STEP)*tcf_mean(K,kr)/dfloat(ko-norg_0HB(kr))
      AUNDER1=AUNDER1+(KINTVL*STEP)*tcf_mean(K,kr)/dfloat(koc(kr)-norg_0HB(kr)) 
   endif
         
15 CONTINUE  

AABOVE=AABOVE1+AABOVE4-AABOVE2-AABOVE3
AUNDER=AUNDER1+AUNDER4-AUNDER2-AUNDER3

AABOVEERR=AABOVE1+AABOVE4+AABOVE2+AABOVE3
AUNDERERR=AUNDER1+AUNDER4+AUNDER2+AUNDER3

!ARITHMETIC AVERAGE OF THE TWO AREAS 

ROT_tau = (AUNDER+AABOVE)/2.0D0   
ERROR_tau = (AABOVEERR-AUNDERERR)/2.0D0  

!   WRITE(n15,89)NNDELS,ROT_tau,ERROR_tau

PSECS=dfloat(NNDELS)*TSTEP
WRITE(nprt,99)PSECS,ROT_tau

IF(NNDELS.LT.NDELS)THEN
   NNDELS = NNDELS+10
   GOTO 70
ELSEIF (NNDELS.GE.NDELS)THEN
   GOTO 75
75 ENDIF  
   
!   89 FORMAT(5X,I9,9X,E11.5,5X,E9.3) 
99 FORMAT(9X,E14.5,7X,E12.5)

return

END SUBROUTINE tcf_solv_tau



SUBROUTINE msd_LSF_Einstein(nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end,a,b,d,r,n)
!Calculate the Einstein diffusion through Linear Squares fit of the MSD
integer,intent(in)                                    :: nroute,nsyst,nf,kr,NDELS,n_LSF_begin,n_LSF_end
real(kind=8),intent(out)                              :: a,b,d,r        !Y = a + b*X; d = standard deviation of the fit ; r = correlation coefficient
integer,intent(out)                                   :: n              !number of points used in the fit: range [n_LSF_begin,n_LSF_end[

!Local
integer                                               :: i,j,k
real(kind=8)                                          :: Xt,Ymsd
!real(kind=8),dimension(NDELS)                         :: X,Y
real(kind=8),dimension(:),allocatable                 :: X,Y

allocate(X(NDELS),Y(NDELS))

if(nsyst==0.and.nroute==2)open(nf,file='MSD_PW_out/msd_mean_bulk.dat',status='old',action='read')                           !Pure water
if(nsyst==0.and.nroute==3)then                                                                                              !Pure water LDL/HDL environments
   if(nf==150)open(nf,file='MSD_LDL_HDL/Bulk.dat',status='old',action='read')                 
   if(nf==160)open(nf,file='MSD_LDL_HDL/LDL.dat',status='old',action='read')
   if(nf==170)open(nf,file='MSD_LDL_HDL/HDL.dat',status='old',action='read')
   if(nf==180)open(nf,file='MSD_LDL_HDL/other.dat',status='old',action='read')
endif
if(nsyst==0.and.nroute==4)then                                                                                              !Pure water Neighbors environments
   if(nf==150)open(nf,file='MSD_PW_out/msd_bulk.dat',status='old',action='read')
   if(nf==160)open(nf,file='MSD_PW_out/msd_4LNbE.dat',status='old',action='read')
   if(nf==170)open(nf,file='MSD_PW_out/msd_3LNbE.dat',status='old',action='read')
   if(nf==180)open(nf,file='MSD_PW_out/msd_5MNbE.dat',status='old',action='read')
endif   

!Solvation
if(nroute==0)then
   if(nsyst==1.and.nf==110)open(nf,file='MSD_env_out/msd_bulk.dat',status='old',action='read')
   if(nsyst==1.and.nf==120)open(nf,file='MSD_env_out/msd_hsh.dat',status='old',action='read')
   if(nsyst==1.and.nf==130)open(nf,file='MSD_env_out/msd_th.dat',status='old',action='read')
   if(nsyst==1.and.nf==140)open(nf,file='MSD_env_out/msd_nth.dat',status='old',action='read')
endif
if(nroute==1)then
   if(nsyst==1.and.nf==110)open(nf,file='MSD_env_out_stat/msd_bulk.dat',status='old',action='read')
   if(nsyst==1.and.nf==120)open(nf,file='MSD_env_out_stat/msd_hsh.dat',status='old',action='read')
   if(nsyst==1.and.nf==130)open(nf,file='MSD_env_out_stat/msd_th.dat',status='old',action='read')
   if(nsyst==1.and.nf==140)open(nf,file='MSD_env_out_stat/msd_nth.dat',status='old',action='read')
endif

n = 0
a = 0.0d0
b = 0.0d0
d = 0.0d0
r = 0.0d0

!check write(*,*)NDELS,n_LSF_begin,n_LSF_end
read(nf,*)
read(nf,*)
do j=1,NDELS
   read(nf,*)Xt,Ymsd                                !read time(ps) and msd(Ang**2/ps) 
   if(j>=n_LSF_begin.and.j<=n_LSF_end)then
      n = n + 1
      X(n) = Xt
      Y(n) = Ymsd
   endif 
end do   

if(n>=10)then
   write(*,*)
   write(*,*)'Number of points for Linear Squares Fit =',n
   write(*,*)
   CALL Least_Square(n,X,Y,NDELS,a,b,d,r)
else
   write(*,*)'Number of msd points for D calculation too low ',n 
endif

close(nf)

!check   write(*,*)a,b,d

deallocate(X,Y)


return

END SUBROUTINE msd_LSF_Einstein



SUBROUTINE Least_Square(n,X,Y,NDELS,a,b,d,r)
!Linear least squares subroutine  
!The input data set is X(m), Y(m)
!The number of data points is n (n must be > 2).
!The returned parameters are: a,b, coefficients of equation Y = a + b X
!r and d, correlation coeficient and standard deviation of the fit.
integer,intent(in)                                                 :: n,NDELS
real(kind=8),dimension(NDELS),intent(in)                           :: X,Y
REAL*8,intent(out)                                                 :: a,b,d,r
  
!Local  
REAL*8                                          :: a1,a2,b0,b1,d1,b2,b_
INTEGER*4                                       :: m
a1 = 0.0
a2 = 0.0
b0 = 0.0
b1 = 0.0
b2 = 0.0
b_ = 0.0
!NG  DO m = 0, n-1
DO m = 1, n
   a1 = a1 + X(m)
   a2 = a2 + X(m) * X(m)
   b0 = b0 + Y(m)
   b1 = b1 + Y(m) * X(m)
!For the calculation of r   
   b2 = b2 + Y(m) * Y(m)
END DO
a1 = a1 / n
a2 = a2 / n
b0 = b0 / n
b1 = b1 / n
!For the calculation of r
b2 = b2 / n

d = a1 * a1 - a2
a = a1 * b1 - a2 * b0
a = a / d
b = a1 * b0 - b1
!For the calculation of r
b_ = b
b = b / d
!calculate r (ref. Wolfram - correlation coefficient)
d = b0 * b0 - b2                  !denominator of b_
b_= b_ / d
r = dsqrt(b * b_)
!Evaluation of standard deviation d (unbiased estimate) 
d = 0
!NG  DO m = 0, n-1
DO m = 1, n
   d1 = Y(m) - a - b * X(m)                      !difference between the data and the fit
    d  = d + d1 * d1
END DO
d = DSQRT(d / (n - 2))

!check   write(*,*)a,b,d

RETURN

END SUBROUTINE Least_Square



!NG     SUBROUTINE VORON3(i,natms,nmolsol,natmsol,nmolwat,nwatsites,nions,x,y,z,MAXVER,MAXCAN,MAXN,cell,rcut)
!NG     !CONSTRUCTION OF VORONOI POLYHEDRON IN 3D.
!NG     integer,intent(in)                                 :: i
!NG     integer,intent(in)                                 :: natms,nmolsol,natmsol,nmolwat,nwatsites,nions
!NG     integer,intent(in)                                 :: MAXVER,MAXCAN,MAXN
!NG     real(kind=4),intent(in)                            :: x(natms),y(natms),z(natms)
!NG     real(kind=4),intent(in)                            :: rcut
!NG     real(kind=8),intent(in)                            :: cell(3) 
!NG     
!NG     !Local Variables
!NG     integer                                            :: jw
!NG     real(kind=8)                                       :: dx,dy,dz,dr,dr2   
!NG     real(kind=4)                                       :: BOX,RCUTSQ, BOXINV
!NG     integer                                            :: NABLST(MAXVER,MAXN), NNAB(MAXN), INAB, JNAB
!NG     real(kind=4)                                       :: PX(MAXCAN), PY(MAXCAN), PZ(MAXCAN), PS(MAXCAN)
!NG     integer                                            :: TAG(MAXCAN), EDGES(MAXCAN)
!NG     
!NG     real(kind=4)                                       :: RXVER(MAXVER), RYVER(MAXVER), RZVER(MAXVER)
!NG     integer                                            :: IVER(MAXVER), JVER(MAXVER), KVER(MAXVER)
!NG     
!NG     integer                                            :: CAN                                                !candidates
!NG     INTEGER                                            :: NCAN, NVER, NEDGE, NFACE, NCOORD
!NG     
!NG     
!NG     BOX = cell(1)
!NG     RCUTSQ = RCUT**2.0
!NG     BOXINV = 1.0 / BOX
!NG     NABLST = 0
!NG     NNAB   = 0
!NG     CAN    = 0
!NG     
!NG     ! ** SELECT CANDIDATES **
!NG     
!NG     do jw = nmolsol*natmsol+1,natms-nions,nwatsites        !Loop over water oxygens
!NG        IF(jw.NE.i)THEN  
!NG     ! COMPONENTS OF VECTOR DISTANCE BETWEEN OXYGEN ATOM i AND OXYGEN ATOMS jw     
!NG           dx = x(i) - x(jw)             
!NG           dy = y(i) - y(jw)
!NG           dz = z(i) - z(jw)
!NG           dx = dx - dnint(dx/cell(1))*cell(1)
!NG           dy = dy - dnint(dy/cell(2))*cell(2)
!NG           dz = dz - dnint(dz/cell(3))*cell(3)
!NG           dr2 = dx**2 + dy**2 + dz**2
!NG     !      dr  = dsqrt(dr2)
!NG           IF ( dr2 .LT. RCUTSQ ) THEN      
!NG              CAN = CAN + 1
!NG              IF ( CAN .GT. MAXCAN ) THEN
!NG                 WRITE(*,'('' TOO MANY CANDIDATES '')')
!NG                 STOP
!NG              ENDIF
!NG     !check ok         write(*,*) i,jw,sqrt(dr2)
!NG     !WARNING!!! CAN takes values 1, 2, 3, etc not the oxygen index (i)
!NG              PX(CAN)  = dx
!NG              PY(CAN)  = dy
!NG              PZ(CAN)  = dz
!NG              PS(CAN)  = dr2
!NG              TAG(CAN) = jw
!NG             
!NG           ENDIF
!NG        ENDIF   
!NG     end do
!NG     !** CANDIDATES HAVE BEEN SELECTED **
!NG     
!NG             NCAN = CAN
!NG     
!NG     !** SORT INTO ASCENDING ORDER OF DISTANCE **
!NG     !** THIS SHOULD IMPROVE EFFICIENCY        **
!NG     
!NG             CALL SORT ( MAXCAN, PX, PY, PZ, PS, TAG, NCAN )
!NG     
!NG     !** PERFORM VORONOI ANALYSIS **
!NG     
!NG     
!NG             CALL WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, NFACE, PX, PY, PZ, PS, EDGES, RXVER, RYVER, RZVER, IVER, JVER, KVER )
!NG     
!NG     
!NG     
!NG     RETURN
!NG     
!NG     END SUBROUTINE VORON3


!     CALL SORT ( MAXCAN, PX, PY, PZ, PS, TAG, NCAN )
SUBROUTINE SORT ( MAXCAN, RX, RY, RZ, RS, TAG, NN, xn, yn, zn, RD )

!    *******************************************************************
!    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
!    **                                                               **
!    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL     **
!    *******************************************************************

INTEGER                               :: MAXCAN, NN
REAL(KIND=4)                          :: RX(MAXCAN), RY(MAXCAN), RZ(MAXCAN), RS(MAXCAN), RD(MAXCAN)
INTEGER                               :: TAG(MAXCAN)
LOGICAL                               :: CHANGE
INTEGER                               :: I, ITOP, I1, TAGI
REAL(KIND=4)                          :: RXI, RYI, RZI, RSI, RDI
REAL(KIND=4)                          :: xn(MAXCAN), yn(MAXCAN), zn(MAXCAN) 
REAL(KIND=4)                          :: xxn, yyn, zzn
!NN - NUMBER OF CANDIDATES = NCAN 
!RS - SQUARED DISTANCE     = PS = dr2

CHANGE = .TRUE.
ITOP = NN - 1

1000    IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

           CHANGE = .FALSE.

           DO 100 I = 1, ITOP    !NG   LOOP OVER [Number of Candidates - 1]

              I1 = I + 1

              IF ( RS(I) .GT. RS(I1) ) THEN       !NG I neighbor farther away than I+1 neighbor

                 RXI  = RX(I)
                 RYI  = RY(I)
                 RZI  = RZ(I)
                 RSI  = RS(I)
                 TAGI = TAG(I)
                 
                 xxn  = xn(I)
                 yyn  = yn(I)
                 zzn  = zn(I)                 
                 RDI  = RD(I)
                 
!exchange particles
                 RX(I)  = RX(I1)
                 RY(I)  = RY(I1)
                 RZ(I)  = RZ(I1)
                 RS(I)  = RS(I1)
                 TAG(I) = TAG(I1)
                 
                 xn(I) = xn(I1)
                 yn(I) = yn(I1)
                 zn(I) = zn(I1)
                 RD(I) = RD(I1)

                 RX(I1)  = RXI
                 RY(I1)  = RYI
                 RZ(I1)  = RZI
                 RS(I1)  = RSI
                 TAG(I1) = TAGI
                 
                 xn(I1) = xxn
                 yn(I1) = yyn
                 zn(I1) = zzn
                 RD(I1) = RDI 
                 
                 CHANGE = .TRUE.

              ENDIF

100        CONTINUE

           ITOP = ITOP - 1
           GOTO 1000

        ENDIF

RETURN
END SUBROUTINE SORT
 
!               ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, NFACE, PX, PY, PZ, PS, EDGES, RXVER, RYVER, RZVER, IVER, JVER, KVER, Area_Vor, Vol_Vor )
!     CALL WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, NFACE, PX, PY, PZ, PS, EDGES, RXVER, RYVER, RZVER, IVER, JVER, KVER )    
SUBROUTINE WORK ( ICP,MAXCAN,MAXV,NN,NV,NE,NF,RX,RY,RZ,RS,EDGES,VX,VY,VZ,IV,JV,KV,Area_t,Vol_t,xi,yi,zi,xn,yn,zn,RD,LEULER )

!    *******************************************************************
!    ** ROUTINE TO PERFORM VORONOI ANALYSIS                           **
!    **                                                               **
!    ** WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                **
!    ** I.E. THE FACES OF THE POLYHEDRON GO THROUGH THE POINTS.       **
!    *******************************************************************

INTEGER              :: ICP, MAXCAN, NN, MAXV, NV, NE, NF
INTEGER              :: EDGES(MAXCAN)
REAL(KIND=4)         :: RX(MAXCAN), RY(MAXCAN), RZ(MAXCAN), RS(MAXCAN), RD(MAXCAN)
REAL(KIND=4)         :: xn(MAXCAN), yn(MAXCAN), zn(MAXCAN)
REAL(KIND=4)         :: VX(MAXV), VY(MAXV), VZ(MAXV)
INTEGER              :: IV(MAXV), JV(MAXV), KV(MAXV)
logical              :: LEULER

real(kind=4)         :: Area_t, Vol_t

LOGICAL              ::  OK, OK_
INTEGER              :: I, J, K, L, NN1, NN2, N, V
REAL(KIND=4)         :: AI, BI, CI, DI, AJ, BJ, CJ, DJ, AK, BK, CK, DK
REAL(KIND=4)         :: AB, BC, CA, DA, DB, DC, DET, DETINV
REAL(KIND=4)         :: VXIJK, VYIJK, VZIJK
REAL(KIND=4)         :: TOL
!PARAMETER ( TOL = 1.E-6 ) 
PARAMETER ( TOL = 1.E-6 )

!NG Variables
integer              :: nedge_face(MAXCAN), nvert_face(MAXCAN)                !Number of edges of each face; Number of vertices of each face
integer              :: nzv, nzf, nm, m, NFV, nmm, nloop, mm, na
real(kind=4)         :: vert_face_x(MAXCAN,MAXV), vert_face_y(MAXCAN,MAXV), vert_face_z(MAXCAN,MAXV)
integer              :: ivv(MAXCAN,MAXV), jvv(MAXCAN,MAXV), kvv(MAXCAN,MAXV)
real(kind=4)         :: vert_face_x_sort, vert_face_y_sort, vert_face_z_sort
integer              :: jvv_sort, kvv_sort
logical              :: change
real(kind=4)         :: Area_f(MAXCAN), Vol_f(MAXCAN)
real(kind=4)         :: Ax_comp, Ay_comp, Az_comp, A_cross,A_cross_dot
real(kind=4)         :: Ax_c, Ay_c, Az_c
real(kind=4)         :: A_dot
real(kind=4)         :: vert_dif_a_x, vert_dif_a_y, vert_dif_a_z, vert_dif_b_x, vert_dif_b_y, vert_dif_b_z 
real(kind=4)         :: vert_dif_c_x, vert_dif_c_y, vert_dif_c_z, norm_c
real(kind=4)         :: norm_a, norm_b, theta, ARG
real(kind=4)         :: R_vert_L, R_vert_IC, R_vert_N1C, R_vert_N2C, R_vert_N3C
real(kind=4)         :: xi, yi, zi
real(kind=4)         :: A_Heron_S, A_Heron, A_Heron_f(MAXCAN), Vcros_dot
real(kind=4)         :: pln_L_V, pln_L_C
real(kind=4)         :: ddi, ddj, ddk
integer              :: MVORON             ! 1: Method 1; 2:Method 2
real(kind=4)         :: vxijk_, vyijk_, vzijk_

!End NG Variables

!NN - NUMBER OF CANDIDATES = NCAN 
!RS - SQUARED DISTANCE     = PS = dr2

vert_face_x = 0; vert_face_y = 0; vert_face_z = 0
Area_f    = 0.0
Area_t    = 0.0
Vol_f     = 0.0
Vol_t     = 0.0
A_Heron_f = 0.0
Vcros_dot = 0.0
MVORON    = 1 
                
!    ** IF THERE ARE LESS THAN 4 POINTS GIVEN **
!    ** WE CANNOT CONSTRUCT A POLYHEDRON      **

IF ( NN .LT. 4 ) THEN

   WRITE(*,'('' LESS THAN 4 POINTS GIVEN TO WORK '',I5)') NN
   STOP

ENDIF

NN1 = NN - 1
NN2 = NN - 2
V = 0

!    ** WE AIM TO EXAMINE EACH POSSIBLE VERTEX  **
!    ** DEFINED BY THE INTERSECTION OF 3 PLANES **
!    ** EACH PLANE IS SPECIFIED BY RX,RY,RZ,RS  **

        DO 400 I = 1, NN2                !NG   I ARE THE NEIGHBOR CANDIDATES - ordered by the lowest distance to the central particle

           AI =  RX(I)                   !NG   RX(I) ARE THE CARTESIAN COORDINATES OF THE VECTOR DISTANCE (I - CENTRAL PARTICLE)
           BI =  RY(I)
           CI =  RZ(I)
           DI = -RS(I)
           
           ddi = RD(I)                   !NG   SECOND METHOD - bisector plane

           DO 300 J = I + 1, NN1

              AJ =  RX(J)
              BJ =  RY(J)
              CJ =  RZ(J)
              DJ = -RS(J)
              
              ddj = RD(J)                !NG   SECOND METHOD
              
!NG   cross product
              AB = AI * BJ - AJ * BI
              BC = BI * CJ - BJ * CI
              CA = CI * AJ - CJ * AI
              DA = DI * AJ - DJ * AI
              DB = DI * BJ - DJ * BI
              DC = DI * CJ - DJ * CI

              DO 200 K = J + 1, NN

                 AK =  RX(K)
                 BK =  RY(K)
                 CK =  RZ(K)
                 DK = -RS(K)
                 
                 ddk = RD(K)             !NG   SECOND METHOD
                 
!NG   THIS IS OK - IT IS THE 3*3 DETERMINANT BUILT WITH Ri (LINE 1), Rj (LINE 2) and Rk (LINE 3)

                 DET = AK * BC + BK * CA + CK * AB

                 IF ( ABS ( DET ) .GT. TOL ) THEN
                    
!                ** THE PLANES INTERSECT **

                    DETINV = 1.0 / DET

                    if(MVORON == 1)then
                       VXIJK = ( - DK * BC + BK * DC - CK * DB ) * DETINV
                       VYIJK = ( - AK * DC - DK * CA + CK * DA ) * DETINV
                       VZIJK = (   AK * DB - BK * DA - DK * AB ) * DETINV
                    elseif(MVORON == 2)then                    
!NG   SECOND METHOD
!Method 1                       ddi=DI
!Method 1                       ddj=DJ
!Method 1                       ddk=DK
                       
                       vxijk = -ddi*BJ*CK - ddj*BK*CI - ddk*BI*CJ + ddk*BJ*CI + ddj*BI*CK + ddi*BK*CJ   
                       vyijk = -ddj*AI*CK - ddi*AK*CJ - ddk*AJ*CI + ddj*AK*CI + ddi*AJ*CK + ddk*AI*CJ                    
                       vzijk = -ddk*AI*BJ - ddi*AJ*BK - ddj*AK*BI + ddi*AK*BJ + ddj*AI*BK + ddk*AJ*BI
                       
                       vxijk =  vxijk*DETINV
                       vyijk =  vyijk*DETINV
                       vzijk =  vzijk*DETINV
                    endif
!NG check code                    
!Method 1                       ddi=DI
!Method 1                       ddj=DJ
!Method 1                       ddk=DK
!                       vxijk_ = -ddi*BJ*CK - ddj*BK*CI - ddk*BI*CJ + ddk*BJ*CI + ddj*BI*CK + ddi*BK*CJ                       
!                       vyijk_ = -ddj*AI*CK - ddi*AK*CJ - ddk*AJ*CI + ddj*AK*CI + ddi*AJ*CK + ddk*AI*CJ                    
!                       vzijk_ = -ddk*AI*BJ - ddi*AJ*BK - ddj*AK*BI + ddi*AK*BJ + ddj*AI*BK + ddk*AJ*BI
!                            
!                       vxijk_ = vxijk_*DETINV
!                       vyijk_ = vyijk_*DETINV
!                       vzijk_ = vzijk_*DETINV
!                        write(*,'(1x,6(1x,F9.3))')VXIJK,vxijk_,VYIJK,vyijk_,VZIJK,vzijk_
! end check code
                       
!NG   END_SECOND METHOD
                    
!NG   DISTANCE OF VERTICE TO THE CENTRAL PARTICLE                    
!                    R_vert_IC = sqrt((VXIJK - xi)**2.0+(VYIJK - yi)**2.0+(VZIJK - zi)**2.0)
!                    R_vert_N1C = sqrt((VXIJK - xn(i))**2.0+(VYIJK - yn(i))**2.0+(VZIJK - zn(i))**2.0)
!                    R_vert_N2C = sqrt((VXIJK - xn(j))**2.0+(VYIJK - yn(j))**2.0+(VZIJK - zn(j))**2.0)
!                    R_vert_N3C = sqrt((VXIJK - xn(k))**2.0+(VYIJK - yn(k))**2.0+(VZIJK - zn(k))**2.0)
!                    write(*,'(3(1x,I4),4(1x,F9.2))')i,j,k,R_vert_IC,R_vert_N1C,R_vert_N2C,R_vert_N3C

!                ** NOW WE TAKE SHOTS AT THE VERTEX **
!                ** USING THE REMAINING PLANES .... **

                    OK = .TRUE.
                    L  = 1

100                 IF ( OK .AND. ( L .LE. NN ) ) THEN
                       IF ( ( L .NE. I ) .AND.( L .NE. J ) .AND. ( L .NE. K ) ) THEN
                          if(MVORON == 1)then
                          
                             OK = ( ( RX(L) * VXIJK + RY(L) * VYIJK + RZ(L) * VZIJK ) .LE. RS(L) )
                             
                          elseif(MVORON == 2)then
!Method 1                             RD(L)=-RS(L)
                             pln_L_V = RX(L)*VXIJK + RY(L)*VYIJK + RZ(L)*VZIJK + RD(L)
                             pln_L_C = RX(L)*xi    + RY(L)*yi    + RZ(L)*zi    + RD(L)

                             OK_ = (pln_L_V*pln_L_C >= 0)
                            
                             OK = ( ( RX(L) * VXIJK + RY(L) * VYIJK + RZ(L) * VZIJK ) .LE. -RD(L) )
                             if(OK_ /= OK)write(*,*)'Failed VERTICE second condition check'
                             
                          endif   
!NG Lab
!NG   DISTANCE OF NEIGHBOR L TO THE VERTICE
!check                          R_vert_L = sqrt((xn(L)-VXIJK)**2.0+(yn(L)-VYIJK)**2.0+(zn(L)-VZIJK)**2.0)
!NG Lab                           
                       ENDIF
                       L = L + 1
                       GOTO 100
                    ENDIF

!                ** IF THE VERTEX MADE IT      **
!                ** ADD IT TO THE HALL OF FAME **
!                ** CONVERT TO CORRECT SCALE   **

                    IF ( OK ) THEN

                       V = V + 1
                      
                       IF ( V .GT. MAXV ) STOP 'TOO MANY VERTICES'

                       IV(V)  = I                                       !NG   NEIGHBOR id - FOR EACH VERTEX WE GET 3 PARTICLES ASSOCIATED WITH THAT VERTEX
                       JV(V)  = J                                       !NG   VERTICE 1, 2, 3, 4,... ASSOCIATED WITH ANY THREE OXYGEN ATOMS
                       KV(V)  = K
                       if(MVORON == 1)then
                          VX(V) = 0.5 * VXIJK                           !NG   THIS IS BECAUSE PLANES WERE NOT DEFINED ON THE BISECTOR
                          VY(V) = 0.5 * VYIJK
                          VZ(V) = 0.5 * VZIJK
                       elseif (MVORON == 2)then
                          VX(V) = VXIJK
                          VY(V) = VYIJK
                          VZ(V) = VZIJK
                       endif 
                       
!CHECK-1                       write(*,'(3(2x,I5),3(2x,F9.3))') IV(V),JV(V),KV(V),VX(V),VY(V),VZ(V)
                    ENDIF

                 ENDIF

200           CONTINUE

300        CONTINUE

400     CONTINUE
!CHECK-1        write(*,*)

        NV = V
        
!        write(*,*)MVORON,NV
!        stop

        IF ( NV .LT. 4 ) THEN
           WRITE(*,'('' LESS THAN 4 VERTICES FOUND IN WORK '',I5)') NV
           STOP
        ENDIF

!    ** IDENTIFY NEIGHBOURING POINTS **

        DO 500 N = 1, NN
           EDGES(N) = 0
!NG   
           nedge_face(N) = 0                                 !NG   NUMBER OF EDGES OF EACH FACE - TO BE USED TO DETERMINE THE AREA OF EACH FACE
!End NG           

500     CONTINUE

        DO 600 V = 1, NV

           EDGES(IV(V)) = EDGES(IV(V)) + 1         !NG   SUM NUMBER OF TIMES AN INDEX SHOWS UP EITHER AS I,J OR K - EACH IJK DEFINES A VERTEX
           EDGES(JV(V)) = EDGES(JV(V)) + 1         !NG   THE NUMB OF TIMES AN INDEX APPEARS EITHER AS I, J OR K IS EQUAL TO THE NUMBER OF EDGES/VERTICES OF A FACE
           EDGES(KV(V)) = EDGES(KV(V)) + 1         !NG   THE SUM OF ALL EDGES WILL GIVE 2*NUMBER OF EDGES = 2*NE - each edge is shared by two faces

600     CONTINUE

!    ** POINTS WITH NONZERO EDGES ARE NEIGHBOURS **

!    ** CHECK EULER RELATION **

!                                                            !NG   Euler Relation Nv - Ne + Nf = 2          !see Ruocco 1991
        NF = 0
        NE = 0

        DO 700 N = 1, NN

           IF ( EDGES(N) .GT. 0 ) NF = NF + 1                 !NG   FOR EACH EDGE ASSOCIATED WITH AN INDEX (DIFFERENT FROM ZERO) SUM A FACE
           
!NG   modified by NG
           IF ( EDGES(N) .GT. 0 ) nedge_face(N)= EDGES(N)
!NG   end modified by NG
           
           NE = NE + EDGES(N)

700     CONTINUE

!NG   INTRODUCED by NG

!IDENTIFY THE VERTICES BELONGING TO EACH FACE 
        nvert_face = nedge_face                    !NG   Numb of vertices of each face = Number of edges of each face 
        nzf = 0
        do N = 1, NN
           nzv = 0                                 !NG   index of the vertice number
           if(nvert_face(N) > 0) then              !NG   The face N exists and has nvert_face(N) vertices
              nzf = nzf + 1                        !NG   not used
!CHECK-2              write(*,*) 'Face', N, nvert_face(N)
!Check which vertices belong to each face              
              do V = 1, NV                         !NG   Loop over the total number of vertices - V is different from nzv (vertice index for each face)  
! V ranges from 1 to the largest vertice number (not necessarily eq to the total number of vertices - some numbers may be skipped)
! nzv ranges between 1 and the number of of vertices for each face
                 if( (IV(V) == N ) .or. ( JV(V) == N ) .or. ( KV(V) == N)) then
                    nzv = nzv + 1                                                   !Number of the vertice
                    vert_face_x(N,nzv) = VX(V)                                      !COMPONENT X OF THE VERTICE NZV BELONGING TO FACE N
                    vert_face_y(N,nzv) = VY(V)
                    vert_face_z(N,nzv) = VZ(V)                    
!                    vert_face_x(nzf,nzv) = VX(V)
!                    vert_face_y(nzf,nzv) = VY(V)
!                    vert_face_z(nzf,nzv) = VZ(V)
                    ivv(N,nzv) = IV(V)
                    jvv(N,nzv) = JV(V)
                    kvv(N,nzv) = KV(V)
!CHECK-3                    write(*,'(a7, 2x,2(1x,i5),3(3x,I5))') 'Vertice', nzv, V, IV(V), JV(V), KV(V) 
!CHECK-3                    write(*,'(i5,3(3x,F9.3))') nzv, vert_face_x(N,nzv), vert_face_y(N,nzv), vert_face_z(N,nzv)
                 endif
              end do   
           endif   
        end do
        
!CHECK-4        write(*,*)
!CHECK-4        write(*,*) 'sort'
!CHECK-4        write(*,*)

! Sort the vertices of each face to calculate the area and volume based on the external product definition   

        do N = 1,NN
           m = 1
           mm = 2
           nloop = 0
1000       if(nvert_face(N) > 0) then
              nloop = nloop + 1
              change = .false.
              NFV = nvert_face(N) 
              do nm = mm, NFV
!CHECK-5                 write(*,*)'INDEX', nm
!CHECK-6                 write(*,'(a7,1x,7(1x,i7))')'***Face', nloop, nm, m, jvv(N,nm), kvv(N,nm), jvv(N,m), kvv(N,m)
                 if ( (jvv(N,nm) == jvv(N,m)) .or. (kvv(N,nm) == kvv(N,m)) .or. (jvv(N,nm) == kvv(N,m)) .or. (kvv(N,nm) == jvv(N,m)) ) then
                    nmm = nm - m
                    if(nmm > 1)then
!CHECK-7                    write(*,*)'**Face', N, nm, jvv(N,nm), kvv(N,nm), jvv(N,m), kvv(N,m)
                       vert_face_x_sort = vert_face_x(N,nm)
                       vert_face_y_sort = vert_face_y(N,nm)
                       vert_face_z_sort = vert_face_z(N,nm)
                       
                       jvv_sort = jvv(N,nm)
                       kvv_sort = kvv(N,nm)
                       
                       vert_face_x(N,nm) = vert_face_x(N,m+1) 
                       vert_face_y(N,nm) = vert_face_y(N,m+1) 
                       vert_face_z(N,nm) = vert_face_z(N,m+1) 
                       
                       jvv(N,nm) = jvv(N,m+1)
                       kvv(N,nm) = kvv(N,m+1)                
!Exchange                    
                       vert_face_x(N,m+1) = vert_face_x_sort
                       vert_face_y(N,m+1) = vert_face_y_sort
                       vert_face_z(N,m+1) = vert_face_z_sort
                       
                       jvv(N,m+1) = jvv_sort 
                       kvv(N,m+1) = kvv_sort
                    
!CHECK-8                    write(*,'(a5,7I5)')'*Face', N, nm, jvv(N,nm), kvv(N,nm), m, jvv(N,m), kvv(N,m) 
                       change = .true.
                    endif
                    m = m + 1
                    mm = mm + 1
                 endif
              end do
              if (change) goto 1000
           endif
        end do 
        
!check sorting     !CHECK_9

!CHECK_9        do N = 1,NN
!CHECK_9           if(nvert_face(N) > 0) then
!CHECK_9              write(*,*) 'Face', N, nvert_face(N)
!CHECK_9              do nm = 1, nvert_face(N)
!CHECK_9                 write(*,'(3(2x,I5),3(2x,F9.3))') ivv(N,nm), jvv(N,nm), kvv(N,nm), vert_face_x(N,nm), vert_face_y(N,nm), vert_face_z(N,nm)
!CHECK_9              end do
!CHECK_9           endif
!CHECK_9        end do
       
!        stop

!NG   CALCULATE AREAS - calculate the area of each face and the total area
        na = 0                                  !NG   Numb of faces counter
        do N = 1,NN                             !NG   Warning - Loop over possible canditates NOT faces - most do not correspond to a face
           if(nvert_face(N) > 0) then
              na = na + 1                                                 !NG   sum number of faces counted for checking
!              write(*,*)
!              write(*,*)'Face',na
              do nm = 2, nvert_face(N) - 1                                !NG   loop over number of vertices
                 vert_dif_a_x = vert_face_x(N,nm) - vert_face_x(N,1)
                 vert_dif_a_y = vert_face_y(N,nm) - vert_face_y(N,1)
                 vert_dif_a_z = vert_face_z(N,nm) - vert_face_z(N,1)
                 norm_a = sqrt( vert_dif_a_x**2.0 + vert_dif_a_y**2.0 + vert_dif_a_z**2.0 )
!                 
                 vert_dif_b_x = vert_face_x(N,nm + 1) - vert_face_x(N,1)
                 vert_dif_b_y = vert_face_y(N,nm + 1) - vert_face_y(N,1)
                 vert_dif_b_z = vert_face_z(N,nm + 1) - vert_face_z(N,1)
                 norm_b = sqrt( vert_dif_b_x**2.0 + vert_dif_b_y**2.0 + vert_dif_b_z**2.0 )
!check purposes only - Heron's formula                 
                 vert_dif_c_x = vert_face_x(N,nm + 1) - vert_face_x(N,nm)
                 vert_dif_c_y = vert_face_y(N,nm + 1) - vert_face_y(N,nm)
                 vert_dif_c_z = vert_face_z(N,nm + 1) - vert_face_z(N,nm)
                 norm_c = sqrt( vert_dif_c_x**2.0 + vert_dif_c_y**2.0 + vert_dif_c_z**2.0 )                 
!Heron's formula
                 A_Heron_S = (norm_a + norm_b + norm_c)/2.0
                 A_Heron   = sqrt(A_Heron_S*(A_Heron_S-norm_a)*(A_Heron_S-norm_b)*(A_Heron_S-norm_c))
                 
!check ok                 write(*,*) vert_face_x(N,nm),vert_face_y(N,nm), vert_face_z(N,nm)
!
!Dot Product - for check purposes                 
                 A_dot = vert_dif_a_x*vert_dif_b_x + vert_dif_a_y*vert_dif_b_y + vert_dif_a_z*vert_dif_b_z
                 ARG = A_dot/(norm_a*norm_b)           
                 theta = acos(ARG)                               !acos(x) x in radians
!Cross Product - parallelogram area - METHOD USED HERE                 
                 Ax_comp = vert_dif_a_y * vert_dif_b_z - vert_dif_a_z * vert_dif_b_y
                 Ay_comp = vert_dif_a_z * vert_dif_b_x - vert_dif_a_x * vert_dif_b_z
                 Az_comp = vert_dif_a_x * vert_dif_b_y - vert_dif_a_y * vert_dif_b_x  
!Magnitude                 
                 A_cross    = sqrt( Ax_comp**2.0 + Ay_comp**2.0 + Az_comp**2.0 )
                 A_cross_dot   = norm_a*norm_b*abs(sin(theta))                    !check purposes: ok A_cross and A_cross_dot give the same results
!check ok                 write(*,*) N, na, A_cross, A_cross_dot, A_Heron*2.0              !check purposes: ok Heron's formula gives the same area                                        
                 Area_f(na) = Area_f(na) + A_cross                                        !cross product parallelogram areas = 2*area
                 A_Heron_f(na)  = A_Heron_f(na) + A_Heron                                 !Heron's area is exact - no need to divide by 2
              end do
              Area_f(na) = Area_f(na)/2.0
               
!              write(*,*) na, Area_f(na)
           endif
        end do
        
!        stop
        
        do N = 1, NF
           Area_t = Area_t + Area_f(N) 
!           write(*,*)'Face =', N, 'Area =', Area_f(N)
        end do
!        write(*,*)
!         write(*,*)'AREA', Area_t
        
        if(na /= NF)then
           write(*,*)'Number of faces is wrong in area calculation'
           stop
        endif
        
!NG   CALCULATE VOLUME
        na = 0
        do N = 1,NN                          !NG   Warning - Loop over possible canditates NOT faces - some do not correspond to a face
           if(nvert_face(N) > 0) then
              na = na + 1
              do nm = 2, nvert_face(N) - 1
                 Ax_c = vert_face_y(N,1) * vert_face_z(N,nm) - vert_face_z(N,1) * vert_face_y(N,nm)
                 Ay_c = vert_face_z(N,1) * vert_face_x(N,nm) - vert_face_x(N,1) * vert_face_z(N,nm)
                 Az_c = vert_face_x(N,1) * vert_face_y(N,nm) - vert_face_y(N,1) * vert_face_x(N,nm)
                 Ax_comp = Ax_c * vert_face_x(N,nm+1)
                 Ay_comp = Ay_c * vert_face_y(N,nm+1)
                 Az_comp = Az_c * vert_face_z(N,nm+1)
                 Vcros_dot = Ax_comp + Ay_comp + Az_comp
                 Vol_f(na) = Vol_f(na) + abs(Vcros_dot)                                                   !Volume partial contributions              
              end do
           endif
        end do
        
        do N = 1, NF
           Vol_t = Vol_t + Vol_f(N) 
!           write(*,*)Vol_f(N) 
        end do 
        
        Vol_t = Vol_t/6.0
!        write(*,*)'VOLUME',Vol_t
        if(na /= NF)then
           write(*,*)'Number of faces is wrong in volume calculation'
           stop
        endif
        
!check        WRITE(*,*)ICP,NF,NV,NE
!Uncomment the following stop for tests for a single particle        
!        stop

!NG   end modified by NG        
        
        IF ( MOD ( NE, 2 ) .NE. 0 ) THEN
           write(*,*)
           WRITE(*,'('' NONINTEGER NUMBER OF EDGES'',I5)') NE
           WRITE(*,'('' WATER MOLECULE INDEX'',I5)')ICP
           WRITE(*,'('' NUMBER OF VERTICES'',I5)')NV
           WRITE(*,'('' NUMBER OF FACES'',I5)')NF
           write(*,*)
           LEULER = .false.
!           STOP

        ENDIF

        NE = NE / 2

        IF ( ( NV - NE + NF ) .NE. 2 ) THEN

           WRITE(*,'('' **** EULER ERROR: DEGENERACY ? **** '')')
           WRITE(*,'('' WATER MOLECULE INDEX'',I5)')ICP
           WRITE(*,'('' NUMBER OF EDGES'',I5)')NE*2
           WRITE(*,'('' NUMBER OF VERTICES'',I5)')NV
           WRITE(*,'('' NUMBER OF FACES'',I5)')NF
           write(*,*)
           LEULER = .false.

        ENDIF

        RETURN
        
END SUBROUTINE WORK
        

END MODULE water_solv
