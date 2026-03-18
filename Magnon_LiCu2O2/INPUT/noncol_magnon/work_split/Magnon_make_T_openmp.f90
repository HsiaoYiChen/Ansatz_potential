program make_Tab
        use omp_lib
        implicit none
        integer :: io, ik,nk,iwan,jwan,nwan,is,i_dim,j_dim,H_dim,nspin,i_mag,i,j,ikq,js
        integer :: period_x,period_y,period_z,N_mag,n_super,i_super,n_site,i_site,j_site
        integer :: n1,n2,n3,n4,w1,w2,w3,w4,j_super,ix,cal_q
        integer :: i_ud,j_ud,ud1,ud2,ud3,ud4
        complex,allocatable     :: Unmk(:,:,:)
        integer,allocatable     :: mag_wan_ind(:,:),K_ind(:,:),ind(:,:)
        real,allocatable        :: band_E(:,:),k_list(:,:),w_list(:)
        real    ::read3(3),Ef,delta!,dE
        character(len=200)      :: tot_HR_path,Hk_com_path,super_H_path,useless_string
        character(len=20)       ::q_char,mag_char,eig_char,i_mag_char,j_mag_char
        complex,allocatable :: Uk(:,:),Ukq(:,:)
        complex,allocatable :: Kab(:,:,:)
        complex,allocatable :: Kab_L(:,:,:,:,:),KL(:,:,:,:)
        complex,allocatable :: Kab_R(:,:,:,:,:),KR(:,:,:,:)
        complex,allocatable :: Kab_LR(:,:,:,:,:,:,:)
        !complex,allocatable :: sum_tmp(:),sum_tmp_R(:),sum_tmp_L(:)
        !complex,allocatable :: sum_tmp_LR(:),G_tmp(:)
        complex,allocatable :: Id_mat(:,:)
        complex(kind(0d0)),allocatable :: tmp_rank2(:,:),inv_rank2(:,:)
        complex(kind(0d0)),allocatable  :: Wread(:,:,:,:,:,:),Wk_rank2(:,:)
        complex,allocatable     :: T_rank2(:,:),T_rank2_tmp(:,:)
        complex,allocatable     ::Wkab(:,:,:)!,Wkab_tmp(:,:)
        complex,allocatable     :: R_ij(:,:,:,:,:),R_ij_local(:,:,:,:)
        complex,allocatable     :: R_ij_Tr(:,:,:)
        complex,allocatable     :: TK(:,:,:,:),KTK(:,:,:,:,:,:)
        complex,allocatable        :: R_mat(:,:),UR_tmp(:,:)
        real,allocatable :: eig_R(:)
        real,allocatable    :: Wab(:,:,:)
        real,allocatable    :: project_vec(:,:),project_component(:)
        complex ::cI,Pauli(2,2,3),Pij!,ctmp
        complex,allocatable ::Pauli_rot(:,:,:,:)
        integer       :: iE,nE,NWF,nEw,red_M_dim,j_mag,n_mag_unit,red_dim_unit!,nwan_mag
        integer         :: print_iE
        real,allocatable :: Ek(:),Ekq(:)
        !integer,allocatable     :: pool(:,:)
        !integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        !integer                 :: npool ,ipool
        !character(len=20)       :: i_mag_char,j_mag_char
        real    :: Emin,Emax,W_lambda
        integer,allocatable     :: cal_unit_wan(:,:),red_ind(:,:),n_gp_mag(:)
        integer             :: n_Wgp,i_Wgp,j_Wgp,n_Wgp_unit,W_gp_nwan_max,i_eig,prt
        integer             :: read_int2(2)
        integer,allocatable :: W_ind(:,:),W_ind_unit(:,:)
        integer,allocatable :: Wgp_mag_nwan(:,:),nwan_per_mag_unit(:)
        logical ::file_exist
        integer :: thread_id

!print*,'start'

        Pauli(:,:,:)= cmplx(0.0d0,0.0d0)
        Pauli(1,2,1)= cmplx(1.0d0,0.0d0)
        Pauli(2,1,1)= cmplx(1.0d0,0.0d0)
        Pauli(1,2,2)=-cmplx(0.0d0,1.0d0)
        Pauli(2,1,2)= cmplx(0.0d0,1.0d0)
        Pauli(1,1,3)= cmplx(1.0d0,0.0d0)
        Pauli(2,2,3)=-cmplx(1.0d0,0.0d0)





        cI=cmplx(0.0,1.0)
        tot_HR_path="../tot_HR"
        Hk_com_path='../Hk_commense'
        super_H_path="../super_H"

open(unit=101,file='../../THE.in',action='read',access='sequential',status='old')
read(101,*)!'(A9,A200)')useless_string,func_dir
read(101,*)!'(A9,A200)')useless_string,nscf_dir
read(101,*)!'(A8,A200)')useless_string,THE_dir
read(101,*)!'(A7,I3)')  useless_string,n_func
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I2)') useless_string, nspin
read(101,*)!'(A3,I3)') useless_string, nx
read(101,*)!'(A3,I3)') useless_string, ny
read(101,*)!'(A3,I3)') useless_string, nz
read(101,*)!'(A6,F8.4)')  useless_string,axial
read(101,'(A9,I3)') useless_string,period_x
read(101,'(A9,I3)') useless_string, period_y
read(101,'(A9,I3)') useless_string, period_z
close(101)


H_dim=nwan*2/nspin*period_x*period_y*period_z
n_super=period_x*period_y*period_z


open(unit=102,file='./mag.in',access='sequential',action='read',status='old')
  read(102,'(A3,F8.4)') useless_string, Ef
  read(102,'(A6,F8.4)') useless_string, delta
  read(102,'(A5,F12.6)') useless_string, Emin
  read(102,'(A5,F12.6)') useless_string, Emax
  read(102,'(A3,I4)') useless_string, nE
  read(102,'(A6,I5)') useless_string,cal_q
  read(102,*)!------------------------------------------
  read(102,'(A11,I3)') useless_string, n_mag_unit
  allocate(nwan_per_mag_unit(n_mag_unit))
  nwan_per_mag_unit(:)=0
  do i_mag=1,N_mag_unit
        read(102,*) nwan_per_mag_unit(i_mag)
        !print*,i_mag,nwan_per_mag_unit(i_mag)
  enddo
  allocate(cal_unit_wan(N_mag_unit,maxval(nwan_per_mag_unit(:))))
  read(102,*)!------------------------------------------
  do i_mag=1,N_mag_unit
  do iwan=1,nwan_per_mag_unit(i_mag)
        read(102,*)  cal_unit_wan(i_mag,iwan)
        !print*,i_mag,iwan,cal_unit_wan(i_mag,iwan)
  enddo
  enddo
close(102)

  call int2char(cal_q,q_char)
  red_dim_unit=0
  do i_mag=1,n_mag_unit
        red_dim_unit=red_dim_unit+nwan_per_mag_unit(i_mag)
  enddo





open(unit=101,file='./W.in',action='read',access='sequential',status='old')
        read(101,'(A9,F12.8)') useless_string, W_lambda
        read(101,'(A6,I3)') useless_string, n_Wgp_unit
        allocate(n_gp_mag(n_Wgp_unit))
        read(101,*) useless_string, n_gp_mag(:)
        if (.not. sum(n_gp_mag(:))  .eq.  n_mag_unit) then
                print*,"n_Wgp*n_gp_mag should eq n_mag_unit"
                stop
        endif
        allocate(Wgp_mag_nwan(n_Wgp_unit,maxval(n_gp_mag(:))))
        do i=1,sum(n_gp_mag(:))
                read(101,*) i_Wgp,i_mag,j  !i_Wgp_th,i_mag_th,j_number_of_wan
                Wgp_mag_nwan(i_Wgp,i_mag)=j
        enddo
        read(101,*)


        allocate(W_ind_unit(red_dim_unit,2))
        i_dim=1
        do i_mag=1,n_mag_unit
        do iwan=1,nwan_per_mag_unit(i_mag)
                read(101,*)  read_int2(:)
                W_ind_unit(i_dim,:)=read_int2(:)
                i_dim=i_dim+1
        enddo
        enddo
close(101)

n_mag=n_super*N_mag_unit
n_Wgp=n_super*n_Wgp_unit

allocate(W_ind(n_super*red_dim_unit,2))


do i_super=1,n_super
do i_dim=1,red_dim_unit
        W_ind(i_dim+(i_super-1)*red_dim_unit,1)=W_ind_unit(i_dim,1)+(i_super-1)*n_Wgp_unit
        W_ind(i_dim+(i_super-1)*red_dim_unit,2)=W_ind_unit(i_dim,2)
enddo
enddo

  W_gp_nwan_max=0
  do i_Wgp=1,n_Wgp_unit
          if (sum(Wgp_mag_nwan(i_Wgp,:)) .gt. W_gp_nwan_max) W_gp_nwan_max=sum(Wgp_mag_nwan(i_Wgp,:))
  enddo

 ! print*,'W_gp_nwan_max = ', W_gp_nwan_max

  allocate(Wread(W_gp_nwan_max*2,W_gp_nwan_max*2,&
                &W_gp_nwan_max*2,W_gp_nwan_max*2,1,n_Wgp_unit))


  call setup()
  call find_axis()

!  print*,'read The.in done'


  allocate(red_ind(n_mag,maxval(nwan_per_mag_unit(:))))
  do i_super=1,n_super
        do i_mag=1,n_mag_unit
                i_site=i_mag+(i_super-1)*n_mag_unit
                do iwan=1,nwan_per_mag_unit(i_mag)
                        red_ind(i_site,iwan)=cal_unit_wan(i_mag,iwan)+(i_super-1)*nwan
                enddo
        enddo
  enddo

  red_M_dim=0
  do i_mag=1,n_mag_unit
        red_M_dim=red_M_dim+nwan_per_mag_unit(i_mag)**2*2*2
  enddo
  red_M_dim=red_M_dim*n_super

 !print*,'red_M_dim=',red_M_dim
  
 allocate(K_ind(red_M_dim,5))
  K_ind(:,:)=0
  i_dim=1
  do i_super=1,n_super
  do i_mag=1,n_mag_unit
  do i_ud=1,2
  do iwan=1,nwan_per_mag_unit(i_mag)
  do j_ud=1,2
  do jwan=1,nwan_per_mag_unit(i_mag)
        K_ind(i_dim,:)=(/i_mag+(i_super-1)*n_mag_unit,i_ud,iwan,j_ud,jwan/)
        i_dim=i_dim+1
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo



  allocate(Wab(red_M_dim,red_M_dim,1))
  allocate(Id_mat(red_M_dim,red_M_dim))
  allocate(tmp_rank2(red_M_dim,red_M_dim))
  allocate(inv_rank2(red_M_dim,red_M_dim))
  allocate(Wk_rank2(red_M_dim,red_M_dim))
  allocate(T_rank2(red_M_dim,red_M_dim))
  allocate(T_rank2_tmp(red_M_dim,red_M_dim))

 ! print*,'index done'
  iE=0
  open(unit=103,file='./data/w_list', action='read',access='sequential',status='old')
  do
          read(103,'(E16.8)',iostat=io) read3(1)
          if (io/=0) exit
          iE=iE+1
  enddo
  close(103)

  nE=iE
  allocate(w_list(nE))
  allocate(R_ij(n_mag,n_mag,3,3,nE))
  allocate(R_ij_local(n_mag,n_mag,3,3))
  allocate(R_mat(n_mag,n_mag))
  allocate(UR_tmp(n_mag,n_mag))
  allocate(eig_R(n_mag))
  allocate(R_ij_Tr(3,3,nE))
  iE=1
  open(unit=103,file='./data/w_list' ,action='read',access='sequential',status='old')
  do
          read(103,'(E16.8)',iostat=io) read3(1)
          if (io/=0) exit
          w_list(iE)=read3(1)
          iE=iE+1
  enddo
  close(103)


        allocate(Kab   (    red_M_dim,red_M_dim,nE))
        allocate(Kab_L (2,2,n_mag,    red_M_dim,nE),KL(2,2,n_mag,    red_M_dim))
        allocate(Kab_R (2,2,red_M_dim,    n_mag,nE),KR(2,2,red_M_dim,    n_mag))
        allocate(Kab_LR(2,2,2,2,n_mag,        n_mag,nE))
        allocate(Wkab  (    red_M_dim,red_M_dim,nE))
   !     print*,'allocate done'

  Kab   (:,:,:)         =cmplx(0.0,0.0)
  Kab_L (:,:,:,:,:)     =cmplx(0.0,0.0)
  Kab_R (:,:,:,:,:)     =cmplx(0.0,0.0)
  Kab_LR(:,:,:,:,:,:,:) =cmplx(0.0,0.0)

        open(unit=102,file='data/Kab_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab
        close(102)



        open(unit=102,file='data/Kab_L_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_L
        close(102)

        open(unit=102,file='data/Kab_R_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_R
        close(102)

        open(unit=102,file='data/Kab_LR_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_LR
        close(102)

        print*,'read done'





 !       open(unit=102 ,file='cmap11',action='write',access='sequential',status='replace')
 !       do iE=1,nE
 !               write(102,'(F12.6,16E18.6)') w_list(iE),real(Kab(1,1,:,:,iE))*6.28*2.03
 !       enddo
 !       close(102)

 !       open(unit=102 ,file='cmap12',action='write',access='sequential',status='replace')
 !       do iE=1,nE
 !               write(102,'(F12.6,16E18.6)') w_list(iE),real(Kab(1,2,:,:,iE))*6.28*2.03
 !       enddo
 !       close(102)
 !       open(unit=102 ,file='cmap21',action='write',access='sequential',status='replace')
 !       do iE=1,nE
 !               write(102,'(F12.6,16E18.6)') w_list(iE),real(Kab(2,1,:,:,iE))*6.28*2.03
 !       enddo
 !       close(102)
 !       open(unit=102 ,file='cmap22',action='write',access='sequential',status='replace')
 !       do iE=1,nE
 !               write(102,'(F12.6,16E18.6)') w_list(iE),real(Kab(2,2,:,:,iE))*6.28*2.03
 !       enddo
 !       close(102)




        R_ij(:,:,:,:,:)=0.0
        do i=1,2
        do j=1,2
        do ud1=1,2
        do ud2=1,2
        do ud3=1,2
        do ud4=1,2
        do i_mag=1,n_mag
        do j_mag=1,n_mag
cycle

        R_ij(i_mag,j_mag,i,j,:)=R_ij(i_mag,j_mag,i,j,:)+&
                &Pauli(ud1,ud2,i)*Kab_LR(ud1,ud2,ud3,ud4,i_mag,j_mag,:)*Pauli(ud4,ud3,j)



!R_ij(i_mag,j_mag,i,j,:)=R_ij(i_mag,j_mag,i,j,:)+Pauli(js,is,j)*Pauli(is,js,i)*&
!                        &Kab_LR(is,js,i_mag,j_mag,:)

 
 !
           !     if (i_mag .eq.1 .and. j_mag .eq. 1) then
!                print*,Pauli(js,is,j)*Pauli(is,js,i)*&
!                        &Kab_LR(is,js,i_mag,j_mag,138),Pauli(js,is,j)*Pauli(is,js,i)*&
!                        &Kab_LR(is,js,i_mag,j_mag,1864)

            !    print*,Kab_LR(is,js,i_mag,j_mag,138),Kab_LR(is,js,i_mag,j_mag,1864)


             !   endif


        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

!goto 500

!stop

        Id_mat(:,:)=cmplx(0.0,0.0)
        do iwan=1,red_M_dim
                Id_mat(iwan,iwan)=cmplx(1.0,0.0)
        enddo


        inquire( file='./data/Wab', exist=file_exist )
        if (file_exist) then
                print*,'Wab found, read it'
                open(unit=102,file='./data/Wab',action='read'&
                        &,access='sequential',status='old',form='unformatted')
                        read(102) Wab(:,:,:)
                close(102)
        else


        do  i_Wgp=1,N_Wgp_unit
                call int2char(i_Wgp,mag_char)
                open(unit=102,file='../dir-intW_'//trim(mag_char)//'/dat.WR0',action='read'&
                                &,access='sequential',status='old',form='unformatted')

                                read(102) NWF,nEw
                                !if(.not. NWF .eq. 2*nwan)  then
                                !        print*,'wrong nwan', NWF,nwan
                                !        stop
                                !endif
                                if(nEw .gt. 1 ) then
                                        print*,'only static screening is fine'
                                        stop
                                endif
                                read(102)
                                read(102) Wread(1:NWF,1:NWF,1:NWF,1:NWF,1,i_Wgp)
                        close(102)
        enddo
 

                        
     !   print*,'convert to Wab'
     !   Wab(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
        Wab(:,:,:)=cmplx(0.0d0, 0.0d0, kind=kind(1.0d0))
        do i_dim=1,red_M_dim
        do j_dim=1,red_M_dim

                
                if(.not. K_ind(i_dim,2) .eq. K_ind(j_dim,2) ) cycle
                if(.not. K_ind(i_dim,4) .eq. K_ind(j_dim,4) ) cycle
                if (.not. (K_ind(i_dim,1)-1)/n_mag_unit .eq.(K_ind(j_dim,1)-1)/n_mag_unit ) cycle


                i_mag=mod(K_ind(i_dim,1)-1,n_mag_unit)+1
                j_mag=mod(K_ind(j_dim,1)-1,n_mag_unit)+1




                if (i_mag .eq. 1) then
                        n1=K_ind(i_dim,3)
                        n4=K_ind(i_dim,5)
                elseif (i_mag .eq. 2) then
                        n1=K_ind(i_dim,3)+nwan_per_mag_unit(1)
                        n4=K_ind(i_dim,5)+nwan_per_mag_unit(1)
                else
                        n1=K_ind(i_dim,3)+sum(nwan_per_mag_unit(1:(i_mag-1)))
                        n4=K_ind(i_dim,5)+sum(nwan_per_mag_unit(1:(i_mag-1)))
                endif


                if (j_mag .eq. 1) then
                        n2=K_ind(j_dim,3)
                        n3=K_ind(j_dim,5)
                elseif (j_mag .eq. 2) then
                        n2=K_ind(j_dim,3)+nwan_per_mag_unit(1)
                        n3=K_ind(j_dim,5)+nwan_per_mag_unit(1)
                else
                        n2=K_ind(j_dim,3)+sum(nwan_per_mag_unit(1:(j_mag-1)))
                        n3=K_ind(j_dim,5)+sum(nwan_per_mag_unit(1:(j_mag-1)))
                endif


                i_Wgp=W_ind_unit(n1,1)
                j_Wgp=W_ind_unit(n2,1)

                if (.not. i_Wgp .eq. j_Wgp) cycle

               ! do is=1,2
               ! do js=1,2

                        !print*,n1,n2,n3,n4

                        w1=K_ind(i_dim,2)+2*(W_ind_unit(n1,2)-1)
                        w2=K_ind(j_dim,2)+2*(W_ind_unit(n2,2)-1)
                        w3=K_ind(j_dim,4)+2*(W_ind_unit(n3,2)-1)
                        w4=K_ind(i_dim,4)+2*(W_ind_unit(n4,2)-1)
 
                        
 !                       print*,i_dim,j_dim,real(Wread(w1,w2,w4,w3,:,i_Wgp))

                        Wab(i_dim,j_dim,:)=real(Wread(w1,w2,w4,w3,:,i_Wgp))

                         if (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 5.7) then
                                Wab(i_dim,j_dim,:)=5.7977!real(Wread(w1,w2,w4,w3,:,i_Wgp)) d-d
                        elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 4.0) then 
                                Wab(i_dim,j_dim,:)=4.21 !p-p
                        elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 2.0) then
                                Wab(i_dim,j_dim,:)=2.61 !px-py
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 1.008) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!1.013 p-d
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.7) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.745  p1-p2
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.335) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.3425 p1-p3
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.313) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.321 d1-d2 p1-p4
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.28) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.284 d1-p2
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.23) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.235 d-d-d-p
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.175) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.181
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.09) then
                       !         !print*,w1,w2,w3,w4
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.0937 !d-p-d-p 1515 1526
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.057) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.023
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.049) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!0.015
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.027) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.0188) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!-0.015
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.01) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!-0.023
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. -0.04) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!-0.0937
                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. -0.138) then
                       !         Wab(is,js,i_dim,j_dim,:)=0.0!-0.181



                        


                       ! elseif (real(Wread(w1,w2,w4,w3,1,i_Wgp)) .gt. 0.0 .and. real(Wread(w1,w2,w4,w3,1,i_Wgp)).lt. 0.2 ) then
                        else
                        Wab(i_dim,j_dim,:)=0.0!-0.235

                        endif




                !enddo !js
                !enddo !is
        enddo !j_dim
        enddo !i_dim


                open(unit=102,file='./data/Wab',action='write'&
                        &,access='append',status='replace',form='unformatted')
                        write(102) Wab(:,:,:)
                close(102)



        endif

   !     print*,'read done'

                        Wab(:,:,:)=Wab(:,:,:)*W_lambda
        print*,'start parallel'
        !$OMP PARALLEL PRIVATE(thread_id,iE,i_dim,j_dim, TK,KTK,&
        !$OMP& is,js,i,j,i_mag,j_mag,Pij,KL,KR,Wk_rank2,tmp_rank2,inv_rank2,T_rank2_tmp,R_ij_local) &
        !$OMP& SHARED(Kab_L,Kab_R,&
        !$OMP& red_M_dim,nE,n_mag,N_mag_unit,cI,n_super,Wab,Kab,Id_mat,Pauli,R_ij)&
        !$OMP& DEFAULT(NONE)
        
        allocate(TK(red_M_dim,n_mag,2,2),KTK(2,2,n_mag,n_mag,2,2))
        
        !$OMP DO
        do iE=1,nE
        !        allocate(TK(red_M_dim,n_mag),KTK(n_mag,n_mag))
                R_ij_local(:,:,:,:)=cmplx(0.0,0.0)
                !print*,"Thread", thread_id, "iE", iE
                thread_id = omp_get_thread_num()
                !print*,"Thread", thread_id, "iE", iE
                ! if (mod(iE,10).eq. 1)     print *, "Thread", thread_id, "iE", iE
                !print*,is,js
                KL(:,:,:,:)=Kab_L(:,:,:,:,iE)
                KR(:,:,:,:)=Kab_R(:,:,:,:,iE)

                Wk_rank2(:,:)=matmul(Wab(:,:,1),Kab(:,:,iE))
                !print*,sum(Kab(:,:,iE))
                tmp_rank2(:,:)=Id_mat(:,:)-Wk_rank2(:,:)
                call invA(tmp_rank2, red_M_dim, inv_rank2)
                T_rank2_tmp(:,:)=matmul(inv_rank2(:,:),Wk_rank2(:,:))

                do ud3=1,2
                do ud4=1,2
                TK(:,:,ud3,ud4)=matmul(T_rank2_tmp(:,:),KR(ud3,ud4,:,:))
                enddo
                enddo
                
                
                do ud1=1,2
                do ud2=1,2
                do ud3=1,2
                do ud4=1,2
                        KTK(ud1,ud2,:,:,ud3,ud4)=matmul(KL(ud1,ud2,:,:),TK(:,:,ud3,ud4))
                enddo
                enddo
                enddo
                enddo




                do i=1,2
                do j=1,2
                        do ud1=1,2
                        do ud2=1,2
                        do ud3=1,2
                        do ud4=1,2
                                Pij=Pauli(ud1,ud2,i)*Pauli(ud4,ud3,j)
                                R_ij_local(:,:,i,j)=R_ij_local(:,:,i,j)+KTK(ud1,ud2,:,:,ud3,ud4)*Pij                               
                        enddo
                        enddo
                        enddo
                        enddo        
                enddo
                enddo
         !       deallocate(TK, KTK)
                R_ij(:,:,:,:,iE)=R_ij_local(:,:,:,:)
        enddo
        !$OMP END  DO
        !$OMP END PARALLEL




500 continue 

print_iE=156

allocate(project_vec(1,n_mag))
allocate(project_component(1))


project_vec(1,:)=1
!project_vec(1,1:n_mag/2)=(/1,1,1,-1,-1,-1,1,1,1,-1,-1,-1/)
!project_vec(1,n_mag/2+1:n_mag)=1*(/1,1,1,-1,-1,-1,1,1,1,-1,-1,-1/)


!project_vec(2,1:n_mag/2)=(/-1,-1,-1,1,1,1,1,1,1,-1,-1,-1/)
!project_vec(2,n_mag/2+1:n_mag)=(/-1,-1,-1,1,1,1,1,1,1,-1,-1,-1/)


!project_vec(3,1:n_mag/2)=(/1,1,1,1,1,1,-1,-1,-1,-1,-1,-1/)
!project_vec(3,n_mag/2+1:n_mag)=(/1,1,1,1,1,1,-1,-1,-1,-1,-1,-1/)


!project_vec(4,1:n_mag/2)=(/1,1,1,1,1,1,1,1,1,1,1,1/)
!project_vec(4,n_mag/2+1:n_mag)=(/1,1,1,1,1,1,1,1,1,1,1,1/)


!project_vec(5,1:n_mag/2)=(/1,1,1,-1,-1,-1,1,1,1,-1,-1,-1/)
!project_vec(5,n_mag/2+1:n_mag)=-1*(/1,1,1,-1,-1,-1,1,1,1,-1,-1,-1/)


!project_vec(6,1:n_mag/2)=(/-1,-1,-1,1,1,1,1,1,1,-1,-1,-1/)
!project_vec(6,n_mag/2+1:n_mag)=-1*(/-1,-1,-1,1,1,1,1,1,1,-1,-1,-1/)


!project_vec(7,1:n_mag/2)=(/1,1,1,1,1,1,-1,-1,-1,-1,-1,-1/)
!project_vec(7,n_mag/2+1:n_mag)=-1*(/1,1,1,1,1,1,-1,-1,-1,-1,-1,-1/)


!project_vec(8,1:n_mag/2)=(/1,1,1,1,1,1,1,1,1,1,1,1/)
!project_vec(8,n_mag/2+1:n_mag)=-1*(/1,1,1,1,1,1,1,1,1,1,1,1/)





        R_ij_Tr(:,:,:)=cmplx(0.0,0.0)
        do i_mag=1,n_mag
        do j_mag=1,n_mag

                call int2char(i_mag,i_mag_char)
                call int2char(j_mag,j_mag_char)

         !       open(unit=102,file='./Rij_'//trim(i_mag_char)//'_'//trim(j_mag_char)//'_q'//trim(q_char),action='write'&
         !               &,access='append',status='replace')
         !       do i=1,nE
         !    write(102,'(F12.8,2E16.8)') w_list(i),&
         !    &real (R_ij(i_mag,j_mag,1,1,i)+R_ij(i_mag,j_mag,2,2,i)+cI*(R_ij(i_mag,j_mag,1,2,i)-R_ij(i_mag,j_mag,2,1,i))),&
         !    &aimag(R_ij(i_mag,j_mag,1,1,i)+R_ij(i_mag,j_mag,2,2,i)+cI*(R_ij(i_mag,j_mag,1,2,i)-R_ij(i_mag,j_mag,2,1,i)))
         !       enddo
         !       close(102)

                R_mat(i_mag,j_mag)=(R_ij(i_mag,j_mag,1,1,print_iE)+R_ij(i_mag,j_mag,2,2,print_iE)&
                        &+cI*(R_ij(i_mag,j_mag,1,2,print_iE)-R_ij(i_mag,j_mag,2,1,print_iE)))

        enddo
                R_ij_Tr(:,:,:)=R_ij_Tr(:,:,:)+R_ij(i_mag,i_mag,:,:,:)

        enddo
        
        R_mat(:,:)=R_mat(:,:)/cI
        R_mat(:,:)=0.5*(R_mat(:,:)+transpose(conjg(R_mat(:,:))))






        prt=0
        do i_eig=n_mag,n_mag-10,-1
                call int2char(i_eig,eig_char)

                open(unit=102,file='./Rij_trace_q'//trim(q_char)//'_'//trim(eig_char),action='write'&
                        &,access='append',status='replace')
                do i=1,nE

    !                    R_mat(:,:)=(R_ij(:,:,1,1,i)+R_ij(:,:,2,2,i)&
    !                    &+cI*(R_ij(:,:,1,2,i)-R_ij(:,:,2,1,i)))                

    !                    R_mat(:,:)=R_mat(:,:)/cI
    !                    R_mat(:,:)=0.5*(R_mat(:,:)+transpose(conjg(R_mat(:,:))))


                
     !                   call diasym(R_mat,eig_R,n_mag,UR_tmp)

                        project_component(1)=abs(sum(UR_tmp(:,i_eig)*project_vec(1,:)))
                        !project_component(2)=abs(sum(UR_tmp(:,i_eig)*project_vec(2,:)))
                        !project_component(3)=abs(sum(UR_tmp(:,i_eig)*project_vec(3,:)))
                        !project_component(4)=abs(sum(UR_tmp(:,i_eig)*project_vec(4,:)))
                        !project_component(5)=abs(sum(UR_tmp(:,i_eig)*project_vec(5,:)))
                        !project_component(6)=abs(sum(UR_tmp(:,i_eig)*project_vec(6,:)))
                        !project_component(7)=abs(sum(UR_tmp(:,i_eig)*project_vec(7,:)))
                        !project_component(8)=abs(sum(UR_tmp(:,i_eig)*project_vec(8,:)))



                        write(102,'(F12.8,2E16.8,2F16.8)') w_list(i),&
                                &real (R_ij_Tr(1,1,i)+R_ij_Tr(2,2,i)+cI*(R_ij_Tr(1,2,i)-R_ij_Tr(2,1,i))),&
                                &aimag(R_ij_Tr(1,1,i)+R_ij_Tr(2,2,i)+cI*(R_ij_Tr(1,2,i)-R_ij_Tr(2,1,i))),&
                                &eig_R(i_eig),project_component(1)
                        if (aimag(R_ij_Tr(1,1,i)+R_ij_Tr(2,2,i)+cI*(R_ij_Tr(1,2,i)-R_ij_Tr(2,1,i))) .gt. 200 .and. prt .eq. 0) then
                                prt=1
                                print*,W_lambda
                        endif
                enddo
                close(102)
        enddo


        call diasym(R_mat,eig_R,n_mag,UR_tmp)

        print*,eig_R

!print*,''
!        print*,real(UR_tmp(:,1))
!print*,''
!        print*,real(UR_tmp(:,2))
!print*,''

!print*,real(UR_tmp(:,3))
!print*,''

!print*,''
print'(8F8.3)',(UR_tmp(1+3*(i-1),12),i=1,4)
!print*,''
!print'(8F8.3)',(UR_tmp(1+3*(i-1),139),i=1,48)
!print*,''
!print'(8F8.3)',(UR_tmp(1+3*(i-1),140),i=1,48)
!print*,''
!print'(8F8.3)',(UR_tmp(1+3*(i-1),141),i=1,48)
!print*,''
!print'(8F8.3)',(UR_tmp(1+3*(i-1),142),i=1,48)
!print*,''
!print'(8F8.3)',(UR_tmp(1+3*(i-1),143),i=1,48)
!print*,''
!print'(8F8.3)',(UR_tmp(1+3*(i-1),144),i=1,48)
!
!        print*,real(UR_tmp(:,11))
!print*,''
!        print*,real(UR_tmp(:,12))


contains



        subroutine setup()
                integer ::read3_int(3)



        ik=0
        open(unit=102,file='../Hk_commense/k_list',action='read',access='sequential',status='old')
        do
                read(102,'(3F12.8)',iostat=io) read3(:)
                if (io/=0) exit
                ik=ik+1
        enddo
        close(102)
        nk=ik
        allocate(k_list(nk,3))

        ik=1
        open(unit=102,file='../Hk_commense/k_list',action='read',access='sequential',status='old')
        do
                read(102,'(3F12.8)',iostat=io) read3(:)
                if (io/=0) exit
                k_list(ik,:)=read3(:)
                ik=ik+1
        enddo
        close(102)

        allocate(Unmk(H_dim,H_dim,nk))
        allocate(Uk(H_dim,H_dim))
        allocate(Ukq(H_dim,H_dim))
        allocate(band_E(H_dim,nk))
        allocate(Ek(H_dim),Ekq(H_dim))

        open(unit=102,file='../tot_HR/U_mat',action='read',access='sequential',status='old',form='unformatted')
                read(102) read3_int(:)
                if (.not. read3_int(1) .eq. H_dim)  then
                        print*,'Umat dim 1 wrong'
                        stop
                endif
                if (.not. read3_int(2) .eq. H_dim)  then
                        print*,'Umat dim 2 wrong'
                        stop
                endif
                if (.not. read3_int(3) .eq. nk)  then
                        print*,'Umat dim 3 wrong'
                        stop
                endif
                read(102) Unmk
        close(102)

        open(unit=102,file='../tot_HR/band_E',action='read',access='sequential',status='old')
        do
                read(102,'(2I8,F16.8)',iostat=io) read3_int(1:2), read3(1)
                if(io/=0) exit
                band_E(read3_int(2),read3_int(1))=read3(1)
        enddo
        close(102)

        end subroutine setup



subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )    WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char




subroutine invA(a,n,inv_a)
 implicit none
 integer,intent(in) :: n
 complex(kind(0d0)),intent(in)  :: a(n,n)
 complex(kind(0d0))    :: a_tmp(n,n)
 complex(kind(0d0)),intent(out)  :: inv_a(n,n)
 complex(kind(0d0)), allocatable :: work(:)
 integer, dimension(n) :: ipiv
 integer        :: info,lwork

if (n <= 0) then
        print *, "Error: n must be positive"
        stop
    end if
 a_tmp(:,:)=a(:,:)
! print*,'call zgetrf'
 call zgetrf(n, n, a_tmp, n, ipiv, info)
! print*,'call zgetrf done'
    if (info /= 0) then
        print *, "Error in LU factorization: ", info
 !       print*,a(:,:)
        stop
    end if
 lwork = -1
! print*,'allocate work1'
 allocate(work(1))
! print*,'call zgetri'
 call zgetri(n, a_tmp, n, ipiv, work, lwork, info)
 if (info /= 0) then
        print *, "Error in workspace query: ", info
        stop
    end if
 lwork = int(work(1))  ! Optimal workspace size
! print*,lwork
 deallocate(work)
 !print *, "Optimal workspace size: ", lwork
 allocate(work(lwork))
 call zgetri(n, a_tmp, n, ipiv, work, lwork, info)
    if (info /= 0) then
        print *, "Error in matrix inversion: ", info
        stop
    end if
 inv_a(:,:)=a_tmp(:,:)
 deallocate(work)
end subroutine invA



subroutine diasym(a,eig,n,D)
 implicit none

 integer,intent(in)  :: n
 complex,intent(in)  :: a(n,n)
 complex              :: a_tmp(n,n)
 complex,intent(out)  :: D(n,n)

 integer  :: l,inf
 complex :: Lwork(2*n-1)
 real,intent(out)      :: eig(n)
 real                :: Rwork(3*n-2)

 a_tmp(:,:)=a(:,:)

 l=2*n-1
 call cheev('V','U',n,a_tmp,n,eig,Lwork,l,Rwork,inf)
 D(:,:)=a_tmp(:,:)



end subroutine diasym





subroutine find_axis()
        integer :: irept,iwan_sp
        integer ::read3_int(3)
        complex(kind(0d0)) :: H_R0(H_dim,H_dim)
        complex::       M_tmp(2,2)
        character(len=20) :: irept_char
        real    :: st,ct,sp,cp,s1,s2,s3!,eig2(2)
        complex :: local_U(2,2)!,Ux(2,2)
        real    :: Mz(3,3),My(3,3),rot_M(3,3)

        allocate(Pauli_rot(N_mag,2,2,3))

        irept=1
        open(unit=101,file='../super_H/ir_vec_HYC',action='read',access='sequential')
        do
        read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                if (all(read3_int(:) .eq. (/0,0,0/))) then
                        exit
                endif
                irept=irept+1
        enddo
        close(101)

      !  print*, irept
        call int2char(irept,irept_char)
        open(unit=101,file='../H_sk/data/Hsk_'//trim(irept_char),action='read'&
                        &,access='sequential',form='unformatted')
                        read(101) H_R0(:,:)
        close(101)

        !print*,sum(abs(H_R0(:,:)))

!        M_tmp(:,:)=Pauli(:,:,3)
!        call diasym(M_tmp(:,:),eig2,2,Uz(:,:))


        do i_super=1,n_super
        do i_mag=1,n_mag_unit
                M_tmp(:,:)=cmplx(0.0,0.0)
                local_U(:,:)=cmplx(0.0,0.0)
                do iwan=1,nwan_per_mag_unit(i_mag)
                        iwan_sp=cal_unit_wan(i_mag,iwan)+(i_super-1)*nwan
                        !print*,2*iwan_sp-1,2*iwan_sp, H_dim 
                        M_tmp(:,:)= M_tmp(:,:)+H_R0(2*iwan_sp-1:2*iwan_sp,2*iwan_sp-1:2*iwan_sp)
        !                print*,i_mag,H_R0(2*iwan_sp-1:2*iwan_sp,2*iwan_sp-1:2*iwan_sp)
                enddo
                !print*,i_mag,M_tmp(:,:)

                M_tmp(:,:)=-M_tmp(:,:)
                M_tmp(1,1)=(M_tmp(1,1)-M_tmp(2,2))/2.0
                M_tmp(2,2)=-M_tmp(1,1)
                M_tmp(1,2)=conjg(M_tmp(2,1))


                s3=M_tmp(1,1)
                s1=real(M_tmp(2,1))
                s2=aimag(M_tmp(2,1))

                st=((s1**2+s2**2)/(s1**2+s2**2+s3**2))**(0.5)
                ct=s3/(s1**2+s2**2+s3**2)**(0.5)

                sp=s2/(s1**2+s2**2)**(0.5)
                cp=s1/(s1**2+s2**2)**(0.5)


                My(1,:)=(/ct ,0.0, st/)
                My(2,:)=(/0.0,1.0,0.0/)
                My(3,:)=(/-st,0.0, ct/)
                
                Mz(1,:)=(/cp ,-sp,0.0/)
                Mz(2,:)=(/sp,  cp,0.0/)
                Mz(3,:)=(/0.0,0.0,1.0/)


                rot_M(:,:)=matmul(Mz(:,:),My(:,:))

 !               print*,sin(10.0*3.14159*2/360)

                !print*,s1,s2,s3
 !               print*,ct,st,cp,sp
                !print'(3F12.4)',rot_M(:,1)
                !print'(3F12.4)',rot_M(:,2)
                !print'(3F12.4)',rot_M(:,3)

                do ix=1,3
                        Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)=&
                                &rot_M(1,ix)*Pauli(:,:,1)+rot_M(2,ix)*Pauli(:,:,2)+rot_M(3,ix)*Pauli(:,:,3)
 !                       print'(I8,8F12.4)',i_mag,Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)
                enddo


      !          call diasym(M_tmp,eig2,2,local_U(:,:))

     !           print'(I8,10F12.4)',i_mag, M_tmp(:,:)
     !           print'(I8,10F12.4)',i_mag,local_U(:,1),local_U(:,2),eig2(:)
     !           do ix=1,3
     !                   M_tmp(:,:)=cmplx(0.0,0.0)        
     !                   M_tmp(:,:)=matmul(Pauli(:,:,ix),transpose(conjg(local_U(:,:))))
     !                   Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)=matmul(local_U(:,:), M_tmp(:,:))
                
     !                   print'(I8,8F12.4)',i_mag,Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)
     !           enddo



        
        enddo
        enddo



do i_super=1,n_super
        do i_mag=1,n_mag_unit

 !               Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,:)=Pauli_rot(1,:,:,:)

enddo
        enddo



end subroutine find_axis




      end program
