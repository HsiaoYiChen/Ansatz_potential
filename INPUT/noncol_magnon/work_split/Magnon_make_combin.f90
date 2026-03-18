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
        character(len=20)       :: i_char,j_char
        complex,allocatable :: Uk(:,:),Ukq(:,:)
        complex,allocatable :: Kab(:,:,:)
        complex,allocatable :: Kab_L(:,:,:,:,:),KL(:,:,:,:)
        complex,allocatable :: Kab_R(:,:,:,:,:),KR(:,:,:,:)
        complex,allocatable :: Kab_LR(:,:,:,:,:,:,:)
        complex,allocatable :: Kab_read(:,:,:)
        complex,allocatable :: Kab_L_read(:,:,:,:,:)
        complex,allocatable :: Kab_R_read(:,:,:,:,:)
        complex,allocatable :: Kab_LR_read(:,:,:,:,:,:,:)

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


        call EXECUTE_COMMAND_LINE("mkdir -p  data")
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
  open(unit=103,file='../work_split_1_4/data/w_list', action='read',access='sequential',status='old')
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
  open(unit=103,file='../work_split_1_4/data/w_list' ,action='read',access='sequential',status='old')
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
        allocate(Kab_read   (    red_M_dim,red_M_dim,nE))
        allocate(Kab_L_read (2,2,n_mag,    red_M_dim,nE))
        allocate(Kab_R_read (2,2,red_M_dim,    n_mag,nE))
        allocate(Kab_LR_read(2,2,2,2,n_mag,        n_mag,nE))




  Kab   (:,:,:)         =cmplx(0.0,0.0)
  Kab_L (:,:,:,:,:)     =cmplx(0.0,0.0)
  Kab_R (:,:,:,:,:)     =cmplx(0.0,0.0)
  Kab_LR(:,:,:,:,:,:,:) =cmplx(0.0,0.0)


        do i =1,168,4
                j=i+3
                call int2char(i,i_char)
                call int2char(j,j_char)


        open(unit=102,file='../work_split_'//trim(i_char)//'_'//trim(j_char)//'/data/Kab_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_read
        close(102)



        open(unit=102,file='../work_split_'//trim(i_char)//'_'//trim(j_char)//'/data/Kab_L_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_L_read
        close(102)

        open(unit=102,file='../work_split_'//trim(i_char)//'_'//trim(j_char)//'/data/Kab_R_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_R_read
        close(102)

        open(unit=102,file='../work_split_'//trim(i_char)//'_'//trim(j_char)//'/data/Kab_LR_q'//trim(q_char),action='read'&
                &,access='sequential',status='old',form='unformatted')
        read(102) Kab_LR_read
        close(102)

        print*,'read done',i,j

        Kab   (:,:,:) =Kab   (:,:,:) +Kab_read   (:,:,:) 
        Kab_L (:,:,:,:,:)     =Kab_L (:,:,:,:,:)     +Kab_L_read (:,:,:,:,:)     
        Kab_R (:,:,:,:,:)=Kab_R (:,:,:,:,:)+Kab_R_read (:,:,:,:,:)
        Kab_LR(:,:,:,:,:,:,:)=Kab_LR(:,:,:,:,:,:,:)+Kab_LR_read(:,:,:,:,:,:,:)
        enddo



        open(unit=102,file='data/Kab_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab
        close(102)

        open(unit=102,file='data/Kab_L_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab_L
        close(102)
        open(unit=102,file='data/Kab_R_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab_R
        close(102)

        open(unit=102,file='data/Kab_LR_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab_LR
        close(102)


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

       ! allocate(Unmk(H_dim,H_dim,nk))
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
        !        read(102) Unmk
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
