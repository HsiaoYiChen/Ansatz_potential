program wan_interpolat
        use mpi

        implicit none
        complex,allocatable  :: H_wan(:,:,:),H_k(:,:),H_buff(:,:),eig_vec(:,:,:)
        complex,allocatable  :: eig_vec_local(:,:,:),AAR(:,:,:,:)
        complex,allocatable  ::unmk(:,:),mat_tmp(:,:),unmk_Hermit(:,:)
        complex,allocatable  :: delHH(:,:,:),AAq_tmp(:,:,:),AA(:,:,:),D_h(:,:,:)
        integer :: ik,jk,nk,nkx,nky,nkz,irept,nrept,io,read3_int(3)
        integer :: i,j,wan_x,wan_y,wan_z,nstep,ix,iy,iz,iwan,jwan
        integer :: jx,jy,jz,ikxy,H_dim,ibnd,jbnd,i_tmp,j_tmp
        integer :: nspin, period_x, period_y, period_z,super_size
        integer :: c_ix,c_iy,c_iz,c_jx,c_jy,c_jz,n1,n2,i_dim,j_dim,is,js
        integer :: super_i,super_j,c_i,c_j,npt
        integer :: nE,iE,i1,i2,i3,i4,f1,f2,f3,f4,i_dir
        integer,allocatable ::super_deg(:)
        real,allocatable     :: E_list(:),eig_list(:,:),spin_list(:,:),spin_tmp(:),proj_tmp(:)
        real,allocatable     :: proj_list(:,:),x_tmp(:),y_tmp(:),z_tmp(:)
        real,allocatable     :: x_list(:,:),y_list(:,:),z_list(:,:)
        complex,allocatable  :: sig_H(:,:,:),sig_AH(:,:,:)
        complex,allocatable  :: sig_H_local(:,:,:),sig_AH_local(:,:,:)
        real,allocatable     ::k_list(:,:),eigen(:,:),eigen_local(:,:)
        real,allocatable     :: cr_vec(:,:),eig_tmp(:),Sig_kxy(:)
        real,allocatable     :: nb(:,:),R_list(:,:),pt_k(:)
        real,allocatable     :: fold_tmp(:,:),fold(:,:,:)
        real :: a0,a(3,3),Ef,Sig_k,dE,fac,pt_ini_k(3),pt_fin_k(3)
        real :: length,vec_tmp(3),b(3,3)
        integer :: n_cell,nbnd,nwan,i_cell
        integer,allocatable     :: ir_vec(:,:)
        character(len=20)       :: irept_char
        complex      :: cI,sig_xy_local,phase,phase2,ctmp1,ctmp2,ctmp3
        complex      :: I2pi,s(2)
        complex,allocatable:: phase0(:,:),phase1(:,:),phase_k(:)
        complex,allocatable:: H_wan_tmp(:,:),RH_tmp(:,:,:),AA_tmp(:,:,:)
        real         :: pi,E_max,E_min       
        logical, allocatable    :: rept_TF(:)
        integer,allocatable     :: deg(:),inv_k(:)
        character(len=200)      :: useless_string,nscf_dir,THE_dir,path,title
        character(len=20)       :: kz_char,kxy_char,i_char,j_char
        character(len=20)       :: super_size_char
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool,v_max,nline,iline
        real    :: t1,t2,t3,t4,eta,w,celldm3
        complex      ::A2_tmp(2,2)
        complex,allocatable  :: wTw(:,:,:),wSw(:,:,:)
        complex,allocatable  :: uTu(:,:),uSu(:,:),V_tmp(:)
        complex,allocatable  :: uTu_list(:,:,:),uSu_list(:,:,:)
        complex,allocatable  :: rot_U(:,:),rot_Uk(:,:)



open(unit=101,file='../THE.in',action='read',access='sequential',status='old')
read(101,*)!useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,*)!'(A7,I3)')  useless_string,n_func
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I2)') useless_string, nspin
read(101,'(A3,I3)') useless_string, nkx
read(101,'(A3,I3)') useless_string, nky
read(101,'(A3,I3)') useless_string, nkz
read(101,*) ! A6 axial
read(101,'(A9,I3)') useless_string,period_x
read(101,'(A9,I3)') useless_string, period_y
read(101,'(A9,I3)') useless_string, period_z
read(101,*)! useless_string, n_electron
read(101,'(A4,F8.4)') useless_string, eta
read(101,'(A4,F16.8)') useless_string, fac
read(101,*)!'(A6,I3)') useless_string,b_min
read(101,*)!'(A6,I3)') useless_string,b_max
read(101,*)!'(A6,I3)') useless_string,wan_min
read(101,*)!'(A6,I3)') useless_string,wan_max
read(101,'(A3,F16.8)') useless_string, a0
read(101,'(A2,F16.8)') useless_string, celldm3
read(101,'(A3,F16.8)') useless_string, Ef
close(101)

               a(1,:) = (/   1.000000,   0.000000 ,  0.000000 /)
               a(2,:) = (/   0.000000,   0.504533 ,  0.000000 /)
               a(3,:) = (/   0.000000,   0.000000 ,  2.169706 /)

               b(1,:) = (/  1.000000,  0.000000,  0.000000 /)
               b(2,:) = (/  0.000000,  1.982030,  0.000000 /)
               b(3,:) = (/  0.000000,  0.000000,  0.460892 /)

a(:,:)=a(:,:)*a0*0.529177 !Ang

a(1,:)=a(1,:)!*period_x
a(2,:)=a(2,:)!*period_y
a(3,:)=a(3,:)!*period_z


!b(1,:) = b(1,:)/period_x
!b(2,:) = b(2,:)/period_y
!b(3,:) = b(3,:)/period_z


E_max=3.0
E_min=0.0
nE=300

allocate(wTw(nwan,nwan,3))
allocate(wSw(nwan,nwan,3))
allocate(uTu(nwan,3))
allocate(uSu(nwan,3))
allocate(V_tmp(nwan))
allocate(E_list(nE))
allocate(sig_H(2,2,nE))
allocate(sig_AH(2,2,nE))
allocate(sig_H_local(2,2,nE))
allocate(sig_AH_local(2,2,nE))
Do iE=1,nE
        E_list(iE)=E_min+(iE-1)*(E_max-E_min)/nE
ENDDO




cI=cmplx(0.000,1.000)
pi=3.141592653500
I2pi=2*pi*cI
H_dim=nwan*2/nspin*period_x*period_y*period_z




if (.not. mod(nkx,period_x) .eq. 0) then
print*, 'not conmenserate x',nkx,period_x
stop
endif
if (.not. mod(nky,period_y) .eq. 0) then
print*, 'not conmenserate y', nky,period_y
stop
endif
if (.not. mod(nkz,period_z) .eq. 0) then
print*, 'not conmenseratez', nkz,period_z
stop
endif
nkx=nkx/period_x
nky=nky/period_y
nkz=nkz/period_z

allocate(H_k(H_dim,H_dim))
allocate(H_buff(H_dim,H_dim))
!allocate(eig_vec_local(nkx*nky,H_dim,H_dim))
!allocate(eigen_local(nkx*nky,H_dim))
allocate(unmk(H_dim,H_dim))
allocate(eig_tmp(H_dim))
allocate(delHH(H_dim,H_dim,3))
allocate(rot_U(H_dim,H_dim),rot_Uk(H_dim,H_dim))
!allocate(AA(H_dim,H_dim,3))
!allocate(AAq_tmp(H_dim,H_dim,3))
allocate(D_h(H_dim,H_dim,3))
allocate(mat_tmp(H_dim,H_dim))
allocate(proj_tmp(H_dim))
allocate(x_tmp(H_dim))
allocate(y_tmp(H_dim))
allocate(z_tmp(H_dim))
allocate(unmk_Hermit(H_dim,H_dim))
allocate(spin_tmp(H_dim))
irept=0
open(unit=101,file='./super_H/ir_vec_HYC',action='read',access='sequential')
do
        read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                irept=irept+1
        enddo
close(101)
nrept=irept
allocate(ir_vec(nrept,3))
allocate(cr_vec(nrept,3))
allocate(H_wan(H_dim,H_dim,nrept))
allocate(H_wan_tmp(H_dim,H_dim))
allocate(RH_tmp(H_dim,H_dim,3))
allocate(super_deg(nrept))
!allocate(AAR(H_dim,H_dim,nrept,3))
!allocate(AA_tmp(H_dim,H_dim,3))
        !allocate(deg(nrept))
irept=1
open(unit=101,file='./super_H/ir_vec_HYC',action='read',access='sequential')
do
        read(101,'(3I8)',iostat=io) read3_int(:)
        if(io/=0) exit
        ir_vec(irept,:)=read3_int(:)
        irept=irept+1
enddo
close(101)


ir_vec(:,1)=ir_vec(:,1)*period_x
ir_vec(:,2)=ir_vec(:,2)*period_y
ir_vec(:,3)=ir_vec(:,3)*period_z



open(unit=101,file='./super_H/H_deg_HYC',action='read',access='sequential')
do
        read(101,'(2I8)',iostat=io) i,j
        super_deg(i)=j
        if (io/=0)  exit
enddo
close(101)





        !open(unit=102,file='./H_sk/data/rot_U',action='read',access='sequential'&
       ! open(unit=102,file='./rot_U',action='read',access='sequential'&        
       ! &,status='old',form='unformatted')
       ! read(102)      rot_U(:,:)
       ! close(102)




call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,size,ierr)


open(unit=101,file='./plt_path.in',action='read',access='sequential',status='old')
read(101,'(A200)') title
read(101,'(A6,I3)') useless_string,nline  !nline=
read(101,'(A4,I3)') useless_string,npt  !npt=
nk=nline*npt
allocate(pt_k(nk))
allocate(eig_list(nk,H_dim))
allocate(spin_list(nk,H_dim))
allocate(proj_list(nk,H_dim))
allocate(x_list(nk,H_dim))
allocate(y_list(nk,H_dim))
allocate(z_list(nk,H_dim))
allocate(k_list(nk,3))
npool=nk
allocate(pool(npool,1))
if (my_rank .eq. 0 ) then
open(unit=102,file='./path_check.out',action='write',access='append',status='replace')
endif
do iline=1,nline
        read(101,'(6F12.6)') pt_ini_k(:),pt_fin_k(:)
        if (my_rank .eq. 0 ) write(102,'(6F12.6)') pt_ini_k(:),pt_fin_k(:)
        vec_tmp(1)=sum(b(:,1)*(pt_fin_k(:)-pt_ini_k(:)))
        vec_tmp(2)=sum(b(:,2)*(pt_fin_k(:)-pt_ini_k(:)))/6
        vec_tmp(3)=sum(b(:,3)*(pt_fin_k(:)-pt_ini_k(:)))
        length=abs(sum(vec_tmp(:)**2))**0.5/npt
        do ik=npt*(iline-1)+1,npt*iline
                k_list(ik,:)=pt_ini_k(:)+(ik-npt*(iline-1))*(pt_fin_k(:)-pt_ini_k(:))/npt
                pool(ik,1)=ik
                if (ik .eq. 1) then
                pt_k(ik)=length
                else
                pt_k(ik)=pt_k(ik-1)+length
                endif
        enddo
enddo
close(101)
if (my_rank .eq. 0 ) close(102)


k_list(:,1)=k_list(:,1)/period_x
k_list(:,2)=k_list(:,2)/period_y
k_list(:,3)=k_list(:,3)/period_z



super_size=period_x*period_y*period_z
call int2char(super_size,super_size_char)
allocate(nb(super_size,3))
allocate(phase0(super_size,super_size))
allocate(phase1(super_size,super_size))
allocate(phase_k(nrept))
allocate(R_list(super_size,3))
allocate(fold_tmp(H_dim,super_size))
allocate(fold(nk,H_dim,super_size))
do ix=0,period_x-1
do iy=0,period_y-1
do iz=0,period_z-1
        i=1+iz+period_z*iy+period_z*period_y*ix
        nb(i,:)=(/ix*1.000/period_x,iy*1.000/period_y,iz*1.000/period_z/)
        
        R_list(i,1)=sum(a(:,1)*(/ix,iy,iz/))
        R_list(i,2)=sum(a(:,2)*(/ix,iy,iz/))
        R_list(i,3)=sum(a(:,3)*(/ix,iy,iz/))


enddo
enddo
enddo

muti_pool=0
DO
        muti_pool=muti_pool+1
        IF(muti_pool*(size-1) .gt. npool) EXIT
ENDDO
muti_pool=muti_pool-1


if (my_rank.eq.0)  then

        open(unit=101,file='report',action='write',access='append',status='replace')
        write(101,*) 'start reading'
        close(101)
        
        do irept=1,nrept
        !print*,irept
                call int2char(irept, irept_char)
                open(unit=101,file='./tot_HR/hr_'//trim(irept_char),action='read'&
                        &,access='stream',form='unformatted')
                        read(101) H_wan(:,:,irept)
                close(101)
                call MPI_Bcast( H_wan(:,:,irept),H_dim**2 , MPI_complex,0,MPI_COMM_WORLD, ierr)
                call MPI_Barrier(  MPI_COMM_WORLD, ierr)
                open(unit=101,file='report',action='write',access='append',status='old')
                write(101,*) 'read HR irept=,',irept,'out of ',nrept, 'done'
                close(101)



        enddo

        eig_list(:,:) =0.000

        Do iz=1,1
                call int2char(iz-1,kz_char)
                do ipool=1,nk

                        call MPI_RECV(eig_list(ipool,:),H_dim,MPI_REAL,&
                        &MPI_ANY_SOURCE,ipool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        call MPI_RECV(spin_list(ipool,:),H_dim,MPI_complex,&
                        &MPI_ANY_SOURCE,ipool+npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        call MPI_RECV(proj_list(ipool,:),H_dim,MPI_complex,&
                        &MPI_ANY_SOURCE,ipool+2*npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        call MPI_RECV(x_list(ipool,:),H_dim,MPI_complex,&
                        &MPI_ANY_SOURCE,ipool+3*npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        call MPI_RECV(y_list(ipool,:),H_dim,MPI_complex,&
                        &MPI_ANY_SOURCE,ipool+4*npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        call MPI_RECV(z_list(ipool,:),H_dim,MPI_complex,&
                        &MPI_ANY_SOURCE,ipool+5*npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                enddo
                call MPI_Barrier(  MPI_COMM_WORLD, ierr)

        enddo


        open(unit=101,file='./plt_bnd_'//trim(title),action='write',access='append',status='replace')
        do i_dim=1,H_dim
        do ik=1,nk
               ! if (abs(eig_list(ik,i_dim)-Ef) .gt. 0.25) cycle
                write(101,'(I8,10F18.6)') ik, k_list(ik,:),pt_k(ik),eig_list(ik,i_dim)&
                &,proj_list(ik,i_dim),spin_list(ik,i_dim),x_list(ik,i_dim),y_list(ik,i_dim),z_list(ik,i_dim)
        enddo
                write(101,*)
        enddo
        close(101)





endif



if (.not. my_rank.eq.0)  then

        do irept=1,nrept
                call MPI_Bcast( H_wan(:,:,irept),H_dim**2 ,MPI_complex,0,MPI_COMM_WORLD, ierr)
                call MPI_Barrier(  MPI_COMM_WORLD, ierr)
        enddo

        !open(unit=101,file='./data_wTw/wSw',action='read',access='stream',status='old',form='unformatted')
        !read(101) wSw(:,:,:)
        !close(101)
        !open(unit=101,file='./data_wTw/wTw',action='read',access='stream',status='old',form='unformatted')
        !read(101) wTw(:,:,:)
        !close(101)






        Do iz=1,1!HYC_z,HYC_z
        DO imuti_pool=0,muti_pool
                ipool=imuti_pool*(size-1)+my_rank
                if (ipool .gt. npool) cycle
                call cpu_time(t1)
                ikxy=pool(ipool,iz)
                ik=ikxy!1+(iz-1)+(ikxy-1)*nkz


                H_k(:,:)=cmplx(0.000,0.000)
                do irept=1,nrept
                phase_k(irept)=exp(I2pi*sum(ir_vec(irept,:)*k_list(ik,:)))
                H_k(:,:)=H_k(:,:)+H_wan(:,:,irept)*phase_k(irept)/super_deg(irept)
                enddo

               
                !H_k=H_k/super_size**2
                call cpu_time(t2)
                call  diasym(H_k,eig_tmp(:),H_dim,unmk(:,:))
 
        !print*,eig_tmp(:)              
        !                        write(102,'(I8,10E24.10)'),i_dim, abs(rot_Uk(i_dim,1241:1250))**2
        !                enddo
        !                close(102)
        !        endif
                rot_Uk(:,:)=matmul(transpose(conjg(rot_U(:,:))),unmk(:,:))
                !rot_Uk(:,:)=unmk(:,:)



                spin_tmp(:)=0.00
                proj_tmp(:)=0.00
                x_tmp(:)=0.00
                y_tmp(:)=0.00
                z_tmp(:)=0.00
                
                do i_dim=1,H_dim
                do j_dim=1,H_dim
       cycle
                spin_tmp(i_dim)=spin_tmp(i_dim)-(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                enddo
                enddo
                
                
                do i_dim=1,H_dim
                do i=1,1
      !cycle
                !cycle               
                !do j_dim=9,10!H_dim
                 !       if (.not. eig_tmp(j_dim) .ge. 9) cycle
                        j_dim=9+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=10+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=19+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=20+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=29+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=30+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=39+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=40+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        do j=1,8
                        j_dim=87+(j-1)*12+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        j_dim=88+(j-1)*12+(i-1)*128
                        spin_tmp(i_dim)=spin_tmp(i_dim)+(-1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        proj_tmp(i_dim)=proj_tmp(i_dim)+( 1.0)**(j_dim+1)*abs(rot_Uk(j_dim,i_dim))**2
                        enddo
                        
                        
                        
                enddo
                enddo
                
                !spin_tmp(:)=spin_tmp(:)/proj_tmp(:)



                
                call cpu_time(t4)

        if ( my_rank.eq.1)  then
                open(unit=101,file='report',action='write',access='append')
                        write(101,*) 'build', t2-t1, 'diag', t4-t2, 'all',t4-t1
                close(101)
        endif

       ! print*,'build', t2-t1, 'diag', t3-t2, 'dipole', t4-t3, 'all',t4-t1

         !       print*,'build', t2-t1, 'diag', t3-t2, 'dipole', t4-t3, 'all',t4-t1,sig_xy_local



                call MPI_SEND(eig_tmp(:),H_dim&
                                &,MPI_REAL,0,ikxy,MPI_COMM_WORLD,ierr)
                call MPI_SEND(spin_tmp(:),H_dim&
                                &,MPI_REAL,0,ikxy+npool,MPI_COMM_WORLD,ierr)
                call MPI_SEND(proj_tmp(:),H_dim&
                                &,MPI_REAL,0,ikxy+2*npool,MPI_COMM_WORLD,ierr)

                call MPI_SEND(x_tmp(:),H_dim&
                                &,MPI_REAL,0,ikxy+3*npool,MPI_COMM_WORLD,ierr)

                call MPI_SEND(y_tmp(:),H_dim&
                                &,MPI_REAL,0,ikxy+4*npool,MPI_COMM_WORLD,ierr)
                call MPI_SEND(z_tmp(:),H_dim&
                                &,MPI_REAL,0,ikxy+5*npool,MPI_COMM_WORLD,ierr)



        enddo

        call MPI_Barrier(  MPI_COMM_WORLD, ierr)
        enddo


endif



call mpi_finalize(ierr)



contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )    WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )    WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )    WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )   WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) ) WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine diasym(a,eig,n,D)
 implicit none

 integer,intent(in)  :: n
 complex,intent(in)  :: a(n,n)
 complex              :: a_tmp(n,n)
 complex,intent(out)  :: D(n,n)

 integer  :: l,inf
 complex  :: Lwork(2*n-1)
 real,intent(out)      :: eig(n)
 real                 :: Rwork(3*n-2)

 a_tmp(:,:)=a(:,:)

 l=2*n-1
 call cheev('V','U',n,a_tmp,n,eig,Lwork,l,Rwork,inf)
 D(:,:)=a_tmp(:,:)



end subroutine diasym



end program
