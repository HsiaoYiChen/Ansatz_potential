program rearrange_HR

        implicit none
        complex(kind(0d0)),allocatable  ::H_wan(:,:,:),H_k(:,:),H_buff(:,:),eig_vec(:,:,:)
        complex(kind(0d0)),allocatable  ::H_k_super(:,:),H_k_super_uc(:,:),Us(:,:)
        complex(kind(0d0)),allocatable  ::H_k_super_tmp(:,:)
        complex(kind(0d0)),allocatable  ::AAR(:,:,:,:),H_soc(:,:,:),super_H_soc(:,:,:)
        complex(kind(0d0)),allocatable  :: H_FM(:,:,:),super_H_FM(:,:,:)
        complex(kind(0d0)),allocatable  :: super_H_wan(:,:,:),super_AAR(:,:,:,:),Dh(:,:,:)
        complex(kind(0d0)),allocatable  :: H_check(:,:),u_check(:,:),delHH(:,:,:),M_tmp(:,:)
        real(kind(0d0)),allocatable::eig_check(:)
        integer :: ik,jk,nk,nkx,nky,nkz,irept,nrept,io,read3_int(3)
        integer :: super_nrept,super_nx,super_ny,super_nz,H_dim
        integer :: i,j,wan_x,wan_y,wan_z,nstep,ix,iy,iz,iwan,jwan
        integer :: ib,ibx,iby,ibz,ikp,iks,ikx,iky,ikz
        integer :: nscf_nkx,nscf_nky,nscf_nkz
        integer :: i_dim,j_dim,nspin,is,js,ir
        integer :: jx,jy,jz,period_x,period_y,period_z,i_shift,n1,n2,n3
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j
        integer :: super_x,super_y,super_z,super_rept,super_size
        real(kind(0d0)),allocatable     ::k_list(:,:),eigen(:,:)
        real(kind(0d0)),allocatable     ::super_k_list(:,:)
        integer :: n_cell,nbnd,nwan,i_cell,soc_i,soc_j
        integer,allocatable     :: ir_vec(:,:),super_ir_vec(:,:)
        real(kind(0d0)),allocatable     :: nb(:,:)
        character(len=20)       :: irept_char,ks_char,kp_char
        complex(kind(0d0))      :: cI,phase1,phase2
        real(kind(0d0))         :: pi,dE,read3(3),k_tmp(3),m(3)
        logical, allocatable    :: rept_TF(:)
        integer,allocatable     :: deg(:),inv_k(:)
        character(len=200)      :: useless_string,nscf_dir,THE_dir,path
        logical ::SOC_TF




cI=dcmplx(0.0d0,1.0d0)
pi=3.1415926535d0



open(unit=101,file='../THE.in',action='read',access='sequential',status='old')
read(101,*)!useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,*)!'(A7,I3)')  useless_string,n_func
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I3)') useless_string, nspin
read(101,*)!'(A3,I3)') useless_string, nkx
read(101,*)!'(A3,I3)') useless_string, nky
read(101,*)!'(A3,I3)') useless_string, nkz
read(101,*) ! A6 axial
read(101,'(A9,I3)')  useless_string, period_x
read(101,'(A9,I3)')  useless_string, period_y
read(101,'(A9,I3)')  useless_string, period_z
!read(101, useless_string, n_electron
!read(101,'(A4,F8.4)') useless_string, eta
!read(101,'(A4,F16.8)') useless_string, fac
!read(101,'(A6,I3)') useless_string,b_min
!read(101,'(A6,I3)') useless_string,b_max
close(101)



H_dim=nwan*2/nspin
allocate(H_k(H_dim,H_dim))
allocate(H_k_super(H_dim*period_x*period_y*period_z,&
                  &H_dim*period_x*period_y*period_z))

allocate(H_k_super_uc(H_dim*period_x*period_y*period_z,&
                  &H_dim*period_x*period_y*period_z))   
allocate(H_k_super_tmp(H_dim*period_x*period_y*period_z,&
                  &H_dim*period_x*period_y*period_z))          
allocate(Us(H_dim*period_x*period_y*period_z,&
                  &H_dim*period_x*period_y*period_z))

allocate(H_buff(H_dim*period_x*period_y*period_z,&
                  &H_dim*period_x*period_y*period_z))
allocate(eigen(H_dim*period_x*period_y*period_z,1))


k_tmp(:)=9999.0

open(unit=101,file=trim(nscf_dir)//'./k_list',status='old',action='read')
do 
        read(101,'(3F12.8)',iostat=io) read3(:)
        if (io/=0) exit

        if (abs(read3(1)) .lt. k_tmp(1) .and. abs(read3(1)) .gt. 0.000001) k_tmp(1)=abs(read3(1))
        if (abs(read3(2)) .lt. k_tmp(2) .and. abs(read3(2)) .gt. 0.000001) k_tmp(2)=abs(read3(2))
        if (abs(read3(3)) .lt. k_tmp(3) .and. abs(read3(3)) .gt. 0.000001) k_tmp(3)=abs(read3(3))
enddo
close(101)

 nscf_nkx=nint(1.0/k_tmp(1))
 nscf_nky=nint(1.0/k_tmp(2))
 nscf_nkz=nint(1.0/k_tmp(3))
 print*,'nk1,nk2,nk3=', nscf_nkx,nscf_nky,nscf_nkz
 print*,'periodic=',period_x,period_y,period_z
 nk=nscf_nkx*nscf_nky*nscf_nkz
 allocate(k_list(nk,3))

 if (.not. mod(nscf_nkx,period_x) .eq. 0) then
         print*,' x  not commense'
         stop
 endif
 if (.not. mod(nscf_nky,period_y) .eq. 0) then
         print*,' y  not commense'
         stop
 endif
 if (.not. mod(nscf_nkz,period_z) .eq. 0) then
         print*,' z  not commense'
         stop
 endif

ik=1
open(unit=101,file=trim(nscf_dir)//'./k_list',status='old',action='read')
do
        read(101,'(3F12.8)',iostat=io) read3(:)
        if (io/=0) exit
        k_list(ik,:)=read3(:)
        ik=ik+1
enddo
close(101)


super_nx= nscf_nkx/period_x
super_ny= nscf_nky/period_y
super_nz= nscf_nkz/period_z

super_nrept=super_nx*super_ny*super_nz


allocate(super_k_list(super_nrept,3))
do ix=0,super_nx-1
do iy=0,super_ny-1
do iz=0,super_nz-1
        ik=1+iz+super_nz*iy+super_ny*super_nz*ix
        super_k_list(ik,:)=(/ix*1.0/super_nx,iy*1.0/super_ny,iz*1.0/super_nz/)
enddo
enddo
enddo

open(unit=101,file='.//Hk_commense/k_list',action='write',access='append',status='replace')
do ik=1,super_nrept
        write(101,'(3F12.8)') super_k_list(ik,:)
enddo
close(101)

do ikx=0,super_nx-1
do iky=0,super_ny-1
do ikz=0,super_nz-1
        iks=1+ikz+super_nz*iky+super_nz*super_ny*ikx
        print*,iks, super_nrept
        m(1)= super_k_list(iks,1)/period_x
        m(2)= super_k_list(iks,2)/period_y
        m(3)= super_k_list(iks,3)/period_z
        call int2char(iks,ks_char)
        Us(:,:)=dcmplx(0.0d0,0.0d0)
        H_k_super_uc(:,:)=dcmplx(0.0d0,0.0d0)
        do ibx=0,period_x-1
        do iby=0,period_y-1
        do ibz=0,period_z-1
                print*,'doing',ibx,iby,ibz
                jx=ikx+super_nx*ibx
                jy=iky+super_ny*iby
                jz=ikz+super_nz*ibz
                ikp=1+jz+nscf_nkz*jy+nscf_nkz*nscf_nky*jx
                call int2char(ikp,kp_char)
                        open(unit=101,file=trim(nscf_dir)//'./make_full/basis_data/Hkw_k'//trim(kp_char),&
                                &action='read',access='sequential',status='old',form='unformatted')
                                read(101) H_k
                        close(101)

        
          !      call  diasym(H_k(:,:),eigen(1:H_dim,1),H_dim,H_buff(1:H_dim,1:H_dim))        
        
          !      print*,eigen(1:H_dim,1)
        !stop

                ib=1+ibz+period_z*iby+period_z*period_y*ibx
                H_k_super_uc(1+(ib-1)*H_dim:ib*H_dim,1+(ib-1)*H_dim:ib*H_dim)=H_k(:,:)
        


                do iwan=1,H_dim
                do ix=0,period_x-1
                do iy=0,period_y-1
                do iz=0,period_z-1

                        ir=1+iz+period_z*iy+period_z*period_y*ix
                        
                        Us(iwan+(ib-1)*H_dim,iwan+(ir-1)*H_dim)=exp(-2*pi*cI*sum(m(:)*(/ix,iy,iz/)))*&
                              &exp(-2*pi*cI*sum((/ibx*1.0/period_x,iby*1.0/period_y,ibz*1.0/period_z/)*(/ix,iy,iz/)))   !(b,r)
                


                enddo
                enddo
                enddo
                enddo

                
        enddo
        enddo
        enddo

        H_k_super_tmp(:,:)=matmul( H_k_super_uc(:,:), Us(:,:))
        H_k_super(:,:)=matmul(dconjg(transpose(Us(:,:))),H_k_super_tmp(:,:))

        H_k_super(:,:)=H_k_super(:,:)/(period_x*period_y*period_z)

        open(unit=101,file='./Hk_commense/Hks_'//trim(ks_char),action='write',access='append',&
                &status='replace',form='unformatted')
        write(101) H_k_super(:,:)
        close(101)




      !  call  diasym(H_k_super(:,:),eigen(:,1),H_dim*period_x*period_y*period_z,H_buff(:,:))


      !  print*,eigen(:,1)
      !  stop


enddo
enddo
enddo










contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )    WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine diasym(a,eig,n,D)
 implicit none

 integer,intent(in)  :: n
 complex(kind(0d0)),intent(in)  :: a(n,n)
 complex(kind(0d0))              :: a_tmp(n,n)
 complex(kind(0d0)),intent(out)  :: D(n,n)

 integer  :: l,inf
 complex(kind(0d0))  :: Lwork(2*n-1)
 real(kind(0d0)),intent(out)      :: eig(n)
 real(kind(0d0))                 :: Rwork(3*n-2)

 a_tmp(:,:)=a(:,:)

 l=2*n-1
 call zheev('V','U',n,a_tmp,n,eig,Lwork,l,Rwork,inf)
 D(:,:)=a_tmp(:,:)



end subroutine diasym


end program
