program cal_qww
        use mpi
        implicit none
        complex(kind(0d0)),allocatable  ::H_wan(:,:,:),H_k(:,:),H_buff(:,:),eig_vec(:,:,:)
        complex(kind(0d0)),allocatable  ::H_k_super(:,:),H_k_super_uc(:,:),Us(:,:)
        complex(kind(0d0)),allocatable  ::H_k_super_tmp(:,:)
        complex(kind(0d0)),allocatable  :: H_soc(:,:,:),super_H_soc(:,:,:)
        complex(kind(0d0)),allocatable  :: H_FM(:,:,:),super_H_FM(:,:,:)
        complex(kind(0d0)),allocatable  :: super_H_wan(:,:,:)
        complex(kind(0d0)),allocatable  :: H_check(:,:),u_check(:,:),M_tmp(:,:)
        complex(kind(0d0)),allocatable  :: wfn1(:),wfn2(:)
        real(kind(0d0)),allocatable::eig_check(:)
        integer :: ik,jk,nk,nkx,nky,nkz,irept,nrept,io,read3_int(3)
        integer :: super_nrept,super_nx,super_ny,super_nz,H_dim
        integer :: i,j,wan_x,wan_y,wan_z,nstep,ix,iy,iz,iwan,jwan
        integer :: ib,ibx,iby,ibz,ikp,iks,ikx,iky,ikz,ikqs
        integer :: jkx,jky,jkz,iGx,iGy,iGz,iqx,iqy,iqz
        integer :: jbx,jby,jbz,iksx,iksy,iksz,ikq,i_super
        integer :: nscf_nkx,nscf_nky,nscf_nkz,nx,ny,nz
        integer :: i_dim,j_dim,nspin,is,js,ir,n_super,i_mag
        integer :: jx,jy,jz,period_x,period_y,period_z,i_shift,n1,n2,n3
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j
        integer :: super_x,super_y,super_z,super_rept,super_size
        real(kind(0d0)),allocatable     ::k_list(:,:),eigen(:,:)
        real(kind(0d0)),allocatable     ::super_k_list(:,:)
        integer,allocatable     :: super_k_int(:,:),k_ind(:),k_int(:,:)
        integer :: n_cell,nbnd,nwan,i_cell,soc_i,soc_j
        integer,allocatable     :: ir_vec(:,:),super_ir_vec(:,:)
        real(kind(0d0)),allocatable     :: nb(:,:)
        character(len=20)       :: irept_char,ks_char,kp_char
        complex(kind(0d0))      :: cI,phase1,phase2,ctmp
        complex(kind(0d0)),allocatable :: qww(:,:,:,:),phase(:),qww_local(:,:,:)
        complex           ,allocatable :: qww_single(:,:,:,:)
        real(kind(0d0))         :: pi,dE,read3(3),k_tmp(3),m(3),q_vec(3)
        logical, allocatable    :: rept_TF(:)
        integer,allocatable     :: deg(:),inv_k(:),nwan_per_mag_unit(:),cal_unit_wan(:,:)
        character(len=200)      :: useless_string,nscf_dir,THE_dir,path
        character(len=20)      :: k_char,kq_char,iwan_char,jwan_char,q_char
        logical ::SOC_TF
        integer :: cal_q,n_mag_unit
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool
       



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





open(unit=102,file='./mag.in',access='sequential',action='read',status='old')
  read(102,*)!'(A3,F8.4)') useless_string, Ef
  read(102,*)!'(A6,F8.4)') useless_string, delta
  read(102,*)!'(A5,F12.6)') useless_string, Emin
  read(102,*)!'(A5,F12.6)') useless_string, Emax
  read(102,*)!'(A3,I4)') useless_string, nE
  read(102,'(A6,I5)') useless_string,cal_q
  read(102,*)!------------------------------------------
  read(102,'(A11,I3)') useless_string, n_mag_unit
  allocate(nwan_per_mag_unit(n_mag_unit))
  nwan_per_mag_unit(:)=0
  do i_mag=1,N_mag_unit
        read(102,*) nwan_per_mag_unit(i_mag)
  enddo
  allocate(cal_unit_wan(N_mag_unit,maxval(nwan_per_mag_unit(:))))
  read(102,*)!------------------------------------------
  do i_mag=1,N_mag_unit
  do iwan=1,nwan_per_mag_unit(i_mag)
        read(102,*)  cal_unit_wan(i_mag,iwan)
  enddo
  enddo
close(102)


call int2char(cal_q,q_char)


if (.not. nspin .eq. 1) then
        print*, 'only for nspin = 1'
        stop
endif

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
 allocate(k_int(nk,3))
 allocate (k_ind(nk))
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
        k_int(ik,1)=nint(read3(1)*nscf_nkx)
        k_int(ik,2)=nint(read3(2)*nscf_nky)
        k_int(ik,3)=nint(read3(3)*nscf_nkz)
        k_ind(1+k_int(ik,3)+nscf_nkz*k_int(ik,2)+nscf_nkz*nscf_nky*k_int(ik,1))=ik
        ik=ik+1
enddo
close(101)





super_nx= nscf_nkx/period_x
super_ny= nscf_nky/period_y
super_nz= nscf_nkz/period_z
n_super=period_x*period_y*period_z

super_nrept=super_nx*super_ny*super_nz


allocate(super_k_list(super_nrept,3))
allocate(super_k_int (super_nrept,3))
do ix=0,super_nx-1
do iy=0,super_ny-1
do iz=0,super_nz-1
        ik=1+iz+super_nz*iy+super_ny*super_nz*ix
        super_k_list(ik,:)=(/ix*1.0/super_nx,iy*1.0/super_ny,iz*1.0/super_nz/)
        super_k_int(ik,:)=(/ix,iy,iz/)
enddo
enddo
enddo




open(unit=101,file='.//Hk_commense/k_list',action='write',access='append',status='replace')
do ik=1,super_nrept
        write(101,'(3F12.8)') super_k_list(ik,:)
enddo
close(101)




open(unit=101,file=trim(nscf_dir)//'./V_HYC_dp',access='sequential',form='unformatted')
read(101) nx,ny,nz
print*,nx,ny,nz
close(101)


allocate(wfn1(nx*ny*nz))
allocate(wfn2(nx*ny*nz))
allocate(phase(nx*ny*nz))

iqx=super_k_int(cal_q,1)
iqy=super_k_int(cal_q,2)
iqz=super_k_int(cal_q,3)
        
q_vec(1)=iqx*1.0/nscf_nkx
q_vec(2)=iqy*1.0/nscf_nky
q_vec(3)=iqz*1.0/nscf_nkz


npool=super_nx*super_ny*super_nz
allocate(pool(npool,3))

do iksx=0,super_nx-1
do iksy=0,super_ny-1
do iksz=0,super_nz-1
        iks=1+iksz+super_nz*iksy+super_nz*super_ny*iksx
        pool(iks,:)=(/iksx,iksy,iksz/)
enddo
enddo
enddo


  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,size,ierr)

  print*,'mpi ini'
  muti_pool=0
  DO
        muti_pool=muti_pool+1
        IF(muti_pool*(size-1) .gt. npool) EXIT
  ENDDO
  muti_pool=muti_pool-1


if (my_rank.eq.0)  then
        open(unit=101,file='report',action='write',access='append',status='replace')
                write(101,*) 'start'
                close(101)

        call EXECUTE_COMMAND_LINE("mkdir -p  dir-qww")

        allocate(qww(super_nrept,n_super*n_mag_unit,maxval(nwan_per_mag_unit(:)),maxval(nwan_per_mag_unit(:))))
        allocate(qww_single(super_nrept,n_super*n_mag_unit,maxval(nwan_per_mag_unit(:)),maxval(nwan_per_mag_unit(:))))
        !allocate(qww_local(n_super*n_mag_unit,maxval(nwan_per_mag_unit(:)),maxval(nwan_per_mag_unit(:))))

        do ipool=1,npool
        call MPI_RECV(qww(ipool,:,:,:),n_super*n_mag_unit*maxval(nwan_per_mag_unit(:))**2,mpi_double_complex,&
                        &MPI_ANY_SOURCE,ipool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        enddo

        

        qww_single(:,:,:,:)=qww(:,:,:,:)/n_super

        print*,sum(abs(qww_single(:,:,:,:)))

        open(unit=103,file='dir-qww/qww_q'//trim(q_char)//'.dat',action='write'&
                &, access='append', status='replace',form='unformatted')
        write(103) qww_single
        close(103)
endif
if ( .not. my_rank.eq.0)  then
        allocate(qww_local(n_super*n_mag_unit,maxval(nwan_per_mag_unit(:)),maxval(nwan_per_mag_unit(:))))

        DO imuti_pool=0,muti_pool
                ipool=imuti_pool*(size-1)+my_rank
                if (ipool .gt. npool) cycle
                iks=ipool
                iksx=pool(ipool,1)
                iksy=pool(ipool,2)
                iksz=pool(ipool,3)

                qww_local(:,:,:)=dcmplx(0.0d0,0.0d0)



        print*,iks
        do ibx=0,period_x-1
        do iby=0,period_y-1
        do ibz=0,period_z-1
                ikx=iksx+ibx*super_nx
                iky=iksy+iby*super_ny
                ikz=iksz+ibz*super_nz
                ik=k_ind(1+ikz+nscf_nkz*iky+nscf_nkz*nscf_nky*ikx)

         !       if (ik .eq. 0) then
         !               print*,iksx,iksy,iksz
         !               print*,ikx,iky,ikz
         !               print*,period_x,period_y,period_z
         !               print*,ibx,iby,ibz
!
!                print*,1+ikz+nscf_nkz*iky+nscf_nkz*nscf_nky*ikx,ik
!                        stop
!                endif
                jkx=ikx+iqx
                jky=iky+iqy
                jkz=ikz+iqz

                iGx=0
                iGy=0
                iGz=0
                if (jkx .lt.  0) then
                        iGx=-1
                        jkx=jkx+nscf_nkx
                elseif(jkx .ge. nscf_nkx ) then
                        iGx=1
                        jkx=jkx-nscf_nkx
                endif

                if (jky .lt.  0) then
                        iGy=-1
                        jky=jky+nscf_nky
                elseif(jky .ge. nscf_nky ) then
                        iGy=1
                        jky=jky-nscf_nky
                endif

                if (jkz .lt.  0) then
                        iGz=-1
                        jkz=jkz+nscf_nkz
                elseif(jkz .ge. nscf_nkz ) then
                        iGz=1
                        jkz=jkz-nscf_nkz
                endif


                ikq=k_ind(1+jkz+nscf_nkz*jky+nscf_nkz*nscf_nky*jkx)
 
 !               print*,1+jkz+nscf_nkz*jky+nscf_nkz*nscf_nky*jkx,ikq
!

                do ix=0,nx-1
                do iy=0,ny-1
                do iz=0,nz-1
                        phase(1+ix+nx*iy+nx*ny*iz)=exp(-2*3.1415926535d0*cI*&
                        &sum((/iGx,iGy,iGz/)*(/ ix*1.0/nx,iy*1.0/ny,iz*1.0/nz /)))
                enddo
                enddo
                enddo



                do i_mag=1,n_mag_unit
                do iwan=1,nwan_per_mag_unit(i_mag)
                        call int2char(ikq,kq_char)               
                        call int2char(cal_unit_wan(i_mag,iwan),iwan_char)
                        path=trim(nscf_dir)//'/disentangled/WFN/wfn_k'//trim(kq_char)//'_w'//trim(iwan_char)//'_vec_dp'
                        open(unit=101,file=trim(path)&
                                &,action='read',access='sequential',form='unformatted')
                                read(101) wfn1
                        close(101)

                        wfn1(:)=wfn1(:)/sum(abs(wfn1(:))**2)**0.5*phase(:)

                do jwan=1,nwan_per_mag_unit(i_mag)
                        call int2char(ik,k_char)
                        call int2char(cal_unit_wan(i_mag,jwan),jwan_char)
                        path=trim(nscf_dir)//'/disentangled/WFN/wfn_k'//trim(k_char)//'_w'//trim(jwan_char)//'_vec_dp'
                        open(unit=101,file=trim(path)&
                                &,action='read',access='sequential',form='unformatted')
                                read(101) wfn2
                        close(101)

                        wfn2(:)=wfn2(:)/sum(abs(wfn2(:))**2)**0.5

                        
                        ctmp=sum(wfn1(:)*dconjg(wfn2(:)))


                        do ix=0,period_x-1
                        do iy=0,period_y-1
                        do iz=0,period_z-1

                                i_super=1+iz+period_z*iy+period_z*period_y*ix
                                qww_local(i_mag+n_mag_unit*(i_super-1),iwan,jwan)=&
                                        &qww_local(i_mag+n_mag_unit*(i_super-1),iwan,jwan)+&
                                        &exp(-2*3.1415926535d0*cI*sum((/ix,iy,iz/)*q_vec(:)))*ctmp
                        enddo ! iz
                        enddo ! iy
                        enddo ! ix
                enddo ! jwan
                enddo ! iwan
                enddo ! i_mag
        enddo ! ibz
        enddo ! iby
        enddo ! ibx

        call MPI_SEND(qww_local(:,:,:),n_super*n_mag_unit*maxval(nwan_per_mag_unit(:))**2&
                                        &,MPI_double_complex,0,ipool,MPI_COMM_WORLD,ierr)


        enddo ! ipool

endif
call mpi_finalize(ierr)




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
