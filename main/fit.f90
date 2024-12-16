program fitting_B

        !This code fits the ansatz coefficients from the extracted potential 
        !Detail should be changed for different material and ansatz 
        !We consider a linear ansatz between V and the four Spins
        !By changing  act='plot', it can output a predicted potential 
        !with a given spin orientation to be compared using plt_compare.f90 

        implicit none

        integer :: i,ix,iy,iz,nx,ny,nz,j
        integer :: z_box
        integer :: n_file, n_kind,nR,i_file,i_kind
        integer :: i_dir,j_dir,i2_dir
       
        real(kind(0d0)),allocatable  :: read_VA(:,:)
        real(kind(0d0)),allocatable  :: VAs(:,:,:,:)
        real(kind(0d0)),allocatable  :: spin(:,:,:)
        real(kind(0d0)),allocatable  :: r_vec(:,:,:,:)
        real(kind(0d0)),allocatable  :: read_angle(:)

        character(len=20)       :: file_char,kind_char
        character(len=200)      :: path,read_angle_form
        real(kind(0d0)),allocatable  :: VA_ave(:,:),VA_abs_sum(:,:)
        real(kind(0d0)),allocatable  :: V0(:,:),Vn(:,:)
        real(kind(0d0)),allocatable  :: B0(:),BH(:,:),JH(:,:)
        real(kind(0d0)),allocatable  :: B_DM(:,:,:)
        real(kind(0d0)),allocatable  :: B0_aniso(:,:),BH_aniso(:,:,:)
        real(kind(0d0)),allocatable  :: B_general(:,:,:,:)
        real(kind(0d0)) ::pi
        real(kind(0d0)) ::lambda,lambda_H,lambda_DM,lambda_B
        real(kind(0d0)),allocatable  :: DV(:,:),cal_V(:,:)
        integer :: i_iter,n_iter,i_count
        integer :: n_spinor
        real(kind(0d0)),allocatable:: Err(:),r_hat(:,:),sub_skip(:)
        real:: skip,l_rate,wt
        character(len=4) :: act
        integer,allocatable:: one_zero(:)
        integer,allocatable:: sym_BZ(:,:),pt_TF(:)
        integer :: ni


  n_spinor=4
  allocate(Err(n_spinor))
  allocate(r_hat(n_spinor,3))
  allocate(sub_skip(n_spinor))
  allocate(read_angle(3*n_spinor))

  pi=3.1415926d0
  n_file=128
  n_kind=1
  n_iter=24001
  lambda=0.0005d0
  lambda_H=0.0d0
  lambda_B=0.002d0
  lambda_DM=0.00005d0
  act='run'


  open(unit=101,file='./sym_BZ',action='read'&
                &,access='sequential',form='unformatted')
        read(101) ni
        print*, ni
        allocate(read_VA(ni,4))
        allocate(VAs(n_file,n_kind,ni,4))
        allocate(pt_TF(ni))
  close(101)

  pt_TF(:)=1

  allocate(spin(n_file,n_kind,3*n_spinor))

  allocate(r_vec(n_file,n_kind,n_spinor,3))
  allocate(VA_ave(ni,4))

  allocate(V0(ni,3))
  allocate(Vn(ni,8))
  allocate(JH(ni,n_spinor))

  allocate(DV(ni,4))
  allocate(cal_V(ni,4))

  allocate(B_general(ni,3,n_spinor,3)) ! nR, i_dir, j_dir, i2_dir


  if (act .eq. 'plot' )then
        call do_plot()
        stop
  endif


  Do i_file=1,n_file
        print*,i_file
        call int2char(i_file,file_char)
        path='./angle_'//trim(file_char)//'/spin_info'

        open(unit=102,file=path,access='sequential',action='read')
        do i_kind=1,n_kind
                read(102,'(12F12.6)') read_angle(:)
                spin(i_file,i_kind,:)=read_angle(:)
                
                do i =1,n_spinor
                        r_vec(i_file,i_kind,i,1)=spin(i_file,i_kind,(i-1)*3+1)
                        r_vec(i_file,i_kind,i,2)=spin(i_file,i_kind,(i-1)*3+2)
                        r_vec(i_file,i_kind,i,3)=spin(i_file,i_kind,(i-1)*3+3)
                enddo
        enddo
        close(102)
        Do i_kind=1,n_kind
                call int2char(i_kind,kind_char)
                path='./angle_'//trim(file_char)//'/V'//trim(kind_char)
                open(unit=101,file=trim(path)&
                        &,access='sequential',action='read'&
                        &,form='unformatted')
                        read(101) i
                        read(101) read_VA
                close(101)
                VAs(i_file,i_kind,:,:)=read_VA(:,:)
        enddo

  enddo




  VA_ave(:,:)=0.0d0

  Do  i_file=1,n_file
  Do  i_kind=1,n_kind
        VAs(i_file,i_kind,:,:)=VAs(i_file,i_kind,:,:)-VA_ave(:,:)
  enddo
  enddo

  V0(:,:)=0.0001d0
  Vn(:,:)=0.0001d0
  B_general(:,:,:,:)=0.001d0
  goto 300 ! remove if you want to restart from a existing one
  open(unit=102,file='V_func_tmp',action='read',access='sequential',form='unformatted')
        read(102)
        read(102) V0
        read(102) Vn
        read(102) B_general
        read(102) JH
  close(102)

  300 continue


  Do i_iter=1,n_iter
        Err(:)=0.0d0
        i_count=0
        DO i_file=1,n_file
        DO i_kind=1,n_kind
                i_count=i_count+1
                r_hat(:,:)=r_vec(i_file,i_kind,:,:)
                call V_func_value(r_hat(:,:),cal_V)
                DV(:,:)=cal_V(:,:)-VAs(i_file,i_kind,:,:)
                
                Err(1)=Err(1)+sum(abs(DV(:,1)*pt_TF(:))**2)
                Err(2)=Err(2)+sum(abs(DV(:,2)*pt_TF(:))**2)
                Err(3)=Err(3)+sum(abs(DV(:,3)*pt_TF(:))**2)
                Err(4)=Err(4)+sum(abs(DV(:,4)*pt_TF(:))**2)

        
                i=0
                Do j_dir=1,n_spinor
                Do i2_dir=1,3
                        if (j_dir .eq. 1 .and. i2_dir .eq. 1) then
                                if ( mod(i_iter,n_spinor*3*5)  .lt. 5) then
                                do i_dir=1,3
                                         l_rate=1.0*wt
                                        B_general(:,i_dir,j_dir,i2_dir)=&
                                                &B_general(:,i_dir,j_dir,i2_dir)-&
                                                &lambda_DM*DV(:,i_dir+1)*r_hat(j_dir,i2_dir)*l_rate
                                


                                enddo
                                goto 700
                                endif
                        
                        elseif(mod(i_iter,n_spinor*3*5) .ge. i*5 .and. mod(i_iter,n_spinor*3*5)  .lt. (i+1)*5 ) then
                                do i_dir=1,3
                                        l_rate=1.0*wt
                                        B_general(:,i_dir,j_dir,i2_dir)=&
                                                &B_general(:,i_dir,j_dir,i2_dir)-&
                                                &lambda_DM*DV(:,i_dir+1)*r_hat(j_dir,i2_dir)*l_rate
                                enddo
                                goto 700
                        
                        endif


                        i=i+1

                enddo
                enddo

                700 continue




        enddo
        enddo

        Err(:)=Err(:)/i_count*n_file*n_kind
        print'(I8,3E24.10)',i_iter,Err(2)**0.5/(sum(abs(VAs(:,:,:,2))**2))**0.5,&
                &Err(3)**0.5/(sum(abs(VAs(:,:,:,3))**2))**0.5,&
                &Err(4)**0.5/(sum(abs(VAs(:,:,:,4))**2))**0.5
        if (mod(i_iter,n_spinor*3*5*2*5) .eq. n_spinor*3*5*2*5-1) then
        open(unit=102,file='V_func_tmp',action='write',access='append',status='replace',form='unformatted')
        write(102) VA_ave
        write(102) V0
        write(102) Vn
        write(102) B_general
        write(102) JH
        close(102)
        endif


ENDDO




contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )            WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )   WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )   WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )  WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char

subroutine V_func_value(angle_in,V_out)
        implicit none
        real(kind(0d0)), intent(in) :: angle_in(n_spinor,3)
        real(kind(0d0)), intent(out):: V_out(ni,4)
        real(kind(0d0)) ::r(n_spinor,3)
        integer ::i_site,i_dir,j_dir,i2_dir

r(:,:)=angle_in(:,:)



V_out(:,:)=0d0

do i_dir =1,3
do j_dir =1,n_spinor
do i2_dir=1,3
        V_out(:,i_dir+1)=V_out(:,i_dir+1)+B_general(:,i_dir,j_dir,i2_dir)*r(j_dir,i2_dir)
enddo
enddo
enddo


end subroutine V_func_value



subroutine do_plot()
        implicit none
        real(kind(0d0)):: V_plot(ni,4)
        real(kind(0d0))::angle_in(n_spinor,3)

        open(unit=102,file='V_func_tmp',action='read',access='sequential',form='unformatted')
        read(102) VA_ave
        read(102) V0
        read(102) Vn
        read(102) B_general
        close(102)



        !!!set your own, following is a FM case
        angle_in(1,:)=(/0.0d0,0.0d0,7.137257d0/)
        angle_in(2,:)=(/0.0d0,0.0d0,7.137257d0/)
        angle_in(3,:)=(/0.0d0,0.0d0,7.137257d0/)
        angle_in(4,:)=(/0.0d0,0.0d0,7.137257d0/)




        call  V_func_value(angle_in,V_plot)
        V_plot(:,:)= V_plot(:,:)+VA_ave(:,:)
        open(unit=103,file='V_fit',action='write',access='append'&
                &,form='unformatted',status='replace')
        write(103) ni
        write(103) V_plot
        close(103)

        



end subroutine do_plot


subroutine make_block()
        ! Not in use yet
        implicit none
        integer :: i_sym,origin(3),unit_vec(3,3)
        integer :: ix,iy,iz,jx,jy,jz,vec(ni,3)



        open(unit=101,file='./angle_1/V_vrs_dp',access='sequential',action='read'&
                &,form='unformatted')
                read(101) nx,ny,nz
        close(101)
        open(unit=101,file='./R_vec',action='read'&
                &,access='sequential',form='unformatted')
                read(101) vec
        close(101)


        pt_TF(:)=1
        do i_sym=1,6
        !print*,i_sym
        if (i_sym .eq. 1) then
                origin(:)=(/0,0,0/)
                unit_vec(1,:)=(/1,0,0/)
                unit_vec(2,:)=(/0,1,0/)
                unit_vec(3,:)=(/0,0,1/)
        elseif (i_sym .eq. 2) then
                origin(:)=(/nx/4,0,0/)
                unit_vec(1,:)=(/1,1,0/)
                unit_vec(2,:)=(/-1,0,0/)
                unit_vec(3,:)=(/0,0,1/)
        elseif (i_sym .eq. 3) then
                origin(:)=(/nx/2,nx/4,0/)
                unit_vec(1,:)=(/0,1,0/)
                unit_vec(2,:)=(/-1,-1,0/)
                unit_vec(3,:)=(/0,0,1/)
        elseif (i_sym .eq. 4) then
                origin(:)=(/nx/2,ny/2,0/)
                unit_vec(1,:)=(/-1,0,0/)
                unit_vec(2,:)=(/0,-1,0/)
                unit_vec(3,:)=(/0,0,1/)
        elseif (i_sym .eq. 5) then
                origin(:)=(/nx/4,ny/2,0/)
                unit_vec(1,:)=(/-1,-1,0/)
                unit_vec(2,:)=(/1,0,0/)
                unit_vec(3,:)=(/0,0,1/)
        elseif (i_sym .eq. 6) then
                origin(:)=(/0,ny/4,0/)
                unit_vec(1,:)=(/0,-1,0/)
                unit_vec(2,:)=(/1,1,0/)
                unit_vec(3,:)=(/0,0,1/)
        endif

         Do i=1,ni
                jx=mod(vec(i,1),nx/2)
                jy=mod(vec(i,2),ny/2)
                jz=mod(vec(i,3),nz)
                ix=origin(1)+sum((/jx,jy,jz/)*unit_vec(:,1))
                iy=origin(2)+sum((/jx,jy,jz/)*unit_vec(:,2))
                iz=origin(3)+sum((/jx,jy,jz/)*unit_vec(:,3))
                if (ix .le. 0 .or. ix .ge. nx) pt_TF(i)=0
                if (iy .le. 0 .or. ix .ge. ny) pt_TF(i)=0
                if (iz .le. 0 .or. ix .ge. nz) pt_TF(i)=0
        ENDDO

        enddo


end subroutine make_block


end program

