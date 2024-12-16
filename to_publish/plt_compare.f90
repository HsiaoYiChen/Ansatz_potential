program plt_compare

        ! This code plot the predicted potential from fit.f90 
        ! (remember to set act=plot and assign the spin orientation you need)
        ! along with a chosen DFT potential
        ! some redundent variables are for future use 
        implicit none
        integer :: nx,ny,nz,i
        integer :: nx0,ny0,nz0
        integer :: ix,iy,iz,jx,jy,jz
        integer :: i_sym,ni,iz_site,iy_site,i_site
        real(kind(0d0)),allocatable :: read_V(:,:)
        real(kind(0d0)),allocatable :: V_write(:,:,:)
        real(kind(0d0)),allocatable :: V1(:,:)
        real(kind(0d0)),allocatable :: V_read(:,:)
        character(len=20) :: file_char,i_sym_char
        real(kind(0d0))    :: read4(4),pi
        integer :: z_box,io
        real:: r
        real    :: angle_file(12),to_write(16)
        real    :: angle_central(2),angle_near(12),angle_top(2)
        integer,allocatable :: sym_BZ(:,:)
        integer,allocatable :: vec(:,:)





open(unit=101,file='./angle_1/V_vrs_dp',access='sequential',action='read'&
                &,form='unformatted')
        read(101) nx,ny,nz
        print*,'nx ny nz=',nx,ny,nz
close(101)


open(unit=101,file='./V_fit',action='read'&                
        &,access='sequential',form='unformatted')
        read(101) ni
        print*,ni
        allocate(V_write(ni,4,2))
        allocate(vec(ni,3))
        read(101) V_write(:,:,1)
close(101)

!! chose the V_vrs you want to compare
open(unit=101,file='./angle_1/V_vrs_dp',action='read'&                
        &,access='sequential',form='unformatted')
        read(101) ni
        print*,ni
        read(101) V_write(:,:,2)
close(101)




open(unit=101,file='./R_vec',action='read'&
                &,access='sequential',form='unformatted')
        read(101) vec
close(101)


open(unit=102,file='./pt.dat_compare',action='write'&
        &,access='append',status='replace')
do i=1,ni
        write(102,'(6E16.8)')  V_write(i,2:4,1),V_write(i,2:4,2)
enddo
close(102)

contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )            WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000))WRITE(out_char,'(I5)')in_int

end subroutine int2char

end program

