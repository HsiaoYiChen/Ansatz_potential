program rot

        !! This code aims for further symmetry operation
        !! which is still trivial in its current form,
        !! so some part of the code looks meaningless
        !! to be comined with extract.f90

        implicit none
        integer :: nx,ny,nz,n_file,i_file,i_iter,n_iter,i,j
        integer :: ix,iy,iz,jx,jy,jz
        integer :: i_sym,ni,iz_site,iy_site,i_site
        real(kind(0d0)),allocatable :: read_V(:,:)
        real(kind(0d0)),allocatable :: V_write(:,:)
        real(kind(0d0)),allocatable :: V1(:,:)
        character(len=20) :: file_char,i_sym_char
        real(kind(0d0))    :: read4(4),pi
        integer :: z_box,io,v3(3)
        real:: r
        real    :: angle_file(2,12),to_write(4,3)
        real    :: amp(4),theta(4),phi(4),s(4,3)
        integer :: pt_map(4),orig(3)
        real,allocatable        :: spin(:,:,:)

        integer,allocatable :: sym_BZ(:,:),vec(:,:)
        integer,allocatable :: rot_pt(:,:)

        pi=3.1415926d0
        n_file=2
        allocate(spin(n_file,4,3))

        open(unit=202,file='./angles',action='read',access='sequential',status='old')
        do i=1,2
                read(202,'(12F12.6)') angle_file(i,:)
        enddo
        close(202)




        open(unit=101,file='./V_vrs_dp',access='sequential',action='read'&
                &,form='unformatted')
                read(101) nx,ny,nz
                print*,'nx ny nz=',nx,ny,nz
        close(101)

        open(unit=101,file='../sym_BZ',action='read'&
                &,access='sequential',form='unformatted')
                read(101) ni
                allocate(V_write(ni,4))
                allocate(V1(ni,4))
                allocate(vec(ni,3))
                allocate(rot_pt(ni,1))
                read(101) sym_BZ
        close(101)

        open(unit=101,file='../R_vec',action='read'&
                &,access='sequential',status='old',form='unformatted')
                read(101) vec
        close(101)



        rot_pt(:,:)=0
        do i=1,ni
        do i_sym=1,1
                if (i_sym .eq. 1) then
                        v3(:)=vec(i,:)
                endif
                !   To be add
        enddo
        enddo





        do i_file=1,n_file
                call int2char(i_file,file_char)

                open(unit=103,file='./V'//trim(file_char),action='read'&
                        &, access='sequential',status='old',form='unformatted')
                        read(103) ni
                        read(103) V1
                close(103)


                do i=1,4
                        s(i,1)=angle_file(i_file,(i-1)*3+1)*sin(2*pi/360.0*angle_file(i_file,(i-1)*3+2))&
                                &*cos((2*pi/360.0*angle_file(i_file,(i-1)*3+3)))
                        s(i,2)=angle_file(i_file,(i-1)*3+1)*sin(2*pi/360.0*angle_file(i_file,(i-1)*3+2))&
                                &*sin((2*pi/360.0*angle_file(i_file,(i-1)*3+3)))
                        s(i,3)=angle_file(i_file,(i-1)*3+1)*cos(2*pi/360.0*angle_file(i_file,(i-1)*3+2))
                enddo


                
                do i_sym=1,1
                        if (i_sym .eq. 1) then
                                do i=1,4
                                        pt_map(i)=i
                                enddo
                        endif

                        if (i_sym .eq. 1) then
                                spin((i_sym-1)*n_file+i_file,:,:)=s(:,:)
                        endif


                !               if (i_sym .eq. 2) then
                !               to be add.... 

                enddo
        enddo



        open(unit=103,file='./spin_info',action='write'&
                &, access='append',status='replace')
        do i_sym=1,1
        do i_file=1,n_file
                write(103,'(12F12.6)') (spin((i_sym-1)*n_file+i_file,i,:),i=1,4)
        enddo
        enddo
        close(103)


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

