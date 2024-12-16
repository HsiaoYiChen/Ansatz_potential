program extract

        !! This code convert the V_vrs from QE to 
        !! our format in the real-space index of the subspace
        !! chosen by R_sym.f90, with symmetry operation.

        !! Note that the details need to be modify for other symmetry operation
        !! Note that some variables may looks redundant.
        !! They are kept for future extension.

        implicit none
        integer :: nx,ny,nz,n_file,i_file,i_iter,n_iter,i,j
        integer :: ix,iy,iz,jx,jy,jz,n_sym
        integer :: i_sym,ni,iz_site,iy_site,i_site
        real(kind(0d0)),allocatable :: read_V(:,:)
        real(kind(0d0)),allocatable :: V_write(:,:)
        real(kind(0d0)),allocatable :: V1(:,:)
        character(len=20) :: file_char,i_sym_char
        real(kind(0d0))    :: read4(4),pi
        integer :: z_box,io
        real:: r
        real    :: angle_file(12),to_write(4,3)
        real    :: amp(4),theta(4),phi(4)

        integer,allocatable :: sym_BZ(:,:)

  pi=3.1415926d0
  n_sym=2

  open(unit=101,file='./polar',action='read',access='sequential')
        read(101,'(12F12.6)') angle_file(:)
  close(101)

  do i=1,4
        amp(i)=angle_file((i-1)*3+1)
        theta(i)=angle_file((i-1)*3+2)
        phi(i)=angle_file((i-1)*3+3)
  enddo



  open(unit=101,file='./V_vrs_dp',access='sequential',action='read'&
                &,form='unformatted')
        read(101) nx,ny,nz
        print*,'nx ny nz=',nx,ny,nz
        allocate(read_V(nx*ny*nz,4))
        allocate(V1(nx*ny*nz,4))
        read(101) read_V
  close(101)

  open(unit=101,file='../sym_BZ',action='read'&
                &,access='sequential',form='unformatted')
        read(101) ni
        allocate(sym_BZ(ni,n_sym))
        allocate(V_write(ni,4))
        read(101) sym_BZ
  close(101)


open(unit=202,file='./angles',action='write',access='append',status='replace')
i_file=1
do iy_site=0,0
        i_site=1+(iy_site)

        do ix=0,nx-1
        do iy=0,ny-1
        do iz=0,nz-1
                i=1+ix+nx*iy+nx*ny*iz
                j=1+ix+nx*iy+nx*ny*iz
                V1(i,:)=read_V(j,:)
        enddo
        enddo
        enddo

        do i_sym=1,n_sym
                !!
                call int2char(i_file,file_char)

                if (i_sym .eq. 1) then
                do i=1,ni
                        V_write(i,:)=V1(sym_BZ(i,i_sym),:)
                enddo
                elseif(i_sym .eq. 2) then
                do i=1,ni
                        V_write(i,1)= V1(sym_BZ(i,i_sym),1)
                        V_write(i,4)=-V1(sym_BZ(i,i_sym),4)

                        V_write(i,2)=-V1(sym_BZ(i,i_sym),2)*cos(2*pi*60.0/360)+&
                                    &V1(sym_BZ(i,i_sym),3)*cos(2*pi*30.0/360)
                        V_write(i,3)=V1(sym_BZ(i,i_sym),3)*sin(2*pi*30.0/360)+&
                                    &V1(sym_BZ(i,i_sym),2)*sin(2*pi*60.0/360)

                enddo
                endif





                open(unit=103,file='./V'//trim(file_char),action='write'&
                        &, access='append',status='replace',form='unformatted')
                        write(103) ni
                        write(103) V_write
                close(103)


                if (iy_site .eq. 0 .and. i_sym .eq.1) then
                        to_write(1,1)=amp(1)
                        to_write(1,2)=theta(1)
                        to_write(1,3)=phi(1)

                        to_write(2,1)=amp(2)
                        to_write(2,2)=theta(2)
                        to_write(2,3)=phi(2)

                        to_write(3,1)=amp(3)
                        to_write(3,2)=theta(3)
                        to_write(3,3)=phi(3)

                        to_write(4,1)=amp(4)
                        to_write(4,2)=theta(4)
                        to_write(4,3)=phi(4)
                
                elseif(iy_site .eq. 0 .and. i_sym .eq.2) then

                        to_write(1,1)=amp(4)
                        to_write(1,2)=180-theta(4)
                        to_write(1,3)=mod(120.0-phi(4),360.0)

                        to_write(2,1)=amp(2)
                        to_write(2,2)=180-theta(2)
                        to_write(2,3)=mod(120.0-phi(2),360.0)

                        to_write(3,1)=amp(3)
                        to_write(3,2)=180-theta(3)
                        to_write(3,3)=mod(120.0-phi(3),360.0)

                        to_write(4,1)=amp(1)
                        to_write(4,2)=180-theta(1)
                        to_write(4,3)=mod(120.0-phi(1),360.0)

                endif


                write(202,'(12F12.6)') (to_write(i,:),i=1,4)

                i_file=i_file+1
        enddo


enddo



close(202)
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

