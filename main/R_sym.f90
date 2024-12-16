program R_sym
  
        ! This code creates the mapping from position index to 3D vector
        ! crystal symmetry operation can be included,
        ! a rotation symmetry making z -> -z in presented here 

        ! Note: Rather than using full BZ, we just take a relavent subspace
        !       countted as ni-points

        implicit none
        integer :: ix,iy,iz,jx,jy,jz,nx,ny,nz,n_sym
        integer,allocatable :: map(:)
        integer,allocatable :: sym_BZ(:,:),vec(:,:)
        integer :: i,j,ni
        integer :: orig(3)


   n_sym=2


   open(unit=101,file='./angle_1/V_HYC_dp',access='sequential',action='read'&
                &,form='unformatted')
        read(101) nx,ny,nz
        print*,'nx ny nz=',nx,ny,nz
   close(101)



  !! In this simple case, the subspace is chosen by removing the boundary
  !! This has no special meaning. 
  !! We just want to keep a space for future development
  i=0
  do ix=0,nx-1
  do iy=0,ny-1
  do iz=0,nz-1
        if(ix .eq. 0) cycle
        if(iy .eq. 0) cycle
        if(iz .eq. 0) cycle
        i=i+1
  enddo
  enddo
  enddo

  ni=i

  allocate(sym_BZ(ni,2))
  allocate(vec(ni,3))
  print*,'reduce ',nx*ny*nz, 'points', 'to' i

  !! mapping index to vector
  !! index for Id-sym
  i=1
  do ix=0,nx-1
  do iy=0,ny-1
  do iz=0,nz-1
        if(ix .eq. 0) cycle
        if(iy .eq. 0) cycle
        if(iz .eq. 0) cycle
        j=1+ix+nx*iy+nx*ny*iz
        sym_BZ(i,1)=j
        vec(i,:)=(/ix,iy,iz/)
        i=i+1
  enddo
  enddo
  enddo


  !! index for rotation symmetry making z -> -z
  orig(:)=(/ 0, 0 , nz /)
  do i=1,ni
        ix=vec(i,1)
        iy=vec(i,2)
        iz=vec(i,3)
        
        jx= iy   +orig(1)
        jy= ix   +orig(2)
        jz=-iz   +orig(3)

        j=1+jx+nx*jy+nx*ny*jz
        sym_BZ(i,2)=j
  enddo
  





  open(unit=101,file='./sym_BZ',action='write'&
                &,status='replace',form='unformatted')
        write(101) ni
        write(101) sym_BZ
  close(101)


  open(unit=101,file='./R_vec',action='write'&
                &,status='replace',form='unformatted')
        write(101) vec
  close(101)



end program
