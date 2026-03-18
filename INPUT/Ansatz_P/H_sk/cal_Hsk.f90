program spin
        implicit none
        integer :: ix,iy,iz,nx,ny,nz,i,is,j,io,js
        integer :: period_x,period_y,period_z,jx,jy,jz
        integer :: ini_x,ini_y,ini_z,fin_x,fin_y,fin_z
        integer :: iRVR,jRVR,H_dim,i_dim,j_dim,iwan,jwan
        integer :: super_x,super_y,super_z,jrept
        integer,allocatable :: cell(:,:),super_ir_vec(:,:)
        real(kind(0d0)),allocatable        :: x_list(:,:,:),coefficient(:)
        real(kind(0d0)),allocatable        :: spin_vec(:,:,:),off_set(:,:)
        real(kind(0d0))         :: A2B(3),r_tmp(3),Mz
        real(kind(0d0))    :: rad,pi,x,y_com,x_com,theta
        real(kind(0d0))    :: a(3,3),a0 ,avB(3),avA(3),rrBdr(3),rrAdr(3)     
        real(kind(0d0)) :: a_next,b_next,s_abs,celldm3
        real(kind(0d0)) ::vec6(6,3),vec3(3),axial
        complex(kind(0d0)),allocatable  :: rot_U_all(:,:),Hd(:,:),Hp(:,:)
        complex,allocatable  :: rot_U_sp(:,:)
        complex(kind(0d0)):: rot_U(2,2),rot_H_tmp(2,2),Tmp22(2,2)
        real(kind(0d0)) ::rot_eig(2)

        integer :: int3(3),super_nx,super_ny,super_nz,n_func,convert_size_x,convert_size_y,convert_size_z
        integer :: irept,nrept,super_nrept,i_locat,j_locat
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j,d_ext
        integer :: nspin, nwan_sp, nwan, wan_max,wan_min,i_dir,j_dir,k_dir
        integer,allocatable     :: dir_vec(:,:),inv_dir(:)
        character(len=200)      :: useless_string,func_dir,nscf_dir,THE_dir
        character(len=200)      :: RVR_dir
        character(len=20)       :: rept_char
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_POT_ave(:,:),RVR_tmp(:,:)
        complex(kind(0d0)),allocatable  :: H_sk(:,:,:)
        complex(kind(0d0))      :: cI

open(unit=101,file='../../THE.in',action='read',access='sequential',status='old')
read(101,'(A9,A200)')useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,'(A7,I3)')  useless_string,n_func
read(101,'(A15,I2)') useless_string,convert_size_x
read(101,'(A15,I2)') useless_string,convert_size_y
read(101,'(A15,I2)') useless_string,convert_size_z
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I2)') useless_string, nspin
read(101,*)!'(A3,I3)') useless_string, nx
read(101,*)!'(A3,I3)') useless_string, ny
read(101,*)!'(A3,I3)') useless_string, nz
read(101,'(A6,F8.4)')  useless_string,axial
read(101,'(A9,I3)') useless_string,period_x
read(101,'(A9,I3)') useless_string, period_y
read(101,'(A9,I3)') useless_string, period_z
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,'(A8,I3)') useless_string, wan_min
read(101,'(A8,I3)') useless_string, wan_max
read(101,'(A3,F16.8)') useless_string, a0
read(101,'(A2,F16.8)') useless_string, celldm3
close(101)

call make_spin()
call load_super_size()
if  (nspin .eq. 1) then
        nwan_sp=2*(wan_max-wan_min+1)
else
        nwan_sp=(wan_max-wan_min+1)
endif

H_dim=nwan*2/nspin*period_x*period_y*period_z
allocate(H_sk(H_dim,H_dim,super_nrept))
allocate(rot_U_all(H_dim,H_dim))
allocate(rot_U_sp(H_dim,H_dim))

d_ext=(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
RVR_dir=trim(THE_dir)//'./1_make_RVR/'



call read_RVR()


allocate(coefficient(n_func))

H_sk(:,:,:)=dcmplx(0.0d0,0.0d0)

do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
        i=1+iz+nz*iy+ny*nz*ix
       ! print*,i,nx*ny*nz
        coefficient(:)=0.0d0
        do i_dir=1,d_ext
!                if (.not. all(dir_vec(i_dir,:) .eq. 0)) cycle
                ini_x=cell(i,1)+dir_vec(i_dir,1)
                ini_y=cell(i,2)+dir_vec(i_dir,2)
                ini_z=cell(i,3)+dir_vec(i_dir,3)
               ! if (.not. (( ini_y .eq. 1  .and.  ini_x .eq. 0) .or.&
               !         &( ini_y .eq. 0  .and.  ini_x .eq. 1))) cycle
               ! if (.not. ini_x .eq. 0) cycle
  !              if (.not. ini_z .eq. 0) cycle
  !              if (.not. (all(cell(i,:) .eq. (/0,0,0/)) .or.all(cell(i,:) .eq.(/0,3,0/)))) cycle
                if (.not. (ini_x .ge.0 .and. ini_x/period_x .eq. 0)) cycle
                if (.not. (ini_y .ge.0 .and. ini_y/period_y .eq. 0)) cycle
                if (.not. (ini_z .ge.0 .and. ini_z/period_z .eq. 0)) cycle
                i_locat=1+ini_z+period_z*ini_y+period_z*period_y*ini_x

                coefficient(1:2)=spin_vec(i,1:2,1)
                coefficient(3:4)=spin_vec(i,1:2,2)
                coefficient(5:6)=spin_vec(i,1:2,3)
                coefficient(7:8)=spin_vec(i,1:2,4)
                coefficient(9:10)=spin_vec(i,1:2,5)
                coefficient(11:12)=spin_vec(i,1:2,6)
                coefficient(13:14)=spin_vec(i,1:2,7)
                coefficient(15:16)=spin_vec(i,1:2,8)
                coefficient(17:18)=spin_vec(i,1:2,9)
                coefficient(19:20)=spin_vec(i,1:2,10)

                coefficient(21:22)=spin_vec(i,1:2,11)
                coefficient(23:24)=spin_vec(i,1:2,12)
                
                coefficient(25:26)=spin_vec(i,1:2,13)
                coefficient(27:28)=spin_vec(i,1:2,14)
                coefficient(29:30)=spin_vec(i,1:2,15)
                coefficient(31:32)=spin_vec(i,1:2,16)
                coefficient(33:34)=spin_vec(i,1:2,17)
                coefficient(35:36)=spin_vec(i,1:2,18)
                coefficient(37:38)=spin_vec(i,1:2,19)
                coefficient(39:40)=spin_vec(i,1:2,20)
                coefficient(41:42)=spin_vec(i,1:2,21)
                coefficient(43:44)=spin_vec(i,1:2,22)

                coefficient(45:46)=spin_vec(i,1:2,23)
                coefficient(47:48)=spin_vec(i,1:2,24)









!print*,coefficient(1:27)

                
                do j_dir=1,nwan_sp*d_ext
                        RVR_tmp(:,j_dir)=matmul(coefficient(:),RVR(:,:,j_dir))
                enddo


        do j_dir=1,d_ext
!if (.not. all(dir_vec(j_dir,:) .eq. 0)) cycle 
!               if (.not. (dir_vec(j_dir,3) .eq. 0)) cycle
                fin_x=cell(i,1)+dir_vec(j_dir,1)
                fin_y=cell(i,2)+dir_vec(j_dir,2)
                fin_z=cell(i,3)+dir_vec(j_dir,3)
                !print*,fin_x,fin_y,fin_z
                if (fin_x .ge. 0) then
                        super_x=fin_x/period_x
                else
                        super_x=(fin_x+1)/period_x-1  
                endif
                if (fin_y .ge. 0) then
                        super_y=fin_y/period_y
                else
                        super_y=(fin_y+1)/period_y-1
                endif
                if (fin_z .ge. 0) then
                        super_z=fin_z/period_z
                else
                        super_z=(fin_z+1)/period_z-1
                endif
 



                do irept=1,super_nrept
                        if (all((/super_x,super_y,super_z/) .eq. (super_ir_vec(irept,:)))) then
                                jrept=irept
                                exit
                        endif
                enddo


                
                !print*,fin_x,fin_y,fin_z,super_x,super_y,super_z,jrept

                jx=fin_x-super_x*period_x
                jy=fin_y-super_y*period_y
                jz=fin_z-super_z*period_z
                j_locat=1+jz+period_z*jy+period_z*period_y*jx

       ! print*,j_locat

        do iwan=wan_min,wan_max
        do jwan=wan_min,wan_max
        do is=1,2/nspin
        do js=1,2/nspin
                !if (iwan .gt. 28) cycle
                !if (jwan .gt. 28) cycle
                iRVR=1+(i_dir-1)+d_ext*((iwan-wan_min)*2+is-1)
                jRVR=1+(j_dir-1)+d_ext*((jwan-wan_min)*2+js-1)
                i_dim=(is+2/nspin*(iwan-1))+nwan*2/nspin*(i_locat-1)
                j_dim=(js+2/nspin*(jwan-1))+nwan*2/nspin*(j_locat-1)
                
                !if (abs(RVR_tmp(iRVR,jRVR)) .lt. 10.0**(-4)) cycle
            !    if (.not. jrept .eq. 536) cycle
                !print'(5I8,2E16.8)',i_dim,j_dim,ix,iy,iz,RVR_tmp(iRVR,jRVR)
                H_sk(i_dim,j_dim,jrept)=H_sk(i_dim,j_dim,jrept)+RVR_tmp(iRVR,jRVR)

        enddo
        enddo
        enddo
        enddo
        
        enddo!j_dir
        enddo!i_dir
enddo!iz
enddo!iy
enddo!ix




!do i_dim=11,20
!        print*,H_sk(11,i_dim,23),H_sk(11+nwan*2/nspin*3,i_dim+nwan*2/nspin*3,23)
!
!enddo

!stop

allocate(Hd(10,10))
!allocate(Hp(6,6))


do iwan=1,5
!        H_sk(10+1+(iwan-1)*2:10+iwan*2,10+1+(iwan-1)*2:10+iwan*2,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,:)
!        H_sk(20+1+(iwan-1)*2:20+iwan*2,20+1+(iwan-1)*2:20+iwan*2,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,:)
!        H_sk(30+1+(iwan-1)*2:20+iwan*2,30+1+(iwan-1)*2:20+iwan*2,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,:)
enddo



!H_sk(1:40,41:128,:)=0.0
!H_sk(41:128,1:40,:)=0.0
!H_sk(41:128,41:128,:)=0.0



do irept=1,super_nrept


    !    do iwan=41,128
    !    do jwan=41,128
    !            if (abs(H_sk(iwan,jwan,irept)) .lt. 0.01) then
    !                    H_sk(iwan,jwan,irept)=0.0
    !                    H_sk(jwan,iwan,irept)=0.0
    !                    
    !            endif
    !    enddo
    !    enddo


        !if (.not. all (super_ir_vec(irept,:) .eq. (/0,0,0/))) then
        !        H_sk(:,:,irept)=0.0
        !else
              !  Hd(:,:)=H_sk(1:10,1:10,irept)
       !         H_sk(:,:,irept)=0.0
       !         H_sk(1:10,1:10,irept)=Hd(:,:)                
       !         H_sk(11:20,11:20,irept)=Hd(:,:)
       !         H_sk(21:30,21:30,irept)=Hd(:,:)
       !         H_sk(31:40,31:40,irept)=Hd(:,:)

       ! endif


        call int2char(irept,rept_char)
        open(unit=102,file='./data/Hsk_'//trim(rept_char),action='write',access='append'&
                &,status='replace',form='unformatted')
        write(102)      H_sk(:,:,irept)
        close(102)
        !if (.not. mod(irept,10) .eq. 8) cycle
        print*,irept,sum(abs(H_sk(:,:,irept))**2),super_ir_vec(irept,:)
 !       print'(I8,4F12.6)',irept,real(H_sk(1,2,irept)),real(H_sk(11,12,irept)),real(H_sk(21,22,irept)),real(H_sk(31,32,irept))

!        do iwan=1,40
!        print*,H_sk(iwan,iwan,irept)
!        enddo




        rot_U_all(:,:)=dcmplx(0.0d0,0.0d0)
        do iwan=1,2*nwan/nspin*period_x*period_y*period_z
                rot_U_all(iwan,iwan)=dcmplx(1.0d0,0.0d0)

        enddo

        if (all (super_ir_vec(irept,:) .eq. (/0,0,0/))) then
                print'(4E18.6)', H_sk(1,11,irept),H_sk(1,12,irept),H_sk(11,1,irept),H_sk(12,1,irept)
        endif
        if (all (super_ir_vec(irept,:) .eq. (/0,-1,0/))) then
                print'(4E18.6)', H_sk(1,11,irept),H_sk(1,12,irept),H_sk(11,1,irept),H_sk(12,1,irept)
        endif



        if (all (super_ir_vec(irept,:) .eq. (/0,0,0/))) then
        do iwan=1,20,5!nwan/nspin*period_x*period_y*period_z
                rot_H_tmp(:,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,irept)

                
                print'(I8,8E18.6)',iwan, rot_H_tmp(:,:)



                call diasym(rot_H_tmp,rot_eig,2,rot_U)
                rot_U_all(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2)=rot_U(:,:)
                
                
                if(iwan .eq. 20) then

                rot_H_tmp(:,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,irept)

                !print*,rot_eig(1),rot_eig(2)
                !print*,''
                !print*,rot_H_tmp(1,1),rot_H_tmp(1,2)
                !print*,rot_H_tmp(2,1),rot_H_tmp(2,2)
                !print*,''
                !print*,rot_U(1,1),rot_U(1,2)
                !print*,rot_U(2,1),rot_U(2,2)

                endif
                
                
                rot_H_tmp(:,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,irept)

                Tmp22(:,:)=matmul(rot_H_tmp(:,:),rot_U(:,:))
        
                rot_H_tmp(:,:)=matmul(transpose(dconjg(rot_U(:,:))),Tmp22(:,:))
 !               print*,rot_H_tmp(:,:)
        
        enddo
        

        rot_U_sp(:,:)=rot_U_all(:,:)

        open(unit=102,file='./data/rot_U',action='write',access='append'&
                &,status='replace',form='unformatted')
        write(102)      rot_U_sp(:,:)
        close(102)




        endif
        
        
        
enddo





contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )   WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   ) WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  ) WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) ) WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine make_spin()
        complex(kind(0d0)) :: s1(3),s2(3),s3(3),sz(3)
        real(kind(0d0)) :: r_tmp(3),a_super(3,3),Q(3,3)
        real(kind(0d0)) :: mzz,rT,rT2,rT3,rT4
        real(kind(0d0)) :: theta
max_x=-99
max_y=-99
max_z=-99
min_x=99
min_y=99
min_z=99

open(unit=102,file='../super_H/ir_vec_HYC',action='read',access='sequential',status='old')
do 
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        if (int3(1) .gt. max_x) max_x=int3(1)
        if (int3(2) .gt. max_y) max_y=int3(2)
        if (int3(3) .gt. max_z) max_z=int3(3)
        if (int3(1) .lt. min_x) min_x=int3(1)
        if (int3(2) .lt. min_y) min_y=int3(2)
        if (int3(3) .lt. min_z) min_z=int3(3)
enddo
close(102)

print*,max_x,max_y,max_z
print*,min_x,min_y,min_z



nx=(max_x-min_x+1)*period_x
ny=(max_y-min_y+1)*period_y
nz=(max_z-min_z+1)*period_z


print*, nx,ny,nz

pi=3.1415926535d0
cI=dcmplx(0.0d0,1.0d0)
allocate(x_list(nx*ny*nz,3,1))
allocate(cell(nx*ny*nz,3))
allocate(spin_vec(nx*ny*nz,3,24))
allocate(off_set(24,3))




               a(1,:) = (/   1.000000,   0.000000 ,  0.000000 /)
               a(2,:) = (/   0.000000,   0.504533 ,  0.000000 /)
               a(3,:) = (/   0.000000,   0.000000 ,  2.169706 /)


a_super(1,:)=a(1,:)/period_x
a_super(2,:)=a(2,:)/period_y
a_super(3,:)=a(3,:)/period_z

a(:,:)=a(:,:)!*a0*0.529177

i=1
do ix=0,nx-1!min_x*period_x,(max_x+1)*period_x-1
do iy=0,ny-1!min_y*period_y,(max_y+1)*period_y-1
do iz=0,nz-1!min_z*period_z,(max_z+1)*period_z-1
        i=1+iz+nz*iy+nz*ny*ix
        cell(i,:)=(/ix+min_x*period_x,iy+min_y*period_y,iz+min_z*period_z/)
       ! rint*,i,cell(i,:)
enddo
enddo
enddo



off_set(:,:)=0.0



off_set(1,:)=(/0.1219643354       , 0.2500000000-1.0       , 0.1522172689/)
off_set(2,:)=(/0.1219643354       , 0.2500000000+0.0       , 0.1522172689/)
off_set(3,:)=(/0.1219643354       , 0.2500000000+1.0       , 0.1522172689/)

off_set(4,:)=(/0.3780356646       , 0.7500000000-1.0       , 0.3477827311/)
off_set(5,:)=(/0.3780356646       , 0.7500000000+0.0       , 0.3477827311/)
off_set(6,:)=(/0.3780356646       , 0.7500000000+1.0       , 0.3477827311/)

off_set(7,:)=(/0.8780356646       , 0.7500000000-1.0       , 0.6522172689/)
off_set(8,:)=(/0.8780356646       , 0.7500000000+0.0       , 0.6522172689/)
off_set(9,:)=(/0.8780356646       , 0.7500000000+1.0       , 0.6522172689/)

off_set(10,:)=(/0.6219643354       , 0.2500000000-1.0       , 0.8477827311/)
off_set(11,:)=(/0.6219643354       , 0.2500000000+0.0       , 0.8477827311/)
off_set(12,:)=(/0.6219643354      ,  0.2500000000+1.0      ,  0.8477827311/)

off_set(13,:)=(/0.1219643354+1.0       , 0.2500000000-1.0       , 0.1522172689/)
off_set(14,:)=(/0.1219643354+1.0  ,      0.2500000000      ,  0.1522172689/)
off_set(15,:)=(/0.1219643354+1.0  ,      0.2500000000+1.0  ,      0.1522172689/)

off_set(16,:)=(/0.3780356646+1.0       , 0.7500000000-1.0       , 0.3477827311/)
off_set(17,:)=(/0.3780356646+1.0       , 0.7500000000+0.0       , 0.3477827311/)
off_set(18,:)=(/0.3780356646+1.0       , 0.7500000000+1.0       , 0.3477827311/)




off_set(19,:)=(/0.8780356646-1.0       , 0.7500000000-1.0       , 0.6522172689/)
off_set(20,:)=(/0.8780356646-1.0       , 0.7500000000+0.0       , 0.6522172689/)
off_set(21,:)=(/0.8780356646-1.0       , 0.7500000000+1.0       , 0.6522172689/)

off_set(22,:)=(/0.6219643354-1.0       , 0.2500000000-1.0       , 0.8477827311/)
off_set(23,:)=(/0.6219643354-1.0       , 0.2500000000+0.0       , 0.8477827311/)
off_set(24,:)=(/0.6219643354-1.0      ,  0.2500000000+1.0      ,  0.8477827311/)















         
         Q(2,:)=(/ 0.0d0,5.0d0/6,0.0d0/)
         Q(1,:)=(/ 0.5d0,0.0d0  ,0.0d0/)
        !Q(3,:)=(/            0.0d0, 1.0d0,0.0d0/)
        !Q(1,:)=(/1.0d0,0.0d0,0.0d0/)
        !Q(2,:)=(/-0.5d0, 3.0d0**(0.5d0)/2,0.0d0/)
        !Q(3,:)=(/-0.5d0,-3.0d0**(0.5d0)/2,0.0d0/)

Q(2,:)=Q(2,:)*1.982030 !b2



        mzz=-1.0d0
Mz=0.0d0
do i=1,nx*ny*nz
        !if (.not. all((/cell(i,1),cell(i,3)/).eq. (/0,0/)))cycle
        do j=1,24

                r_tmp(1)=sum(a(:,1)*(cell(i,:)+off_set(j,:)))
                r_tmp(2)=sum(a(:,2)*(cell(i,:)+off_set(j,:)))
                r_tmp(3)=sum(a(:,3)*(cell(i,:)+off_set(j,:)))

                
  !              s1(1)=cos(2*pi*sum(r_tmp(:)*Q(1,:))+2*pi*sum(r_tmp(:)*Q(2,:)))
  !              s1(2)=sin(2*pi*sum(r_tmp(:)*Q(1,:))+2*pi*sum(r_tmp(:)*Q(2,:)))
  !              s1(3)=0.0
 



                r_tmp(1)=sum(a(:,1)*(cell(i,:)+off_set(j,:)))+20000.0
                r_tmp(2)=sum(a(:,2)*(cell(i,:)+off_set(j,:)))/0.504533+20000.0
                r_tmp(3)=sum(a(:,3)*(cell(i,:)+off_set(j,:)))/2.169706+2000.0


                rT=75!30.0
                rT2=0.0
                rT3=0.0
                rT4=0.0

                if (mod(r_tmp(2),6.0) .le. 0.5 ) then
                        theta=0.0
                elseif (mod(r_tmp(2),6.0) .le. 1.0 ) then
                        theta=150.0*1
                elseif (mod(r_tmp(2),6.0) .le. 1.5 ) then
                        theta=150.0*2
                elseif (mod(r_tmp(2),6.0) .le. 2.0 ) then
                        theta=150.0*3
                elseif (mod(r_tmp(2),6.0) .le. 2.5 ) then
                        theta=150.0*4
                elseif (mod(r_tmp(2),6.0) .le. 3.0 ) then
                        theta=150.0*5
                elseif (mod(r_tmp(2),6.0) .le. 3.5 ) then
                        theta=150.0*6
                elseif (mod(r_tmp(2),6.0) .le. 4.0 ) then
                        theta=150.0*7
                elseif (mod(r_tmp(2),6.0) .le. 4.5 ) then
                        theta=150.0*8
                elseif (mod(r_tmp(2),6.0) .le. 5.0 ) then
                        theta=150.0*9
                elseif (mod(r_tmp(2),6.0) .le. 5.5 ) then
                        theta=150.0*10
                elseif (mod(r_tmp(2),6.0) .le. 6.0 ) then
                        theta=150.0*11
                endif








                if (mod(r_tmp(1),2.0) .le. 0.5) then
                        theta=theta
                elseif (mod(r_tmp(1),2.0) .le. 1.0) then
                        theta=theta+90
                elseif (mod(r_tmp(1),2.0) .le. 1.5) then
                        theta=theta+180
                elseif (mod(r_tmp(1),2.0) .le. 2.0) then        
                        theta=theta+270
                endif



                !if (r_tmp(3) .gt. 0.5) then
                !        s1(1)=cos(0.0*2*pi/360)
                !        s1(2)=sin(0.0*2*pi/360)
                !else
                !        s1(1)=cos(-0.0*2*pi/360)
                !        s1(2)=sin(-0.0*2*pi/360)
                !endif


                        s1(1)=cos(theta*2*pi/360)
                        s1(2)=sin(theta*2*pi/360)
                        s1(3)=0.0




               ! -( Q(1,2))*cI*exp(cI*2*pi*sum(r_tmp(:)*Q(1,:)))
               ! s1(2)=-(-Q(1,1))*cI*exp(cI*2*pi*sum(r_tmp(:)*Q(1,:)))
               ! s1(3)= mzz*exp(cI*2*pi*sum(r_tmp(:)*Q(1,:)))


             !   s2(1)=-( Q(2,2))*cI*exp(cI*2*pi*sum(r_tmp(:)*Q(2,:)))
             !   s2(2)=-(-Q(2,1))*cI*exp(cI*2*pi*sum(r_tmp(:)*Q(2,:)))
             !   s2(3)= mzz*exp(cI*2*pi*sum(r_tmp(:)*Q(2,:)))

             !   s3(1)=-( Q(3,2))*cI*exp(cI*2*pi*sum(r_tmp(:)*Q(3,:)))
             !   s3(2)=-(-Q(3,1))*cI*exp(cI*2*pi*sum(r_tmp(:)*Q(3,:)))
             !   s3(3)= mzz*exp(cI*2*pi*sum(r_tmp(:)*Q(3,:)))


 !               s1(:)=-s1(:)
 !               s2(:)=-s2(:)
  !              s3(:)=-s3(:)
                
              !  sz(:)=(/1.0d0,0.0d0,axial*0.00/)*1000000
                !sz(:)=0.0
                !sx(:)=0.0
                !sy(:)=0.0

                

                spin_vec(i,1,j)=s1(1)!(s1(1)+dconjg(s1(1))+s2(1)+dconjg(s2(1))+s3(1)+dconjg(s3(1)))/2
                spin_vec(i,2,j)=s1(2)!(s1(2)+dconjg(s1(2))+s2(2)+dconjg(s2(2))+s3(2)+dconjg(s3(2)))/2
                spin_vec(i,3,j)=s1(3)!(s1(3)+dconjg(s1(3))+s2(3)+dconjg(s2(3))+s3(3)+dconjg(s3(3)))/2
               ! spin_vec(i,:,j)=0.0d0
               ! spin_vec(i,:,j)=spin_vec(i,:,j)+sz(:)

 


                s_abs=sum(spin_vec(i,:,j)**2)**0.5

                spin_vec(i,:,j)=spin_vec(i,:,j)/s_abs
                !spin_vec(i,2,j)=spin_vec(i,2,j)/s_abs
                !spin_vec(i,3,j)=spin_vec(i,3,j)/s_abs
                spin_vec(i,:,j)=spin_vec(i,:,j)*0.526d0
       ! print*,spin_vec(i,:,j),sum(r_tmp(:)*Q(1,:)),sum(r_tmp(:)*Q(2,:)),sum(r_tmp(:)*Q(3,:))


                Mz=Mz+spin_vec(i,3,j)



 !               print*,cI*2,pi,(cell(i,1)+off_set(j,1)),period_x
        enddo
enddo


Mz=Mz/(nx*ny*nz)
print*,Mz/2

open(unit=101,file='cmap',action='write',access='append',status='replace')
do i=1,nx*ny*nz
        if (.not. cell(i,3) .eq. 0) cycle
   !      if (.not. cell(i,1) .eq.0) cycle
     !   if (.not. cell(i,3) .le. 10) cycle
        do j=1,24
        write(101,'(6F24.6)') sum(a(:,1)*(cell(i,:)+off_set(j,:)))*15.53313653*0.529177-spin_vec(i,1,j)*1,spin_vec(i,1,j)*2&
        &,sum(a(:,2)*(cell(i,:)+off_set(j,:)))*15.53313653*0.529177-spin_vec(i,2,j)*1,spin_vec(i,2,j)*2&
        &,sum(a(:,3)*(cell(i,:)+off_set(j,:)))*15.53313653*0.529177-spin_vec(i,3,j)*1,spin_vec(i,3,j)*2
        enddo
enddo
close(101)

open(unit=101,file='cmap_rad',action='write',access='append',status='replace')
do i=1,nx*ny*nz
        if (.not. cell(i,3) .eq. 0) cycle
   !      if (.not. cell(i,1) .eq.0) cycle
     !   if (.not. cell(i,3) .le. 10) cycle
        do j=1,4
        !if (abs(sum(a(:,2)*(cell(i,:)+off_set(j,:)))) .lt. 10.0**(-3)) then
        write(101,'(2F24.6)') sum(a(:,1)*(cell(i,:)+off_set(j,:)))**2+sum(a(:,2)*(cell(i,:)+off_set(j,:)))**2,spin_vec(i,3,j)
        !endif        
        enddo
enddo
close(101)


!open(unit=101,file='spin_texture',action='write',access='append',form='unformatted',status='replace')
!write(101) spin_vec(:,:,:)
!close(101)

end subroutine




subroutine load_super_size()
        irept=0
        open(unit=102,file='../super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        irept=irept+1
        enddo
        close(102)
        super_nrept=irept
        allocate(super_ir_vec(super_nrept,3))

        irept=1
        open(unit=102,file='../super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
                super_ir_vec(irept,:)=int3(:)
                irept=irept+1
        enddo
        close(102)






end subroutine


subroutine read_RVR()


        allocate(RVR(n_func,nwan_sp*d_ext,nwan_sp*d_ext))
        allocate(RVR_POT_ave(nwan_sp*d_ext,nwan_sp*d_ext))
        allocate(RVR_tmp(nwan_sp*d_ext,nwan_sp*d_ext))

        open(unit=103,file=trim(RVR_dir)//'./B_general',action='read',access='sequential',form='unformatted')
                read(103) RVR
        close(103)

        !open(unit=103,file=trim(RVR_dir)//'./B_ave',action='read',access='sequential',form='unformatted')
        !        read(103) RVR_POT_ave
        !close(103)
        RVR_POT_ave(:,:)=0.0d0








        allocate(dir_vec((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1),3))
        allocate(inv_dir((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)))
        inv_dir(:)=0
        i=1
        do i_dir=-convert_size_x,convert_size_x
        do j_dir=-convert_size_y,convert_size_y
        do k_dir=-convert_size_z,convert_size_z
        
                dir_vec(i,:)=(/i_dir,j_dir,k_dir/)
                i=i+1
        enddo
        enddo
        enddo


        






        do i=1,(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
        do j=1,(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
                if (all( dir_vec(i,:)+dir_vec(j,:) .eq. (/0,0,0/))) then
                        inv_dir(i)=j
                        exit
                endif
        enddo
                if (inv_dir(i) .eq. 0) then
                        print*, i, "-dir inverse not found"
                        stop
                endif
        enddo




        !do i_dir=1,d_ext
        !do j_dir=1,d_ext
        !if (.not. all(dir_vec(i_dir,:) .eq. 0)) cycle
        !if (.not. (dir_vec(j_dir,3) .eq. 0)) cycle
        
        !do iwan=6,6!wan_min,wan_max
        !do jwan=6,6!wan_min,wan_max
        !do is=2,2!2/nspin
        !do js=2,2!2/nspin
        !        iRVR=1+(i_dir-1)+d_ext*((iwan-wan_min)*2+is-1)
        !        jRVR=1+(j_dir-1)+d_ext*((jwan-wan_min)*2+js-1)
        !print*, abs(RVR(3,jRVR,iRVR)),dir_vec(j_dir,:)
        !enddo
        !enddo
        !enddo
        !enddo
        !enddo
        !enddo

end subroutine read_RVR


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
