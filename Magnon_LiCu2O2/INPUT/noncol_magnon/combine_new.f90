program spin
        implicit none
        integer :: is,io,iks,nks,js
        integer :: period_x,period_y,period_z
        integer :: H_dim,i_dim,j_dim,iwan,jwan,i_super
        integer,allocatable :: super_ir_vec(:,:)
        real(kind(0d0))    :: pi
        integer :: int3(3)
        integer :: irept,super_nrept
        integer :: nspin, nwan_sp, nwan, wan_max,wan_min
        character(len=200)      :: useless_string,func_dir,nscf_dir,THE_dir
        character(len=20)       :: rept_char
        complex(kind(0d0)) :: cI
        complex(kind(0d0)),allocatable  :: H_sk(:,:),read_H(:,:),full_H(:,:)
        complex(kind(0d0)),allocatable  :: read_AAR(:,:,:,:),full_AAR(:,:,:,:)
        complex(kind(0d0)),allocatable  :: H_soc(:,:),H_FM(:,:)
        complex,allocatable  :: full_H_sp(:,:,:),full_AAR_sp(:,:,:,:)
        complex,allocatable  :: rot_U_all(:,:),U_tmp(:,:)
        complex ::rot_H_tmp(2,2),rot_U(2,2)
        complex,allocatable  :: Unmk_super(:,:,:),H_k_super(:,:)
        real,allocatable     :: eig_super(:,:),super_k_list(:,:)
        real    :: rot_eig(2),read3(3)
        logical ::SOC_TF
        integer,allocatable     :: super_deg(:)


pi=3.1415926535d0
cI=dcmplx(0.0d0,1.0d0)

open(unit=101,file='../THE.in',action='read',access='sequential',status='old')
read(101,'(A9,A200)')useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
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
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,'(A8,I3)') useless_string, wan_min
read(101,'(A8,I3)') useless_string, wan_max
close(101)

call load_super_size()
if  (nspin .eq. 1) then
        nwan_sp=2*(wan_max-wan_min+1)
else
        nwan_sp=(wan_max-wan_min+1)
endif

SOC_TF=.false.
If (nspin .eq. 1) then
        inquire( file=trim(nscf_dir)//'./make_soc/basis_data/Hsoc_R1', exist=SOC_TF )
        if (.not.SOC_TF) then
                print*, 'soc not calculated for nspin=1'
        else
                print*, 'soc file found. include in H_full'
        endif
endif


SOC_TF=.true.




H_dim=nwan*2/nspin*period_x*period_y*period_z
allocate(H_sk(H_dim,H_dim))
allocate(full_H(H_dim,H_dim))
if (SOC_TF) then
        allocate(H_soc(H_dim,H_dim))
!        allocate(H_FM(H_dim,H_dim))
endif
allocate(read_H(H_dim*nspin/2,H_dim*nspin/2))
allocate(read_AAR(H_dim*nspin/2,H_dim*nspin/2,super_nrept,3))
allocate(full_AAR(H_dim,H_dim,super_nrept,3))
allocate(full_H_sp(H_dim,H_dim,super_nrept))
allocate(H_k_super(H_dim,H_dim))
allocate(rot_U_all(H_dim,H_dim))
allocate(U_tmp(H_dim,H_dim))
allocate(full_AAR_sp(H_dim,H_dim,super_nrept,2))
allocate(Unmk_super(H_dim,H_dim,nks))
allocate(eig_super(H_dim,nks))
print*,'H_dim=',H_dim
do irept=1,super_nrept
        print*,irept
        full_H(:,:)=dcmplx(0.0d0,0.0d0)
        H_sk(:,:)=dcmplx(0.0d0,0.0d0)
        call int2char(irept,rept_char)
    !    open(unit=102,file='./super_H/hr_HYC_'//trim(rept_char),action='read',access='sequential'&
    !            &,status='old',form='unformatted')
    !    read(102)   read_H(:,:)
    !    close(102)
 

        !open(unit=102,file='./FM_hop/FM_data/HFM_'//trim(rept_char),action='read',access='stream'&
        !        &,status='old',form='unformatted')
        !read(102)   H_sk(:,:)
        !close(102)



        if (SOC_TF) then
        open(unit=102,file='./super_H/hr_SOC_'//trim(rept_char),action='read',access='sequential'&
                &,status='old',form='unformatted')
        read(102)   H_soc(:,:)
        !print*,sum(H_soc(:,:))
        close(102)
!        open(unit=102,file='./super_H/hr_FM_'//trim(rept_char),action='read',access='sequential'&
!                &,status='old',form='unformatted')
!        read(102)   H_FM(:,:)
!        close(102)

        !  print*,sum(H_soc(:,:))
        endif
        
        open(unit=102,file='./H_sk/data/Hsk_'//trim(rept_char),action='read',access='sequential'&
                &,status='old',form='unformatted')
        read(102)   H_sk(:,:)
        close(102)
        

        !open(unit=102,file='./H_sk/data/HFM_'//trim(rept_char),action='read',access='stream'&
        !        &,status='old',form='unformatted')
        !read(102)   H_sk(:,:)
        !close(102)
        



        do iwan=1,H_dim*nspin/2
        do jwan=1,H_dim*nspin/2
        do is=1,2/nspin
        do js=1,2/nspin
                i_dim=is+2/nspin*(iwan-1)
                j_dim=js+2/nspin*(jwan-1)
                if (abs(H_sk(i_dim,j_dim)) .gt.0.008 .and. abs(H_sk(i_dim,j_dim))  .lt. 0.013)then    !<0.05
 !                     print*,irept,iwan,jwan,  H_sk(i_dim,j_dim)
 !                       H_sk(i_dim,j_dim)=-1.11E-002!5.187330016993909E-002
                elseif (abs(H_sk(i_dim,j_dim)) .gt.0.0 .and. abs(H_sk(i_dim,j_dim))  .lt. 0.008)then 
               !         print*,irept,i_dim,j_dim,  H_sk(i_dim,j_dim)
 !                       H_sk(i_dim,j_dim)=0.0

                endif
              !  if (abs(H_sk(i_dim,j_dim)) .lt.0.013)then    !>0.05 
              !          H_sk(i_dim,j_dim)=0.0
              !  endif
                
        enddo
        enddo
        enddo
        enddo




        do iwan=1,H_dim*nspin/2
        do jwan=1,H_dim*nspin/2
        do is=1,2/nspin
                i_dim=is+2/nspin*(iwan-1)
                j_dim=is+2/nspin*(jwan-1)
                full_H(i_dim,j_dim)=read_H(iwan,jwan)
        enddo
        enddo
        enddo
        full_H(:,:)=0.0

       
!if (all (super_ir_vec(irept,:) .eq. (/0,0,0/))) then
!        do iwan=1,H_dim*nspin/2
!        do is=1,2/nspin
!                i_dim=is+2/nspin*(iwan-1)
!                full_H(i_dim,i_dim)=0.2*(-1)**(is)
!        enddo
!        enddo
!endif


        if (all (super_ir_vec(irept,:) .eq. (/0,0,0/))) then
                do i_super=0,period_x*period_y*period_z-1
                H_sk(9 +128*i_super:10+128*i_super, 9+128*i_super:10+128*i_super)=H_sk( 9+128*i_super:10+128*i_super, 9+128*i_super:10+128*i_super)*1.7
                H_sk(19+128*i_super:20+128*i_super,19+128*i_super:20+128*i_super)=H_sk(19+128*i_super:20+128*i_super,19+128*i_super:20+128*i_super)*1.7
                H_sk(29+128*i_super:30+128*i_super,29+128*i_super:30+128*i_super)=H_sk(29+128*i_super:30+128*i_super,29+128*i_super:30+128*i_super)*1.7
                H_sk(39+128*i_super:40+128*i_super,39+128*i_super:40+128*i_super)=H_sk(39+128*i_super:40+128*i_super,39+128*i_super:40+128*i_super)*1.7
                enddo
        endif

 
 
        full_H(:,:)=full_H(:,:)+H_sk(:,:)*13.6
        if (SOC_TF)  full_H(:,:)=full_H(:,:)+H_soc(:,:)!+H_FM(:,:)*13.6

        full_H_sp(:,:,irept)=full_H(:,:)



        open(unit=102,file='./tot_HR/hr_'//trim(rept_char),action='write',access='stream'&
                &,status='replace',form='unformatted')
        write(102)      full_H_sp(:,:,irept)
        !print*,sum(full_H_sp(:,:))
        close(102)
        !if (.not. mod(irept,10) .eq. 8) cycle
       ! print*,irept,sum(abs(H_sk(:,:,irept))**2),super_ir_vec(irept,:)

       
        if (all (super_ir_vec(irept,:) .eq. (/0,0,0/))) then
           !     do iwan=1,H_dim*nspin/2
           !             rot_H_tmp(:,:)=full_H_sp(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,irept)

           !             call diasym(rot_H_tmp,rot_eig,2,rot_U)

           !             rot_U_all(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2)=rot_U(:,:)

               ! rot_H_tmp(:,:)=H_sk(1+(iwan-1)*2:iwan*2,1+(iwan-1)*2:iwan*2,irept)

               ! Tmp22(:,:)=matmul(rot_H_tmp(:,:),rot_U(:,:))

               ! rot_H_tmp(:,:)=matmul(transpose(dconjg(rot_U(:,:))),Tmp22(:,:))
 !               print*,rot_H_tmp(:,:)

           !     enddo
           !     open(unit=102,file='./rot_U',action='write',access='append'&
           !             &,status='replace',form='unformatted')
           !     write(102)      rot_U_all(:,:)
           !     close(102)
        endif
enddo








do  iks=1,nks
        H_k_super(:,:)=cmplx(0.0d0,0.0d0)

        do irept=1,super_nrept
                H_k_super(:,:)=H_k_super(:,:)+exp(2*cI*pi*sum(super_k_list(iks,:)*super_ir_vec(irept,:)))&
                                        &*full_H_sp(:,:,irept)/super_deg(irept)
        enddo

        call  diasym(H_k_super(:,:),eig_super(:,iks),H_dim,U_tmp(:,:))

           !     Unmk_super(:,:,iks)=matmul(transpose(conjg(rot_U_all(:,:))),U_tmp(:,:))  ! iwan,ibnd
        Unmk_super(:,:,iks)=U_tmp(:,:)

       ! print*,iks,sum(abs(Unmk_super(:,:,iks))**2)               
enddo


 open(unit=101,file='./tot_HR/U_mat',action='write',access='append',status='replace',form='unformatted')
        write(101) H_dim,H_dim,nks
        write(101) Unmk_super(:,:,:)
 close(101)

 open(unit=101,file='./tot_HR/band_E',action='write',access='append',status='replace')
        do iks=1,nks
        do i_dim=1,H_dim
                write(101,'(2I8,F16.8)') iks,i_dim,eig_super(i_dim,iks)
        enddo
        enddo
 close(101)





!open(unit=102,file='./super_H/AAR',action='read',access='sequential'&
!                &,status='old',form='unformatted')
!        read(102)   read_AAR(:,:,:,:)
!close(102)


!full_AAR(:,:,:,:)=dcmplx(0.0d0,0.0d0)
!do iwan=1,H_dim*nspin/2
!do jwan=1,H_dim*nspin/2
!do is=1,2/nspin
!        i_dim=is+2/nspin*(iwan-1)
!        j_dim=is+2/nspin*(jwan-1)
!        full_AAR(i_dim,j_dim,:,:)=read_AAR(iwan,jwan,:,:)
!enddo
!enddo
!enddo

!full_AAR_sp(:,:,:,1:3)=full_AAR(:,:,:,1:3)


!do irept=1,super_nrept
!        call int2char(irept,rept_char)

!        open(unit=102,file='./tot_HR/AAR_'//trim(rept_char),action='write',access='stream'&
!                &,status='replace',form='unformatted')
!                write(102)      full_AAR_sp(:,:,irept,1:3)
!        close(102)
!enddo



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




subroutine load_super_size()
        integer :: i,j
        irept=0
        open(unit=102,file='./super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        irept=irept+1
        enddo
        close(102)
        super_nrept=irept
        allocate(super_ir_vec(super_nrept,3))
        allocate(super_deg(super_nrept))

        irept=1
        open(unit=102,file='./super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
                super_ir_vec(irept,:)=int3(:)
                irept=irept+1
        enddo
        close(102)


        iks=0
        open(unit=102,file='./Hk_commense/k_list',action='read',access='sequential',status='old')
        do
                read(102,'(3F12.8)',iostat=io) read3(:)
                if (io/=0) exit
                iks=iks+1
        enddo
        close(102)
        nks=iks
        allocate(super_k_list(nks,3))

        iks=1
        open(unit=102,file='./Hk_commense/k_list',action='read',access='sequential',status='old')
        do
                read(102,'(3F12.8)',iostat=io) read3(:)
                if (io/=0) exit
                super_k_list(iks,:)=read3(:)
                iks=iks+1
        enddo
        close(102)



        open(unit=101,file='./super_H/H_deg_HYC',action='read',access='sequential',status='old')
                do
                read(101,'(2I8)',iostat=io) i,j
                super_deg(i)=j
                if (io/=0)  exit
        enddo
        close(101)





end subroutine




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
