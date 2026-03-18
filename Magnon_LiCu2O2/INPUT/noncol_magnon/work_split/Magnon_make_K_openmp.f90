program make_Kab
        use omp_lib

        implicit none
        integer :: io, ik,nk,iwan,jwan,nwan,is,i_dim,j_dim,H_dim,nspin,i_mag,i,j,ikq,js
        integer :: jk,cal_q
        integer :: period_x,period_y,period_z,N_mag,n_super,i_super,n_site,i_site,j_site
        integer :: n1,n2,n3,n4,w1,w2,w3,w4,nwan_mag,n_mag_unit,j_mag
        integer :: m1,m2,m3,m4,i_ud,j_ud,ud1,ud2,ud3,ud4
        complex,allocatable     :: Unmk(:,:,:)
        integer,allocatable     :: mag_wan_ind(:,:),K_ind(:,:),ind(:,:)
        real,allocatable        :: band_E(:,:),k_list(:,:),w_list(:)
        real    ::read3(3),Ef,dE,delta,t1,t2,Emax,Emin,q_vec(3)
        character(len=200)      :: tot_HR_path,Hk_com_path,super_H_path,useless_string
        character(len=20)       ::q_char,irept_char
        complex,allocatable :: Uk(:,:),Ukq(:,:),Uk_rot(:,:),Ukq_rot(:,:)
        complex,allocatable :: Kab(:,:,:),Kab_local(:,:,:)
        complex,allocatable :: Kab_L(:,:,:,:,:),Kab_L_local(:,:,:,:,:)
        complex,allocatable :: Kab_R(:,:,:,:,:),Kab_R_local(:,:,:,:,:)
        complex,allocatable :: Kab_LR(:,:,:,:,:,:,:),Kab_LR_local(:,:,:,:,:,:,:)
        complex,allocatable :: sum_tmp(:),sum_tmp_R(:),sum_tmp_L(:)
        complex,allocatable :: sum_tmp_LR(:),G_tmp(:),H_wan(:,:,:)
        complex,allocatable :: qww(:,:,:,:),qww_tmp(:,:,:)
        complex ::ctmp,cI,Pauli(2,2,3)
        integer       :: iE,nE,red_M_dim,irept,nrept
        real,allocatable :: Ek(:),Ekq(:),Ek_shift(:,:),Ekq_shift(:,:)
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool,tmp_pool(2),jpool
        integer,allocatable     :: ind_rank2(:,:),nwan_per_mag_unit(:)
        integer,allocatable     :: cal_unit_wan(:,:),red_ind(:,:),ir_vec(:,:),super_deg(:)
        integer :: interpo_x,interpo_y,interpo_z,nk1,nk2,nk3,i_inter,n_interpo
        real    :: Ek_inter,Ekq_inter, c_edge,r
        complex,allocatable ::Pauli_rot(:,:,:,:),rot_local(:,:,:)
        integer :: thread_id,i_count

        cI=cmplx(0.0,1.0)
        delta=0.01
        tot_HR_path="./tot_HR"
        Hk_com_path='./Hk_commense'
        super_H_path="./super_H"
        interpo_x=1
        interpo_y=1
        interpo_z=0
        n_interpo=(1+2*interpo_x)*(1+2*interpo_y)*(1+2*interpo_z)

        Pauli(:,:,:)= cmplx(0.00,0.00)
        Pauli(1,2,1)= cmplx(1.00,0.00)
        Pauli(2,1,1)= cmplx(1.00,0.00)
        Pauli(1,2,2)=-cmplx(0.00,1.00)
        Pauli(2,1,2)= cmplx(0.00,1.00)
        Pauli(1,1,3)= cmplx(1.00,0.00)
        Pauli(2,2,3)=-cmplx(1.00,0.00)




open(unit=101,file='../../THE.in',action='read',access='sequential',status='old')
read(101,*)!'(A9,A200)')useless_string,func_dir
read(101,*)!'(A9,A200)')useless_string,nscf_dir
read(101,*)!'(A8,A200)')useless_string,THE_dir
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
close(101)


H_dim=nwan*2/nspin*period_x*period_y*period_z
n_super=period_x*period_y*period_z
open(unit=102,file='./mag.in',access='sequential',action='read',status='old')
  read(102,'(A3,F8.4)') useless_string, Ef
  read(102,'(A6,F8.4)') useless_string, delta
  read(102,'(A5,F12.6)') useless_string, Emin
  read(102,'(A5,F12.6)') useless_string, Emax
  read(102,'(A3,I4)') useless_string, nE
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


  
  
  
  call setup()
  print*,'setup done'
  n_mag=n_super*N_mag_unit




  allocate(red_ind(n_mag,maxval(nwan_per_mag_unit(:))))
  do i_super=1,n_super
        do i_mag=1,n_mag_unit
                i_site=i_mag+(i_super-1)*n_mag_unit
                do iwan=1,nwan_per_mag_unit(i_mag)
                        red_ind(i_site,iwan)=cal_unit_wan(i_mag,iwan)+(i_super-1)*nwan
                enddo
        enddo
  enddo
  red_M_dim=0
  do i_mag=1,n_mag_unit
        red_M_dim=red_M_dim+nwan_per_mag_unit(i_mag)**2*2*2
  enddo
  red_M_dim=red_M_dim*n_super


  allocate(K_ind(red_M_dim,5))
  K_ind(:,:)=0
  i_dim=1
  do i_super=1,n_super
  do i_mag=1,n_mag_unit
  do i_ud=1,2
  do iwan=1,nwan_per_mag_unit(i_mag)
  do j_ud=1,2
  do jwan=1,nwan_per_mag_unit(i_mag)
        K_ind(i_dim,:)=(/i_mag+(i_super-1)*n_mag_unit,i_ud,iwan,j_ud,jwan/)
        i_dim=i_dim+1
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  !allocate(ind(red_M_dim,2))
  !do i_dim=1,red_M_dim
  !      i_mag=K_ind(i_dim,1)
  !      n1   =K_ind(i_dim,2)
  !      n2   =K_ind(i_dim,3)
  !      ind(i_dim,:)=(/red_ind(i_mag,n1),red_ind(i_mag,n2)/)
  !enddo

  allocate(Kab   (        red_M_dim,red_M_dim,nE))
  allocate(Kab_L (2,2,    n_mag,    red_M_dim,nE))
  allocate(Kab_R (2,2,    red_M_dim,    n_mag,nE))
  allocate(Kab_LR(2,2,2,2,n_mag,    n_mag,    nE))
  allocate(sum_tmp(nE),sum_tmp_R(nE),sum_tmp_L(nE),sum_tmp_LR(nE))

  allocate(qww(nk,n_super*n_mag_unit,maxval(nwan_per_mag_unit(:)),maxval(nwan_per_mag_unit(:))))

  open(unit=103,file='../dir-qww/qww_q'//trim(q_char)//'.dat',action='read'&
                &, access='sequential', status='old',form='unformatted')
        read(103) qww
  close(103)



  Kab   (:,:,:)         =cmplx(0.0,0.0)
  Kab_L (:,:,:,:,:)     =cmplx(0.0,0.0)
  Kab_R (:,:,:,:,:)     =cmplx(0.0,0.0)
  Kab_LR(:,:,:,:,:,:,:) =cmplx(0.0,0.0)




        open(unit=101,file='report',action='write',access='append',status='replace')
                write(101,*) 'start'
                close(101)

        do irept=1,nrept
     !   print*,irept
                call int2char(irept, irept_char)
                open(unit=101,file='../tot_HR/hr_'//trim(irept_char),action='read'&
                        &,access='stream',form='unformatted')
                        read(101) H_wan(:,:,irept)
                close(101)
                open(unit=101,file='report',action='write',access='append',status='old')
                write(101,*) 'read HR irept=,',irept,'out of ',nrept, 'done'
                close(101)

        enddo


        call EXECUTE_COMMAND_LINE("mkdir -p  data")
 


        open(unit=103,file='./data/w_list',action='write',access='append',status='replace')
        do iE=1,nE
                write(103,'(E16.8)') w_list(iE)
        enddo
        close(103)


        call find_axis()
        allocate(qww_tmp(n_super*n_mag_unit,maxval(nwan_per_mag_unit(:)),maxval(nwan_per_mag_unit(:))))

        npool=H_dim**2
        allocate(pool(npool,2))
        ipool=1
        do i_dim=1,H_dim
        do j_dim=1,H_dim
                pool(ipool,:)=(/i_dim,j_dim/)
                ipool=ipool+1
        enddo
        enddo

        do ipool=npool,2,-1
                call random_number(r)
                jpool= int(r*ipool)+1
                tmp_pool(:)=pool(ipool,:)
                pool(ipool,:)=pool(jpool,:)
                pool(jpool,:) =tmp_pool(:)
        enddo



        do ik=HYC_K1,HYC_K2!1,nk


                q_vec(:)=k_list(ik,:)+k_list(cal_q,:)
                q_vec(1)=mod(q_vec(1),1.00)
                q_vec(2)=mod(q_vec(2),1.00)
                q_vec(3)=mod(q_vec(3),1.00)


                ikq=0
                do jk=1,nk
                        if (sum(abs(k_list(jk,:)-q_vec(:))**2) .lt. 10.0**(-5)) then
                                ikq=jk
                                exit
                        endif
                enddo
                if (ikq .eq. 0) then
                        print*,'ikq not found for ik=',ik
                        stop
                endif



                Ek(:)=band_E(:,ik)
                Ekq(:)=band_E(:,ikq)

                c_edge=minval(band_E(124*n_super+1:128*n_super,:))




                Uk(:,:)=Unmk(:,:,ik)
                Ukq(:,:)=Unmk(:,:,ikq)



!                do iwan=121,128
!                        print'(I8,16F16.6)',iwan,Uk(9:10,iwan),Uk(19:20,iwan),Uk(29:30,iwan),Uk(39:40,iwan)

!                enddo

!stop





                do i_super=1,n_super
               ! print*,i_super
                do i_mag=1,n_mag_unit
                do n1=1,nwan_per_mag_unit(i_mag)

                         m1=1+(red_ind(i_mag+(i_super-1)*n_mag_unit,n1)-1)*2/nspin
                         m2=2+(red_ind(i_mag+(i_super-1)*n_mag_unit,n1)-1)*2/nspin


                         Uk_rot(m1,:)=rot_local(i_mag+(i_super-1)*n_mag_unit,1,1)*Uk(m1,:)&
                                     &+rot_local(i_mag+(i_super-1)*n_mag_unit,1,2)*Uk(m2,:)
!
                         Uk_rot(m2,:)=rot_local(i_mag+(i_super-1)*n_mag_unit,2,1)*Uk(m1,:)&
                                     &+rot_local(i_mag+(i_super-1)*n_mag_unit,2,2)*Uk(m2,:)
!
                         Ukq_rot(m1,:)=rot_local(i_mag+(i_super-1)*n_mag_unit,1,1)*Ukq(m1,:)&
                                     &+rot_local(i_mag+(i_super-1)*n_mag_unit,1,2)*Ukq(m2,:)
!
                         Ukq_rot(m2,:)=rot_local(i_mag+(i_super-1)*n_mag_unit,2,1)*Ukq(m1,:)&
                                     &+rot_local(i_mag+(i_super-1)*n_mag_unit,2,2)*Ukq(m2,:)
!
                enddo
                enddo
                enddo




                qww_tmp(:,:,:)=qww(ik,:,:,:)

                call interpolate_E(ik,Ek_shift(:,:))
                call interpolate_E(ikq,Ekq_shift(:,:))


              !  print'(8F16.6)',Ek(121:128)
              !  print'(8F16.6)',Ek_shift(121:128,1)+0.256

!stop

        print*, ik,'start parallel'       
        allocate(G_tmp(nE))



        print*,SIZEOF(Kab),red_M_dim,red_M_dim,nE
        print*,SIZEOF(Kab_L),N_mag,    red_M_dim,nE
        print*,SIZEOF(Kab_R)
        print*,SIZEOF(Kab_LR)




        !$OMP PARALLEL PRIVATE(i_count,thread_id,ipool,i_dim,j_dim,G_tmp,i_inter,Ek_inter,Ekq_inter,dE, &
        !$OMP& is,js,i,j,i_mag,j_mag,i_ud,j_ud,ud1,ud2,ud3,ud4,n1,n2,n3,n4,m1,m2,m3,m4,ctmp,Kab_local,Kab_L_local,Kab_R_local,Kab_LR_local) &
        !$OMP& SHARED(Ek_shift,Ekq_shift,w_list,k_ind,red_ind,Ukq_rot,Uk_rot,Ukq,Uk,Kab,Kab_L,Kab_R,Kab_LR,qww_tmp&
        !$OMP&,red_M_dim,nE,N_mag,H_dim,interpo_x,interpo_y,interpo_z,Ef,cI,n_interpo,Delta,nspin,n_super,npool,pool)&
        !$OMP& DEFAULT(NONE)
        
        allocate(Kab_local   (red_M_dim,red_M_dim,nE))
        allocate(Kab_L_local (2,2,N_mag,    red_M_dim,nE))
        allocate(Kab_R_local (2,2,red_M_dim,    n_mag,nE))
        allocate(Kab_LR_local(2,2,2,2,n_mag,        n_mag,nE))

        Kab_local   (:,:,:)         =cmplx(0.0d0,0.0d0)
        Kab_L_local (:,:,:,:,:)     =cmplx(0.0d0,0.0d0)
        Kab_R_local (:,:,:,:,:)     =cmplx(0.0d0,0.0d0)
        Kab_LR_local(:,:,:,:,:,:,:) =cmplx(0.0d0,0.0d0)
        i_count=0
        !$OMP DO
        do ipool=1,npool
        i_dim=pool(ipool,1)
        j_dim=pool(ipool,2)
        thread_id = omp_get_thread_num()
        i_count=i_count+1
    !    if (mod(i_count,1000).eq. 1)     print *, "Thread", thread_id, "i_count", i_count
    !     if (mod(ipool,1000).eq. 1)     print *, "Thread", thread_id, "ipool", ipool,'/',npool
        !print*,i_dim

                G_tmp(:)=cmplx(0.0)
                do i_inter=1, (1+2*interpo_x)*(1+2*interpo_y)*(1+2*interpo_z)
                        Ek_inter=Ek_shift(j_dim,i_inter)
                                Ekq_inter=Ekq_shift(i_dim,i_inter)
                                if (Ekq_inter .lt. Ef .and. Ek_inter .lt. Ef) cycle
                                if (Ekq_inter .gt. Ef .and. Ek_inter .gt. Ef) cycle
                                dE=Ek_inter-Ekq_inter
                    !            if (dE .gt. 0.0 ) cycle !normal
                                if    (Ek_inter .lt. Ef .and. Ekq_inter .gt. Ef) then
                                        G_tmp(:)=G_tmp(:)-1.0/(w_list(:)+dE+cI*delta)
                                elseif(Ek_inter .gt. Ef .and. Ekq_inter .lt. Ef) then
                                        G_tmp(:)=G_tmp(:)+1.0/(w_list(:)+dE+cI*delta)
                                endif
                enddo



 !               if (Ekq_inter .lt. Ef .and. Ek_inter .lt. Ef) cycle
 !               if (Ekq_inter .gt. Ef .and. Ek_inter .gt. Ef) cycle
 !               if (abs(dE) .lt. 0.943) cycle
 !                if (abs(dE) .gt. 0.945) cycle

 !       print*,i_dim,j_dim

                G_tmp(:)=G_tmp(:)/n_interpo

                if (sum(abs(G_tmp(:))) .lt.  10.0**(-1)) cycle
!print*,i_dim,j_dim
                !do is=1,2
                !do js=1,2
                do i=1,red_M_dim
                        i_mag=k_ind(i,1)
                        ud1  =k_ind(i,2)
                        n1   =k_ind(i,3)  !1, nwan_mag
                        ud2  =k_ind(i,4)
                        n2   =k_ind(i,5)  !1, nwan_mag
                        m1=ud1+(red_ind(i_mag,n1)-1)*2/nspin
                        m2=ud2+(red_ind(i_mag,n2)-1)*2/nspin
!print*,m1,m2,m3,m4
                do j=1,red_M_dim
                        j_mag=k_ind(j,1)
                        ud3  =k_ind(j,2)
                        n3   =k_ind(j,3)
                        ud4  =k_ind(j,4)
                        n4   =k_ind(j,5)
                        m3=ud3+(red_ind(j_mag,n3)-1)*2/nspin
                        m4=ud4+(red_ind(j_mag,n4)-1)*2/nspin

  !                      if (.not. i_mag .eq. j_mag) cycle
  !                      if (.not. n1 .eq. 1) cycle
  !                      if (.not. n2 .eq. 1) cycle
  !                      if (.not. n3 .eq. 1) cycle
  !                      if (.not. n4 .eq. 1) cycle                
           !                     print*,
                        
!
!
!
                        ctmp=conjg(Ukq_rot(m1,i_dim))*Ukq_rot(m3,i_dim)*&
                            &conjg(Uk_rot(m4,j_dim))*Uk_rot(m2,j_dim)
!
                        ctmp=conjg(ctmp) !important
                        !
!
                        if (abs(ctmp) .lt. 0.00001/n_super) cycle

 !                       print'(7I8)',i_mag,i_dim,j_dim,m1,m2,m3,m4
                        !print'(8F16.8)',Ukq_rot(m1,i_dim),Uk_rot(m2,j_dim),Ukq_rot(m3,i_dim),Uk_rot(m4,j_dim)
 !                       print'(8F16.8)',Ukq(m1,i_dim),Uk(m2,j_dim),Ukq(m3,i_dim),Uk(m4,j_dim)
                        !
!                        !do iE=1,nE
                        Kab_local(i,j,:)=Kab_local(i,j,:)+ctmp*G_tmp(:)
!              !   print*,is,js,i_mag,n1,n2,qww_tmp(i_mag,n1,n2)
                        Kab_L_local(ud1,ud2,i_mag,j,:)      =Kab_L_local(ud1,ud2,i_mag,j,:) +ctmp*G_tmp(:)*qww_tmp(i_mag,n1,n2)
                        Kab_R_local(ud3,ud4,i,j_mag,:)      =Kab_R_local(ud3,ud4,i,j_mag,:) +ctmp*G_tmp(:)*conjg(qww_tmp(j_mag,n3,n4))
                        Kab_LR_local (ud1,ud2,ud3,ud4,i_mag,j_mag,:)=Kab_LR_local(ud1,ud2,ud3,ud4,i_mag,j_mag,:) +ctmp*G_tmp(:)*qww_tmp(i_mag,n1,n2)*conjg(qww_tmp(j_mag,n3,n4))
!                        !enddo 
                enddo !j
                enddo !i



                !enddo!js
                !enddo!is
       ! enddo ! i_dim
        enddo ! ipool


        !$OMP END  DO
        !$OMP CRITICAL
                Kab(:,:,:)           =Kab(:,:,:)            +Kab_local(:,:,:)
        !do iE=1,nE
                !print*,sum(Kab_local(:,:,:))
        !enddo
                Kab_L(:,:,:,:,:)     =Kab_L(:,:,:,:,:)      +Kab_L_local(:,:,:,:,:)
                Kab_R(:,:,:,:,:)     =Kab_R(:,:,:,:,:)      +Kab_R_local(:,:,:,:,:)
                Kab_LR(:,:,:,:,:,:,:)=Kab_LR(:,:,:,:,:,:,:) +Kab_LR_local(:,:,:,:,:,:,:)

        !$OMP END CRITICAL
        
        deallocate(Kab_local)
        deallocate(Kab_L_local)
        deallocate(Kab_R_local)
        deallocate(Kab_LR_local)
        !$OMP END PARALLEL
        deallocate(G_tmp)

        enddo !ik

        
        Kab   (:,:,:)=Kab   (:,:,:)/nk




        Kab_L (:,:,:,:,:)=Kab_L (:,:,:,:,:)/nk/(1.0*n_super)**0.5
        Kab_R (:,:,:,:,:)=Kab_R (:,:,:,:,:)/nk/(1.0*n_super)**0.5
        Kab_LR(:,:,:,:,:,:,:)=Kab_LR(:,:,:,:,:,:,:)/nk


!        open(unit=101,file='cmap',action='write',access='append',status='replace')
!        do i_mag=1,n_mag
!        do j_mag=1,n_mag
!                write(101,'(2I8,2E16.6)') i_mag,j_mag,Kab_LR(2,1,i_mag,j_mag,1)
!        enddo
!        enddo
!        close(101)
!print*,sum(abs(Kab_LR(:,:,:,:,:)))


        open(unit=102,file='data/Kab_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab
        close(102)






        open(unit=102,file='data/Kab_L_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab_L
        close(102)
        open(unit=102,file='data/Kab_R_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab_R
        close(102)

        open(unit=102,file='data/Kab_LR_q'//trim(q_char),action='write'&
                &,access='append',status='replace',form='unformatted')
        write(102) Kab_LR
        close(102)




contains
        subroutine setup()
                integer ::read3_int(3)
                real    :: x_min,y_min,z_min


        ik=0
        open(unit=102,file='../Hk_commense/k_list',action='read',access='sequential',status='old')
        do
                read(102,'(3F12.8)',iostat=io) read3(:)
                if (io/=0) exit
                ik=ik+1
        enddo
        close(102)
        nk=ik
        allocate(k_list(nk,3))

        x_min=999.0
        y_min=999.0
        z_min=999.0
        
        ik=1
        open(unit=102,file='../Hk_commense/k_list',action='read',access='sequential',status='old')
        do
                read(102,'(3F12.8)',iostat=io) read3(:)
                if (io/=0) exit
                k_list(ik,:)=read3(:)

                        if (.not. abs(read3(1)) .lt. 0.00001) then
                                if (x_min .gt. abs(read3(1))) x_min=abs(read3(1))
                        endif
                        if (.not. abs(read3(2)) .lt. 0.00001) then
                                if (y_min .gt. abs(read3(2))) y_min=abs(read3(2))
                        endif
                        if (.not. abs(read3(3)) .lt. 0.00001) then
                                if (z_min .gt. abs(read3(3))) z_min=abs(read3(3))
                        endif
                
                ik=ik+1
        enddo
        close(102)


        nk1=nint(1.0/x_min)
        nk2=nint(1.0/y_min)
        nk3=nint(1.0/z_min)

        print*,'nk1,nk2,nk3',nk1,nk2,nk3

        allocate(Unmk(H_dim,H_dim,nk))
        allocate(Uk(H_dim,H_dim))
        allocate(Ukq(H_dim,H_dim))
        allocate(Uk_rot(H_dim,H_dim))
        allocate(Ukq_rot(H_dim,H_dim))
        allocate(band_E(H_dim,nk))
        allocate(Ek(H_dim),Ekq(H_dim))
        allocate(Ek_shift(H_dim,n_interpo),Ekq_shift(H_dim,n_interpo))

       ! print*,'allocate done for U'

        open(unit=102,file='../tot_HR/U_mat',action='read',access='sequential',status='old',form='unformatted')
                read(102) read3_int(:)
                if (.not. read3_int(1) .eq. H_dim)  then
                        print*,'Umat dim 1 wrong'
                        stop
                endif
                if (.not. read3_int(2) .eq. H_dim)  then
                        print*,'Umat dim 2 wrong'
                        stop
                endif
                if (.not. read3_int(3) .eq. nk)  then
                        print*,'Umat dim 3 wrong', read3_int(3), 'vs', nk
                        stop
                endif
                read(102) Unmk
        close(102)

        !print*,'read Umat done'

        open(unit=102,file='../tot_HR/band_E',action='read',access='sequential',status='old')
        do
                read(102,'(2I8,F16.8)',iostat=io) read3_int(1:2), read3(1)
                if(io/=0) exit
                band_E(read3_int(2),read3_int(1))=read3(1)
        enddo
        close(102)

        !print*,'read band_E done'

        allocate(w_list(nE))

        do iE=1,nE
                w_list(iE)=Emin+(iE-1)*(Emax-Emin)/nE
        enddo


        irept=0
        open(unit=101,file='../super_H/ir_vec_HYC',action='read',access='sequential')
        do
                read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                irept=irept+1
        enddo
        close(101)
        nrept=irept
        allocate(ir_vec(nrept,3))
        allocate(super_deg(nrept))
        allocate(H_wan(H_dim,H_dim,nrept))
        irept=1
        open(unit=101,file='../super_H/ir_vec_HYC',action='read',access='sequential')
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

        open(unit=101,file='../super_H/H_deg_HYC',action='read',access='sequential')
        do
                read(101,'(2I8)',iostat=io) i,j
                super_deg(i)=j
                if (io/=0)  exit
        enddo
        close(101)




        end subroutine setup



!subroutine find_axis()
!        integer :: irept,iwan_sp
!        integer ::read3_int(3)
!        complex(kind(0d0)) :: H_R0(H_dim,H_dim)
!        complex::       M_tmp(2,2)
!        character(len=20) :: irept_char 
!        real    :: eig2(2)
!
!        allocate(local_U(N_mag,2,2))
!
!
!        irept=1
!        open(unit=101,file='./super_H/ir_vec_HYC',action='read',access='sequential')
!        do
!        read(101,'(3I8)',iostat=io) read3_int(:)
!                if(io/=0) exit
!                if (all(read3_int(:) .eq. (/0,0,0/))) then
!                        exit
!                endif
!                irept=irept+1
!        enddo
!        close(101)
!        open(unit=101,file='./H_sk/data/Hsk_'//trim(irept_char),action='read'&
!                        &,access='stream',form='unformatted')
!                        read(101) H_R0(:,:)
!        close(101)
!
!
!        do i_super=1,n_super
!        do i_mag=1,n_mag_unit
!                M_tmp(:,:)=cmplx(0.0,0.0)
!                do iwan=1,nwan_per_mag_unit(i_mag)
!                        iwan_sp=cal_unit_wan(i_mag,iwan)+(i_super-1)*nwan
!                        M_tmp(:,:)= M_tmp(:,:)+H_R0(2*iwan_sp-1:2*iwan_sp,2*iwan_sp-1:2*iwan_sp)
!                enddo
!
!                call diasym(M_tmp,eig2,2,local_U(i_mag+(i_super-1)*n_mag_unit,:,:))
!        enddo
!        enddo
!
!end subroutine find_axis




subroutine interpolate_E(ik_in,E_out)
        implicit none
        integer,intent(in) :: ik_in
        real,intent(out):: E_out(H_dim,n_interpo)
        integer :: ix,iy,iz
        real    :: k_tmp(3)
        complex :: phase,Hk_tmp(H_dim,H_dim),U_inter(H_dim,H_dim)
        complex :: I2pi

        I2pi=2*cmplx(0.0,1.0)*3.1415926535

        i_inter=1
        do ix=-interpo_x,interpo_x
        do iy=-interpo_y,interpo_y
        do iz=-interpo_z,interpo_z
 !       print*,'xyz',ix,iy,iz
                k_tmp(1)=k_list(ik_in,1)/period_x
                k_tmp(2)=k_list(ik_in,2)/period_y
                k_tmp(3)=k_list(ik_in,3)/period_z
                
                k_tmp(1)=k_tmp(1)+1.0/(1+2*interpo_x)*ix/period_x
                k_tmp(2)=k_tmp(2)+1.0/(1+2*interpo_y)*iy/period_y
                k_tmp(3)=k_tmp(3)+1.0/(1+2*interpo_z)*iz/period_z

                Hk_tmp(:,:)=cmplx(0.0,0.0)
                do irept=1,nrept
                        phase=exp(I2pi*sum(ir_vec(irept,:)*k_tmp(:)))
                        Hk_tmp(:,:)=Hk_tmp(:,:)+H_wan(:,:,irept)*phase/super_deg(irept)
                enddo
                call  diasym(Hk_tmp(:,:),E_out(:,i_inter),H_dim,U_inter(:,:))
                i_inter=i_inter+1
        enddo
        enddo
        enddo

       ! E_out(125:128,:)= E_out(125:128,:)+0.06

        E_out(:,:)=E_out(:,:)-c_edge+Ef+0.01

        E_out(1:124*n_super,:)= E_out(1:124*n_super,:)-0.00


end subroutine interpolate_E


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

subroutine find_axis()
        integer :: irept,iwan_sp,ix
        integer ::read3_int(3)
        complex(kind(0d0)) :: H_R0(H_dim,H_dim)
        complex::       M_tmp(2,2)
        character(len=20) :: irept_char
        real    :: eig2(2),st,ct,sp,cp,s1,s2,s3
        complex :: local_U(2,2),Ux(2,2)
        real    :: Mz(3,3),My(3,3),rot_M(3,3)

        allocate(Pauli_rot(N_mag,2,2,3))
        allocate(rot_local(N_mag,2,2))
        irept=1
        open(unit=101,file='../super_H/ir_vec_HYC',action='read',access='sequential')
        do
        read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                if (all(read3_int(:) .eq. (/0,0,0/))) then
                        exit
                endif
                irept=irept+1
        enddo
        close(101)

      !  print*, irept
        call int2char(irept,irept_char)
        open(unit=101,file='../H_sk/data/Hsk_'//trim(irept_char),action='read'&
                        &,access='sequential',form='unformatted')
                        read(101) H_R0(:,:)
        close(101)

        !print*,sum(abs(H_R0(:,:)))

!        M_tmp(:,:)=Pauli(:,:,3)
!        call diasym(M_tmp(:,:),eig2,2,Uz(:,:))


        do i_super=1,n_super
        do i_mag=1,n_mag_unit
                M_tmp(:,:)=cmplx(0.0,0.0)
                local_U(:,:)=cmplx(0.0,0.0)
                do iwan=1,nwan_per_mag_unit(i_mag)
                        iwan_sp=cal_unit_wan(i_mag,iwan)+(i_super-1)*nwan
                        !print*,2*iwan_sp-1,2*iwan_sp, H_dim 
                        M_tmp(:,:)= M_tmp(:,:)+H_R0(2*iwan_sp-1:2*iwan_sp,2*iwan_sp-1:2*iwan_sp)
        !                print*,i_mag,H_R0(2*iwan_sp-1:2*iwan_sp,2*iwan_sp-1:2*iwan_sp)
                enddo
                !print*,i_mag,M_tmp(:,:)

                M_tmp(:,:)=M_tmp(:,:)      !HYC check sign
                M_tmp(1,1)=(M_tmp(1,1)-M_tmp(2,2))/2.0
                M_tmp(2,2)=-M_tmp(1,1)
                M_tmp(1,2)=conjg(M_tmp(2,1))


             !   s3=M_tmp(1,1)
             !   s1=real(M_tmp(2,1))
             !   s2=aimag(M_tmp(2,1))

             !   st=((s1**2+s2**2)/(s1**2+s2**2+s3**2))**(0.5)
             !   ct=s3/(s1**2+s2**2+s3**2)**(0.5)

             !   sp=s2/(s1**2+s2**2)**(0.5)
             !   cp=s1/(s1**2+s2**2)**(0.5)


             !   My(1,:)=(/ct ,0.0, st/)
             !   My(2,:)=(/0.0,1.0,0.0/)
             !   My(3,:)=(/-st,0.0, ct/)

             !   Mz(1,:)=(/cp ,-sp,0.0/)
             !   Mz(2,:)=(/sp,  cp,0.0/)
             !   Mz(3,:)=(/0.0,0.0,1.0/)


              !  rot_M(:,:)=matmul(Mz(:,:),My(:,:))

 !               print*,sin(10.0*3.14159*2/360)

                !print*,s1,s2,s3
              !  print*,ct,st,cp,sp
                !print'(3F12.4)',rot_M(:,1)
                !print'(3F12.4)',rot_M(:,2)
                !print'(3F12.4)',rot_M(:,3)

              !  do ix=1,3
              !          Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)=&
              !                  &rot_M(1,ix)*Pauli(:,:,1)+rot_M(2,ix)*Pauli(:,:,2)+rot_M(3,ix)*Pauli(:,:,3)
              !          print'(I8,8F12.4)',i_mag,Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)
              !  enddo


                call diasym(M_tmp,eig2,2,local_U(:,:))
                rot_local(i_mag+(i_super-1)*n_mag_unit,:,:)=conjg(transpose(local_U(:,:)))   !HYC check conjg
      !           print'(I8,10F12.4)',i_mag, M_tmp(:,:)
      !           print'(I8,10F12.4)',i_mag,local_U(:,1),local_U(:,2),eig2(:)
     !           do ix=1,3
     !                   M_tmp(:,:)=cmplx(0.0,0.0)        
     !                   M_tmp(:,:)=matmul(Pauli(:,:,ix),transpose(conjg(local_U(:,:))))
     !                   Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)=matmul(local_U(:,:), M_tmp(:,:))

     !                   print'(I8,8F12.4)',i_mag,Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,ix)
     !           enddo




        enddo
        enddo



!do i_super=1,n_super
!        do i_mag=1,n_mag_unit

 !               Pauli_rot(i_mag+(i_super-1)*n_mag_unit,:,:,:)=Pauli_rot(1,:,:,:)

!enddo
!enddo

end subroutine find_axis




end program
