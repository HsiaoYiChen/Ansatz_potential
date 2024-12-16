!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM post_processing_example
  !-----------------------------------------------------------------------
  !
  ! This code is modified from the sample code provided by QE for PP.
  ! Therefore, there could contain some redundant lines which are not necessary. 
  ! This simple code
  ! 1. reads the data directory of QE, then
  ! 2. compute the exchange-correlation potential in KS-equation
  ! 3. print the potential in differnet formats as well as the k-point list  
  ! BEWARE: The Hubbard U is not included
  !
  ! Input: namelist &inputpp [outdir=...] [prefix=...] / as in QE input
  ! (default values as in QE).
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE io_files,   ONLY : tmp_dir, prefix
  USE mp_global,  ONLY : mp_startup
  USE mp_images,  ONLY : intra_image_comm
  USE mp,         ONLY : mp_bcast
  USE environment,ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  LOGICAL :: needwf = .true.
  INTEGER :: ios
  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  NAMELIST / inputpp / outdir, prefix
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  CALL environment_start ( 'PPEXAMPLE' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, intra_image_comm)
  IF ( ios /= 0) CALL errore ('postproc', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  !
  !   Read xml file, allocate and initialize general variables
  !
  CALL read_file_new ( needwf )
  !
  CALL vrs_print ( )
  !
  CALL environment_end ( 'PPEXAMPLE' )
  !
  CALL stop_pp()
  !
END PROGRAM post_processing_example
!
!-----------------------------------------------------------------------
SUBROUTINE vrs_print ( )
  !-----------------------------------------------------------------------
  !
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE noncollin_module,   ONLY: npol
  USE fft_base,   ONLY : dfftp
  USE klist,      ONLY : xk, nks, nkstot, igk_k, ngk
  USE lsda_mod,   ONLY : nspin, isk, current_spin
  USE io_files,   ONLY : restart_dir
  USE scf,        ONLY : vrs, vltot, v, kedtau
  USE gvecs,      ONLY : doublegrid
  USE uspp,       ONLY : nkb, vkb
  USE uspp_init,  ONLY : init_us_2
  USE cell_base,       ONLY :at
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.h'
  !
  INTEGER :: i,ik
  REAL(kind(0d0))  :: to_write_V_vec_dp(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nspin)
  REAL             :: to_write_V_vec(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nspin)
  REAL(DP) :: xkg(3)
  !

  CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)

  print*,'size of vrs:', size(vrs(:,1)),size(vrs(1,:))
  print*,'nnr=', dfftp%nnr
  print*,'nr1x, nr2x, nr3x,=',dfftp%nr1x, dfftp%nr2x, dfftp%nr3x
  print*,'nr1x x nr2x x nr3x,=',dfftp%nr1x* dfftp%nr2x* dfftp%nr3x


  !write vrs in the same format as QE
  open(file='./V_QE',unit=1002,action='write',access='append',status='replace'&
                                                                        &,form='unformatted')
        write(1002) dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x
        write(1002) vrs
  close(1002)
        
  !write vrs in kind(0d0) format
  to_write_V_vec_dp(:,:)=dble(vrs(1:dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,:))
  open(file='./V_vrs_dp',unit=1002,action='write',access='append',status='replace'&
                                                                        &,form='unformatted')
        write(1002) dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x
        write(1002) to_write_V_vec_dp
  close(1002)

  !write vrs in single precision
  to_write_V_vec(:,:)=vrs(1:dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,:)
  open(file='./V_vrs',unit=1002,action='write',access='append',status='replace'&
                                                                        &,form='unformatted')
        write(1002) dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x
        write(1002) to_write_V_vec
  close(1002)


  !write k in crystal coordinate for later use
  open(file='./k_list',unit=1002,action='write',access='append',status='replace')
  DO ik = 1, nkstot
     DO i = 1, 3
           xkg(i) = at(1,i)*xk(1,ik) + at(2,i)*xk(2,ik) + &
                       at(3,i)*xk(3,ik)
           ! xkg are the component in the crystal RL basis
     ENDDO
     write(1002,'(3F12.8)')  xkg(:) 
  ENDDO
  close(1002)


END SUBROUTINE vrs_print



