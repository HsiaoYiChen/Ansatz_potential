program make_angle

        !! This code makes a bash script to create 
        !! folder and modify the scf.in put for a 
        !! random spin orientation set.

        !! In this exmaple code, there are 4 magnetic atoms.


        implicit none
        integer :: i,j
        character(len=20) :: i_char
        real(kind(0d0))    :: to_write(8)
        real(kind(0d0))    :: angle8(8)
        real(kind(0d0)),allocatable :: angle_list(:,:) 
        integer :: n_file,restart_i

        call RANDOM_SEED()

        restart_i=0
        n_file=128
        allocate(angle_list(n_file,8))
        
        do i=1,n_file,2
                
                print*, i
        
                CALL RANDOM_NUMBER(angle8(:))
        
                do j=1,4
                        angle_list(i,j*2-1)=angle8(j*2-1)*180
                        angle_list(i,j*2  )=angle8(j*2  )*360
                enddo
                !! Here we don't make the set fully random
                !! n and n+1 initial guess z-component is inversed
                do j=1,4
                        angle_list(i+1,j*2-1)=180-angle8(j*2-1)*180
                        angle_list(i+1,j*2  )=    angle8(j*2  )*360
                enddo
        enddo



        open(unit=101,file='run.sh',access='append',action='write',status='replace')
        do i=1,n_file
                to_write(:)=angle_list(i,:)
                call  int2char(i+restart_i,i_char)

                write(101,*) 'mkdir angle_'//trim(i_char)
                write(101,*) 'cd angle_'//trim(i_char)

                write(101,*) 'cp ../IN_files/* ./'

                write(101,9003) i
                write(101,9007) to_write(1)
                write(101,9008) to_write(2)
                write(101,9009) to_write(3)
                write(101,9010) to_write(4)
                write(101,9011) to_write(5)
                write(101,9012) to_write(6)
                write(101,9013) to_write(7)
                write(101,9014) to_write(8)
                write(101,*)"cd ../"
        enddo
        close(101)



9003 FORMAT('sed -i "s/JObname/',I0.3,'/g" sub_scf.job')

9007 FORMAT('sed -i "s/Ang1_1/',F0.6,'/g" scf.in')
9008 FORMAT('sed -i "s/Ang1_2/',F0.6,'/g" scf.in')
9009 FORMAT('sed -i "s/Ang2_1/',F0.6,'/g" scf.in')
9010 FORMAT('sed -i "s/Ang2_2/',F0.6,'/g" scf.in')

9011 FORMAT('sed -i "s/Ang3_1/',F0.6,'/g" scf.in')
9012 FORMAT('sed -i "s/Ang3_2/',F0.6,'/g" scf.in')
9013 FORMAT('sed -i "s/Ang4_1/',F0.6,'/g" scf.in')
9014 FORMAT('sed -i "s/Ang4_2/',F0.6,'/g" scf.in')

9015 FORMAT('sed -i "s/Ang5_1/',F0.6,'/g" scf.in')
9016 FORMAT('sed -i "s/Ang5_2/',F0.6,'/g" scf.in')
9017 FORMAT('sed -i "s/Ang6_1/',F0.6,'/g" scf.in')
9018 FORMAT('sed -i "s/Ang6_2/',F0.6,'/g" scf.in')

9019 FORMAT('sed -i "s/Ang7_1/',F0.6,'/g" scf.in')
9020 FORMAT('sed -i "s/Ang7_2/',F0.6,'/g" scf.in')
9021 FORMAT('sed -i "s/Ang8_1/',F0.6,'/g" scf.in')
9022 FORMAT('sed -i "s/Ang8_2/',F0.6,'/g" scf.in')




contains 

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )      WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )   WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )   WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )  WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char


end program
