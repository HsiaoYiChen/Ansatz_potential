for i in {1..128} 
do
	cp extract.f90 angle_${i}
	cp extract2.f90 angle_${i}
	cd angle_${i}
	ifort extract.f90  -O2
	./a.out
	ifort extract2.f90 -O2
	./a.out
	cd ../
	echo ${i}
done
