for i in 1 769 775 781 817 823 829 835 841  787 793 799 805 811 #97 #817 823 829 835 841 847 853 859 865  7 13 19 25 31 37  43 49 55 61 67 73 79 85 91 193 385 577 #1 769 775 781 787 793 799 805 811
do
	cp -r work_1/W.in  work_${i}
	#cp -r no_O/work_1/Magnon_make_T_openmp.f90  work_${i}
	cd work_${i}
#		rm data/Wab
	ifort   -qopenmp Magnon_make_T_openmp.f90  -O2 -mkl -o mT_op

	./mT_op
	cd ../
done
