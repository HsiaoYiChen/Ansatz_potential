for i in {1..168..4}
do
	let j=i+3
	echo ${i} ${j}
	cp -r work_split work_split_${i}_${j}
	cd work_split_${i}_${j}
	sed -i "s/HYC_K1/${i}/g" Magnon_make_K_openmp.f90 
	sed -i "s/HYC_K2/${j}/g" Magnon_make_K_openmp.f90
	sed -i "s/HYC_K1/${i}/g" sub_mag.job
        sed -i "s/HYC_K2/${j}/g" sub_mag.job
	qsub sub_mag.job
	cd ../

	sleep 60 
done
