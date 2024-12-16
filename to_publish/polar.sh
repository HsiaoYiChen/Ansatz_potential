#extract the spin polar angle from scf.out
for i in {1..128} 
do
	cd angle_${i}
	echo ${i}
	grep 'polar coord.'  scf.out  > polar
	cut -c 42-  polar  > polar2
	tail -n 12 polar2 > polar
	head -4 polar > polar2
	(readarray -t ARRAY < polar2; IFS=''; echo "${ARRAY[*]}">polar3)
	mv polar3 polar
	rm polar2
	cd ../
done
