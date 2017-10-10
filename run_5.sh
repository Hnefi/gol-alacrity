#array=("randtrack_element_lock")
#array=("-g"  "-O2" "-O3" "-Os" '-g -fprofile-arcs -ftest-coverage' '-g -pg')
array=("./hw5src/gol")

for j in "${array[@]}"

do 
	#echo "making with $j compile option"
	#sed -i "s/^OPT_FLAGS =.*/OPT_FLAGS = -O3 -finline-limit=$j /g" Makefile
	#rm time.out
	#for I in 1  
	#do
	#	make clean > /dev/null
#		/usr/bin/time -q -a -o time.out make -j 4 &> /dev/null
		#make clean
	#done 

	echo "Running $j" 
	rm time.out
	for J in 1 2 3 4 5
	do
		#echo "/usr/bin/time -q -a -o time.out $j 4 50 &> /dev/null"	
		/usr/bin/time -q -a -o time.out $j 10000 inputs/1k.pbm outputs/1k.pbm > output 
	done 
	perl avg_times.pl 
done 

