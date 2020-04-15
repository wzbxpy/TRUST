n=$1
# nn=`expr $n - 1`
# echo $nn
for((i=0;i< $n ;i++)) 
do
    for((j=0;j< $n ;j++)) 
    do
        echo partition$i'_'$j
        mkdir partition$i'_'$j
    done
done
./partition $n #n, total partition number=n*n
 
#1 is dataset folder
#2 soc-LiveJournal1.txt dataset file
#3 n partition number n*n parts need n*n*n worker(gpu)
