cd $1
g++ TransformTxt2Binary.cpp -O3 -o TransformTxt2Binary
g++ preprocess_out.cpp -O3 -o preprocess
g++ partition.cpp -O3 -o partition

filename=$2
edgenum=$4
./TransformTxt2Binary $filename $4
# ./fromDirectToUndirece $filename #directgraph, edge list
./preprocess 
n=$3
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
