wget http://snap.stanford.edu/data/cit-Patents.txt.gz
gzip -d cit-Patents.txt.gz
cp ../../Preprocess/fromDirectToUndirect .
cp ../../Preprocess/preprocess .
cp ../../Preprocess/partition .
./fromDirectToUndirect cit-Patents.txt
./preprocess