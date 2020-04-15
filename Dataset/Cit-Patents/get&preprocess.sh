wget http://snap.stanford.edu/data/cit-Patents.txt.gz
gzip -d cit-Patents.txt.gz
cp ../../Preprocess/fromDirectToUndirece .
cp ../../Preprocess/preprocess .
cp ../../Preprocess/partition .
./fromDirectToUndirece cit-Patents.txt
./preprocess