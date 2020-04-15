for ext in .properties .graph .md5sums; 
do
	    wget -c http://data.law.di.unimi.it/webdata/$1/$1$ext
done

java -cp "../depend/*" it.unimi.dsi.webgraph.BVGraph -o -O -L $1
java -cp "../depend/*" it.unimi.dsi.webgraph.ArcListASCIIGraph  -g BVGraph $1 $1.edgelist


# $1 filename