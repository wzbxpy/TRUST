//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


#include "graph.h"
#include <fstream>
#include <omp.h>

#define FILE_NOT_EXIST	1
#define FILE_EXIST	0

#define CPU_id GPU_NUM
using namespace std;

graph::graph(
	string jsonfile)//,
{
//	cout<<"read from folder "<<jsonfile<<endl;
	
	string s_begin = jsonfile+"/begin.bin";
	string s_adj = jsonfile+"/adjacent.bin";
	string s_degree = jsonfile+"/edge";

	char* begin_file = const_cast<char*>(s_begin.c_str());
	char* adj_file = const_cast<char*>(s_adj.c_str());
	char* degree_file = const_cast<char*>(s_degree.c_str());

	vert_count = fsize(begin_file)/sizeof(index_t) - 1;
	edge_count = fsize(adj_file)/sizeof(vertex_t);

	// cout<<"vert:"<< vert_count<<"  edge: "<<edge_count<<endl;


	FILE *pFile= fopen(adj_file,"rb");
	adj_list = (vertex_t *)malloc(fsize(adj_file));
//	printf("adj_file size: %d\n", fsize(adj_file));
	fread(adj_list,sizeof(vertex_t),edge_count,pFile);
	fclose(pFile);


	FILE *pFile3 = fopen(begin_file,"rb");
	beg_pos = (index_t *)malloc(fsize(begin_file));
	fread(beg_pos,sizeof(index_t),vert_count+1,pFile3);
	fclose(pFile3);

// 	count = new index_t[256];
// //	valid = (int *)malloc(vert_count*sizeof(int));
// 	gdata = new GPU_data[GPU_NUM];
// 	for(int i=0; i<GPU_NUM; i++){
// 		gdata[i].id = i;
// 		gdata[i].EdgeBuffer = new Edge* [2];
// 	}

}
graph::~graph()
{
	free(adj_list);
	// free(edge_list);
	free(beg_pos);
}

