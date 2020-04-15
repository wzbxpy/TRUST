//graph.h
//Graph format: Json based format: [src_id, src_weigh,[[connected_ver_0, edge_weight],[connected_ver_1, edge_weight],[connected_ver_2, edge_weight]]]
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
#ifndef	GRAPH_H
#define	GRAPH_H

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "comm.h"

class graph{
	
	//variable
public:
	vertex_t 	vert_count;
	vertex_t	*adj_list;
	vertex_t	*head_list;
	vertex_t		*edge_list;
	index_t		*beg_pos;
	//after sort
//	vertex_t	*upperAdj;
//	vertex_t	*upperHead;
	
	Edge		*OrientedEdge;

	index_t	*upperBegin;
	index_t	*upperDegree;
	index_t	upperEdgeCount;
	
	index_t		edge_count;
	index_t 	edge_list_count;
	//after partition
	vertex_t**	partAdj;
	vertex_t**	partHead;
	index_t**	partBegin;
	index_t*	partEdgeCount;
	
	index_t		*count;	
	int 		*valid;
//dynaic scheduling
	index_t		ChunkNum;	// = roof(upperEdgeCount/BufferSize), the number of chunks for workload edge list

	index_t		*ds_count;	//ds_count[P][i] is the count of partition P chunk i	
	index_t		*ds_complete;	//ds_progress[D]: the ID of first chunk need worker
	index_t		*ds_status;	// how many work remain for [P][C]
	index_t		*ds_help;	//ds_last[P]; initiate with the last chunk it has to work with
	index_t		*ds_last;	//ds_last[P]; initiate with the last chunk it has to work with

	//gpu data
	GPU_data *gdata;

	//constructor
public:
	graph() {};
	graph(	std::string filename);//,
	~graph();
	void preproc();
	void reduceResult();

	void initDevice(int GPU_id,int Part_id);
	void DeviceCompute(int GPU_id,index_t Chunk_id);
	void gpuReduce(int GPU_id);
	void gpuProc(int GPU_id);

	void cpuProc();
	void cpuCompute(int Part_id, index_t Chunk_id);

};

//#include "graph.c"
#include "graph.hpp"
//#include "kernel.cu"
#endif
