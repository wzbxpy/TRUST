#include <iostream>
#include "graph.h"
#include "wtime.h"
#include <queue>
#include <set>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include "herror.h"
#include <math.h>
#include "TC.cuh"
#include <assert.h>
#include <sstream> 
#include <cstring>
#include <string>
#include <fstream>
#include "comm.cuh"

// #define dynamic
#define static
int BUCKET_SIZE=100;
#define shared_BUCKET_SIZE 6
#define SUM_SIZE 1
#define USE_CTA 100
#define USE_WARP 2
#define without_combination 0
#define use_static 0

using namespace std;

__device__ 
void d_display(int *a, int column,int row,int start)
{
		printf("\n");
		for(int i=0; i<row; i++)
		{
			for(int j=0;j<column;j++)
			{
				printf("%d\t",a[i*column+j+start]);
			}
			printf("\n");
		}		  
}

__device__
void kogge_sum(int *A,int *B, int len, int WARP_TID,index_t *beg_pos)
{
    /* We require enough threads for this method */
    int step=log2f(len)+1;
    int i=WARP_TID;
	while(i<len)
	{
		A[i]=beg_pos[B[i]+1]-beg_pos[B[i]];
		i+=32;
	}
	
	__syncwarp();
    for(i=0;i<step;i++)
    {
        int pos=powf(2,i);
        int j=len-WARP_TID;
        while(j-pos>=0)
        {
			int temp=A[j-pos];
            A[j]+=temp;
            //printf("Write:%d , Read:%d , Written: %d\n",j,j-pos,A[j]);
            j-=32;
        }
        //if(threadIdx.x==0){printf("\n\n");}
		__syncwarp();
	}
}

__device__
int linear_search(int neighbor,int *shared_partition,int *partition, int *bin_count, int bin, int BIN_START)
{
	
	for(;;)
	{
		int i = bin;
		int len = bin_count[i];
		int step=0;
		int nowlen;
		if (len<shared_BUCKET_SIZE) nowlen=len; else nowlen=shared_BUCKET_SIZE;
		while(step<nowlen)
		{
			if(shared_partition[i]==neighbor)
			{
				return 1;
			}
			i+=blockDim.x;
			step+=1;
		}
		
		len-=shared_BUCKET_SIZE;
		i = bin+BIN_START;
		step=0;
		while(step<len)
		{	
			if(partition[i]==neighbor)
			{
				return 1;
			}
			i+=blockDim.x;
			step+=1;
		}
		if (len+shared_BUCKET_SIZE<99) break;
		bin++; 
	}
	return 0;
}
__device__
int merge(int *A, int *B, int ai, int bi, int l1_e, int l2_e,int steps)
{
	/*Reminder: As the partition is coalesced, accessing next element in each partition would require: next element --> prev + Warpsize */
	int WARPSIZE=64;
	int count = 0;
	int steps_count=0;
	while((ai<=l1_e) && (bi<=l2_e))
	{
		steps_count+=1;
		//printf("\nAI: %d, value: %d \t",ai, A[ai]);
		//printf("BI: %d, value: %d \t",bi, B[bi]);
		if(A[ai]>B[bi])
		{
			bi+=WARPSIZE;	
		}
		else if(A[ai]<B[bi])
		{
			ai+=WARPSIZE;
		}
		else
		{
			count+=1;
			ai+=WARPSIZE;
			bi+=WARPSIZE;
		}
		//printf("\n");
		__syncthreads();
	}
	//printf("Thread: %d, count: %d \n",threadIdx.x,count);
	return count;
}

__device__
int device_binary_search(int *arr,int value,int start,int end)
{
	int l=start,r=end;
	while (l<r-1)
	{
		int mid=(l+r)>>1;
		if (arr[mid]<=value) l=mid;
			else r=mid;
	}
	// if (arr[r]<=value) return r;
	if (arr[l]>value) return -1;
	return l;
}

int binary_search(int start,int end,int value, int *arr)
{
    //printf("low:%d,high:%d,value:%f\n",start,end,value); 
    int low=start;
    int high=end;
    int index=start;
    while (low<=high)
    {
	index=((low+high)/2);
        if (value<arr[index])
		{
            //set high to index-1
            high= index-1;
	    //printf("high:%d\n",high);
        }
        else if (value>arr[index])
        {
            // set low to index+1
            low = index+1;
            //printf("low:%d\n",low);

	}
        else
        {
            break;
        } 
	}
	//printf("Vaue: %d,Found: %d\n",value,arr[index]);
    return index;
}

int my_binary_search(int len,int val, index_t *beg)
{
	int l=0,r=len;
	while (l<r-1)
	{
		int mid=(l+r)/2;
		if (beg[mid+1]-beg[mid]>val) l=mid;
		else r=mid;
	}
	if (beg[l+1]-beg[l]<=val) return -1;
	return l;
}
__device__
int max_count(int *bin_count,int start,int end,int len)
{
	int max_count =bin_count[start];
	int min_count=bin_count[start];
	int zero_count=0;
	for (int i=start;i<end;i++)
	{
		if(bin_count[i]>max_count)
		{
			max_count=bin_count[i];			
		}
		if(bin_count[i]<min_count)
		{
			min_count=bin_count[i];
		}
		if(bin_count[i]==0)
		{
			zero_count+=1;
		}
	}
	//printf("%d,%d,%d\n",zero_count,max_count,len);
	return max_count-1;
}

// void graph_reordering(graph *graph_temp)
// {

// }

// __global__ void
// warp_hash_count(vertex_t* adj_list, index_t* beg_pos, vertex_t* edge_list, int edge_count, int vertex_count,int edge_list_count, int *partition,unsigned long long *GLOBAL_COUNT,int E_START, int E_END, int device, int BUCKETS, int G_BUCKET_SIZE, int T_Group, int *G_INDEX, int CHUNK_SIZE)
// {
// 	// Uncomment the lines below and change partition to Gpartition for using shared version
// 	int *part;
// 	int tid=threadIdx.x+blockIdx.x*blockDim.x;
// 	int WARPSIZE=T_Group;
// 	int PER_BLOCK_WARP=blockDim.x/WARPSIZE;
// 	int G_WARPID= tid/WARPSIZE;
// 	int WARPID = threadIdx.x/WARPSIZE;
// 	int WARP_TID=threadIdx.x%WARPSIZE;
// 	int __shared__ bin_count[32*4];
// 	//int __shared__ partition[160*4];
// 	int o=0, i =0;
// 	unsigned long long __shared__ G_counter;
// 	// int __shared__ warp_counter[4];
// 	// int __shared__ index[4];
// 	if (threadIdx.x==0)
// 	{
// 		G_counter=0;
// 	}
// 	unsigned long long P_counter=0;
// 	int BINsize = BUCKETS*G_BUCKET_SIZE;
// 	//int BINsize = BUCKETS*5;
// 	int BIN_START = G_WARPID*BINsize;
// 	//int BIN_START = WARPID*BINsize;
// 	int total_WARPS=gridDim.x* blockDim.x/32;
// 	int STOP = E_END - total_WARPS;
// 	//int i=G_WARPID+E_START;
// 	int RANGE= E_END-E_START;
// 	int BIN_OFFSET= WARPID*BUCKETS;

// 	int MODULO=BUCKETS-1;
// 	// int DIVIDE=(vertex_count+BUCKETS)/BUCKETS;

// 	double start_time;

// 	int SubSize=4;
// 	int SubNumber=WARPSIZE/SubSize; // Total number of sub groups in warp
// 	int Sub_Tid= WARP_TID%SubSize;	// Thread ID within a warp
// 	int SUB_ID = WARP_TID / SubSize;// ID of sub threads
// 	__shared__ int prefix_sum[SUM_SIZE*4];
// 	// #pragma unroll 5

// 	int START_INDEX= 0;
// 	//E_END=200;

// #ifdef dynamic
// 	while (*G_INDEX<E_END)
// #endif
// 	{
// #ifdef dynamic
		
// 		if(WARP_TID==0)
// 		{
// 			if (START_INDEX==0)
// 			{
// 				START_INDEX=(blockIdx.x*4+WARPID)*CHUNK_SIZE;
// 			}
// 			else
// 			{
// 				START_INDEX= atomicAdd(G_INDEX, CHUNK_SIZE);
// 			}
// 		}
// 		__syncwarp();
// 		START_INDEX= __shfl_sync(0xffffffff,START_INDEX,0);
// 		int i_end= START_INDEX + CHUNK_SIZE;
// 		if (i_end>E_END) i_end=E_END;
// #endif
		
// #ifdef static
// 		int i_end=E_END;
// 		START_INDEX=blockIdx.x*4+WARPID;
// #endif
// 		while((START_INDEX < i_end))
// 		{
			
// 			// if (WARP_TID==0)
// 			// 	printf("%d\n",START_INDEX);
// 			// N2 is for hashing and N1 is lookup
// 			//if(threadIdx.x==0){printf("I: %d,G: %d\n",i,G_WARPID);}
			
// 			int N2_start=beg_pos[START_INDEX];
// 			int N2_end= beg_pos[START_INDEX+1];	
// 			int L2= N2_end-N2_start;

// 			//-------------------- If L2 is equal to 0, continue-------------------
// 			if(L2<=1)
// 			{
// #ifdef dynamic
// 				START_INDEX+=1;
// #endif
// #ifdef static
// 				START_INDEX+=gridDim.x*PER_BLOCK_WARP;
// #endif
// 				continue;
// 			}
// 			//------------------------Clear bin counter--------------------------------
// 			//if(WARP_TID==0){printf("\n%d,%d,",START_INDEX,L2);}
// 			int id=WARP_TID+BIN_OFFSET;
// 			int end = BIN_OFFSET+BUCKETS;
// 			while(id<(end))
// 			{
// 				bin_count[id]=0;
// 				//printf("BIN: %d\n",id);
// 				id+=WARPSIZE;
// 			}
// 			__syncwarp();
// 			//--------------------------Hash source node------------------------------
// 			start_time = clock64();
// 			int start=WARP_TID + N2_start;
// 			// BIN_OFFSET is for count of number of element of each bin for all 4 warps
// 			// Hash one list 
// 			while(start<N2_end)
// 			{
// 				//if(threadIdx.x==0){printf("I: %d, Hashing: %d\n",i,L2);}
// 				int temp= adj_list[start];
// 				int bin=temp & MODULO;
// 				int index=atomicAdd(&bin_count[bin+BIN_OFFSET],1);
// 				partition[index*BUCKETS+ bin + BIN_START]=temp;
// 				//{printf("thread: %d,warp:%d, write: %d bin %d, index %d  at: %d\n",threadIdx.x,WARPID,temp,bin,index,(index*WARPSIZE+bin+BIN_START));}	
// 				start+=WARPSIZE;
// 			}
// 			__syncwarp();
// 			double hash_time=clock64();
// 			hash_time-=start_time;
// 			double prefix_sum_time;
// 			double with_combination_time;
// 			double without_combination_time;

// 			//for large degree neighbor
// 			if (0)
// 			{
// 				int N_start=N2_start;
// 				// printf("%d %d\n",edge_list[START_INDEX],START_INDEX);
// 				while (N_start<N2_start+L2)
// 				{
// 					int v=adj_list[N_start];
// 					int s=beg_pos[v]+WARP_TID;
// 					int e=beg_pos[v+1];
// 					while (s<e)
// 					{
// 						int neighbor=adj_list[s];
// 						int bin=neighbor & MODULO;
// 						P_counter+=linear_search(neighbor,partition,bin_count,bin,BIN_OFFSET,BIN_START,BUCKETS);
// 						s+=32;
// 					}
// 					N_start++;
// 				}

// 				N2_start=N2_start+L2;
// 				L2=N2_end-N2_start;
// 				if (L2<0) printf("%d\n",beg_pos[START_INDEX+1]-beg_pos[START_INDEX],edge_list[START_INDEX]);
// 			}
// 			// for short degree neighbor
// 			if (L2>0)
// 			{
// 				if (L2<SUM_SIZE)
// 				{
// 					__syncwarp();

// 					start_time = clock64();
// 					kogge_sum(&prefix_sum[WARPID*SUM_SIZE],&adj_list[N2_start],L2,WARP_TID,beg_pos);

// 					__syncwarp();
				
// 					prefix_sum_time=clock64();
// 					prefix_sum_time-=start_time;
// 					start_time = clock64();


// 					int N_start=WARP_TID;
// 					while (N_start<prefix_sum[WARPID*SUM_SIZE+L2-1])
// 					{
// 						int v=device_binary_search(prefix_sum,N_start,WARPID*SUM_SIZE,WARPID*SUM_SIZE+L2-1);
// 						int u=N_start;
// 						if (v>=0) 
// 						{
// 							u-=prefix_sum[v];
// 							v-=WARPID*SUM_SIZE;
// 						}
// 						v++;
// 						v=adj_list[N2_start+v];
// 						int neighbor=adj_list[beg_pos[v]+u];
						
// 						// printf("%d %d %d\n",u,v,neighbor);
// 						int bin=neighbor & MODULO;
// 						P_counter+=linear_search(neighbor,partition,bin_count,bin,BIN_OFFSET,BIN_START,BUCKETS);
// 						N_start+=32;
// 					}
// 					__syncwarp();
				
// 					with_combination_time=clock64();
// 					with_combination_time-=start_time;
// 					start_time = clock64();
// 				}
// 				else
// 				{
// 					int Nstart= N2_start+ SUB_ID;
// 					while (Nstart<N2_end)
// 					{
// 						int neighbor= adj_list[Nstart];
// 						////if(threadIdx.x==0){printf("Intersecting: %d\n",neighbor);}
// 						int N_start=beg_pos[neighbor];
// 						int N_end= beg_pos[neighbor+1];	
// 						int L1= N_end-N_start;
// 						//if(WARP_TID==0){printf("%d,",L1);}
// 						start=Sub_Tid + N_start;
// 						while(start<N_end)
// 						{
// 							int neighbor=adj_list[start];
// 							int bin=neighbor & MODULO;
// 							P_counter+=linear_search(neighbor,partition,bin_count,bin,BIN_OFFSET,BIN_START,BUCKETS);
// 							start+=SubSize;
// 							//printf("Tid: %d, Search:%d\n",threadIdx.x,neighbor);
// 						}
// 						Nstart+=SubNumber;
// 					}
					
// 					__syncwarp();
// 					without_combination_time=clock64();
// 					without_combination_time-=start_time;
// 				}
// 			}
// 			// if (WARP_TID==0)
// 			// 	printf("%d,%d,%d,%lf,%lf,%lf,%lf\n",START_INDEX,L2,prefix_sum[WARPID*SUM_SIZE+L2-1],hash_time,prefix_sum_time,with_combination_time,without_combination_time);
// 			// unsigned long long int stop_time= clock64();
// 			// unsigned long long int search_time = (stop_time-start_time);
// 			//if(WARP_TID==0){printf("%.d\n",search_time);}
			

// #ifdef dynamic
// 			START_INDEX+=1;
// #endif
// #ifdef static
// 			START_INDEX+=gridDim.x*PER_BLOCK_WARP;
// #endif
// 			// if (WARP_TID==0)
// 			// 	printf("startindex%d\n",START_INDEX);
// 			__syncwarp();
// 		}	
// 	}
// 	//unsigned long long int stop_time_warp=clock64();
// 	//unsigned long long int warp_time= stop_time_warp -start_time_warp;
// 	// if(WARP_TID==0)
// 	// 	{
// 	// 		printf("%d, %d\n",G_WARPID,warp_time);
// 	// 	}	
// 	atomicAdd(&G_counter,P_counter);
// 	__syncthreads();
// 	if(threadIdx.x==0)
// 	{
// 		atomicAdd(&GLOBAL_COUNT[0],G_counter);
// 		// printf("%lld\n",GLOBAL_COUNT[0]);
// 	}
	
// 	//if(tid==0){printf("Device: %d, Count:%d\n",device,GLOBAL_COUNT[0]);}
// }
   
// __global__ void
// CTA_hash_count(vertex_t* adj_list, index_t* beg_pos, vertex_t* edge_list, int edge_count, int vertex_count,int edge_list_count, int *partition,unsigned long long *GLOBAL_COUNT,int E_START, int E_END, int device, int BUCKETS, int G_BUCKET_SIZE,int T_Group)
// {
// 	// Uncomment the lines below and change partition to Gpartition for using shared version
// 	int *part;
// 	int S_BUCKET_SIZE=320;
// 	int tid=threadIdx.x+blockIdx.x*blockDim.x;
// 	int WARPSIZE=blockDim.x;
// 	int PER_BLOCK_WARP=1;
// 	int G_WARPID= blockIdx.x;
// 	int WARPID = 0;
// 	int WARP_TID=threadIdx.x;
// 	int __shared__ bin_count[256*4];
// 	//int __shared__ partition[160*4];
// 	int __shared__ G_counter;
// 	int __shared__ warp_counter[4];
// 	int __shared__ index[4];
// 	G_counter=0;
// 	int P_counter=0;
// 	int BINsize = BUCKETS*G_BUCKET_SIZE;
// 	//int BINsize = BUCKETS*5;
// 	int BIN_START = G_WARPID*BINsize;
// 	//int BIN_START = WARPID*BINsize;
	
// 	// if (WARP_TID==0)
// 	// {	
// 	// 	index[WARPID]= atomicAdd(&G_INDEX[0], 1);
// 	// }
// 	// __syncwarp();
// 	//int i=index[WARPID];
// 	//int i=0,o=0;
// 	int total_WARPS=gridDim.x* blockDim.x/32;
// 	int i=G_WARPID+E_START;
// 	int RANGE= E_END-E_START;
// 	int BIN_OFFSET= WARPID*BUCKETS;
// 	int count;
// 	int MODULO=BUCKETS-1;
// 	//TODO: Static assignment to dynamic assignment of vertices 
// 	#pragma unroll 5
// 	//if(threadIdx.x==0){printf("I: %d\n",RANGE);}
// 	clock_t start_time_warp=clock();
// 	while(i<( E_END))
// 	{
// 		// N2 is for hashing and N1 is lookup
// 		//if(threadIdx.x==0){printf("I: %d\n",i);}
// 		int N2_start=beg_pos[i];
// 		int N2_end= beg_pos[i+1];	
// 		int L2= N2_end-N2_start;
// 		//-------------------- If L2 is equal to 0, continue-------------------
// 		if(L2==0 || L2==1)
// 		{
// 			i+=gridDim.x*PER_BLOCK_WARP;
			
// 			continue;
// 			}
// 		//------------------------Clear bin counter--------------------------------
// 		int id=WARP_TID+BIN_OFFSET;
// 		int end = BIN_OFFSET+BUCKETS;
// 		while(id<(end))
// 		{
// 			bin_count[id]=0;
// 			id+=WARPSIZE;
// 		}
// 		__syncthreads();
// 		//--------------------------Hash source node------------------------------
// 		// clock_t start_time = clock();
// 		int start=WARP_TID + N2_start;
// 		// BIN_OFFSET is for count of number of element of each bin for all 4 warps
// 		// Hash one list 
// 		while(start<N2_end)
// 		{
// 			//if(threadIdx.x==0){printf("I: %d, Hashing: %d\n",i,L2);}
// 			int temp= adj_list[start];
// 			int bin=temp & MODULO;
// 			int index=atomicAdd(&bin_count[bin+BIN_OFFSET],1);
// 			partition[index*BUCKETS+ bin + BIN_START]=temp;
// 			//{printf("thread: %d,warp:%d, write: %d bin %d, index %d  at: %d\n",threadIdx.x,WARPID,temp,bin,index,(index*WARPSIZE+bin+BIN_START));}	
// 			start+=WARPSIZE;
// 		}
		
// 		__syncthreads();
// 		// clock_t stop_time = clock(); 
// 		// float hash_time = float(stop_time-start_time);
// 		// int max_len_collision= max_count(bin_count,BIN_OFFSET,BIN_OFFSET+BUCKETS,L2);
// 		//--------------------------Loop over the neighbors of the vertex--------------------------------
// 		start_time=clock();
// 		int Nstart= N2_start;
// 		while (Nstart<N2_end)
// 		{
// 			int neighbor= adj_list[Nstart];
// 			////if(threadIdx.x==0){printf("Intersecting: %d\n",neighbor);}
// 			int N_start=beg_pos[neighbor];
// 			int N_end= beg_pos[neighbor+1];	
// 			// int L1= N_end-N_start;
// 			//-----------------------------
// 			// if(L1==0)
// 			// {
// 			// 	Nstart+=1;
// 			// 	continue;
// 			// }
// 			//---------------------------
// 			start=WARP_TID + N_start;
// 			while(start<N_end)
// 			{
// 				count=0;
// 				int neighbor=adj_list[start];
// 				int bin=neighbor&MODULO;
// 				count=linear_search(neighbor,partition,bin_count,bin,BIN_OFFSET,BIN_START,BUCKETS);
// 				P_counter+=count;
// 				start+=WARPSIZE;
// 				//printf("Tid: %d, Search:%d\n",threadIdx.x,neighbor);
// 			}
// 			Nstart+=1;
			
// 			__syncthreads();
// 		}
// 		stop_time=clock();
// 		float search_time= float(stop_time-start_time);
		
// 		// if(WARP_TID==0)
// 		// {
// 		// 	printf("%d, %d, %d,%.2f, %.2f\n",i,L2,max_len_collision ,hash_time, search_time);
// 		// }
// 		i+=gridDim.x*PER_BLOCK_WARP;
// 	}
// 	clock_t stop_time_warp=clock();
// 	float warp_time= float(stop_time_warp -start_time_warp);
// 	if(WARP_TID==0)
// 		{
// 			//printf("%d, %.2f\n",G_WARPID,warp_time);
// 		}
// 	atomicAdd(&G_counter,P_counter);
// 	__syncthreads();
// 	if(threadIdx.x==0){atomicAdd(&GLOBAL_COUNT[0],G_counter);}
	
// 	//if(threadIdx.x==0){printf("Device: %d, Count:%d\n",device,GLOBAL_COUNT[0]);}
// }

__global__ void
dynamic_assign(vertex_t* adj_list_HT, index_t* beg_pos_HT, vertex_t* adj_list_intersection, index_t* beg_pos_intersection, vertex_t* adj_list_neighbor, index_t* beg_pos_neighbor, int edge_count, int vertex_count,int edge_list_count, int *partition,unsigned long long *GLOBAL_COUNT, int rank, int total_process,  int BUCKET_SIZE, int T_Group, int *G_INDEX, int CHUNK_SIZE,int warpfirstvertex,int nothreadfirstvertex, int* vertexmap, unsigned long long* gettime,unsigned long long* maxcollision)
{
	
	// printf("thread%d",threadIdx.x);
	int tid=threadIdx.x+blockIdx.x*blockDim.x;
	__shared__ int bin_count[1024];
	__shared__ int shared_partition[1024*shared_BUCKET_SIZE+1];
	unsigned long long __shared__ G_counter;
	int WARPSIZE=32;
	if (threadIdx.x==0)
	{
		G_counter=0;
	}
	//timetest
	// unsigned long long TT=0,HT=0,IT=0;
	// unsigned long long __shared__ G_TT,G_HT,G_IT;
	// G_TT=0,G_HT=0,G_IT=0;

	__syncthreads();
	unsigned long long P_counter=0;
	
	// unsigned long long start_time;
	
	// start_time = clock64();
	//CTA for large degree vertex
	int vertex=(blockIdx.x*total_process+rank)*CHUNK_SIZE;
	int vertex_end=vertex+CHUNK_SIZE;
	__shared__ int ver;
	while (vertex<warpfirstvertex)
	{
		int vertexID=vertexmap[vertex];
		int BINsize = blockDim.x*BUCKET_SIZE;
		int BIN_START = blockIdx.x*BINsize;
		// if (degree<=USE_CTA) break;
		int start=beg_pos_HT[vertexID];
		int end=beg_pos_HT[vertexID+1];
		int degree=end-start;
		int now=threadIdx.x + start;
		int MODULO=blockDim.x-1;
		// int divide=(vert_count/blockDim.x);
		int BIN_OFFSET=0;
		//clean bin_count
		bin_count[threadIdx.x]=0;
		__syncthreads();
		
		
		// start_time = clock64();
		//count hash bin
		while(now<end)
		{
			int temp= adj_list_HT[now];
			int bin=temp & MODULO;
			int index;
			for(;;)
			{
				index=atomicAdd(&bin_count[bin],1);
				if (index<shared_BUCKET_SIZE)
				{
					shared_partition[index*blockDim.x+ bin]=temp;
					break;
				}
				else if (index<BUCKET_SIZE)
				{
					index=index-shared_BUCKET_SIZE;
					partition[index*blockDim.x+ bin + BIN_START]=temp;
					break;
				}
				break;
				index=atomicAdd(&bin_count[bin],-1);
				bin=(bin+1)%blockDim.x;
			}
			now+=blockDim.x;
		}
		__syncthreads();
		
		// unsigned long long hash_time=clock64()-start_time;
		// start_time = clock64();
		//list intersection
		now=beg_pos_intersection[vertexID];
		end=beg_pos_intersection[vertexID+1];
		if (without_combination)
		{
			while (now<end)
			{
				int neighbor=adj_list_intersection[now];
				int neighbor_start=beg_pos_neighbor[neighbor];
				int neighbor_end=beg_pos_neighbor[neighbor+1];
				int neighbor_now=neighbor_start+threadIdx.x;
				while (neighbor_now<neighbor_end)
				{
					int temp=adj_list_neighbor[neighbor_now];
					int bin=temp & MODULO;
					P_counter+=linear_search(temp,shared_partition,partition,bin_count,bin+BIN_OFFSET,BIN_START);
					neighbor_now+=blockDim.x;
				}
				now++;
			}
		}
		else
		{
			int superwarp_ID=threadIdx.x/64;
			int superwarp_TID=threadIdx.x%64;
			int workid=superwarp_TID;
			now=now+superwarp_ID;
			while (now<end)
			{
				int neighbor=adj_list_intersection[now];
				int neighbor_start=beg_pos_neighbor[neighbor];
				int neighbor_end=beg_pos_neighbor[neighbor+1];
				while (now<end && workid-(neighbor_end-neighbor_start)>=0)
				{
					now+=16;
					workid-=(neighbor_end-neighbor_start);
					neighbor=adj_list_intersection[now];
					neighbor_start=beg_pos_neighbor[neighbor];
					neighbor_end=beg_pos_neighbor[neighbor+1];
				}
				if (now<end) 
				{				
					int temp=adj_list_neighbor[neighbor_start+workid];
					int bin=temp & MODULO;
					P_counter+=linear_search(temp,shared_partition,partition,bin_count,bin+BIN_OFFSET,BIN_START);
				}
				workid+=64;
			}
		}
		
		// unsigned long long intersection_time=clock64()-start_time;
		// if (threadIdx.x==0 &&degree>3000)
		// {
		// 	int max_len_collision= max_count(bin_count,0,blockDim.x,1);
		// 	printf("%d %d %d %d %lld %lld\n",degree,vertexID,blockIdx.x,max_len_collision,hash_time,intersection_time);
		// }

		__syncthreads();
		// if (vertex>1) break;
		if (use_static)
		{
			vertex+=gridDim.x*total_process;
		}
		else
		{
			vertex++;
			if (vertex==vertex_end)
			{
				if (threadIdx.x==0)
				{
					ver= atomicAdd(&G_INDEX[1], CHUNK_SIZE*total_process);
				}
				__syncthreads();
				vertex=ver;
				vertex_end=vertex+CHUNK_SIZE;
			}
		}
		__syncthreads();
	}
	__syncthreads();
	// unsigned long long CTA_time=clock64()-start_time;
	// start_time = clock64();

	//warp method
	int WARPID = threadIdx.x/WARPSIZE;
	int WARP_TID=threadIdx.x%WARPSIZE;
	int WARPDIM=blockDim.x*gridDim.x/WARPSIZE;
	vertex=warpfirstvertex+((WARPID+blockIdx.x*blockDim.x/WARPSIZE)*total_process+rank)*CHUNK_SIZE;
	vertex_end=vertex+CHUNK_SIZE;
	// vertex++;
	// while (vertex<warpfirstvertex+32768*2)
	while (vertex<nothreadfirstvertex)
	{
		int vertexID=vertexmap[vertex];
		// if (vertex==warpfirstvertex+32768*2-1 && WARP_TID==0)
		// 	printf("%d %d ok\n",vertex,vertexID);
		// unsigned long long start_time=clock64();
		int BINsize = blockDim.x*BUCKET_SIZE;
		int BIN_START = blockIdx.x*BINsize;
		int start=beg_pos_HT[vertexID];
		int end=beg_pos_HT[vertexID+1];
		int degree=end-start;
		int now=WARP_TID + start;
		int MODULO=WARPSIZE-1;
		int BIN_OFFSET=WARPID*WARPSIZE;

		// if (vertex==warpfirstvertex+32768*2-1 && WARP_TID==0)
		// 	printf("%d %d\n",start,end);
		//clean bin_count
		// unsigned long long hash_start=clock64();
		bin_count[threadIdx.x]=0;
		__syncwarp();
		
		//count hash bin
		while(now<end)
		{
			int temp= adj_list_HT[now];
			int bin=temp & MODULO;
			bin+=BIN_OFFSET;
			int index;
			for(;;)
			{
				index=atomicAdd(&bin_count[bin],1);
				if (index<shared_BUCKET_SIZE)
				{
					shared_partition[index*blockDim.x+ bin]=temp;
					break;
				}
				else if (index<BUCKET_SIZE)
				{
					index=index-shared_BUCKET_SIZE;
					partition[index*blockDim.x+ bin + BIN_START]=temp;
					break;
				}
				break;
				index=atomicAdd(&bin_count[bin],-1);
				bin++;
				if (bin-BIN_OFFSET==32) bin=BIN_OFFSET;
			}
			now+=WARPSIZE;
		}
		__syncwarp();
		
		// if (vertex==warpfirstvertex+32768*2-1 && WARP_TID==0)
		// {
		// 	printf("\n");
		// 	for (int i=0;i<32;i++)
		// 		printf("%d ",bin_count[i+BIN_OFFSET]);
		// }
		// return ;
		// unsigned long long hash_time=clock64()-hash_start;
		// unsigned long long intersection_start=clock64();
		//list intersection
		now=beg_pos_intersection[vertexID];
		end=beg_pos_intersection[vertexID+1];
		
		if (without_combination)
		{
			while (now<end)
			{
				int neighbor=adj_list_intersection[now];
				int neighbor_start=beg_pos_neighbor[neighbor];
				int neighbor_end=beg_pos_neighbor[neighbor+1];
				int neighbor_now=neighbor_start+WARP_TID;
				while (neighbor_now<neighbor_end)
				{
					int temp=adj_list_neighbor[neighbor_now];
					int bin=temp & MODULO;
					P_counter+=linear_search(temp,shared_partition,partition,bin_count,bin+BIN_OFFSET,BIN_START);
					neighbor_now+=WARPSIZE;
				}
				now++;
			}
		}
		else
		{
			int workid=WARP_TID;
			while (now<end)
			{
				int neighbor=adj_list_intersection[now];
				int neighbor_start=beg_pos_neighbor[neighbor];
				int neighbor_end=beg_pos_neighbor[neighbor+1];
				while (now<end && workid-(neighbor_end-neighbor_start)>=0)
				{
					now++;
					workid-=(neighbor_end-neighbor_start);
					neighbor=adj_list_intersection[now];
					neighbor_start=beg_pos_neighbor[neighbor];
					neighbor_end=beg_pos_neighbor[neighbor+1];
				}
				__syncwarp();
				if (now<end)
				{
					int temp=adj_list_neighbor[neighbor_start+workid];
					int bin=temp & MODULO;
					P_counter+=linear_search(temp,shared_partition,partition,bin_count,bin+BIN_OFFSET,BIN_START);
				}
				__syncwarp();
				now=__shfl_sync(0xffffffff,now,31);
				workid=__shfl_sync(0xffffffff,workid,31);
				
				workid+=WARP_TID+1;
				
				// workid+=WARPSIZE;
			}
		}
		__syncwarp();
		// unsigned long long intersection_time=clock64()-intersection_start;
		// unsigned long long total_time=clock64()-start_time;
		// if(threadIdx.x%32==0){
		// 	// printf("%d %d %d\n",total_time, hash_time, intersection_time);
		// 	// TT+=total_time;
		// 	// HT+=hash_time;
		// 	// IT+=intersection_time;
		// 	gettime[vertexID]=total_time;
		// 	maxcollision[vertex]=max_count(bin_count,BIN_OFFSET,BIN_OFFSET+WARPSIZE,0);
		// }
		// if(threadIdx.x%32==0){
		// 	gettime[vertex]=1;}
		// __syncwarp();
		// if (vertex>1) break;
		if (use_static)
		{
			vertex+=WARPDIM*total_process;
		}
		else
		{
			vertex++;
			if (vertex==vertex_end)
			{
				if (WARP_TID==0)
				{
					vertex= atomicAdd(&G_INDEX[2], CHUNK_SIZE*total_process);
				}
				__syncwarp();
				vertex= __shfl_sync(0xffffffff,vertex,0);
				vertex_end=vertex+CHUNK_SIZE;
			}
		}

		
	}

	
	
	// unsigned long long warp_time=clock64()-start_time;

	// if (threadIdx.x==0)
	// {
	// 	printf("%d %lld %lld\n",blockIdx.x,CTA_time,warp_time);
	// }
	atomicAdd(&G_counter,P_counter);
	// atomicAdd(&G_HT,HT);
	// atomicAdd(&G_TT,TT);
	// atomicAdd(&G_IT,IT);
	
	__syncthreads();
	if(threadIdx.x==0)
	{
		// printf("%d\n",G_TT);
		atomicAdd(&GLOBAL_COUNT[0],G_counter);
		// atomicAdd(&GLOBAL_COUNT[1],G_TT);
		// atomicAdd(&GLOBAL_COUNT[2],G_HT);
		// atomicAdd(&GLOBAL_COUNT[3],G_IT);
	}
}


struct arguments Triangle_count(int rank, char name[100], struct arguments args, int total_process,int n_threads , int n_blocks, int chunk_size,int partition_num)
{

	// int partition_num=2;// should be n
	int T_Group= 32;
	int PER_BLOCK_WARP= n_threads/T_Group;
	int total=n_blocks*PER_BLOCK_WARP*32*BUCKET_SIZE;
    unsigned long long *counter=(unsigned long long *)malloc(sizeof(unsigned long long)*10);
	string json_file 	= name;
	stringstream ii;
	stringstream jj;
	stringstream kk;
	ii << rank%partition_num;    
	jj << rank/partition_num%partition_num;    

	kk << rank/partition_num/partition_num%partition_num;    
	graph *graph_HT = new graph	(json_file+"/partition"+ii.str()+"_"+jj.str());
	graph *graph_intersection = new graph	(json_file+"/partition"+ii.str()+"_"+kk.str());
	graph *graph_neighbor = new graph	(json_file+"/partition"+kk.str()+"_"+jj.str());
	
	json_file=json_file+"/partition"+ii.str()+"_"+kk.str();
	string loadfile=json_file+"/division";
    fstream inFile(loadfile.c_str(), ios::in);
	int warpfirstvertex,nothreadfirstvertex;
	inFile>>warpfirstvertex>>nothreadfirstvertex;
	// cout<<warpfirstvertex<<' '<<nothreadfirstvertex<<endl;
	// warpfirstvertex=nothreadfirstvertex;
	// warpfirstvertex--;
	// nothreadfirstvertex++;
	// warpfirstvertex=0;
	inFile.close();
	loadfile=json_file+"/map";

	inFile.open(loadfile.c_str(), ios::in|ios::binary);
	int *vertexmap=new int[nothreadfirstvertex];
	inFile.read((char*)vertexmap,nothreadfirstvertex*sizeof(int));
	inFile.close();

	

	int deviceCount;
	HRR(cudaGetDeviceCount(&deviceCount));
	HRR(cudaSetDevice(rank%deviceCount));


	// float memory_req = (sizeof(int)*total + sizeof(index_t)*(vertex_count+1)+ sizeof(vertex_t)*(edge_count)+sizeof(vertex_t)*(edge_list_count))/(1024*1024);
	//fprintf(stderr,"-------------------GPU: %d, Memory required: %f MB\n",rank,memory_req);
	// printf("%f\n",memory_req);
	rank=rank/partition_num/partition_num/partition_num;
	total_process=total_process/partition_num/partition_num/partition_num;
	int *hash,* BIN_MEM;
	unsigned long long *GLOBAL_COUNT,*g_gettime, *g_maxcollision;
	int *G_INDEX, *g_vertexmap;

	// cout<<graph_HT-> vert_count<<' '<<graph_HT-> edge_count<<endl;
	index_t vertex_count=	graph_HT-> vert_count;
	index_t edge_count= graph_HT-> edge_count;
	index_t edge_list_count= graph_HT-> edge_list_count;
	index_t edges= graph_HT-> edge_count;

	index_t *d_beg_pos_HT;
	vertex_t *d_adj_list_HT;

	// cout<<vertex_count<<' '<<edge_count<<endl;
	HRR(cudaMalloc((void **) &d_beg_pos_HT,sizeof(index_t)*(vertex_count+1)));
	HRR(cudaMalloc((void **) &d_adj_list_HT,sizeof(vertex_t)*(edge_count)));

	HRR(cudaMemcpy(d_beg_pos_HT,graph_HT->beg_pos,sizeof(index_t)*(vertex_count+1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_adj_list_HT,graph_HT->adj_list,sizeof(vertex_t)*edge_count, cudaMemcpyHostToDevice));




	vertex_count=	graph_intersection-> vert_count;
	edge_count= graph_intersection-> edge_count;
	edge_list_count= graph_intersection-> edge_list_count;
	edges= graph_intersection-> edge_count;

	index_t *d_beg_pos_intersection;
	vertex_t *d_adj_list_intersection;

	HRR(cudaMalloc((void **) &d_beg_pos_intersection,sizeof(index_t)*(vertex_count+1)));
	HRR(cudaMalloc((void **) &d_adj_list_intersection,sizeof(vertex_t)*(edge_count)));

	HRR(cudaMemcpy(d_beg_pos_intersection,graph_intersection->beg_pos,sizeof(index_t)*(vertex_count+1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_adj_list_intersection,graph_intersection->adj_list,sizeof(vertex_t)*edge_count, cudaMemcpyHostToDevice));



	vertex_count=	graph_neighbor-> vert_count;
	edge_count= graph_neighbor-> edge_count;
	edge_list_count= graph_neighbor-> edge_list_count;
	edges= graph_neighbor-> edge_count;

	index_t *d_beg_pos_neighbor;
	vertex_t *d_adj_list_neighbor;

	HRR(cudaMalloc((void **) &d_beg_pos_neighbor,sizeof(index_t)*(vertex_count+1)));
	HRR(cudaMalloc((void **) &d_adj_list_neighbor,sizeof(vertex_t)*(edge_count)));

	HRR(cudaMemcpy(d_beg_pos_neighbor,graph_neighbor->beg_pos,sizeof(index_t)*(vertex_count+1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_adj_list_neighbor,graph_neighbor->adj_list,sizeof(vertex_t)*edge_count, cudaMemcpyHostToDevice));
	// cout<<"total:"<<total<<endl;

	if (1)
	{	
		HRR(cudaMalloc((void **) &GLOBAL_COUNT,sizeof(unsigned long long)*10));
		HRR(cudaMalloc((void **) &g_gettime,sizeof(unsigned long long)*(vertex_count+1)));
		HRR(cudaMalloc((void **) &g_maxcollision,sizeof(unsigned long long)*(vertex_count+1)));
		HRR(cudaMalloc((void **) &G_INDEX,sizeof(int)*3));
		HRR(cudaMalloc((void **) &BIN_MEM,sizeof(int)*total));
		HRR(cudaMalloc((void **) &g_vertexmap,sizeof(int)*(nothreadfirstvertex)));
		
		int nowindex[3];
		nowindex[0]=chunk_size*n_blocks*n_threads/T_Group;
		nowindex[1]=chunk_size*(n_blocks*total_process+rank);
		nowindex[2]=warpfirstvertex+chunk_size*(n_blocks*n_threads/T_Group*total_process+rank);
		// unsigned long long cou=0;
		// int nowindex=0;

		HRR(cudaMemcpy(G_INDEX, nowindex, sizeof(int)*3, cudaMemcpyHostToDevice));
		HRR(cudaMemcpy(g_vertexmap, vertexmap, sizeof(int)*(nothreadfirstvertex), cudaMemcpyHostToDevice));
	}


	double t1=wtime();
	double cmp_time;

	if (1)
	{
		double time_start=wtime();
		dynamic_assign<<<n_blocks,n_threads>>>(d_adj_list_HT, d_beg_pos_HT, d_adj_list_intersection, d_beg_pos_intersection, d_adj_list_neighbor, d_beg_pos_neighbor, edge_count, vertex_count,edge_list_count, BIN_MEM,GLOBAL_COUNT,rank,total_process, BUCKET_SIZE, T_Group, G_INDEX, chunk_size,warpfirstvertex,nothreadfirstvertex,g_vertexmap, g_gettime,g_maxcollision);
		HRR(cudaDeviceSynchronize());
    	cmp_time = wtime()-time_start;
	}


	HRR(cudaMemcpy(counter,GLOBAL_COUNT,sizeof(unsigned long long)*10, cudaMemcpyDeviceToHost));
	// unsigned long long *gettime= new unsigned long long[vertex_count+1];
	// HRR(cudaMemcpy(gettime,g_gettime,sizeof(unsigned long long)*(vertex_count+1), cudaMemcpyDeviceToHost));
	// unsigned long long *maxcollision= new unsigned long long[vertex_count+1];
	// HRR(cudaMemcpy(maxcollision,g_maxcollision,sizeof(unsigned long long)*(vertex_count+1), cudaMemcpyDeviceToHost));


	HRR(cudaFree(GLOBAL_COUNT));
	HRR(cudaFree(g_gettime));
	HRR(cudaFree(g_maxcollision));
	HRR(cudaFree(G_INDEX));
	HRR(cudaFree(BIN_MEM));
	HRR(cudaFree(g_vertexmap));
	HRR(cudaFree(d_beg_pos_HT));
	HRR(cudaFree(d_adj_list_HT));
	HRR(cudaFree(d_beg_pos_intersection));
	HRR(cudaFree(d_adj_list_intersection));
	HRR(cudaFree(d_beg_pos_neighbor));
	HRR(cudaFree(d_adj_list_neighbor));
	// free(counter);
	
	delete graph_HT;
	delete graph_intersection;
	delete graph_neighbor;
	delete vertexmap;
	args.time=cmp_time;
	args.count=counter[0];
	// cout<<counter[0]<<endl;
	// printf("%lld\n",args.count);
	// cout<<counter[1]<<' '<<counter[2]<<' '<<counter[3]<<endl;
	// for (int i=0;i<vertex_count;i++)
	// 	if (gettime[i]>0)
	// 		cout<<gettime[i]<<' '<<graph_d->beg_pos[i+1]-graph_d->beg_pos[i]<<' '<<maxcollision[i]+1<<endl;

	// for (int i=0;i<vertex_count;i++)
	// 	if (gettime[i]==0&&i%total_process==rank || gettime[i]==1&&i%total_process!=rank)
	// 	{
	// 		cout<<i<<endl;
	// 		break;
	// 	}
	args.edge_count=edges;
	// args.degree= SIZE;
	args.vertices= vertex_count-1;
	return args;
}    
