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

// #define dynamic
#define static
int BUCKET_SIZE = 100;
#define shared_BUCKET_SIZE 6
#define SUM_SIZE 1
#define USE_CTA 100
#define USE_WARP 2
#define without_combination 0
#define use_static 0

#define block_bucketnum 1024
#define warp_bucketnum 32

using namespace std;

__device__ void d_display(int *a, int column, int row, int start)
{
	printf("\n");
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			printf("%d\t", a[i * column + j + start]);
		}
		printf("\n");
	}
}

__device__ void kogge_sum(int *A, int *B, int len, int WARP_TID, index_t *beg_pos)
{
	/* We require enough threads for this method */
	int step = log2f(len) + 1;
	int i = WARP_TID;
	while (i < len)
	{
		A[i] = beg_pos[B[i] + 1] - beg_pos[B[i]];
		i += 32;
	}

	__syncwarp();
	for (i = 0; i < step; i++)
	{
		int pos = powf(2, i);
		int j = len - WARP_TID;
		while (j - pos >= 0)
		{
			int temp = A[j - pos];
			A[j] += temp;
			// printf("Write:%d , Read:%d , Written: %d\n",j,j-pos,A[j]);
			j -= 32;
		}
		// if(threadIdx.x==0){printf("\n\n");}
		__syncwarp();
	}
}

__device__ int linear_search(int neighbor, int *shared_partition, int *partition, int *bin_count, int bin, int BIN_START)
{

	for (;;)
	{
		int i = bin;
		int len = bin_count[i];
		int step = 0;
		int nowlen;
		if (len < shared_BUCKET_SIZE)
			nowlen = len;
		else
			nowlen = shared_BUCKET_SIZE;
		while (step < nowlen)
		{
			if (shared_partition[i] == neighbor)
			{
				return 1;
			}
			i += block_bucketnum;
			step += 1;
		}

		len -= shared_BUCKET_SIZE;
		i = bin + BIN_START;
		step = 0;
		while (step < len)
		{
			if (partition[i] == neighbor)
			{
				return 1;
			}
			i += block_bucketnum;
			step += 1;
		}
		if (len + shared_BUCKET_SIZE < 99)
			break;
		bin++;
	}
	return 0;
}

__device__ int merge(int *A, int *B, int ai, int bi, int l1_e, int l2_e, int steps)
{
	/*Reminder: As the partition is coalesced, accessing next element in each partition would require: next element --> prev + Warpsize */
	int WARPSIZE = 64;
	int count = 0;
	int steps_count = 0;
	while ((ai <= l1_e) && (bi <= l2_e))
	{
		steps_count += 1;
		// printf("\nAI: %d, value: %d \t",ai, A[ai]);
		// printf("BI: %d, value: %d \t",bi, B[bi]);
		if (A[ai] > B[bi])
		{
			bi += WARPSIZE;
		}
		else if (A[ai] < B[bi])
		{
			ai += WARPSIZE;
		}
		else
		{
			count += 1;
			ai += WARPSIZE;
			bi += WARPSIZE;
		}
		// printf("\n");
		__syncthreads();
	}
	// printf("Thread: %d, count: %d \n",threadIdx.x,count);
	return count;
}

__device__ int device_binary_search(int *arr, int value, int start, int end)
{
	int l = start, r = end;
	while (l < r - 1)
	{
		int mid = (l + r) >> 1;
		if (arr[mid] <= value)
			l = mid;
		else
			r = mid;
	}
	// if (arr[r]<=value) return r;
	if (arr[l] > value)
		return -1;
	return l;
}

int binary_search(int start, int end, int value, int *arr)
{
	// printf("low:%d,high:%d,value:%f\n",start,end,value);
	int low = start;
	int high = end;
	int index = start;
	while (low <= high)
	{
		index = ((low + high) / 2);
		if (value < arr[index])
		{
			// set high to index-1
			high = index - 1;
			// printf("high:%d\n",high);
		}
		else if (value > arr[index])
		{
			// set low to index+1
			low = index + 1;
			// printf("low:%d\n",low);
		}
		else
		{
			break;
		}
	}
	// printf("Vaue: %d,Found: %d\n",value,arr[index]);
	return index;
}

int my_binary_search(int len, int val, index_t *beg)
{
	int l = 0, r = len;
	while (l < r - 1)
	{
		int mid = (l + r) / 2;
		if (beg[mid + 1] - beg[mid] > val)
			l = mid;
		else
			r = mid;
	}
	if (beg[l + 1] - beg[l] <= val)
		return -1;
	return l;
}
__device__ int max_count(int *bin_count, int start, int end, int len)
{
	int max_count = bin_count[start];
	int min_count = bin_count[start];
	int zero_count = 0;
	for (int i = start; i < end; i++)
	{
		if (bin_count[i] > max_count)
		{
			max_count = bin_count[i];
		}
		if (bin_count[i] < min_count)
		{
			min_count = bin_count[i];
		}
		if (bin_count[i] == 0)
		{
			zero_count += 1;
		}
	}
	// printf("%d,%d,%d\n",zero_count,max_count,len);
	return max_count - 1;
}

__global__ void
dynamic_assign(vertex_t *adj_list, index_t *beg_pos, int edge_count, int vertex_count, int *partition, unsigned long long *GLOBAL_COUNT, int rank, int total_process, int BUCKET_SIZE, int T_Group, int *G_INDEX, int CHUNK_SIZE, int warpfirstvertex, unsigned long long *gettime, unsigned long long *maxcollision)
{

	// int tid=threadIdx.x+blockIdx.x*blockDim.x;
	__shared__ int bin_count[block_bucketnum];
	__shared__ int shared_partition[block_bucketnum * shared_BUCKET_SIZE + 1];
	// __shared__ int shared_now,shared_workid;
	// __shared__ int useless[1024*9];
	// useless[threadIdx.x]=1;
	unsigned long long __shared__ G_counter;
	int WARPSIZE = 32;
	if (threadIdx.x == 0)
	{
		G_counter = 0;
	}
	// timetest
	unsigned long long TT = 0, HT = 0, IT = 0;
	unsigned long long __shared__ G_TT, G_HT, G_IT;
	G_TT = 0, G_HT = 0, G_IT = 0;

	int BIN_START = blockIdx.x * block_bucketnum * BUCKET_SIZE;
	// __syncthreads();
	unsigned long long P_counter = 0;

	// unsigned long long start_time;

	// start_time = clock64();
	// CTA for large degree vertex
	int vertex = (blockIdx.x * total_process + rank) * CHUNK_SIZE;
	int vertex_end = vertex + CHUNK_SIZE;
	__shared__ int ver;
	while (vertex < warpfirstvertex)
	{

		int degree = beg_pos[vertex + 1] - beg_pos[vertex];
		// if (degree<=USE_CTA) break;
		int start = beg_pos[vertex];
		int end = beg_pos[vertex + 1];
		int now = threadIdx.x + start;
		int MODULO = block_bucketnum - 1;
		// int divide=(vert_count/blockDim.x);
		int BIN_OFFSET = 0;
		// clean bin_count
		for (int i = threadIdx.x; i < block_bucketnum; i += blockDim.x)
			bin_count[i] = 0;
		__syncthreads();

		// start_time = clock64();
		// count hash bin
		while (now < end)
		{
			int temp = adj_list[now];
			int bin = temp & MODULO;
			int index;
			for (;;)
			{
				index = atomicAdd(&bin_count[bin], 1);
				if (index < shared_BUCKET_SIZE)
				{
					shared_partition[index * block_bucketnum + bin] = temp;
					break;
				}
				else if (index < BUCKET_SIZE)
				{
					index = index - shared_BUCKET_SIZE;
					partition[index * block_bucketnum + bin + BIN_START] = temp;
					break;
				}
				break;
				index = atomicAdd(&bin_count[bin], -1);
				bin = (bin + 1) % blockDim.x;
			}
			now += blockDim.x;
		}
		__syncthreads();

		// unsigned long long hash_time=clock64()-start_time;
		// start_time = clock64();
		// list intersection
		now = beg_pos[vertex];
		end = beg_pos[vertex + 1];
		if (without_combination)
		{
			while (now < end)
			{
				int neighbor = adj_list[now];
				int neighbor_start = beg_pos[neighbor];
				int neighbor_end = beg_pos[neighbor + 1];
				int neighbor_now = neighbor_start + threadIdx.x;
				while (neighbor_now < neighbor_end)
				{
					int temp = adj_list[neighbor_now];
					int bin = temp & MODULO;
					P_counter += linear_search(temp, shared_partition, partition, bin_count, bin + BIN_OFFSET, BIN_START);
					neighbor_now += blockDim.x;
				}
				now++;
			}
		}
		else
		{
			int superwarp_ID = threadIdx.x / 64;
			int superwarp_TID = threadIdx.x % 64;
			int workid = superwarp_TID;
			now = now + superwarp_ID;
			int neighbor = adj_list[now];
			int neighbor_start = beg_pos[neighbor];
			int neighbor_degree = beg_pos[neighbor + 1] - neighbor_start;
			while (now < end)
			// while (0)
			{
				while (now < end && workid >= neighbor_degree)
				{
					now += 16;
					workid -= neighbor_degree;
					neighbor = adj_list[now];
					neighbor_start = beg_pos[neighbor];
					neighbor_degree = beg_pos[neighbor + 1] - neighbor_start;
				}
				if (now < end)
				{
					int temp = adj_list[neighbor_start + workid];
					int bin = temp & MODULO;
					P_counter += linear_search(temp, shared_partition, partition, bin_count, bin + BIN_OFFSET, BIN_START);
				}
				// __syncthreads();
				workid += 64;
			}
		}
		if (0)
		{
			int workid = threadIdx.x;
			while (now < end)
			// while (0)
			{
				int neighbor = adj_list[now];
				int neighbor_start = beg_pos[neighbor];
				int neighbor_end = beg_pos[neighbor + 1];
				while (now < end && workid - (neighbor_end - neighbor_start) >= 0)
				{
					now++;
					workid -= (neighbor_end - neighbor_start);
					neighbor = adj_list[now];
					neighbor_start = beg_pos[neighbor];
					neighbor_end = beg_pos[neighbor + 1];
				}

				// if (threadIdx.x==0)
				// {
				// 	shared_now=now;
				// 	shared_workid=workid;
				// }
				// __syncthreads();
				if (now == end)
					break;
				int temp = adj_list[neighbor_start + workid];
				int bin = temp & MODULO;
				P_counter += linear_search(temp, shared_partition, partition, bin_count, bin + BIN_OFFSET, BIN_START);
				// __syncthreads();
				workid += blockDim.x;
				// workid=shared_workid+threadIdx.x+1;
				// now=shared_now;
			}
		}

		// unsigned long long intersection_time=clock64()-start_time;
		// if (threadIdx.x==0 &&degree>3000)
		// {
		// 	int max_len_collision= max_count(bin_count,0,blockDim.x,1);
		// 	printf("%d %d %d %d %lld %lld\n",degree,vertex,blockIdx.x,max_len_collision,hash_time,intersection_time);
		// }

		__syncthreads();
		// if (vertex>1) break;
		if (use_static)
		{
			vertex += gridDim.x * total_process;
		}
		else
		{
			vertex++;
			if (vertex == vertex_end)
			{
				if (threadIdx.x == 0)
				{
					ver = atomicAdd(&G_INDEX[1], CHUNK_SIZE * total_process);
				}
				__syncthreads();
				vertex = ver;
				vertex_end = vertex + CHUNK_SIZE;
			}
		}
		// __syncthreads();
	}
	// __syncthreads();
	// unsigned long long CTA_time=clock64()-start_time;
	// start_time = clock64();

	// warp method
	int WARPID = threadIdx.x / WARPSIZE;
	int WARP_TID = threadIdx.x % WARPSIZE;
	int WARPDIM = blockDim.x * gridDim.x / WARPSIZE;
	vertex = warpfirstvertex + ((WARPID + blockIdx.x * blockDim.x / WARPSIZE) * total_process + rank) * CHUNK_SIZE;
	vertex_end = vertex + CHUNK_SIZE;
	while (vertex < vertex_count)
	// while (0)
	{
		unsigned long long start_time = clock64();
		int degree = beg_pos[vertex + 1] - beg_pos[vertex];
		if (degree < USE_WARP)
			break;
		int start = beg_pos[vertex];
		int end = beg_pos[vertex + 1];
		int now = WARP_TID + start;
		int MODULO = warp_bucketnum - 1;
		int BIN_OFFSET = WARPID * warp_bucketnum;
		// clean bin_count
		unsigned long long hash_start = clock64();

		for (int i = BIN_OFFSET + WARP_TID; i < BIN_OFFSET + warp_bucketnum; i += WARPSIZE)
			bin_count[i] = 0;
		// bin_count[threadIdx.x]=0;
		__syncwarp();

		// count hash bin
		while (now < end)
		{
			int temp = adj_list[now];
			int bin = temp & MODULO;
			bin += BIN_OFFSET;
			int index;
			for (;;)
			{
				index = atomicAdd(&bin_count[bin], 1);
				if (index < shared_BUCKET_SIZE)
				{
					shared_partition[index * block_bucketnum + bin] = temp;
					break;
				}
				else if (index < BUCKET_SIZE)
				{
					index = index - shared_BUCKET_SIZE;
					partition[index * block_bucketnum + bin + BIN_START] = temp;
					break;
				}
				break;
				index = atomicAdd(&bin_count[bin], -1);
				bin++;
				if (bin - BIN_OFFSET == 32)
					bin = BIN_OFFSET;
			}
			now += WARPSIZE;
		}
		__syncwarp();

		// unsigned long long hash_time=clock64()-hash_start;
		// unsigned long long intersection_start=clock64();
		// list intersection
		now = beg_pos[vertex];
		end = beg_pos[vertex + 1];

		if (without_combination)
		{
			while (now < end)
			{
				int neighbor = adj_list[now];
				int neighbor_start = beg_pos[neighbor];
				int neighbor_end = beg_pos[neighbor + 1];
				int neighbor_now = neighbor_start + WARP_TID;
				while (neighbor_now < neighbor_end)
				{
					int temp = adj_list[neighbor_now];
					int bin = temp & MODULO;
					P_counter += linear_search(temp, shared_partition, partition, bin_count, bin + BIN_OFFSET, BIN_START);
					neighbor_now += WARPSIZE;
				}
				now++;
			}
		}
		else
		{
			int workid = WARP_TID;
			while (now < end)
			{
				int neighbor = adj_list[now];
				int neighbor_start = beg_pos[neighbor];
				int neighbor_degree = beg_pos[neighbor + 1] - neighbor_start;
				// neighbor=__shfl_sync(0xffffffff,neighbor,31);
				// neighbor_start=__shfl_sync(0xffffffff,neighbor_start,31);
				// neighbor_degree=__shfl_sync(0xffffffff,neighbor_degree,31);

				while (now < end && workid >= neighbor_degree)
				{
					now++;
					workid -= neighbor_degree;
					neighbor = adj_list[now];
					neighbor_start = beg_pos[neighbor];
					neighbor_degree = beg_pos[neighbor + 1] - neighbor_start;
				}
				if (now < end)
				{
					int temp = adj_list[neighbor_start + workid];
					int bin = temp & MODULO;
					P_counter += linear_search(temp, shared_partition, partition, bin_count, bin + BIN_OFFSET, BIN_START);
				}
				__syncwarp();
				now = __shfl_sync(0xffffffff, now, 31);
				workid = __shfl_sync(0xffffffff, workid, 31);

				workid += WARP_TID + 1;

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
		// 	gettime[vertex]=total_time;
		// 	maxcollision[vertex]=max_count(bin_count,BIN_OFFSET,BIN_OFFSET+WARPSIZE,0);
		// }
		// if(threadIdx.x%32==0){
		// 	gettime[vertex]=1;}
		__syncwarp();
		// if (vertex>1) break;
		if (use_static)
		{
			vertex += WARPDIM * total_process;
		}
		else
		{
			vertex++;
			if (vertex == vertex_end)
			{
				if (WARP_TID == 0)
				{
					vertex = atomicAdd(&G_INDEX[2], CHUNK_SIZE * total_process);
				}
				__syncwarp();
				vertex = __shfl_sync(0xffffffff, vertex, 0);
				vertex_end = vertex + CHUNK_SIZE;
			}
		}
	}

	// unsigned long long warp_time=clock64()-start_time;

	// if (threadIdx.x==0)
	// {
	// 	printf("%d %lld %lld\n",blockIdx.x,CTA_time,warp_time);
	// }
	atomicAdd(&G_counter, P_counter);
	// atomicAdd(&G_HT,HT);
	// atomicAdd(&G_TT,TT);
	// atomicAdd(&G_IT,IT);

	__syncthreads();
	if (threadIdx.x == 0)
	{
		// printf("%d\n",G_TT);
		atomicAdd(&GLOBAL_COUNT[0], G_counter);
		// atomicAdd(&GLOBAL_COUNT[1],G_TT);
		// atomicAdd(&GLOBAL_COUNT[2],G_HT);
		// atomicAdd(&GLOBAL_COUNT[3],G_IT);
	}
}

struct arguments Triangle_count(int rank, char name[100], struct arguments args, int total_process, int n_threads, int n_blocks, int chunk_size)
{

	// fprintf(stderr,"---------------Here----------------");
	int T_Group = 32;
	int PER_BLOCK_WARP = n_threads / T_Group;
	int total = n_blocks * block_bucketnum * BUCKET_SIZE;
	unsigned long long *counter = (unsigned long long *)malloc(sizeof(unsigned long long) * 10);
	string json_file = name;
	graph *graph_d = new graph(json_file);

	// printf("Graph Adj Read: %d",graph_d->adj_list[10]);
	// int N_GPUS=argv[1];
	int deviceCount;
	HRR(cudaGetDeviceCount(&deviceCount));
	// fprintf(stderr,"----------------Device count: %d\n",deviceCount);
	//  HRR(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
	//  HRR(cudaFuncSetAttribute(dynamic_assign,cudaFuncAttributePreferredSharedMemoryCarveout,16));
	// cudaSetDevice();
	HRR(cudaSetDevice((rank + 1) % deviceCount));
	// cudaDeviceProp devProp;
	// HRR(cudaGetDeviceProperties(&devProp, rank));
	index_t vertex_count = graph_d->vert_count;
	index_t edge_count = graph_d->edge_count;
	index_t edge_list_count = graph_d->edge_list_count;
	index_t edges = graph_d->edge_count;
	/* Preprocessing Step to calculate the ratio */
	int *prefix = (int *)malloc(sizeof(int) * vertex_count);
	int temp;

	int warpfirstvertex = my_binary_search(vertex_count, USE_CTA, graph_d->beg_pos) + 1;

	// cout<<my_binary_search(vertex_count,1,graph_d->beg_pos)<<' '<<my_binary_search(vertex_count,0,graph_d->beg_pos)<<endl;
	// printf("divide position%d %d %d %d\n",graph_d->beg_pos[warpfirstvertex+1]-graph_d->beg_pos[warpfirstvertex],warpfirstvertex,graph_d->beg_pos[warpfirstvertex]-graph_d->beg_pos[warpfirstvertex-1],graph_d->beg_pos[warpfirstvertex-1]-graph_d->beg_pos[warpfirstvertex-2]);

	// fprintf(stderr,"Rank: %d, Devicecount: %d,  Start: %d, End: %d, Selected: %d\n",rank,deviceCount,E_START,E_END,(rank%deviceCount));

	// cout<<edge_count<<' '<<rank<<' '<<total_process<<endl;
	// for (int i=graph_d->beg_pos[9631];i<graph_d->beg_pos[9632];i++)
	// {
	// 	cout<<graph_d->adj_list[i]<<endl;
	// }
	// cout<<vertex_count<<' '<<edge_count<<endl;
	// cout<<graph_d->beg_pos[vertex_count-100]<<' '<<graph_d->beg_pos[vertex_count-200]<<endl;

	int *hash, *BIN_MEM;
	unsigned long long *GLOBAL_COUNT, *g_gettime, *g_maxcollision;
	int *G_INDEX;
	index_t *d_beg_pos;
	vertex_t *d_adj_list, *d_edge_list;
	float memory_req = (sizeof(int) * total + sizeof(index_t) * (vertex_count + 1) + sizeof(vertex_t) * (edge_count) + sizeof(vertex_t) * (edge_list_count)) / (1024 * 1024);
	// fprintf(stderr,"-------------------GPU: %d, Memory required: %f MB\n",rank,memory_req);
	//  printf("%f\n",memory_req);
	HRR(cudaMalloc((void **)&GLOBAL_COUNT, sizeof(unsigned long long) * 10));
	HRR(cudaMalloc((void **)&g_gettime, sizeof(unsigned long long) * (vertex_count + 1)));
	HRR(cudaMalloc((void **)&g_maxcollision, sizeof(unsigned long long) * (vertex_count + 1)));
	HRR(cudaMalloc((void **)&G_INDEX, sizeof(int) * 3));
	HRR(cudaMalloc((void **)&d_beg_pos, sizeof(index_t) * (vertex_count + 1)));
	HRR(cudaMalloc((void **)&d_adj_list, sizeof(vertex_t) * (edge_count)));
	// HRR(cudaMalloc((void **) &d_edge_list,sizeof(vertex_t)*(vertex_count+1)));
	// Swap edge list count with Eend - estart; --> gives error; may add some more

	// fprintf(stderr,">>>>>>>>>>>>>>>>>Malloc:adj_list[10]: %d\n",graph_d->adj_list[10]);

	int nowindex[3];
	nowindex[0] = chunk_size * n_blocks * n_threads / T_Group;
	nowindex[1] = chunk_size * (n_blocks * total_process + rank);
	nowindex[2] = warpfirstvertex + chunk_size * (n_blocks * n_threads / T_Group * total_process + rank);
	// unsigned long long cou=0;
	// int nowindex=0;

	HRR(cudaMemcpy(G_INDEX, &nowindex, sizeof(int) * 3, cudaMemcpyHostToDevice));
	// HRR(cudaMemcpy(GLOBAL_COUNT, &cou, sizeof(unsigned long long), cudaMemcpyHostToDevice));
	// HRR(cudaMemcpy(d_edge_list,graph_d->edge_list,sizeof(vertex_t)*(vertex_count+1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_beg_pos, graph_d->beg_pos, sizeof(index_t) * (vertex_count + 1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_adj_list, graph_d->adj_list, sizeof(vertex_t) * edge_count, cudaMemcpyHostToDevice));
	// fprintf(stderr,">>>>>>>>>>>>>>>>>>>Memcopy completed");
	// HRR(cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128*1024*1024));
	double t1 = wtime();
	double cmp_time;
	HRR(cudaMalloc((void **)&BIN_MEM, sizeof(int) * total));

	if (1)
	{
		double time_start = clock();
		// HRR(cudaMalloc((void **) &BIN_MEM,sizeof(int)*total));
		dynamic_assign<<<n_blocks, n_threads>>>(d_adj_list, d_beg_pos, edge_count, vertex_count, BIN_MEM, GLOBAL_COUNT, rank, total_process, BUCKET_SIZE, T_Group, G_INDEX, chunk_size, warpfirstvertex, g_gettime, g_maxcollision);
		HRR(cudaDeviceSynchronize());
		// HRR(cudaFree(BIN_MEM));
		cmp_time = clock() - time_start;
		cmp_time = cmp_time / CLOCKS_PER_SEC;
	}
	HRR(cudaFree(BIN_MEM));

	HRR(cudaMemcpy(counter, GLOBAL_COUNT, sizeof(unsigned long long) * 10, cudaMemcpyDeviceToHost));
	// unsigned long long *gettime= new unsigned long long[vertex_count+1];
	// HRR(cudaMemcpy(gettime,g_gettime,sizeof(unsigned long long)*(vertex_count+1), cudaMemcpyDeviceToHost));
	// unsigned long long *maxcollision= new unsigned long long[vertex_count+1];
	// HRR(cudaMemcpy(maxcollision,g_maxcollision,sizeof(unsigned long long)*(vertex_count+1), cudaMemcpyDeviceToHost));

	HRR(cudaFree(GLOBAL_COUNT));
	HRR(cudaFree(g_gettime));
	HRR(cudaFree(g_maxcollision));
	HRR(cudaFree(G_INDEX));
	HRR(cudaFree(d_beg_pos));
	HRR(cudaFree(d_adj_list));
	// HRR(cudaFree(d_edge_list));
	// free(counter);
	free(prefix);
	delete graph_d;
	args.time = cmp_time;
	args.count = counter[0];
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
	args.edge_count = edges;
	args.degree = edges / vertex_count;
	args.vertices = vertex_count;
	return args;
}
