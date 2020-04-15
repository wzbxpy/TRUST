#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "TC.cuh"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <errno.h>
#include <netdb.h>
using namespace std;
int get_cpu_id()
{
    /* Get the the current process' stat file from the proc filesystem */
    FILE* procfile = fopen("/proc/self/stat", "r");
    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof(char), to_read, procfile);
    fclose(procfile);

    // Field with index 38 (zero-based counting) is the one we want
    char* line = strtok(buffer, " ");
    for (int i = 1; i < 38; i++)
    {
        line = strtok(NULL, " ");
    }

    line = strtok(NULL, " ");
    int cpu_id = atoi(line);
    return cpu_id;
}


void check_host_name(int hostname) { //This function returns host name for
  
   if (hostname == -1) {
      perror("gethostname");
      exit(1);
   }
}
void check_host_entry(struct hostent * hostentry) { //find host info from

   if (hostentry == NULL) {
      perror("gethostbyname");
      exit(1);
   }
}
void IP_formatter(char *IPbuffer) { //convert IP string to dotted decimal
   if (NULL == IPbuffer) {
      perror("inet_ntoa");
      exit(1);
   }
}
void IP(int rank, int cpu_id)
{
    char host[256];
   char *IP;
   struct hostent *host_entry;
   int hostname;
   hostname = gethostname(host, sizeof(host)); //find the host name
   check_host_name(hostname);
   host_entry = gethostbyname(host); //find host information
   check_host_entry(host_entry);
   IP = inet_ntoa(*((struct in_addr*) host_entry->h_addr_list[0]));
   //Convert into IP string
   //printf("Current Host Name: %s, Host IP: %s, Rank: %d, CPU_ID: %d\n", host,IP,rank,cpu_id);
   
}


int main(int argc, char *argv[])
{
   int myrank,numtasks;
   int total_rank=atoi(argv[2]);
   int N_THREADS= atoi(argv[3]);    
   int N_BLOCKS=atoi(argv[4]);
   int chunk_size=atoi(argv[5]);
   int partition_num=atoi(argv[6]);
//   	assert(1>3);

   //fprintf(stderr," %s \n", argv[1]);
   //[For p2p08]:
   //int stat[16]= {1313,2451,3704,4682,6005,7118,8562,9878,11177,12506,13800,15148,16410,17716,19200,20777};
   // For Facebook combined:
   //int stat[16]={6660,12003,17400,23208,29515,34223,38163,41771,45113,48885,52745,56454,59673,65882,74390,88234};
   MPI_Status status;
   struct arguments args;
   /* Initialize the MPI library */
   MPI_Init(&argc, &argv);
   /* Determine unique id of the calling process of all processes participating
      in this MPI program. This id is usually called MPI rank. */
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

   long long global_sum,now_sum=0;
   double global_min_time,global_max_time,now_min_time=9999,now_max_time=0;
   // myrank=705;
   while (myrank<total_rank)
   {
      // cout<<myrank<<"start"<<endl;
      args=Triangle_count(myrank,argv[1],args, total_rank,N_THREADS, N_BLOCKS, chunk_size, partition_num);
      // cout<<myrank<<"end"<<endl;
      myrank+=numtasks;
      if (now_min_time>args.time) now_min_time=args.time;
      if (now_max_time<args.time) now_max_time=args.time;
      now_sum+=args.count;
   }
   int cpu_id = get_cpu_id();
   //printf("Rank %d, CPU: %d\n",myrank,cpu_id);
   //IP(myrank, cpu_id);
   //printf("%s,GPU: %d,%d, %d, %f\n",argv[1],myrank,args.edge_count,args.degree,args.time);
   //MPI_Barrier(MPI_COMM_WORLD);
   // cout<<args.time<<endl;
   MPI_Reduce(&now_sum, &global_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&now_max_time, &global_max_time, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Reduce(&now_min_time, &global_min_time, 1, MPI_DOUBLE,MPI_MIN, 0, MPI_COMM_WORLD);
   if(myrank%numtasks==0)
   {
      printf("\n%s,triangle:%lld,%f,%f,%f \n",argv[1],global_sum,global_min_time,global_max_time,(args.edge_count/global_max_time/1000000000));
   }
   MPI_Finalize();
   return 0;
}
