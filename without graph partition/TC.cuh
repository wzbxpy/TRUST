
struct arguments{
    int edge_count;
    long long count;
    double time;
    int degree;
    int vertices;
};

struct arguments Triangle_count(int rank, char input[100], struct arguments args,int total_process, int threads, int blocks , int chunk_size, int partition_num);