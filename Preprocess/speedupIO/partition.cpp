#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream> 
#include <cstring>
#include <string>
#define bounder 100

using namespace std;
typedef 	long long 	index_t;
typedef 	int 	vertex_t;
long long times;
long long vertex_count, edge_count;
typedef struct edge_list
{
    int vertexID;
    vector<int> edge;
    int initaldegree;
};
vertex_t * adj_list;
index_t * beg_pos;
vector<edge_list> vertex;
int endofprocess=-1;
bool cmp1(edge_list a, edge_list b)
{
    return a.edge.size()<b.edge.size() ;
}
bool cmp2(edge_list a, edge_list b)
{
    return a.edge.size()>b.edge.size() ;
}
long long file_size(char* filename)  
{  
    FILE *fp=fopen(filename,"r");  
    if(!fp) return -1;  
    fseek(fp,0L,SEEK_END);  
    long long size=ftell(fp);  
    cout<<size<<endl;
    fclose(fp);  
      
    return size;  
}
void loadgraph()
{
    string s_begin = "begin.bin";
	string s_adj = "adjacent.bin";

	char* begin_file = const_cast<char*>(s_begin.c_str());
	char* adj_file = const_cast<char*>(s_adj.c_str());
    cout<<1<<endl;

	vertex_count = file_size(begin_file)/sizeof(index_t) - 1;
	edge_count = file_size(adj_file)/sizeof(vertex_t);

	cout<<"vert:"<< vertex_count<<"  edge: "<<edge_count<<endl;


	FILE *pFile= fopen(adj_file,"rb");
	adj_list = (vertex_t *)malloc(file_size(adj_file));
	// printf("adj_file size: %d\n", file_size(adj_file));
	fread(adj_list,sizeof(vertex_t),edge_count,pFile);
	fclose(pFile);


	FILE *pFile3 = fopen(begin_file,"rb");
	beg_pos = (index_t *)malloc(file_size(begin_file));
	fread(beg_pos,sizeof(index_t),vertex_count+1,pFile3);
	fclose(pFile3);
}
void partition(int k)
{
    string s_begin[100][100];
	string s_adj[100][100];
    ofstream beginFile[100];
    ofstream adjFile[100];
    index_t zero=0;
    for (int i=0;i<k;i++)
    {
        for (int j=0;j<k;j++)
        {
            stringstream ii;
            ii << i;    
            stringstream jj;
            jj << j;    
            s_begin[i][j]="partition"+ii.str()+"_"+jj.str()+"/begin.bin";
            s_adj[i][j]="partition"+ii.str()+"_"+jj.str()+"/adjacent.bin";
            // cout<<s_begin[i][j]<<' '<<s_adj[i][j]<<endl;
            
            beginFile[i*k+j].open(s_begin[i][j].c_str(), ios::out|ios::binary);
            adjFile[i*k+j].open(s_adj[i][j].c_str(), ios::out|ios::binary);
            beginFile[i*k+j].write((char*)&zero,sizeof(index_t));
        }
    }
    index_t beg_count[100][100];

    for (int i=0;i<k;i++)
        for (int j=0;j<k;j++)
            beg_count[i][j]=0;
    for (int i=0;i<vertex_count;i++)
    {
        index_t count[100];
        for (int j=0;j<k;j++)
            count[j]=0;
        int ii=i%k;
        if (endofprocess==-1 && beg_pos[i+1]-beg_pos[i]<2) endofprocess=i;
        for (index_t j=beg_pos[i];j<beg_pos[i+1];j++)
        {
            int jj=adj_list[j]%k;
            count[jj]++;
            int res=adj_list[j]/k;
            adjFile[ii*k+jj].write((char*)&res,sizeof(vertex_t));
        }
        for (int j=0;j<k;j++)
        {
            beg_count[ii][j]+=count[j];
            beginFile[ii*k+j].write((char*)&beg_count[ii][j],sizeof(index_t));
        }
    }
}
typedef struct vertexwithdegree
{
    int vertexid;
    int degree;
};
bool cmpfordegree(vertexwithdegree a, vertexwithdegree b)
{
    return a.degree>b.degree;
}
int mapvertex(string foldname,int k)
{
    string s_begin = foldname+"/begin.bin";
	char* begin_file = const_cast<char*>(s_begin.c_str());
	FILE *pFile3 = fopen(begin_file,"rb");
	fread(beg_pos,sizeof(index_t),endofprocess/k+2,pFile3);
    vertexwithdegree *ver=new vertexwithdegree[endofprocess/k+1];
    for (int i=0;i<endofprocess/k+1;i++)
    {
        ver[i].vertexid=i;
        ver[i].degree=beg_pos[i+1]-beg_pos[i];
    }
    sort(ver,ver+endofprocess/k+1,cmpfordegree);
    s_begin = foldname+"/map";
    fstream outFile(s_begin.c_str(), ios::out|ios::binary);
    int firstwarpvertex=-1,firstnothreadvertex=-1;
    int sum=0;
    // cout<<"assd"<<endl;
    for (int i=0;i<endofprocess/k+1;i++)
    {
        // cout<<ver[i].degree<<endl;
        if (firstwarpvertex==-1 && ver[i].degree<=bounder) firstwarpvertex=i;
        if (firstnothreadvertex==-1 && ver[i].degree<1) {firstnothreadvertex=i;break;}
        outFile.write((char*)&ver[i].vertexid,sizeof(int));
        // cout<<ver[i].vertexid<<endl;
        sum++;
    }
    if (firstnothreadvertex==-1) firstnothreadvertex=endofprocess/k;
    cout<<sum<<' '<<firstnothreadvertex<<endl;
    s_begin = foldname+"/division";
    outFile.close();
    outFile.open(s_begin.c_str(), ios::out);
    outFile<<firstwarpvertex<<' '<<firstnothreadvertex<<endl;

}
void puttime(string s)
{
    times=clock()-times;
    cout<<s<<times/CLOCKS_PER_SEC<<"s"<<endl;
    times=clock();
}
int main(int argc,char *argv[])
{
    times=clock();
    int k=atoi(argv[1]);
    loadgraph();
    puttime("load");
    partition(k); 
    puttime("partition");   
    // k--;
    for (int i=0;i<k;i++)
    {
        for (int j=0;j<k;j++)
        {            
            stringstream ii;
            ii << i;    
            stringstream jj;
            jj << j;    
            mapvertex("partition"+ii.str()+"_"+jj.str(),k);
        }
    }
    
    puttime("map");   
    system("pause");
    return 0;
} 