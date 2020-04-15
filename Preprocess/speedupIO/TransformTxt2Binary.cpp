#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream> 

#define bounder 100

using namespace std;
typedef 	long long 	index_t;
double times;
int vertex_count, edge_count;
int maxvertex;
typedef struct edge
{
    int u,v;
};
typedef struct edge_list
{
    int vertexID;
    vector<int> edge;
    int newid;
};
vector<edge_list> vertex;
vector<edge_list> vertexb;
vector<edge> edgelist;
bool cmp(int a,int b)
{
    return a>b;
}
bool cmp1(edge_list a, edge_list b)
{
    return a.edge.size()<b.edge.size() ;
}
bool cmp2(edge_list a, edge_list b)
{
    return a.edge.size()>b.edge.size() ;
}
void puttime(string s)
{
    times=clock()-times;
    cout<<s<<times/CLOCKS_PER_SEC<<"s"<<endl;
    times=clock();
}

void loadgraph()
{
    ifstream inFile("1.mmio", ios::in);
    if(!inFile) {
        cout << "error" <<endl;
        // return 0;
    }
    int x;
    int p=0;
    string line;
    stringstream ss;
    while(getline( inFile, line ))
	{
		if(line[0] < '0' || line[0] > '9')
			continue;
        ss.str("");
        ss.clear();
        ss << line;
        if (p==0)
        {
            ss>>vertex_count>>x>>edge_count;
            p=1;
            vertex.resize(vertex_count);
            for (int i=0;i<vertex_count;i++)
            {
                vertex[i].vertexID=i;
            }
            continue;
        }
        int u,v;
        ss>>u>>v>>x;
        u--;
        v--;
        vertex[u].edge.push_back(v);
        vertex[v].edge.push_back(u);
    }
}

void newloadgraph(string filename)
{
    FILE *fp=fopen(filename.c_str(),"r");
    edge e;
    while(fscanf(fp,"%d%d",&e.u,&e.v)!=EOF)
	{
        edgelist.push_back(e);
        maxvertex=max(maxvertex,max(e.u,e.v));
    }
    fclose(fp);
}

void selectVertex()
{
    int *a=new int[maxvertex+10];
    int *b=new int[maxvertex+10];
    int *c=new int[maxvertex+10];
    for (int i=0;i<=maxvertex;i++)
    {
        a[i]=0;
    }
    for (long long i=0;i<edgelist.size();i++)
    {
        a[edgelist[i].u]++;
        a[edgelist[i].v]++;
    }
    int k=0;
    for (int i=0;i<=maxvertex;i++)
    {
        if (a[i]) 
        {
            b[k]=a[i];
            a[i]=k;
            k++;
        }
    }
    
    ofstream graph_prop("prop.txt", ios::out);
    graph_prop<<edgelist.size()<<' '<<k;
    ofstream F_degree("degree", ios::out|ios::binary);
    F_degree.write((char*)&b[0],sizeof(int)*k);


    ofstream F_EL("BinEdgelist", ios::out|ios::binary);
    puttime("check time");
    vertex_count=k;
    vertex.resize(k);
    int *writeresult=new int[10002];
    int p=0;
    for (long long i=0;i<edgelist.size();i++)
    {
        writeresult[p]=a[edgelist[i].u];
        p++;
        writeresult[p]=a[edgelist[i].v];
        p++;
        if (p==10000)
        {
            F_EL.write((char*)&writeresult[0],sizeof(int)*p);
            p=0;
        }
        
    }
    if (p>0)
        F_EL.write((char*)&writeresult[0],sizeof(int)*p);
    free(a);
    free(b);
    free(c);
}

int main(int argc, char* argv[])
{
    times=clock();
    string Infilename="test.txt";
    string Outfilename="1.mmio";
    if (argc>1)
    {
        Infilename=argv[1];
    }
    if (argc>2)
    {
        long long edge_size=atoll(argv[2]);
        // cout<<edge_size<<endl;
        edgelist.reserve(edge_size);
        cout<<edgelist.capacity()<<endl;
    }
    if (argc>3)
    {
        Outfilename=argv[3];
    }
    newloadgraph(Infilename);

    puttime("loadgraph");
    selectVertex();
    puttime("write time");

    system("pause");
    return 0;
} 