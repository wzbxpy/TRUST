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
    for (int i=0;i<=maxvertex;i++)
    {
        a[i]=0;
    }
    for (int i=0;i<edgelist.size();i++)
    {
        a[edgelist[i].u]=1;
        a[edgelist[i].v]=1;
    }
    int k=0;
    for (int i=0;i<=maxvertex;i++)
    {
        if (a[i]) 
        {
            a[i]=k;
            k++;
        }
    }
    puttime("check time");
    vertex_count=k;
    vertex.resize(k);
    for (int i=0;i<edgelist.size();i++)
    {
        int u=a[edgelist[i].u];
        int v=a[edgelist[i].v];
        vertex[u].edge.push_back(v);
        vertex[v].edge.push_back(u);
    }
    free(a);
    free(b);
}
void deleteedge()
{
    edge_count=0;
    for (int i = 0; i < vertex.size(); i++)
    {
        vector <int> x;
        x.swap(vertex[i].edge);
        sort(x.begin(),x.end(),cmp);
        int p=0,v=-1;
        while (!x.empty())
        {
            while (!x.empty()&&v==x.back())
            {
                x.pop_back();
            }
            if (x.empty()) break;
            v=x.back();
            x.pop_back();
            vertex[i].edge.push_back(v);
        }
        vertex[i].vertexID=i;
        edge_count+=vertex[i].edge.size();
    }
    vector<edge>().swap(edgelist);
}
void orientation()
{
    int *a =new int[vertex_count];
    for (int i = 0; i < vertex_count; i++)
    {
        a[vertex[i].vertexID]=i;
    }
    
    for (int i = 0; i < vertex_count; i++)
    {
        vector<int> x(vertex[i].edge);
        vertex[i].edge.clear();
        while (!x.empty())
        {
            int v=x.back();
            x.pop_back();
            if (a[v]>i) vertex[i].edge.push_back(v);
        }
    }
    free(a);
}

void reassignID()
{
    int k1=0,k2=-1,k3=-1;
    for (int i = 0; i < vertex_count; i++)
    {
        vertex[i].newid=-1;
        if (k2==-1 && vertex[i].edge.size()<=bounder)
            k2=i;
            
        if (k3==-1 && vertex[i].edge.size()<2)
            k3=i;
    }
    cout<<k2<<' '<<k3<<endl;
    int s1=k1,s2=k2,s3=k3;
    for (int i = 0; i < vertex_count; i++)
    {
        if (vertex[i].edge.size()<=2) break;
        for (int j = 0; j < vertex[i].edge.size(); j++)
        {
            int v=vertex[i].edge[j];
            if (vertex[v].newid==-1)
            {
                if (v>=s3) 
                {
                    vertex[v].newid=k3;
                    k3++;
                }
                else if (v>=s2) 
                {
                    vertex[v].newid=k2;
                    k2++;
                }
                else
                {
                    vertex[v].newid=k1;
                    k1++;
                }
            }
        }
    }
    for (int i = 0; i < vertex_count; i++)
    {
        int u=vertex[i].newid;
        if (u==-1)
        {
            if (i>=s3) 
            {
                vertex[i].newid=k3;
                k3++;
            }
            else if (i>=s2) 
            {
                vertex[i].newid=k2;
                k2++;
            }
            else
            {
                vertex[i].newid=k1;
                k1++;
            }
        }
    }
    vertexb.swap(vertex);
    vertex.resize(vertex_count);
    
    for (int i = 0; i < vertex_count; i++)
    {
        int u=vertexb[i].newid;
        
        for (int j = 0; j < vertexb[i].edge.size(); j++)
        {
            int v=vertexb[i].edge[j];
            v=vertexb[v].newid;
            // cout<<u<<' '<<v<<endl;
            vertex[u].edge.push_back(v);
        }
        vector<int>().swap(vertexb[i].edge);
    }
    
    // for (int i = 0; i < 10; i++)
    // {
    //     for (int j = 0; j < vertexb[i].edge.size(); j++)
    //         cout<<vertexb[vertexb[i].edge[j]].newid<<' ';
    //     cout<<endl;
    // }
    
}

void computeCSR()
{
    long long time=clock();
    int *a =new int[vertex_count];
    for (int i = 0; i < vertex_count; i++)
    {
        a[vertex[i].vertexID]=i;
    }
    for (int i = 0; i < vertex_count; i++)
    {
        for (int j = 0; j < vertex[i].edge.size(); j++)
        {
            vertex[i].edge[j]=a[vertex[i].edge[j]];
        }
        vertex[i].vertexID=i;
    }
    
    reassignID();
    time=clock()-time;
    cout<<"reassign time"<<time/CLOCKS_PER_SEC<<"s"<<endl;
    // time=clock();
    ofstream beginFile("begin.bin", ios::out|ios::binary);
    ofstream adjFile("adjacent.bin", ios::out|ios::binary);
    // ofstream edgeFile("edge", ios::out|ios::binary);
    // FILE *beginFile=fopen("begin.bin","wb");
    // FILE *adjFile=fopen("begin.bin","wb");
    long long sum=0;
    // vector<long long> beg;
    // vector<int> adj;

    for (int i = 0; i < vertex_count; i++)
    {
        beginFile.write((char*)&sum,sizeof(long long));
        // beg.push_back(sum);
        // fwrite(&sum,sizeof(long long),1,beginFile);
        sum+=vertex[i].edge.size();
        // for (int j=0;j<vertex[i].edge.size();j++)
        //     adj.push_back(vertex[i].edge[j]);
        adjFile.write((char*)&vertex[i].edge[0],sizeof(int)*vertex[i].edge.size());
        // fwrite(&vertex[i].edge[0],sizeof(int),vertex[i].edge.size(),beginFile);
    }
    // beg.push_back(sum);
    beginFile.write((char*)&sum,sizeof(long long));
    // edgeFile.write((char*)&sum,sizeof(long long));

    // fwrite(&beg,sizeof(long long),beg.size(),beginFile);
    // fwrite(&adj,sizeof(long long),adj.size(),adjFile);

    // time=clock()-time;
    // cout<<"write time"<<time/CLOCKS_PER_SEC<<"s"<<endl;
    // beginFile.write((char*)&sum,sizeof(long long));
    free(a);
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
        Outfilename=argv[2];
    }
    newloadgraph(Infilename);

    puttime("loadgraph");
    selectVertex();
    puttime("select time");
    deleteedge();

    // cout<<vertex_count<<" ???"<<endl;
    // for (int i = 0; i < vertex_count; i++)
    // {
    //     cout<<vertex[i].edge.size()<<endl;
    // }
    puttime("delete time");
    int *tst=new int[vertex_count];
    for (int i=0;i<vertex_count;i++)
        tst[i]=vertex[i].edge.size();
    sort(tst,tst+vertex_count);
    puttime("withoutvector");

    sort(vertex.begin(),vertex.end(),cmp1);
    
    puttime("sort time");
    //print
    // for (int i = 0; i < vertex_count; i++)
    // {
    //     cout<<vertex[i].edge.size()<<endl;
    // }
   
    orientation();
    
        //print
    
    puttime("orientation time");
    sort(vertex.begin(),vertex.end(),cmp2);

    
    puttime("sort time");
    cout<<vertex[0].edge.size()<<endl;

    // cout<<vertex[k].edge.size()<<endl;
    // cout<<vertex[k-1].edge.size()<<endl;
    // cout<<k<<endl;
    
    computeCSR();
    // for (int i = 0; i < k; i++)
    // {
    //     cout<<vertex[i].edge.size()<<endl;
    // }
    
    puttime("store time");

    system("pause");
    return 0;
} 