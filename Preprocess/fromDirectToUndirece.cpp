#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream> 
#include <cmath>

using namespace std;
typedef 	long long 	index_t;

typedef struct edge
{
    int u,v;
};
typedef struct edge_list
{
    int vertexID;
    vector<int> edge;
};
vector<edge_list> vertex;
vector<edge> edgelist;
int maxvertex=0;
int vertex_count,edge_count;
bool cmp(int a,int b)
{
    return a>b;
}
void loadgraph(string filename)
{
    ifstream inFile(filename.c_str(), ios::in);
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
        edge e;
        ss>>e.u>>e.v>>x;
        edgelist.push_back(e);
        maxvertex=max(maxvertex,max(e.u,e.v));
    }
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
            k++;
            a[i]=k;
        }
    }
    vertex_count=k;
    vertex.resize(k+1);
    for (int i=0;i<edgelist.size();i++)
    {
        int u=a[edgelist[i].u];
        int v=a[edgelist[i].v];
        if (u>v) swap(u,v);
        vertex[u].edge.push_back(v);
    }
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
        edge_count+=vertex[i].edge.size();
    }
    
}
void writeback(string filename)
{
    fstream outFile(filename.c_str(), ios::out);
    outFile<<vertex_count<<' '<<vertex_count<<' '<<edge_count<<endl;
    for (int i=0;i<vertex.size();i++)
    {
        for (int j=0;j<vertex[i].edge.size();j++)
            outFile<<i<<' '<<vertex[i].edge[j]<<endl;
    }
}

int main(int argc, char* argv[])
{
    string Infilename="1.mmio";
    string Outfilename="1.mmio";
    if (argc>1)
    {
        Infilename=argv[1];
    }
    if (argc>2)
    {
        Outfilename=argv[2];
    }
    loadgraph(Infilename);
    cout<<"loadok"<<endl;
    selectVertex();
    cout<<"selectok"<<endl;
    deleteedge();
    cout<<"deleteok"<<endl;
    writeback(Outfilename);
    cout<<"writebackok"<<endl;
    system("pause");
}