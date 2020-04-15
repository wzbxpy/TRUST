#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream> 
#define bounder 100
using namespace std;
typedef 	long long 	index_t;

int vertex_count, edge_count;
typedef struct edge_list
{
    int vertexID;
    vector<int> edge;
};
vector<edge_list> vertex;
vector<edge_list> vertexb;
bool cmp1(edge_list a, edge_list b)
{
    return a.edge.size()<b.edge.size() ;
}
bool cmp2(edge_list a, edge_list b)
{
    return a.edge.size()>b.edge.size() ;
}

int binary_search(int value)
{
	int l=0,r=vertex_count-1;
	while (l<r-1)
	{
		int mid=(l+r)>>1;
		if (vertex[mid].edge.size()>=value) l=mid;
			else r=mid;
	}
	// if (arr[r]<=value) return r;
	return l;
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
            ss>>x>>vertex_count>>edge_count;
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
}
void computeCSR(int k)
{
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
    
    ofstream beginFile("begin.bin", ios::out|ios::binary);
    ofstream adjFile("adjacent.bin", ios::out|ios::binary);
    ofstream edgeFile("edge", ios::out|ios::binary);
    long long sum=0;
    for (int i = 0; i < vertex_count; i++)
    {
        beginFile.write((char*)&sum,sizeof(long long));
        sum+=vertex[i].edge.size();
        // sort(vertex[i].edge.begin(),vertex[i].edge.end());
        vector<int>::iterator upp=upper_bound(vertex[i].edge.begin(),vertex[i].edge.end(),k);
        int divide=upp-vertex[i].edge.begin();
        // cout<<divide<<' '<<vertex[i].edge.size()<<endl;
        edgeFile.write((char*)&divide,sizeof(int));
        adjFile.write((char*)&vertex[i].edge[0],sizeof(int)*vertex[i].edge.size());
    }
    
    beginFile.write((char*)&sum,sizeof(long long));
}
void divideforresources()
{
    vertexb.swap(vertex);
    for (int i=0;i<vertex_count;i++)
    {
        if (vertexb[i].edge.size()>bounder)
        {
            vertex.push_back(vertexb[i]);
        }
    }
    for (int i=0;i<vertex_count;i++)
    {
        if (vertexb[i].edge.size()>=2&&vertexb[i].edge.size()<=bounder)
        {
            vertex.push_back(vertexb[i]);
        }
    }
    for (int i=0;i<vertex_count;i++)
    {
        if (vertexb[i].edge.size()<2)
        {
            vertex.push_back(vertexb[i]);
        }
    }
    vertexb.clear();
}
int main()
{
    loadgraph();
    sort(vertex.begin(),vertex.end(),cmp1);

    //print
    // for (int i = 0; i < vertex_count; i++)
    // {
    //     cout<<vertex[i].edge.size()<<endl;
    // }
  
    orientation();
    
    // sort(vertex.begin(),vertex.end(),cmp2);
    divideforresources();
    // for (int i = 0; i < vertex_count; i++)
    // {
    //     cout<<vertex[i].edge.size()<<endl;
    // }
  
    cout<<vertex[0].edge.size()<<endl;
    int k=binary_search(32);
    // cout<<k<<endl;
    // cout<<vertex[k].edge.size()<<endl;
    // cout<<vertex[k-1].edge.size()<<endl;
    // cout<<k<<endl;
    computeCSR(k);



    system("pause");
    return 0;
} 