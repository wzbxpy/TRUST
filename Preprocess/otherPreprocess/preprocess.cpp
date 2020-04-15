#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream> 
#include <cstring>

using namespace std;
typedef 	long long 	index_t;
int vertex_count, edge_count;
typedef struct edge_list
{
    int vertexID;
    vector<int> edge;
    int initaldegree;
};
vector<edge_list> vertex;
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
        //for H-index
        for (int j=0;j<vertex[i].edge.size();j++)
        {
            edgeFile.write((char*)&i,sizeof(int));
            edgeFile.write((char*)&vertex[i].edge[j],sizeof(int));
        }
        adjFile.write((char*)&vertex[i].edge[0],sizeof(int)*vertex[i].edge.size());
    }
    
    beginFile.write((char*)&sum,sizeof(long long));
}
int averagedegree(double result[],int num[])
{
    memset(result,0,sizeof(result));
    memset(num,0,sizeof(num));
    for (int i=0;i<vertex_count;i++)
    {
        result[vertex[i].initaldegree]+=vertex[i].edge.size();
        num[vertex[i].initaldegree]++;
    }
    int k=1;
    while (num[k]>0)
    {
        result[k]/=num[k];
        k++;
    }
    return k;
}
int main()
{
    loadgraph();
    sort(vertex.begin(),vertex.end(),cmp1);

    double* result=new double[vertex_count+1];
    int* num=new int[vertex_count+1];
    for (int i=0;i<vertex_count;i++)
    {
        vertex[i].initaldegree=vertex[i].edge.size();
    }
    //print

    orientation();
    
    int maxdegree=averagedegree(result,num);
    cout<<"after orientation"<<endl;
    
    // for (int i = 1; i < maxdegree; i++)
    // {
    //     cout<<result[i]<<endl;
    // }
    // cout<<endl;
    // cout<<endl;
    // cout<<endl;
    sort(vertex.begin(),vertex.end(),cmp2);
    
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