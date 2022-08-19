#include<cstdlib>
#include <ctime>
#include <map>
#include <vector>
#include <iostream>
#include <queue>
#include <fstream>

using namespace std;

vector<int> bfs(int vert, int &branch_min, float weights[], bool visited[], const vector<int> neighborhoods[], int pre[])
{
    //std::cout << "finding branch min for " << vert << std::endl;

	clock_t d_start = clock();
    	
	queue<int> q;
	vector<int> visited_list;
	visited_list.clear();
	visited_list.push_back(vert);
	pre[vert] = -1;
	q.push(vert);

	//cout << "entering queue loop" << endl;
	//branch_min = vert;
	branch_min = vert;

	int calls = 0;

	while (!q.empty())
	{
		calls++;
		int current_vert = q.front();
		q.pop();
	
		bool is_visited = visited[current_vert];	
		if (!is_visited)
		{
			clock_t d_start = clock();
			visited[current_vert] = true;
			clock_t d_end = clock();
			//cout << "discovered: " << (d_end - d_start) / float(CLOCKS_PER_SEC) << endl;
		}

		//std::cout << "current min: " << branch_min << " has weight of " << weights[branch_min] << std::endl;

		float current_val = weights[current_vert];
		if (current_val < weights[branch_min])
		{
			clock_t m_start = clock();
			branch_min = current_vert;
			clock_t m_end = clock();
            //std::cout << "new min: " << branch_min << " has weight of " << weights[branch_min] << std::endl;
			//cout << "min: " << (m_end - m_start) / float(CLOCKS_PER_SEC) << endl;
		}
		else if (current_val == weights[branch_min] && current_vert < branch_min)
		{
			branch_min = current_vert;
            //std::cout << "new min: " << branch_min << " has weight of " << weights[branch_min] << std::endl;

		}

		vector<int> neighborhood = neighborhoods[current_vert];

		for (int i = 0; i < neighborhood.size(); i++)
		{
			//cout << "working on neighbor " << i << endl;
			int neighbor = neighborhood[i];
			if (!visited[neighbor])
			{
			clock_t m_start = clock();
				visited[neighbor] = true;
				visited_list.push_back(neighbor);
				pre[neighbor] = current_vert;
				q.push(neighbor);
			clock_t m_end = clock();
			//cout << "min: " << (m_end - m_start) / float(CLOCKS_PER_SEC) << endl;
			}
		}
	}
	for (int i = 0; i < visited_list.size(); i++)
	{
		visited[visited_list[i]] = false;
	}
	//cout << "number of calls: " << calls << endl;
			clock_t d_end = clock();
			//cout << "whole call: " << (d_end - d_start) / float(CLOCKS_PER_SEC) << endl;
	return visited_list;
}

void retrieve_path(int vert, vector<int> &vPath, int pre[]){
    vPath.clear();
    //cout << "starting at vert " << vert << endl;
    while(pre[vert] != -1){
        vPath.push_back(vert);
        int vpre = pre[vert];
	//cout << vert << " " << vpre << endl;
        vert = vpre;
    }
    vPath.push_back(vert);
}

int main(int argc, char* argv[])
{
	string weights_filename = argv[1];
	string input_edge_filename = argv[2];
	float persistence_threshold = atof(argv[3]);
    //float weight_thresh = atof(argv[4]);
	string output_dir = argv[4];

    ifstream fin;
    
    vector<float> weights;
    weights.clear();
    float r;
    fin.open(weights_filename.c_str());
    cout << "reading in weights from: " << weights_filename << endl;
    while (fin >> r)
    {
        //std::cout << "weight: " << r << std::endl;
        weights.push_back(r);
    }
    fin.close();
    int num_verts = weights.size();
    std::cout << "read in " << num_verts << " verts" << std::endl;

    float* weight_arr = new float[num_verts];
    for (int i = 0; i < num_verts; i++)
    {
        weight_arr[i] = weights[i];
    }

	vector<vector<int> > edges;
	edges.clear();
    std::vector<float> edge_persistence;
    edge_persistence.clear();
	int u, v;
    string type;
    float i;
    fin.open(input_edge_filename.c_str());
    cout << "reading in edges from: " << input_edge_filename << endl;
	int cnt = 0;
    while (fin >> u >> v >> type >> i)
    {
	    cnt++;
        
        //std::cout << "edge: " << u << " " << v << " " << type << " " << i << std::endl;

        vector<int> edge;
    	edge.clear();
	    edge.push_back(u);
	    edge.push_back(v);
        if (type == "0")
        {
            i = -i - 1;
        }
	
        if (type == "2")
        {
            i = 999;
            cout << "persistent loop in domain" << endl;
        }
		edge_persistence.push_back(i);
		edges.push_back(edge);
	}
    fin.close();
	std::cout << "read in " << edges.size() << " edges." << std::endl;

    cout << "Computing vector field" << endl;
    vector<int>* neighborhoods = new vector<int>[num_verts];
    for (int i = 0; i < num_verts; i++)
    {
        vector<int> neighborhood;
        neighborhood.clear();
        neighborhoods[i] = neighborhood;
    }

	std::cout << "initialized neighbors" << std::endl;

	vector<vector<int> > vf_edges;
    vf_edges.clear();
    int ve_in_vf = 0;
    for (int i = 0; i < edges.size(); i++)
    {
		vector<int> edge = edges[i];
        float persistence = edge_persistence[i];
		if (persistence < 0)
		{
			persistence = -(persistence + 1);
		}
		else
		{
			continue;
		}
		if (persistence > persistence_threshold)
		{
			continue;
		}
		

        ve_in_vf++;
        int v0 = edge[0];
        int v1 = edge[1];
        vector<int> field_edge;
        field_edge.clear();
        field_edge.push_back(v0);
        field_edge.push_back(v1);
        vf_edges.push_back(field_edge);
        neighborhoods[v0].push_back(v1);
        neighborhoods[v1].push_back(v0);
    }


	cout << "edges in vector field: " << ve_in_vf << endl;

	cout << "Computing manifold" << endl;
    vector<int> min_computed;
    min_computed.clear();

    bool* visit = new bool[num_verts];

    int* next_in_path = new int[num_verts];

    vector<vector<int> > manifold;
    manifold.clear();

    vector<int> manifold_type;
    manifold_type.clear();

    for (int i = 0; i < num_verts; i++)
    {
        min_computed.push_back(-1);
        visit[i] = false;
        next_in_path[i] = -1;
    }


	int know_min = 0;
    int not_know_min = 0;
    int critical_count = 0;

    int total_pos = 0;
    int total_neg = 0;

    for (int i = 0; i < edges.size(); i++)
    {
        if (i % 10000 == 0)
        {
	       //cout << "working on edge " << i << " out of " << edges.size() << endl;
	    }

    	vector<int> edge = edges[i];
        float persistence = edge_persistence[i];
        bool pos;
    	if (persistence < 0)
        {
            pos = false;
            persistence = -(persistence + 1);
        }
        else
        {
            pos = true;
        }
        if (persistence <= persistence_threshold)
        {
		  continue;
        }

        int v0 = edge[0];
        int v1 = edge[1];
	/*
        //float weight_thresh = 29;
        if (weights[v0] > weight_thresh or weights[v1] > weight_thresh)
        {
            std::cout << "skipping due to low density" << std::endl;
            continue;
        }
	*/

        if (pos)
        {
            total_pos++;
            std::cout << "positive: " << persistence << std::endl;
        }
        else
        {
            total_neg++;
            //continue;
            std::cout << "negative: " << persistence << std::endl;
        }

        
        critical_count++;

        vector<int> critical_edge;
        critical_edge.clear();
        critical_edge.push_back(edge[0]);
        critical_edge.push_back(edge[1]);
        manifold.push_back(critical_edge);
        if (pos)
        {
            manifold_type.push_back(2);
        }
        else
        {
            manifold_type.push_back(1);
        }

    	for (int j = 0; j < 2; j++)
        {
            int v = edge[j];
            vector<int> vPath;
            vPath.clear();
            if (min_computed[v] == -1)
            {
                not_know_min++;
                int branch_min;
                vector<int> component = bfs(v, branch_min, weight_arr, visit, neighborhoods, next_in_path);
		        for (int k = 0; k < component.size(); k++)
                {
                    min_computed[component[k]] = branch_min;
                }
                bfs(branch_min, branch_min, weight_arr, visit, neighborhoods, next_in_path);
                retrieve_path(v, vPath, next_in_path);
            }
            else
            {
                know_min++;
                retrieve_path(v, vPath, next_in_path);
            }
            manifold.push_back(vPath);
            manifold_type.push_back(-1);
	   }
    }

    std::cout << "Positive: " << total_pos << std::endl; 
    std::cout << "Negative: " << total_neg << std::endl; 

	delete[] neighborhoods;
	//delete[] vert_arr;

	cout << "outputting..." << endl;
    vector<vector<int> > output_edges;
    output_edges.clear();

	for (int i = 0; i < manifold.size(); i++)
    {
        vector<int> component = manifold[i];
        for (int j = 0; j < component.size() - 1; j++)
        {
                        //cout << "beginning of i loop: " << output_index << endl;
            int v0 = component[j];
            int v1 = component[j + 1];

            vector<int> edge;
            edge.clear();
            edge.push_back(v0);
            edge.push_back(v1);
            edge.push_back(manifold_type[i]);
            output_edges.push_back(edge);
        }
    }

    cout << "writing files" << endl;

    string edge_filename = output_dir + "dimo_edge.txt";
    ofstream eFile(edge_filename.c_str());
    for (int i = 0; i < output_edges.size(); i++)
    {
        vector<int> edge = output_edges[i];
        eFile << edge[0] << " " << edge[1] << " " << edge[2] << endl;
    }

    return 0;
}
