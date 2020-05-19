#include <iostream>
#include <vector>
#include <queue>
#include <limits>

using std::cin;
using std::cout;
using std::vector;
using std::queue;

class Graph {
public:
	Graph() {}
	Graph(int vertex_count) {
		graph.assign(vertex_count, vector<int>());
	}
	Graph(int vertex_count, int edges) {
		graph.assign(vertex_count, vector<int>());
		this->edges = edges;
	}
	~Graph() = default;
	void AddEdge(int head, int tail) {
		graph[head].push_back(tail);
		graph[tail].push_back(head);
	}
	bool HasEdge(int head, int tail) const {
		for (int i = 0; i < graph[head].size(); i++)
			if (graph[head][i] == tail)
				return true;
	}
	int Edges() const {
		return edges;
	}
	int size() const {
		return graph.size();
	}
	const vector <int>& GetNext(int parent) const {
		return graph[parent];
	}
private:
	vector <vector <int>> graph;
	int edges;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

void InitNums(int* num_1, int* num_2, int* num_3, int* num_4, int* num_5) {
	cin >> *num_1;
	cin >> *num_2;
	cin >> *num_3;
	*num_3 -= 1;
	cin >> *num_4;
	*num_4 -= 1;
	cin >> *num_5;
	*num_5 -= 1;
}

void ReadGraph(Graph &graph) {
	int head, tail;
	for (int i = 0; i < graph.Edges(); i++)  {	
		cin >> head >> tail;
		graph.AddEdge(head - 1, tail - 1);
	}
}

void BFS(const Graph &graph, int start, int* distance) {
	queue <int> jail;
	vector <int> children;
	int v = start;
	for (int i = 0; i < graph.size(); ++i) {
		distance[i] = -1;
	}
	distance[start] = 0;
	bool* visited = new bool [graph.size()];

	for (int i = 0; i < graph.size(); i++) 
		visited[i] = false;

	jail.push(start);
	while (!jail.empty())
	{
		v = jail.front();
		jail.pop();
		visited[v] = true;
		children = graph.GetNext(v);
		for (int i = 0; i < children.size(); i++)
		{
			if (distance[ children[i] ] == -1)
			{
				distance[ children[i] ] = distance[v] + 1;
				jail.push(children[i]);
			}
		}
	}
	delete[] visited;
}

void PrintAnswer(const int &answer) {
	cout << answer;
}

int main() {
	InOutOptimizaton();
	int n, m, leon, matilda, milk;
	InitNums(&n, &m, &leon, &matilda, &milk);
	Graph graph(n, m);
	ReadGraph(graph);

	int* leon_distance = new int[n];
	BFS(graph, leon, leon_distance);
	int* matilda_distance = new int[n];
	BFS(graph, matilda, matilda_distance);
	int* milk_distance = new int[n];
	BFS(graph, milk, milk_distance);

	int min_distance = std::numeric_limits<int>::max();
	for (int i = 0; i < n; i++)
		if ( (milk_distance[i] + matilda_distance[i] + leon_distance[i]) < min_distance)
			min_distance = milk_distance[i] + matilda_distance[i] + leon_distance[i];

	delete[] leon_distance;
	delete[] matilda_distance;
	delete[] milk_distance;
	PrintAnswer(min_distance);

	return 0;
}
