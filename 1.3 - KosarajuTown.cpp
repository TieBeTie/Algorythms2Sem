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
	int GetNextSize(int parent) const {
		return graph[parent].size();
	}
private:
	vector <vector <int>> graph;
	int edges;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

void ReadGraphWithInvert(Graph& graph, Graph& graph_invert) {
	int head, tail;
	for (int i = 0; i < graph.Edges(); i++) {
		cin >> head >> tail;
		graph.AddEdge(head - 1, tail - 1);
		graph_invert.AddEdge(tail - 1, head - 1);
	}
}

void InitParam(int* n, int* m) {
	cin >> *n >> *m;
}

void PrintAnswer(int answer) {
	cout << answer;
}

void DFSCLock(const Graph& graph, vector <char>& colors, int v, vector <int>& out_time) {
	colors[v] = 'b';
	for (int u : graph.GetNext(v)) {
		if (colors[u] == 'w') {
			DFSCLock(graph, colors, u, out_time);
		}
	}
	out_time.push_back(v);
}

void DFSTimeMachine(const Graph& graph, vector <char>& colors, int v, vector <int>& components, int component_num) {
	colors[v] = 'b';
	components[v] = component_num;
	for (int u : graph.GetNext(v)) {
		if (colors[u] == 'w') {
			DFSTimeMachine(graph, colors, u, components, component_num);
		}
	}
}

void Kosaraju(const Graph& graph, Graph& graph_invert, vector <int>& components, int& component_num) {
	vector <char> colors;
	vector <int> out_time;
	colors.assign(graph.size(), 'w');

	for (int vertex = 0; vertex < graph.size(); vertex++) {
		if (colors[vertex] == 'w') {
			DFSCLock(graph, colors, vertex, out_time);
		}
	}

	colors.assign(graph.size(), 'w');
	components.assign(graph.size(), 0);
	component_num = 0;

	for (int i = out_time.size() - 1; i >= 0; --i) {
		if (colors[ out_time[i] ] == 'w') {
			DFSTimeMachine(graph_invert, colors, out_time[i], components, component_num);
			component_num++;
		}
	}
}

void CondenseGraph(const Graph& graph, Graph& condensed_graph, vector <int>& components) {
	for (int parent = 0; parent < graph.size(); ++parent) {
		for (int children : graph.GetNext(parent)) {
			if (components[parent] != components[children]) {
				condensed_graph.AddEdge(components[parent], components[children]);
			}
		}
	}
}

void InvertGraph(const Graph& graph, Graph& inverted_graph) {
	for (int parent = 0; parent < graph.size(); ++parent) {
		for (int children : graph.GetNext(parent)) {
			inverted_graph.AddEdge(children, parent);
		}
	}
}
int MaxInOutOfVertex(Graph& graph) {
	if (graph.size() == 1) {
		return 0;
	}

	int in = 0, out = 0;
	vector <int> children;
	Graph inverted_graph(graph.size(), graph.Edges());
	InvertGraph(graph, inverted_graph);
	
	for (int parent = 0; parent < graph.size(); ++parent) {
		if (graph.GetNextSize(parent) == 0) {
			out++;
		}
	}
	for (int parent = 0; parent < inverted_graph.size(); ++parent) {
		if (inverted_graph.GetNextSize(parent) == 0) {
			in++;
		}
	}
	return in > out ? in : out;
}

int main() {
	InOutOptimizaton();
	int n, m;
	InitParam(&n, &m);
	Graph roads(n, m);
	Graph roads_invert(n, m);
	ReadGraphWithInvert(roads, roads_invert);
	
	vector <int> map;
	int districts_count;
	Kosaraju(roads, roads_invert, map, districts_count);
	Graph districts(districts_count);
	CondenseGraph(roads, districts, map);
	PrintAnswer(MaxInOutOfVertex(districts));
}
