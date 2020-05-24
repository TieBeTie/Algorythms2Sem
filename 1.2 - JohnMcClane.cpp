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
	const vector <int> &GetNext(int parent) const {
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

void ReadGraph(Graph& graph) {
	int head, tail;
	for (int i = 0; i < graph.Edges(); i++) {
		cin >> head >> tail;
		graph.AddEdge(head, tail);
	}
}

void InitParam(int* n, int* m) {
	cin >> *n >> *m;
}

void PrintAnswer(bool answer) {
	cout << answer;
}

void DFS(const Graph &graph, vector <char> &colors, int v, 
vector <int> &linear_ordering, bool &is_cyclic) {

	colors[v] = 'g';
	for (int u : graph.GetNext(v)) {
		if (colors[u] == 'w') {
			DFS(graph, colors, u, linear_ordering, is_cyclic);
		} else if (colors[u] == 'g') {
			is_cyclic = true;
			return;
		}
	}
	colors[v] = 'b';
	linear_ordering.push_back(v);
}

void TopologicalSort(const Graph& graph, bool& is_cyclic, vector <int> &linear_ordering) {
	vector <char> colors(graph.size(), 'w');
	for (int i = 0; i < graph.size(); ++i) {
		if (colors[i] == 'w') {
			DFS(graph, colors, i, linear_ordering, is_cyclic);
		}
	}
}

void PrintReport(bool anxiety, const vector <int>& call_order) {
	if (anxiety) {
		cout << "NO";
	} else {
		cout << "YES" << "\n";
		for (int i = call_order.size() - 1; i >= 0; i--) {
			cout << call_order[i] << " ";
		}
	}
}

int main() {
	InOutOptimizaton();
	int n, m;
	InitParam(&n, &m);
	Graph policemen (n, m);
	ReadGraph(policemen);

	bool anxiety = false;
	vector <int> call_order;
	TopologicalSort(policemen, anxiety, call_order);

	PrintReport(anxiety, call_order);
}
