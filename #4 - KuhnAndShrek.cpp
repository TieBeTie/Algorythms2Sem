#include <iostream>
#include <vector>
#include <limits>
#include <queue>

using std::cin;
using std::cout;
using std::priority_queue;
using std::vector;
using std::pair;

class Graph {
public:
	Graph() {}
	Graph(size_t vertex_count) {
		graph.assign(vertex_count, vector< size_t >());
	}
	Graph(size_t vertex_count, size_t edges) {
		graph.assign(vertex_count, vector< size_t >());
		this->edges = edges;
	}
	~Graph() = default;
	void AddEdge(size_t head, size_t tail) {
		graph[head].emplace_back(tail); //.push_back(std::pair<size_t, size_t> (tail, cost))
	}
	bool HasEdge(size_t head, size_t tail) const {
		for (size_t i = 0; i < graph[head].size(); i++)
			if (graph[head][i] == tail)
				return true;
	}
	size_t Edges() const {
		return edges;
	}
	size_t size() const {
		return graph.size();
	}
	const vector< size_t >& GetNext(size_t parent) const {
		return graph[parent];
	}
	size_t GetNextSize(size_t parent) {
		return graph[parent].size();
	}
private:
	vector <vector< size_t >> graph;
	size_t edges;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

struct Comparator {
	bool operator() (const pair <size_t, size_t>& a, 
                     const pair <size_t, size_t>& b) const {
		return a.second > b.second;
	}
};

void PrintAnswer(size_t answer) {
	cout << answer;
}

void ReadMap(vector <char>& map, int n, int m, int& holes, vector <bool>& black) {
	char c;
	int chess = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			cin >> c;
			map.push_back(c);
			if ((i + j) % 2 == 1) {
				black[chess] = true;
			}
			if (c == '*') {
				holes += 1;
			}
			chess++;
		}
	}
}

void AddEdgeToWhite(vector <char>& map, int i,
	Graph& graph, int n, int m) {
	if (i % m != 0 && map[i - 1] == '*') { //left
		graph.AddEdge(i, i - 1);
	}

	if (i % m != m - 1 && map[i + 1] == '*') { //right
		graph.AddEdge(i, i + 1);
	}

	if (i >= m  && map[i - m] == '*') { //down
		graph.AddEdge(i, i - m);
	}

	if (i < (n - 1) * m && map[i + m] == '*') { //up
		graph.AddEdge(i, i + m);
	}
}

void CreateChessGraph(vector <char>& map, Graph& graph, int n, int m) {
	for (int i = 0; i < n * m; ++i) {
		if (map[i] == '*') {
			AddEdgeToWhite(map, i, graph, n, m);
		}
	}
}

bool DFS(const Graph& graph, vector <char>& color, int v, 
         vector <int>& matching) {

	if (color[v] == 'b') {
		return false;
	}

	color[v] = 'b';

	for (int u : graph.GetNext(v)) {
		if (matching[u] == -1 || DFS(graph, color, matching[u], matching)) {
			matching[u] = v;
			return true;
		}
	}

	return false;
}

void Kuhn(Graph& graph, vector <int>& matching, int& answer, vector <bool>& black) {
	for (size_t i = 0; i < matching.size(); ++i) {
		matching[i] = -1;
	}

	vector <char> color (graph.size());
	for (size_t i = 0; i < graph.size(); ++i) {
		for (size_t j = 0; j < graph.size(); ++j) {
			color[j] = 'w';
		}
		if (black[i]) {
			answer += DFS(graph, color, i, matching);
		}
	}
}


int main() {
	int n, m, a, b;
	cin >> n >> m >> a >> b;

	vector <char> map;
	vector <bool> black (m * n, false);

	int holes = 0;
	ReadMap(map, n, m, holes, black);
	if (a < 2 * b) {
		Graph graph(n * m);
		CreateChessGraph(map, graph, n, m);

		int answer = 0;
		vector <int> matching(n * m);
		Kuhn(graph, matching, answer, black);
		cout << (answer * a + (holes - 2 * answer) * b);
	}
	else {
		cout << holes  * b;
	}
}
