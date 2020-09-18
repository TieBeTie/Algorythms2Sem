#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>

using std::cin;
using std::cout;
using std::vector;
using std::queue;


class Graph {
public:
	Graph() {}
	Graph(size_t vertex_count) {
		graph.assign(vertex_count, vector<size_t>());
	}
	Graph(size_t vertex_count, size_t edges) {
		graph.assign(vertex_count, vector<size_t>());
		this->edges = edges;
	}
	~Graph() = default;
	void AddEdge(size_t head, size_t tail) {
		graph[head].push_back(tail);
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
	const vector <size_t>& GetNext(size_t parent) const {
		return graph[parent];
	}
	const size_t GetNextSize(size_t parent) const {
		return graph[parent].size();
	}
private:
	vector <vector <size_t>> graph;
	size_t edges;
};

void DFS(const Graph& graph, const size_t v, const size_t p,
	vector <size_t>& tin, vector <size_t>& tout, vector < vector<size_t> >& up,
	size_t& timer, const size_t h) {

	up[v][0] = p;
	tin[v] = ++timer;

	for (size_t i = 1; i <= h; ++i) {
		up[v][i] = up[up[v][i - 1]][i - 1];
	}
	
	for (size_t i = 0; i < graph.GetNextSize(v); ++i) {
		size_t u = graph.GetNext(v)[i];
		if (u != p) {
			DFS(graph, u, v, tin, tout, up, timer, h);
		}
	}
	tout[v] = ++timer;
}

bool Upper(const size_t a, const size_t b, 
	const vector <size_t>& tin, const vector <size_t>& tout) {

	return tin[a] <= tin[b] && tout[a] >= tout[b];
}

size_t LCA(size_t a, size_t b, const size_t h, 
	const vector <size_t>& tin, const vector <size_t>& tout,
	const vector < vector<size_t> >& up) {

	if (a == b) {
		return a;
	}

	if (Upper(a, b, tin, tout)) {
		return a;
	}
	if (Upper(b, a, tin, tout)) {
		return b;
	}

	for (size_t i = h; i <= h; --i) {
		if (!Upper(up[a][i], b, tin, tout)) {
			a = up[a][i];
		}
	}
	return up[a][0];
}

int main() {

	size_t n, m;
	cin >> n >> m;
	
	Graph graph(n);
	size_t parent;
	for (size_t son = 1; son < n; ++son) {
		cin >> parent;
		graph.AddEdge(parent, son);
	}

	vector <size_t> tin(n, 0);
	vector <size_t> tout(n, 0);
	vector < vector<size_t> > up(n);

	size_t h = 1;
	while ((1<<h) <= n)  ++h;
	for (size_t i = 0; i < n; ++i) {
		up[i].resize(h + 1);
	}

	size_t timer = 0;
	DFS(graph, 0, 0, tin, tout, up, timer, h);

	size_t a, b, x, y, z;
	cin >> a >> b >> x >> y >> z;

	size_t result = LCA(a, b, h, tin, tout, up);
	size_t last = result;

	size_t a_tmp = (x * a + y * b + z) % n;
	size_t b_tmp = (x * b + y * a_tmp + z) % n;

	a = a_tmp;
	b = b_tmp;

	for (size_t i = 2; i <= m; ++i) {

		last = LCA((a + last) % n, b, h,
			tin, tout, up);
		result += last;
		a_tmp = (x * a + y * b + z) % n;
		b_tmp = (x * b + y * a_tmp + z) % n;
		
		a = a_tmp;
		b = b_tmp;
	}
	cout << result;
}