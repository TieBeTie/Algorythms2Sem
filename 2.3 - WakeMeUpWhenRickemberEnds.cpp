#include <iostream>
#include <vector>
#include <limits>
#include <queue>

using std::cin;
using std::cout;
using std::priority_queue;
using std::vector;
using std::pair;

struct Node {
	size_t vertex;
	double chance;
	Node(size_t vertex, double chance) {
		this->vertex = vertex;
		this->chance = chance;
	}
};

class WDGraph {
public:
	WDGraph() {}
	WDGraph(size_t vertex_count) {
		graph.assign(vertex_count, vector< pair <size_t, double> >());
	}
	WDGraph(size_t vertex_count, size_t edges) {
		graph.assign(vertex_count, vector< pair <size_t, double> >());
		this->edges = edges;
	}
	~WDGraph() = default;
	void AddEdge(size_t head, size_t tail, double cost) {
		graph[head].emplace_back(tail, cost); //.push_back(std::pair<size_t, double> (tail, cost))
	}
	bool HasEdge(size_t head, size_t tail) const {
		for (size_t i = 0; i < graph[head].size(); i++)
			if (graph[head][i].first == tail)
				return true;
	}
	size_t Edges() const {
		return edges;
	}
	size_t size() const {
		return graph.size();
	}
	const vector <pair <size_t, double> >& GetNext(size_t parent) const {
		return graph[parent];
	}
	size_t GetNextSize(size_t parent) const {
		return graph[parent].size();
	}
private:
	vector <vector< pair <size_t, double> > > graph;
	size_t edges;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

struct Comparator {
	bool operator() (const Node& a, const Node& b) const {
		return a.chance > b.chance;
	}
};

void Dijkstra(const WDGraph& graph, size_t s, vector <double>& d) {
	d[s] = 0;
	priority_queue < Node, vector < Node >, Comparator> q;

	Node tmp(s, d[s]);
	q.push(tmp);
	size_t to, v;
	double chance, chance_v;

	while (!q.empty()) {
		v = q.top().vertex;
		chance_v = q.top().chance;

		q.pop();
		if (chance_v > d[v]) { continue; }

		for (size_t i = 0; i < graph.GetNextSize(v); ++i) {
			to = graph.GetNext(v)[i].first;
			chance = graph.GetNext(v)[i].second;

			if (d[v] + chance - d[v] * chance < d[to]) {
				d[to] = d[v] + chance - d[v] * chance;
				tmp.vertex = to;
				tmp.chance = d[to];
				q.push(tmp);
			}
		}
	}
}

void CreateMultiUniverseGraph(WDGraph& multi_universe) {
	size_t n = multi_universe.Edges();
	size_t head, tail;
	double chance;
	for (size_t i = 0; i < n; ++i) {
		cin >> head >> tail >> chance;
		chance /= 100.0;
		multi_universe.AddEdge(head - 1, tail - 1, chance);
		multi_universe.AddEdge(tail - 1, head - 1, chance);
	}
}

void Read4Nums(size_t& n1, size_t& n2, size_t& n3, size_t& n4) {
	cin >> n1 >> n2 >> n3 >> n4;
}

void PrintAnswer(double answer) {
	cout << answer;
}


int main() {
	InOutOptimizaton();
	size_t count_universes, hyperjumps, start, finish;
	Read4Nums(count_universes, hyperjumps, start, finish);
	WDGraph multi_universe(count_universes, hyperjumps);
	start -= 1; finish -= 1;
	CreateMultiUniverseGraph(multi_universe);
	vector <double> chances(multi_universe.size(), std::numeric_limits<double>::max());
	Dijkstra(multi_universe, start, chances);
	PrintAnswer(chances[finish]);
}
