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
	size_t cost; 
	size_t step;
	Node(size_t vertex, size_t cost, size_t step) {
		this->vertex = vertex;
		this->cost = cost;
		this->step = step;
	}
};

class WDGraph {
public:
	WDGraph() {}
	WDGraph(size_t vertex_count) {
		graph.assign(vertex_count, vector< pair <size_t, size_t> >());
	}
	WDGraph(size_t vertex_count, size_t edges) {
		graph.assign(vertex_count, vector< pair <size_t, size_t> >());
		this->edges = edges;
	}
	~WDGraph() = default;
	void AddEdge(size_t head, size_t tail, size_t cost) {
		graph[head].emplace_back(tail, cost); //.push_back(std::pair<size_t, size_t> (tail, cost))
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
	const vector <pair <size_t, size_t> >& GetNext(size_t parent) const {
		return graph[parent];
	}
	size_t GetNextSize(size_t parent) const {
		return graph[parent].size();
	}
private:
	vector <vector< pair <size_t, size_t> > > graph;
	size_t edges;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

struct Comparator {
	bool operator() (const Node & a, const Node & b) const {
		return a.cost > b.cost;
	}
};

void Dijkstra(const WDGraph& graph, size_t s, 
              size_t step_limit, vector <size_t>& d) {

	d[s] = 0;
	priority_queue < Node, vector < Node >, Comparator> q;
	
	Node tmp(s, d[s], 0);
	q.push(tmp);
	size_t to, cost, v, cost_v, step_v;

	while (!q.empty()) {
		v = q.top().vertex;
		cost_v = q.top().cost;
		step_v = q.top().step;

		q.pop();
		if (cost_v > d[v] || step_v >= step_limit) { continue; }

		for (size_t i = 0; i < graph.GetNextSize(v); ++i) {
			to = graph.GetNext(v)[i].first;
			cost = graph.GetNext(v)[i].second;

			if (d[v] + cost < d[to]) {
				d[to] = d[v] + cost;
				tmp.vertex = to;
				tmp.cost = d[to];
				tmp.step = step_v + 1;
				q.push(tmp);
			}
		}
	}
}

void CreateMultiUniverseGraph(WDGraph& multi_universe) {
	size_t n = multi_universe.Edges();
	size_t head, tail, cost;
	for (size_t i = 0; i < n; ++i) {
		cin >> head >> tail >> cost;
		multi_universe.AddEdge(head - 1, tail - 1, cost);
	}
}

void Read5Nums(size_t& n1, size_t& n2, size_t& n3, size_t& n4, size_t& n5) {
	cin >> n1 >> n2 >> n3 >> n4 >> n5;
}

void PrintAnswer(long long answer) {
	cout << answer;
}


int main() {
	InOutOptimizaton();
	size_t count_universes, hyperjumps, max_hyperjumps, start, finish;
	Read5Nums(count_universes, hyperjumps, max_hyperjumps, start, finish);
	start -= 1; finish -= 1;
	WDGraph multi_universe(count_universes, hyperjumps);

	CreateMultiUniverseGraph(multi_universe);
	vector <size_t> distance(multi_universe.size(), std::numeric_limits<size_t>::max());
	Dijkstra(multi_universe, start, max_hyperjumps, distance);
	if (distance[finish] == std::numeric_limits<size_t>::max()) {
		PrintAnswer(-1);
	} else {
		PrintAnswer(distance[finish]);
	}
}
