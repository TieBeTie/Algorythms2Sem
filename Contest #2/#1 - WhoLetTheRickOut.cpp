#include <iostream>
#include <vector>
#include <limits>
#include <queue>

using std::cin;
using std::cout;
using std::priority_queue;
using std::vector;
using std::pair;

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
	const vector <pair <size_t, size_t>>& GetNext(size_t parent) const {
		return graph[parent];
	}
	size_t GetNextSize(size_t parent) const {
		return graph[parent].size();
	}
private:
	vector <vector< pair<size_t, size_t>>> graph;
	size_t edges;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

struct Comparator {
	bool operator() (const pair <size_t, size_t>& a, const pair <size_t, size_t>& b) const {
		return a.second > b.second;
	}
};

void Dijkstra(const WDGraph& graph, size_t s, 
              vector <size_t>& d, size_t pay_1, size_t pay_2) {

	d[s] = 0;
	priority_queue < pair <size_t, size_t>, vector < pair <size_t, size_t> >, Comparator> q;

	q.emplace(s, d[s]);
	size_t to, cost, v, cost_v;

	while (!q.empty()) {
		v = q.top().first;
		cost_v = q.top().second;
		q.pop();
		if (cost_v > d[v]) { continue; }

		for (size_t i = 0; i < graph.GetNextSize(v); ++i) {
			to = graph.GetNext(v)[i].first;
			cost = graph.GetNext(v)[i].second;

			if (d[v] + cost < d[to]) {
				d[to] = d[v] + cost;
				q.emplace(to, d[to]);
			}
		}
	}
}

void CreateMultiUniverseGraph(WDGraph& multi_universe, size_t pay_1, size_t pay_2, size_t start, size_t finish) {
	for (size_t i = 0; i < multi_universe.size(); ++i) {
		multi_universe.AddEdge(i, (i + 1) % multi_universe.size(), pay_1);
	}

	for (size_t i = 0; i < multi_universe.size(); ++i) {
		multi_universe.AddEdge(i, (i * i + 1) % multi_universe.size(), pay_2);
	}
}

void Read5Nums(size_t& n1, size_t& n2, size_t& n3, size_t& n4, size_t& n5) {
	cin >> n1 >> n2 >> n3 >> n4 >> n5;
}

void PrintAnswer(size_t answer) {
	cout << answer;
}


int main() {
	InOutOptimizaton();
	size_t pay_1, pay_2, count_universes, start, finish;
	Read5Nums(pay_1, pay_2, count_universes, start, finish);
	WDGraph multi_universe(count_universes);

	CreateMultiUniverseGraph(multi_universe, pay_1, pay_2, start, finish);
	vector <size_t> distance(multi_universe.size(), std::numeric_limits<size_t>::max());
	Dijkstra(multi_universe, start, distance, pay_1, pay_2);
	PrintAnswer(distance[finish]);
}
