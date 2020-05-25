#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <set>

using std::cin;
using std::cout;
using std::set;
using std::vector;
using std::pair;

class WGraph {
public:
	WGraph() {}
	WGraph(size_t vertex_count) {
		graph.assign(vertex_count, vector< pair <size_t, size_t> >());
	}
	WGraph(size_t vertex_count, size_t edges) {
		graph.assign(vertex_count, vector< pair <size_t, size_t> >());
		this->edges = edges;
	}
	~WGraph() = default;
	void AddEdge(size_t head, size_t tail, size_t cost) {
		graph[head].emplace_back(tail, cost);
		graph[tail].emplace_back(head, cost);
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

struct Comparator {
	bool operator() (const pair <size_t, size_t>& a, 
                     const pair <size_t, size_t>& b) const {
		return a.first == b.first;
	}
};

void PrimFindMST(const WGraph& graph, vector <size_t>& min_cost,
	vector <long long>& min_way) {

	set < pair <size_t, size_t>> q;
	size_t to, v, cost;
	min_cost[0] = 0;
	vector <bool> used(graph.size(), false);
	q.emplace(0, 0);

	

	while (!q.empty()) {
		v = q.begin()->second;
		q.erase(q.begin());

		for (size_t j = 0; j < graph.GetNextSize(v); ++j) {
			used[v] = true;
			to = graph.GetNext(v)[j].first;
			cost = graph.GetNext(v)[j].second;

			if (cost < min_cost[to] && !used[to]) {
				q.erase(pair <size_t, size_t>(min_cost[to], to));
				min_cost[to] = cost;
				min_way[j] = to;
				q.emplace(cost, to);
			}
		}
	}
}


void ReadGraph(WGraph& graph) {
	size_t head, tail, cost;
	for (size_t i = 0; i < graph.Edges(); ++i) {
		cin >> head >> tail >> cost;
		graph.AddEdge(head - 1, tail - 1, cost);
	}
}



int main() {
	size_t n, m;
	cin >> n >> m;
	WGraph graph(n, m);
	ReadGraph(graph);

	vector <size_t> min_cost(n, std::numeric_limits<size_t>::max());
	vector <long long> min_way(n, -1);
	PrimFindMST(graph, min_cost, min_way);
	
	size_t summ = 0;
	for (size_t i = 0; i < n; ++i) {
		summ += min_cost[i];
	}

	cout << summ;

}