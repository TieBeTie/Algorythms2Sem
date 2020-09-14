#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <set>
#include <algorithm> 

using std::cin;
using std::cout;
using std::set;
using std::vector;
using std::pair;

class WEGraph {
public:
	WEGraph() {}
	WEGraph(size_t vertex_count) {
		this->vertex_count = vertex_count;
	}
	~WEGraph() = default;
	void AddEdge(size_t head, size_t tail, size_t cost) {
		pair < size_t, pair<size_t, size_t> > 
			edge(cost, pair<size_t, size_t>(head, tail));
		graph.push_back(edge);
	}
	size_t Edges() const {
		return graph.size();
	}
	size_t VertexCount() const {
		return vertex_count;
	}
	size_t GetHead(size_t i) const {
		return graph[i].second.first;
	}
	size_t GetTail(size_t i) const {
		return graph[i].second.second;
	}
	size_t GetCost(size_t i) const {
		return graph[i].first;
	} 
	pair < size_t, pair<size_t, size_t> > GetEdge(size_t i) const {
		return graph[i];
	}
private:
	vector < pair < size_t, pair<size_t, size_t> > > graph;
	size_t vertex_count;
};

bool Comparator(const pair < size_t, pair<size_t, size_t> > a,
				const pair < size_t, pair<size_t, size_t> > b) {
	return a.first < b.first;
}

void KruskalFindMST(const WEGraph& graph, WEGraph& result) {

	vector < pair < size_t, pair<size_t, size_t> > > e;
	size_t n = graph.VertexCount();
	size_t m = graph.Edges();
	for (size_t i = 0; i < m; ++i) {
		e.push_back(graph.GetEdge(i));
	}
	sort(e.begin(), e.end());

	vector<size_t> color(n);
	for (size_t i = 0; i < n; ++i) {
		color[i] = i;
	}

	size_t head, tail, cost, old_color, new_color;
	for (size_t i = 0; i < m; ++i) {
		cost = e[i].first;
		tail = e[i].second.first;
		head = e[i].second.second;

		if (color[head] != color[tail]) {
			new_color = color[head];
			old_color = color[tail];
			result.AddEdge(head, tail, cost);

			for (size_t i = 0; i < n; ++i) {
				if (color[i] == old_color) {
					color[i] = new_color;
				}
			}
		}
	}
}


void ReadGraph(WEGraph& graph, size_t edges) {
	size_t head, tail, cost;
	for (size_t i = 0; i < edges; ++i) {
		cin >> head >> tail >> cost;
		graph.AddEdge(head - 1, tail - 1, cost);
	}
}

int main() {
	size_t n, m;
	cin >> n >> m;
	WEGraph graph(n);
	ReadGraph(graph, m);

	WEGraph result;
	KruskalFindMST(graph, result);

	size_t summ = 0;
	for (size_t i = 0; i < result.Edges(); ++i) {
		summ += result.GetCost(i);
	}

	cout << summ;


}