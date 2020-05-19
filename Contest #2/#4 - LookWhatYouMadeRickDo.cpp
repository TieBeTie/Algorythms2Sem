#include <iostream>
#include <vector>
#include <limits>
#include <bitset>
#include <string>

using std::cin;
using std::cout;
using std::vector;
using std::pair;
using std::bitset;
using std::string;

struct Node {
	size_t vertex;
	double chance;
	Node(size_t vertex, double chance) {
		this->vertex = vertex;
		this->chance = chance;
	}
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}
void CreateMultiUniverseGraph(int size, vector <bitset<1005>>& multi_universe) {
	size_t n = size;
	string mask;
	for (size_t i = 0; i < n; ++i) {
		cin >> mask;
		multi_universe.push_back(bitset<1005>(mask));
	}
}

void FloydWarshall(vector <bitset<1005>>& graph) {
	for (int i = 0; i < graph.size(); ++i) {
		for (int j = 0; j < graph.size(); ++j) {
			if (graph[j][i]) {
				graph[j] = graph[j] | graph[graph.size() - i - 1];
			}
		}
	}
}

void COutGraph(vector <bitset<1005>> graph) {
	for (int i = 0; i < graph.size(); ++i) {
		for (int j = graph.size() - 1; j >= 0; --j) {
			cout << graph[i][j];
		}
		cout << "\n";
	}
}

int main() {
	InOutOptimizaton();
	size_t count_universes;
	cin >> count_universes;	
	vector <bitset<1005>> multi_universe;
	CreateMultiUniverseGraph(count_universes, multi_universe);
	FloydWarshall(multi_universe);
	COutGraph(multi_universe);
}
