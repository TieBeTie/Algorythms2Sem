#include <iostream>
#include <vector>
#include <limits>
#include <queue>

using std::cin;
using std::cout;
using std::priority_queue;
using std::vector;
using std::pair;

struct Edge {
	Edge() {}
	Edge(size_t from, size_t to, double cost, double commision) {
		this->from = from;
		this->to = to;
		this->cost = cost;
		this->commision = commision;
	}
	size_t from;
	size_t to;
	double cost;
	double commision;
};

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

void PrintAnswer(bool answer) {
	if (answer) {
		cout << "YES";
	} else {
		cout << "NO";
	}
}

void ReadNums(size_t& s1, size_t& s2, size_t& s3, double& d1) {
	cin >> s1 >> s2 >> s3 >> d1;
}

void InitExchangers(vector <Edge>& exchangers, size_t& exchangers_count) {
	size_t from, to;
	double from_cost, from_commision, to_cost, to_commision;
	for (size_t i = 0; i < exchangers_count; ++i) {
		cin >> from >> to >> from_cost >> from_commision >> to_cost >> to_commision;
		exchangers.push_back(Edge(from - 1, to - 1, from_cost, from_commision));
		exchangers.push_back(Edge(to - 1, from - 1, to_cost, to_commision));
	}
}

bool Floyd_Bellman(const vector <Edge>& edges, size_t vertex_count, 
				   size_t start, double start_money) {
	vector <double> d (vertex_count, 0);
	d[start] = start_money;
	size_t edges_count = edges.size();
	bool relaxation;
	for (size_t i = 0; i < vertex_count; ++i) {
		relaxation = false;
		for (size_t j = 0; j < edges_count; ++j) {
			if (d[edges[j].to] < (d[edges[j].from] - edges[j].commision) * edges[j].cost) {
				(d[edges[j].to] = (d[edges[j].from] - edges[j].commision) * edges[j].cost);
				relaxation = true;
			}
		}
	}
	return relaxation;
}
int main() {
	InOutOptimizaton();
	size_t currencies, exchangers_count, rick_currency;
	double rick_money;
	ReadNums(currencies, exchangers_count, rick_currency, rick_money);
	rick_currency--;
	vector <Edge> exchangers;
	InitExchangers(exchangers, exchangers_count);
	

	bool answer;
	answer = Floyd_Bellman(exchangers, currencies, rick_currency, rick_money);
	PrintAnswer(answer);
}
