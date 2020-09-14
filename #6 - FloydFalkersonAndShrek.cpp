#include <iostream>
#include <vector>
#include <limits>
#include <queue>

using std::cin;
using std::cout;
using std::vector;
using std::pair;

void InOutOptimizaton() {
	std::ios::sync_with_stdio(0);
	cin.tie(0);
}

void StoerWagner(vector < vector <int> >& g, int n, 
	vector <int>& best_cut, int& best_cost) {

	vector < vector <int> > compression(n);
	vector <int> cost(n);
	vector <bool> exist(n);
	vector <bool> in_a(n);
	int max_cost, pre;

	for (int i = 0; i < n; ++i) {
		exist[i] = true;
	}

	for (int i = 0; i < n; ++i) {
		compression[i].assign(1, i);
	}

	for (int j = 0; j < n - 1; ++j) {

		for (int i = 0; i < n; ++i) {
			in_a[i] = false;
		}


		for (int i = 0; i < n; ++i) {
			cost[i] = 0;
		}

		for (int x = 0; x < n - j; ++x) {
			max_cost = -1;

			for (int i = 0; i < n; ++i) {

				if (exist[i] && !in_a[i] && 
					(max_cost == -1 
						|| cost[i] > cost[max_cost])) {
					max_cost = i;
				}

			}

			if (x == n - j - 1) {

				if (cost[max_cost] < best_cost) {
					best_cost = cost[max_cost]; 
					best_cut = compression[max_cost];
				}

				compression[pre].insert
				(compression[pre].end(), 
					compression[max_cost].begin(), 
					compression[max_cost].end());


				for (int i = 0; i < n; ++i) {
					g[i][pre] += g[max_cost][i];
					g[pre][i] = g[i][pre];
				}

				exist[max_cost] = false;
			} else {
				in_a[max_cost] = true;

				for (int i = 0; i < n; ++i) {
					cost[i] += g[max_cost][i];
				}

				pre = max_cost;
			}
		}
	}
}

void ReadMatrix(vector < vector <int> >& matrix, int n) {
	char c;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cin >> c;
			if (c == '0') {
				matrix[i].push_back(0);
			} else {
				matrix[i].push_back(1);
			}
		}
	}
}

void PrintLobes(vector <int>& best_cut, int vertex_count) {
	vector <bool> lobe(vertex_count, false);
	for (int i = 0; i < best_cut.size(); ++i) {
		lobe[best_cut[i]] = true;
	}
	for (int i = 0; i < vertex_count; ++i) {
		if (!lobe[i]) {
			cout << i + 1 << ' ';
		}
	}
	cout << '\n';
	for (int i = 0; i < vertex_count; ++i) {
		if (lobe[i]) {
			cout << i + 1 << ' ';
		}
	}
}

int main() {
	int n;
	cin >> n;
	vector < vector <int> > graph (n);

	ReadMatrix(graph, n);

	vector <int> best_cut;
	int best_cost = std::numeric_limits <int>::max();
	StoerWagner(graph, n, best_cut, best_cost);

	PrintLobes(best_cut, n);
}