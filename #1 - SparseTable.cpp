#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using std::cin;
using std::cout;
using std::vector;
using std::min;
using std::pow;

void ReadVector(vector <int>& v, int n) {
	int c;
	for (int i = 0; i < n; ++i) {
		cin >> c;
		v.push_back(c);
	}
}

struct IndVal {
	IndVal() = default;
	IndVal(int ind, int val) {
		this->ind = ind;
		this->val = val;
	}
	int ind;
	int val;
};

void CreateSparseTable(vector <int>& v, vector <vector <IndVal> >& t) {
	t.assign(v.size(), vector<IndVal>());
	for (int i = 0; i < v.size(); ++i) {
		t[i].push_back(IndVal(i, v[i]));
	}

	for (int j = 1; j <= floor(log2(t.size())); ++j) {
		for (int i = 0; i < t.size(); ++i) {

			if (i + pow(2, j) > t.size()) {
				continue;
			}

			if (t[i][j - 1].val < t[i + pow(2, j - 1)][j - 1].val) {
				t[i].push_back(IndVal(t[i][j - 1].ind, t[i][j - 1].val));
			}
			else {
				t[i].push_back(IndVal(
					t[i + pow(2, j - 1)][j - 1].ind,
					t[i + pow(2, j - 1)][j - 1].val));
			}
		}
	}
}


IndVal RMQMin(vector <vector <IndVal> >& t, int head, int tale) {
	if (head == tale) {
		return t[head][0];
	}

	int len = floor(log2(tale - head + 1));
	if (t[head][len].val < t[tale - pow(2, len) + 1][len].val) {
		return t[head][len];
	}
	else {
		return t[tale - pow(2, len) + 1][len];
	}
}

void PrintVector(vector <int>& v) {
	for (int i = 0; i < v.size(); ++i) {
		cout << v[i] << '\n';
	}
}

int MinOnSegments(vector <vector <IndVal> >& sparse_table,
	int head, IndVal hole, int tale) {

	if (head > hole.ind - 1) {
		return RMQMin(sparse_table, hole.ind + 1, tale).val;
	}

	else if (hole.ind + 1 > tale) {
		return RMQMin(sparse_table, head, hole.ind - 1).val;
	}

	else if (RMQMin(sparse_table, head, hole.ind - 1).val <
		RMQMin(sparse_table, hole.ind + 1, tale).val) {

		return RMQMin(sparse_table, head, hole.ind - 1).val;
	}
	else {
		return RMQMin(sparse_table, hole.ind + 1, tale).val;
	}
}

void Print2OrderStat(vector < vector<IndVal> > sparse_table, int requests) {
	int head, tale;
	IndVal first_ord;

	for (int i = 0; i < requests; ++i) {
		cin >> head >> tale;
		head -= 1;
		tale -= 1;

		first_ord = RMQMin(sparse_table, head, tale);
		cout << MinOnSegments(sparse_table, head, first_ord, tale) << '\n';
	}
}

int main() {
	int n, m;
	cin >> n >> m;
	vector <int> sequence;
	ReadVector(sequence, n);

	vector < vector<IndVal> > sparse_table;
	CreateSparseTable(sequence, sparse_table);

	Print2OrderStat(sparse_table, m);
}