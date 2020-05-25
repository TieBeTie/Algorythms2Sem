#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using std::cin;
using std::cout;
using std::vector;
using std::min;
using std::max;
using std::pow;

class SegmentTree {
public:
	SegmentTree(int size) {
		tree = vector <int>(4 * size);
		this->size = size;
	}
	SegmentTree(vector <int>& root) {
		tree = vector <int>(4 * root.size(), -1);
		this->size = root.size();
		Build(root, 1, 0, size);
		updates.assign(size * 4, -1);
	}

	int MaxInquiry(int l, int r) {
		return Max(1, 0, size, l, r);
	}
	
	void UpdateTree(int l, int r, int val) {
		Update(1, 0, size, l, r, val);
	}

private:
	void Build(vector <int>& root, int node, int l, int r) {
		if (l + 1 == r) {
			tree[node] = root[l];
		}
		else {
			int m = (l + r) / 2;
			Build(root, node * 2, l, m);
			Build(root, node * 2 + 1, m, r);
			tree[node] = max(tree[node * 2], tree[node * 2 + 1]);
		}
	}

	void Push(int v) {
		if (updates[v] != -1) {
			tree[v * 2] = tree[v];
			tree[v * 2 + 1] = tree[v];
		}
		if (updates[v] != -1) {
			updates[v * 2] = updates[v];
			updates[v * 2 + 1] = updates[v];
			updates[v] = -1;
		}
	}

	int Max(int node, int l, int r, int l_in, int r_in) {
		if (l_in >= r_in) {
			return -1;
		}

		if (l == l_in && r == r_in) {
			return tree[node];
		}

		Push(node);
		int m = (l + r) / 2;
		return max(Max(node * 2, l, m, l_in, min(r_in, m)),
			Max(node * 2 + 1, m, r, max(l_in, m), r_in));
	}
	
	void Update(int node, int l, int r, int l_in, int r_in, int val) {
		if (l_in >= r_in) {
			return;
		}
		
		if (l == l_in && r == r_in) {
			tree[node] += val;
			updates[node] = tree[node];
		}
		else {
			Push(node);
			int m = (l + r) / 2;
			Update(node * 2, l, m, l_in, min(m, r_in), val);
			Update(node * 2 + 1, m, r, max(l_in, m), r_in, val);
			tree[node] = max(tree[node * 2], tree[node * 2 + 1]);
		}
	}

	vector <int> tree;
	vector <int> updates;
	int size;
};

int main() {
	int stations;
	cin >> stations;

	vector <int> people_count(stations - 1, 0);
	int tickets;
	for (int i = 0; i < stations - 1; ++i) {
		cin >> tickets;
		people_count[i] = tickets;
	}
	SegmentTree tree(people_count);

	int train_capacity;
	cin >> train_capacity;

	int inquiry_count;
	cin >> inquiry_count;

	int head, tale, people;
	for (int i = 0; i < inquiry_count; ++i) {
		cin >> head >> tale >> people;

		if (train_capacity < tree.MaxInquiry(head, tale) + people) {
			cout << i << ' ';
		}
		else {
			tree.UpdateTree(head, tale, people);
		}
	}
}