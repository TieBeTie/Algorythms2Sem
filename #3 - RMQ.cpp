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
		tree = vector <int>(4 * size, std::numeric_limits<int>::max());
		this->size = size;
	}
	SegmentTree(vector <int>& root) {
		tree = vector <int>(4 * root.size(), std::numeric_limits<int>::max());
		this->size = root.size();
		Build(root, 1, 0, size);
		updates.assign(size * 4, std::numeric_limits<int>::max());
	}
	~SegmentTree() = default;
	int MinInquiry(int l, int r) {
		return Min(1, 0, size, l, r);
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
			tree[node] = min(tree[node * 2], tree[node * 2 + 1]);
		}
	}

	void Push(int v) {
		if (updates[v] != std::numeric_limits<int>::max()) {
			tree[v * 2] = tree[v];
			tree[v * 2 + 1] = tree[v];
		}
		if (updates[v] != std::numeric_limits<int>::max()) {
			updates[v * 2] = updates[v];
			updates[v * 2 + 1] = updates[v];
			updates[v] = std::numeric_limits<int>::max();
		}
	}

	int Min(int node, int l, int r, int l_in, int r_in) {
		if (l_in >= r_in) {
			return std::numeric_limits<int>::max();
		}

		if (l == l_in && r == r_in) {
			return tree[node];
		}

		Push(node);
		int m = (l + r) / 2;
		return min(Min(node * 2, l, m, l_in, min(r_in, m)),
			Min(node * 2 + 1, m, r, max(l_in, m), r_in));
	}
	void Update(int node, int l, int r, int l_in, int r_in, int val) {
		if (l_in >= r_in) {
			return;
		}

		if (l == l_in && r == r_in) {
			tree[node] = val;
			updates[node] = val;
		}
		else {
			Push(node);
			int m = (l + r) / 2;
			Update(node * 2, l, m, l_in, min(m, r_in), val);
			Update(node * 2 + 1, m, r, max(l_in, m), r_in, val);
			tree[node] = min(tree[node * 2], tree[node * 2 + 1]);
		}
	}

	vector <int> tree;
	vector <int> updates;
	int size;
};

int main() {
	int len;
	cin >> len;

	vector <int> segment(len);
	int r, g, b;
	for (int i = 0; i < len; ++i) {
		cin >> r >> g >> b;
		segment[i] = r + g + b;
	}
	SegmentTree tree(segment);

	int art_processes;
	cin >> art_processes;

	int head, tale;
	for (int i = 0; i < art_processes; ++i) {
		cin >> head >> tale >> r >> g >> b;

		tree.UpdateTree(head, tale + 1, r + g + b);
		
		cin >> head >> tale;
		cout << tree.MinInquiry(head, tale + 1) << ' ';
	}
}