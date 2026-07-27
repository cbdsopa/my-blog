+++
date = '2026-07-27T14:46:04+08:00'
draft = false
title = 'Test Post'
categories = ["测试"]
tags = ["测试"]
+++

本页面用于测试功能。

# markdowm测试

$$\phi(n) = \sum_{i = 1}^{n} [\gcd(i, n) = 1]$$

你好？

时间复杂度为 $\mathcal{O}(n)$。我们考虑使用一个超级长的行来测试他同一段内文字的这个间距是否符合一个舒适的观看条件，这要求我们的文本需要特别特别的长来帮助我们完成这一点。

this is a question.

what i wanna tell you is that...

# H1标题

## H2标题

### H3标题

#### H4标题

##### H5标题

###### H6标题



# 图片测试

![test图片](https://cdn.luogu.com.cn/upload/image_hosting/4epse4bs.png)

[旧博客链接](https://www.cnblogs.com/cbdsopa)

# 代码测试

```cpp
#include <bits/stdc++.h>
using namespace std;
using ull = unsigned long long;
using ll = long long;
using ldb = long double;
using db = double;

template <typename T>
auto makeVector(int n, T val) {
	return vector<T>(n, val);
}
template <typename... Args>
auto makeVector(int n, Args... args) {
	auto val = makeVector(args...);
	return vector<decltype(val)>(n, val);
}

const int mod = 998244353;
int inc(int x, int y, int p = mod) { return (x += y) >= p ? x - p : x; }
int dec(int x, int y, int p = mod) { return (x -= y) < 0 ? x + p : x; }
int mul(int x, int y, int p = mod) { return 1ll * x * y % mod; }
int qpow(int a, ll b, int p = mod) {
    int res = 1;
    while(b) {
        if(b & 1) res = mul(res, a, p);
        a = mul(a, a, p);
        b >>= 1;
    }
    return res;
}


struct Matrix {
	int n;
	decltype(makeVector(0, 0, (int)0)) c;
	Matrix(int n, int val = 0) : n(n), c(makeVector(n, n, val)) {

	}
	void clearE() {
		for (int i = 0; i < n; ++i) {
			c[i][i] = 1;
		}
	}

	Matrix operator * (const Matrix &rhs) const {
		Matrix res(n);
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < n; ++k) {
				if (!c[i][k]) continue;
				for (int j = 0; j < n; ++j) {
					res.c[i][j] = inc(
						res.c[i][j],
						mul(c[i][k], rhs.c[k][j])
					);
				}
			}
		}
		return res;
	}
};

vector<int> operator * (const vector<int> &v, const Matrix &m) {
	int n = v.size();
	vector<int> res(n);
	for (int k = 0; k < n; ++k) {
		if (!v[k]) continue;
		for (int j = 0; j < n; ++j) {
			res[j] = inc(
				res[j],
				mul(v[k], m.c[k][j])
			);
		}
	}
	return res;
}


vector<int> bm(const vector<int> &s) {
	vector<int> C = {1}, B = {1};
	int L = 0, m = 1, b = 1;
	for (int i = 0; i < s.size(); ++i) {
		int d = 0;
		for (int j = 0; j <= L; ++j) {
			d = inc(d, mul(s[i - j], C[j]));
		}
		if (d == 0) {
			++m;
			continue;
		}
		vector D = C;
		int c = mul(mod - d, qpow(b, mod - 2));
		if (C.size() < B.size() + m) C.resize(B.size() + m);
		for (int j = 0; j < B.size(); ++j) {
			C[j + m] = inc(C[j + m], mul(c, B[j]));
		}
		if (2 * L <= i) {
			L = i - L + 1;
			B = D;
			m = 1;
			b = d;
		} else {
			++m;
		}
	}
	C.resize(L + 1);
	return C;
}

vector<int> polyMul(const vector<int> &a, const vector<int> &b, const vector<int> &bm_seq) {
	int n = a.size() + b.size() - 1;
	int L = bm_seq.size() - 1;
	vector<int> res(n);
	for (int i = 0; i < a.size(); ++i) {
		for (int j = 0; j < b.size(); ++j) {
			res[i + j] = inc(res[i + j], mul(a[i], b[j])); 
		}
	}
	for (int i = res.size() - 1; i >= L; --i) {
		if (!res[i]) continue;
		for (int j = 1; j <= L; ++j) {
			res[i - j] = dec(res[i - j], mul(bm_seq[j], res[i]));
		}
	}
	res.resize(min<int>(L, res.size()));
	return res;
}
vector<int> polyQpow(ll b, const vector<int> &bm_seq) {
	vector<int> res = {1}, a = {0, 1};
	while (b) {
		if (b & 1) res = polyMul(res, a, bm_seq);
		a = polyMul(a, a, bm_seq);
		b >>= 1;
	}
	return res;
}


void solve() {
	int n, m;
	cin >> n >> m;
	if (n == 1) {
		cout << mul(3, qpow(2, m - 1) ) << "\n";
		return;
	}
	int k = 1 << (n - 2);
	Matrix trans(k);

	for (int s = 0; s < k; ++s) {
		vector<int> sa(n);
		sa[0] = 0;
		sa[1] = 1;
		for (int i = 2; i < n; ++i) {
			sa[i] = (s >> i - 2 & 1) ? sa[i - 2] : 3 - sa[i - 1] - sa[i - 2];
		}

		for (int x = 0; x < 3; ++x)
		for (int y = 0; y < 3; ++y)
		for (int t = 0; t < k; ++t) {
			if (x == y) continue;
			vector<int> ta(n);
			ta[0] = x;
			ta[1] = y;
			for (int i = 2; i < n; ++i) {
				ta[i] = (t >> i - 2 & 1) ? ta[i - 2] : 3 - ta[i - 1] - ta[i - 2];
			}

			bool flag = 1;
			for (int i = 0; i < n; ++i) {
				if (ta[i] == sa[i]) { 
					flag = 0;
					break;
				}
			}
			trans.c[s][t] += flag;
		}
	}

	int lim = min(128 * 2, m);
	vector<int> ans(k, 6);
	vector<int> seq(lim);
	for (int i = 0; i < lim; ++i) {
		for (int j = 0; j < k; ++j) seq[i] = inc(seq[i], ans[j]);
		
		if (i == lim - 1)  break;
		
		ans = ans * trans;
	}

	if (m - 1 < lim) {
		cout << seq[m - 1] << "\n";
		return;
	}
	vector<int> bm_seq = bm(seq);
	vector<int> R = polyQpow(m - 1, bm_seq);
	int result = 0;
	for (int i = 0; i < R.size(); ++i) {
		result = inc(result, mul(R[i], seq[i]));
	}
	cout << result << "\n";
}
int main(){
	ios::sync_with_stdio(false); cin.tie(0);
#ifdef LOCAL
	freopen("a.in", "r", stdin);
	freopen("a.out", "w", stdout);
#endif
	int T = 1; 
	// cin >> T;
	while (T--) {
	    solve();
	}
	return 0;
}
```
