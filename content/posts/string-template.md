+++
date = '2026-08-17T02:44:17+08:00'
draft = false
title = '算法竞赛字符串性质与模板'
categories = ["字符串"]
tags = ["Hash", "KMP", "Z函数", "Manacher", "AC自动机", "后缀数组", "后缀自动机", "回文自动机", "最小表示法", "Lyndon分解", "Runs"]
+++

# 字符串

## hash

### 数字工具

小质数：101, 103, 1009, 10007, 10009, 100003, 1000003, 1145141

大质数：998244353, 1000000007, 1000000009

## kmp 和 border

### 性质

fail 树子树里的位置为该 border 的出现位置。

> **性质1.** border 的 border 还是 border，并且每次取最大的 border 可以遍历原串所有 border。

> **性质2. （弱周期引理）** 对于周期 $p, q$，如果 $p + q \le n$，则 $gcd(p, q)$ 也是周期。

> **性质3.** 对于字符串 $s, t$，$s$ 为 $t$ 的前缀，$t$ 有周期 $a$，$s$ 有周期 $b$，满足 $a \le |s|$ 且 $b \mid a$，则 $b$ 也为 $t$ 的周期。

> **性质4.** 对于字符串 $s, t$，若 $2|s| \ge |t|$，则 $s$ 在 $t$ 上的匹配位置形成等差数列。

> **性质5.** 一个字符串长为 $n$，其长度 $\ge \lceil \frac{n}{2} \rceil$ 的 border 构成等差数列。具体的，如果最长 border 为 $n - p$，那么公差为 $p$。

> **性质6.** 一个字符串其 border 形成 $O(log n)$ 个等差数列。

### 代码

```cpp
// 均摊 O(n) 版
void kmp(const string &s) {
	int n = s.size() - 1;
	vector<int> fail(n + 1);
	for (int i = 2, j = 0; i <= n; ++i) {
		while (j && s[j + 1] != s[i]) j = fail[j];
		if (s[j + 1] == s[i]) ++j;
		fail[i] = j;
	}
}
// 单次 O(log n) 版
void kmp(const string &s) {
	int n = s.size() - 1;
	vector<int> fail(n + 1);
	for (int i = 2, j = 0; i <= n; ++i) {
		while(j && s[j + 1] != s[i]){
			if(j + 1 < i && fail[j + 1] * 2 > j + 1) 
				j = (j - 1) % (j - fail[j]) + 1;
			else j = fail[j];
		}
		if (s[j + 1] == s[i]) ++j;
		fail[i] = j;
	}
}
```

## z 函数（扩展kmp）

### 代码

```cpp
void z_func(const string &s) {
	int n = s.size() - 1;
	vector<int> z(n + 1);
	z[1] = n;
	for (int i = 2, l = 0, r = -1; i <= n; ++i) {
		int k = 0;
		if (i <= r) k = min(z[i - l + 1], r - i + 1);
		// s[i + k] 改成 t[i + k]，可以求 t 的后缀与 s 的LCP
		while (i + k <= n && s[i + k] == s[1 + k])
			++k;
		z[i] = k;
		if (i + k - 1 > r) {
			l = i;
			r = i + k - 1;
		}
	}
}
```

## manacher

### 代码

```cpp
void manacher(const string &s) {
	int n = s.size() - 1;
	vector<int> d1(n + 1, 1);
	// len = 2 * d1[i] - 1
	for (int i = 1, l = 0, r = -1; i <= n; ++i) {
		int k = 1;
		if (i <= r) k = min(d1[l + r - i], r - i + 1);
		while (i + k <= n && i - k >= 1 && s[i + k] == s[i - k]) 
			++k;
		d1[i] = k--;
		if (i + k > r) {
			l = i - k;
			r = i + k; 
		}
	}
	vector<int> d2(n + 1, 0);
	// len = 2 * d2[i]
	for (int i = 1, l = 0, r = -1; i <= n; ++i) {
		int k = 0;
		if (i <= r) k = min(d2[l + r - i + 1], r - i + 1);
		while (i + k <= n && i - k - 1 >= 1 && s[i + k] == s[i - k - 1]) 
			++k;
		d2[i] = k--;
		if (i + k > r) {
			l = i - k - 1;
			r = i + k; 
		}
	}
}
// 统一处理版
void manacher(const string &t) {
	int n = t.size() - 1;
	string s = " ";
	for (int i = 1; i <= n; ++i) {
		s += '#';
		s += t[i]; 
	}
	s += '#';
	n = 2 * n + 1;
	vector<int> d(n + 1, 1);
	// len = d[i]
	for (int i = 1, l = 0, r = -1; i <= n; ++i) {
		int k = 1;
		if (i <= r) k = min(d[l + r - i], r - i + 1);
		while (i + k <= n && i - k >= 1 && s[i + k] == s[i - k]) 
			++k;
		d[i] = --k;
		if (i + k > r) {
			l = i - k;
			r = i + k; 
		}
	}
}
```

## AC自动机

### 性质

fail 指向 border，up 指向 fail 链中上一个真实具有字符串的节点。

> **性质1.** AC自动机为 DAG，fail 形成树，有时可以采用拓扑排序或者子树 dfs 优化。

> **性质2.** 字符串有总长限制时，不同长度的真实串只有 $O(\sqrt S)$ 个，遍历 up 链时可利用。

### 代码

```cpp
struct ACAM {
	#define fail(x) t[x].fail
	#define trans(x) t[x].trans 
	struct node {
		bool end_pos = 0;
		int fail = 0;
		int up = 0;
		int trans[26] = {};
		int len = 0;
	};
	vector<node> t;
	ACAM(int n = 0) {
		t.reserve(n + 1);
		newnode();
	}
	int newnode() {
		t.push_back(node());
		return t.size() - 1;
	}
	int ins(const string &s) {
		int p = 0;
		for (char cha : s) {
			int c = cha - 'a';
			if (!t[p].trans[c]) {
				t[p].trans[c] = newnode();
				t[t[p].trans[c]].len = t[p].len + 1;
			} 
			p = t[p].trans[c];
		}
		t[p].end_pos = 1;
		return p;
	}
	void build() {
		queue<int> q;
		for (int i = 0; i < 26; ++i) {
			if (t[0].trans[i]) {
				q.push(t[0].trans[i]);
			}
		}
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (int i = 0; i < 26; ++i) {
				int &v = t[u].trans[i];
				if (v) {
					t[v].fail = t[t[u].fail].trans[i];
					if (t[t[v].fail].end_pos) {
						t[v].up = t[v].fail;
					} else {
						t[v].up = t[t[v].fail].up;
					}
					q.push(v);
				} else {
					v = t[t[u].fail].trans[i];
				}
			}
		}
	}
};
```

## 后缀数组

### 性质

$high_i$ 为 $sa_i$ 和 $sa_{i - 1}$ 的 LCP。

> **性质1.** 按 high 最小值分治得到的本质是简化后缀树。

> **性质2.** $high_{rk_i} \ge high_{rk_{i-1}} - 1$

> **性质3.** 两个后缀的 LCP 为对应 high 的区间 min。

### 应用

#### 求子串出现次数

对于子串 $s[l,r]$ ，求出前缀为该子串的字典序最小最大的后缀的排名 $x,y$，答案 $y-x+1$。

等价求有多少个后缀与 $s[l,n]$ 的 LCP $\ge r-l+1$。

#### 比较子串大小

对于 $A=[a,b],B=[c,d]$，若 $LCP(a,c)\ge \min\{|A|,|B|\}$ ，那么比较 $|A|,|B|$ 大小，否则比较 $rk[a],rk[c]$。

#### 本质不同子串个数

$\frac{n(n+1)}{2}-\sum_{i=2}^nhigh[i]$

计算所有新增的 LCP 即可。

#### 最小表示法

求一个串循环同构中字典序最小的那个。

倍长原串，其同构转为后缀后字典序相对大小不改变，于是可得。

为什么相对大小不变？因为发现多出来的部分是个前缀，你原来比它字典序小，那之后仍然是小的。特殊情况是两个相同的循环同构，但是由于两者相同求出来的结果不变。

#### 最长公共子串

等价于两（或多个）个串拼起来求最大的后缀 LCP 。两个串间添加不同的不在字符集中的元素作为分隔符。

我们考虑二分答案后，检验是否可达。具体就是找到一段连续段 $[l,r]$ 满足 $x\in[l,r],high[x]\ge ans$ ，且对于 $x\in[l-1,r]$ ，使得每个串都有至少一个元素在里面。

### 代码

```cpp
struct SA{
    int n;
	vector<int> rk, id, cnt, sz, lark, sa;
	vector<int> high, lg2;
	vector<vector<pair<int, int>>> mn;
    string s;

	SA(const string &c, int _n) : n(_n), s(c), 
		rk(_n + 1), id(_n + 1), cnt(max(128, _n) + 1), sz(_n + 1), 
		lark(_n + 1), sa(_n + 1), high(_n + 1), lg2(_n + 1),
		mn(make_vector(_n + 1, __lg(max(1, _n)) + 2, pair<int, int>()))
	{
        int m = max(128, n);
        for (int i = 0; i <= m; ++i) cnt[i] = 0;
        for (int i = 1; i <= n; ++i) ++cnt[rk[i] = c[i] ], s[i] = c[i];
        for (int i = 1; i <= m; ++i) cnt[i] += cnt[i - 1];
        for (int i = n; i >= 1; --i) sa[cnt[rk[i] ]--] = i;
        for (int w = 1; w < n; w <<= 1){
            int p = 0;
            for (int i = n; i > n - w; --i) id[++p] = i;
            for (int i = 1; i <= n; ++i) 
                if (sa[i] > w) id[++p] = sa[i] - w;
            for (int i = 0; i <= m; ++i) cnt[i] = 0;
            for (int i = 1; i <= n; ++i) ++cnt[rk[id[i] ] ];
            for (int i = 1; i <= m; ++i) cnt[i] += cnt[i - 1];
            for (int i = n; i >= 1; --i) sa[cnt[rk[id[i] ] ]--] = id[i];
            for (int i = 1; i <= n; ++i) lark[i] = rk[i];
            m = 0;
            for (int i = 1; i <= n; ++i){
				int x = sa[i] + w <= n ? lark[sa[i] + w] : 0;
				int y = sa[i - 1] + w <= n ? lark[sa[i - 1] + w] : 0;
				if (
					i == 1 ||
					lark[sa[i]] != lark[sa[i - 1]] ||
					x != y
				) ++m;
                rk[sa[i] ] = m;
            }
            if (m == n) break;
        }
        get_high();
        init_st();
	}
    void get_high(){
        int t = 0;
        for(int i = 1; i <= n; ++i){
			if (rk[i] == 1) continue;
            if (t) --t;
            while (
				i + t <= n &&
       			sa[rk[i] - 1] + t <= n &&
	   			s[i + t] == s[sa[rk[i] - 1] + t]
			) ++t;
            high[rk[i] ] = t;
        }
    }
    void init_st(){
        for (int i = 1; i <= n; ++i) mn[i][0] = {high[i], i};
        lg2[1] = 0;
        for (int i = 2; i <= n; ++i){
            lg2[i] = lg2[i >> 1] + 1;
        }
        for (int j = 1; (1 << j) <= n; ++j){
            for (int i = 1; i + (1 << j) - 1 <= n; ++i){
                mn[i][j] = min(mn[i][j - 1], mn[i + (1 << j - 1)][j - 1]);
            }
        }
    }
    pair<int, int> query(int l, int r){
        if (l > r) return {-1, -1};
        int s = lg2[r - l + 1];
        return min(mn[l][s], mn[r - (1 << s) + 1][s]);
    }
};
```

## 后缀自动机

### 性质

endpos 集合相同的点，我们认为他们等价。后缀自动机上一个点表示一个上述等价关系的等价类。定义后缀链接 link 指向 endpos 包含当前点的长度最长的等价类。

> **性质1.** 等价类内部长度连续。

> **性质2.** 两个子串的 endpos 要么没有交集，要么满足包含关系。

> **性质3.** $len_{min}(x) = len_{max}(link(x)) + 1$

> **性质4.** 点数最多 $2n - 1$，边数最多 $3n - 4$。

> **性质5.** 只有 link 树的叶子的 endpos 个数为 1，其他为 0。

如果需要统计等价类出现次数，只需统计子树内 endpos 个数即可。也就是叶子个数。

> **性质6.** link树内，一个点的 len 等于子树的所有 endpos 最长公共后缀长度。

> **性质7.** 反串的 SAM 为其正串的简化后缀树。

倍长之后建 SAM，一直走最小的字符，走 $|s|$ 个可以得到 s 的最小表示法。

对于广义后缀自动机，在我们每次设置的 lst 的位置维护信息即可。

### 构造

考虑怎么求出 $link(cur)$。

令 $p = lst$， 在 $p$ 的 link 链上找到第一个 $p'$ 满足 $q = \delta(p', c) \neq \emptyset$。尝试使用 $q$ 拆出新节点。 

如果没有 $q$ 直接连源点了。如果存在 $q$，查看是否满足性质3，满足的话直接 $link(cur) = q$。否则考虑把 $q$ 拆了，把出现位置增加的那部分后缀拆出来成 $q_2$，也就是 p 那边过来的部分，原本部分为 $q_1$，然后让原本的 $q_1$ 和 $cur$ 的 link 都连向 $q_2$ 即可。

### 代码

```cpp
struct SAM {
	struct node {
		int len = 0;
		int link = 0;
		int trans[26] = {};
	};
	int lst = 1;
	vector<node> t;
	SAM(int n) {
		t.reserve(2 * n + 2);
		newnode();
		newnode();
	}
	int newnode() {
		t.push_back(node());
		return t.size() - 1;
	}
	void ins(int c) {
		int cur = newnode(), p = lst;
		t[cur].len = t[lst].len + 1;
		lst = cur;
		while (p && !t[p].trans[c]) {
			t[p].trans[c] = cur;
			p = t[p].link;
		}
		if (!p) t[cur].link = 1;
		else {
			int q1 = t[p].trans[c];
			if (t[p].len + 1 == t[q1].len) t[cur].link = q1;
			else {
				int q2 = newnode();
				t[q2] = t[q1];
				t[q2].len = t[p].len + 1;
				t[q1].link = t[cur].link = q2;
				while (p && t[p].trans[c] == q1) {
					t[p].trans[c] = q2;
					p = t[p].link;
				} 
			}
		}
	}
};

struct GSAM {
	struct node {
		int len = 0;
		int link = 0;
		int trans[26] = {};
	};
	int lst = 1;
	vector<node> t;
	GSAM(int n) {
		t.reserve(2 * n + 2);
		newnode();
		newnode();
	}
	int newnode() {
		t.push_back(node());
		return t.size() - 1;
	}
	void next_str() {
		lst = 1;
	}
	void ins(int c) {
		int p = lst;
		if (t[p].trans[c]) {
			int q1 = t[p].trans[c];
			if(t[q1].len == t[p].len + 1) { 
				lst = q1;
				return;
			}
			int q2 = newnode();
			t[q2] = t[q1];
			t[q2].len = t[p].len + 1;
			t[q1].link = q2;
			while (p && t[p].trans[c] == q1) {
				t[p].trans[c] = q2;
				p = t[p].link;
			}
			lst = q2;
			return;
		}
		int cur = newnode();
		t[cur].len = t[lst].len + 1;
		lst = cur;
		while (p && !t[p].trans[c]) {
			t[p].trans[c] = cur;
			p = t[p].link;
		}
		if (!p) t[cur].link = 1;
		else {
			int q1 = t[p].trans[c];
			if (t[p].len + 1 == t[q1].len) t[cur].link = q1;
			else {
				int q2 = newnode();
				t[q2] = t[q1];
				t[q2].len = t[p].len + 1;
				t[q1].link = t[cur].link = q2;
				while (p && t[p].trans[c] == q1) {
					t[p].trans[c] = q2;
					p = t[p].link;
				} 
			}
		}
	}
};
```

## 回文自动机

### 性质

fail 是最长回文后缀。

可以注意到回文自动机左插和右插其实是一样的，因为回文串的最长回文后缀也是前缀，只需要额外维护 last_l 即可。

子串出现次数是 fail 链上访问次数和，可以拓扑排序一下。

> **性质1.** 回文串的 border 是回文串。

> **性质2.** 一个长为 $n$ 的字符串具有长度 $\ge \lceil \frac{n}{2} \rceil$ 的 border，则其为回文串。

> **性质3.** 回文串 $x$ ，$y$ 为 $x$ 的最长回文真后缀，$z$ 为 $y$ 的最长回文真后缀， $u,v$ 满足 $x=uy,y=vz$，则：

1.$|u|\ge |v|$。

2.若 $|u|\gt|v|$ ，则 $|u|\gt|z|$ 。

3.若 $|u|=|v|$，则 $u=v$。

> **性质4.** 字符串 $s$ 的所有回文后缀按长度排序后可划分为 $O(\log |s|)$ 组等差数列。

> **性质5.** 一个长为 $n$ 的字符串具有长度 $\ge \lceil \frac{n}{2} \rceil$ 的 border，则该 border 的只会作为前缀和后缀出现这两次。

优化回文划分的 $O(n ^ 2)$ 的 dp 时可以采用，只需要添加每个等差数列的新的值即可。

### 代码

```cpp
struct PAM {
    int lst = 0;
    string s;
    struct node {
        int len = 0;
        int fail = 0;
        int trans[26] = {};
		int dif = 0;
		int slink = 0;
    };
    vector<node> t;
    PAM(int n = 0) {
        t.reserve(n + 2);
        newnode(0);
        newnode(-1);
        t[1].fail = t[0].fail = 1;
    }
    int newnode(int L = 0) {
        t.push_back(node());
        t.back().len = L;
        return t.size() - 1;
    }
    void ins(char ch) {
        s += ch;
        int c = ch - 'a';
        int tip = s.size() - 1;
        auto extend = [&](int x) {
            while (tip - t[x].len - 1 < 0 || s[tip] != s[tip - t[x].len - 1]) 
                x = t[x].fail;
            return x;
        };
        int x = extend(lst);
        if (!t[x].trans[c]) {
            int y = newnode(t[x].len + 2);
            t[y].fail = t[extend(t[x].fail)].trans[c];
            t[x].trans[c] = y;
			t[y].dif = t[y].len - t[t[y].fail].len;
			if (t[y].dif == t[t[y].fail].dif) {
				t[y].slink = t[t[y].fail].slink;
			} else {
				t[y].slink = t[y].fail;
			}
        }
        lst = t[x].trans[c];
    }
};
```

## 最小表示法

```cpp
int min_cyclic_string() {
	int k = 0, i = 0, j = 1;
	while (k < n && i < n && j < n) {
		if (s[(i + k) % n] == s[(j + k) % n]) {
			k++;
		} else {
			s[(i + k) % n] > s[(j + k) % n] ? i = i + k + 1 : j = j + k + 1;
			if (i == j) i++;
			k = 0;
		}
	}
	return min(i, j);
}
```

## Lyndon 分解

### 性质

Lyndon 串指自己比所有后缀都小的串。

Lyndon 分解指将串分为多个 Lyndon 串，且字典序非严格单调降。不嫌麻烦可以 SA 每次取最小后缀分解得到。

动态维护的话，可以求出每个前缀的最小后缀，最小后缀为 Lyndon 分解中最后一个 Lyndon 串。

> **性质1.** Lyndon 分解必定存在且唯一。

> **性质2.** Lyndon 串没有 border。也就是说，Lyndon 串不应该具有周期。

> **性质3.** 对于 Lyndon 串 $s = ab$，有 $a < b$。

> **性质4.** 对于 Lyndon 串 $b$ 和另一个串 $a$, 有 $a < b \Leftrightarrow ab < b$。

> **性质5.** Lyndon 串可以被分为小 Lyndon 串按字典序严格单调增拼接。

> **性质6.** 若字符串 $s$ 和字符 $\overline x$ 满足 $s\overline x$ 是某个 Lyndon 串的前缀，则对于 $\overline y > \overline x$，$s\overline y$ 是 Ly。

如果一个串 $t=w^k\overline{w}$，其中 $w$ 为 Lyndon 串，$\overline{w}$ 是 $w$ 的严格前缀，那么称 $t$ 为近似 Lyndon 串。

### Duval 算法

duval 算法维护这个近似 Lyndon 串：

当新字符相匹配上了 $w$ 的下一位，那么周期不变；

当新字符比 $w$ 的下一位大，那么根据 `性质6` 形成 Lyndon 串，该 Lyndon 相较于 $w$ 字典序更大，不满足要求，需要把前面的 Lyndon 全部根据 `性质4` 合并；

当新字符比 $w$ 的下一位小，根据定义的字典序非严格单调降，前面的 $w$ 全部变成 Lyndon 分解里固定的 Lyndon 串。

### 代码

```cpp
vector<string> duval(string const &s) {
	int n = s.size(), i = 0;
	vector<string> factor;
	vector<int> min_suf(n);
	while (i < n) {
		min_suf[i] = i;
		int j = i + 1, k = i;
		while (j < n && s[k] <= s[j]) {
			if (s[k] < s[j])
				min_suf[j] = k = i;
			else
				min_suf[j] = min_suf[k] + j - k, ++k;
			++j;
		}
		while (i <= k) {
			factor.push_back(s.substr(i, j - k));
			i += j - k;
		}
	}
	return factor;
}
```

## Runs

