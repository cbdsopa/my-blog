+++
date = '2026-09-05T19:13:40+08:00'
draft = false
title = '算法竞赛数学公式与模板'
categories = ["数学"]
tags = ["多项式", "NTT", "FWT", "拉格朗日插值", "拉格朗日反演", "Berlekamp-Massey", "生成函数", "卡特兰数", "斐波那契", "欧拉数", "Polya计数", "牛顿迭代", "伯努利数", "分拆数", "数论"]
+++

# 数学

## 多项式

### 多项式基本运算板子

```cpp
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
namespace Poly {
	using poly = vector<int>;
	const int G = 3, invG = qpow(G, mod - 2);
	const int Imgunit = 911660635; // sqrt(-1)
	vector<int> rev, w, inv;
	void get_rev(int n) {
		rev.resize(n, 0);
		for(int i = 1; i < n; ++i){
			rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
		}
	}
	void get_inv(int n) {
		if (n < inv.size()) return; 
		int last = inv.size();
		inv.resize(n + 1);
		inv[1] = 1;
		for(int i = last; i <= n; ++i) {
			inv[i] = mul(mod - mod / i, inv[mod % i]);
		}
	}
	void NTT(poly &f, int n, bool I) {
		for(int i = 0; i < n; ++i){
			if(i < rev[i]) swap(f[i], f[rev[i] ]);
		}
		w.resize(n, 0);
		for(int p = 2; p <= n; p <<= 1){
			int len = p >> 1;
			w[0] = 1; w[1] = qpow(I ? invG : G, (mod - 1) / p);
			for(int i = 2; i < len; ++i){
				w[i] = mul(w[i - 1], w[1]);
			}
			for(int k = 0; k < n; k += p){
				for(int l = k; l < k + len; ++l){
					int tmp = mul(f[l + len], w[l - k]);
					f[l + len] = dec(f[l], tmp);
					f[l] = inc(f[l], tmp); 
				}
			}
		}
		if(I){
			int invn = qpow(n, mod - 2);
			for(int i = 0; i < n; ++i){
				f[i] = mul(f[i], invn);
			}
		}
	}
	poly poly_mod(poly f, int n) {
		f.resize(n, 0);
		return move(f);
	}
	poly operator + (poly f, poly g) {
		if(f.size() < g.size() ) swap(f, g);
		for(int i = 0; i < g.size(); ++i){
			f[i] = inc(f[i], g[i]);
		}
		return move(f);
	}
	poly operator - (poly f, poly g) {
		bool flag = 0;
		if(f.size() < g.size() ) swap(f, g), flag = 1;
		for(int i = 0; i < g.size(); ++i){
			f[i] = dec(f[i], g[i]);
		}
		if(flag)
			for(int i = 0; i < f.size(); ++i){
				f[i] = dec(0, f[i]);
			}
		return move(f);
	}
	poly operator * (poly f, poly g) {
		int n = f.size() - 1, m = g.size() - 1;
		for(m = n + m, n = 1; n <= m; n <<= 1);
		get_rev(n);
		f = poly_mod(f, n); 
		g = poly_mod(g, n);
		NTT(f, n, 0); NTT(g, n, 0);
		for(int i = 0; i < n; ++i) f[i] = mul(f[i], g[i]);
		NTT(f, n, 1);
		f = poly_mod(f, m + 1);
		return move(f);
	}
	poly poly_inv(const poly &f) {
		int n = f.size();
		poly f_inv(1, qpow(f[0], mod - 2) );
		for(int p = 2; (p >> 1) < n; p <<= 1){
			int len = p << 1;
			get_rev(len);
			f_inv = poly_mod(f_inv, len); 
			poly g = poly_mod(f, p);
			g = poly_mod(g, len);
			NTT(g, len, 0); NTT(f_inv, len, 0);
			for(int i = 0; i < len; ++i){
				f_inv[i] = mul(f_inv[i], dec(2, mul(f_inv[i], g[i]) ) );
			}
			NTT(f_inv, len, 1);
			f_inv = poly_mod(f_inv, p);
		}
		f_inv = poly_mod(f_inv, n);
		return move(f_inv);
	}
	poly poly_d(const poly &f) {
		int n = f.size() - 1;
		poly g(n + 1);
		for(int i = 0; i + 1 <= n; ++i) {
			g[i] = mul(f[i + 1], i + 1);
		}
		g[n] = 0;
		return move(g);
	}
	poly poly_i(const poly &f) {
		int n = f.size();
		poly g(n + 1);
		get_inv(n);
		for(int i = 1; i <= n; ++i) {
			g[i] = mul(f[i - 1], inv[i]);
		}
		g[0] = 0;
		return move(g);
	}
	poly poly_ln(const poly &f) {
		if(f[0] == 1) {
			int n = f.size();
			poly fp = poly_d(f);
			poly fi = poly_inv(f);
			fp = fp * fi;
			fp = poly_i(fp);
			fp = poly_mod(fp, n);
			return move(fp);
		}
		return {-1};
	}
	poly poly_exp(const poly &f) {
		if(f[0] == 0) {
			int n = f.size();
			poly g(1, 1);
			for(int t = 2; (t >> 1) < n; t <<= 1) {
				g = poly_mod(g, t);
				poly exp_t = poly_ln(g);
				exp_t[0] = dec(inc(f[0], 1), exp_t[0]);
				for(int i = 1; i < t; ++i) {
					if(i < n) exp_t[i] = dec(f[i], exp_t[i]);
					else exp_t[i] = dec(0, exp_t[i]);
				}
				g = g * exp_t;
			}
			g = poly_mod(g, n);
			return move(g);
		}
		return {-1};
	}
	poly poly_pow(const poly &f, int k, int k2 = 0, int k3 = 0) {
        // k2 是指数意义下的 k，k3 是 k 的一段足够大的前缀用于判断 x^ik
		poly g = f;
		if(f[0] == 1) {
			g = poly_ln(g);
			for(int i = 0; i < g.size(); ++i) {
				g[i] = mul(g[i], k);
			}
			g = poly_exp(g);
			return move(g);
		} 
		
		pair<int, int> p = {-1, -1};
		for(int i = 0; i < g.size(); ++i) {
			if(g[i] != 0) {
				p = {i, g[i]};
				break;
			}
		}
		if(p == make_pair(-1, -1) ) {
			return move(g);
		}
		int inv_f = qpow(p.second, mod - 2);
		for(int i = p.first; i < g.size(); ++i) {
			g[i - p.first] = mul(g[i], inv_f);
			if(p.first != 0) g[i] = 0;
		}

		g = poly_ln(g);
		for(int i = 0; i < g.size(); ++i) {
			g[i] = mul(g[i], k);
		}
		g = poly_exp(g);

		if(1ll * p.first * k3 > g.size() ) {
			p.first = g.size();
		} else p.first *= k3;
		p.second = qpow(p.second, k2);
		for(int i = g.size() - 1; i >= 0; --i) {
			if(i + p.first < g.size() ) g[i + p.first] = mul(g[i], p.second);
			if(p.first != 0) g[i] = 0;
		}
		return move(g);
	}
	pair<poly, poly> poly_div(poly f, poly g) {
		int n = f.size() - 1, m = g.size() - 1;
		if (n < m) return make_pair(poly{0}, f);
		int L = n - m + 1;
		poly q = g, t = f;
		reverse(q.begin(), q.end());
		reverse(t.begin(), t.end());
		q = poly_inv(poly_mod(q, L));
		q = poly_mod(q * t, L);
		reverse(q.begin(), q.end());
		g = poly_mod(g * q, m);
		for (int i = 0; i < m; ++i) {
			g[i] = dec(f[i], g[i]);
		}
		return make_pair(move(q), move(g));
	}
	poly poly_sqrt(poly f) {
		int n = f.size(), m = 1;
		while (m < n) m <<= 1;
		f = poly_mod(f, m);
		poly b = {1}, b2;
		for (int p = 2; (p >> 1) < n; p <<= 1) {
			int len = p >> 1;
			b2 = poly_mod(b2, p);
			for (int i = 0; i < len; ++i) b2[i] = mul(b[i], 2);
			b2 = poly_inv(b2);
			b = poly_mod(b * b, p);
			for (int i = 0; i < p; ++i) b[i] = inc(f[i], b[i]);
			b = poly_mod(b * b2, p);
		}
		b = poly_mod(b, n);
		return move(b);
	}
};
```

如果需要任意模数 NTT，这里提供三个模数 $998244353, 1004535809,469762049$，他们原根都是 $3$，最后 CRT 合并答案，然后取模即可。

### 多项式牛顿迭代

对于 $G(F(x))=0$，求解 $F(x) ()$

$$F(x) = F_*(x) - \frac{G(F_*(x))}{G'(F_*(x))}$$

### FWT 快速沃尔什变换

注意系数可以不是数，经典的子集卷积就是系数为多项式，然后外层做 Or，内层多项式维护长度。当然也可以把多项式提出到外层。这也算是不依赖单位根的好处。

$k$ 进制下，And (min) 的系数为上三角，Or (max) 的系数为下三角，Xor (不进位加法) 的系数是 $w^{ij}_k$（这种如果 $k > 2$ 就依赖单位根了）。

逆矩阵，前两个为差分的在对角边上的斜线设为 $-1$，后续的斜线为 $0$。Xor 的是改成 $w^{-ij}$，外面乘一个 $\frac{1}{k}$ 的系数。

```cpp
const int inv2 = qpow(2, mod - 2);
const int
	c_or[2][2] = {{1, 0}, {1, 1}},
	ic_or[2][2] = {{1, 0}, {mod - 1, 1}},
	c_and[2][2] = {{1, 1}, {0, 1}}, 
	ic_and[2][2] = {{1, mod - 1}, {0, 1}},
	c_xor[2][2] = {{1, 1}, {1, mod - 1}}, 
	ic_xor[2][2] = {{inv2, inv2}, {inv2, dec(mod, inv2)}};
void FWT(vector<int> &f, int n, const int c[2][2]) {
	for (int p = 2; p <= n; p <<= 1) {
		int len = (p >> 1);
		for (int k = 0; k < n; k += p) {
			for (int l = k; l < k + len; ++l) {
				int tl = f[l], tr = f[l + len];
				f[l] = inc(mul(c[0][0], tl), mul(c[0][1], tr));
				f[l + len] = inc(mul(c[1][0], tl), mul(c[1][1], tr));
			}
		}
	}
}
```

### 拉格朗日插值

$n$ 次多项式使用 $n + 1$ 个点值可以完成插值，复杂度 $O(n)$。

$$f(x) = \sum_{i = 0}^n y_i\prod_{j \ne i} \frac{x - x_j}{x_i - x_j}$$

如果点值连续，有：

$$f(x) = \sum_{i = 0}^n y_i \frac{\prod_{j \ne i} x - x_j}{(-1)^{n - i}(n - i)!i!}$$

如果点值呈现公差为 $d$ 的等差数列，有：

$$f(x) = \sum_{i = 0}^n y_i \frac{\prod_{j \ne i} x - x_j}{d^n(-1)^{n - i}(n - i)!i!}$$

### 拉格朗日反演

如果有 $G(F(x)) = x$，两者互为复合逆（也就是下面的式子是可交换的），有：

$$[x^n]F(x)=\frac{1}{n}[x^{-1}]G^{-m}(x)$$

$$[x^n]F(x)=\frac{1}{n}[x^{n - 1}](\frac{x}{G(x)})^n$$

$$[x^n]H(F(x)) = \frac{1}{n}[x^{n - 1}]H'(x)(\frac{x}{G(x)})^n$$

### 常系数其次线性递推

#### Berlekamp-Massey 算法

如果转移依赖项有 $k$ 个，一般拿出 $2k$ 项由 $bm$ 算法提取特征多项式，然后快速幂计算取模乘法即可。

```cpp
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
vector<int> poly_modmul(const vector<int> &a, const vector<int> &b, const vector<int> &bm_seq) {
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
vector<int> poly_modqpow(ll b, const vector<int> &bm_seq) {
	vector<int> res = {1}, a = {0, 1};
	while (b) {
		if (b & 1) res = poly_modmul(res, a, bm_seq);
		a = poly_modmul(a, a, bm_seq);
		b >>= 1;
	}
	return res;
}
```

使用例：

```cpp
vector<int> bm_seq = bm(seq);
vector<int> R = poly_modqpow(m - 1, bm_seq);
int result = 0;
for (int i = 0; i < R.size(); ++i) {
	result = inc(result, mul(R[i], seq[i]));
}
```

### 普通生成函数

#### 常用生成函数

乘 $\frac{1}{1-x}$ 可以得到系数前缀和的生成函数，相应的 $1-x$ 为差分。

$$\sum_{n = 0}^{\infty} x^n=\frac{1}{1-x}$$

$$\sum_{n = 0}^{\infty} (n + 1)x^n = \frac{1}{(1 - x)^2}$$

$$\sum_{n = 0}^{\infty} \binom{m}{n}x^n = (1 + x)^m$$

牛顿二项式定理：

扩展组合数： $\binom{m}{n} = \frac{m^{\underline n}}{n!}, m \in C, n \in N$

$$(1 + x)^{\alpha} = \sum_{n = 0}^{\infty} \binom{\alpha}{n}x^n$$ 

$$\sum_{n = 0}^{\infty} \binom{m + n}{n} x^n = \frac{1}{(1 - x)^{m + 1}}$$


### 指数生成函数

#### 定义

$$\sum_{n = 0}^{\infty} \frac{x^n}{n!} = e^x$$

#### 组合意义

##### 乘法

两个代表标号组件的指数型生成函数当成 $OGF$ 相乘会自动完成标号分配。

$$\hat{F}(x)\hat{G}(x) = \sum_{n = 0}^{\infty} \left [ \sum_{i = 0} ^ n \binom{n}{i} F_iG_{n -i} \right ] \frac{x^n}{n!}$$

就是总共 $n$ 个标号，我选出 $i$ 个标号给到组件 $F_i$，剩下的给到组件 $G_{n-i}$。

##### 幂次

多个相同的组件组合。 $\hat{F}^k(x)$ 相当于把 $n$ 个标号，分配给 $k$ 个有顺序的区别的相同类型的组件。

##### exp

$$\exp{(F(x))} = \sum_{k = 0}^{\infty} \frac{\hat{F}^k(x)}{k!}$$

正好把组件直接的顺序去掉了。

##### 求导和积分

求导：标定或移除某个特定的元素。

积分：加入一个元素。

### ！！！概率生成函数



### ！！！狄利克雷生成函数



## 组合

### 特殊数

#### 斐波那契

递推公式：$F_n = F_{n - 1} + F_{n - 2}$

斐波那契生成函数：

$$\frac{x}{1-x-x^2}= \sum_{n=0}^{\infty}\frac{(\frac{1 + \sqrt 5} 2)^n - (\frac{(1 - \sqrt 5)}2)^n}{\sqrt 5}$$

斐波那契前缀和 $S_n = f_{n + 2} - 1$

> **卡西尼性质** $F_{n - 1}F_{n + 1} - F_n^2= (-1)^n$

> **附加性质** $F_{n + k} = F_k F_{n + 1} + F_{k - 1}F_n$

> **整除性质** $F_a|F_b \Leftrightarrow a | b$ 

> **GCD性质** $F_{\gcd(a, b)} = \gcd(F_a, F_b)$

模 $m$ 意义下的周期不超过 $6m$

#### 卡特兰数

合法括号匹配数。

$$
C_n = \begin{cases}
1, & n = 0 \\
\sum_{i = 0}^nC_iC_{n - i}, & else \\
\end{cases}
$$

通项公式：

$$C_n = \frac{\binom{2n}{n}}{n + 1} = \binom{2n}{n} - \binom{2n}{n - 1}$$

递推公式：

$$C_n = \frac{4n - 2}{n + 1}C_{n - 1}$$

生成函数：

$$C^2(x) = \sum_{n=0}^{\infty}\sum_{i = 0}^n C_{i} C_{n-i}x^n$$

$$xC^2(x) + 1 = C(x)$$

所以有：

$$C(x) = \frac{1 - \sqrt{1 - 4x}}{2x}$$

#### 欧拉数

$$\sum_{n = 0}^{\infty} n ^k x^n = \frac{A_k(x)}{(1-x)^{k + 1}}$$

其中 $A_k(x)$ 为欧拉多项式，其为：

$$A_n(x) = \sum_{i = 1} ^n\left\langle \begin{matrix} k \\ m \end{matrix} \right\rangle x^n$$

$\left\langle \begin{matrix} n \\ m \end{matrix} \right\rangle$ 为欧拉数（Eulerian Number），指 $1$ 到 $n$ 的排列中，有恰好 $m$ 个数比前一个数大（认为开头也是上升），的排列的数量。 

递推公式：

$$\left\langle \begin{matrix} n \\ m \end{matrix} \right\rangle = m \left\langle \begin{matrix} n - 1 \\ m \end{matrix} \right\rangle + (n - m + 1) \left\langle \begin{matrix} n - 1 \\ m - 1 \end{matrix} \right\rangle$$

##### 自然数幂和

$$
\sum_{i=1}^n i^2=\frac{n(n+1)(2n+1)}{6}
$$

$$
\sum_{i=1}^n i^3=(\sum_{i=1}^ni)^2 
$$


#### 伯努利数

研究自然数幂和多项式的系数。

已知 $\{1, c, c^2, ...\}$ 的 EGF 为 $G_c(x) = e^{cx}$，那我们通过这个构造自然数幂和的 EGF。

$$S_n(x) = \sum_{c = 0}^{n - 1} G_c(x) = \frac{1 - e^{nx}}{1 - e^x}$$

这个东西的意义是，第 $i$ 项系数是前 $n$ 个数的每个幂 $i$ 下的和。

分离出伯努利数的 EGF，$B(x) = \frac{x}{e^x - 1}$，则：

$$S_n(x) = B(x) \frac{e^{nx} - 1}{x}$$

令：

$$G_n(x) = \frac{e^{nx} - 1}{x} = \sum_{i = 0}^{\infty} \frac{n^{i + 1}}{(i + 1)!}x^i$$


所以有：

$$\sum_{i = 0} ^ {n - 1}i^k = [\frac{x^k}{k!}]S_n(x) = k!\sum_{i + j = k} [x^i]B(x) [x^j]G_n(x) = k!\sum_{i = 0}^k\frac{(n + 1)^i}{(i + 1)!}B_{k - i}$$

对于伯努利数，直接多项式求逆求解即可 $B(x) = \frac{x}{e^x - 1}$。

#### 分拆数

正整数拆分的方案数。

$$P(x) = \prod_{k = 1}^{\infty} \sum_{i = 0}^{\infty} x^{ik} = \prod_{k = 1}^{\infty} \frac{1}{1 - x^k}$$

也就是先枚举每个数，然后枚举选几次。

乘法的直接取 $Ln$ 然后再 $Exp$ 回去。

$$P(x) = \prod_{k = 1}^{\infty}  \frac{1}{1 - x^k} = \exp\sum_{k = 1} ^{\infty} \ln\frac{1}{1 - x^k} = \exp\sum_{k = 1} ^{\infty} \sum_{i = 1}^{\infty} \frac{x^{ik}}{i}$$

### 斯特林数

#### 第一类斯特林数

##### 定义

$$\begin{bmatrix} n \\ k \end{bmatrix}$$

$n$ 个不同元素分为 $k$ 个互不区分的非空轮换的方案数。

##### 递推公式

$$\begin{bmatrix} n \\ k \end{bmatrix} = \begin{bmatrix} n - 1 \\ k - 1 \end{bmatrix} + (n - 1)\begin{bmatrix} n - 1 \\ k \end{bmatrix}$$

##### 同行计算

根据递推公式，我们需要求这样的一个多项式：

$$\prod_{i = 1}^n(x + i - 1)$$

分治 NTT 即可。

同行和：

$$\sum_{i =0}^n\begin{bmatrix} n \\ i \end{bmatrix} = n!$$

##### 同列计算

给出 EGF：

$$\frac{1}{k!}(- \ln(1 - x))^k$$

#### 第二类斯特林数

##### 定义

$$\begin{Bmatrix} n \\ k \end{Bmatrix}$$

$n$ 个不同元素分为 $k$ 个互不区分的非空集合的方案数。

##### 递推公式

$$\begin{Bmatrix} n \\ k \end{Bmatrix} = \begin{Bmatrix} n - 1 \\ k - 1 \end{Bmatrix} + k\begin{Bmatrix} n - 1 \\ k \end{Bmatrix}$$

##### 通项公式

$$
\begin{Bmatrix} n \\ k \end{Bmatrix}=\sum_{i=0}^k\frac{(-1)^{k-i}i^n}{i!(k-i)!}
$$

##### 同行计算

发现是一个卷积的形式，$f(x)=\sum_{i=0}^n\frac{(-1)^i}{i!}x^i$ 卷积 $g(x)=\sum_{i=0}^n\frac{i^n}{i!}$ 即可 $O(n\log n)$ 的时间内求出。

同行和，其中 $B_n$ 为贝尔数：

$$\sum_{i =0}^n\begin{Bmatrix} n \\ i \end{Bmatrix} = B_n$$

##### 同列计算

给出 EGF：

$$\frac{1}{k!}(e^x - 1)^k$$

同时，这个 $e^x - 1$ 的 exp 是贝尔数的生成函数。

### 组合恒等式与反演

#### 上升幂与下降幂

##### 上升幂-普通幂

$$
x^{\overline{n}}=\sum_{i=0}^n\begin{bmatrix} n \\ i \end{bmatrix}x^i
$$

$$
x^n=\sum_{i=0}^n\begin{Bmatrix} n \\ i \end{Bmatrix}(-1)^{n-i}x^{\overline i}
$$

##### 下降幂-普通幂

$$
x^{\underline n}=\sum_{i=0}^n\begin{bmatrix} n \\ k \end{bmatrix} (-1)^{n-i}x^i
$$

$$
x^n=\sum_{i=0}^n\begin{Bmatrix} n \\ i \end{Bmatrix}x^{\underline i}
$$

##### 牛顿（下降幂）多项式系数与点值表示

假设有：

$$f(x) = \sum_{i = 0}^{n} b_i x^{\underline i}$$

和一组点值：

$$(i, a_i), i \in [0, n]$$

则：

$$\frac{a_k}{k!} = \sum_{i = 0}^k\frac{b_i}{(k - i)!}$$

卷积可以完成插值。

#### 二项式反演和斯特林反演

##### 二项式反演

$$
f(n)=\sum_{i=0}^n\binom{n}{i}g(i)\Leftrightarrow g(n)=\sum_{i=0}^n(-1)^{n-i}\binom{n}{i}f(i)
$$


$$
f(n)=\sum_{i=0}^n(-1)^i\binom{n}{i}g(i)\Leftrightarrow g(n)=\sum_{i=0}^n(-1)^i\binom{n}{i}f(i)
$$

对于其 EGF，有：

$$f * e^x = g, g * e^{-x} = f$$

##### 斯特林反演

下式的一二类斯特林和 $-1$ 可以任意交换。

$$
f(n)=\sum_{i=0}^n\begin{Bmatrix} n \\ i \end{Bmatrix}g(i)\Leftrightarrow g(N)=\sum_{i=0}^n(-1)^{n-i}\begin{bmatrix} n \\ i \end{bmatrix}f(i)
$$


$$
f(n)=\sum_{i=0}^n(-1)^i\begin{Bmatrix} n \\ i \end{Bmatrix}g(i)\Leftrightarrow g(N)=\sum_{i=0}^n(-1)^{i}\begin{bmatrix} n \\ i \end{bmatrix}f(i)
$$

##### 反转公式

$$
\sum_{i=m}^n(-1)^{n-i}\begin{Bmatrix} n \\ i \end{Bmatrix}\begin{bmatrix} i \\ m \end{bmatrix}=[m=n]\tag{1} 
$$
$$
\sum_{i=m}^n(-1)^{m-i}\begin{bmatrix} n \\ i \end{bmatrix}\begin{Bmatrix} i \\ m \end{Bmatrix}=[m=n]\tag{2}
$$

互为反演的都有上述类似形式的性质。

#### min-max 容斥

对于 $n$ 长全序序列 $\{x_i\}$ ，$S=\{1,2,...,n\}$，有：
$$
\text{kth}\max_{i\in S}{x_i} =\sum_{T\subseteq S}(-1)^{|T|-k}\binom{|T|-1}{k-1}\min_{j\subseteq T} x_j
$$
$$
\text{kth}\min_{i\in S}{x_i}=\sum_{T\subseteq S}(-1)^{|T|-k}\binom{|T|-1}{k-1}\max_{j\subseteq T} x_j
$$

#### 单位根反演

$$[n | a] = \frac{1}{n} \sum_{i = 0}^{n - 1} w_n^{ia}$$

#### 二项式恒等式

##### 插板法

$n$ 个相同元素分为 $k$ 个有区别的非空集合的方案数为 $\binom{n - 1}{k - 1}$。

##### 乘式

$$\binom{n}{r}\binom{r}{k} = \binom{n}{k}\binom{n - k}{r - k}$$

##### 上指标反转

$$\binom{r}{k} = (-1)^k \binom{k - r - 1}{k} $$

##### 变上项求和

$$\sum_{l = 0}^n\binom{l}{k} = \binom{n + 1}{k + 1}$$

##### 范德蒙德卷积

$$\sum_{k = 0}^{r}\binom{m}{k}\binom{n}{r - k} = \binom{n + m}{r}$$

变式，只要能弄出底下和为定值，上面两个定值就可以直接套。

$$\sum_{k = 0}^r\binom{k}{a}\binom{r - k}{b} = \binom{r + 1}{a + b + 1}$$

变式，只要能弄出上面和为定值，底下两个定值就可以直接套。

##### 斐波那契恒等式

$$\sum_{i = 0}^{n} \binom{n - i}{i} = F_{n + 1}$$

### polya 计数

#### Burnside 引理

给定群 $G$ 在集合 $X$ 上的作用，则所有不同的轨道的数目

$$
|X/G| = \frac{1}{|G|}\sum_{g \in G} |X^g|
$$

其中 $X^g = \{x \in X|gx = x\}$ ，表示变换 $g$ 的不动点集合。

应用时列出所有对称操作，然后求出各个操作的不动点个数即可。

#### polya 定理

如果将 $G$ 看作对于原本可染色部分的结构的一种置换，那么对于 $g \in G$，他的可以独立染色的部分为其轮换表示时其的轮换个数，我们记作 $c(g)$。 

如果可用颜色有 $m$ 种，那么其不动点个数为 $m ^ {c(g)}$。

回代 Burnside 引理中可以得到：

$$
|X/G| = \frac{1}{|G|}\sum_{g \in G} m ^ {c(g)}
$$

## 数论



## 线性代数

