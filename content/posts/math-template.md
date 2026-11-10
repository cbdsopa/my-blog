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
		int last = max<int>(2, inv.size());
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
		poly b = {1}, b2; // b[0] 是 f[0] 的模意义下开根
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

$k$ 进制下，And (min) 的系数为上三角，Or (max) 的系数为下三角，Xor (不进位加法) 的系数是 $w^{ij}_k$（这种依赖单位根了）。

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

### FFT 快速傅里叶变换

```cpp
struct cp{
	double x,y;
	cp(double _x=0,double _y=0){
		x=_x;y=_y;
	}
	inline friend cp operator+(cp a,cp b){
		return cp(a.x+b.x,a.y+b.y);
	} 
	inline friend cp operator-(cp a,cp b){
		return cp(a.x-b.x,a.y-b.y);
	}
	inline friend cp operator*(cp a,cp b){
		return cp(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
	}
};
const double pi = acos(-1);
namespace Poly {
	using poly = vector<cp>;
	vector<int> rev;
    poly w;
	void get_rev(ll n) {
		rev.resize(n, 0);
		for(ll i = 1; i < n; ++i){
			rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
		}
	}
	void FFT(poly &f, ll n, bool I) {
		for(ll i = 0; i < n; ++i){
			if(i < rev[i]) swap(f[i], f[rev[i] ]);
		}
		w.resize(n);
		for(ll p = 2; p <= n; p <<= 1){
			ll len = p >> 1;
			w[0] = 1; w[1] = cp(cos(2 * pi / p), sin(2 * pi / p));
            if (I) w[1].y *= -1;
			for(ll i = 2; i < len; ++i){
				w[i] = w[i - 1] * w[1];
			}
			for(ll k = 0; k < n; k += p){
				for(ll l = k; l < k + len; ++l){
					cp tmp = f[l + len] * w[l - k];
					f[l + len] = f[l] - tmp;
					f[l] = f[l] + tmp; 
				}
			}
		}
		if(I){
			for(ll i = 0; i < n; ++i){
				f[i].x /= n;
			}
		}
	}
	poly poly_mod(poly f, ll n) {
		f.resize(n, 0);
		return move(f);
	}
	poly operator * (poly f, poly g) {
		ll n = f.size() - 1, m = g.size() - 1;
		for(m = n + m, n = 1; n <= m; n <<= 1);
		get_rev(n);
		f = poly_mod(f, n); 
		g = poly_mod(g, n);
		FFT(f, n, 0); FFT(g, n, 0);
		for(ll i = 0; i < n; ++i) f[i] = f[i] * g[i];
		FFT(f, n, 1);
		f = poly_mod(f, m + 1);
		return move(f);
    }
};
```

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

### 概率生成函数

对于随机变量 $X$，有 $P(X = i) = f_i$，则：

$$F(x) = \sum_{i = 0}^n f_i x^i$$

有：
 
$$F(1) = 1$$

$$E(X) = F'(1)$$

$$E(X^{\underline k}) = F^{(n)}(1)$$

$$E(X^2) = F'(1) + F''(1)$$

$$D(X) = F'^2(1) - F'(1) - F''(1)$$

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

#### 杨表

同行单调递增，行长单调不增。

插入时，在当前行二分找到最小的比他大的替换，如果是最大的则无需替换，放在末尾直接结束；否则拿着替换得到的数去下一行接着做。

##### 勾长公式

记 $\operatorname{hook}(x)$ 表示方格 $x$ 同列下方、同行右边、自身格子的数量。

记 $\operatorname{dim} \pi_\lambda$ 表示杨表 $\pi_\lambda$ 中填入 $[1,n]$ 中的整数，且行列全都单调递增的方案数。

有：

$$\operatorname{dim}\pi_\lambda=\frac{n!}{\prod \operatorname{hook}(x)}$$

一个排列对应一对杨表，一个结果表，一个时间表，所以对于一个杨表的形状，代入勾长公式的平方即为对应排列数。

##### LIS

LIS 长度不超过 $k$ 的子序列的最长长度，为杨表前 $k$ 列的长度和。

序列翻转，得到的杨表是原表的转置。

比较方式反转，得到的杨表是原表的转置。

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
\sum_{i=m}^n(-1)^{n-i}\begin{Bmatrix} n \\ i \end{Bmatrix}\begin{bmatrix} i \\ m \end{bmatrix}=[m=n]
$$
$$
\sum_{i=m}^n(-1)^{m-i}\begin{bmatrix} n \\ i \end{bmatrix}\begin{Bmatrix} i \\ m \end{Bmatrix}=[m=n]
$$

互为反演的都有上述类似形式的性质。

#### min-max 容斥

对于 $n$ 长全序序列 $\{x_i\}$ ，$S=\{1,2,...,n\}$，有：
$$
\operatorname{kth}\max_{i\in S}{x_i} =\sum_{T\subseteq S}(-1)^{|T|-k}\binom{|T|-1}{k-1}\min_{j\subseteq T} x_j
$$
$$
\operatorname{kth}\min_{i\in S}{x_i}=\sum_{T\subseteq S}(-1)^{|T|-k}\binom{|T|-1}{k-1}\max_{j\subseteq T} x_j
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

## 抽象代数

### 群论

#### Lagrange 定理

同一子群的不同陪集大小都相同，都等于对应子群的大小。子群的全体陪集构成群的一个分划，有限群的阶必然是子群的阶的整数倍。也就是 $|G| = [G:H]|H|$

$[G:H]$ 为 $H$ 在 $G$ 中的陪集个数，称为群 $G$ 中 $H$ 的指数。


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
#### 有限生成 Abel 群

如果群 $G$ 有一个有限的生成子集，则称群 $G$ 是有限生成的。

##### 直积

群 $(G,\cdot_G)$ 和群 $(H,\cdot_H)$ 的直积，为群 $(G\times H, \cdot)$，且满足 $(g_1, h_1)\cdot(g_2, h_2)=(g_1\cdot_G g_2, h_1\cdot_H h_2) $。

对于 $G$ 的子群 $H_1, H_2$，$G\cong H_1 \times H_2$ 当且仅当 $H_1, H_2$ 都是正规子群，$H_1\cup H_2 = \{e\}$ 且 $G=H_1H_2$

也就是满足：全覆盖，零重叠，零耦合。

##### 分类定理

对于有限生成 Abel 群 $G$，实现如下拆分（$r \ge 0, n_i \ge 2$）：

$$G \cong C_\infty^r \times C_{n_1} \times \cdots \times C_{n_s}$$

其中 $r$ 为 $G$ 的秩。

拆分方法：

一、质数幂次分解。此时称 $C_{n_i}$ 为初等因子。

二、倍数分解。前一个 $n_i$ 是 $n_{i + 1}$ 的因子。此时称 $C_{n_i}$ 为不变因子。

### 环论

#### 整环

可交换，有幺元，无零因子的环。

##### 欧几里得整环

定义一个范数 $N : R \setminus \{0\} \rightarrow N$，满足可以 $a, b\in R, b \neq0, \exist q, r\in R$ 满足 $a = qb + r$，且 $r = 0$ 或 $N(r) < N(b)$。

可以带余除法以及辗转相除，存在最大公约数，满足裴蜀定理。

##### 主理想整环

欧几里得整环弱化，满足所有理想都是主理想。

存在最大公约数，满足裴蜀定理。

##### 唯一分解整环

主理想整环弱化，满足算术基本定理。

任何非零且不可逆的元素 $r$ 可以唯一分解为多个不可约分元素的乘积。

此时对于非零元素，不可约分等价于其为素元。

存在最大公约数。

#### 多项式环

##### Guass 引理

对于唯一分解整环 $R$ 及其分式域 $F$，次数大于 $0$ 的多项式 $f(x) \in R[x]$ 在 $R[x]$ 中不可约，那么它在 $F[x]$ 中不可约。 

当 $f(x)$ 为本原多项式时其系数互素，此时上述推导关系变为等价关系。

##### Eisenstein 判别法

对于 $n$ 次整系数多项式 $f(x) = a_0+a_1x+...+a_nx^n$，如果存在质数 $p$ 满足 $p \mid a_i$ 对 $0, 1, ... , n - 1$ 成立，且 $p \nmid a_n$， $p^2 \nmid a_0$，则 $f(x)$ 在有理域 $Q$ 上 不可约分。如果其还是本原多项式，则在整数域 $Z$ 上也不可约分。 

#### 中国剩余定理

称理想 $I_1, I_2 \le R$ 互素，当 $I_1 + I_2 = R$。

现有非零交换幺环 $R$，以及其理想 $I_1, I_2, ..., I_n$，如果它们两两互素，则 $\varphi : R \rightarrow I_1 \times I_2  ...\times I_n$ 是满射，$\operatorname{ker} \varphi = I_1\cap I_2 \cap... \cap I_n=I_1I_2...I_n$ 

有：

$R / (I_1I_2...I_n)= R / (I_1\cap I_2 \cap... \cap I_n) \cong R/I_1 \times R/I_2 \times ... \times R/ I_n$

##### 整数同余类的乘法群

整数模 $n$ 乘法群，是 $Z/nZ$ 中取出所有可逆元得到的群，记为 $(\mathbf Z/n\mathbf Z)^{\times}$。阶为 $\varphi(n)$，$n$ 可分为：

$$n=p_1^{\alpha_1}\cdots p_s^{\alpha_s}$$

$$\mathbf Z/n\mathbf Z\cong\mathbf Z/p_1^{\alpha_1}\mathbf Z\times\cdots\times\mathbf Z/p_s^{\alpha_s}\mathbf Z$$

$$(\mathbf Z/n\mathbf Z)^\times\cong(\mathbf Z/p_1^{\alpha_1}\mathbf Z)^\times\times\cdots\times(\mathbf Z/p_s^{\alpha_s}\mathbf Z)^\times$$

其中：

$(\mathbf Z/2\mathbf Z)^{\times}\cong C_1$，$(\mathbf Z/2^2\mathbf Z)^{\times}\cong C_2$，对于 $k \ge 2$ 有 $(\mathbf Z/2^k\mathbf Z)^{\times}\cong C_2 \times C_{2^{k - 2}}$。对于奇素数 $p$，$(\mathbf Z/p^k\mathbf Z)^{\times} \cong C_{\varphi(p^k)}$

当 $n = 1, 2, 4, p^k, 2p^k$ 时，没有子群 $C_2\times C_2$ ，此时为循环群，存在原根。

## ！数论

### 素性测试和 Pollard-rho

```cpp
std::mt19937_64 rd(time(0)*1000+clock() );
bool Fermat(ll p) {
	if (p < 3) return p == 2;
	for (int i = 1; i <= 10; ++i){
		int a = rd() % (n - 2) + 2;
		if (qpow(a, p - 1, p) != 1) return 0;
	}
	return 1;
}
int ts[8] = {2, 3, 5, 7, 13, 19, 61, 233};
bool MillerRabin(ll p) {
	if (p < 3) return p == 2;
	if (!(p & 1)) return 0; 
	ll t = p - 1, k = 0;
	while (!(t & 1)) { t >>= 1; ++k; }
	for (int i = 0; i < 8; ++i){
		ll a = qpow(ts[i], t, p);
		for (int j = 0; j < k; ++j){
			ll b = mul(a, a, p);
			if (b == 1 && a != 1 && a != p - 1) return 0;
			a = b;
		}
		if (a != 1) return 0;
	}
	return 1;
}
ll PollardRho(ll p) {
	static ll sav[128 + 3];
	int top = 0;
	ll c = rd() % (p - 1) + 1;
	ll x, y = 0, buf = 1;
	for (int st = 1; ; st<<=1){
		x = y;
		for (int i = 0; i < st; ++i){
			y = inc(mul(y, y, p), c, p);
			buf = mul(buf,dec(y, x, p), p);
			sav[++top] = dec(y, x, p);
			if (top == 128){
				if (gcd(buf, p) > 1) break;
				top = 0; buf = 1;
			}
		}
		if (top == 128) break;
	}
	for (int i = 1; i <= top; ++i){
		buf = gcd(sav[i], p);
		if(buf > 1) return buf;
	}
	return p;
}
void divide(ll p) {
	if (p < 2) return;
	if (MillerRabin(p)){
		// 回收因子
		return;
	}
	ll q = PollardRho(p);
	while (p == q) q = PollardRho(p);
	divide(p / q); divide(q);
}
```

### 扩展欧几里得

```cpp
int exgcd(int a, int b, int &x, int &y) {
	if (!b) {
		x = 1; y = 0;
		return a;
	}
	int d = exgcd(b, a % b, x, y);
	int temp = x;
	x = y;
	y = temp - (a / b) * y;
	return d;
}
```

### 欧拉定理

$$
a^k \equiv \begin{cases}
a^{k \bmod \varphi(m)},                &\gcd(a,m) =  1,                   \\
a^k,                                   &\gcd(a,m)\ne 1, k <   \varphi(m), \\
a^{(k \bmod \varphi(m)) + \varphi(m)}, &\gcd(a,m)\ne 1, k \ge \varphi(m).
\end{cases} \pmod m
$$

### 中国剩余定理

$n_i$ 两两互质。

$$
\begin{cases}
x &\equiv a_1 \pmod {n_1} \\
x &\equiv a_2 \pmod {n_2} \\
  &\vdots \\
x &\equiv a_k \pmod {n_k} \\
\end{cases}
$$

每一维构造仅这一维度为 $1$，其他维度为 $0$。

取 $m_i = \frac{n}{n_i}$，并计算其在模 $n_i$ 下逆元。

然后计算 $c_i = m_i m_i^{-1}$，不需要对 $n_i$ 取模，这个数在模 $n_i$ 意义下为 $1$，而在其他的模数 $n_j$ 意义下，$m_i \equiv 0 \pmod{n_j}$，所以成功构造单位元素。

后续按比例乘积求和即可。

#### 扩展中国剩余定理

考虑 $x\equiv a_1 \pmod{m_1} $ 和 $x\equiv a_2 \pmod{m_2}$，有 $x=m_1 p_1 +a_1=m_2p_2+a_2$，得 $m_1p_1-m_2p_2=a_2-a_1$ 用 `exgcd` 可解得 $p_1, p_2$（要求 $a_2-a_1$ 被 $\gcd(m_1, m_2) $ 整除）。

得到新方程 $x\equiv m_1p_1+a_1 \pmod{ {\operatorname{lcm} (m_1, m_2)}}$

### 阶乘取模

#### 升幂引理

##### 引理1

素数 $p$ 和所有满足 $(n, p) = 1$ 的 $n$：

若 $ p \mid x - y$，则：

$$
\nu_p(x^n-y^n)=\nu_p(x-y)
$$

若 $p \mid x + y$，则对于奇数 $n$：

$$
\nu_p(x^n+y^n)=\nu_p(x+y)
$$

##### 引理2

若 $p$ 是奇素数：

若 $ p \mid x - y$，则：

$$
\nu_p(x^n-y^n)=\nu_p(x-y)+\nu_p(n)
$$

若 $p \mid x + y$，则对于奇数 $n$：

$$
\nu_p(x^n+y^n)=\nu_p(x+y)+\nu_p(n)
$$

##### 引理3

对于 $p = 2$，且 $ p \mid x - y$，则：

对奇数 $n$：
$$
\nu_p(x^n-y^n)=\nu_p(x-y)
$$

对偶数 $n$：

$$
\nu_p(x^n+y^n)=\nu_p(x-y)+\nu_p(x+y)+\nu_p(n)-1
$$

若 $4 \mid x - y$：

$$
\nu_2(x + y) = 1\\
\nu_p(x^n+y^n)=\nu_p(x-y)+\nu_p(n)
$$

#### 威尔逊定理

当 $n$ 为素数，有：

$$(n - 1)!\equiv -1 \pmod n$$

##### 扩展

对于自然数 $m$，有：

$$
\prod_{1\le k < m, k\perp m = 1} k = \pm 1 \pmod m
$$

取 $-1$ 当且仅当原根存在。

#### Legendre 公式

$$
\nu_p(n!) = \sum_{i=1}^{\infty} \left\lfloor \dfrac{n}{p^i} \right\rfloor = \dfrac{n-S_p(n)}{p-1} 
$$

$S_p(n)$ 为 $n$ 在 $p$ 进制下各个数位之和。

#### Kummer 定理

$$
\nu_p\left(\dbinom{m}{n}\right)=\frac{S_p(n)+S_p(m-n)-S_p(m)}{p-1}.
$$

预处理阶乘之后递归计算。

#### 阶乘余数

对于素数 $p$：

$$
(n!)_p \equiv (-1)^{\left\lfloor n/p\right\rfloor}\cdot (n\bmod p)!\cdot\left(\left\lfloor n/p\right\rfloor!\right)_p\pmod{p}.
$$


对于素数 $p$ 和正整数 $\alpha, n$

$$
(n!)_{p} \equiv (\pm 1)^{\lfloor n/p^\alpha\rfloor}\cdot\left(\prod_{1\le j\le (n\bmod p^\alpha),\ j\perp p}j\right)\cdot(\lfloor n/p\rfloor!)_p\pmod{p^\alpha}.
$$

$\pm 1$ 取值同威尔逊定理。

### 卢卡斯定理

对于素数 $p$，有：

$$
\binom{n}{k}\equiv \binom{\lfloor \frac{n}{p}\rfloor}{\lfloor \frac{k}{p}\rfloor}\binom{n\bmod p}{k\bmod p}\pmod p.
$$

#### 扩展卢卡斯定理

考虑模数为 $p^\alpha$，因此，组合数可以写作：

$$
\binom{n}{k} = p^{\nu_p(n!)-\nu_p(k!)-\nu_p((n-k)!)}\dfrac{(n!)_p}{(k!)_p((n-k)!)_p}.
$$

后续 `CRT` 合并即可。

### 原根与单位根

#### 阶

##### 性质1

$$
\delta_m(a^k) = \dfrac{\delta_m(a)}{(\delta_m(a),k)}.
$$

##### 性质2

对于 $a,b\in\mathbf Z,m\in\mathbf N_+$ 且$a,b\perp m$，有：

$$
\delta_m(ab) = \delta_m(a)\delta_m(b) \iff \delta_m(a)\perp\delta_m(b)
$$

##### 性质3

对于 $a,b\in\mathbf Z,m\in\mathbf N_+$ 且 $a,b\perp m$，总是存在 $c\in\mathbf Z$ 且 $c\perp m$ 使得：

$$
\delta_m(c) = [\delta_m(a),\delta_m(b)]
$$

#### 原根

##### 原根判定定理

对于整数 $m\ge 3$ 和 $g\perp m$，$g$ 是模 $m$ 的原根，当且仅当对于 $\varphi(m)$ 的每个素因数 $p$，都有 

$$
g^{\frac{\varphi(m)}{p}}\not\equiv 1 \pmod m
$$

##### 原根个数

如果正整数 $m$ 有原根 $g$，那么，当且仅当 $d\mid\varphi(m)$ 时，模 $m$ 的 $d$ 阶元素存在，且恰有 $\varphi(d)$ 个。

##### 原根存在定理

模 $m$ 的原根存在，当且仅当 $m=1,2,4,p^e,2p^e$，其中，$p$ 是奇素数且 $e\in\mathbf N_+$．

### 离散对数

如果存在 $x^2\equiv a\pmod p$ ，则称 $x$ 为模 $p$ 的二次剩余。

先声明 `Legendre` 记号:

对奇素数 $p$，
$$
\left(\frac{a}{p}\right)=\begin{cases}
    0,  & p\mid a,\\
    1,  & (p\nmid a) \land ((\exists x\in\mathbf{Z}),~~a\equiv x^2\pmod p),\\
    -1, & \text{otherwise}.\\
\end{cases}
$$

#### Euler 叛别准则

若奇素数 $p$ 满足 $p\nmid a$，则：

$$
a^{\frac{p-1}{2} }\equiv \left(\frac{x}{p}\right) 
$$

`Legendre` 记号为完全积性函数。

小于 $p$ 的正数中二次剩余的数量与非二次剩余数量相等。

可以扩展 $n$ 次剩余，解数为 $n$。

#### Cipolla算法

考虑求 $x^2\equiv n\pmod p$。

先 `Euler判别准则` 判无解，再特判 $0$。

如果有解则随意选一个数 $a$ 开始。这个 $a$ 需要满足，$a^2-n$ 为非二次剩余。

这个 $a$ 我们随机算取，由上面的内容我们知道概率为 $\frac{1}{2}$。

然后设 $i^2=(a^2-n)$ ，然后答案就为 $\pm(a+i)^{\frac{p+1}{2} }$。这个 $i$ 我们视作虚数，因为实际上 $a^2-n$ 是非二次剩余，这个 $i$ 实际上并不存在实际的意义。

因为我们知道:

$$
(a+i)^{p+1}=(a^p+i^p)(a+i)=(a+i)(a-i)=n
$$

第三步的理由是 $i^p=i^{p-1}i=(a^2-n)^{\frac{p-1}{2}}i=-i$ ，由 `Euler判别准则`。

#### BSGS

解关于 $x$ 的同余方程 $a^x\equiv n\pmod p$

写成这种形式查哈希表即可。

$$
a^{x\times B}\equiv na^{-y} \pmod p
$$

如果 $(a, p) \neq 1$ 可以一直除 $\gcd(a, p)$。 

### 数论函数

#### 狄利克雷生成函数

定义 Dirichlet 级数：

$$
F(s) = \sum_{n=1}^{\infty}\dfrac{f(n)}{n^s}
$$

另一种形式是：

$$
F(x) = \prod_{p\in \mathbf P}(\sum_{i\ge 0} \frac{f(p^i)}{p^{ix}})
$$

##### 常见狄利克雷生成函数

$$
E(s) = \sum_{n=1}^{\infty}\dfrac{\varepsilon(n)}{n^s} = 1.
$$ 

$$
I(s) = \sum_{n=1}^{\infty}\dfrac{1}{n^s} = \prod_{p\in\mathbf P}\dfrac{1}{1-p^{-s}} = \zeta(s).
$$

$$
M(s) = \sum_{n=1}^{\infty}\dfrac{\mu(n)}{n^s} = \prod_{p\in\mathbf P}(1-p^{-s}) = \dfrac{1}{\zeta(s)}.
$$

$$
I_k(s) = \sum_{n=1}^{\infty}\dfrac{n^k}{n^s} = \prod_{p\in\mathbf P}\dfrac{1}{1-p^{k-s}} = \zeta(s-k).
$$

$$
\begin{aligned}
\Phi(s) &= \prod_{p\in\mathbf P}\left(1+\dfrac{p-1}{p^s} + \dfrac{p(p-1)}{p^{2s}} + \dfrac{p^2(p-1)}{p^{3s}}+\cdots\right)\\
&= \prod_{p\in\mathbf P}\left(\dfrac{1}{1-p^{1-s}}-\dfrac{1}{p^s}\dfrac{1}{1-p^{1-s}}\right) = \prod_{p\in\mathbf P}\dfrac{1-p^{-s}}{1-p^{1-s}} = \dfrac{\zeta(s-1)}{\zeta(s)}.
\end{aligned}
$$

$$
\begin{aligned}
\Sigma_k(s) &= \prod_{p\in\mathbf P}\left(1+\dfrac{1+p^k}{p^s}+\dfrac{1+p^k+p^{2k}}{p^{2s}}+\dfrac{1+p^k+p^{2k}+p^{3k}}{p^{3s}}+\cdots\right) \\
&= \prod_{p\in\mathbf P}\dfrac{1}{1-p^k}\left((1-p^k)+\dfrac{1-p^{2k}}{p^s}+\dfrac{1-p^{3k}}{p^{2s}}+\dfrac{1-p^{4k}}{p^{3s}}+\cdots\right)\\
&= \prod_{p\in\mathbf P}\dfrac{1}{1-p^k}\left(\dfrac{1}{1-p^{-s}} - \dfrac{p^k}{1-p^{k-s}}\right)\\
&= \prod_{p\in\mathbf P}\dfrac{1}{(1-p^{-s})(1-p^{k-s})} = \zeta(s-k)\zeta(s).
\end{aligned}
$$

### 筛法

#### 杜教筛

假设 $f\ast g = h$, $g, h$ 已知，求 $f$ 的前缀和 $S(n) = \sum_{i = 1} ^ n f(i)$。


$$
\sum_{i = 1}^n h(i) = \sum_{i = 1}^n \sum_{d\mid i}f(d)g(\frac{i}{d})\\
=\sum_{d = 1}^n g(d)\sum_{j = 1}^{\lfloor\frac{n}{d}\rfloor}f(j)\\
=g(1)S(n)+\sum_{d=2}^n g(d)S(\lfloor\frac{n}{d}\rfloor)
$$

得：

$$
g(1)S(n) = \sum_{i = 1}^n h(i) - \sum_{d=2}^n g(d)S(\lfloor\frac{n}{d}\rfloor)
$$

采用记忆化后复杂度 $O(n^{\frac{3}{4}})$。如果预处理前 $n^{\frac{2}{3}}$ 个，则复杂度 $O(n^{\frac{2}{3}})$。

## 数值算法

### 自适应辛普森积分

```cpp
double intergral(double l, double r, double eps, int t){
	auto simpson(double l, double r) -> double {
		return (r - l) / 6 * (func(l) + 4 * func((l + r) / 2) + func(r) );
	}
	double mid = (l + r) / 2;
	double f = simpson(l, r);
	double fl = simpson(l, mid);
	double fr = simpson(mid, r);
	if(fabs(f - fl - fr) <= eps && t <= 0)
		return f;
	return intergral(l, mid, eps / 2, t - 1) + intergral(mid, r, eps / 2, t - 1);
}
```

## 博弈论

### SG 定理

对于公平组合游戏 $G_1, G_2, ..., G_n$：

$$SG(x) = \operatorname{mex} \{SG(x'):x'\in x\}$$

$$SG(G_1+G_2+...+G_n)=\operatorname{SG}(G_1)\oplus \operatorname{SG}(G_2)\oplus\cdots\oplus\operatorname{SG}(G_n).$$

其中公平游戏 $G$ 中的一个状态 $x$ 是先手必胜状态，当且仅当 
$\operatorname{SG}(x)\neq 0$。

