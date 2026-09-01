+++
date = '2026-09-01T16:24:32+08:00'
draft = false
title = '算法竞赛计算几何常用公式与模板'

categories = ["计算几何"]
tags = ["凸包", "半平面交", "旋转卡壳", "平面最近点对", "闵可夫斯基和", "反演变换", "极角排序", "动态凸包", "李超线段树", "Pick定理"]
+++

# 计算几何与浮点精度

## 浮点处理

### 符号判定

```cpp
const db eps = 1e-9;
int sign(db x) { return x > eps ? 1 : x > -eps ? 0 : -1; }
```

### 数乘压缩

当多个数字乘积很大时，采用 $\log$ 对数字进行压缩，将乘积转为加和，最后再使用 $\exp$ 反解出值。

经典例子：压缩组合数计算中的阶乘，从而压缩组合数。

### 对数换底公式

$$log_ba = \frac{\log_c a}{\log_c b} = \frac{\ln a}{\ln b}$$

### 高精度计算圆周率

```cpp
const db pi = acos(-1);
```

## 点线相关

### 向量类

#### 性质

$$ \vec a \cdot \vec b = |\vec a||\vec b|\cos\theta$$

$$ \vec a \times \vec b = |\vec a||\vec b|\sin\theta$$

叉积大于 $0$ 说明 $b$ 在 $a$ 的逆时针。

#### 代码

```cpp
const db eps = 1e-9;
int sign(db x) { return x > eps ? 1 : x > -eps ? 0 : -1; }
struct vec {
	db x = 0, y = 0;
	vec() {}
	vec(db _x, db _y) : x(_x), y(_y) {}
	vec operator + (const vec &rhs) const { return {x + rhs.x, y + rhs.y}; }
	vec operator - (const vec &rhs) const { return {x - rhs.x, y - rhs.y}; }
	vec operator * (const db &k) const { return {k * x, k * y}; }
	vec operator / (const db &k) const { return {x / k, y / k}; }
	db operator * (const vec &rhs) const { return x * rhs.y - y * rhs.x; }
	db operator / (const vec &rhs) const { return x * rhs.x + y * rhs.y; }
	bool operator == (const vec &rhs) const { return !sign(x - rhs.x) && !sign(y - rhs.y); }
	db length() const { return sqrt(x * x + y * y); }
	db length2() const { return x * x + y * y; }
	vec normalized() const { return *this / length(); } // notice (0, 0)
	bool is_zero() const { return length() < eps; }
	vec rotate90() const { return {-y, x}; }
	vec rotate(db ang) const { // radian
		db s = sin(ang), c = cos(ang);
		return {x * c - y * s, x * s + y * c};
	}
};
```
### 极角排序

最方便的是直接按 $atan2$ 算出角度来排序，但是他是按 $x$ 轴负向给隔开，角度范围 $(-\pi, \pi]$，而且精度差一点。

如果采用叉积，当所有向量位于一个半平面内时，直接按叉积排序即可。否则需要先按象限排序。

按实际情况的需要选择是否添加 $sign$ 函数。

```cpp
bool cmp1(vec x, vec y) {
	db cross = x * y;
	if (sign(cross) == 0) return x.length2() < y.length2();
	return sign(cross) > 0;
}
bool cmp2(vec x, vec y) {
	bool cx, cy;
	cx = sign(x.y) < 0 || (sign(x.y) == 0 && sign(x.x) < 0);
	cy = sign(y.y) < 0 || (sign(y.y) == 0 && sign(y.x) < 0);
	if (cx != cy) return cx < cy;
	db cross = x * y;
	if (sign(cross) == 0) return x.length2() < y.length2();
	return sign(cross) > 0;
}
```

### 跨立测试

如果线段 $P_1P_2$ 和 $Q_1Q_2$ 相交，首先两者和坐标轴平行的矩形（包围盒）必须相交，其次 $(\overrightarrow{Q_1 P_1} \times \overrightarrow{Q_1 Q_2})(\overrightarrow{Q_1 P_2} \times \overrightarrow{Q_1 Q_2}) \le 0$ 并且 $(\overrightarrow{P_1 Q_1} \times \overrightarrow{P_1 P_2})(\overrightarrow{P_1 Q_2} \times \overrightarrow{P_1 P_2}) \le 0$，说明两者相交。

如果是直线则无需判定包围盒，且只需判定一次叉积。

### 直线求交点

如果直线 $P_1P_2$ 和 $Q_1Q_2$ 相交，画出四个点连出的四边形（不是凸的也没事），其被直线 $Q_1Q_2$ 分成两个三角形，如果将 $P_1P_2$ 看做底，由于等高，其面积比为底边长度比，按比例乘即可得到交点。

```cpp
vec get_inter(vec P1, vec P2, vec Q1, vec Q2) {
	db ls = (Q2 - Q1) * (P1 - Q1);
	db rs = (Q2 - Q1) * (P2 - Q1);
	return P1 + (P2 - P1) * ls / (ls - rs);
}
```

### 点到线的距离

如果是直线，叉积得到面积除底边即可。

如果是线段，看与端点形成的夹角是否存在一个超过 $90 \degree$，如果是那就是到该端点为最短，采用叉积判断。否则就是直接同直线计算。

## 三角形相关

### 三角形的五心

记三个顶点 $A, B, C$，对应对边为 $a, b, c$。

#### 重心

##### 性质和计算

$$G = A + B + C$$

均质三角形的质心（几何中心）为重心。

三条中线的交点，是中线的靠近对边的三等分点。

三点距离的平方和最小。

##### 重心坐标

$$P = \alpha A + \beta B + \gamma C$$

其中 $\alpha = \frac{S_{\triangle PBC}}{S}$，以此类推。

三个参数和为 $1$。

在内部时 $\alpha > 0, \beta > 0, \gamma > 0$。

在边 $a$ 上 $\alpha = 0$，其余两 $> 0$。

在外面时，存在一个参数为负数。

可以利用叉积求解。

```cpp
void get_barycentric(vec A, vec B, vec C, vec P, 
	db &alpha, db &beta, db &gamma) {
	auto cross = [](vec A, vec B, vec C) -> db {
		return (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
	};
    db area_ABC = cross(A, B, C);
    alpha = cross(P, B, C) / area_ABC;
    beta  = cross(A, P, C) / area_ABC;
    gamma = 1.0 - alpha - beta;
}
```

#### 外心

重心坐标参数比例 $(\sin 2A, \sin 2B, \sin 2C)$。

三角形外接圆圆心，到三点距离相等。

$$r = \frac{abc}{4S}$$

锐角在内，直角在斜边中点，钝角在外。

```cpp
vec get_circle(vec a, vec b, vec c) {
	vec bc = c - b, ca = a - c, ab = b - a;
	return (b + c - bc.rotate90() * (ca / ab) / (ca * ab) ) / 2;
}
```

#### 内心

重心坐标参数比例 $(a, b, c)$。

$$I = \frac{aA+bB+cC}{a+b+c}$$

三条角平分线的交点，到三边距离相等且切于三边。

$$r = \frac{S}{p}$$ 

其中 $p$ 为半周长。

#### 垂心

重心坐标参数比例 $(\tan a, \tan b, \tan c)$。

三条垂线的交点。

由欧拉线 $OH = 3OG$ 可得:

$$H = 3G - 2O$$

#### 旁心

旁心是以一条边和另外两条边的延长线为切线的圆（旁切圆）的圆心。每个三角形有三个。

重心坐标参数比例 $(-a, b, c)$。

$$I_A = \frac{-aA+bB+cC}{-a+b+c}$$

$$r_a = \frac{S}{p - a}$$

### 正弦定理和余弦定理

$$\frac{a}{\sin A} = \frac{b}{\sin B} = \frac{c}{\sin C}$$

$$c^2 = a^2 + b^2 - 2ab\cos C$$

## 圆相关

### 求直线与圆的交点

如果直线与圆相离则无交点，若相切则可以利用垂线转化为求两直线交点。

若有两交点，先求出垂线长度，得到中点，然后加减半弦长得到交点。

### 求两圆交点

如果外离或内含则无交点。

如果相切，可以算出两圆心连线的方向向量，然后利用两圆半径计算出平移距离，最后将圆心沿这个方向向量进行平移即可。

如果两圆相交，则必有两个交点，并且关于两圆心连线对称。

设两个向量分别为 $\vec v_1, \vec v_2$，则有：

$$
|\vec v_1| = r_1 \\

|\vec v_2| = r_2 \\

|v_1 - v_2| = |O_1O_2|
$$

最底下式子改为：

$$(v_1 - v_2) ^2 = r_1^2 +  r_2^2 - 2r_1r_2 \cos \theta = |O_1O_2|^2$$

得：

$$\cos \theta = \frac{r_1^2+r_2^2-d^2}{2r_1r_2}$$

反解出 $\theta$，旋转连线即可。

### 求过点与圆的切线

根据 $d, r$ 算出夹角，旋转连线得到切线，然后得到垂线，然后平移 $r$ 即可。

### 求两圆公切线

把切线移动到一者的圆心，等效于一个点到一个半径 $r_1+r_2$ 的圆的切线（这里是内公切线，如果是外公切线那就是$r_1-r_2$）。求出后平移回去即可。

### 反演变换

#### 定义

用于处理圆线相切问题。

设反演中心为 $O$，反演半径为 $R$，有：

$$|OP||OP'|=R^2$$

转换写起来可能有点搞，画个图推一下半径和圆心等特殊信息就好。

#### 性质

> **保角性** 反演后两曲线交点的夹角不变，原本相切的依旧相切。

> **自反性** 对同一反演中心和反演半径反演两次可以得到原图形。

> **线圆性质** 过反演中心的圆反演后变为一条直线。

## 凸包与多边形相关

### 凸包

采用 Graham 算法。

```cpp
void get_convex(vector<vec> &s) {
	int n = s.size();
	for(int i = 1; i < n; ++i) {
		if(s[0].x > s[i].x || (s[0].x == s[i].x && s[0].y > s[i].y) ) {
			swap(s[0], s[i]);
		}
	}
	vec base = s[0];
	for(int i = 0; i < n; ++i) {
		s[i] = s[i] - base;
	}
	sort(s.begin() + 1, s.end(), [](const vec &x, const vec &y) {
		db res = x * y;
		if(res == 0) return x.length2() < y.length2();
		return res > 0;
	});
	vector<int> sta(s.size() + 3);
	int top;
	sta[top = 1] = 0;
	for(int i = 1; i < n; ++i) {
		while(top >= 2 && (s[i] - s[sta[top - 1] ]) * (s[sta[top] ] - s[sta[top - 1] ]) >= 0)
			--top;
		sta[++top] = i;
	}
	for(int i = 0; i < top; ++i) {
		s[i] = s[sta[i + 1] ] + base;
	} 
	s.resize(top);
}
```

### 闵可夫斯基和

多边形 $S, T$ 内部所有点相加得到的新多边形。

当为凸包相加时，等价为所有边加在一起极角排序，可以采用归并排序合并。

```cpp
void minkowsik(const vector<vec> &s1, const vector<vec> &s2, vector<vec> &t) {
	int n = s1.size(), m = s2.size();
	t.clear();
	t.push_back(s1[0] + s2[0]);
	vector<vec> d1(n), d2(m);
	for(int i = 0; i < n; ++i) d1[i] = s1[(i + 1) % n] - s1[i];
	for(int i = 0; i < m; ++i) d2[i] = s2[(i + 1) % m] - s2[i];
	int p1 = 0, p2 = 0;
	while(p1 < n && p2 < m) {
		t.push_back(t.back() + (d1[p1] * d2[p2] >= 0 ? d1[p1++] : d2[p2++]) );
	}
	while(p1 < n) t.push_back(t.back() + d1[p1++]);
	while(p2 < m) t.push_back(t.back() + d2[p2++]);
	t.pop_back();
}
```

### 查询点是否在凸包内

减去 base 之后在凸包上二分所在三角形即可。

```cpp
bool query(const vec &q, const vector<vec> &t) { // 均已经减去 base
	if(q * t[0] > 0 || t.back() * q > 0) return 0;
	auto it = lower_bound(t.begin(), t.end(), q, [](const vec &x, const vec &y) {
		db res = x * y;
		if(res == 0) return x.length2() < y.length2();
		return res > 0;
	});
	if(it == t.end() ) return 0;
	if(it == t.begin() ) return q.length2() <= t.begin()->length2();
	--it;
	vec lst = *it, nxt = *(it + 1);
	return (nxt - lst) * (q - lst) >= 0;
}
```

### 动态凸包

采用 set 维护即可：

```cpp
struct ConvexSet {
	db len = 0;
	set<vec> s;
	set<vec>::iterator pre(vec x){
		set<vec>::iterator p = s.lower_bound(x);
		if(p == s.begin()) p = s.end();
		--p;
		return p;
	}
	set<vec>::iterator nxt(vec x){
		set<vec>::iterator p = s.upper_bound(x);
		if(p == s.end()) return s.begin();
		return p;
	}
	void insert(vec x){
		if(s.size() < 3){
			s.insert(x);
			return;
		}
		if(query(x)) return;
		len -= dist(*pre(x), *nxt(x));
		set<vec>::iterator a, b;
		a = pre(x); b = pre(*a);
		while(s.size() >= 3 and (*b - x) * (*a - x) >= 0){
			s.erase(a);
			len -= dist(*a, *b);
			a = b; b = pre(*a);
		}
		len += dist(*pre(x), x);
		a = nxt(x); b = nxt(*a);
		while(s.size() >= 3 and (*b - x) * (*a - x) <= 0){
			s.erase(a);
			len -= dist(*a, *b);
			a = b; b = nxt(*a);
		}
		len += dist(*nxt(x), x);
		s.insert(x);
	}
	bool query(vec x){
		if(s.count(x)) return true;
		set<vec>::iterator pr, nx;
		pr = pre(x); nx = nxt(x);
		if((*nx - *pr) * (x - *pr) > 0) return false;
		return true;
	}
};
```

### pick 定理

$$A=i + \frac{b}{2} - 1$$

其中 $A$ 为面积，$i$ 为内部整点数，$b$ 为边界整点数。

本质和欧拉公式：

$$V - E + F = 2$$ 

等价，其中 $V, E, F$ 分别为顶点数，边数，和面数。

### 三角剖分

以任意一点 $P$ 为剖分中心，然后依次连向多边形上相邻的两点，以叉积正负标定三角形方向。

如果希望叉积都是正的，剖分中心需要选在多边形的核内部。

#### 求多边形面积

直接以原点作为剖分中心，叉积直接计算有向面积即可。

```cpp
db get_area(const vector<vec> &s) {
	db res = 0;
	int n = s.size();
	if (n < 3) return 0; 
	for (int i = 0; i < n; ++i) {
		res += s[i] * s[(i + 1) % n];
	}
	return res / 2.0;
}
```

#### 求多边形重心

三角剖分后，对于每个三角形算出重心（几何中心），然后根据面积大小为权做加权平均即可。类似求面积。

### 凸壳

按 $x$ 为第一关键字，$y$ 为第二关键字，排序后依次加入，叉积弹出凹的点即可。

正做一遍得到上凸壳，反做一遍得到下凸壳。两个拼起来做可以得到凸包。

本质是 Andrew 算法。

### 李超线段树

每次把中点处劣的线段下传，如果这个点没有永久化标记那就停止，直接把线段插入。

然后具体下传到哪边，就是哪边的端点处我比原标记大我下哪边。复杂度 $O(n\log^2 n)$

如果要做线段树合并，每个点合并完还是得把劣的下放。复杂度 $O(n\log n)$

```cpp
const db eps = 1e-9;
int sign(db x) { return x > eps ? 1 : x > -eps ? 0 : -1; }
struct SegTree {
	struct line {
		db k, b;
		line(int sx, int sy, int fx, int fy) {
			if (sx == fx) {
				k = 0;
				b = max(sy, fy);
				return;
			}
			k = 1.0 * (fy - sy) / (fx - sx);
			b = sy - sx * k;
		}
		db calc(int x) { return k * x + b; }
	};
	int V;
	vector<line> p;
	vector<int> id; 
	SegTree(int n, int _V) : V(_V), id(V * 4 + 3) {
		p.reserve(n + 1);
		p.push_back(line(0, -1e9, 0, -1e9));
	}
	bool cmp(int x, int y, int k) {
		return 
			p[x].calc(k) < p[y].calc(k) ||
			(sign(p[x].calc(k) - p[y].calc(k)) == 0 && x > y) ||
			x == 0;
	}
	void modify_ins(int x, int l, int r, int k) {
		if (!id[x]) {
			id[x] = k;
			return;
		}
		int mid = l + r >> 1;
		if (cmp(id[x], k, mid)) swap(id[x], k);
		if (l == r) {
			return;
		}
		bool fl = cmp(id[x], k, l), fr = cmp(id[x], k, r);
		if (fl) modify_ins(x << 1, l, mid, k);
		if (fr) modify_ins(x << 1 | 1, mid + 1, r, k);
	}
	void modify_insdiv(int x, int l, int r, int L, int R, int k) {
		if (L <= l && r <= R) {
			modify_ins(x, l, r, k);
			return;
		};
		int mid = l + r >> 1;
		if (L <= mid) modify_insdiv(x << 1, l, mid, L, R, k);
		if (mid + 1 <= R) modify_insdiv(x << 1 | 1, mid + 1, r, L, R, k);
	}
	int query(int x, int l, int r, int k) {
		int ans = id[x];
		int mid = l + r >> 1;
		if (l == r) return ans;
		int res = 0;
		if (k <= mid) {
			res = query(x << 1, l, mid, k);
		} else {
			res = query(x << 1 | 1, mid + 1, r, k);
		}
		if (cmp(ans, res, k)) ans = res;
		return ans;
	}
	void insert(int sx, int sy, int fx, int fy) {
		p.emplace_back(sx, sy, fx, fy);
		modify_insdiv(1, 1, V, min(sx, fx), max(sx, fx), p.size() - 1);
	}
	int query(int x) {
		return query(1, 1, V, x);
	}
};
```

### 旋转卡壳

枚举边，再维护一个指针找到距离这个边最远的点即可。

如果还需要左右的垂线，则选取点积最大和最小的两个即可。（投影长度最长）

### 平面最近点对

#### KDT 做法

点插入 KD tree，在上面 dfs 剪枝求解，复杂度不稳定，最坏是 $O(k n ^ 2)$ 的。

#### 分治做法

二维做法，总体先按 $x$ 排序，然后分治，按 $y$ 归并，假设已经求出底下分治的答案 $d$，每次取出与中间分界距离小于 $d$ 的，然后这些点中，枚举每一个点 $i$，要求 $j$ 和 $i$ 的 $y$ 上的距离小于 $d$，得到左右各一个矩形，然后只在这个矩形里暴力枚举即可。因此也可以做其他的一些最小化的内容。

每个 $d\times d$ 的矩形分成四个 $\frac{d}{2} \times \frac{d}{2} $ 的小矩形，可以知道每个小矩形里只有一个元素，所以整个区域（左右两个矩形）只有 $8$ 个点，暴力复杂度正确。

多维需要反复的取这个中间的平面，实现降维。但是这个点数增长很快，大概是 $O(k^3)$ 的，（就是三维就有 $27$ 个点了），所以高维常数有点大。

#### 随机划分网格做法

期望复杂度 $O(n)$。

把点随机打乱，然后依次加入。如果当前的答案为 $d$，把整个空间划分为 $\frac{d}{\sqrt k}$ 的网格，每个格子里只有一个点，对于新加入的点查询周围的网格即可。网格中的点用 map 之类的东西存储。

### 半平面交

#### 多边形的核

如果一个点集中的点与多边形上任意一点的连线与多边形没有其他交点，那么这个点集被称为多边形的核。

把多边形的每条边看成是首尾相连的向量，那么这些向量在多边形内部方向的半平面交就是多边形的核。

#### 单调队列求解半平面交

极角排序后，新加入的线可能会弹出队首和队尾的元素，判断条件是看上一个交点在这个线的内侧还是外侧，记录交点然后叉积判断即可。

结束时用队首重新把队尾多余的线弹出，因为我们之前只用新线弹出过队首，但是没用队首检查新加入的线。（相当于闭合一下）

注意需要先删尾后删头，比如只有两条线时，交点位于新线外侧，按理说应该先弹出尾部的。这是因为这两个线共用一个交点了。平行线最好做一下去重，只保留最左侧的。

复杂度 $O(n \log n)$（瓶颈是排序）。

```cpp
vector<vec> HPI(vector<line> s) {
	sort(s.begin(), s.end(), [](line lx, line ly) {
		vec x = lx.b - lx.a;
		vec y = ly.b - ly.a;
		bool cx, cy;
		cx = sign(x.y) < 0 || (sign(x.y) == 0 && sign(x.x) < 0);
		cy = sign(y.y) < 0 || (sign(y.y) == 0 && sign(y.x) < 0);
		if (cx != cy) return cx < cy;
		db cross = x * y;
		if (sign(cross) == 0) return sign((lx.b - lx.a) * (ly.b - lx.a)) < 0;
		return sign(cross) > 0;
	});

	int m = 0;
	for (int i = 0; i < s.size(); ++i) {
		if (i > 0) {
			vec v1 = s[i].b - s[i].a;
			vec v2 = s[i - 1].b - s[i - 1].a;
			if (sign(v1 * v2) == 0 && sign(v1 / v2) > 0) {
				continue;
			}
		}
		s[m++] = s[i];
	}
	s.resize(m);

	vector<vec> t(s.size() + 1);
	vector<line> q(s.size());
	int l = 0, r = -1;
	for (int i = 0; i < s.size(); ++i) {
		while (r - l + 1 > 1 && sign((s[i].b - s[i].a) * (t[r] - s[i].a)) <= 0) --r;
		while (r - l + 1 > 1 && sign((s[i].b - s[i].a) * (t[l + 1] - s[i].a)) <= 0) ++l;
		q[++r] = s[i];
		if (r - l + 1 > 1) t[r] = get_inter(q[r], q[r - 1]);
	}
	while (r - l + 1 > 1 && sign((q[l].b - q[l].a) * (t[r] - q[l].a)) <= 0) --r;
	t[r + 1] = get_inter(q[l], q[r]);
	for (int i = l + 1; i <= r + 1; ++i) {
		t[i - l - 1] = t[i];
	}
	t.resize(r + 1 - l);
	return t;
}
```

#### 暴力拆多边形求解半平面交

给定多边形顶点，枚举每条边，在线内侧的保留，跨立的边加入交点即可。

要求原本就有一个多边形（有封闭的顶点），不然有点麻烦。如果要的就是切割多边形那种效果，只能使用这种，这种严格意义上不是半平面交。否则一般由于复杂度比较劣也不考虑采用。

若对一个 $n$ 个点的多边形进行 $k$ 次切割，复杂度 $O(k n)$。

```cpp
vector<vec> HPI(const vector<vec> &p, const vector<line> &sp) {
	vector<vec> s = p;
	vector<vec> res;
	res.reserve(s.size() + p.size());
	for (line l : sp) {
		res.clear();
		for (int i = 0; i < s.size(); ++i) {
			vec a = s[i], b = s[(i + 1) % s.size()];
			int x = sign((a - l.a) * (l.b - l.a));
			int y = sign((b - l.a) * (l.b - l.a));
			if (x < 0) res.push_back(a);
			if (x * y < 0) res.push_back(get_inter(line(a, b), l));
		}
		s = res;
	}
	return s;
}
```

