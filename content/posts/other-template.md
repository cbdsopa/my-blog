+++
date = '2026-09-05T19:00:03+08:00'
draft = false
title = '算法竞赛杂项模板'

tags = ["自适应辛普森积分", "Pollard-rho", "素性测试"]
+++

# 杂项

## 自适应辛普森积分

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

## 素性测试和 Pollard-rho

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
