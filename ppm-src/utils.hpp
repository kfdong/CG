#ifndef UTILS_H
#define UTILS_H

#include "config.h"
#include <bits/stdc++.h>
using std::min;
using std::max;
using std::vector;
using std::pair;

template <class T> T sqr(T x) { return x * x; }
template <class T> T pow5(T x) { return sqr(sqr(x)) * x; }
int sign(real_t x) { return x < 0 ? -1 : (x > 0 ? 1 : 0); }
real_t clip(real_t x, real_t l, real_t r) {
	return x < l ? l : (x > r ? r : x);
}

class vect {
public:
	real_t e[3];
	vect() {}
	vect(real_t e0, real_t e1, real_t e2) { e[0] = e0, e[1] = e1, e[2] = e2; }
	vect(real_t *_e) { e[0] = _e[0], e[1] = _e[1], e[2] = _e[2]; }

	inline real_t operator [] (const int &p) const { return e[p]; }
	inline real_t& operator [] (const int &p) { return e[p]; }

    inline vect& operator += (const vect &v); 
    inline vect& operator -= (const vect &v); 
    inline vect& operator *= (const vect &v); 
    inline vect& operator *= (const real_t &v); 
    inline vect& operator /= (const real_t &v); 
    // inline vect& operator /= (const vect &v); 
    inline vect operator + (const vect &v) const;
    inline vect operator - (const vect &v) const;
    inline vect operator * (const vect &v) const; 
    inline vect operator * (const real_t &v) const;  
    // inline vect operator / (const vect &v) const;  
    inline vect operator / (const real_t &v) const; 

	inline real_t abs() const;
	inline real_t abs2() const;
	inline vect norm() const;

	void read(FILE *f = stdin) {
		fscanf(f, "%lf%lf%lf", &e[0], &e[1], &e[2]);
	}
};

class ray {
	/* 
	 * direction vector v is normalized;
	 */
public:
	vect s, v;
	ray() {}
	ray(vect _s, vect _v) { s = _s, v = _v.norm(); }
	inline vect operator [] (const real_t &x) const;
};

inline real_t vect::abs2() const {
	return sqr(e[0]) + sqr(e[1]) + sqr(e[2]);
}
inline real_t vect::abs() const {
	return sqrt(abs2());
}
inline vect cross(const vect &a, const vect &b) {
	return vect(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
inline real_t dot(const vect &a, const vect &b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline vect vect::norm() const {
	real_t r = 1.0 / abs();
	return *this * r;
}

inline vect& vect::operator += (const vect &v) {
	e[0] += v[0], e[1] += v[1], e[2] += v[2];
	return *this;
}

inline vect& vect::operator -= (const vect &v) {
	e[0] -= v[0], e[1] -= v[1], e[2] -= v[2];
	return *this;
}

inline vect& vect::operator *= (const vect &v) {
	e[0] *= v[0], e[1] *= v[1], e[2] *= v[2];
	return *this;
}

inline vect& vect::operator *= (const real_t &v) {
	e[0] *= v, e[1] *= v, e[2] *= v;
	return *this;
}

inline vect& vect::operator /= (const real_t &v) {
	e[0] /= v, e[1] /= v, e[2] /= v;
	return *this;
}

inline vect vect::operator + (const vect &v) const {
	return vect(e[0] + v[0], e[1] + v[1], e[2] + v[2]);
}

inline vect vect::operator - (const vect &v) const {
	return vect(e[0] - v[0], e[1] - v[1], e[2] - v[2]);
}

inline vect vect::operator * (const vect &v) const {
	return vect(e[0] * v[0], e[1] * v[1], e[2] * v[2]);
}

inline vect vect::operator * (const real_t &v) const {
	return vect(e[0] * v, e[1] * v, e[2] * v);
}

inline vect operator * (const real_t &x, const vect &v) {
	return vect(v[0] * x, v[1] * x, v[2] * x); 
}

inline vect vect::operator / (const real_t &v) const {
	return vect(e[0] / v, e[1] / v, e[2] / v);
}

inline std::istream& operator >> (std::istream &is, vect &t) {
	is >> t[0] >> t[1] >> t[2];
	return is;
}

inline std::ostream& operator << (std::ostream &os, const vect &t) {
	os << '(' << t[0] << ", " << t[1] << ", " << t[2] << ')';
	return os;
}

inline vect ray::operator[] (const real_t &x) const {
	return s + x * v;
}

vect random_in_unit_ball() {
	vect ret;
	do {
		ret = 2.0 * vect(rand_unit(gen), rand_unit(gen), rand_unit(gen)) - vect(1, 1, 1);
	} while (ret.abs2() >= 1.0);
	return ret;
}

vect random_in_unit_disk() {
	vect ret;
	do {
		ret = 2.0 * vect(rand_unit(gen), rand_unit(gen), 0) - vect(1, 1, 0);
	} while (ret.abs2() >= 1.0);
	return ret;
}

vect reflect(const vect &v, const vect &norm) {
	return v - 2 * dot(v, norm) * norm;
}

class Progress_bar {
	int total, current;
	int barWidth = 70;
	std::string header;
public:

	Progress_bar() {}
	Progress_bar(std::string s) { header = s; }
	void set_total(int _total) {
		current = 0;
		total = _total;
	}
	void update(int part) {
		current += part;
		real_t progress = std::min((real_t)current / total, (real_t)1.0);

		if (header.size()) printf("%s: ", header.c_str());
		putchar('[');
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) putchar('=');
			else if (i == pos) putchar('>');
			else putchar(' ');
		}
		printf("] %d %c", int(progress * 100.0), '%');
		putchar('\r');
		fflush(stdout);
	}
	void finish() {
		puts("");
	}
};

/* ppm */
class recording_point {
public:
	vect x; 			/* position */
	vect attenuation; 	/* color attenuation */
	vect norm;
	vect flux;
	real_t rad2; 		/* square of radius */
	int n; 				/* N / alpha */
	int pixx, pixy; 	/* pixel id */
};


struct kd_tree {
	kd_tree *son[2];
	recording_point *key;
	vect xmin, xmax;
	int dim;
};

vect min(const vect &a, const vect &b) {
	return vect(min(a.e[0], b.e[0]), min(a.e[1], b.e[1]), min(a.e[2], b.e[2]));
}
vect max(const vect &a, const vect &b) {
	return vect(max(a.e[0], b.e[0]), max(a.e[1], b.e[1]), max(a.e[2], b.e[2]));
}

int kd_tree_dim_id;
class KDTree {
public:
	kd_tree *root;
	vector <recording_point*> h;
	kd_tree **Q;
	real_t r2;

	static bool kd_tree_compare(recording_point *a, recording_point *b) {
		return a->x[kd_tree_dim_id] < b->x[kd_tree_dim_id];
	}

	kd_tree *build_tree(int l, int r, int dim) {
		if (l > r) return NULL;
		kd_tree *now = new kd_tree;
		now->dim = dim;
		now->key = h[l];
		now->xmin = now->xmax = now->key->x;
		kd_tree_dim_id = dim;
		int mid = (l + r) >> 1;
		std::nth_element(h.begin() + l, h.begin() + mid, h.begin() + r + 1, kd_tree_compare);
		now->son[0] = build_tree(l, mid - 1, (dim + 1) % 3);
		now->son[1] = build_tree(mid + 1, r, (dim + 1) % 3);
		if (now->son[0]) {
			now->xmin = min(now->xmin, now->son[0]->xmin);
			now->xmax = max(now->xmax, now->son[0]->xmax);
		}
		if (now->son[1]) {
			now->xmin = min(now->xmin, now->son[1]->xmin);
			now->xmax = max(now->xmax, now->son[1]->xmax);
		}
		now->key = h[mid];
		return now;
	}
	void build(const vector <recording_point*> &_h, real_t _r2) {
		h = _h;
		root = build_tree(0, (int)h.size() - 1, 0);
		Q = new kd_tree*[h.size()];
		r2 = _r2;
	}

	vector <recording_point* > rec;
	vector <recording_point* >* query(vect x) {
		rec.clear();
		int l = 0, r = 0;
		Q[r++] = root;
		while (l != r) {
			kd_tree *now = Q[l++];

			real_t dist2 = 0;
			dist2 += sqr(x.e[0] > now->xmax.e[0] ? x.e[0] - now->xmax.e[0] : 
						 (x.e[0] < now->xmin.e[0] ? now->xmin.e[0] - x.e[0] : 0));
			dist2 += sqr(x.e[1] > now->xmax.e[1] ? x.e[1] - now->xmax.e[1] : 
						 (x.e[1] < now->xmin.e[1] ? now->xmin.e[1] - x.e[1] : 0));
			dist2 += sqr(x.e[2] > now->xmax.e[2] ? x.e[2] - now->xmax.e[2] : 
						 (x.e[2] < now->xmin.e[2] ? now->xmin.e[2] - x.e[2] : 0));

			if (dist2 >= r2) continue;
			if ((now->key->x - x).abs2() < r2) rec.push_back(now->key);
			if (now->son[0] != NULL) Q[r++] = now->son[0];
			if (now->son[1] != NULL) Q[r++] = now->son[1];
		}
		return &rec;
	}
};

int cmod(int x, int v) {
	x = x % v;
	return x < 0 ? x + v : x;
}
struct hash_entry {
	long long key;
	recording_point *value;
	bool operator < (const hash_entry &x) const { return key < x.key; }
};
class hash_table {
public:
	static const int mod = 8000009;
	static const long long base = 12345703, base2 = base * base;
	hash_entry *hash;
	vector <pair<long long, vector <recording_point*>* > > start[mod];
	real_t r;
	int n;

	long long get_key(real_t x, real_t y, real_t z) {
		return (long long)floor(x / r) * base2 + (long long)floor(y / r) * base + (long long)floor(z / r);
	}

	void build(const vector <recording_point*> &h, real_t _r2) {
		r = 2 * sqrt(_r2) + EPS;
		// duplicate input for hash_grid
		hash = new hash_entry[h.size() * 8];
		n = h.size();
		for (int i = 0; i < n; ++i) {
			real_t dx = (fmod(fmod(h[i]->x.e[0], r) + r, r) < r / 2 ? -r : r);
			real_t dy = (fmod(fmod(h[i]->x.e[1], r) + r, r) < r / 2 ? -r : r);
			real_t dz = (fmod(fmod(h[i]->x.e[2], r) + r, r) < r / 2 ? -r : r);
			hash[i * 8 + 0].key = get_key(h[i]->x.e[0], h[i]->x.e[1], h[i]->x.e[2]);
			hash[i * 8 + 1].key = get_key(h[i]->x.e[0], h[i]->x.e[1], h[i]->x.e[2] + dz);
			hash[i * 8 + 2].key = get_key(h[i]->x.e[0], h[i]->x.e[1] + dy, h[i]->x.e[2]);
			hash[i * 8 + 3].key = get_key(h[i]->x.e[0], h[i]->x.e[1] + dy, h[i]->x.e[2] + dz);
			hash[i * 8 + 4].key = get_key(h[i]->x.e[0] + dx, h[i]->x.e[1], h[i]->x.e[2]);
			hash[i * 8 + 5].key = get_key(h[i]->x.e[0] + dx, h[i]->x.e[1], h[i]->x.e[2] + dz);
			hash[i * 8 + 6].key = get_key(h[i]->x.e[0] + dx, h[i]->x.e[1] + dy, h[i]->x.e[2]);
			hash[i * 8 + 7].key = get_key(h[i]->x.e[0] + dx, h[i]->x.e[1] + dy, h[i]->x.e[2] + dz);
			for (int j = i * 8; j < i * 8 + 8; ++j)
				hash[j].value = h[i];
		}
		n = n * 8;
		std::sort(hash, hash + n);

		vector <recording_point*> *rec = NULL;
		for (int i = 0; i < n; ++i) {
			if (i == 0 || hash[i].key != hash[i - 1].key) 
				rec = new vector <recording_point*>;
			rec->push_back(hash[i].value);
			if (i + 1 == n || hash[i + 1].key != hash[i].key)
				start[cmod(hash[i].key, mod)].push_back(std::make_pair(hash[i].key, rec));
		}
	}

	vector <recording_point*>* query(const vect &x) {
		long long key = get_key(x.e[0], x.e[1], x.e[2]);
		for (auto p : start[cmod(key, mod)])
			if (p.first == key) return p.second;
		return NULL;
	}

	~hash_table() {
		delete[] hash;
	}
};

class Bezier { /* Bezier curve */
public:
	static const int MAXN = 7;
	int n; /* 0 .. n */
	real_t cox[MAXN], coy[MAXN];
	real_t cox2[MAXN * 2], coy2[MAXN * 2];
	real_t fac[MAXN];

	Bezier() {}
	Bezier(int _n, real_t *px, real_t *py) {
		n = _n;
		fac[0] = 1;
		for (int i = 1; i <= n; ++i) fac[i] = fac[i - 1] * i;
		memset(cox, 0, sizeof(cox));
		memset(coy, 0, sizeof(coy));
		for (int i = 0; i <= n; ++i) {
			real_t c = binomial(n, i);
			for (int j = 0; j <= n - i; ++j) {
				cox[i + j] += c * px[i] * binomial(n - i, j) * (j % 2 == 0 ? 1 : -1);
				coy[i + j] += c * py[i] * binomial(n - i, j) * (j % 2 == 0 ? 1 : -1);
			}
		}
		memset(cox2, 0, sizeof(cox2));
		memset(coy2, 0, sizeof(coy2));
		for (int i = 0; i <= n; ++i)
			for (int j = 0; j <= n; ++j) {
				cox2[i + j] += cox[i] * cox[j];
				coy2[i + j] += coy[i] * coy[j];
			}
	}

	real_t binomial(int n, int m) { return fac[n] / (fac[m] * fac[n - m]); }
	void eval(real_t t, real_t &x, real_t &y) {
		x = y = 0;
		real_t pw = 1;
		for (int i = 0; i <= n; ++i, pw *= t) {
			x += cox[i] * pw;
			y += coy[i] * pw;
		}
	}
	void eval_deriv(real_t t, real_t &x, real_t &y) {
		x = y = 0;
		real_t pw = 1;
		for (int i = 1; i <= n; ++i, pw *= t) {
			x += cox[i] * pw * i;
			y += coy[i] * pw * i;
		}
	}
	void eval_jacob(real_t t, real_t &x, real_t &y) {
		x = y = 0;
		real_t pw = 1;
		for (int i = 2; i <= n; ++i, pw *= t) {
			x += cox[i] * pw * i * (i - 1);
			y += coy[i] * pw * i * (i - 1);
		}
	}
};

#endif
