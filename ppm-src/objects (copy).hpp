#ifndef OBJECTS_H
#define OBJECTS_H

#include "config.h"
#include "utils.hpp"

class material;
class hit_record;
class bounding_box;

enum material_type {LAMBERTIAN, METAL, GLASS};
class material {
public:
	virtual bool scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const = 0;
	material_type type;
};

class hit_record {
public:
	bool hit;
	real_t t;
	vect p;
	vect norm; // ! normalized 
	material *mat;

	bool operator < (const hit_record &h) const {
		if (hit && h.hit) return t < h.t;
		else return hit;
	}

	hit_record() { hit = false; }
};

class object {
public:
	/* directed face, only intersect when dot(norm, r.v) > 0. */
	material *mat;
	virtual hit_record hit(const ray &r) const { 
		/* the first intersection from r.s along r.v */
		return hit_record();
	}
};

class sphere : public object { /* directed sphere */
public:
	vect c;
	real_t radius;

	sphere() {}
	sphere(vect _c, real_t _r, material *_mat) : c(_c), radius(_r) {
		mat = _mat;
	}
	virtual hit_record hit(const ray &r) const;
};

hit_record sphere::hit(const ray &r) const {
	hit_record ret;
	real_t dist2 = (c - r.s).abs2();
	real_t dotp = dot(c - r.s, r.v);
	real_t l2 = dist2 - sqr(dotp);
	real_t delta = sqr(radius) - l2;
	if (delta > 0) {
		real_t tp = sqrt(delta);
		real_t tmp = dotp - sign(radius) * tp;
		if (tmp >= EPS) { // EPS to avoid self-intersection
			ret.hit = true;
			ret.t = tmp;
			ret.p = r[tmp];
			ret.norm = (sign(radius) * (ret.p - c)).norm();
			ret.mat = mat;
		}
	}
	return ret;
}

class hemisphere : public object { /* directed sphere */
public:
	vect c;
	real_t radius;
	vect norm; /* only for dot(p - c, norm) > 0 */

	hemisphere() {}
	hemisphere(vect _c, real_t _r, vect _norm, material *_mat) : c(_c), radius(_r), norm(_norm) {
		mat = _mat;
	}
	virtual hit_record hit(const ray &r) const;
};

hit_record hemisphere::hit(const ray &r) const {
	hit_record ret;
	real_t dist2 = (c - r.s).abs2();
	real_t dotp = dot(c - r.s, r.v);
	real_t l2 = dist2 - sqr(dotp);
	real_t delta = sqr(radius) - l2;
	if (delta > 0) {
		real_t tp = sqrt(delta);
		real_t tmp = dotp - sign(radius) * tp;
		if (tmp >= EPS && dot(r[tmp] - c, norm) > 0) { // EPS to avoid self-intersection
			ret.hit = true;
			ret.t = tmp;
			ret.p = r[tmp];
			ret.norm = (sign(radius) * (ret.p - c)).norm();
			ret.mat = mat;
		}
	}
	return ret;
}

class triangle : public object { /* directed triangle */
public:
	vect a, b, c; /* a -> b -> c, righthand direction */
	vect norm;

	triangle() {}
	triangle(vect _a, vect _b, vect _c, material *_mat) {
		a = _a, b = _b, c = _c;
		mat = _mat;
		norm = cross(b - a, c - a).norm();
	}
	virtual hit_record hit(const ray &r) const;
};

hit_record triangle::hit(const ray &r) const {
	hit_record ret;
	real_t dt = dot(r.v, norm);
	if (fabs(dt) > EPS) {
		real_t t = -dot(r.s - a, norm) / dt;
		if (t >= EPS) {
			vect p = r[t];
			if (dot(cross(b - a, p - a), norm) >= 0 && dot(cross(c - b, p - b), norm) >= 0
					&& dot(cross(a - c, p - c), norm) >= 0) {
				ret.p = p;
				ret.hit = true;
				ret.t = t;
				ret.norm = norm;
				ret.mat = mat;
			}
		}
	}
	return ret;
}

class polygon : public object { /* directed triangle */
public:
	vector <vect> v; /* v[0], v[1], ..., v[n - 1] */
	int n;
	vect norm;

	polygon() {}
	polygon(vector <vect> _v, material *_mat) {
		v = _v;
		n = _v.size();
		assert(n >= 3);
		mat = _mat;
		norm = cross(v[2] - v[1], v[0] - v[1]).norm();
	}
	virtual hit_record hit(const ray &r) const;
};

hit_record polygon::hit(const ray &r) const {
	hit_record ret;
	real_t dt = dot(r.v, norm);
	if (fabs(dt) > EPS) {
		real_t t = -dot(r.s - v[0], norm) / dt;
		if (t >= EPS) {
			vect p = r[t];
			bool flag = 1;
			for (int i = 0; i < n; ++i) {
				int j = (i + 1) % n;
				if (dot(cross(v[j] - v[i], p - v[i]), norm) < 0) {
					flag = 0;
					break;
				}
			}
			if (flag) {
				ret.p = p;
				ret.hit = true;
				ret.t = t;
				ret.norm = norm;
				ret.mat = mat;
			}
		}
	}
	return ret;
}

class object_list : public object {
public:
	std::vector < object* > olist;

	object_list() {}

	size_t size() const { return olist.size(); }
	void add(object &o);
	virtual hit_record hit(const ray &r) const;
};

void object_list::add(object &o) {
	olist.push_back(&o);
}

hit_record object_list::hit(const ray &r) const {
	hit_record ret;
	for (auto o : olist) 
		ret = std::min(ret, o->hit(r));
	return ret;
}

/* material */
class lambertian : public material {
public:
	vect albedo;
	lambertian(const vect &_albedo) : albedo(_albedo) {
		type = LAMBERTIAN;
	}
	bool scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const; 
};

bool lambertian::scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const {
	vect v = h.norm + random_in_unit_ball();
	out = ray(h.p, v);
	attenuation = albedo;
	return true;
}

class metal : public material {
public:
	vect albedo;
	real_t fuzzy; // 0 <= fuzzy <= 1
	metal(const vect &_albedo, real_t fuzzy = 0) : albedo(_albedo), fuzzy(fuzzy) {
		type = METAL;
	}
	bool scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const; 
};

bool metal::scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const {
	vect v = reflect(in.v, h.norm) + fuzzy * random_in_unit_ball();
	out = ray(h.p, v);
	attenuation = albedo;
	return true;
}

class glass : public material {
public:
	real_t refractive_in, refractive_out; /* in : incoming light, out : outcoming light */
	real_t obsorb;
	glass(real_t _refractive_in, real_t _refractive_out, real_t _obsorb = 1.0) : 
		refractive_in(_refractive_in), refractive_out(_refractive_out), obsorb(_obsorb) {
		type = GLASS;
	}
	bool scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const; 
	bool scatter(const ray &in, const hit_record &h, vect &attenuation, ray &refrac, ray &reflac, real_t &refrac_p) const; 
};

// real_t schlick_approximation(real_t cos, real_t n1, real_t n2) {
	// real_t r0 = sqr((n1 - n2) / (n1 + n2));
	// return r0 + (1 - r0) * pow5(1 - cos);

real_t schlick_approximation(real_t cos, real_t n2, real_t n1) {
	assert(cos >= 0);
	real_t sin = sqrt(1.0 - cos * cos);
	real_t ct = sqrt(1 - sqr(n1 / n2 * sin));
	real_t rs = sqr((n1 * cos - n2 * ct) / (n1 * cos + n2 * ct));
	real_t rp = sqr((n1 * ct - n2 * cos) / (n1 * ct + n2 * cos));
	return (rs + rp) * 0.5;
}

bool refract(vect v, const vect &norm, real_t n1, real_t n2, vect &out) {
	v = v.norm();
	real_t tmp = n2 / n1;
	real_t cos = -dot(v, norm);
	real_t delta = 1 - sqr(tmp) * (1 - sqr(cos));
	if (delta <= 0) return false;
	out = tmp * v + (cos * tmp - sqrt(delta)) * norm;
	return true;
}

bool glass::scatter(const ray &in, const hit_record &h, vect &attenuation, ray &out) const {
	vect v_refrac;
	vect v_reflac = reflect(in.v, h.norm);
	if (refract(in.v, h.norm, refractive_in, refractive_out, v_refrac)) { 
		real_t cos = -dot(in.v, h.norm);
		real_t refrac_prob = 1 - schlick_approximation(cos, refractive_in, refractive_out);
		if (rand_unit(gen) < refrac_prob)
			out = ray(h.p, v_refrac);
		else 
			out = ray(h.p, v_reflac);
	} else {
		out = ray(h.p, v_reflac);
	}
	attenuation = vect(1, 1, 1) * obsorb;
	return true;
}

bool glass::scatter(const ray &in, const hit_record &h, vect &attenuation, ray &refrac, ray &reflac, real_t &refrac_p) const {
	vect v_refrac;
	vect v_reflac = reflect(in.v, h.norm);
	if (refract(in.v, h.norm, refractive_in, refractive_out, v_refrac)) { 
		real_t cos = -dot(in.v, h.norm);
		real_t refrac_prob = 1 - schlick_approximation(cos, refractive_in, refractive_out);
		refrac_p = refrac_prob;
		refrac = ray(h.p, v_refrac);
		reflac = ray(h.p, v_reflac);
	} else {
		refrac_p = 0;
		reflac = ray(h.p, v_reflac);
	}
	attenuation = vect(1, 1, 1) * obsorb;
	return true;
}

/* lights */
class lights {
public:
	virtual ray emit(vect &col_in) { assert(false); }
};

class point_light : public lights {
public:
	vect x; /* position */
	vect norm; /* outgoing norm of the hemisphere */
	vect color; /* color of the light */

	point_light() {}
	point_light(const vect &_x, const vect &_norm, const vect &_color) {
		x = _x;
		norm = _norm.norm();
		color = _color;
	}

	ray emit(vect &col_in) {
		col_in = color;
		ray ret(x, norm + random_in_unit_ball());
		return ret;
	}
};

class hemisphere_light: public lights {
public:
	vect x; /* position */
	vect norm; /* outgoing norm of the hemisphere */
	vect color; /* color of the light */

	hemisphere_light() {}
	hemisphere_light(const vect &_x, const vect &_norm, const vect &_color) {
		x = _x;
		norm = _norm.norm();
		color = _color;
	}

	ray emit(vect &col_in) {
		col_in = color;
		vect v = random_in_unit_ball();
		if (dot(v, norm) < 0) v -= 2 * norm * dot(v, norm);
		ray ret(x, v);
		return ret;
	}
};

class bounding_box {
public:
	object *obj;
	bool hit(const ray &r) const {
		if ((obj->hit(r)).hit) return true;
		else return false;
	}
};

/* rotating objects along y-axis */
class rotation : public object {
public:
	real_t dir; /* when dir = 1, norm direction is from rotation axis to Bezier function */
	vect axis;
	Bezier *f;
	vector <real_t> iy, ix, it;
	bounding_box *bd;

	rotation() {}
	rotation(real_t _dir, vect _axis, Bezier *_f, material *_mat) : dir(_dir), axis(_axis), f(_f) {
		mat = _mat;
		for (int i = 0; i <= BEZIER_SAMPLE; ++i) {
			real_t x, y;
			f->eval(1.0 * i / BEZIER_SAMPLE, x, y);
			iy.push_back(y);
			ix.push_back(x);
			it.push_back((i + 0.5) / BEZIER_SAMPLE);
		}
	}
	virtual hit_record hit(const ray &r) const;
	void get_init(const ray &r, vector <real_t> &init) const;
};

void rotation::get_init(const ray &r, vector <real_t> &init) const {
	real_t pa = sqr(r.v.e[0]) + sqr(r.v.e[2]);
	real_t pb = 2 * (r.v.e[0] * (r.s.e[0] - axis.e[0]) + r.v.e[2] * (r.s.e[2] - axis.e[2]));
	real_t pc = sqr(r.s.e[0] - axis.e[0]) + sqr(r.s.e[2] - axis.e[2]);
	real_t a = pa / sqr(r.v.e[1]);
	real_t b = pb / r.v.e[1] - 2 * pa * r.s.e[1] / sqr(r.v.e[1]);
	real_t c = pa * sqr(r.s.e[1] / r.v.e[1]) - pb * r.s.e[1] / r.v.e[1] + pc;

	real_t dist = -b * b / (4 * a) + c;
	if (fabs(r.v.e[1]) > EPS) {
		for (int i = 0; i < BEZIER_SAMPLE; ++i) 
			if (max(ix[i], ix[i + 1]) >= dist) {
				real_t ky = iy[i + 1] - iy[i], by = iy[i];
				real_t kx = ix[i + 1] - ix[i], bx = ix[i];
				real_t ta = ky * ky * a - kx * kx;
				real_t tb = 2 * a * ky * by + ky * b - 2 * kx * bx;
				real_t tc = a * by * by + b * by + c - bx * bx;
				real_t delta = tb * tb - 4 * ta * tc;
				if (delta >= 0) {
					real_t t = (-tb + sqrt(delta)) / (2 * ta);
					if (t >= 0 && t <= 1) init.push_back(t / BEZIER_SAMPLE + it[i]);
					t = (-tb - sqrt(delta)) / (2 * ta);
					if (t >= 0 && t <= 1) init.push_back(t / BEZIER_SAMPLE + it[i]);
				}
			}
	}
	// for (auto t : init) cerr << t << endl;
}

hit_record rotation::hit(const ray &r) const {
	hit_record ret;
	vector <real_t> init;
	init.reserve(NEWTON_INIT);

	if (bd != NULL && !bd->hit(r)) return ret;
	real_t t = 0;
	real_t x, y, dx, dy;
	bool flag = 0;
	for (int step = 0; step < NEWTON_LIMIT; ++step) {
		f->eval(t, x, y);
		f->eval_deriv(t, dx, dy);
		real_t f = y - r.s.e[1];
		if (fabs(f) < NEWTON_EPS) { flag = 1; break; }
		real_t df = dy;
		real_t d = -f / df;
		d = max(d, -t);
		d = min(d, 1 - t);
		t += d;
	}
	if (flag) {
		if (fabs(r.v.e[1]) > EPS) {
			init.push_back(t);
		} else {
			real_t a = sqr(r.v.e[0]) + sqr(r.v.e[2]);
			real_t b = 2 * (r.v.e[0] * r.s.e[0] + r.v.e[2] * r.s.e[2]);
			real_t c = sqr(r.s.e[0]) + sqr(r.s.e[2]) - sqr(x);
			real_t delta = b * b - 4 * a * c;
			if (delta > 0) {
				real_t k = (-b + sqrt(delta)) / (2 * a);
				if (k >= EPS && (ret.hit == false || (ret.hit == true && ret.t > k))) {
					vect norm;
					norm = vect(ret.p.e[0] - axis.e[0], 0, ret.p.e[2] - axis.e[2]).norm();
					norm *= dir;
					if (-dot(norm, r.v) > 0) {
						ret.hit = true;
						ret.t = k;
						ret.p = r[ret.t];
						ret.mat = mat;
						ret.norm = norm;
					}
				}
				k = (-b - sqrt(delta)) / (2 * a);
				if (k >= EPS && (ret.hit == false || (ret.hit == true && ret.t > k))) {
					vect norm;
					norm = vect(ret.p.e[0] - axis.e[0], 0, ret.p.e[2] - axis.e[2]).norm();
					norm *= dir;
					if (-dot(norm, r.v) > EPS) {
						ret.hit = true;
						ret.t = k;
						ret.p = r[ret.t];
						ret.mat = mat;
						ret.norm = norm;
					}
				}
			}
			return ret;
		}
	}

	if (fabs(r.v.e[1]) > EPS) {
		real_t pa = sqr(r.v.e[0]) + sqr(r.v.e[2]);
		real_t pb = 2 * (r.v.e[0] * (r.s.e[0] - axis.e[0]) + r.v.e[2] * (r.s.e[2] - axis.e[2]));
		real_t pc = sqr(r.s.e[0] - axis.e[0]) + sqr(r.s.e[2] - axis.e[2]);
		real_t a = pa / sqr(r.v.e[1]);
		real_t b = pb / r.v.e[1] - 2 * pa * r.s.e[1] / sqr(r.v.e[1]);
		real_t c = pa * sqr(r.s.e[1] / r.v.e[1]) - pb * r.s.e[1] / r.v.e[1] + pc;

		real_t x, y;
		for (real_t t = -0.0; t <= 0.2; t += 0.001) {
			f->eval(t, x, y);
			real_t f = ((a * y) + b) * y + c - x * x;
			cerr << "{" << t << ' ' << f << "}" << endl;
		}


		for (int i = 0; i < NEWTON_INIT; ++i) init.push_back((0.5 + i) / NEWTON_INIT);
		get_init(r, init);

		for (auto t : init) {
			real_t x, y, dx, dy;
			bool flag = 0;
			for (int step = 0; step < NEWTON_LIMIT; ++step) {
				f->eval(t, x, y);
				f->eval_deriv(t, dx, dy);
				real_t f = ((a * y) + b) * y + c - x * x;
				if (fabs(f) < NEWTON_EPS) { flag = 1; break; }
				real_t df = 2 * a * y * dy + b * dy - 2 * x * dx;
				real_t d = -f / df;
				// d = max(d, -t);
				// d = min(d, 1.0 - t);
				//cerr << f << ' ' << x << ' ' << y << endl;
				t += d;
				while (t > 1 || t < 0) {
					if (t > 1) t = 1 - alpha * (t - 1);
					if (t < 0) t = -alpha * t;
				}
			}
			if (flag) {
				real_t k = (y - r.s.e[1]) / r.v.e[1];
				// cerr << k << ' ' << r[k] << endl;
				if (k >= EPS && (ret.hit == false || (ret.hit == true && ret.t > k))) {
					vect norm;
					ret.p = r[k];
					norm = vect(ret.p.e[0] - axis.e[0], 0, ret.p.e[2] - axis.e[2]).norm();
					norm.e[1] = -dx / dy;
					norm *= dir;
					norm = norm.norm();
					if (-dot(norm, r.v) > EPS) {
						ret.hit = true;
						ret.t = k;
						ret.mat = mat;
						ret.norm = norm;
					}
					// cerr << dx / dy << ' ' << norm << endl;
				}
			}
		}
	}
	return ret;
}

#endif
