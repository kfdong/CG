#ifndef CAMERA_H
#define CAMERA_H
#include "config.h"
#include "objects.hpp"
#include "utils.hpp"
#include <bits/stdc++.h>
#include <random>
#include <png++/png.hpp>
#include <omp.h>

/*
vect path_tracing(ray r, object_list *obj, int depth) {
	auto h = obj -> hit(r);
	if (!h.hit) {
		vect unit_direction = (r.v).norm();
		real_t t = 0.5 * (unit_direction.e[1] + 1.0);
		return (1.0 - t) * vect(1, 1, 1) + t * vect(0.5, 0.7, 1.0);
	} else {
		if ((h.mat->is_lambertian()) && (((lambertian*)h.mat)->albedo).e[0] == -1) {
			real_t dist = (h.p - vect(0, 3, 0)).abs();
			real_t il = 1.0 / (1.0 + 0.1 * dist);
			return il * vect(1, 1, 1);
		}
		ray out;
		vect attenuation;
		if (depth < DEPTH_LIMIT && h.mat->scatter(r, h, attenuation, out)) {
			return attenuation * path_tracing(out, obj, depth + 1);
		} else {
			return vect(0, 0, 0);
		}
	}
}
*/

class camera {
public:
	int nx, ny; /* number of pixels */

	vect lookfrom;
	vect lookat;
	vect up;	/* up direction (i.e., dy) of the photot */
	real_t fov; /* fov is left to right in radian, always less than PI */
	real_t init_radius;
	real_t energy;

	camera() {}
	camera(int _nx, int _ny, vect _lookfrom, vect _lookat, vect _up, real_t _fov, real_t _radius, 
			real_t _energy) {
		nx = _nx, ny = _ny, nx = _nx;
		lookfrom = _lookfrom;
		lookat = _lookat;
		up = _up;
		fov = _fov;
		init_radius = _radius;
		energy = _energy;
	}

	object_list *obj;
	vector <recording_point*> rec;
	// KDTree *tree;
	hash_table *tree;
	void render(png::image <png::rgb_pixel> &image, object_list *_obj, lights *lit);
	void photon_pass(const ray &r, const vect &att_in, int depth);
	void hitpoint_pass(const ray &r, const vect &att_in, int depth, int pixx, int pixy);
};

void camera::render(png::image <png::rgb_pixel> &image, object_list *_obj, lights *lit) {
	omp_set_num_threads(nr_thread);
	obj = _obj;
	real_t half_width = tan(fov / 2);
	real_t half_height = half_width * ny / nx;
	vect w = (lookat - lookfrom).norm();
	vect horizontal = cross(w, up).norm();
	vect vertical = cross(horizontal, w).norm();
	vect corner = lookfrom - half_width * horizontal - half_height * vertical;

	vect *col = new vect[ny * nx];

	for (int round = 0; round < nr_round; ++round) {
		std::cerr << "running round " << round << std::endl;

		rec.clear();
		vector <pair <int, int> > hitpoints;
		for (int y = ny - 1; y >= 0; --y)
			for (int x = 0; x < nx; ++x) hitpoints.push_back(std::make_pair(x, y));
		random_shuffle(hitpoints.begin(), hitpoints.end());

		int nr_batch = 1000;
		Progress_bar progress("hitpoint pass");
		progress.set_total(nr_batch);
		int n = hitpoints.size();
#pragma omp parallel for 
		for (int t = 0; t < n; ++t) {
			int x = hitpoints[t].first, y = hitpoints[t].second;
			ray r(lookfrom, corner - lookfrom + w + 
				2 * half_width * (x + rand_unit(gen)) / nx * horizontal + 
				2 * half_height * (y + rand_unit(gen)) / ny * vertical);
			/* ray r(lookfrom, corner - lookfrom + w + 
				2 * half_width * x / nx * horizontal + 
				2 * half_height * y / ny * vertical); */
			int pixx = x, pixy = ny - y - 1;
			// if (pixx == 662 && pixy == 408) cerr << r.s << ' ' << r.v << endl;
			// if (pixx == 662 && pixy == 407) cerr << r.s << ' ' << r.v << endl;
			hitpoint_pass(r, vect(1, 1, 1), 0, pixx, pixy);
			if (t % (n / nr_batch) == 0) 
#pragma omp critical 
				progress.update(1);
		}
		for (auto hp : rec) hp->rad2 = sqr(init_radius);
		progress.finish();
		// tree = new KDTree;
		tree = new hash_table;
		tree->build(rec, sqr(init_radius));

		progress = Progress_bar("photon pass");
		progress.set_total(nr_batch);
		for (int i = 0; i < nr_batch; ++i) {
#pragma omp parallel for 
			for (int j = 0; j < photon / nr_batch; ++j) {
				vect col_in;
				ray r = lit->emit(col_in);
				photon_pass(r, col_in, 0);
			}
			progress.update(1);
		}
		progress.finish();
		
		for (auto hp : rec) {
			int x = hp->pixx, y = hp->pixy;
			col[y * nx + x] += hp->flux / (PI * hp->rad2 * photon);
		}
		for (auto hp : rec) delete hp;
		delete tree;
	}
	for (int y = 0; y < ny; ++y)
		for (int x = 0; x < nx; ++x) {
			int id = y * nx + x;
			col[id] *= energy / nr_round;
			for (int t = 0; t < 3; ++t) col[id].e[t] = pow(1.0 - exp(-col[id].e[t]), 1.0 / 2.2) * 255 + 0.5;
			image[y][x] = png::rgb_pixel(int(col[id].e[0]), int(col[id].e[1]), int(col[id].e[2]));
			// std::cerr << col[id] << std::endl;
		}
}

void camera::hitpoint_pass(const ray &r, const vect &att_in, int depth, int pixx, int pixy) {
	if (att_in.abs() < EPS) return;
	// cerr << "pass " << depth << endl;
	auto h = obj -> hit(r);
	// cerr << "pass " << h.p << ' ' << att_in << ' ' << r.s << ' ' << r.v << ' ' << h.t << ' ' << depth << endl;
	//if (att_in.e[0] + att_in.e[1] + att_in.e[2] < EPS) return;
	if (!h.hit) return;
	ray out;
	vect attenuation;
	if (h.mat->type == LAMBERTIAN) {
		h.mat->scatter(r, h, attenuation, out);
		attenuation *= att_in;
		recording_point *now = new recording_point;
		now->x = h.p;
		now->attenuation = attenuation;
		// now->attenuation = attenuation * -dot(h.norm, r.v);
		now->norm = h.norm;
		now->flux = vect(0, 0, 0);
		now->n = 0;
		now->pixx = pixx; now->pixy = pixy;
#pragma omp critical 
		rec.push_back(now);
		return;
	} else if (h.mat->type == METAL) {
		h.mat->scatter(r, h, attenuation, out);
		attenuation *= att_in;
		hitpoint_pass(out, attenuation, depth + 1, pixx, pixy);
	} else if (h.mat->type == GLASS) {
		ray refrac, reflac;
		real_t refrac_p;
		((glass*)h.mat)->scatter(r, h, attenuation, refrac, reflac, refrac_p);
		attenuation *= att_in;
		if (refrac_p != 0)
			hitpoint_pass(refrac, attenuation * refrac_p, depth + 1, pixx, pixy);
		hitpoint_pass(reflac, attenuation * (1.0 - refrac_p), depth + 1, pixx, pixy);
	}
}

void camera::photon_pass(const ray &r, const vect &att_in, int depth) {
	if (depth >= DEPTH_LIMIT) return;
	auto h = obj -> hit(r);
	if (!h.hit) return;
	ray out;
	vect attenuation;
	h.mat->scatter(r, h, attenuation, out);
	if (h.mat->type == LAMBERTIAN) {
		auto rec = tree->query(h.p);
		if (rec != NULL) {
#pragma omp critical 
			for (auto hp : *rec) {
				if (dot(hp->norm, h.norm) > EPS && (hp->x - h.p).abs2() < hp->rad2) {
					real_t factor = (hp->n + alpha) / (hp->n + 1.0);
					hp->rad2 *= factor;
					hp->flux = (hp->flux + hp->attenuation * att_in) * factor;
					hp->n++;
					// std::cerr << att_in << std::endl;
				}
			}
		}
	}
	attenuation *= att_in;
	real_t p = std::max(attenuation.e[0], std::max(attenuation.e[1], attenuation.e[2]));
	if (rand_unit(gen) < p) photon_pass(out, attenuation / p, depth + 1);
	// photon_pass(out, attenuation, depth + 1);
}

#endif
