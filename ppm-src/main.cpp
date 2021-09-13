#include <bits/stdc++.h>
#include <random>
#include <png++/png.hpp>
#include "config.h"
#include "utils.hpp"
#include "objects.hpp"
#include "camera.hpp"

object_list obj;

Bezier *read_bezier(FILE *model, real_t scale, real_t y_shift) {
	int n;
	real_t x[20], y[20];
	fscanf(model, "%d", &n);
	for (int i = 0; i <= n; ++i) {
		fscanf(model, "%lf%lf", x + i, y + i);
		x[i] *= scale, y[i] *= scale;
		y[i] += y_shift;
	}
	Bezier *B = new Bezier(n, x, y);
	return B;
}
bounding_box *read_bounding_box(FILE *model, real_t scale, real_t y_shift) {
	int nr_v, nr_s;
	if (fscanf(model, "%d%d", &nr_v, &nr_s) == 2) {
		bounding_box *ret = new bounding_box();
		object_list *bd = new object_list;
		vector <vect> vs;
		for (int i = 0; i < nr_v; ++i) {
			vect v;
			v.read(model);
			vs.push_back(v * scale + vect(0, y_shift, 0));
		}
		for (int i = 0; i < nr_s; ++i) {
			int nt = 0;
			fscanf(model, "%d", &nt);
			vector <vect> poly;
			for (int j = 0; j < nt; ++j) {
				int p;
				fscanf(model, "%d", &p);
				poly.push_back(vs[p]);
			}
			bd->add(*new polygon(poly, NULL));
			std::reverse(poly.begin(), poly.end());
			bd->add(*new polygon(poly, NULL));
		}
		ret->obj = (object*)bd;
		return ret;
	} else return NULL;
}

void build_test_scene(png::image<png::rgb_pixel> &image, int nx, int ny) {
	bool raindrop = false;
	bool winecup = true;

	/*
	real_t half_width = 4;
	real_t half_height = 3;
	real_t half_depth = 5;*/
	real_t half_width = 40;
	real_t half_height = 3;
	real_t half_depth = 50;
	real_t back_depth = 20;
	/* box */
	vect p[8] = {
		vect(half_width, -half_height, -back_depth),
		vect(half_width, -half_height, half_depth),
		vect(-half_width, -half_height, half_depth),
		vect(-half_width, -half_height, -back_depth),
		vect(half_width, half_height, -back_depth),
		vect(half_width, half_height, half_depth),
		vect(-half_width, half_height, half_depth),
		vect(-half_width, half_height, -back_depth),
	};

	polygon *poly;
	// poly = new polygon({p[0], p[1], p[5], p[4]}, new lambertian(vect(0.8, 0.3, 0.3))); obj.add(*poly); // left
	// poly = new polygon({p[3], p[7], p[6], p[2]}, new lambertian(vect(0.3, 0.3, 0.8))); obj.add(*poly); // right
	// poly = new polygon({p[4], p[5], p[6], p[7]}, new lambertian(vect(0.95, 0.95, 0.95))); obj.add(*poly); // up
	// poly = new polygon({p[0], p[3], p[2], p[1]}, new lambertian_mesh(vect(0.95, 0.95, 0.95), vect(0.1, 0.1, 0.1))); obj.add(*poly); // down
	poly = new polygon({p[0], p[3], p[2], p[1]}, new lambertian_texture("texture/wood_texture.jpg", vect(100, 0, 0), vect(0, 0, 100))); obj.add(*poly); // down
	// poly = new polygon({p[0], p[3], p[2], p[1]}, new lambertian(vect(0.95, 0.95, 0.95))); obj.add(*poly); // down
	// poly = new polygon({p[1], p[2], p[6], p[5]}, new lambertian(vect(0.0, 0.95, 0.95))); obj.add(*poly); // front
	/*
	s = new sphere({r + half_width, 0, 0}, r, new lambertian(vect(0.8, 0.3, 0.3))); obj.add(*s); // left
	s = new sphere({-r - half_width, 0, 0}, r, new lambertian(vect(0.3, 0.3, 0.8))); obj.add(*s); // right
	s = new sphere({0, r + half_height, 0}, r, new lambertian(vect(0.95, 0.95, 0.95))); obj.add(*s); // up
	s = new sphere({0,-r - half_height, 0}, r, new lambertian(vect(0.95, 0.95, 0.95))); obj.add(*s); // down
	s = new sphere({0, 0, r + half_depth}, r, new lambertian(vect(0.95, 0.95, 0.95))); obj.add(*s); // front
	s = new sphere({0, 0,-r - back_depth}, r, new lambertian(vect(0, 0, 0))); obj.add(*s); // back
	*/
	/*
	triangle *t;
	t = new triangle({sqrt(3), -1, half_depth - EPS}, {-sqrt(3), -1, half_depth - EPS}, {0, 2, half_depth - EPS}, new lambertian(vect(0.95, 0.6, 0.95))); obj.add(*t); // back

	sphere *s;
	s = new sphere({0.5, -2.4, 0}, 0.6, new lambertian(vect(0.4, 0.5, 0.4))); obj.add(*s);
	s = new sphere({2, -1.5, -2}, 1.5, new metal(vect(0.3, 0.7, 0.2), 0.0)); obj.add(*s);
	// s = new sphere({-1, 0, -1}, 0.5, new metal(vect(0.8, 0.8, 0.8), 0.3)); obj.add(*s);
	s = new sphere({-2, -1.5, -1}, 1.5, new glass(1.0, 1.0 / 1.5, 1.0)); obj.add(*s);
	s = new sphere({-2, -1.5, -1}, -1.5, new glass(1.0 / 1.5, 1.0, 1.0)); obj.add(*s);
	*/
	/* raindrop */
	if (raindrop) {
		FILE *model = fopen("raindrop.data", "r");

		real_t y_shift = 0.5;
		real_t scale = 1.0;
		Bezier *B = read_bezier(model, scale, y_shift);
		rotation *r;
		r = new rotation(1.0, {0, 0, 0}, B, new glass(1.0, 1.0 / 1.5, 1.0));
		// r = new rotation(1.0, {-1, 0, -2}, B, new metal(vect(0.8, 0.8, 0.8)));
		r->bd = read_bounding_box(model, scale, y_shift);
		obj.add(*r);

		rotation *rev_r = new rotation;
		*rev_r = *r;
		rev_r->dir = -1.0;
		rev_r->mat = new glass(1.0 / 1.5, 1.0, 1.0);
		obj.add(*rev_r);

		hemisphere *s;
		s = new hemisphere({0, 0 + 0.5 * scale + y_shift, 0}, 0.5 * scale, {0, -1, 0}, new glass(1.0, 1.0 / 1.5, 1.0)); obj.add(*s);
		s = new hemisphere({0, 0 + 0.5 * scale + y_shift, 0},-0.5 * scale, {0, -1, 0}, new glass(1.0 / 1.5, 1.0, 1.0)); obj.add(*s);
	}
	if (winecup) {
		real_t y_shift = -3.0;
		real_t scale = 3.0 / 35;
		{
			FILE *model = fopen("winecup_1.data", "r");

			Bezier *B = read_bezier(model, scale, y_shift);
			rotation *r;
			r = new rotation(1.0, {0, 0, 0}, B, new glass(1.0, 1.0 / 1.5, 0.9));
			// r = new rotation(1.0, {0, 0, 0}, B, new lambertian(vect(0.8, 0.8, 0.8)));
			r->bd = read_bounding_box(model, scale, y_shift);
			obj.add(*r);

			rotation *rev_r = new rotation;
			*rev_r = *r;
			rev_r->dir = -1.0;
			rev_r->mat = new glass(1.0 / 1.5, 1.0, 0.9);
			// rev_r->mat = new lambertian(vect(0, 0, 0));
			obj.add(*rev_r);
		}

		{
			FILE *model = fopen("winecup_2.data", "r");

			Bezier *B = read_bezier(model, scale, y_shift);
			rotation *r;
			r = new rotation(-1.0, {0, 0, 0}, B, new glass(1.0, 1.0 / 1.5, 0.9));
			// r = new rotation(-1.0, {0, 0, 0}, B, new lambertian(vect(0.8, 0.8, 0.8)));
			r->bd = read_bounding_box(model, scale, y_shift);
			obj.add(*r);

			rotation *rev_r = new rotation;
			*rev_r = *r;
			rev_r->dir = 1.0;
			rev_r->mat = new glass(1.0 / 1.5, 1.0, 0.9);
			// rev_r->mat = new lambertian(vect(0, 0, 0));
			obj.add(*rev_r);
		}
	}

	camera *cam;
	// vect lookfrom(0.0, 6, -back_depth + 1.0);
	// vect lookat(0.0, 0.5, 0.0);
	vect lookfrom(-2, 10, -back_depth - 2);
	vect lookat(-2, 0.5, -4.5);
	vect up(0.0, 1.0, 0.0);
	real_t fov = 0.15 * PI;
	real_t init_radius = 0.01;
	// lights *lit = new point_light(vect(3, 3, 4), vect(-1, -1, -1), vect(1, 1, 1), 0.5);
	// lights *lit = new point_light(vect(6, 6, 8), vect(-1, -1, -1), vect(1, 1, 1), 0.7);
	lights *lit = new gaussian_light(vect(3, 3, 4), vect(-1, -2, -1), vect(1, 0.9, 0.9), 0.7, 0.05);
	// lights *lit = new point_light(vect(4, 3, 5), vect(-1, -0.5, -1), vect(1, 1, 1));
	// lights *lit = new hemisphere_light(vect(0, 3, 0), vect(0, -1, 0), vect(1, 1, 1));
	cam = new camera(nx, ny, lookfrom, lookat, up, fov, init_radius, 20);
	cam->render(image, &obj, lit);

	/*
	// ray r({0, 1, -19}, {0, -0.0586749, 0.998277});
	ray r({0, 1, -19}, {-0.0335479, -0.0564816, 0.99784});
	// ray r({0, 1, -19}, {-0.0335475, -0.0559415, 0.99787});

	cam->obj = &obj;
	cam->hitpoint_pass(r, vect(1, 1, 1), 0);
	cerr << cam->rec.size() << endl;
	for (auto hp : cam->rec) {
		cerr << hp->x << ' ' << hp->norm << ' ' << hp->attenuation << endl;
	}
	auto h = obj.hit(r);
	std::cerr << h.hit << std::endl;
	std::cerr << h.p << std::endl;
	std::cerr << h.norm << std::endl;
	cerr << -dot(h.norm, r.v) << endl;
	cerr << h.t << endl;
	return;
*/
}

int main(int argc, char **argv) {
	int nx = 1200, ny = 900;
	png::image< png::rgb_pixel > image(nx, ny);
	// build_test_scene();
	build_test_scene(image, nx, ny);
	image.write("cornel_box.png");

/*
	FILE *fout = fopen("ball.ppm", "w");
	fprintf(fout, "P3\n%d %d\n255\n", nx, ny);
	for (int j = 0; j < ny; ++j)
		for (int i = 0; i < nx; ++i)
			fprintf(fout, "%d %d %d\n", image[j][i].red, image[j][i].green, image[j][i].blue);
*/
	return 0;
}
