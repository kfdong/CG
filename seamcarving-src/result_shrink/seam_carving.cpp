#include <stdint.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <bits/stdc++.h>
#define real_t double
using namespace std;
const real_t EPS = 1e-4;
int sqr(int x) { return x * x; }
int iabs(int x) { return x < 0 ? -x : x; }
int sign(real_t x) { return x < -EPS ? -1 : (x > EPS ? 1 : 0); }
int cmp(real_t x, real_t y) { return sign(x - y); }

class Progress_bar {
	int total, current;
	int barWidth = 70;
public:

	Progress_bar() {}
	void set_total(int _total) {
		current = 0;
		total = _total;
	}
	void update(int part) {
		current += part;
		real_t progress = min((real_t)current / total, (real_t)1.0);

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress * 100.0) << " %\r";
		std::cout.flush();
	}
	void finish() {
		std::cout << std::endl;
	}
} progress;

real_t alpha = 0.5;
const int MaxN = 3000, CHANNEL = 3, INF = 1 << 30;
int image[MaxN][MaxN][CHANNEL];
bool flag[MaxN][MaxN];
real_t ene[MaxN][MaxN];
real_t energy[MaxN][MaxN];
real_t dp[MaxN][MaxN];
int width, height, bpp;
#ifdef SHOWSEAM
pair <int, int> id[MaxN][MaxN];
#endif

void compute_energy(real_t alpha) {
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j) {
			energy[i][j] = 0;
			int il = i - 1 < 0 ? i : i - 1;
			int ih = i + 1 == height ? i : i + 1;
			int jl = j - 1 < 0 ? j : j - 1;
			int jh = j + 1 == width ? j : j + 1;
			for (int k = 0; k < 3; ++k) {
				energy[i][j] += 2 * (alpha * iabs(image[i][jh][k] - image[i][jl][k]) + (1 - alpha) * iabs(image[ih][j][k] - image[il][j][k]));
				energy[i][j] += (alpha * iabs(image[il][jh][k] - image[il][jl][k]) + (1 - alpha) * iabs(image[ih][jl][k] - image[il][jl][k]));
				energy[i][j] += (alpha * iabs(image[ih][jh][k] - image[ih][jl][k]) + (1 - alpha) * iabs(image[ih][jh][k] - image[il][jh][k]));
			}
		}
}

void find_column_seam(int width, int height, vector<vector<int> > &paths) {
	paths.clear();
	memset(flag, 0, sizeof(flag));
	compute_energy(alpha);
	/* dp */
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			dp[i][j] = energy[i][j];

	for (int i = 1; i < height; ++i) {
		dp[i][0] += min(dp[i - 1][0], dp[i - 1][1]);
		for (int j = 1; j < width - 1; ++j) {
			dp[i][j] += min(dp[i - 1][j - 1], min(dp[i - 1][j], dp[i - 1][j + 1]));
		}
		dp[i][width - 1] += min(dp[i - 1][width - 2], dp[i - 1][width - 1]);
	}
	/* find path */
	static pair <int, int> pt[MaxN];
	for (int j = 0; j < width; ++j) 
		pt[j] = make_pair(dp[height - 1][j], j);
	sort(pt, pt + width);

	real_t min_cost = INF;
	for (int j = 0; j < width; ++j) {
		vector <int> p;
		int now = pt[j].second;
		p.push_back(now);
		for (int i = height - 2; i >= 0; --i) {
			if (cmp(dp[i][now], dp[i + 1][now] - energy[i + 1][now]) == 0)
				now = now + 0;
			else if (now - 1 >= 0 && cmp(dp[i][now - 1], dp[i + 1][now] - energy[i + 1][now]) == 0)
				now = now - 1;
			else if (now + 1 < width && cmp(dp[i][now + 1], dp[i + 1][now] - energy[i + 1][now]) == 0)
				now = now + 1;
			if (flag[i][now]) break;
			p.push_back(now);
		}
		if ((int)p.size() == height) {
			real_t current_cost = dp[height - 1][p[0]];
			min_cost = min(min_cost, current_cost);
			reverse(p.begin(), p.end());
			for (int i = 0; i < height; ++i)
				flag[i][p[i]] = 1;
			paths.push_back(p);
			if (current_cost > min_cost * 2) break;
		}
	}
}

void find_row_seam(int width, int height, vector<vector<int> > &paths) {
	paths.clear();
	memset(flag, 0, sizeof(flag));
	compute_energy(1 - alpha);
	/* dp */
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			dp[i][j] = energy[i][j];

	for (int j = 1; j < width; ++j) {
		dp[0][j] += min(dp[0][j - 1], dp[1][j - 1]);
		for (int i = 1; i < height - 1; ++i) {
			dp[i][j] += min(dp[i - 1][j - 1], min(dp[i][j - 1], dp[i + 1][j - 1]));
		}
		dp[height - 1][j] += min(dp[height - 1][j - 1], dp[height - 2][j - 1]);
	}
	/* find path */
	static pair <int, int> pt[MaxN];
	for (int i = 0; i < height; ++i)
		pt[i] = make_pair(dp[i][width - 1], i);
	sort(pt, pt + height);

	real_t min_cost = INF;
	for (int i = 0; i < height; ++i) {
		vector <int> p;
		int now = pt[i].second;
		p.push_back(now);
		for (int j = width - 2; j >= 0; --j) {
			if (now - 1 >= 0 && cmp(dp[now - 1][j], dp[now][j + 1] - energy[now][j + 1]) == 0)
				now = now - 1;
			else if (now + 1 < height && cmp(dp[now + 1][j], dp[now][j + 1] - energy[now][j + 1]) == 0)
				now = now + 1;
			if (flag[now][j]) break;
			p.push_back(now);
		}
		if ((int)p.size() == width) {
			real_t current_cost = dp[p[0]][width - 1];
			min_cost = min(min_cost, current_cost);
			reverse(p.begin(), p.end());
			for (int j = 0; j < width; ++j)
				flag[p[j]][j] = 1;
			paths.push_back(p);
			if (current_cost > min_cost * 2) break;
		}
	}
}

void seam_carving(int h_rem, int w_rem) {
	vector<vector<int> > paths;

	/* remove column */
	cout << "removing colulms" << endl;
	progress.set_total(w_rem);

	while (w_rem > 0) {
		/* padding */
		for (int i = 0; i < height; ++i)
			for (int k = 0; k < 3; ++k)
				image[i][width][k] = 2 * image[i][width - 1][k] - image[i][width - 2][k];
		for (int j = 0; j < width; ++j)
			for (int k = 0; k < 3; ++k)
				image[height][j][k] = 2 * image[height - 1][j][k] - image[height - 2][j][k];

		find_column_seam(width, height, paths);
		int del = 0;
		if (w_rem < (int)paths.size()) {
			/* remove fewer paths than found */
			memset(flag, 0, sizeof(flag));
			for (auto &p : paths) {
				for (int i = 0; i < height; ++i)
					flag[i][p[i]] = 1;
				++del;
				if (del == w_rem) break;
			}
		} else del = paths.size();

		for (int i = 0; i < height; ++i) {
			int pos = 0;
			for (int j = 0; j < width; ++j)
				if (!flag[i][j]) {
					memcpy(image[i][pos], image[i][j], sizeof(image[i][j]));
#ifdef SHOWSEAM
					id[i][pos] = id[i][j];
#endif
					++pos;
				}
		}
		w_rem -= del;
		width -= del;
		progress.update(del);
	}
	progress.finish();

	/* remove row */
	cout << "removing rows" << endl;
	progress.set_total(h_rem);
	while (h_rem > 0) {
		/* padding */
		for (int i = 0; i < height; ++i)
			for (int k = 0; k < 3; ++k)
				image[i][width][k] = 2 * image[i][width - 1][k] - image[i][width - 2][k];
		for (int j = 0; j < width; ++j)
			for (int k = 0; k < 3; ++k)
				image[height][j][k] = 2 * image[height - 1][j][k] - image[height - 2][j][k];

		find_row_seam(width, height, paths);
		int del = 0;
		if (h_rem < (int)paths.size()) {
			/* remove fewer paths than found */
			memset(flag, 0, sizeof(flag));
			for (auto &p : paths) {
				for (int j = 0; j < width; ++j)
					flag[p[j]][j] = 1;
				++del;
				if (del == h_rem) break;
			}
		} else del = paths.size();

		for (int j = 0; j < width; ++j) {
			int pos = 0;
			for (int i = 0; i < height; ++i) 
				if (!flag[i][j]) {
					memcpy(image[pos][j], image[i][j], sizeof(image[i][j]));
#ifdef SHOWSEAM
					id[pos][j] = id[i][j];
#endif
					++pos;
				}
		}
		h_rem -= del;
		height -= del;
		progress.update(del);
	}
	progress.finish();
}

int main(int argc, char **argv) {
	if (argc <= 1) {
		printf("usage %s <filename> alpha\n", argv[0]);
		return 0;
	}
	alpha = stod(argv[2]);

    uint8_t* rgb_image = stbi_load(argv[1], &width, &height, &bpp, 3);
	printf("image size : (%d %d)\n", width, height);
	if (width < 0 || width >= MaxN || height < 0 || height >= MaxN) {
		printf("size error\n");
		return 0;
	}
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			for (int k = 0; k < 3; ++k)
				image[i][j][k] = rgb_image[(i * width + j) * CHANNEL + k];
	/* padding */
	for (int i = 0; i < height; ++i)
		for (int k = 0; k < 3; ++k)
			image[i][width][k] = 2 * image[i][width - 1][k] - image[i][width - 2][k];
	for (int j = 0; j < width; ++j)
		for (int k = 0; k < 3; ++k)
			image[height][j][k] = 2 * image[height - 1][j][k] - image[height - 2][j][k];

#ifdef SHOWSEAM
	int w = width, h = height;
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			id[i][j] = make_pair(i, j);
	static int img[MaxN][MaxN][3];
	memcpy(img, image, sizeof(img));
#endif
	stbi_image_free(rgb_image);

	seam_carving(int(0.2 * height), int(0.2 * width));

	FILE *fout = fopen("output.ppm", "w");
	fprintf(fout, "P3\n%d %d\n255\n", width, height);
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			fprintf(fout, "%d %d %d\n", image[i][j][0], image[i][j][1], image[i][j][2]);
#ifdef SHOWSEAM
	FILE *fshow = fopen("output_showseam.ppm", "w");
	fprintf(fshow, "P3\n%d %d\n255\n", w, h);
	memset(flag, 0, sizeof(flag));
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j) {
			int x = id[i][j].first, y = id[i][j].second;
			flag[x][y] = 1;
		}
	int count = 0;
	for (int i = 0; i < h; ++i)
		for (int j = 0; j < w; ++j) {
			if (!flag[i][j]) img[i][j][0] = 255;
			else ++count;
			fprintf(fshow, "%d %d %d\n", img[i][j][0], img[i][j][1], img[i][j][2]);
		}
	cout << count - width * height << endl;
#endif
    return 0;
}

