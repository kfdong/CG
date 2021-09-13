#include <stdint.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <bits/stdc++.h>
#define real_t double
using namespace std;
const real_t EPS = 1e-8;
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

const real_t alpha = 0.5;
const int SEPERATE = 3;
const int MaxN = 4000, CHANNEL = 3, INF = 1 << 30;
int image[MaxN][MaxN][CHANNEL];
bool flag[MaxN][MaxN];
real_t ene[MaxN][MaxN];
real_t energy[MaxN][MaxN];
real_t dp[MaxN][MaxN];
int width, height, bpp;
pair <int, int> id[MaxN][MaxN];

void compute_energy(real_t alpha) {
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j) {
			energy[i][j] = 0;
#ifdef OPERATOR1
			for (int k = 0; k < 3; ++k) 
				energy[i][j] += (alpha * iabs(image[i][j + 1][k] - image[i][j][k]) + (1 - alpha) * iabs(image[i + 1][j][k] - image[i][j][k]));
			energy[i][j] = sqrt(energy[i][j]);
#else
			for (int k = 0; k < 3; ++k) 
				energy[i][j] += (alpha * sqr(image[i][j + 1][k] - image[i][j][k]) + (1 - alpha) * sqr(image[i + 1][j][k] - image[i][j][k]));
			energy[i][j] = sqrt(energy[i][j]);
			energy[i][j] = sqrt(energy[i][j]);
#endif
			if (flag[i][j]) energy[i][j] = INF;
		}
}

void find_column_seam(int width, int height, vector<vector<int> > &paths) {
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
			if (now - 1 >= 0 && cmp(dp[i][now - 1], dp[i + 1][now] - energy[i + 1][now]) == 0)
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
		}
	}
}

void find_row_seam(int width, int height, vector<vector<int> > &paths) {
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
			if (now - 1 >= 0 && dp[now - 1][j] == dp[now][j + 1] - energy[now][j + 1])
				now = now - 1;
			else if (now + 1 < height && dp[now + 1][j] == dp[now][j + 1] - energy[now][j + 1])
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
		}
	}
}

void seam_carving(int h_rem, int w_rem) {
	vector<vector<int> > paths;

	/* inserting column */
	cout << "inserting colulms" << endl;
	progress.set_total(w_rem);

	int STEP_SIZE = w_rem * 0.2;
	while (w_rem > 0) {
		paths.clear();
		memset(flag, 0, sizeof(flag));
		while ((int)paths.size() < STEP_SIZE) {
			find_column_seam(width, height, paths);
		}
		int add = 0;
		memset(flag, 0, sizeof(flag));
		for (auto &p : paths) {
			for (int i = 0; i < height; ++i)
				flag[i][p[i]] = 1;
			++add ;
			if (add == w_rem) break;
		}
		
		static int from[MaxN][MaxN];
		memset(from, -1, sizeof(from));
		for (int i = 0; i < height; ++i) {
			int pos = 0;
			for (int j = 0; j < width; ++j) {
				if (flag[i][j]) ++pos;
				from[i][pos] = j;
				++pos;
			}
		}
		for (int i = 0; i < height; ++i)
			for (int j = width + add - 1; j >= 0; --j)
				if (from[i][j] != -1) {
					memcpy(image[i][j], image[i][from[i][j]], sizeof(image[i][j]));
					id[i][j] = id[i][from[i][j]];
				} else id[i][j] = make_pair(-1, -1);
		for (int i = 0; i < height; ++i)
			for (int j = width + add - 1; j >= 0; --j)
				if (from[i][j] == -1) {
					for (int k = 0; k < 3; ++k) {
						if (j - 1 < 0) image[i][j][k] = image[i][j + 1][k];
						else if (j + 1 >= width + add) image[i][j][k] = image[i][j - 1][k];
						else image[i][j][k] = (image[i][j - 1][k] + image[i][j + 1][k]) / 2;
					}
				}

		w_rem -= add;
		width += add;
		progress.update(add);

		/* padding */
		for (int i = 0; i < height; ++i)
			for (int k = 0; k < 3; ++k)
				image[i][width][k] = 2 * image[i][width - 1][k] - image[i][width - 2][k];
		for (int j = 0; j < width; ++j)
			for (int k = 0; k < 3; ++k)
				image[height][j][k] = 2 * image[height - 1][j][k] - image[height - 2][j][k];
	}
	progress.finish();

	/* inserting row */
	cout << "inserting rows" << endl;
	STEP_SIZE = h_rem * 0.2;
	progress.set_total(h_rem);
	while (h_rem > 0) {
		paths.clear();
		memset(flag, 0, sizeof(flag));
		while ((int)paths.size() < STEP_SIZE) {
			find_row_seam(width, height, paths);
		}
		int add = 0;
		memset(flag, 0, sizeof(flag));
		for (auto &p : paths) {
			for (int j = 0; j < width; ++j)
				flag[p[j]][j] = 1;
			++add;
			if (add == h_rem) break;
		}

		static int from[MaxN][MaxN];
		memset(from, -1, sizeof(from));
		for (int j = 0; j < width; ++j) {
			int pos = 0;
			for (int i = 0; i < height; ++i) {
				if (flag[i][j]) ++pos;
				from[pos][j] = i;
				++pos;
			}
		}
		for (int j = 0; j < width; ++j) 
			for (int i = height + add - 1; i >= 0; --i)
				if (from[i][j] != -1) {
					memcpy(image[i][j], image[from[i][j]][j], sizeof(image[i][j]));
					id[i][j] = id[from[i][j]][j];
				} else id[i][j] = make_pair(-1, -1);
		for (int j = 0; j < width; ++j) 
			for (int i = height + add - 1; i >= 0; --i)
				if (from[i][j] == -1) {
					for (int k = 0; k < 3; ++k) {
						if (i - 1 < 0) image[i][j][k] = image[i + 1][j][k];
						else if (i + 1 >= height + add) image[i][j][k] = image[i - 1][j][k];
						else image[i][j][k] = (image[i - 1][j][k] + image[i + 1][j][k]) / 2;
					}
				}
		h_rem -= add;
		height += add;
		progress.update(add);

		/* padding */
		for (int i = 0; i < height; ++i)
			for (int k = 0; k < 3; ++k)
				image[i][width][k] = 2 * image[i][width - 1][k] - image[i][width - 2][k];
		for (int j = 0; j < width; ++j)
			for (int k = 0; k < 3; ++k)
				image[height][j][k] = 2 * image[height - 1][j][k] - image[height - 2][j][k];
	}
	progress.finish();
}

int main(int argc, char **argv) {
	if (argc <= 1) {
		printf("usage %s <filename>\n", argv[0]);
		return 0;
	}

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

	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			id[i][j] = make_pair(i, j);

	stbi_image_free(rgb_image);

	seam_carving(int(0.0 * height), int(0.5 * width));

	FILE *fout = fopen("output.ppm", "w");
	fprintf(fout, "P3\n%d %d\n255\n", width, height);
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			fprintf(fout, "%d %d %d\n", image[i][j][0], image[i][j][1], image[i][j][2]);
#ifdef SHOWSEAM
	FILE *fshow = fopen("output_showseam.ppm", "w");
	fprintf(fshow, "P3\n%d %d\n255\n", width, height);
	memset(flag, 0, sizeof(flag));
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j) {
			if (id[i][j].first == -1) image[i][j][0] = 255;
			fprintf(fshow, "%d %d %d\n", image[i][j][0], image[i][j][1], image[i][j][2]);
		}
#endif
    return 0;
}

