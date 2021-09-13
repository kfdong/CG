#ifndef CONFIG_H
#define CONFIG_H

#include <bits/stdc++.h>
#include <random>
typedef std::mt19937 random_generator;
typedef double real_t;

using std::cerr;
using std::cout;
using std::endl;

const int nr_thread = 2;
const real_t EPS = 1e-6;
const real_t NEWTON_EPS = 1e-9;
const real_t INF = 1e15;
const real_t PI = acosl(-1);
const real_t alpha = 0.7;
const int nr_round = 4;
const int photon = 5000000;
const int DEPTH_LIMIT = 20;
const int NEWTON_INIT = 3;
const int BEZIER_SAMPLE = 20;
const int NEWTON_LIMIT = 20;
thread_local std::random_device rd;
thread_local random_generator gen(rd());
thread_local std::uniform_real_distribution<> rand_unit(0.0, 1.0);

#endif
