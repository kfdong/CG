# g++ main.cpp `libpng-config --ldflags` -std=c++11 -O2 -Wall -pg
g++ -fopenmp main.cpp `libpng-config --ldflags` -std=c++11 -O2 -Wall 
