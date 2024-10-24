// Tobias Rubel 2024
// UMD Parallel Algorithms Group
// This program generates a collection of n d-dimensional vectors on the unit sphere.
// It then perturbs the vectors in a manner inspired by the Watts-Strogatz random graph model:
// for each vector, for each of k vec

#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "../ParlayANN_utils/euclidian_point.h"
#include "../ParlayANN_utils/point_range.h"

PointRange<float, Euclidian_Point<float>> generate_random_points(int n, int d, float delta){


}


int main(int argc, char* argv[]) {
  commandLine P(argc,argv,
  "[-out <outfile>] [-n <n>] [-d <d>] [-delta <d>] [-interpolation <ip>] [-k <k>]");

  char* oFile = P.getOptionValue("-out");
  char* ipc = P.getOptionValue("-interpolation");
  int d = P.getOptionIntValue("-d", 128);
  float delta = P.getOptionFloatValue("-delta", 0.1);
  int n = P.getOptionIntValue("-n", 1000);
  int k = P.getOptionIntValue("-k", 10);


  std::string ip = std::string(ipc);
  if(ip != "LERP" && ip != "SLERP"){
    std::cout << "Error: invalid interpolation method. LERP or SLERP" << std::endl;
    abort();
  }

  std::cout << "Generating " << n << " points with dimension " << d << " using " << ip << " interpolation and delta " << delta << std::endl;



  parlay::sequence<parlay::sequence<pid>> answers;
  std::string base = std::string(bFile);
  std::string query = std::string(qFile);


  if(ip == "LERP"){
  } else if(ip == "SLERP"){

  }
  

  return 0;
}

