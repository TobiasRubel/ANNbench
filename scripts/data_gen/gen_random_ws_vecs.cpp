// Tobias Rubel 2024
// UMD Parallel Algorithms Group
// This program generates a collection of n d-dimensional vectors on the unit sphere.
// It then perturbs the vectors in a manner inspired by the Watts-Strogatz random graph model:
// for each vector, for each of k vec

// Include necessary headers
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/internal/get_time.h"
#include "../bench/parse_command_line.h"
#include "types.h"
#include "../ParlayANN_utils/Euclidian_Point.h"
#include "../ParlayANN_utils/PointRange.h"

// Heat function
float compute_heat(const float* x, const parlay::sequence<Euclidian_Point<float>>& points,
                   int d, float sigma, float hi) {
  float total_heat = 0.0f;
  float sigma_sq = sigma * sigma;
  for (const auto& point : points) {
    const float* p_vals = point.values;
    float distance_sq = 0.0f;
    for (int j = 0; j < d; ++j) {
      float diff = x[j] - p_vals[j];
      distance_sq += diff * diff;
    }
    float heat_contribution = hi * std::exp(-distance_sq / (2.0f * sigma_sq));
    total_heat += heat_contribution;
  }
  return total_heat;
}

// Metropolis-Hastings sampling function
float* metropolis_hastings_sampling(const parlay::sequence<Euclidian_Point<float>>& points,
                                    int d, float sigma, float hi,
                                    int burn_in, int thinning, std::mt19937& rng) {
  // Initialize starting point within the unit cube
  float* x_current = new float[d];
  std::uniform_real_distribution<float> uniform_dist(-1.0f, 1.0f);
  for (int j = 0; j < d; ++j) {
    x_current[j] = uniform_dist(rng);
  }
  
  float heat_current = compute_heat(x_current, points, d, sigma, hi);
  
  // Proposal standard deviation
  float proposal_std = sigma * 0.5f;
  std::normal_distribution<float> normal_dist(0.0f, proposal_std);
  
  int total_iterations = burn_in + thinning;
  float* x_proposal = new float[d];
  
  for (int iteration = 0; iteration < total_iterations; ++iteration) {
    // Propose a new point
    for (int j = 0; j < d; ++j) {
      x_proposal[j] = x_current[j] + normal_dist(rng);
      // Reflect back into unit cube if outside
      if (x_proposal[j] > 1.0f) x_proposal[j] = 1.0f - (x_proposal[j] - 1.0f);
      if (x_proposal[j] < -1.0f) x_proposal[j] = -1.0f - (x_proposal[j] + 1.0f);
      // Ensure within [-1, 1]
      x_proposal[j] = std::max(std::min(x_proposal[j], 1.0f), -1.0f);
    }
    
    float heat_proposal = compute_heat(x_proposal, points, d, sigma, hi);
    
    // Compute acceptance probability
    float acceptance_ratio = heat_proposal / (heat_current + 1e-8f); // Avoid division by zero
    float acceptance_probability = std::min(1.0f, acceptance_ratio);
    
    // Accept or reject the proposal
    std::uniform_real_distribution<float> accept_dist(0.0f, 1.0f);
    if (accept_dist(rng) < acceptance_probability) {
      // Accept proposal
      std::swap(x_current, x_proposal);
      heat_current = heat_proposal;
    }
    
    // No need to store samples for burn-in and thinning since we only need one sample
  }
  
  delete[] x_proposal;
  
  return x_current; // Caller is responsible for deleting x_current
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
    "[-out <outfile>] [-n <n>] [-d <d>] [-m <m>] [-sigma <sigma>] [-hi <hi>] [-burn_in <burn_in>] [-thinning <thinning>]");
  
  char* oFile = P.getOptionValue("-out");
  int d = P.getOptionIntValue("-d", 128);
  size_t n = P.getOptionLongValue("-n", 1000);
  size_t m = P.getOptionLongValue("-m", 10);
  
  float sigma = P.getOptionFloatValue("-sigma", 0.1f);
  float hi = P.getOptionFloatValue("-hi", 1.0f);
  int burn_in = P.getOptionIntValue("-burn_in", 20);
  int thinning = P.getOptionIntValue("-thinning", 5);
  
  // Step 1: Initialize m initial points
  using T = float; // Data type for point coordinates
  using Point = Euclidian_Point<T>;
  using PointRangeType = PointRange<T, Point>;
  
  // Prepare data structures
  parlay::sequence<Point> points;
  points.reserve(n); // Reserve space for n points
  
  // Random number generators
  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<T> uniform_dist(-1.0f, 1.0f);
  
  // Initialize m initial points
  for (size_t i = 0; i < m; ++i) {
    T* coords = new T[d];
    for (int j = 0; j < d; ++j) {
      coords[j] = uniform_dist(rng);
    }
    points.emplace_back(coords, i, Point::parameters(d));
  }
  
  // Step 2: Generate remaining points using Metropolis-Hastings
  for (size_t i = m; i < n; ++i) {
    T* new_coords = metropolis_hastings_sampling(points, d, sigma, hi, burn_in, thinning, rng);
    points.emplace_back(new_coords, i, Point::parameters(d));
  }
  
  // (Optional) Output the points to a file
  if (oFile != nullptr) {
    std::ofstream outFile(oFile, std::ios::binary);
    if (!outFile) {
      std::cerr << "Error opening output file: " << oFile << std::endl;
      return 1;
    }
    // Write number of points and dimensions
    unsigned int num_points = static_cast<unsigned int>(n);
    unsigned int dims = static_cast<unsigned int>(d);
    outFile.write(reinterpret_cast<char*>(&num_points), sizeof(unsigned int));
    outFile.write(reinterpret_cast<char*>(&dims), sizeof(unsigned int));
    // Write point data
    for (const auto& point : points) {
      outFile.write(reinterpret_cast<const char*>(point.values), d * sizeof(T));
    }
    outFile.close();
  }
  
  // Clean up dynamically allocated memory
  for (auto& point : points) {
    delete[] point.values;
  }
  
  return 0;
}

