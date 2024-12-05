// Tobias Rubel 2024
// UMD Parallel Algorithms Group
// generate random points in a unit cube using a heat kernel
// Include necessary headers
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/random.h"
#include "parlay/internal/get_time.h"
#include "../../ParlayANN/algorithms/utils/types.h"
#include "../../ParlayANN/algorithms/utils/euclidian_point.h"
#include "../../ParlayANN/algorithms/utils/point_range.h"
#include "../../ParlayANN/algorithms/bench/parse_command_line.h"

/*
    // a refresher on how to generate random vectors...
    parlay::sequence<parlay::sequence<float>> random_vectors; 
    for (size_t i = 0; i < num_vectors; i++) {
      parlay::random_generator gen(i+(repl*num_vectors));
      //std::uniform_real_distribution<float> dis(-1.0, 1.0);
      std::normal_distribution<float> dis(0.0, 1.0);
        auto vec = parlay::tabulate(d, [&](size_t i) {
        auto r = gen[i];
        return dis(r);
        });
      random_vectors.push_back(vec);  
    }
*/


// Heat function
float compute_heat(const Euclidian_Point<float>& x, const parlay::sequence<Euclidian_Point<float>>& points, float sigma, float hi) {
  // float total_heat = 0.0f;
  float sigma_sq = sigma * sigma;
  // for (const auto& point : points) {
  //   // Use the optimized squared distance function
  //   float distance_sq = std::pow(x.distance(point),2);
  //   float heat_contribution = hi * std::exp(-distance_sq / (2.0f * sigma_sq));
  //   total_heat += heat_contribution;
  // }
  float total_heat = parlay::reduce(parlay::map(points, [&](const Euclidian_Point<float>& point) {
    // Use the optimized squared distance function
    float distance_sq = std::pow(x.distance(point),2);
    float heat_contribution = hi * std::exp(-distance_sq / (2.0f * sigma_sq));
    return heat_contribution;
  }), parlay::addm<float>());
  return total_heat;
}

// Metropolis-Hastings sampling function
parlay::sequence<Euclidian_Point<float>> metropolis_hastings_sampling(const parlay::sequence<Euclidian_Point<float>>& points, int d, float sigma, float hi,
                                                                      int burn_in, int thinning, size_t npit) {
  // // Initialize random number generator with a unique seed
  // std::mt19937 rng();
  // parlay::random_generator gen(rng);

  // std::uniform_real_distribution<float> uniform_dist(-1.0f, 1.0f);

    // Initialize random number generator with a unique seed
  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  std::uniform_real_distribution<float> uniform_dist(-1.0f, 1.0f);

  // Generate initial point coordinates
  // auto x_current_coords = parlay::tabulate(d, [&](size_t i) {
  //   auto r = gen[i];
  //   auto val = uniform_dist(r);
  //   return std::min(val,1.0)
  // });
  std::vector<float> x_current_coords(d);
  for (int i = 0; i < d; ++i) {
    x_current_coords[i] = uniform_dist(rng);
  }

  Euclidian_Point<float> x_current(x_current_coords.data(), -1, Euclidian_Point<float>::parameters(d));

  float heat_current = compute_heat(x_current, points, sigma, hi);

  // Proposal standard deviation
  float proposal_std = sigma * 0.5f;
  std::normal_distribution<float> normal_dist(0.0f, proposal_std);

  // Burn-in phase
  for (int iteration = 0; iteration < burn_in; ++iteration) {
    // Propose a new point
    std::vector<float> x_proposal_coords(d);
    for (int i = 0; i < d; ++i) {
      float proposal = x_current.get_values()[i] + normal_dist(rng);
      // Reflect back into unit cube if outside
      if (proposal > 1.0f) proposal = 1.0f - (proposal - 1.0f);
      if (proposal < -1.0f) proposal = -1.0f - (proposal + 1.0f);
      // Ensure within [-1, 1]
      proposal = std::max(std::min(proposal, 1.0f), -1.0f);
      x_proposal_coords[i] = proposal;
    }

    Euclidian_Point<float> x_proposal(x_proposal_coords.data(), -1, Euclidian_Point<float>::parameters(d));

    float heat_proposal = compute_heat(x_proposal, points, sigma, hi);

    // Compute acceptance probability
    float acceptance_ratio = heat_proposal / (heat_current + 1e-8f); // Avoid division by zero
    float acceptance_probability = std::min(1.0f, acceptance_ratio);

    // Accept or reject the proposal
    std::uniform_real_distribution<float> accept_dist(0.0f, 1.0f);
    if (accept_dist(rng) < acceptance_probability) {
      // Accept proposal
      x_current = x_proposal;
      heat_current = heat_proposal;
    }
    // Else, x_current and heat_current remain unchanged
  }

  // Sampling phase
  parlay::sequence<Euclidian_Point<float>> samples;
  int samples_collected = 0;
  int iteration_since_last_sample = 0;

  while (samples_collected < npit) {
    // Propose a new point
    std::vector<float> x_proposal_coords(d);
    for (int i = 0; i < d; ++i) {
      float proposal = x_current.get_values()[i] + normal_dist(rng);
      // Reflect back into unit cube if outside
      if (proposal > 1.0f) proposal = 1.0f - (proposal - 1.0f);
      if (proposal < -1.0f) proposal = -1.0f - (proposal + 1.0f);
      // Ensure within [-1, 1]
      proposal = std::max(std::min(proposal, 1.0f), -1.0f);
      x_proposal_coords[i] = proposal;
    }

    Euclidian_Point<float> x_proposal(x_proposal_coords.data(), -1, Euclidian_Point<float>::parameters(d));

    float heat_proposal = compute_heat(x_proposal, points, sigma, hi);

    // Compute acceptance probability
    float acceptance_ratio = heat_proposal / (heat_current + 1e-8f); // Avoid division by zero
    float acceptance_probability = std::min(1.0f, acceptance_ratio);

    // Accept or reject the proposal
    std::uniform_real_distribution<float> accept_dist(0.0f, 1.0f);
    if (accept_dist(rng) < acceptance_probability) {
      // Accept proposal
      x_current = x_proposal;
      heat_current = heat_proposal;
    }
    // Else, x_current and heat_current remain unchanged

    iteration_since_last_sample++;

    if (iteration_since_last_sample >= thinning) {
      // Collect sample
      samples.push_back(x_current);
      samples_collected++;
      iteration_since_last_sample = 0;
    }
  }

  return samples;
}

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
    "[-out <outfile>] [-n <n>] [-d <d>] [-m <m>] [-sigma <sigma>] [-hi <hi>] [-burn_in <burn_in>] [-thinning <thinning>] [-npit <npit>]");

  char* oFile = P.getOptionValue("-out");
  int d = P.getOptionIntValue("-d", 128);
  size_t n = P.getOptionLongValue("-n", 1000);
  size_t m = P.getOptionLongValue("-m", 10);

  // Parse float parameters using getOptionFloatValue
  float sigma = P.getOptionFloatValue("-sigma", 0.1f);
  float hi = P.getOptionFloatValue("-hi", 1.0f);
  int burn_in = P.getOptionIntValue("-burn_in", 40);
  int thinning = P.getOptionIntValue("-thinning", 10);
  size_t npit = P.getOptionLongValue("-npit", 1); // Default npit is 1

  // Prepare data structures
  parlay::sequence<Euclidian_Point<float>> points;
  points.reserve(n); // Reserve space for n points

  // Random number generator
  std::mt19937 rng(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<float> uniform_dist(-1.0f, 1.0f);

  // Initialize m initial points
  for (size_t i = 0; i < m; ++i) {
    std::vector<float> coords(d);
    for (int j = 0; j < d; ++j) {
      coords[j] = uniform_dist(rng);
    }

    Euclidian_Point<float> point(coords.data(), i, Euclidian_Point<float>::parameters(d));
    points.push_back(point);
  }

  // Generate remaining points using Metropolis-Hastings
  size_t i = m;
  while (i < n) {
    size_t points_to_generate = std::min(npit, n - i);
    if (i % 1000 == 0) std::cout << "Generated " << i << " points" << std::endl;
    parlay::sequence<Euclidian_Point<float>> new_points = metropolis_hastings_sampling(points, d, sigma, hi, burn_in, thinning, points_to_generate);
    for (size_t j = 0; j < new_points.size(); ++j) {
      new_points[j].set_id(i);
      points.push_back(new_points[j]);
      ++i;
    }
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
      outFile.write(reinterpret_cast<const char*>(point.get_values()), d * sizeof(float)); // Use get_values() accessor
    }
    outFile.close();
  }

  // No explicit memory cleanup is needed if Euclidian_Point manages its own memory properly

  return 0;
}
