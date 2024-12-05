#include <iostream>
#include <fstream>
#include <random>
#include <cstring>
#include <cmath>
#include <memory>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "../../ParlayANN/algorithms/utils/euclidian_point.h"
#include "../../ParlayANN/algorithms/utils/point_range.h"
#include "../../ParlayANN/algorithms/bench/parse_command_line.h"

int main(int argc, char* argv[]) {
    commandLine P(argc, argv,
        "-input_file <input> -output_file <output> -np <np> -sigma <sigma>");

    char* input_file = P.getOptionValue("-input_file");
    char* output_file = P.getOptionValue("-output_file");
    int np = P.getOptionIntValue("-np", 1);
    double sigma = P.getOptionDoubleValue("-sigma", 1.0);

    if (!input_file || !output_file) {
        std::cout << "Usage: " << argv[0] << " -input_file <input> -output_file <output> -np <np> -sigma <sigma>" << std::endl;
        return -1;
    }

    std::cout << "Generating " << np << " new points per original point with sigma = " << sigma << std::endl;

    // Read the original PointRange from the input file
    PointRange<float, Euclidian_Point<float>> original_points(input_file);
    size_t dimension = original_points.dimension();
    size_t n = original_points.size();
    size_t total_new_points = n * np;

    size_t aligned_dimension = original_points.aligned_dimension();
    size_t num_bytes = total_new_points * aligned_dimension * sizeof(float);
    float* ptr = (float*) aligned_alloc(64, num_bytes);
    if (!ptr) {
        std::cerr << "Memory allocation failed!" << std::endl;
        return -1;
    }
    madvise(ptr, num_bytes, MADV_HUGEPAGE);
    std::shared_ptr<float[]> new_values(ptr, std::free);

    // Generate new points by adding Gaussian noise
    parlay::parallel_for(0, total_new_points, [&](size_t idx) {
        size_t p_idx = idx / np; // Index of the original point
        Euclidian_Point<float> p = original_points[p_idx];
        float* new_point_values = new_values.get() + idx * aligned_dimension;

        // Seed the RNG differently for each point to ensure randomness
        std::mt19937 local_rng(idx + std::random_device{}());
        std::normal_distribution<float> local_dist(0.0, sigma);

        for (size_t d = 0; d < dimension; ++d) {
            float noise = local_dist(local_rng);
            new_point_values[d] = p[d] + noise;
        }
        // Zero out any padding
        for (size_t d = dimension; d < aligned_dimension; ++d) {
            new_point_values[d] = 0.0f;
        }
    });

    // Write the new points directly to the output file
    std::ofstream writer(output_file, std::ios::binary);
    if (!writer.is_open()) {
        std::cerr << "Failed to open output file: " << output_file << std::endl;
        return -1;
    }

    // Write header (number of points and dimension)
    unsigned int num_points = static_cast<unsigned int>(total_new_points);
    unsigned int dim = static_cast<unsigned int>(dimension);
    writer.write(reinterpret_cast<char*>(&num_points), sizeof(unsigned int));
    writer.write(reinterpret_cast<char*>(&dim), sizeof(unsigned int));

    // Write point data in blocks to manage memory usage
    size_t BLOCK_SIZE = 1000000;
    size_t index = 0;
    while (index < total_new_points) {
        size_t floor = index;
        size_t ceiling = std::min(index + BLOCK_SIZE, total_new_points);
        size_t block_size = (ceiling - floor) * dimension * sizeof(float);

        // Create a temporary buffer to hold the data to be written
        std::unique_ptr<float[]> temp_buffer(new float[(ceiling - floor) * dimension]);

        parlay::parallel_for(floor, ceiling, [&](size_t i) {
            float* source = new_values.get() + i * aligned_dimension;
            float* dest = temp_buffer.get() + (i - floor) * dimension;
            std::memcpy(dest, source, dimension * sizeof(float));
        });

        writer.write(reinterpret_cast<char*>(temp_buffer.get()), block_size);
        index = ceiling;
    }

    writer.close();

    std::cout << "New PointRange with " << total_new_points << " points saved to " << output_file << std::endl;

    return 0;
}
