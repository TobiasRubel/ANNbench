#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>
#include <utility>

#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

using ind_type = uint32_t;
using val_type = float;

inline val_type euclidean_distance(const val_type* a, const val_type* b, uint32_t num_dims) {
    val_type dist = 0;
    for (uint32_t i = 0; i < num_dims; i++) {
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return std::sqrt(dist);
}

parlay::sequence<val_type> sorted_dists(const parlay::sequence<val_type>& points, ind_type point_idx, uint32_t num_vecs, uint32_t num_dims) {
    const val_type* query = &points[point_idx * num_dims];

    auto distances = parlay::tabulate<val_type>(num_vecs, [&] (size_t i) {
        const val_type* neighbor = &points[i * num_dims];
        return euclidean_distance(query, neighbor, num_dims);
    });

    std::sort(distances.begin(), distances.end());
    return distances;
}

int main(int argc, const char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];

    int sample_size = 100;

    std::ifstream reader(input_file);
    if (!reader.is_open()) {
        std::cerr << "Unable to open file for reading: " << input_file << std::endl;
        return 0;
    }

    uint32_t num_vecs, num_dims;
    reader.read((char*)(&num_vecs), sizeof(uint32_t));
    reader.read((char*)(&num_dims), sizeof(uint32_t));

    std::cout << "Vecs: " << num_vecs << "\tDims: " << num_dims << std::endl;

    auto point_set = parlay::sequence<val_type>::uninitialized(num_vecs * num_dims);
    reader.read((char*)(&point_set[0]), num_vecs * num_dims * sizeof(val_type));

    auto sampled_indices = parlay::sequence<ind_type>::uninitialized(sample_size);
    for (int i = 0; i < sample_size; i++) {
        sampled_indices[i] = rand() % num_vecs;
    }

    auto distances = parlay::sequence<parlay::sequence<float>>::uninitialized(sample_size);
    parlay::parallel_for(0, sample_size, [&] (size_t i) {
        distances[i] = sorted_dists(point_set, sampled_indices[i], num_vecs, num_dims);
    }, 1);

    std::ofstream writer(output_file);
    if (!writer.is_open()) {
        std::cerr << "Unable to open file for writing: " << output_file << std::endl;
        return 1;
    }

    writer.write((char*)(&sample_size), sizeof(int));
    for (const auto& dist : distances) {
        uint32_t count = dist.size();
        writer.write((char*)(&count), sizeof(uint32_t));
        writer.write((char*)(&dist[0]), count * sizeof(val_type));
    }

    std::cout << "Clustering coefficients written to " << output_file << std::endl;
    return 0;
}
