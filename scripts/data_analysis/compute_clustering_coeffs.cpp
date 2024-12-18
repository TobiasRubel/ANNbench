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

inline val_type euclidean_distance_sq(const val_type* a, const val_type* b, uint32_t num_dims) {
    val_type dist = 0;
    for (uint32_t i = 0; i < num_dims; i++) {
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return dist;
}

parlay::sequence<uint32_t> find_knn(const parlay::sequence<val_type>& points, int point_idx, int k, uint32_t num_vecs, uint32_t num_dims) {
    const val_type* query = &points[point_idx * num_dims];

    auto distances = parlay::tabulate<std::pair<uint32_t, val_type>>(num_vecs, [&] (size_t i) {
        const val_type* neighbor = &points[i * num_dims];
        val_type dist = euclidean_distance_sq(query, neighbor, num_dims);
        return std::make_pair(i, dist);
    });

    std::partial_sort(distances.begin(), distances.begin() + k, distances.end(), [&] (auto a, auto b) {
        return a.second < b.second;
    });
    auto neighbors = parlay::tabulate<ind_type>(k, [&] (size_t i) {
        return distances[i + 1].first;
    });
    return neighbors;
}

float compute_local_clustering(const parlay::sequence<val_type>& points, uint32_t index, uint32_t num_vecs, uint32_t num_dims, int k) {
    int connections = 0;

    auto neighbors = find_knn(points, index, k, num_vecs, num_dims);
    neighbors.push_back(index);
    std::unordered_set<uint32_t> neighbor_set(neighbors.begin(), neighbors.end());

    for (int i = 0; i < neighbors.size(); i++) {
        auto other_neighbors = find_knn(points, neighbors[i], k, num_vecs, num_dims);
        for (int other_neighbor : other_neighbors) {
            if (neighbor_set.find(other_neighbor) != neighbor_set.end()) {
                connections++;
            }
        }
    }

    return static_cast<float>(connections) / (k * k);
}

int main(int argc, const char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];

    int sample_size = 100;
    int k = 10;

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

    auto sampled_indices = parlay::sequence<uint32_t>::uninitialized(sample_size);
    for (int i = 0; i < sample_size; i++) {
        sampled_indices[i] = rand() % num_vecs;
    }

    auto clustering_coeffs = parlay::sequence<float>::uninitialized(sample_size);
    parlay::parallel_for(0, sample_size, [&] (size_t i) {
        clustering_coeffs[i] = compute_local_clustering(point_set, sampled_indices[i], num_vecs, num_dims, k);
    }, 1);

    float total_clustering_coeff = 0;
    for (int i = 0; i < sample_size; i++) {
        total_clustering_coeff += clustering_coeffs[i];
    }
    std::cout << "Average clustering coefficient: " << total_clustering_coeff / sample_size << std::endl;

    std::ofstream writer(output_file);
    if (!writer.is_open()) {
        std::cerr << "Unable to open file for writing: " << output_file << std::endl;
        return 1;
    }

    for (val_type coeff : clustering_coeffs) {
        writer << coeff << "\n";
    }

    std::cout << "Clustering coefficients written to " << output_file << std::endl;
    return 0;
}
