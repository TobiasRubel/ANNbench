import numpy as np
import matplotlib.pyplot as plt

def read_distances(file_path):
    with open(file_path, 'rb') as f:
        sample_size = np.frombuffer(f.read(4), dtype=np.int32)[0]
        all_distances = []
        for _ in range(sample_size):
            count = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            distances = np.frombuffer(f.read(count * 4), dtype=np.float32)
            all_distances.append(distances)
    return all_distances

def plot_distances(all_distances, file_path):
    plt.figure(figsize=(10, 6))
    for distances in all_distances:
        plt.plot(distances, alpha=0.6, linewidth=0.8)
    plt.xlabel("Index (sorted distances)")
    plt.ylabel("Distance")
    plt.title("Distribution of Distances for Sampled Points")
    plt.grid(True)
    plt.savefig(file_path, dpi=300)

if __name__ == "__main__":
    input_path = "distance_distr.bin"
    output_path = "distance_distr.png"
    all_distances = read_distances(input_path)
    plot_distances(all_distances, output_path)
