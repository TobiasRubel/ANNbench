import matplotlib.pyplot as plt
import sys


def main(argv):
    input_file = argv[1]
    output_plot = argv[2]

    with open(input_file, "r") as f:
        coefficients = [float(line.strip()) for line in f]

    coefficients.sort()

    plt.figure(figsize=(10, 6))
    plt.plot(coefficients, label="Clustering Coefficients", linewidth=2)
    plt.title("Distribution of Clustering Coefficients", fontsize=16)
    plt.xlabel("Sample Index (sorted by value)", fontsize=14)
    plt.ylabel("Clustering Coefficient", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=12)

    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()

    print(f"Plot saved to {output_plot}")

if __name__ == "__main__":
    main(sys.argv)