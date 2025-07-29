import matplotlib.pyplot as plt

data = {
    "RTX3070": {
        20: {
            "symmetric": [2.15,  7.07, 13.86, 21.67, 35.67, 53.99, 73.37, 103.12, 123.33],
            "triangle": [2.29, 7.81, 15.33, 23.39, 38.38, 58.57, 81.96, 108.60, 128.81],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        25: {
            "symmetric": [5.47, 21.75, 42.73, 64.94, 103.48, 126.70, 176.40, 215.27, 241.96],
            "triangle": [5.95, 24.17, 47.84, 71.73, 105.07, 144.22, 191.22, 230.54, 261.30],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        30: {
            "symmetric": [15.36, 67.69, 133.87, 200.91, 278.85, 359.94, 431.70, 518.88, 588.33],
            "triangle": [17.11, 75.12, 149.01, 223.31, 306.21, 400.06, 485.84, 562.88, 642.06],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        }
    },
    "RTX4090": {
        20: {
            "symmetric": [1.14, 2.50, 4.53, 8.69, 13.44, 20.16, 29.56, 36.88, 49.77],
            "triangle": [1.19, 2.81, 4.86, 9.32, 14.96, 21.77, 32.11, 42.59, 55.41],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        25: {
            "symmetric": [2.30, 6.32, 11.71, 23.07, 37.79, 55.75, 77.46, 105.49, 126.36],
            "triangle": [2.48, 6.98, 12.80, 25.77, 41.74, 61.83, 93.31, 109.58, 132.27],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        30: {
            "symmetric": [5.42, 17.97, 34.78, 59.67, 102.40, 114.91, 143.51, 194.50, 241.79],
            "triangle": [6.04, 20.15, 38.65, 68.93, 105.89, 135.73, 175.56, 216.12, 257.25],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        }
    },
    "RTX5090": {
        20: {
            "symmetric": [0.69, 1.73, 3.24, 4.96, 8.10, 11.95, 17.09, 20.23, 21.56],
            "triangle": [0.72, 1.85, 3.49, 5.41, 8.77, 13.01, 18.18, 21.51, 23.41],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        25: {
            "symmetric": [1.42, 4.43, 8.53, 13.37, 22.55, 34.32, 48.23, 55.75, 61.29],
            "triangle": [1.53, 4.94, 9.42, 14.84, 24.93, 37.21, 52.62, 61.00, 67.14],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        30: {
            "symmetric": [3.942, 12.797, 25.308, 38.713, 60.070, 97.792, 118.17, 137.84, 150.44],
            "triangle": [3.90, 14.32, 27.85, 42.50, 68.53, 100.93, 130.72, 147.70, 158.46],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        }
    },
    "A100": {
        20: {
            "symmetric": [0.824, 1.024, 1.459, 1.962, 2.488, 2.817, 3.976, 6.511],
            "triangle": [0.762, 1.104, 1.637, 2.160, 2.720, 3.373, 4.684, 8.201],
            "x": [1, 5, 10, 15, 20, 25, 30, 40]
        },
        25: {
            "symmetric": [1.432, 2.195, 3.342, 4.748, 6.122, 7.417, 9.723, 12.700, 16.441],
            "triangle": [1.503, 2.644, 4.427, 6.228, 7.880, 9.738, 12.348, 16.640, 21.555],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        },
        30: {
            "symmetric": [3.055, 5.638, 9.862, 13.859, 17.846, 22.471, 27.364, 35.880, 47.497],
            "triangle": [3.3968, 7.9428, 14.1236, 20.1772, 26.5501, 33.931, 41.010, 53.699, 66.586],
            "x": [1, 5, 10, 15, 20, 25, 30, 35, 40]
        }
    }
}

def plot_speedup(gpu: str, water_size: int, baseline: float):
    vals_half = data[gpu][water_size]["symmetric"]
    vals_triangle = data[gpu][water_size]["triangle"]
    x = data[gpu][water_size]["x"]

    speedup_half = [baseline / (v / n) for v, n in zip(vals_half, x)]
    speedup_triangle = [baseline / (v / n) for v, n in zip(vals_triangle, x)]

    plt.figure(figsize=(10, 6))
    plt.plot(x, speedup_half, marker='o', linestyle='-', label='Symmetric Rectangular Block')
    plt.plot(x, speedup_triangle, marker='s', linestyle='--', label='Triangle Block')

    for xi, val in zip(x, speedup_half):
        plt.text(xi, val + 0.3, f'{val:.2f}', ha='center', va='bottom', fontsize=8)
    for xi, val in zip(x, speedup_triangle):
        plt.text(xi, val + 0.3, f'{val:.2f}', ha='center', va='bottom', fontsize=8)

    plt.title(f"Speedup for {gpu} at water_{water_size}")
    plt.xlabel("Number of Concurrent Simulations")
    plt.ylabel("Speedup over Original Serial")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_speedup_subplots(gpu: str):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for idx, water_size in enumerate([20, 25, 30]):
        ax = axs[idx]
        vals_half = data[gpu][water_size]["symmetric"]
        x = data[gpu][water_size]["x"]
        baseline = baseline_times[gpu][water_size]
        speedup_half = [baseline / (v / n) for v, n in zip(vals_half, x)]

        vals_triangle = data[gpu][water_size]["triangle"]
        speedup_triangle = [baseline / (v / n) for v, n in zip(vals_triangle, x)]

        ax.plot(x, speedup_half, marker='o', linestyle='-', label='Symmetric Rectangular Block')
        ax.plot(x, speedup_triangle, marker='s', linestyle='--', label='Triangle Block')

        ax.set_title(f"{gpu} water_{water_size}")
        ax.set_xlabel("Concurrent Simulations")
        if idx == 0:
            ax.set_ylabel("Speedup")
        ax.grid(True)
        ax.legend(fontsize=8)

        for xi, val in zip(x, speedup_half):
            ax.text(xi, val + 0.3, f'{val:.1f}', ha='center', va='bottom', fontsize=7)
        for xi, val in zip(x, speedup_triangle):
            ax.text(xi, val + 0.3, f'{val:.1f}', ha='center', va='bottom', fontsize=7)

    plt.suptitle(f"Speedup of {gpu} across water ball sizes", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()



baseline_times = {
    "RTX3070": {20: 17.7, 25: 59.26, 30: 173.59},
    "RTX4090": {20: 4.96, 25: 14.96, 30: 45.05},
    "RTX5090": {20: 3.68, 25: 11.47, 30: 33.79},
    "A100": {20: 2.52, 25: 5.25, 30: 12.25}
}
# plot_speedup_subplots("RTX3070")
# plot_speedup_subplots("RTX4090")
# plot_speedup_subplots("RTX5090")
# plot_speedup_subplots("A100")

def plot_all_gpus_all_waterballs():
    gpus = ["RTX3070", "RTX4090", "RTX5090", "A100"]
    water_sizes = [20, 25, 30]

    fig, axs = plt.subplots(len(gpus), len(water_sizes), figsize=(18, 12), sharey=True)

    for row, gpu in enumerate(gpus):
        for col, water_size in enumerate(water_sizes):
            ax = axs[row, col]
            vals_half = data[gpu][water_size]["symmetric"]
            # vals_triangle = data[gpu][water_size]["triangle"]
            x = data[gpu][water_size]["x"]
            baseline = baseline_times[gpu][water_size]

            speedup_half = [baseline / (v / n) for v, n in zip(vals_half, x)]
            # speedup_triangle = [baseline / (v / n) for v, n in zip(vals_triangle, x)]

            ax.plot(x, speedup_half, marker='o', linestyle='-', label='Symmetric')
            # ax.plot(x, speedup_triangle, marker='s', linestyle='--', label='Triangle')

            ax.set_title(f"{gpu} water_{water_size}")

            if col == 0:
                ax.set_ylabel("Speedup")
            if row == len(gpus) - 1:
                ax.set_xlabel("Concurrent Simulations")

            ax.grid(True)
            ax.legend(fontsize=7)

            for xi, val in zip(x, speedup_half):
                ax.text(xi, val + 0.3, f'{val:.1f}', ha='center', va='bottom', fontsize=6)
            # for xi, val in zip(x, speedup_triangle):
            #     ax.text(xi, val + 0.3, f'{val:.1f}', ha='center', va='bottom', fontsize=6)

    plt.suptitle("Speedup of All GPUs at Different Water Ball Sizes", fontsize=18)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

plot_all_gpus_all_waterballs()
# plot_speedup_by_gpu("RTX3070")


rtx3070_baseline_time_20 = 17.7
rtx3070_baseline_time_25 = 59.26
rtx3070_baseline_time_30 = 173.59
# plot_speedup("RTX3070", 20, rtx3070_baseline_time_20)
# plot_speedup("RTX3070", 25, rtx3070_baseline_time_25)
# plot_speedup("RTX3070", 30, rtx3070_baseline_time_30)

# rtx4090_baseline_time_20 = 4.96
# rtx4090_baseline_time_25 = 14.96
# rtx4090_baseline_time_30 = 45.05
# plot_speedup("RTX4090", 20, rtx4090_baseline_time_20)
# plot_speedup("RTX4090", 25, rtx4090_baseline_time_25)
# plot_speedup("RTX4090", 30, rtx4090_baseline_time_30)



# rtx5090_baseline_time_20 = 3.68
# rtx5090_baseline_time_25 = 11.47
# rtx5090_baseline_time_30 = 33.79
# plot_speedup("RTX5090", 20, rtx5090_baseline_time_20)
# plot_speedup("RTX5090", 25, rtx5090_baseline_time_25)
# plot_speedup("RTX5090", 30, rtx5090_baseline_time_30)



# rtxa100_baseline_time_20 = 2.52
# rtxa100_baseline_time_25 = 5.25
# rtxa100_baseline_time_30 = 12.25

# plot_speedup("A100", 20, rtxa100_baseline_time_20)
# plot_speedup("A100", 25, rtxa100_baseline_time_25)
# plot_speedup("A100", 30, rtxa100_baseline_time_30)

 