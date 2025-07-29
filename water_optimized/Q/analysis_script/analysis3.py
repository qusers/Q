# Re-import required libraries after kernel reset
import matplotlib.pyplot as plt

# Re-define the data
vals = [318, 5 * 60 + 39, 425, 8 * 60 + 23, 9 * 60 + 34, 11 * 60 + 4,
        11 * 60 + 58, 13 * 60 + 39, 15 * 60 + 2, 16 * 60 + 45,
        17 * 60 + 49, 19 * 60 + 39, 20 * 60 + 41, 22 * 60 + 17,
        23 * 60 + 31, 30 * 60 + 35]

# vals2 = [33.968, 45.208, 56.08, 60 + 8.412, 60 + 19.428, 60 + 29.122,
#          60 + 42.016, 120 + 21.236, 180 + 21.772, 240 + 25.501]

vals2 = [8.24, 9.01, 9.07, 9.67, 10.24, 11.36, 12.02, 14.59, 19.62, 24.88, 28.17, 39.76, 65.11]

x = [i for i in range(1, 16)]
x.append(20)

x2 = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20, 25, 30, 40]

# baseline times
old = 20 * 60 + 1  # for vals
old2 = 25.2       # for vals2

# Convert times to speedup ratios
vals = [old / (v / n) for v, n in zip(vals, x)]
vals2 = [old2 / (v / n) for v, n in zip(vals2, x2)]

plt.figure(figsize=(10, 6))
plt.plot(x, vals, marker='o', linestyle='-', label='Symmetric Rectangular Block')
plt.plot(x2, vals2, marker='s', linestyle='--', label='Triangle Block')

for xi, val in zip(x, vals):
    plt.text(xi, val + 0.3, f'{val:.2f}', ha='center', va='bottom', fontsize=8)
for xi, val in zip(x2, vals2):
    plt.text(xi, val + 0.3, f'{val:.2f}', ha='center', va='bottom', fontsize=8)

plt.title("Speedup vs. Number of Concurrent Simulations (MPS)")
plt.xlabel("Number of Concurrent Simulations")
plt.ylabel("Speedup over Original Serial")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
