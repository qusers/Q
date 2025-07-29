import matplotlib.pyplot as plt
import numpy as np

# Data
gpus = ['A100', '5090', '4090', '3070']
new_performance = [49.87, 17.01, 15.12, 3.8]
old_performance = [8, 1.88, 1.37, 0.34]

# reverse
gpus.reverse()
new_performance.reverse()
old_performance.reverse()

# Bar width
bar_width = 0.35

index = np.arange(len(gpus))


# Plot with values on top of bars
fig, ax = plt.subplots(figsize=(8, 5))

# Plot bars
bar1 = ax.bar(index, new_performance, bar_width, label='New')
bar2 = ax.bar(index + bar_width, old_performance, bar_width, label='Old')

# Labels and title
ax.set_xlabel('GPU Models')
ax.set_ylabel('ns/day')
ax.set_title('Performance Comparison of GPUs (New vs Old)')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(gpus)
ax.legend()

# Display values on top of bars
for bar in bar1:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, yval + 0.5, round(yval, 2), ha='center', va='bottom')

for bar in bar2:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, yval + 0.5, round(yval, 2), ha='center', va='bottom')

# Show plot
plt.tight_layout()
plt.show()
