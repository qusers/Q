import matplotlib.pyplot as plt
import numpy as np

float_times = [406.00, 147.13, 108.78, 184.99]
double_times = [1578, 542, 396, 308]
devices = ['RTX 3070', 'RTX 4090', 'RTX 5090', 'A100']

x = np.arange(len(devices)) 
width = 0.35 

fig, ax = plt.subplots()

rects1 = ax.bar(x - width/2, float_times, width, label='Float', alpha=0.8)
rects2 = ax.bar(x + width/2, double_times, width, label='Double', alpha=0.8)

ax.set_ylabel('Time (s)')
ax.set_xlabel('GPU Devices')
ax.set_title('Kernel Performance Comparison (Float vs Double)')
ax.set_xticks(x)
ax.set_xticklabels(devices)
ax.legend()

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.1f}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 上移 3
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

plt.tight_layout()
plt.show()
