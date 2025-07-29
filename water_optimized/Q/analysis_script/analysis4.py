import matplotlib.pyplot as plt

methods = ['RTX 3070', 'RTX 4090', 'RTX 5090', 'A100']

times_sec = [1612, 7*60+32, 5*60+40, 4*60+21]  # seconds

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(methods, times_sec)

for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 20, f'{int(yval)}s', ha='center', va='bottom')

ax.set_ylabel('Time (seconds)')
ax.set_title(f'Different GPU Performance Comparison\n')

plt.tight_layout()
plt.show()