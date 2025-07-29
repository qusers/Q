import matplotlib.pyplot as plt

methods = ['Old Method', 'Triangular Method', 'Rectangular Method(Optimized)']
old = 20*60 + 1
times_sec = [20*60+1, 5*60+40, 4*60+20]  # seconds

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(methods, times_sec)

for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 20, f'{int(yval)}s, x{old / yval:.1f}', ha='center', va='bottom')

ax.set_ylabel('Time (seconds)')
ax.set_title(f'Performance Comparison\n')

plt.tight_layout()
plt.show()