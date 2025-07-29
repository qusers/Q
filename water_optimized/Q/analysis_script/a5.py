import pandas as pd
import matplotlib.pyplot as plt


# Configuring matplotlib to use a specific font
data = {
    'GPU': ['3070', '3070', '4090', '4090', '5090', '5090', 'A100', 'A100'],
    'Version': ['old', 'new'] * 4,
    '20A': [177, 26.848, 49.565, 12.563, 36.832, 11.279, 25.248, 11.844],
    '25A': [592.59, 58.505, 149.554, 22.538, 114.748, 16.786, 52.452, 14.729],
    '30A': [1735.94, 157.85, 450.534, 54.231, 337.916, 39.607, 122.507, 30.776]
}

df = pd.DataFrame(data)

# Use GPU and Version as a multi-index for better organization
df['Label'] = df['GPU'] + ' ' + df['Version']


fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
sizes = ['20A', '25A', '30A']
x_labels = ['3070 old', '3070 new', '4090 old', '4090 new', '5090 old', '5090 new', 'A100 old', 'A100 new']
bar_width = 0.5
colors = ['#FF9999', '#FF3333', '#9999FF', '#3333FF', '#99FF99', '#33CC33', '#FFD700', '#FFA500']

for idx, size in enumerate(sizes):
    ax = axes[idx]
    values = df[size]
    bars = ax.bar(x_labels, values, color=colors, width=bar_width)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 10,
                f'{val:.1f}', ha='center', va='bottom', fontsize=8)

    ax.set_title(f'Runtime for {size} (10000 steps)')
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.set_ylabel('Runtime (s)' if idx == 0 else "")
    ax.grid(True)

plt.tight_layout()
plt.show()# 
# fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
# sizes = ['20A', '25A', '30A']
# x_labels = ['3070 old', '3070 new', '4090 old', '4090 new', '5090 old', '5090 new', 'A100 old', 'A100 new']
# bar_width = 0.35
# colors = ['#FF9999', '#FF3333', '#9999FF', '#3333FF', '#99FF99', '#33CC33', '#FFD700', '#FFA500']

# for idx, size in enumerate(sizes):
#     ax = axes[idx]
#     values = df[size]
#     ax.bar(x_labels, values, color=colors)
#     ax.set_title(f'Runtime for {size} (10000 steps)')
#     ax.set_xticklabels(x_labels, rotation=45, ha='right')
#     ax.set_ylabel('Runtime (s)' if idx == 0 else "")
#     ax.grid(True)

# plt.tight_layout()
# plt.show()
