import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from pathlib import Path

# Set font and size
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 24

# Data extracted from the table
molecules = ["DG", "FC3", "I1P"]

crest_data = {
    "DG": [375, None, None, None],
    "FC3": [118, None, None, None],
    "I1P": [2167, None, None, None],
}

conquest_data = {
    "DG": [1386, 400, 412, 311],
    "FC3": [150, 128, 134, 83],
    "I1P": [3443, 2230, 617, 637],
}

nodes = [1, 2, 4, 8]

# Calculating the speedup
speedup_data = {}

for molecule in molecules:
    crest_time_1_node = crest_data[molecule][0]
    speedup_data[molecule] = [
        crest_time_1_node / conquest_time if conquest_time else None
        for conquest_time in conquest_data[molecule]
    ]

fig, ax = plt.subplots(figsize=(10, 8))

# Add a black horizontal line through 1
ax.axhline(y=1, color="black", linewidth=3)

# Set color palette
palette = sns.color_palette("Dark2", len(molecules))

for molecule, color in zip(molecules, palette):
    plt.plot(
        nodes,
        speedup_data[molecule],
        "-o",
        label=f"{molecule}",
        linewidth=3,
        color=color,
        markersize=10,
    )

# Setting up the plot
ax.set_xlabel("Number of Nodes")
ax.set_ylabel("Speedup")
ax.legend()
ax.grid(True, which="both", linewidth=1)
ax.set_yscale("log")  # Using a logarithmic scale for the y-axis
ax.set_yticks([0.2, 0.5, 1, 2, 3])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xscale("log")
ax.set_xticks(nodes)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

plt.tight_layout()
plt.savefig(Path("graphs/speedup.pdf"), dpi=300, bbox_inches="tight")
