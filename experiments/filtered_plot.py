import matplotlib
import math
import sys
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap

res_datasets_path = sys.argv[1]
csv_separator = ";"

input_file = open(res_datasets_path, "r")
lines = [l.split(csv_separator) for l in input_file.read().splitlines()]
input_file.close()


def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier + 0.5) / multiplier


labels = ["\n".join(wrap(l[0], 15)) for l in lines[1:]]
print(labels)
label_method = lines[0][1:]
n = len(label_method)
res = []
for i in range(n):
    res.append([round_half_up(int(l[i + 1]) / 10 ** 6, 2) for l in lines[1:]])
    print(label_method[i] + ":", res[i])

x = np.arange(len(labels))  # the label locations
width = 0.9  # the width of the bars
fig, ax = plt.subplots()

rects = []
for i in range(n):
    rects.append(
        ax.bar(x + (i - n // 2) * (width / n), res[i], width / n, label=label_method[i])
    )
# rects2 = ax.bar(x, synthetic_EBWT_RLO_no_dollar, width / 3, label="RLO+EBWT no $")
# rects3 = ax.bar(x + width / 3, synthetic_XBWT, width / 3, label="XBWT")


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel("Millions of runs")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate(
            "{}".format(height),
            xy=(rect.get_x() + rect.get_width() / 2, height),
            xytext=(0, 3),  # 3 points vertical offset
            fontsize=6,
            textcoords="offset points",
            ha="center",
            va="bottom",
        )


for r in rects:
    autolabel(r)

fig.tight_layout()
fig.savefig("fig_" + res_datasets_path.split(".")[0] + ".png")
# plt.show()
