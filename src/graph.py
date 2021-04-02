import sys
import json

import matplotlib.pyplot as plt

jsn_plt = None
try:
    with open("plot_points.json", 'r') as inp_f:
        jsn_plt = json.load(inp_f)
except:
    print("Error opening file, might need running 'to_plot' first")
    exit(-1)

r = jsn_plt["r"]
ptncl_l = [jsn_plt["AA"], jsn_plt["AB"], jsn_plt["BB"]]
lbl_l = ["A-A", "A-B", "B-B"]

ptncl_min_list = [min(zip(r, u), key= lambda pair: pair[1]) for u in ptncl_l]

fig, ax = plt.subplots(1, 3, figsize=(18, 4))
for ax_i, lbl, u, u_min in zip(ax, lbl_l, ptncl_l, ptncl_min_list):
    ax_i.set_title(lbl)
    ax_i.set_xlim(u_min[0] - 1, u_min[0] + 5)
    ax_i.set_ylim(u_min[1] - 1, u_min[1] + 5)
    ax_i.set_xlabel("r[A]")
    ax_i.set_ylabel("u(r) [eV]")
    ax_i.grid(True)
    ax_i.plot(r, u)

plt.savefig("plots.png")
plt.show()

