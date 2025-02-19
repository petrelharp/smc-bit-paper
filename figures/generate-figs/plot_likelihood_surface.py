# Data generated in the likelihood_surface.ipynb notebook

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

REC = 1.25e-8

df_rat = pd.read_csv("likelihood_surface.csv").set_index("rho")

fig, ax = plt.subplots(1, 1)

colours = {
    "full_arg_smc": "tab:red",
    "unary_coal": "tab:blue",
    "fully_simplified": "tab:orange",
    "full_arg_hudson": "tab:purple",
}
titles = {
    "full_arg_smc": "Full ARG (SMC)",
    "unary_coal": "RE nodes removed",
    "fully_simplified": "Fully simplified",
    "full_arg_hudson": "Full ARG (HKY)",
}

line_styles = {
    "0": "solid",
    "0.05": "dashed",
    "0.25": "dotted",
}

ls_patches = [
    mlines.Line2D([0], [0], color="black", lw=1, ls=ls, label=f"{float(label):.0%}")
    for label, ls in line_styles.items()
]

x = df_rat.index
for col in df_rat:
    if col.startswith("poly"):
        splits = col.split("_")
        arg = "_".join(splits[1:-1])
        poly = splits[-1]
    else:
        arg = col
        poly = "0"

    y = df_rat[col]  # / df_rat[arg]
    line = ax.semilogx(
        x,
        y,
        color=colours[arg],
        ls=line_styles[poly],
        label="_nolegend_" if poly != "0" else titles[arg],
    )

legend1 = ax.legend()
ax.legend(handles=ls_patches, loc="lower left", title="Polytomies")
ax.annotate(
    "True parameter", xy=(REC, 1.1), xytext=(-80, 0), textcoords="offset points"
)
ax.add_artist(legend1)

ax.set_ylabel("Relative log likelihood")
ax.axvline(REC, linestyle="dashed", color="grey")
ax.set_xlabel("Recombination rate")

plt.savefig("../likelihood_surface.pdf")
