# Data generated in the likelihood_surface.ipynb notebook

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

REC = 1.25e-8


def plot_single(ax, df, title, include_legend=False):
    print(list(df))

    x = df.index
    for col in df:
        if col.startswith("poly"):
            splits = col.split("_")
            arg = "_".join(splits[1:-1])
            poly = splits[-1]
        else:
            arg = col
            if arg == "full_arg_hudson":
                continue
            poly = "0"

        y = df[col]
        line = ax.semilogx(x, y, label=f"{float(poly):.0%}")

    if include_legend:
        ax.legend(title="Polytomies")

    ax.set_title(title)

    ax.set_ylabel("Log likelihood ratio")
    ax.axvline(REC, linestyle="dashed", color="grey")


def plot(df_rat, path):

    fig, ax = plt.subplots(1, 1)

    colours = {
        "full_arg_smc": "tab:red",
        "unary_coal": "tab:blue",
        "fully_simplified": "tab:orange",
        "tsinfer_tsdate": "tab:olive",
        "tsinfer_fs_tsdate": "tab:cyan",
        "full_arg_hudson": "tab:purple",
    }
    titles = {
        "full_arg_smc": "Full ARG",
        "unary_coal": "RE nodes removed",
        "fully_simplified": "Fully simplified",
        "full_arg_hudson": "Full ARG (HKY)",
        "tsinfer_tsdate": "tsinfer+tsdate",
        "tsinfer_fs_tsdate": "tsinfer+tsdate (FS)",
    }

    line_styles = {
        "0": "solid",
        "0.05": "dashed",
        "0.25": "dotted",
        "0.5": "dashdot",
        "0.75": "-",
    }

    used_line_styles = set()

    x = df_rat.index
    for col in df_rat:
        if col.startswith("poly"):
            splits = col.split("_")
            arg = "_".join(splits[1:-1])
            poly = splits[-1]
        else:
            arg = col
            if arg == "full_arg_hudson":
                continue
            poly = "0"
        marker = None
        if arg == "full_arg_smc":
            marker = "o"

        used_line_styles.add(poly)

        y = df_rat[col]  # / df_rat[arg]
        line = ax.semilogx(
            x,
            y,
            color=colours[arg],
            ls=line_styles[poly],
            marker=marker,
            markersize=2,
            label="_nolegend_" if poly != "0" else titles[arg],
        )

    ls_patches = [
        mlines.Line2D(
            [0],
            [0],
            color="black",
            lw=1,
            ls=line_styles[label],
            label=f"{float(label):.0%}",
        )
        for label in sorted(used_line_styles)
    ]

    legend1 = ax.legend(loc="upper right")
    ax.legend(handles=ls_patches, loc="lower left", title="Polytomies")
    ax.annotate(
        "True parameter", xy=(REC, 18000), xytext=(-80, 0), textcoords="offset points"
    )
    ax.add_artist(legend1)

    ax.set_ylabel("Log likelihood ratio")
    ax.axvline(REC, linestyle="dashed", color="grey")
    ax.set_xlabel("Recombination rate")
    # ax.set_ylim(-3000, 25000)

    plt.savefig(path)


df = pd.read_csv("likelihood_surface.csv").set_index("rho")

df_rat = df.iloc[1:].copy()
for col in list(df_rat):
    df_rat[col] -= df.iloc[0][col]
    df_rat[col] *= -1
    # df_rat[col] /= df.iloc[0][col]

df1 = df_rat[[x for x in list(df_rat) if not (x.split("_")[-1] in ("0.5", "0.75"))]]

plot(df1, "../likelihood_surface.pdf")
plot(df1, "../likelihood_surface.png")

# Drop the 0.05 line as it's so close to full.
df_rat = df_rat[[x for x in list(df_rat) if not (x.split("_")[-1] in ["0.05"])]]

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(4, 8))
df2 = df_rat[[x for x in list(df_rat) if "fully" in x]]
plot_single(axes[2], df2, "Fully simplified")
df2 = df_rat[[x for x in list(df_rat) if "full_arg_smc" in x]]
plot_single(axes[0], df2, "Full ARG", include_legend=True)
df2 = df_rat[[x for x in list(df_rat) if "unary" in x]]
plot_single(axes[1], df2, "RE nodes removed")
axes[2].set_xlabel("Recombination rate")

plt.tight_layout()
# plt.savefig("../likelihood_surface_suppl.png")
plt.savefig("../likelihood_surface_suppl.pdf")
