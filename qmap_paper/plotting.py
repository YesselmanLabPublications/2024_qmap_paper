import matplotlib.pyplot as plt
import numpy as np

from rna_secstruct_design.selection import get_selection, SecStruct
from q_dms_ttr_paper.titration import normalize_data, normalized_hill_equation


def colors_for_sequence(seq: str):
    """
    get colors for a sequence for DMS reactivity
    :param seq: sequence
    :return: list of colors
    """
    colors = []
    for e in seq:
        if e == "A":
            colors.append("red")
        elif e == "C":
            colors.append("blue")
        elif e == "G":
            colors.append("orange")
        else:
            colors.append("green")
    return colors


def find_stretches(nums):
    """
    Identify stretches of consecutive numbers
    :param nums: list of numbers
    :return: list of stretches
    """
    if len(nums) == 0:
        return []
    nums.sort()
    stretches = []
    start = end = nums[0]

    for num in nums[1:]:
        if num == end + 1:
            end = num
        else:
            stretches.append([start, end])
            start = end = num

    stretches.append([start, end])
    return stretches


def fill_between(ax, color, x, y, alpha=0.15, **kwargs):
    ax.fill_between(x, y, color=color, alpha=0.15, zorder=-1)


def trim(content, prime_5, prime_3):
    """
    trims a string or list from the 5' and 3' end
    :param content: string or list
    :param prime_5: number of bases to trim from 5' end
    :param prime_3: number of bases to trim from 3' end
    :return: trimmed string or list
    """
    if isinstance(content, str):
        return content[prime_5 : -prime_3 or None]
    elif isinstance(content, list):
        return content[prime_5 : -prime_3 or len(content)]
    else:
        return "Invalid content type. Please provide a string or a list."


def plot_pop_avg(
    seq,
    ss,
    reactivities,
    ax=None,
    axis="sequence_structure",
    trim_5p=0,
    trim_3p=0,
    highlights=None,
):
    """
    plot dms reactivity for a sequence and secondary structure
    :param seq: sequence
    :param ss: secondary structure
    :param reactivities: list of reactivities
    :param ax: matplotlib axis
    :param axis: axis to plot on
    :param trim_5p: trim 5' end
    :param trim_3p: trim 3' end

    """
    seq = trim(seq, trim_5p, trim_3p)
    ss = trim(ss, trim_5p, trim_3p)
    reactivities = trim(reactivities, trim_5p, trim_3p)
    highlight_bounds = []
    if highlights is None:
        highlights = []
    for h in highlights:
        selection = get_selection(SecStruct(seq, ss), h)
        for bounds in find_stretches(selection):
            highlight_bounds.append(bounds)
    colors = colors_for_sequence(seq)
    x = list(range(len(seq)))
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(20, 4))
    ax.bar(range(0, len(reactivities)), reactivities, color=colors)
    ax.set_xticks(x)
    for bounds in highlight_bounds:
        fill_between(ax, "gray", bounds, [0, 10])
    if axis == "sequence_structure":
        ax.set_xticklabels([f"{s}\n{nt}" for s, nt in zip(seq, ss)])
    elif axis == "sequence":
        ax.set_xticklabels([f"{s}" for s in seq])
    elif axis == "structure":
        ax.set_xticklabels([f"{s}" for s in ss])
    else:
        pass
    return ax


def plot_pop_avg_from_row(row, data_col="data", ax=None):
    return plot_pop_avg(row["sequence"], row["structure"], row[data_col], ax)


def plot_pop_avg_all(df, data_col="data", axis="sequence_structure", **kwargs):
    fig, axes = plt.subplots(len(df), 1, **kwargs)
    j = 0
    for i, row in df.iterrows():
        colors = colors_for_sequence(row["sequence"])
        axes[j].bar(range(0, len(row[data_col])), row[data_col], color=colors)
        axes[j].set_title(row["rna_name"])
        j += 1
    plot_pop_avg_from_row(df.iloc[-1], ax=axes[-1], axis=axis)
    return fig


def plot_pop_avg_titration(df, titration_col, highlights=None, **kwargs):
    fig, axes = plt.subplots(len(df), 1, **kwargs)
    j = 0
    secstruct = SecStruct(df.iloc[0]["sequence"], df.iloc[0]["structure"])
    highlight_bounds = []
    if highlights is None:
        highlights = []
    else:
        for h in highlights:
            selection = get_selection(secstruct, h)
            for bounds in find_stretches(selection):
                highlight_bounds.append(bounds)
    for i, row in df.iterrows():
        colors = colors_for_sequence(row["sequence"])
        axes[j].bar(range(0, len(row["data"])), row["data"], color=colors)
        axes[j].set_title(str(row[titration_col]) + " mM")
        axes[j].set_ylim([0, 0.1])
        axes[j].set_xlim([-0.1, len(row["data"]) + 0.1])
        axes[j].set_xticks([])
        for bounds in highlight_bounds:
            fill_between(axes[j], "gray", bounds, [0, 10])
        j += 1
    plot_pop_avg_from_row(df.iloc[-1], ax=axes[-1])
    return fig


def plot_mg_titration_fit(x, y, mg_1_2, n, max_val, **kwargs):
    norm_data = -normalize_data(np.array(y)) + 1
    fig, ax = plt.subplots(1, 1, **kwargs)
    ax.scatter(x, norm_data, s=100)
    xs, ys = [], []
    for j in np.arange(0, 45, 0.25):
        y = normalized_hill_equation(j, mg_1_2, n, max_val)
        xs.append(j)
        ys.append(y)
    plt.plot(xs, ys, lw=3)
    # plt.fill_between(xs, ys - r[1][0], ys + r[1][0], alpha=0.2, lw=0)
    # plt.ylim(-0.05, 1.1)
    publication_style_ax(ax)
    return ax


# style functions #############################################################


def publication_style_ax(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(width=2)
    fsize = 24
    ax.xaxis.label.set_fontsize(fsize)
    ax.yaxis.label.set_fontsize(fsize)
    ax.tick_params(axis="both", which="major", labelsize=fsize - 2)
