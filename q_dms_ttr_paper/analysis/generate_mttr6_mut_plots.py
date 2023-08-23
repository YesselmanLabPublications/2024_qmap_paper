import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


def main():
    """
    main function for script
    """


def plot_tlr_res_against_mg2(df, bounds=[3, 10]):
    pal = sns.color_palette("tab10", 11)[::-1]
    colors = pal
    plt.title(df.iloc[0]["name"] + " " + str(round(df.iloc[0]["mg_1_2"], 2)))
    for i in range(bounds[0], bounds[1]):
        seq = df.iloc[0]["aligned_seq"][i - 1]
        # not DMS active
        if seq == "U" or seq == "G":
            continue
        plt.plot(
            df["mg_conc"],
            df["tlr"].apply(lambda x: x[i - 1]),
            lw=2,
            marker="o",
            markersize=7,
            color=colors[i],
            label=str(i) + " " + seq,
        )
    plt.legend()
    plt.xlabel("mg_conc")
    plt.ylabel("normalized_reactivity")


def plot_pos_against_org_uucg(df, df_org, df_uucg, pos):
    plt.title(df.iloc[0]["name"] + " " + str(round(df.iloc[0]["mg_1_2"], 2)))
    plt.plot(
        df_org["mg_conc"],
        df_org["tlr"].apply(lambda x: x[pos - 1]),
        marker="o",
        label="wt",
    )
    plt.plot(
        df_uucg["mg_conc"],
        df_uucg["tlr"].apply(lambda x: x[pos - 1]),
        marker="o",
        label="uucg",
    )
    plt.plot(
        df["mg_conc"],
        df["tlr"].apply(lambda x: x[pos - 1]),
        marker="o",
        label=df.iloc[0]["name"],
    )
    plt.legend()


def plot_all_res(df):
    df_names = df["name"].unique()
    with PdfPages("all_res.pdf") as pdf:
        for name in df_names:
            group = df[df["name"] == name]
            group = group.sort_values(by=["mg_conc"])
            plot_tlr_res_against_mg2(group)
            pdf.savefig()
            plt.clf()


def main():
    DATA_PATH = "data"
    df = pd.read_json(f"{DATA_PATH}/mttr6_data_full.json")
    df = df[
        ~df["name"].isin(
            [
                "CAUGA_UCUAAA",
                "UAUGG_CUUAAC",
                "UACGG_CCUACA",
                "CACGG_CCUCAC",
                "CAUGC_GCUCAA",
                "CAUGC_GCUGAA",
            ]
        )
    ]
    df_org = pd.read_json(f"{DATA_PATH}/wt_mg_titra.json")
    df_org = df_org[df_org["exp_name"] == "2022_07_27_C0117_50mM_NaC_Mg2+_titra_CM"]
    df_uucg = pd.read_json(f"{DATA_PATH}/uucg_mg_titra.json")
    df = df[df["mg_conc"] != 5.0]
    df = df.sort_values(by=["mg_1_2"])
    df_names = df["name"].unique()
    pos = 8
    with PdfPages("output.pdf") as pdf:
        for name in df_names:
            group = df[df["name"] == name]
            group = group.sort_values(by=["mg_conc"])
            if group.iloc[0]["aligned_seq"][pos - 1] != "A":
                continue
            plot_pos_against_org_uucg(group, df_org, df_uucg, pos)
            pdf.savefig()
            plt.clf()


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
