import pandas as pd
import matplotlib.pyplot as plt
from q_dms_ttr_paper.paths import RESOURCES_PATH


def main():
    """
    main function for script
    """
    df = pd.read_csv("q_dms_ttr_paper/resources/csvs/ttr_mutation_dgs_subset.csv")
    df_rna_map = pd.read_csv(
        "/Users/jyesselman2/Dropbox/projects/construct_design/2022/2022_01_27_steve_ttrs/raw_data/tectorna_results_tertcontacts.csv"
    )
    df_rna_map = df_rna_map[df_rna_map["length"] == 10]
    df_rna_map = df_rna_map[df_rna_map["junction_seq"] == "_"]
    df_rna_map = df_rna_map[df_rna_map["b_name"] == "normal"]
    df_rna_map.to_csv("test.csv", index=False)
    data = []
    for i, row in df.iterrows():
        df_sub = df_rna_map[df_rna_map["r_seq"] == row["name"]]
        data.append(
            [
                row["name"],
                row["dg"],
                df_sub["dG_Mut2_GAAA"].mean(),
                df_sub["dGerr_Mut2_GAAA"].min(),
            ]
        )
    df_data = pd.DataFrame(data, columns=["name", "dg", "rna_map_dg", "rna_map_dg_err"])
    plt.errorbar(
        df_data["dg"],
        df_data["rna_map_dg"],
        yerr=df_data["rna_map_dg_err"],
        xerr=df_data["rna_map_dg_err"],
        fmt="o",
        capsize=5,
    )
    # df_data.to_csv("q_dms_ttr_paper/resources/csvs/rna_map_dg.csv", index=False)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
