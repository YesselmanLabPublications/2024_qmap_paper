"""
attempts to refine and simplify mutations clasification of steves ttrs
"""
import pandas as pd
import numpy as np
import json


def get_raw_data():
    dropbox_path = "/Users/jyesselm/Dropbox"
    runs = [
        "2022_08_25_mtt6_set4_1st3_seq",
        "2022_08_26_mtt6_set4_2nd3_seq",
        "2022_08_29_mtt6_seq",
        "2022_09_01_mtt6_set4_seq",
        "2022_08_31_mtt6_set4_seq",
        "2022_08_26_mtt6_set1-3_MgTitra_KU_seq",
    ]
    path = f"{dropbox_path}/data/sequencing_analysis/summary/"
    dfs = []
    for run in runs:
        df = pd.read_json(f"{path}{run}.json")
        dfs.append(df)
    df = pd.concat(dfs)
    return df


def get_existing_mut_data():
    dropbox_path = "/Users/jyesselm/Dropbox"
    df_ref = pd.read_json(
        f"{dropbox_path}/projects/other/2022_01_27_steve_ttrs/data/all_sets_wt_muts.json"
    )
    df_ref.rename({"r_seq": "name", "dg_gaaa": "dg"}, axis=1, inplace=True)
    return df_ref


def assign_new_classifications(row, pos):
    d = {
        "name": row["name"],
        "pos": pos,
        "dg": row["dg"],
        "aligned_seq": row["aligned_seq"],
        "act_seq" : row["act_seq"],
        "act_ss": row["act_ss"],
        "insertions" : row["inserts"],
        "mutations" : row["mut_pos"],
        "deletions" : row["deletes"],
        "count_ins": len(row["inserts"]),
        "count_dels": len(row["deletes"]),
        "count_muts": len(row["muts"]),
        "has_bp_muts": 0,
        "bp_mut_1": None,
        "bp_mut_2": None,
        "bp_mut_3": None,
    }
    for pos in range(1, 12):
        d[f"mut_{pos}"] = None
        d[f"ins_{pos}"] = None
        d[f"dels_{pos}"] = None

    for mut in row["muts"]:
        new_seq = mut[0]
        pos = int(mut[1:-1])
        d[f"mut_{pos}"] = new_seq
    for insert in row["inserts"]:
        new_seq = insert[-1]
        pos = int(insert[:-1])
        d[f"ins_{pos}"] = new_seq
    # double check basepair mutations
    bp_pos = [[1, 11], [2, 10], [6, 7]]
    for i, bp in enumerate(bp_pos):
        if d[f"mut_{bp[0]}"] is not None or d[f"mut_{bp[1]}"] is not None:
            d[f"bp_mut_{i + 1}"] = (
                row["aligned_seq"][bp[0] - 1] + row["aligned_seq"][bp[1] - 1]
            )
    if (
        d["bp_mut_1"] is not None
        or d["bp_mut_2"] is not None
        or d["bp_mut_3"] is not None
    ):
        d["has_bp_muts"] = 1

    return d


def main():
    """
    main function for script
    """
    df_ref = get_existing_mut_data()
    ds = []
    for i, row in df_ref.iterrows():
        d = assign_new_classifications(row, i + 1)
        ds.append(d)
    df = pd.DataFrame(ds)
    df.to_json("classified_muts.json", orient="records")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
