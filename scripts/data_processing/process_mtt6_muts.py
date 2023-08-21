from pathlib import Path
import pandas as pd
import yaml
import numpy as np

from seq_tools import has_5p_sequence, to_rna
from seq_tools import trim as seq_ss_trim
from seq_tools import SequenceStructure
from seq_tools.structure import find as seq_ss_find
from rna_secstruct import SecStruct, MotifSearchParams

from get_data import get_preprocessed_data
from process import trim_p5_and_p3

# TODO use the rna_secstruct_design code to do the selections maybe even move that into
# its own package
# TODO also look at the recent code for parsing that I wrote for the mtt6 data


def get_motif_data_in_dataframe(df, params):
    params = MotifSearchParams(**params)
    all_data = []
    for i, row in df.iterrows():
        ss = SecStruct(row["sequence"], row["structure"])
        motifs = ss.get_motifs(params)
        if len(motifs) != 1:
            raise ValueError("More than one motif found")
        data = []
        for s in motifs[0].strands:
            for e in s:
                data.append(row["data"][e])
        all_data.append(data)
    return all_data


def get_sec_struct_data_in_dataframe(df, sub_seq_struct, start=None, end=None):
    all_data = []
    for i, row in df.iterrows():
        ss = SequenceStructure(row["sequence"], row["structure"])
        r = seq_ss_find(ss, sub_seq_struct, start, end)
        if len(r) != 1:
            raise ValueError("More than one segment found")
        pos = []
        bounds = r[0]
        for r in bounds:
            pos.extend(list(range(r[0], r[1])))
        all_data.append([row["data"][p] for p in pos])
    return all_data


def print_out_duplicates(df):
    for i, g in df.groupby(["name", "mg_conc"]):
        unique_runs = g["run_name"].unique()
        unique_rna = g["rna_name"].unique()
        if len(unique_runs) == 1:
            continue
        print(i, unique_runs)
        print(unique_rna)


def get_existing_mut_data():
    dropbox_path = "/Users/jyesselman2/Dropbox"
    df_ref = pd.read_json(
        f"{dropbox_path}/projects/other/2022_01_27_steve_ttrs/data/all_sets_wt_muts.json"
    )
    df_ref.rename({"r_seq": "name", "dg_gaaa": "dg"}, axis=1, inplace=True)
    return df_ref


def clean_up_dataset(df):
    # ensure data is rna
    df = to_rna(df)
    # fix some errors in the data
    # 0.20 was set instead of 0.25 for some constructs
    df.loc[df["mg_conc"] == 0.20, "mg_conc"] = 0.25
    # remove constructs not related to this analysis
    exclude = ["PurRe4L5", "PtrRe2L5", "PurRe2L4"]
    df = df.query("name not in @exclude")
    # duplicate data not sure right now which is better
    q1 = df["run_name"] == "2022_08_29_mtt6_seq"
    q2 = df["rna_name"] == "mtt6_mutations_set_2_50mM-NaC_0.1-mM-Mg2+"
    df = df[~(q1 & q2)]
    q1 = df["run_name"] == "2022_08_29_mtt6_seq"
    q3 = df["rna_name"] == "mtt6_mutations_set_1_50mM-NaC_5-mM-Mg2+"
    df = df[~(q1 & q3)]
    # check there are no duplicates
    print_out_duplicates(df)
    # round data
    # for the data column round each value to 5 decimal places
    df["data"] = df["data"].apply(lambda x: [round(y, 5) for y in x])
    return df


def assign_tlr_reactivity(row, data):
    # print(data)
    # print(len(data))
    # print(json.dumps(row.to_dict(), indent=4))
    # print(row['aligned_seq'])
    # print(len(row['aligned_seq']))
    # seq = row["act_seq"].split("&")[1][:-1] + row["act_seq"].split("&")[0][1:]
    # print(len(seq))
    new_cols = {}
    insert_pos = []
    for insert in row["inserts"]:
        pos = int(insert[:-1])
        insert_pos.append(pos)
    mapped_pos = 1
    pos = 0
    new_cols["tlr"] = []
    new_cols["tlr_norm"] = []
    while pos < len(data):
        if mapped_pos in insert_pos:
            new_cols["tlr"].append(data[pos])
            new_cols["tlr_norm"].append(data[pos] / row["ref_hp"])
            new_cols[f"tlr_in_{mapped_pos}"] = data[pos + 1] / row["ref_hp"]
            pos += 1
        else:
            new_cols["tlr"].append(data[pos])
            new_cols["tlr_norm"].append(data[pos] / row["ref_hp"])
        mapped_pos += 1
        pos += 1
    return new_cols


def get_tlr_reactivities(df) -> pd.DataFrame:
    all_data = []
    for i, row in df.iterrows():
        full_ss = SequenceStructure(row["sequence"], row["structure"])
        ss = SequenceStructure(row["act_seq"], row["act_ss"])
        bounds = seq_ss_find(full_ss, ss)[0]
        # trims off last basepair of 5' and first basepair of 3'
        # flips the orientation to match aligend sequence and structure
        data = (
            row["data"][bounds[1][0] : bounds[1][1] - 1]
            + row["data"][bounds[0][0] + 1 : bounds[0][1]]
        )
        data = assign_tlr_reactivity(row, data)

        all_data.append(data)
    df_m = pd.DataFrame(all_data)
    df = pd.concat([df, df_m], axis=1)
    return df


def main():
    """
    main function for script
    """
    runs = [
        "2022_08_25_mtt6_set4_1st3_seq",
        "2022_08_26_mtt6_set4_2nd3_seq",
        "2022_08_26_mtt6_set1-3_MgTitra_KU_seq",
        "2022_08_29_mtt6_seq",
        "2022_08_30_mtt6_set4_seq",
        "2022_08_31_mtt6_set4_seq",
        "2022_09_01_mtt6_set4_seq",
    ]
    # data to remove
    # CAUGA_UCUAAA <- messed up tetraloop
    # UACGG_CCUACA <- messed up tetraloop
    # UACGG_CCUACA <- odd tetraloop reactivity

    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    df = clean_up_dataset(df)
    # trim off extra sequences from p5 and p3 and also for the data
    df = trim_p5_and_p3(df)
    df_ref = pd.read_csv("data/rna_map/chemical_mapped.csv")
    df_ref = df_ref[["name", "dg", "act_seq", "act_ss"]]
    df = df.merge(df_ref, on="name")
    # add in motif data
    df["gaaa"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GGAAAC", "(....)")
    )
    df["gaaa_avg"] = df["gaaa"].apply(lambda x: np.mean(x[2:-1]))
    df["ref_hp"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), end=50
    )
    df["ref_hp"] = df["ref_hp"].apply(lambda x: np.mean([x[2], x[5]]))
    df_muts = pd.read_json("data/mttr6_mut_charactization.json")
    df_muts.rename({"r_seq": "name"}, axis=1, inplace=True)
    df_muts = df_muts[["name", "aligned_seq", "inserts", "deletes", "mut_pos", "muts"]]
    df = df.merge(df_muts, on="name")
    df_mg = pd.read_csv("data/mg_1_2.csv")
    df_mg.rename({"k": "mg_1_2"}, axis=1, inplace=True)
    df = df.merge(df_mg, on="name")
    df = get_tlr_reactivities(df)
    df.to_json("data/mttr6_data_full.json", orient="records")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
