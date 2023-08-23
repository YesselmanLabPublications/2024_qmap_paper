from pathlib import Path
import numpy as np

from seq_tools import to_rna
from seq_tools.structure import SequenceStructure
from rna_secstruct import SecStruct, MotifSearchParams

from get_data import get_preprocessed_data
from process import trim_p5_and_p3
from segment_data import (
    get_sec_struct_data_in_dataframe,
    get_motif_data_in_dataframe,
    get_wt_tlr_data_in_dataframe,
)


def print_out_duplicates(df):
    for i, g in df.groupby(["name", "mg_conc", "buffer_conc"]):
        unique_runs = g["run_name"].unique()
        unique_rna = g["rna_name"].unique()
        if len(unique_runs) == 1:
            continue
        print(i, unique_runs)
        print(unique_rna)


def clean_up_dataset(df):
    # ensure data is rna
    df = to_rna(df)
    # fix some errors in the data
    # remove constructs not related to this analysis
    include = ["minittr-6-2HP-ref"]
    df = df.query("name  in @include")
    # round data
    # for the data column round each value to 5 decimal places
    # df["data"] = df["data"].apply(lambda x: [round(y, 5) for y in x])
    return df


def main():
    """
    main function for script
    """
    runs = [
        "2022_07_26_minittr-6-2HP-ref_buffer_seq",
        "2022_07_20_minittr_Hepes-titra_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    df = clean_up_dataset(df)
    df = trim_p5_and_p3(df)
    # get motif data
    df["gaaa"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GGAAAC", "(....)")
    )
    df["gaaa_avg"] = df["gaaa"].apply(lambda x: np.mean(x[2:-1]))
    df["tlr"] = get_wt_tlr_data_in_dataframe(df)
    df["ref_hp_1"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), end=50
    )
    df["ref_hp_1_avg"] = df["ref_hp_1"].apply(lambda x: np.mean([x[2], x[5]]))
    df["ref_hp_2"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), start=50
    )
    df["ref_hp_2_avg"] = df["ref_hp_2"].apply(lambda x: np.mean([x[2], x[5]]))
    df["ires"] = get_motif_data_in_dataframe(
        df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
    )
    df["kink_turn"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
    )
    df["3x3_motif"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GAACA&UACCC", "(...(&)...)")
    )
    df.to_json("data/wt_buffer_titra.json", orient="records")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
