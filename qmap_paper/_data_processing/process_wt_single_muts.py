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


def standard_mttr6_processing(df):
    df = to_rna(df)
    df = trim_p5_and_p3(df)
    # get motif data
    df["tlr"] = get_wt_tlr_data_in_dataframe(df)
    df["gaaa"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GGAAAC", "(....)")
    )
    df["gaaa_avg"] = df["gaaa"].apply(lambda x: np.mean(x[2:-1]))
    df["ref_hp_1"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), end=50
    )
    df["ref_hp_1_avg"] = df["ref_hp_1"].apply(lambda x: np.mean([x[2], x[5]]))
    df["ires"] = get_motif_data_in_dataframe(
        df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
    )
    df["kink_turn"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
    )
    df["3x3_motif"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GAACA&UACCC", "(...(&)...)")
    )
    return df


def process_no_tlr():
    runs = [
        "2023_03_13_no_tlr_Mg_titra_redo_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    df = to_rna(df)
    df = trim_p5_and_p3(df)
    df["gaaa"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GGAAAC", "(....)")
    )
    df["gaaa_avg"] = df["gaaa"].apply(lambda x: np.mean(x[2:-1]))
    df["ref_hp_1"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), end=50
    )
    df["ref_hp_1_avg"] = df["ref_hp_1"].apply(lambda x: np.mean([x[2], x[5]]))
    df["ires"] = get_motif_data_in_dataframe(
        df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
    )
    df["kink_turn"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
    )
    df["3x3_motif"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GAACA&UACCC", "(...(&)...)")
    )
    df.to_json("data/no_tlr_mg.json", orient="records")


def process_uucg():
    runs = [
        "2023_02_02_minittr_6_uucg_Mg_titra_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    df = to_rna(df)
    df = trim_p5_and_p3(df)
    df["tlr"] = get_wt_tlr_data_in_dataframe(df)
    df["ref_hp_1"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), end=50
    )
    df["ref_hp_1_avg"] = df["ref_hp_1"].apply(lambda x: np.mean([x[2], x[5]]))
    df["ires"] = get_motif_data_in_dataframe(
        df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
    )
    df["kink_turn"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
    )
    df["3x3_motif"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GAACA&UACCC", "(...(&)...)")
    )
    df.to_json("data/uucg_mg_titra.json", orient="records")


def process_no_3x3():
    runs = [
        "2023_02_17_no_3_3_junc_Mg_titr_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    df = to_rna(df)
    df = trim_p5_and_p3(df)
    df["tlr"] = get_wt_tlr_data_in_dataframe(df)
    df["gaaa"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("GGAAAC", "(....)")
    )
    df["gaaa_avg"] = df["gaaa"].apply(lambda x: np.mean(x[2:-1]))
    df["ref_hp_1"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CGAGUAG", "(.....)"), end=50
    )
    df["ref_hp_1_avg"] = df["ref_hp_1"].apply(lambda x: np.mean([x[2], x[5]]))
    df["ires"] = get_motif_data_in_dataframe(
        df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
    )
    df["kink_turn"] = get_sec_struct_data_in_dataframe(
        df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
    )
    df.to_json("data/no_3x3_mg_titra.json", orient="records")


def process_h1_3bp_longer():
    runs = [
        "2023_03_10_h1_3bp_longer_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    df = df[df["name"] == "minittr_6_h1_3bp_longer_fixed"]
    df = standard_mttr6_processing(df)
    df.to_json("data/h1_3bp_longer_mg_titra.json", orient="records")


def process_h2_3bp_longer():
    runs = [
        "2023_03_15_h2_3bp_longer_Mg_titra_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    # df = df[df["name"] == "minittr_6_h1_3bp_longer_fixed"]
    df = standard_mttr6_processing(df)
    df.to_json("data/h2_3bp_longer_mg_titra.json", orient="records")


def process_h3_3bp_longer():
    runs = [
        "2023_03_22_h3_3bp_longer_Mg_titra_seq",
    ]
    path = Path("/Users/jyesselman2/Dropbox/data/sequencing_analysis")
    df = get_preprocessed_data(path, runs)
    # df = df[df["name"] == "minittr_6_h1_3bp_longer_fixed"]
    df = standard_mttr6_processing(df)
    df.to_json("data/h3_3bp_longer_mg_titra.json", orient="records")


def main():
    """
    main function for script
    """
    process_h3_3bp_longer()


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
