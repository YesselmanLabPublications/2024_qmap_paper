from pathlib import Path
import pandas as pd
import numpy as np

from seq_tools import has_5p_sequence, to_rna
from seq_tools import trim as seq_ss_trim
from seq_tools import SequenceStructure
from seq_tools.structure import find as seq_ss_find
from rna_secstruct import SecStruct, MotifSearchParams


def trim(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """
    trims the dataframe to the given start and end
    :param df: a dataframe with data
    :param start: the start index
    :param end: the end index
    :return: a trimmed dataframe
    """
    df = seq_ss_trim(df, start, end)
    df["data"] = df["data"].apply(lambda x: x[start:-end])
    return df


def trim_p5_and_p3(df):
    """
    trims the 5' and 3' ends of the data
    :param df: a dataframe with data
    :return: a trimmed dataframe
    """
    df_p5 = pd.read_csv(f"resources/p5_sequences.csv")
    df_p5 = to_rna(df_p5)
    common_p5_seq = ""
    for p5_seq in df_p5["sequence"]:
        if has_5p_sequence(df, p5_seq):
            common_p5_seq = p5_seq
    if len(common_p5_seq) == 0:
        raise ValueError("No common p5 sequence found")
    print(f"common p5 sequence: {common_p5_seq}")
    return trim(df, len(common_p5_seq), 20)
