"""
Functions to generate the constructs used in this study
"""

import pandas as pd
import vienna
from rna_lib_design.structure import rna_structure
from rna_secstruct_design.selection import get_selection, SecStruct
from rna_secstruct_design.helix_randomizer import HelixRandomizer


def split_dataframe(df, chunk_size=10000):
    """
    Splits a DataFrame into smaller chunks.

    This function takes a pandas DataFrame and a chunk size as input, and returns a list
    of smaller DataFrames, each with a maximum length of chunk_size. If the length of
    the DataFrame is not a multiple of chunk_size, the last chunk will be smaller.

    Args:
        df (pandas.DataFrame): The DataFrame to split.
        chunk_size (int): The maximum size of each chunk.

    Returns:
        list of pandas.DataFrame: A list of DataFrame chunks.

    Raises:
        ZeroDivisionError: If chunk_size is zero.

    Example:
    >>> df = pd.DataFrame(range(50))
    >>> chunks = split_dataframe(df, 10)
    >>> len(chunks)
    5
    >>> len(chunks[0])
    10
    """
    if chunk_size == 0:
        raise ZeroDivisionError("chunk_size cannot be zero")
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i * chunk_size : (i + 1) * chunk_size])
    return chunks


def indentify_scaffold(seq) -> str:
    """
    Find the scaffold used to compute the average dG from Bonilla et al.

    The function checks if the given RNA sequence contains any of the known scaffolds:
    S27, S28, S33, S36, and S37. If a match is found, the corresponding scaffold name is returned.
    If none of the known scaffolds are found, "OTHER" is returned.

    Args:
        seq (str): The sequence of the RNA.

    Returns:
        str: The scaffold name or "OTHER" if not found.

    Example:
    >>> indentify_scaffold("AACUGAGUCG")
    'S27'
    """
    if "AACUGAGUCG" in seq:
        return "S27"
    elif "AAUGCACAGG" in seq:
        return "S28"
    elif "AGGGAUCUUGGGAACAAGAUCCCU" in seq:
        return "S33"
    elif "AAGCCGGUCGGGAACGACCAGGCU" in seq:
        return "S36"
    elif "AAGCCGGUCGGGAACGACCGUGGC" in seq:
        return "S37"
    else:
        return "OTHER"


def get_full_tlr_seq_struct(seq, ss, r_seq):
    """
    Takes the sequence and secondary structure of the RNA and the TLR sequence
    to determine the full sequence and secondary structure of the TLR sequence. For
    some reason the TLR sequence does not have the flanking pairs and doesn't say what
    the secondary structure is.

    Args:
        seq (str): The sequence of the RNA.
        ss (str): The secondary structure of the RNA.
        r_seq (str): The TLR sequence.

    Returns:
        tuple: A tuple containing the full sequence and secondary structure of the TLR sequence.

    """
    ttr = r_seq.split("_")
    pos_1 = seq.find(ttr[0]) - 1
    pos_2 = seq.find(ttr[1])
    ttr_seq_1 = seq[pos_1 : pos_1 + len(ttr[0]) + 1]
    ttr_ss_1 = ss[pos_1 : pos_1 + len(ttr[0]) + 1]
    ttr_seq_2 = seq[pos_2 : pos_2 + len(ttr[1]) + 1]
    ttr_ss_2 = ss[pos_2 : pos_2 + len(ttr[1]) + 1]
    ttr_seq = ttr_seq_1 + "&" + ttr_seq_2
    ttr_ss = ttr_ss_1 + "&" + ttr_ss_2
    return ttr_seq, ttr_ss


def add_act_seq_and_ss(df):
    """
    Adds the full sequence and secondary structure of the TLR sequence to the dataframe.

    Args:
        df (pandas.DataFrame): The dataframe with the RNA sequences.

    Returns:
        pandas.DataFrame: The dataframe with the actual sequence and secondary structure of the TLR.

    """
    act_seq, act_ss, ac_counts = [], [], []
    for _, row in df.iterrows():
        ss = vienna.fold(row["seq"]).dot_bracket
        ttr_seq, ttr_ss = get_full_tlr_seq_struct(row["seq"], ss, row["name"])
        act_seq.append(ttr_seq)
        act_ss.append(ttr_ss)
        ac_count = 0
        for e, d in zip(ttr_seq, ttr_ss):
            if (e == "A" or e == "C") and d == ".":
                ac_count += 1
        ac_counts.append(ac_count)
    df["act_seq"] = act_seq
    df["act_ss"] = act_ss
    df["ac_count"] = ac_counts
    return df


def get_average_dg_dataframe():
    """
    takes the raw data from Bonilla et al. and computes the average dG for each. Thanks
    again for Steve Bonilla for sharing this raw data.

    Args:
        None

    Returns:
        pandas.DataFrame: The average dG for each TLR mutant.

    Example:
    >>> df = get_average_dg_dataframe()
    >>> df.head()
    """
    df = pd.read_csv(
        "data/construct_design/bonilla_et_al_pnas_2021/individual_measurements.csv"
    )
    # remove all that dont have measurements
    df = df[~df["dG_Mut2_GAAA"].isna()]
    # take only the values we need
    df = df[
        "index,b_name,loop_name,r_name,r_seq,seq,dG_Mut2_GAAA,dGerr_Mut2_GAAA".split(
            ","
        )
    ]
    # assign the scaffold name to each sequence
    df["scaffold"] = df["seq"].apply(indentify_scaffold)
    # remove all the sequences that are not in the 5 scaffolds
    df = df[df["scaffold"] != "OTHER"]
    # take only the value with normal flank sequences except S28 which doesnt seem to
    # have a normal sequence
    keep = []
    for i, row in df.iterrows():
        if row["scaffold"] == "S28" and row["b_name"] == "GC":
            keep.append(True)
        elif row["b_name"] == "normal":
            keep.append(True)
        else:
            keep.append(False)
    df = df[keep]
    groups = df.groupby("r_seq")
    data = []
    # take the average dG value and average dG error for each TLR mutant
    for _, g in groups:
        if len(g) == 0:
            continue
        row = g.iloc[0]
        avg_gaaa = g["dG_Mut2_GAAA"].mean()
        avg_gaaa_err = g["dGerr_Mut2_GAAA"].mean()
        data.append(
            [
                row["r_seq"],
                row["r_name"],
                row["seq"],
                avg_gaaa,
                avg_gaaa_err,
            ]
        )
    df = pd.DataFrame(
        data,
        columns=[
            "name",
            "name_common",  # name that was given to the RNA in the paper
            "seq",
            "dg",
            "dg_err",
        ],
    )
    return df


def flip_ss(db):
    """
    Flips the parentheses in a given string. Fixing a secondary structure when the strands
    of a motif are flipped.

    Args:
        db (str): The input string containing parentheses.

    Returns:
        str: The string with flipped parentheses.

    Example:
        >>> flip_ss(')))&(((')
        '(((&)))'
    """
    new_db = ""
    for s in db:
        if s == "(":
            new_db += ")"
        elif s == ")":
            new_db += "("
        else:
            new_db += s
    return new_db


def insert_motif_into_mttr6_scaffold(m_seq, m_ss):
    """
    Inserts a motif into the mttr6 scaffold sequence.

    Args:
        m_seq (str): The sequence of the motif.
        m_ss (str): The secondary structure of the motif.

    Returns:
        str: The new RNA structure with the motif inserted.

    """
    seq = "GAGCCUAUGGCUGCCACCCGAGCCCUUGAACUACAGGGAACACUGGAAACAGUACCCCCUGCAAGGGCGUUUGACGGUGGCAGCCUAAGGGCUC"
    ss = "((((((..((((((((((((((((((((.....(((((...((((....))))...))))))))))))..)))..))))))))))...))))))"
    struct = rna_structure(seq, ss)
    sub_1 = struct[:4]
    sub_2 = struct[10:-11]
    sub_3 = struct[-4:]
    seqs = m_seq.split("&")
    sss = m_ss.split("&")
    m_struct_1 = rna_structure(seqs[0], sss[0])
    m_struct_2 = rna_structure(seqs[1], sss[1])
    new_struct = sub_1 + m_struct_1 + sub_2 + m_struct_2 + sub_3
    return new_struct


def generate_mttr6_scaffold_library(df, path):
    """
    Generates a scaffold library for MTTR6 based on the given DataFrame.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data for generating the scaffold library.
        path (str): The path to the output file where the scaffold library will be written.

    Returns:
        None
    """
    f = open(path, "w")
    f.write("name,sequence,structure,dg\n")
    for i, row in df.iterrows():
        struct = insert_motif_into_mttr6_scaffold(row["act_seq"], row["act_ss"])
        folded_ss = vienna.fold(str(struct.sequence)).dot_bracket
        if struct.dot_bracket != folded_ss:
            print("failed")
        f.write(f"{row['name']},{struct.sequence},{struct.dot_bracket},{row['dg']}\n")
    f.close()


def generate_mttr6_randomized_helices_library(df, path):
    """
    Generates a randomized helices library based on the given dataframe and saves the results to a CSV file.

    Args:
        df (pandas.DataFrame): The input dataframe containing sequence and structure information.
        path (str): The file path to save the generated library.

    Returns:
        None
    """
    # these are the residues that should not change indenity
    params = {
        "motif_1": {"name": "gaaa_tetraloop"},
        "motif_two_way_junctions": {"m_type": "JUNCTION"},
        "seq_struct_kink_turn": {
            "sequence": "CCGAG&CGUUUGACG",
            "structure": "(((((&)..)))..)",
        },
        "range": "1-10,84-94",
    }
    hr = HelixRandomizer()
    data = []
    for row_i, row in df.iterrows():
        try:
            ss = SecStruct(row["sequence"], row["structure"])
        except:
            continue
        sele = get_selection(ss, params)
        results = hr.run(ss, exclude=sele)
        data.append([row["name"], results[1], row["structure"], results[0]])
    df_results = pd.DataFrame(
        data, columns=["name", "sequence", "structure", "ens_defect"]
    )
    df_results.to_csv(path, index=False)
