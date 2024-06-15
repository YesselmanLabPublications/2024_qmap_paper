import pandas as pd
from typing import List, Generator, Union, Tuple, Any
from Bio import pairwise2


# helper functions ####################################################################


def find_chars_in_str(s: str, ch: str) -> Generator[int, None, None]:
    """Finds the indices of a given character in a string.

    Args:
        s (str): The input string.
        ch (str): The character to search for.

    Yields:
        int: The indices of the character in the string.

    """
    for i, ltr in enumerate(s):
        if ltr == ch:
            yield i


def remove_char_from_str(s: str, i: int) -> tuple:
    """
    Removes a character from a string at the specified index.

    Args:
        s (str): The input string.
        i (int): The index of the character to remove.

    Returns:
        tuple: A tuple containing the removed character and the modified string.

    Example:
        >>> remove_char_from_str('hello', 2)
        ('l', 'helo')
    """
    return s[i], s[:i] + s[i + 1 :]


def is_bp(bp_name: str) -> bool:
    """
    Checks if a given base pair name is valid.

    Args:
        bp_name (str): The name of the base pair to check.

    Returns:
        bool: True if the base pair name is valid, False otherwise.
    """
    names = "AU,UA,GC,CG,GU,UG".split(",")
    if bp_name in names:
        return True
    else:
        return False


def is_bp_mutation(wt: str, mut: str, p1: int, p2: int) -> Union[List[str], None]:
    """
    Determines if a given mutation is a base pair mutation.

    Args:
        wt: The wildtype sequence.
        mut: The mutated sequence.
        p1: The index of the first position to compare.
        p2: The index of the second position to compare.

    Returns:
        If the mutation is a base pair mutation, returns a list containing the original base pair name and the mutated base pair name. Otherwise, returns None.
    """
    if wt[p1] == mut[p1]:
        return None
    if wt[p2] == mut[p2]:
        return None
    bp_name = mut[p1] + mut[p2]
    org_bp_name = wt[p1] + wt[p2]
    if is_bp(bp_name):
        return [org_bp_name, bp_name]
    else:
        return None


# main functions ####################################################################


def get_insertions(wt_seqs: List[str], r_seq: str) -> Tuple[List[str], str]:
    """
    Get the insertions in the mutated sequence compared to the wild-type sequences.

    Args:
        wt_seqs (list): List of wild-type sequences.
        r_seq (str): Mutated sequence.

    Returns:
        tuple: A tuple containing:
            - seq_inserts (list): List of insertions in the mutated sequence.
            - new_str (str): The mutated sequence with insertions removed.
    """
    spl = r_seq.split("_")
    a1 = pairwise2.align.globalms(wt_seqs[0], spl[0], 2, -1, -10, -10)
    a2 = pairwise2.align.globalms(wt_seqs[1], spl[1], 2, -1, -10, -10)
    wt_seq_rev_aligned = a1[0].seqA + a2[0].seqA
    mut_seq_rev_aligned = a1[0].seqB + a2[0].seqB
    inserts = list(find_chars_in_str(wt_seq_rev_aligned, "-"))[::-1]
    seq_inserts = []
    new_str = mut_seq_rev_aligned
    for ins in inserts:
        substr = wt_seq_rev_aligned[:ins]
        char, new_str = remove_char_from_str(new_str, ins)
        seq_inserts.append(f"{ins - substr.count('-') + 1}{char}")
    return seq_inserts, new_str


def get_deletions(wt_seqs: List[str], r_seq: str) -> Tuple[List[str], str]:
    """
    Get deletions in the mutated sequence compared to the wild-type sequences.

    Args:
        wt_seqs (List[str]): List of wild-type sequences.
        r_seq (str): Mutated sequence.

    Returns:
        Tuple[List[str], str]: A tuple containing the list of deletions in the mutated sequence
        and the mutated sequence aligned to the wild-type sequence.
    """
    spl = r_seq.split("_")
    a1 = pairwise2.align.globalms(spl[0], wt_seqs[0], 2, -1, -10, -10)
    a2 = pairwise2.align.globalms(spl[1], wt_seqs[1], 2, -1, -10, -10)
    mut_seq_rev_aligned = a1[0].seqA + a2[0].seqA
    wt_seq_rev_aligned = a1[0].seqB + a2[0].seqB
    deletes = list(find_chars_in_str(mut_seq_rev_aligned, "-"))[::-1]
    seq_deletes = []
    for delete in deletes:
        seq_deletes.append(f"{delete + 1}{wt_seq_rev_aligned[delete]}")
    return seq_deletes, mut_seq_rev_aligned


def get_diff_muts(list1: List[str], list2: List[str]) -> int:
    """
    Calculates the number of different mutations between two lists.

    Args:
        list1 (List[str]): The first list of mutations.
        list2 (List[str]): The second list of mutations.

    Returns:
        int: The number of different mutations between the two lists.
    """
    diff = set(list1).symmetric_difference(list2)
    diff_nums = []
    for d in diff:
        diff_nums.append(int(d[:-1]))
    diff_nums = set(diff_nums)
    return len(diff_nums)


def get_diff(row1: List[Any], row2: List[Any]) -> int:
    """
    Calculates the number of differences between two rows.

    Args:
        row1 (List[Any]): The first row.
        row2 (List[Any]): The second row.

    Returns:
        int: The number of differences between the two rows.
    """
    seq_1 = row1[-6]
    seq_2 = row2[-6]
    count = 0
    for i, (e1, e2) in enumerate(zip(seq_1, seq_2)):
        if e1 == e2:
            continue
        count += 1
    # add insertions
    count += get_diff_muts(row1[-5], row2[-5])
    # add deletions
    count += get_diff_muts(row1[-4], row2[-4])
    return count


def get_aligned_sequence_and_ins_del(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each sequence in the dataframe, get the aligned sequence and the insertions
    and deletions.

    Args:
        df (pandas.DataFrame): The input dataframe containing the sequences.

    Returns:
        pandas.DataFrame: The modified dataframe with additional columns for aligned sequence,
        insertions, deletions, number of insertions, number of deletions, and size.

    Notes:
        This function assumes that larger sequences have insertions and smaller sequences have deletions.
        It is possible that sequences can contain both insertions and deletions.
    """
    # wild type information
    wt_seq = "CCUAAG_UAUGG"
    wt_seqs = ["CCUAAG", "UAUGG"]
    df["aligned_seq"] = ["" for _ in range(len(df))]
    df["insertions"] = [[] for _ in range(len(df))]
    df["deletions"] = [[] for _ in range(len(df))]
    df["num_ins"] = [0 for _ in range(len(df))]
    df["num_dels"] = [0 for _ in range(len(df))]
    df["size"] = [len(s) for s in df["r_seq"]]
    for ii, row in df.iterrows():
        if len(row["name"]) > len(wt_seq):
            inserts, aligned_seq = get_insertions(wt_seqs, row["r_seq"])
            df.at[ii, "insertions"] = inserts
            df.at[ii, "aligned_seq"] = aligned_seq
            df.at[ii, "num_ins"] = len(inserts)
        elif len(row["name"]) == len(wt_seq):
            spl = row["r_seq"].split("_")
            df.at[ii, "aligned_seq"] = spl[0] + spl[1]
        else:
            deletes, aligned_seq = get_deletions(wt_seqs, row["name"])
            df.at[ii, "deletions"] = deletes
            df.at[ii, "aligned_seq"] = aligned_seq
            df.at[ii, "num_dels"] = len(deletes)
    return df


def characterize_mutations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Characterizes mutations in a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing mutation data.

    Returns:
        pd.DataFrame: The DataFrame with additional columns characterizing the mutations.

    Raises:
        None
    """
    df["r_seq"] = df["name"]
    # bonilla scheme for numbering
    wt_seq_rev = "CCUAAGUAUGG"
    df = get_aligned_sequence_and_ins_del(df)
    df = df[df["aligned_seq"] != ""]
    # remove constructs where I couldnt align
    df_sub = df[df["aligned_seq"] != ""]
    count = 0
    all_bp_pos = [[0, -1], [1, -2], [5, 6]]
    data = []
    for _, row in df_sub.iterrows():
        seq_rev = row["aligned_seq"]
        bp_mut_pos = []
        bp_muts = []
        for i, bp_pos in enumerate(all_bp_pos):
            r = is_bp_mutation(wt_seq_rev, seq_rev, *bp_pos)
            if r is None:
                continue
            bp_muts.append(f"{r[0]}{i + 1}{r[1]}")
            bp_mut_pos.extend(bp_pos)
        mut_pos = []
        muts = []
        for i, (wt, new) in enumerate(zip(seq_rev, wt_seq_rev)):
            if wt == new:
                continue
            mut = f"{wt}{i + 1}{new}"
            muts.append(mut)
            mut_pos.append(i + 1)
        data.append([row["name"], bp_muts, mut_pos, muts])
        count += 1
    df_muts = pd.DataFrame(data, columns="name,bp_muts,mut_pos,mutations".split(","))
    df_sub = df_sub.merge(df_muts, on="name")
    return df_sub
