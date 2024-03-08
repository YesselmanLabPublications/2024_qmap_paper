import pandas as pd
import numpy as np
import re
import warnings
from Bio import pairwise2
import click


# helper functions ####################################################################


def find_chars_in_str(s, ch):
    for i, ltr in enumerate(s):
        if ltr == ch:
            yield i


def remove_char_from_str(s, i):
    return s[i], s[:i] + s[i + 1 :]


def is_bp(bp_name):
    names = "AU,UA,GC,CG,GU,UG".split(",")
    if bp_name in names:
        return True
    else:
        return False


def is_bp_mutation(wt, mut, p1, p2):
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


def get_insertions(wt_seqs, r_seq):
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


def get_deletions(wt_seqs, r_seq):
    spl = r_seq.split("_")
    a1 = pairwise2.align.globalms(spl[0], wt_seqs[0], 2, -1, -10, -10)
    a2 = pairwise2.align.globalms(spl[1], wt_seqs[1], 2, -1, -10, -10)
    mut_seq_rev_aligned = a1[0].seqA + a2[0].seqA
    wt_seq_rev_aligned = a1[0].seqB + a2[0].seqB
    deletes = list(find_chars_in_str(mut_seq_rev_aligned, "-"))[::-1]
    seq_deletes = []
    new_str = wt_seq_rev_aligned
    for delete in deletes:
        seq_deletes.append(f"{delete + 1}{wt_seq_rev_aligned[delete]}")
    return seq_deletes, mut_seq_rev_aligned


def get_diff_muts(list1, list2):
    diff = set(list1).symmetric_difference(list2)
    diff_nums = []
    for d in diff:
        diff_nums.append(int(d[:-1]))
    diff_nums = set(diff_nums)
    return len(diff_nums)


def get_diff(row1, row2):
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


def get_aligned_sequence_and_ins_del(df):
    """
    for each sequence in the dataframe, get the aligned sequence and the insertions
    and deletions. Right now this is still naive as it only considers that larger
    sequences have insertions and smaller sequences have deletions. It is possible
    that sequences can contain both insertions and deletions.
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
            print(row["name"])
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


def characterize_mutations(df):
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

    # print(len(df), count)
