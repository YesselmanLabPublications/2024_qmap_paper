import pandas as pd
import os
from biopandas.pdb import PandasPdb
import shutil

fasta_txt = """
>loop
ggaaac

>receptor_1
g{seq1}

>receptor_2
{seq2}c
"""


def get_dataframe_from_pdb(pdb_path):
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_path)
    return ppdb.df["ATOM"]


def write_pdb(df, pdb_path):
    ppdb = PandasPdb()
    ppdb.df["ATOM"] = df
    ppdb.to_pdb(path=pdb_path, records=["ATOM"], gz=False, append_newline=True)


def setup_farfar_modeling_runs(df):
    output_path = "test_output/farfar_modeling"
    os.makedirs(output_path, exist_ok=True)
    for i, row in df.iterrows():
        df_pdb = get_dataframe_from_pdb("start.pdb")
        new_name = row["act_seq"].replace("&", "_")
        dir_name = f"runs/{new_name}"
        os.makedirs(dir_name, exist_ok=True)
        seqs = row["act_seq"].lower().split("&")
        fasta_str = fasta_txt.format(seq1=seqs[0], seq2=seqs[1])
        f = open(f"{dir_name}/test.fasta", "w")
        f.write(fasta_str)
        f.close()
        f = open(f"{dir_name}/test.secstruct_file", "w")
        f.write(f"(....)+({row['act_ss'].replace('&', '+')})\n")
        f.write(f"ggaaac+g{seqs[0]}+{seqs[1]}c\n")
        f.close()
        pos = 6
        res_1_num = pos + 1
        res_2_num = len(seqs[0]) + len(seqs[1]) + pos + 2
        df_pdb.loc[df_pdb["residue_number"] == 21, "residue_number"] = res_2_num
        write_pdb(df_pdb, f"{dir_name}/start.pdb")
        os.chdir(dir_name)
        os.system(
            "rna_denovo.macosclangrelease -fasta test.fasta -secstruct_file "
            "test.secstruct_file -s start.pdb -minimize_rna true -nstruct 1"
        )
        os.chdir("../..")
