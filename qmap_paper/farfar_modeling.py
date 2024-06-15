import os
from biopandas.pdb import PandasPdb

from qmap_paper.logger import get_logger

script_name = os.path.splitext(os.path.basename(__file__))[0]
log = get_logger(script_name)

fasta_txt = """
>loop
ggaaac

>receptor_1
g{seq1}

>receptor_2
{seq2}c
"""


def get_dataframe_from_pdb(pdb_path):
    """
    Reads a PDB file and returns the 'ATOM' records as a pandas DataFrame.

    Args:
        pdb_path (str): The path to the PDB file.

    Returns:
        pandas.DataFrame: The 'ATOM' records from the PDB file.

    """
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_path)
    return ppdb.df["ATOM"]


def write_pdb(df, pdb_path):
    """
    Writes the given DataFrame to a PDB file.

    Args:
        df (pandas.DataFrame): The DataFrame containing the ATOM records.
        pdb_path (str): The path to the output PDB file.

    Returns:
        None
    """
    ppdb = PandasPdb()
    ppdb.df["ATOM"] = df
    ppdb.to_pdb(path=pdb_path, records=["ATOM"], gz=False, append_newline=True)


def setup_farfar_modeling_runs(df):
    """
    Sets up the FARFAR modeling runs based on the given dataframe.

    Args:
        df (pandas.DataFrame): The dataframe containing the necessary data for setting up the modeling runs.

    Returns:
        None
    """
    log.info("Setting up FARFAR modeling runs...")
    output_path = "farfar_runs/"
    log.info(f"Output path: {output_path}")
    os.makedirs(output_path, exist_ok=True)
    pdb_path = "data/pdbs/farfar2_modeling/"
    for i, row in df.iterrows():
        df_pdb = get_dataframe_from_pdb(f"{pdb_path}/start.pdb")
        new_name = row["act_seq"].replace("&", "_")
        dir_name = f"{output_path}/{new_name}"
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
        # os.chdir(dir_name)
        # os.system(
        #    "rna_denovo.macosclangrelease -fasta test.fasta -secstruct_file "
        #    "test.secstruct_file -s start.pdb -minimize_rna true -nstruct 1"
        # )
        # os.chdir("../..")
