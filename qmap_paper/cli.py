import click
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import warnings

# Suppress BiopythonDeprecationWarning
warnings.filterwarnings("ignore", message="Bio.pairwise2 has been deprecated")

from qmap_paper.logger import setup_applevel_logger, get_logger
from qmap_paper.paths import RESOURCES_PATH, DATA_PATH
from qmap_paper.data_processing import (
    MTTR6BufferTitrationDataProcessor,
    MTTR6MgTitrationDataProcessor,
    TTRMutsDataProcessor,
    MTTR6MutsDataProcessor,
    compute_all_mg_1_2,
    compute_mg_1_2,
    trim,
)
from qmap_paper import mutation_characterize
from qmap_paper.plotting import (
    plot_pop_avg_titration,
    plot_mg_titration_fit,
    plot_pop_avg_from_row,
)
from qmap_paper.construct_design import (
    get_average_dg_dataframe,
    split_dataframe,
    add_act_seq_and_ss,
    generate_mttr6_scaffold_library,
    generate_mttr6_randomized_helices_library,
)
from qmap_paper.sasa import compute_solvent_accessability
from qmap_paper.farfar_modeling import setup_farfar_modeling_runs


log = get_logger("CLI")


def get_motif_data(df_wt, df_uucg, motif_name, motif_seq, motif_ss):
    """
    Retrieves motif data from the given dataframes.

    Args:
        df_wt (pandas.DataFrame): DataFrame containing wild-type data.
        df_uucg (pandas.DataFrame): DataFrame containing UUCG data.
        motif_name (str): Name of the motif.
        motif_seq (list): List of motif sequences.
        motif_ss (list): List of motif secondary structures.

    Returns:
        pandas.DataFrame: DataFrame containing the retrieved motif data.
    """
    data = []
    for i, row in df_wt.iterrows():
        count = 0
        uucg_row = df_uucg[df_uucg["mg_conc"] == row["mg_conc"]].iloc[0]
        for seq, ss, wt_val, uucg_val in zip(
            motif_seq, motif_ss, row[motif_name], uucg_row[motif_name]
        ):
            count += 1
            if seq != "A" and seq != "C":
                continue
            data.append(
                {
                    "name": str(count) + seq,
                    "seq": seq,
                    "ss": ss,
                    "wt_val": wt_val,
                    "wt_norm": wt_val / row["ref_hp_1_as"],
                    "uucg_val": uucg_val,
                    "uucg_norm": uucg_val / uucg_row["ref_hp_1_as"],
                    "mg_conc": row["mg_conc"],
                }
            )
    df = pd.DataFrame(data)
    return df


# cli commands ####################################################################


@click.group()
def cli():
    """cli for qmap_paper"""
    pass


# construct generation ################################################################
@cli.command()
def generate_sets():
    """
    Selects the TLR mutants from the Bonilla et. al paper to be used in this study.

    This function reads a CSV file containing average dG values and filters the data based on certain criteria.
    The filtered data is then split into groups and saved as separate CSV files.

    Returns:
        None
    """
    setup_applevel_logger()
    df = get_average_dg_dataframe()
    # this is the table from Bonilla et. al
    # only contains class A which behave like the wild-type
    path = "data/construct_design/bonilla_et_al_pnas_2021/table_s3_dG_avgs.csv"
    df_summary = pd.read_csv(path)
    df_summary["name"] = [r["seq2"] + "_" + r["seq1"] for _, r in df_summary.iterrows()]
    df = df.merge(df_summary, on="name")
    df["dg"] = df["dg_gaaa"]
    df.drop(
        columns=[
            "dg_gaaa",
            "dg_gaua",
            "dg_diff",
            "natural",
            "seq1",
            "seq2",
            "type",
            "class",
            "subgroup",
        ],
        inplace=True,
    )
    df = add_act_seq_and_ss(df)
    df = df.sort_values("dg")
    df = df[df["ac_count"] > 3]
    df = df[df["act_ss"].str.contains("((&))", regex=False)]
    # fix name to be consistent with paper
    for i, row in df.iterrows():
        df.at[i, "name"] = row["name"].split("_")[1] + "_" + row["name"].split("_")[0]
    df_groups = split_dataframe(df, 25)
    dfs = []
    output_path = "data/construct_design/sets"
    os.makedirs(output_path, exist_ok=True)
    for i, df in enumerate(df_groups):
        df.to_csv(f"{output_path}/set_{i + 1}.csv", index=False)
        df = pd.DataFrame(df)
        df["set"] = i + 1
        dfs.append(df)
    df = pd.concat(dfs)
    df.to_csv(f"{output_path}/all_sets.csv", index=False)


@cli.command()
def generate_scaffold_sets():
    """
    Generates scaffold sets based on input CSV files.

    This function reads CSV files from the 'data/construct_design/sets' directory,
    generates scaffold sets, and saves them as CSV files in the 'data/construct_design/scaffold_sets' directory.

    Args:
        None

    Returns:
        None
    """
    setup_applevel_logger()
    input_path = "data/construct_design/sets"
    output_path = "data/construct_design/scaffold_sets"
    os.makedirs(output_path, exist_ok=True)
    for i in range(1, 5):
        df = pd.read_csv(f"{input_path}/set_{i}.csv")
        generate_mttr6_scaffold_library(df, f"{output_path}/mtt6_mutations_set_{i}.csv")


@cli.command()
def randomize_helices():
    """
    Randomizes helices in the input CSV files and saves the randomized sets to the output directory.

    Args:
        None

    Returns:
        None
    """
    setup_applevel_logger()
    input_path = "data/construct_design/scaffold_sets"
    output_path = "data/construct_design/randomized_sets"
    os.makedirs(output_path, exist_ok=True)
    for i in range(1, 5):
        df = pd.read_csv(f"{input_path}/mtt6_mutations_set_{i}.csv")
        generate_mttr6_randomized_helices_library(
            df, f"{output_path}/mtt6_mutations_set_{i}.csv"
        )


@cli.command()
def barcode_libraries():
    """
    Generate barcode libraries for mutated sets.

    This function generates barcode libraries for mutated sets by executing the 'rld barcode2' command
    with specific parameters. The input files are located in the 'data/construct_design/randomized_sets'
    directory, and the output files will be saved in the 'data/construct_design/final_sets' directory.

    The function creates the output directory if it doesn't exist and then iterates over a range of values
    from 1 to 5. For each value, it constructs a command string using f-strings and executes it using the
    'os.system' function.

    Note: This function assumes that the 'rld' command and the required input files are available in the
    system's PATH. rld is from https://github.com/jyesselm/rna_lib_design

    Args:
        None

    Returns:
        None
    """
    setup_applevel_logger()
    input_path = "data/construct_design/randomized_sets"
    output_path = "data/construct_design/final_sets"
    os.makedirs(output_path, exist_ok=True)
    for i in range(1, 5):
        cmd = (
            f"rld barcode2 --trim-p5 4 --trim-p3 4 -o {output_path}/mtt6_mutations_set_{i} "
            f"{input_path}/mtt6_mutations_set_{i}.csv --param-file resources/barcode.yml"
        )
        os.system(cmd)


@cli.command()
def generate_order():
    """
    Generates an order by merging and processing data from multiple CSV files.

    This function reads multiple CSV files containing mutation data and merges them into a single DataFrame.
    It then performs additional processing by merging the merged DataFrame with another DataFrame.
    Finally, it saves the resulting DataFrame to a CSV file.

    Args:
        None

    Returns:
        None
    """
    setup_applevel_logger()
    intput_path = "data/construct_design/final_sets"
    dfs = []
    for i in range(1, 5):
        df = pd.read_csv(f"{intput_path}/mtt6_mutations_set_{i}/results-rna.csv")
        dfs.append(df)
    df = pd.concat(dfs)
    df_org = pd.read_csv("data/construct_design/sets/all_sets.csv")
    df_org = df_org[["name", "dg", "dg_err", "act_seq", "act_ss"]]
    df = df.merge(df_org, on="name")
    df.to_csv("data/construct_design/constructs.csv", index=False)


# data processing ####################################################################


@cli.command()
def get_sequencing_data():
    """
    Moves data into the 'data/sequencing_runs/raw' folder.
    Note: This function requires access to all Yesselman lab data.

    Returns:
        None

    Raises:
        FileNotFoundError: If the specified file is not found.
        Exception: If there is an error reading the CSV file.

    """
    setup_applevel_logger()
    local_data_path = "data/"
    sequencing_runs_path = os.path.join(local_data_path, "sequencing_runs", "raw")
    os.makedirs(sequencing_runs_path, exist_ok=True)
    log.info("DATA_PATH: " + local_data_path)
    # put this in a param file
    path = "/Users/jyesselman2/Dropbox/data/sequencing_analysis/"
    try:
        df = pd.read_csv(os.path.join(RESOURCES_PATH, "csvs", "sequencing_runs.csv"))
    except Exception as e:
        log.error(f"Failed to read CSV file: {e}")
        return
    for _, row in df.iterrows():
        full_path = os.path.join(path, row["run_name"], "analysis", "summary.json")
        if not os.path.isfile(full_path):
            log.warning(f"File {row['run_name']} is not found at {full_path}")
            continue
        else:
            log.info(f"File {row['run_name']} is found at {full_path} copying")
        destination_path = os.path.join(sequencing_runs_path, f"{row['run_name']}.json")
        shutil.copy(full_path, destination_path)


@cli.command()
def process_sequencing_data():
    """
    Analyzes the sequencing data and saves the processed data in
    data/sequencing_runs/processed.

    This function sets up the application-level logger, creates the necessary
    directory for storing the processed data, and processes the sequencing data
    using a list of data processors. Each data processor is responsible for
    loading, cleaning, and processing a specific type of data.

    Data Processors:
    - MTTR6BufferTitrationDataProcessor: Processes MTTR6 buffer titration data.
    - MTTR6MgTitrationDataProcessor: Processes MTTR6 Mg titration data.
    - MTTR6MutsDataProcessor: Processes MTTR6 mutations data.
    - TTRMutsDataProcessor: Processes TTR mutations data.
    """
    setup_applevel_logger()
    os.makedirs("data/sequencing_runs/processed", exist_ok=True)
    data_processors = [
        MTTR6BufferTitrationDataProcessor(),
        MTTR6MgTitrationDataProcessor(),
        MTTR6MutsDataProcessor(),
        TTRMutsDataProcessor(),
    ]
    for p in data_processors:
        log.info("starting to process: " + p.name)
        p.load_data()
        p.clean_data()
        p.process_data()


@cli.command()
def characterize_mutations():
    """
    characterize steve's TLR mutations
    """
    setup_applevel_logger()
    path = "data/construct_design/sets/all_sets.csv"
    df = pd.read_csv(path)
    df_results = mutation_characterize.characterize_mutations(df)
    df_results.to_json("data/mttr6_mut_charactization.json")


@cli.command()
def compute_sasa():
    state_1_pdb_path = "data/pdbs/sasa/no_mg_no_gaaa_state_1.pdb"
    # resi 6 = A4, resi 7 = A5, resi 18 = A8
    df_sasa_1 = compute_solvent_accessability(state_1_pdb_path)
    keep = [6, 7, 18]
    df_sasa_1 = df_sasa_1[df_sasa_1["Nt_num"].isin(keep)]
    df_sasa_1["state"] = 1
    # renumber the residues
    df_sasa_1["Nt_num"] = df_sasa_1["Nt_num"].map({6: 4, 7: 5, 18: 8})
    state_2_pdb_path = "data/pdbs/sasa/mg_no_gaaa_state_2.pdb"
    df_sasa_2 = compute_solvent_accessability(state_2_pdb_path)
    keep = [226, 225, 248]
    df_sasa_2 = df_sasa_2[df_sasa_2["Nt_num"].isin(keep)]
    df_sasa_2["state"] = 2
    # renumber the residues
    df_sasa_2["Nt_num"] = df_sasa_2["Nt_num"].map({226: 4, 225: 5, 248: 8})
    state_3_pdb_path = "data/pdbs/sasa/mg_gaaa_state_3.pdb"
    df_sasa_3 = compute_solvent_accessability(state_3_pdb_path)
    df_sasa_3 = df_sasa_3[df_sasa_3["Nt_num"].isin(keep)]
    df_sasa_3["state"] = 3
    df_sasa_3["Nt_num"] = df_sasa_3["Nt_num"].map({226: 4, 225: 5, 248: 8})
    df_sasa = pd.concat([df_sasa_1, df_sasa_2, df_sasa_3])
    df_sasa = df_sasa.rename(columns={"Nt_num": "resi"})
    df_sasa.to_csv("data/tables/tlr_sasa.csv", index=False)


@cli.command()
def process_mg_1_2():
    """
    Compute mg 1/2 for ttr mutants
    """
    setup_applevel_logger()
    path = "q_dms_ttr_paper/data/processed/mttr6_data_full.json"
    if not os.path.isfile(path):
        log.error("File not found: " + path)
        return
    df = pd.read_json(path)
    df_results = compute_all_mg_1_2(df)
    df_results.to_csv(
        "q_dms_ttr_paper/data/processed/mtt6_data_mg_1_2.csv", index=False
    )


@cli.command()
def generate_3d_modeling_runs():
    pass


# plot generation ####################################################################


@cli.command()
def plot_raw_mut_fractions():
    """
    plot raw mut fractions. Shows the full pop avg for eac titration point
    """
    setup_applevel_logger()
    os.makedirs("plots/titration_full_construct", exist_ok=True)
    # mg2+ titrations
    json_path = "q_dms_ttr_paper/data/processed/wt_mg_titra.json"
    df = pd.read_json(json_path)
    highlights = []
    highlights.append({"motif": {"name": "gaaa_tetraloop"}})
    highlights.append({"motif": {"name": "tlr"}})
    path = f"plots/titration_full_construct/wt_mg_titra/"
    os.makedirs(path, exist_ok=True)
    for (name, exp_name), g in df.groupby(["name", "exp_name"]):
        os.makedirs(f"{path}/{exp_name}", exist_ok=True)
        g = trim(g, 20, 20)
        plot_pop_avg_titration(g, "mg_conc", highlights, figsize=(12, 14))
        plt.tight_layout()
        plt.savefig(f"{path}/{exp_name}/{name}.png", dpi=300)
        plt.clf()
    # buffer titrations
    json_path = "q_dms_ttr_paper/data/processed/wt_buffer_titra.json"
    df = pd.read_json(json_path)
    path = f"plots/titration_full_construct/wt_buffer_titra/"
    os.makedirs(path, exist_ok=True)
    for (name, buffer, exp_name), g in df.groupby(["name", "buffer", "exp_name"]):
        os.makedirs(f"{path}/{exp_name}_buffer_{buffer}", exist_ok=True)
        g = trim(g, 20, 20)
        plot_pop_avg_titration(g, "buffer_conc", highlights, figsize=(12, 14))
        plt.tight_layout()
        plt.savefig(f"{path}/{exp_name}_buffer_{buffer}/{name}.png", dpi=300)
        plt.clf()
    # mttr6 muts
    json_path = "q_dms_ttr_paper/data/processed/mttr6_muts_titra.json"
    df = pd.read_json(json_path)
    path = f"plots/titration_full_construct/mttr6_muts_titra/"
    os.makedirs(path, exist_ok=True)
    for (name, exp_name), g in df.groupby(["name", "exp_name"]):
        highlights = []
        if name != "minittr_6_uucg_fixed":
            highlights.append({"motif": {"name": "gaaa_tetraloop"}})
        if name != "minittr_6_no_tlr_fixed":
            highlights.append({"motif": {"name": "tlr"}})
        g = trim(g, 20, 0)
        plot_pop_avg_titration(g, "mg_conc", highlights, figsize=(12, 14))
        plt.tight_layout()
        plt.savefig(f"{path}/{name}.png", dpi=300)
        plt.clf()
    # tlr muts
    highlights = []
    highlights.append({"motif": {"name": "gaaa_tetraloop"}})
    json_path = "q_dms_ttr_paper/data/processed/mttr6_data_full.json"
    df = pd.read_json(json_path)
    path = f"plots/titration_full_construct/mttr6_data_full/"
    os.makedirs(path, exist_ok=True)
    for name, g in df.groupby(["name"]):
        # g = trim(g, 20, 0)
        plot_pop_avg_titration(g, "mg_conc", highlights, figsize=(12, 14))
        plt.tight_layout()
        plt.savefig(f"{path}/{name}.png", dpi=300)
        plt.clf()


@cli.command()
def plot_titrations():
    """
    plot titration as a function of mg2+ taking the average of gaaa tetraloop avg
    """
    setup_applevel_logger()
    os.makedirs("plots/titration_fits", exist_ok=True)
    json_path = "q_dms_ttr_paper/data/processed/wt_mg_titra.json"
    df = pd.read_json(json_path)
    path = f"plots/titration_fits/wt_mg_titra/"
    os.makedirs(path, exist_ok=True)
    for (name, exp_name), g in df.groupby(["name", "exp_name"]):
        os.makedirs(f"{path}/{exp_name}", exist_ok=True)
        pfit, perr = compute_mg_1_2(g["mg_conc"], g["gaaa_avg"])
        plot_mg_titration_fit(g["mg_conc"], g["gaaa_avg"], pfit[0], pfit[1], pfit[2])
        plt.title(f"{name} - mg_1_2: {round(pfit[0], 3)} - n: {round(pfit[1], 2)}")
        plt.savefig(f"{path}/{exp_name}/{name}_fit.png", dpi=300)
        plt.clf()
    # mttr6 muts
    json_path = "q_dms_ttr_paper/data/processed/mttr6_muts_titra.json"
    df = pd.read_json(json_path)
    path = f"plots/titration_fits/mttr6_muts_titra/"
    os.makedirs(path, exist_ok=True)
    for (name, exp_name), g in df.groupby(["name", "exp_name"]):
        pfit, perr = compute_mg_1_2(g["mg_conc"], g["gaaa_avg"])
        plot_mg_titration_fit(g["mg_conc"], g["gaaa_avg"], pfit[0], pfit[1], pfit[2])
        plt.title(f"{name} - mg_1_2: {round(pfit[0], 3)} - n: {round(pfit[1], 2)}")
        plt.savefig(f"{path}/{name}_fit.png", dpi=300)
        plt.clf()
    # tlr muts
    json_path = "q_dms_ttr_paper/data/processed/mttr6_data_full.json"
    df = pd.read_json(json_path)
    # todo this is only for the first set
    df = df[df["mg_conc"] != 5.0]
    path = f"plots/titration_fits/mttr6_data_full/"
    os.makedirs(path, exist_ok=True)
    for name, g in df.groupby(["name"]):
        pfit, perr = compute_mg_1_2(g["mg_conc"], g["gaaa_avg"])
        plot_mg_titration_fit(g["mg_conc"], g["gaaa_avg"], pfit[0], pfit[1], pfit[2])
        plt.title(f"{name} - mg_1_2: {round(pfit[0], 3)} - n: {round(pfit[1], 2)}")
        plt.savefig(f"{path}/{name}_fit.png", dpi=300)
        plt.clf()


@cli.command()
def plot_other_motif_diff():
    def _plot(df, m_name):
        for i, (name, g) in enumerate(df.groupby("name")):
            plt.title(name, fontsize=20)
            plt.scatter(g["wt_val"], g["uucg_val"])
            plt.xlabel("WT Mut. Frac.", fontsize=16)
            plt.ylabel("TL-Knockout Mut. Frac.", fontsize=16)
            plt.axis("square")
            plt.savefig(f"plots/other_motifs/{m_name}-{name}.png", dpi=300)
            plt.clf()

    setup_applevel_logger()
    os.makedirs("plots/other_motifs", exist_ok=True)
    json_path = "q_dms_ttr_paper/data/processed/wt_mg_titra.json"
    df_wt = pd.read_json(json_path)
    df_wt = df_wt[df_wt["exp_name"] == "2022_07_27_C0117_50mM_NaC_Mg2+_titra_CM"]
    json_path = "q_dms_ttr_paper/data/processed/mttr6_muts_titra.json"
    df_muts = pd.read_json(json_path)
    df_uucg = df_muts[df_muts["name"] == "minittr_6_uucg_fixed"]
    df_uucg = df_uucg.sort_values("mg_conc")
    # ires
    sequence = "GAACUACGC"
    structure = "(.....())"
    df = get_motif_data(df_wt, df_uucg, "ires", sequence, structure)
    _plot(df, "ires")
    sequence = "CCGAG&CGUUUGACG"
    structure = "(((((&)..)))..)"
    df = get_motif_data(df_wt, df_uucg, "kink_turn", sequence, structure)
    _plot(df, "kink_turn")
    sequence = "GAACA&UACCC"
    structure = "(...(&)...)"
    df = get_motif_data(df_wt, df_uucg, "3x3_motif", sequence, structure)
    _plot(df, "3x3_motif")


@cli.command()
@click.option("-n", "--name", default="minittr_6_uucg_fixed")
@click.option("-m", "--mg_conc", default=10.0)
def plot_mut_histo(name, mg_conc):
    data_files = [
        "mttr6_data_full.json",
        "mttr6_muts_titra.json",
        "wt_buffer_titra.json",
        "wt_mg_titra.json",
    ]
    dfs = []
    for f in data_files:
        path = os.path.join(DATA_PATH, "processed", f)
        df = pd.read_json(path)
        dfs.append(df)
    df = pd.concat(dfs)
    df = df[df["name"] == name]
    df = df[df["mg_conc"] == mg_conc]
    plot_pop_avg_from_row(df.iloc[0])
    # plt.ylim(0, 0.08)
    plt.show()


if __name__ == "__main__":
    cli()
