import click
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt


from q_dms_ttr_paper.logger import setup_applevel_logger, get_logger
from q_dms_ttr_paper.paths import RESOURCES_PATH, DATA_PATH
from q_dms_ttr_paper.data_processing import (
    MTTR6BufferTitrationDataProcessor,
    MTTR6MgTitrationDataProcessor,
    TTRMutsDataProcessor,
    MTTR6MutsDataProcessor,
    compute_all_mg_1_2,
    compute_mg_1_2,
    trim,
)
from q_dms_ttr_paper import mutation_characterize
from q_dms_ttr_paper.plotting import plot_pop_avg_titration, plot_mg_titration_fit

log = get_logger("CLI")


@click.group()
def cli():
    """ """
    pass


@cli.command()
def get_sequencing_data():
    """
    move data into this folder data/raw/sequencing_runs. Note will not work without
    access to all Yesselman lab data.
    """
    setup_applevel_logger()
    local_data_path = "q_dms_ttr_paper/data/"
    os.makedirs(local_data_path + "/raw", exist_ok=True)
    os.makedirs(local_data_path + "/raw/sequencing_runs", exist_ok=True)
    log.info("DATA_PATH: " + local_data_path)
    # put this in a param file
    path = "/Users/jyesselman2/Dropbox/data/sequencing_analysis/"
    df = pd.read_csv(RESOURCES_PATH + "/sequencing_runs.csv")
    for _, row in df.iterrows():
        full_path = path + row["run_name"] + "/analysis/summary.json"
        if not os.path.isfile(full_path):
            log.warning(f"File {row['run_name']} is not found at {full_path}")
            continue
        else:
            log.info(f"File {row['run_name']} is found at {full_path} copying")
        os.system(
            f"cp {full_path} {local_data_path}/raw/sequencing_runs/{row['run_name']}.json"
        )


@cli.command()
def process_sequencing_data():
    """
    process sequencing data
    """
    setup_applevel_logger()
    data_processors = [
        # MTTR6BufferTitrationDataProcessor(),
        # MTTR6MgTitrationDataProcessor(),
        MTTR6MutsDataProcessor(),
        # TTRMutsDataProcessor()
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
    path = "q_dms_ttr_paper/data/processed/rna_map/chemical_mapped.csv"
    df = pd.read_csv(path)
    df_results = mutation_characterize.characterize_mutations(df)
    df_results.to_json(
        "q_dms_ttr_paper/data/processed/mutations/mttr6_mut_charactization.json"
    )


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
    setup_applevel_logger()
    os.makedirs("plots/other_motifs", exist_ok=True)
    json_path = "q_dms_ttr_paper/data/processed/wt_mg_titra.json"
    df = pd.read_json(json_path)
    print(df["exp_name"].unique())


if __name__ == "__main__":
    cli()
