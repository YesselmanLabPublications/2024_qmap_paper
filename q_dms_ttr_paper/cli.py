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
    trim,
)
from q_dms_ttr_paper.plotting import plot_pop_avg_titration

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
    p = MTTR6BufferTitrationDataProcessor()
    p.load_data()
    p.clean_data()
    p.process_data()


@cli.command()
def plot_raw_mut_fractions():
    """
    plot raw mut fractions
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


if __name__ == "__main__":
    cli()
