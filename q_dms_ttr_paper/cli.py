import click
import os
import pandas as pd


from q_dms_ttr_paper.logger import setup_applevel_logger, get_logger
from q_dms_ttr_paper.paths import RESOURCES_PATH

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
    os.makedirs("data/raw/sequencing_runs", exist_ok=True)
    path = "/Users/jyesselman2/Dropbox/data/sequencing_analysis/"
    df = pd.read_csv(RESOURCES_PATH + "/sequencing_runs.csv")
    for _, row in df.iterrows():
        full_path = path + row["run_name"] + "/analysis/summary.json"
        if not os.path.isfile(full_path):
            log.warn(f"File {row['run_name']} is not found at {full_path}")
            continue
        else:
            log.info(f"File {row['run_name']} is found at {full_path} copying")
        os.system(f"cp {full_path} data/raw/sequencing_runs/{row['run_name']}.json")


@cli.command()
def process_sequencing_data():
    """
    process sequencing data
    """
    setup_applevel_logger()


if __name__ == "__main__":
    cli()
