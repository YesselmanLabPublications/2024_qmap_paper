import click
import os


from q_dms_ttr_paper.logger import setup_applevel_logger, get_logger


@click.group()
def cli():
    """ """
    pass


def get_sequencing_data():
    """
    move data into this folder data/raw/sequencing_runs. Note will not work without access
    to all Yesselman lab data.
    """
    setup_applevel_logger()


if __name__ == "__main__":
    cli()
