"""
generates a txt file that has the dms values that will be read in by pymol to color the structure
"""
import click
import pandas as pd

from seq_tools.dataframe import trim


def get_data_dataframe(runs):
    dropbox_path = "/Users/jyesselm/Dropbox"
    runs = [run for run in runs.split(",")]
    path = f"{dropbox_path}/data/sequencing_analysis/summary/"
    dfs = []
    for run in runs:
        print(f"{path}{run}.json")
        df = pd.read_json(f"{path}{run}.json")
        dfs.append(df)
    df = pd.concat(dfs)
    return df

def generate_pymol_file(row, pymol_data_file):
    """
    generate the pymol file for a single row
    """
    data = row["data"]
    with open(pymol_data_file, "w") as f:
        for i, d in enumerate(data):
            d = d*25
            f.write(f"{i+1} {d} {d}\n")

@click.command()
@click.argument("runs", type=str)
@click.option("-b", "--buffer", type=str, default=None)
@click.option("-e", "--exp", type=str, default=None)
@click.option("-mg", "--mg-conc", type=float, default=None)
@click.option("--trim-5p", type=int, default=0)
@click.option("--trim-3p", type=int, default=0)
def main(runs, exp, buffer, mg_conc, trim_5p, trim_3p):
    """
    main function for script
    """
    df = get_data_dataframe(runs)
    if exp is not None:
        df = df[df["exp_name"] == exp]
    if buffer is not None:
        df = df[df["buffer"] == buffer]
    if mg_conc is not None:
        df = df[df["mg_conc"] == mg_conc]
    df = trim(df, trim_5p, trim_3p)
    if trim_5p is not None:
        df["data"] = df["data"].apply(lambda x: x[trim_5p:])
    if trim_3p is not None:
        df["data"] = df["data"].apply(lambda x: x[:-trim_3p])
    if len(df) == 0:
        raise ValueError("no data found")
    print(f"there are {len(df)} rows in the dataframe")
    row = df.iloc[0]
    generate_pymol_file(row, "pymol_data.txt")


# pylint: disable=no-value-for-parameter
if __name__ == '__main__':
    main()
