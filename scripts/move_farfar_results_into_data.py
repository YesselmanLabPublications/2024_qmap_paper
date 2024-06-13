import pandas as pd
import os
import shutil
import glob


def main():
    """
    main function for script
    """
    path = "data/construct_design/constructs.csv"
    df = pd.read_csv(path)
    org_path = "analysis/2022_08_10_rebuild_ttr/cluster_runs/"
    data_path = "data/farfar2_models/"
    files = ["default.out", "test.fasta", "test.secstruct_file", "start.pdb"]
    for i, row in df.iterrows():
        act_seq_name = row["act_seq"].replace("&", "_")
        org_dir = f"{org_path}{act_seq_name}"
        sbatch_file = glob.glob(f"{org_dir}/*.sbatch")[0]
        os.makedirs(f"{data_path}{row['name']}", exist_ok=True)
        shutil.copy(sbatch_file, f"{data_path}{row['name']}/{row['name']}.sbatch")
        for file in files:
            shutil.copy(f"{org_dir}/{file}", f"{data_path}{row['name']}/{file}")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
