"""
moves the org sets into the data directory so the final constructs will be in the right place
"""

import pandas as pd
from seq_tools import add, to_dna_template, to_fasta, to_dna


def main():
    """
    main function for script
    """
    path = "/Users/jyesselman2/Dropbox/projects/construct_design/2022/2022_01_27_steve_ttrs/results"
    output_path = "data/construct_design/final_sets/"
    for i in range(1, 5):
        df = pd.read_csv(f"{path}/mtt6_mutations_set_{i}/results-all.csv")
        # update names to be in the correct format
        for j, row in df.iterrows():
            df.at[j, "name"] = (
                row["name"].split("_")[1] + "_" + row["name"].split("_")[0]
            )
        # put files into data directory
        # save results-all.csv
        df.to_csv(f"{output_path}/mtt6_mutations_set_{i}/results-all.csv", index=False)
        df_input = df[
            ["name", "org_sequence", "org_structure", "org_ens_defect"]
        ].copy()
        # fix input need to add back on first stem
        df_input.rename(
            columns={
                "org_sequence": "sequence",
                "org_structure": "structure",
                "org_ens_defect": "ens_defect",
            },
            inplace=True,
        )
        df_input = add(df_input, "GAGC", "GCUC")
        df_input.to_csv(f"{output_path}/mtt6_mutations_set_{i}/input.csv", index=False)
        # create results-rna.csv
        df_rna = df[["name", "sequence", "structure", "ens_defect"]].copy()
        df_rna.to_csv(
            f"{output_path}/mtt6_mutations_set_{i}/results-rna.csv", index=False
        )
        # create results-dna.csv
        df_dna = df[["name", "sequence"]].copy()
        df_dna = to_dna_template(df_dna)
        df_dna.to_csv(
            f"{output_path}/mtt6_mutations_set_{i}/results-dna.csv", index=False
        )
        # create fasta files
        to_fasta(to_dna(df_rna), f"{output_path}/mtt6_mutations_set_{i}/results.fasta")
        # create opool files
        df_opool = df_dna[["name", "sequence"]].copy()
        df_opool["name"] = f"mtt6_mutations_set_{i}"
        df_opool.to_csv(
            f"{output_path}/mtt6_mutations_set_{i}/results-opool.csv", index=False
        )
        df_opool.to_excel(
            f"{output_path}/mtt6_mutations_set_{i}/results-opool.xlsx", index=False
        )


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
