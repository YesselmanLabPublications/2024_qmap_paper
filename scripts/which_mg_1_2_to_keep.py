import pandas as pd


def main():
    """
    main function for script
    """
    supp_df = pd.read_csv("data/tables/supplemental_table_1.csv")
    supp_df.drop(columns=["mg_1_2"], inplace=True)
    existing_df = pd.read_csv("final_mg_1_2.csv")
    existing_df["name"] = existing_df["name"].apply(
        lambda x: x.split("_")[1] + "_" + x.split("_")[0]
    )
    merged = existing_df[["name", "mg_1_2"]].merge(supp_df, on="name", how="left")
    merged.to_csv("data/mg_1_2_fits/mtt6_data_mg_1_2_final.csv", index=False)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
