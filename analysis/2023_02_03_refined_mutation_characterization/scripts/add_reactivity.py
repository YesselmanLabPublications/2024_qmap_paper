"""
adds reactivity to classifed mut dataframe.
"""
import json
import pandas as pd

from seq_tools import structure

def get_dataframe():
    dropbox_path = "/Users/jyesselm/Dropbox"
    runs = [
        "2022_08_25_mtt6_set4_1st3_seq",
        "2022_08_26_mtt6_set4_2nd3_seq",
        "2022_08_29_mtt6_seq",
        "2022_09_01_mtt6_set4_seq",
        "2022_08_31_mtt6_set4_seq",
        "2022_08_26_mtt6_set1-3_MgTitra_KU_seq",
    ]
    path = f"{dropbox_path}/data/sequencing_analysis/summary/"
    dfs = []
    for run in runs:
        df = pd.read_json(f"{path}{run}.json")
        dfs.append(df)
    df = pd.concat(dfs)
    df_ref = pd.read_json("classified_muts.json")
    df = df.merge(df_ref)
    return df


def assign_residue_reactivities(row, data):
    # print(data)
    # print(len(data))
    # print(json.dumps(row.to_dict(), indent=4))
    #print(row['aligned_seq'])
    #print(len(row['aligned_seq']))
    #seq = row["act_seq"].split("&")[1][:-1] + row["act_seq"].split("&")[0][1:]
    #print(len(seq))
    res_data = {}
    insert_pos = []
    for insert in row["insertions"]:
        pos = int(insert[:-1])
        insert_pos.append(pos)
    mapped_pos = 1
    pos = 0
    while pos < len(data):
        if mapped_pos in insert_pos:
            res_data[f"data_{mapped_pos}_ins"] = data[pos]  / row['ref_hp_as']
            res_data[f"data_{mapped_pos}"] = data[pos+1]  / row['ref_hp_as']
            pos += 1
        else:
            res_data[f"data_{mapped_pos}"] = data[pos] / row['ref_hp_as']
        mapped_pos += 1
        pos += 1
    return res_data


def updata_dataframe(df) -> pd.DataFrame:
    df.drop(
        ["m_data", "m_sequence", "m_structure", "m_ind"], axis=1, inplace=True
    )
    all_data = []
    for i, row in df.iterrows():
        full_ss = structure.SequenceStructure(row["sequence"], row["structure"])
        ss = structure.SequenceStructure(row["act_seq"], row["act_ss"])
        bounds = structure.find(full_ss, ss)[0]
        # trims off last basepair of 5' and first basepair of 3'
        # flips the orientation to match aligend sequence and structure
        data = (
            row["data"][bounds[1][0]: bounds[1][1] -1] +
            row["data"][bounds[0][0] + 1 : bounds[0][1]]
        )
        res_data = assign_residue_reactivities(row, data)
        res_data["m_data"] = data
        all_data.append(res_data)
    df_m = pd.DataFrame(all_data)
    df = pd.concat([df, df_m], axis=1)
    return df


def main():
    """
    main function for script
    """
    df = get_dataframe()
    df = updata_dataframe(df)
    df.to_json("classified_muts_w_reactivity_normed.json", orient="records")

# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
