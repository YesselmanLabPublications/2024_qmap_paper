import pandas as pd
import numpy as np

from seq_tools import structure


def get_data():
    runs = [
        "2022_08_11_minittr_0.3M_NaC_Mg_titra_seq",
        "2022_08_10_minittr_0.25M_NaC_titr_seq",
        "2022_08_09_minittr_0.2M_NaC_Mg_titra_seq",
        "2022_07_29_minittr_0.15M_NaC_Mg_titr_seq",
        "2022_07_28_minittr_0.1M_NaC_Mg_titra_seq",
        "2022_07_27_minittr_50mM_NaC_Mg_titra_seq",
        "2022_07_26_minittr-6-2HP-ref_buffer_seq",
        "2022_07_20_minittr_Hepes-titra_seq",
        "2022_06_29_minittr-6-2HP-ref_seq",
    ]
    path = "/Users/jyesselm/Dropbox/data/sequencing_analysis/summary/"
    dfs = []
    for run in runs:
        df = pd.read_json(f"{path}{run}.json")
        dfs.append(df)
    df = pd.concat(dfs)
    return df


def assign_reactivity_in_tlr(full_ss, row):
    # TODO okay this is a problem have to search in a specfic direction.
    # m_ss = structure.SequenceStructure("CCUAAG&UAUGG", "((...(&)..))")
    m_ss = structure.SequenceStructure("UAUGG&CCUAAG", "(..((&))...)")
    try:
        bounds = structure.find(full_ss, m_ss)[0]
    except:
        return [-1 for _ in range(11)]
    data = (
        row["data"][bounds[1][0] : bounds[1][1]]
        + row["data"][bounds[0][0] : bounds[0][1]]
    )
    return data


def assign_reactivity_in_kink_turn(full_ss, row):
    m_ss = structure.SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
    try:
        bounds = structure.find(full_ss, m_ss)[0]
    except:
        return [-1 for _ in range(11)]
    data = (
        row["data"][bounds[0][0] : bounds[0][1]]
        + row["data"][bounds[1][0] : bounds[1][1]]
    )
    return data


def assign_ires_reactivity(full_ss, row):
    # TODO think more about this can I be really sure that the )) is the correct position?
    m_ss = structure.SequenceStructure("GAACUAC&GC", "(.....(&))")
    try:
        bounds = structure.find(full_ss, m_ss)[0]
    except:
        return [-1 for _ in range(11)]
    data = (
        row["data"][bounds[0][0] : bounds[0][1]]
        + row["data"][bounds[1][0] : bounds[1][1]]
    )
    return data


def assign_3_x_3_reactivity(full_ss, row):
    m_ss = structure.SequenceStructure("GAACA&UACCC", "(...(&)...)")
    try:
        bounds = structure.find(full_ss, m_ss)[0]
    except:
        return [-1 for _ in range(11)]
    data = (
        row["data"][bounds[0][0] : bounds[0][1]]
        + row["data"][bounds[1][0] : bounds[1][1]]
    )
    return data


def update_dataframe(df) -> pd.DataFrame:
    df.drop(
        ["m_data", "m_sequence", "m_structure", "m_ind"], axis=1, inplace=True
    )
    df["tlr_data"] = [[] for _ in range(len(df))]
    df["kink_turn_data"] = [[] for _ in range(len(df))]
    df["ires_data"] = [[] for _ in range(len(df))]
    df["3_x_3_data"] = [[] for _ in range(len(df))]
    for i, row in df.iterrows():
        full_ss = structure.SequenceStructure(row["sequence"], row["structure"])
        data = assign_reactivity_in_tlr(full_ss, row)
        df.at[i, "tlr_data"] = np.array(data)
        df.at[i, "kink_turn_data"] = np.array(
            assign_reactivity_in_kink_turn(full_ss, row)
        )
        df.at[i, "ires_data"] = np.array(assign_ires_reactivity(full_ss, row))
        df.at[i, "3_x_3_data"] = np.array(assign_3_x_3_reactivity(full_ss, row))
    return df


def main():
    """
    main function for script
    """
    df = get_data()
    new_index = pd.RangeIndex(start=0, stop=len(df), step=1)
    df = df.set_index(new_index)
    df = update_dataframe(df)
    df.to_json("C0117_titration_data.json", orient="records")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
