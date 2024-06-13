import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict

from seq_tools import has_5p_sequence, to_rna
from seq_tools import trim as seq_ss_trim
from seq_tools import SequenceStructure
from seq_tools.structure import find as seq_ss_find
from rna_secstruct import SecStruct, MotifSearchParams


from qmap_paper.logger import get_logger
from qmap_paper.titration import compute_mg_1_2


log = get_logger("DATA-PROCESSING")

DATA_PATH = "data"
RESOURCES_PATH = "qmap_paper/resources"


def get_data(path: str, sets: list) -> pd.DataFrame:
    """
    Reads multiple JSON files and concatenates them into a single DataFrame.

    Args:
        path (str): The directory path where the JSON files are located.
        sets (list): A list of filenames (without extensions) to read and concatenate.

    Returns:
        pd.DataFrame: A concatenated DataFrame containing the data from all specified JSON files.

    Raises:
        FileNotFoundError: If any of the specified JSON files do not exist in the given path.

    Example:
        >>> path = "/data/json_files"
        >>> sets = ["file1", "file2", "file3"]
        >>> df = get_data(path, sets)
        >>> print(df.head())
    """
    dfs = []
    for run_name in sets:
        full_path = Path(path) / f"{run_name}.json"
        if not full_path.exists():
            raise FileNotFoundError(f"File {full_path} does not exist")
        df = pd.read_json(full_path)
        dfs.append(df)
    return pd.concat(dfs)


def trim(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """
    Trims the 'sequence', 'structure', and 'data' columns of the DataFrame to the given start and end indices.

    Args:
        df (pd.DataFrame): A DataFrame with 'sequence', 'structure', and 'data' columns, where 'data' contains lists of numbers.
        start (int): The start index for trimming.
        end (int): The end index for trimming.

    Returns:
        pd.DataFrame: A trimmed DataFrame with the 'sequence', 'structure', and 'data' columns adjusted to the specified indices.

    Example:
        >>> df = pd.DataFrame({
        ...     "sequence": ["ABCDEFG", "HIJKLMN", "OPQRSTU"],
        ...     "structure": ["1234567", "2345678", "3456789"],
        ...     "data": [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]
        ... })
        >>> trimmed_df = trim(df, 1, 2)
        >>> print(trimmed_df)
          sequence structure         data
        0     BCDEF    23456     [2, 3, 4]
        1     IJKLM    34567     [7, 8, 9]
        2     PQRST    45678  [12, 13, 14]
    """
    df = seq_ss_trim(df, start, end)
    if start == 0 and end != 0:
        df["data"] = df["data"].apply(lambda x: x[:-end])
    elif end == 0 and start != 0:
        df["data"] = df["data"].apply(lambda x: x[start:])
    elif start == 0 and end == 0:
        df["data"] = df["data"].apply(lambda x: x)
    else:
        df["data"] = df["data"].apply(lambda x: x[start:-end])
    return df


def trim_p5_and_p3(df: pd.DataFrame) -> pd.DataFrame:
    """
    Trims the 5' and 3' ends of the data in the DataFrame.

    This function reads a CSV file containing p5 sequences, converts these sequences to RNA,
    checks for a common p5 sequence in the given DataFrame, and trims the DataFrame based on
    the length of this common p5 sequence and a fixed 3' end length.

    Args:
        df (pd.DataFrame): A DataFrame with a 'data' column containing sequences as strings.

    Returns:
        pd.DataFrame: A trimmed DataFrame with the 5' and 3' ends trimmed.

    Raises:
        ValueError: If no common p5 sequence is found or the sequence is not registered in the CSV file.

    Example:
        >>> df = pd.DataFrame({"data": ["GGAAGATCGAGTAGATCAAAGCATGC", "GGAAGATCGAGTAGATCAAAGCATGC", "GGAAGATCGAGTAGATCAAAGCATGC"]})
        >>> trimmed_df = trim_p5_and_p3(df)
        >>> print(trimmed_df)
           data
        0  GCATGCAT
        1  GCATGCAT
        2  GCATGCAT
    """
    df_p5 = pd.read_csv(f"{RESOURCES_PATH}/csvs/p5_sequences.csv")
    df_p5 = to_rna(
        df_p5
    )  # Ensure to_rna is defined elsewhere or import it if it's external
    common_p5_seq = ""
    for p5_seq in df_p5["sequence"]:
        if has_5p_sequence(
            df, p5_seq
        ):  # Ensure has_5p_sequence is defined elsewhere or import it if it's external
            common_p5_seq = p5_seq
    if len(common_p5_seq) == 0:
        raise ValueError("No common p5 sequence found")
    log.info(f"common p5 sequence: {common_p5_seq}")
    return trim(
        df, len(common_p5_seq), 20
    )  # Ensure trim is defined elsewhere or import it if it's external


def get_dms_reactivity_for_motif(
    df: pd.DataFrame, params: Dict, error: bool = True
) -> List[List[int]]:
    """
    Finds a specific sequence and structure motif in each construct and returns the corresponding DMS reactivity data.

    Args:
        df (pd.DataFrame): A DataFrame containing 'sequence', 'structure', and 'data' columns.
        params (dict): A dictionary of parameters for motif searching. Should be compatible with `MotifSearchParams`.
        error (bool): If True, raises an error when more than one motif is found. Default is True.

    Returns:
        List[List[int]]: A list of lists, where each inner list contains the DMS reactivity data for the found motif.
                         If no motif is found, the list contains -1 for each expected data point.

    Raises:
        ValueError: If more than one motif is found in a construct when `error` is True.

    Example:
        >>> df = pd.DataFrame({
             "sequence": ["ACGUACGUACGU", "GCUAGCUAGCUA"],
             "structure": ["(((...)))...", ".((...))..."],
             "data": [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
                  [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]]
            })
        >>> params = {
            "sequence": "ACGU&CGUA",
             "structure": "(((...)))"
            }
        >>> reactivity_data = get_dms_reactivity_for_motif(df, params, false)
        >>> print(reactivity_data)
        [[0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9], [-1, -1, -1, -1, -1, -1, -1]]

    """
    params = MotifSearchParams(**params)
    data_len = len(params.sequence) - params.sequence.count("&")
    all_data = []
    for _, row in df.iterrows():
        ss = SecStruct(
            row["sequence"], row["structure"]
        )  # Ensure SecStruct is defined elsewhere or import it if it's external
        motifs = ss.get_motifs(params)  # Ensure get_motifs is a method of SecStruct
        if len(motifs) != 1 and error:
            raise ValueError("More than one motif found")
        elif len(motifs) == 0:
            all_data.append([-1] * data_len)
            continue
        data = []
        for s in motifs[0].strands:
            for e in s:
                data.append(row["data"][e])
        all_data.append(data)
    return all_data


def get_dms_reactivity_for_sub_structure(
    df: pd.DataFrame,
    sub_seq_struct: SequenceStructure,
    start=None,
    end=None,
    error=True,
) -> List[List[float]]:
    """
    Retrieves the DMS reactivity data for a given sub-sequence and structure from a DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing the data.
        sub_seq_struct (SequenceStructure): The sub-sequence and structure to search for.
        start (int, optional): The starting position of the sub-sequence. Defaults to None.
        end (int, optional): The ending position of the sub-sequence. Defaults to None.
        error (bool, optional): Whether to raise an error if the sub-sequence and structure are not found. Defaults to True.
        If no error is raised it will return -1 for each expected data point.

    Returns:
        List[List[float]]: A list of lists containing the DMS reactivity data for each position in the sub-sequence.
    """
    all_data = []
    data_len = len(sub_seq_struct.sequence) - sub_seq_struct.sequence.count("&")
    for i, row in df.iterrows():
        ss = SequenceStructure(row["sequence"], row["structure"])
        r = seq_ss_find(ss, sub_seq_struct, start, end)
        if len(r) == 0 and error:
            msg = (
                f"Could not find seq:{sub_seq_struct.sequence} "
                f"ss:{sub_seq_struct.structure} in {row['name']} with seq "
                f"{row['sequence']} and ss {row['structure']}"
            )
            raise ValueError(msg)
        elif len(r) == 0:
            all_data.append([-1] * data_len)
            continue
        elif len(r) > 1:
            msg = (
                f"multiple copies of seq:{sub_seq_struct.sequence} "
                f"ss:{sub_seq_struct.structure} in {row['name']} with seq "
                f"{row['sequence']} and ss {row['structure']}"
            )
            raise ValueError(msg)
        pos = []
        bounds = r[0]
        for r in bounds:
            pos.extend(list(range(r[0], r[1])))
        all_data.append([row["data"][p] for p in pos])
    return all_data


def get_dms_reactivity_for_wt_tlr(df, start=None, end=None, error=True):
    """
    Retrieves the DMS reactivity data for the wild-type tetraloop receptor (wt_tlr) from a DataFrame.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
        start (int, optional): The starting position for searching the tetraloop receptor. Defaults to None.
        end (int, optional): The ending position for searching the tetraloop receptor. Defaults to None.
        error (bool, optional): Whether to raise an error if the tetraloop receptor is not found. Defaults to True.

    Returns:
        list: A list of lists containing the DMS reactivity data for each row in the DataFrame. Each inner list represents the reactivity data for the wt_tlr.

    Raises:
        ValueError: If the tetraloop receptor is not found in a row and error is set to True.
        ValueError: If multiple copies of the tetraloop receptor are found in a row.

    """
    # NOTE okay this is a problem have to search in a specfic direction.
    # m_ss = structure.SequenceStructure("CCUAAG&UAUGG", "((...(&)..))")
    m_ss = SequenceStructure("UAUGG&CCUAAG", "(..((&))...)")
    data_len = len(m_ss.sequence) - m_ss.sequence.count("&")
    all_data = []
    for i, row in df.iterrows():
        ss = SequenceStructure(row["sequence"], row["structure"])
        r = seq_ss_find(ss, m_ss, start, end)
        if len(r) == 0 and error:
            msg = (
                f"Could not find the tetraloop receptor "
                f"in {row['name']} with seq "
                f"{row['sequence']} and ss {row['structure']}"
            )
            raise ValueError(msg)
        elif len(r) == 0:
            all_data.append([-1] * data_len)
            continue
        elif len(r) > 1:
            msg = (
                f"multiple copies of the tetraloop receptor where found  "
                f"in {row['name']} with seq "
                f"{row['sequence']} and ss {row['structure']}"
            )
            raise ValueError(msg)
        bounds = r[0]
        data = (
            row["data"][bounds[1][0] : bounds[1][1]]
            + row["data"][bounds[0][0] : bounds[0][1]]
        )
        all_data.append(data)
    return all_data


# reused processing functions #########################################################


def get_gaaa_data(df, error=True):
    """
    generates two new columns in the dataframe one with the gaaa reactivity data
    and one with the average gaaa reactivity data
    :params df: a dataframe with data
    """
    df["gaaa"] = get_dms_reactivity_for_sub_structure(
        df, SequenceStructure("GGAAAC", "(....)"), error=error
    )
    df["gaaa_avg"] = df["gaaa"].apply(lambda x: np.mean(x[2:-1]))


# specific library processing #########################################################
# for chemical mapping sequencing runs ################################################


class DataProcessor:
    """A class for processing data."""

    def load_data(self):
        """
        Get data from the raw folder and load it into a dataframe.

        Returns:
            None
        """
        pass

    def clean_data(self):
        """
        Make sure we have the correct data by removing any duplicates and data that
        is not related to the project.

        Returns:
            None
        """
        pass

    def process_data(self):
        """
        Perform any processing for the data to be useful after. This includes
        getting data on motifs, mutational analysis, titration analysis, etc.

        Returns:
            None
        """
        pass


"""
MTTR6BufferTitrationDataProcessor - covers wild-type sequence with different buffer 
conditions

MTTR6MgTitrationDataProcessor - covers wild-type sequence with different mg conditions

MTTR6MutsDataProcessor - covers point mutants of the wild-type sequence such as uucg / 
no-tlr and their mg titrations

TTRMutsDataProcessor - covers the library of mutations from steves paper
"""


class MTTR6BufferTitrationDataProcessor(DataProcessor):
    """
    This class represents a data processor for MTTR6 buffer titration data.
    It covers the original data tests to find the ideal buffer conditions.

    Attributes:
        name (str): The name of the data processor.
        df (pd.DataFrame): The DataFrame containing the data.

    Methods:
        load_data: Loads the data from the specified runs.
        clean_data: Cleans the loaded data by removing irrelevant entries.
        process_data: Processes the cleaned data by extracting specific features.
    """

    def __init__(self):
        self.name = "MTTR6BufferTitrationDataProcessor"
        self.df: pd.DataFrame = None

    def load_data(self):
        """
        Loads the data from the specified runs.
        """
        runs = [
            "2022_07_26_minittr-6-2HP-ref_buffer_seq",
            "2022_07_20_minittr_Hepes-titra_seq",
        ]
        self.df = get_data(DATA_PATH + "/sequencing_runs/raw", runs)

    def clean_data(self):
        """
        Cleans the loaded data by removing irrelevant entries.

        Removes irrelevant entries from the loaded data based on specific criteria.
        The cleaning process includes the following steps:
        1. Converts the data to RNA format.
        2. Removes entries with names not included in the 'include' list.
        3. Removes entries with experiment names included in the 'exclude_exps' list.
        4. Removes entries with a specific experiment name.

        Note: This method modifies the 'df' attribute of the object.

        Args:
            None

        Returns:
            None
        """
        # ensure data is rna
        self.df = to_rna(self.df)
        # remove p5 and p3 sequence
        # remove constructs not related to this analysis
        include = ["minittr-6-2HP-ref"]
        self.df = self.df[self.df["name"].isin(include)]
        # these experiments dont seem to be related to this analysis
        exclude_exps = [
            "2022_07_18_C0117_100mM_buffer_Mg2+_titra_CM_BL"
            "2022_07_18_minittr-6-2HP-ref_NaC_titra_No_MgCl2_CM_RG",
            "2022_07_19_C0117_Hepes_titra_CM_BL",
        ]
        self.df = self.df[~self.df["exp_name"].isin(exclude_exps)]
        # why do I need this extra line??
        self.df = self.df[
            ~(
                self.df["exp_name"]
                == "2022_07_18_minittr-6-2HP-ref_NaC_titra_No_MgCl2_CM_RG"
            )
        ]

    def process_data(self):
        """
        Process the data by performing various operations on the DataFrame.

        This method performs the following operations:
        1. Removes common p5 and p3 sequences from the DataFrame.
        2. Retrieves GAAA tetraloop reactivity data.
        3. Calculates the reactivity for the wild-type TLR (tlr).
        4. Retrieves the reactivity for the first reference hairpin at the 5' end (ref_hp_1).
        5. Calculates the average reactivity for the first reference hairpin at the 5' end (ref_hp_1_as).
        6. Retrieves the reactivity for the second reference hairpin at the 3' end (ref_hp_2).
        7. Calculates the average reactivity for the second reference hairpin at the 3' end (ref_hp_2_as).
        8. Retrieves the reactivity for the ires motif.
        9. Retrieves the reactivity for the kinke turn motif.
        10. Retrieves the reactivity for the 3x3 motif.
        11. Saves the processed DataFrame to a JSON file.

        Returns:
            None
        """
        # remove common p5 and p3 sequences
        self.df = trim_p5_and_p3(self.df)
        # get GAAA tetraloop reactivity data
        get_gaaa_data(self.df)
        self.df["tlr"] = get_dms_reactivity_for_wt_tlr(self.df)
        # get the first reference hairpin at the 5' end
        self.df["ref_hp_1"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CGAGUAG", "(.....)"), end=50
        )
        self.df["ref_hp_1_as"] = self.df["ref_hp_1"].apply(
            lambda x: np.mean([x[2], x[5]])
        )
        # get the second reference hairpin at the 3' end
        self.df["ref_hp_2"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CGAGUAG", "(.....)"), start=50
        )
        self.df["ref_hp_2_as"] = self.df["ref_hp_2"].apply(
            lambda x: np.mean([x[2], x[5]])
        )
        # get the ires motif
        self.df["ires"] = get_dms_reactivity_for_motif(
            self.df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
        )
        # get kinke turn motif
        self.df["kink_turn"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
        )
        # get the 3x3 motif
        self.df["3x3_motif"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("GAACA&UACCC", "(...(&)...)")
        )

        self.df.to_json(
            "data/sequencing_runs/processed/wt_buffer_titra.json", orient="records"
        )


class MTTR6MgTitrationDataProcessor(DataProcessor):
    """Processes MTTR6 Mg Titration Data.

    This class is responsible for loading, cleaning, and processing MTTR6 Mg titration data.
    It inherits from the DataProcessor base class.

    Attributes:
        name (str): The name of the data processor.
        df (pd.DataFrame): The data frame containing the loaded data.

    """

    def __init__(self):
        self.name = "MTTR6MgTitrationDataProcessor"
        self.df: pd.DataFrame = None

    def load_data(self):
        """Loads the data from sequencing runs."""
        runs = [
            "2022_07_27_minittr_50mM_NaC_Mg_titra_seq",
            "2022_07_28_minittr_0.1M_NaC_Mg_titra_seq",
            "2022_07_29_minittr_0.15M_NaC_Mg_titr_seq",
            "2022_08_09_minittr_0.2M_NaC_Mg_titra_seq",
            "2022_08_10_minittr_0.25M_Mg_titr_seq",
            "2022_08_11_minittr_0.3M_NaC_Mg_titra_seq",
        ]
        self.df = get_data(DATA_PATH + "/sequencing_runs/raw/", runs)

    def clean_data(self):
        """Cleans the loaded data.

        This method performs data cleaning operations on the loaded data. It ensures that the data is in RNA format by
        calling the `to_rna` function. It also removes the p5 and p3 sequences and any constructs that are not related to
        the current analysis.

        Args:
            self: The instance of the class.

        Returns:
            None
        """
        # ensure data is rna
        self.df = to_rna(self.df)
        # remove p5 and p3 sequence
        # remove constructs not related to this analysis
        include = ["minittr-6-2HP-ref"]
        self.df = self.df[self.df["name"].isin(include)]

    def process_data(self):
        """Processes the cleaned data.

        This method performs various data processing steps on the cleaned data.
        It trims the P5 and P3 regions, retrieves motif data, calculates reactivity
        for different hairpin structures, and saves the processed data to a JSON file.

        Args:
            self (object): The instance of the class.

        Returns:
            None
        """
        self.df = trim_p5_and_p3(self.df)
        # get motif data
        # get GAAA tetraloop reactivity data
        get_gaaa_data(self.df)
        self.df["tlr"] = get_dms_reactivity_for_wt_tlr(self.df)
        # get the first reference hairpin at the 5' end
        self.df["ref_hp_1"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CGAGUAG", "(.....)"), end=50
        )
        self.df["ref_hp_1_as"] = self.df["ref_hp_1"].apply(
            lambda x: np.mean([x[2], x[5]])
        )
        # get the second reference hairpin at the 3' end
        self.df["ref_hp_2"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CGAGUAG", "(.....)"), start=50
        )
        self.df["ref_hp_2_as"] = self.df["ref_hp_2"].apply(
            lambda x: np.mean([x[2], x[5]])
        )
        # get the ires motif
        self.df["ires"] = get_dms_reactivity_for_motif(
            self.df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}
        )
        # get kinke turn motif
        self.df["kink_turn"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
        )
        # get the 3x3 motif
        self.df["3x3_motif"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("GAACA&UACCC", "(...(&)...)")
        )
        self.df.to_json(
            "data/sequencing_runs/processed/wt_mg_titra.json", orient="records"
        )


class MTTR6MutsDataProcessor(DataProcessor):
    """
    This class represents a data processor for the MTTR6Muts dataset.
    It covers the point mutants for the wild-type.

    Attributes:
        name (str): The name of the data processor.
        df (pd.DataFrame): The DataFrame to store the processed data.
    """

    def __init__(self):
        self.name = "MTTR6MutsDataProcessor"
        self.df: pd.DataFrame = None

    def load_data(self):
        """
        Loads the data for the MTTR6Muts dataset.
        """
        runs = [
            "2023_03_13_no_tlr_Mg_titra_redo_seq",
            "2023_02_02_minittr_6_uucg_Mg_titra_seq",
            "2023_02_17_no_3_3_junc_Mg_titr_seq",
            "2023_03_10_h1_3bp_longer_seq",
            "2023_03_14_no_ires_Mg_titr_redo_seq",
            "2023_03_15_h2_3bp_longer_Mg_titra_seq",
            "2023_03_22_h3_3bp_longer_Mg_titra_seq",
        ]
        self.df = get_data(DATA_PATH + "/sequencing_runs/raw", runs)
        # remove anything that does not have minittr in the name
        self.df = self.df[self.df["name"].str.contains("minittr")]

    def clean_data(self):
        """
        Cleans the loaded data.
        """
        pass

    def process_data(self):
        """
        Processes the cleaned data.

        This method performs various data processing steps on the cleaned data.
        It removes common p5 and p3 sequences, retrieves GAAA tetraloop reactivity data,
        calculates reactivity for specific sub-structures, and saves the processed data
        to a JSON file.

        Args:
            None

        Returns:
            None
        """
        # remove common p5 and p3 sequences
        self.df = trim_p5_and_p3(self.df)
        # get GAAA tetraloop reactivity data
        get_gaaa_data(self.df, error=False)
        self.df["tlr"] = get_dms_reactivity_for_wt_tlr(self.df, error=False)
        self.df["ref_hp_1"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CGAGUAG", "(.....)"), end=50
        )
        self.df["ref_hp_1_as"] = self.df["ref_hp_1"].apply(
            lambda x: np.mean([x[2], x[5]])
        )
        self.df["ires"] = get_dms_reactivity_for_motif(
            self.df, {"sequence": "GAACUAC&GC", "structure": "(.....(&))"}, error=False
        )
        self.df["3x3_motif"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("GAACA&UACCC", "(...(&)...)"), error=False
        )
        self.df["kink_turn"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CCGAG&CGUUUGACG", "(((((&)..)))..)")
        )
        self.df.to_json(
            "data/sequencing_runs/processed/mttr6_muts_titra.json", orient="records"
        )


class TTRMutsDataProcessor(DataProcessor):
    """
    This class represents a data processor for TTR mutations.

    Attributes:
        name (str): The name of the data processor.
        df (pd.DataFrame): The DataFrame containing the data.
    """

    def __init__(self):
        self.name = "TTRMutsDataProcessor"
        self.df: pd.DataFrame = None

    def load_data(self):
        """
        Loads the data into the DataFrame.
        """
        runs = [
            "2022_08_25_mtt6_set4_1st3_seq",
            "2022_08_26_mtt6_set4_2nd3_seq",
            "2022_08_26_mtt6_set1-3_MgTitra_KU_seq",
            "2022_08_29_mtt6_seq",
            "2022_08_30_mtt6_set4_seq",
            "2022_08_31_mtt6_set4_seq",
            "2022_09_01_mtt6_set4_seq",
        ]
        self.df = get_data(DATA_PATH + "/sequencing_runs/raw", runs)

    def clean_data(self):
        """
        Cleans the data by applying various transformations and filters.
        """
        # ensure data is rna
        self.df = to_rna(self.df)
        # fix some errors in the data
        # 0.20 was set instead of 0.25 for some constructs
        self.df.loc[self.df["mg_conc"] == 0.20, "mg_conc"] = 0.25
        # remove constructs not related to this analysis
        exclude = ["PurRe4L5", "PtrRe2L5", "PurRe2L4"]
        self.df = self.df.query("name not in @exclude")
        # duplicate data not sure right now which is better
        q1 = self.df["run_name"] == "2022_08_29_mtt6_seq"
        q2 = self.df["rna_name"] == "mtt6_mutations_set_2_50mM-NaC_0.1-mM-Mg2+"
        self.df = self.df[~(q1 & q2)]
        q1 = self.df["run_name"] == "2022_08_29_mtt6_seq"
        q3 = self.df["rna_name"] == "mtt6_mutations_set_1_50mM-NaC_5-mM-Mg2+"
        self.df = self.df[~(q1 & q3)]
        # check there are no duplicates
        duplicates = self.__get_duplicates(self.df)
        if len(duplicates) > 0:
            log.error("There are duplicates in the data")
            for d in duplicates:
                log.error(d)
            exit()
        # round data
        # for the data column round each value to 5 decimal places
        self.df["data"] = self.df["data"].apply(lambda x: [round(y, 5) for y in x])

    def process_data(self):
        """
        Processes the data by applying various transformations and calculations.
        """
        # remove common p5 and p3 sequences
        self.df = trim_p5_and_p3(self.df)
        # fix naming
        self.df["name"] = self.df["name"].apply(
            lambda x: x.split("_")[1] + "_" + x.split("_")[0]
        )
        # grab original rna_map experimental data from steves paper
        df_ref = pd.read_csv(f"data/construct_design/sets/all_sets.csv")
        df_ref = df_ref[["name", "dg", "dg_err", "act_seq", "act_ss"]]
        self.df = self.df.merge(df_ref, on="name")
        # get GAAA tetraloop reactivity data
        self.df["gaaa"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("GGAAAC", "(....)")
        )
        # get averaged value over the 3As in the tetraloop
        self.df["gaaa_avg"] = self.df["gaaa"].apply(lambda x: np.mean(x[2:-1]))
        # get reference hairpin at the 5' end of each construct
        self.df["ref_hp"] = get_dms_reactivity_for_sub_structure(
            self.df, SequenceStructure("CGAGUAG", "(.....)"), end=50
        )
        # get averaged value of the 2 As in the hairpin
        self.df["ref_hp_as"] = self.df["ref_hp"].apply(lambda x: np.mean([x[2], x[5]]))
        self.df = self.__get_tlr_reactivities(self.df)
        # make sure I save locally
        self.df.to_json(
            "data/sequencing_runs/processed/mttr6_data_full.json", orient="records"
        )

    def __get_duplicates(self, df):
        """
        Helper method to find duplicates in the DataFrame.

        Args:
            df (pd.DataFrame): The DataFrame to check for duplicates.

        Returns:
            list: A list of tuples containing the duplicate information.
        """
        duplicates = []
        for i, g in df.groupby(["name", "mg_conc"]):
            unique_runs = g["run_name"].unique()
            unique_rna = g["rna_name"].unique()
            if len(unique_runs) == 1:
                continue
            duplicates.append((i, unique_runs, unique_rna))
        return duplicates

    def __get_tlr_reactivities(self, df) -> pd.DataFrame:
        """
        Helper method to calculate TLR reactivities.

        Args:
            df (pd.DataFrame): The DataFrame containing the data.

        Returns:
            pd.DataFrame: The DataFrame with TLR reactivities added.
        """
        all_data = []
        for i, row in df.iterrows():
            full_ss = SequenceStructure(row["sequence"], row["structure"])
            ss = SequenceStructure(row["act_seq"], row["act_ss"])
            bounds = seq_ss_find(full_ss, ss)[0]
            # trims off last basepair of 5' and first basepair of 3'
            # flips the orientation to match aligend sequence and structure
            data = (
                row["data"][bounds[1][0] : bounds[1][1] - 1]
                + row["data"][bounds[0][0] + 1 : bounds[0][1]]
            )
            all_data.append(data)
        df["tlr"] = all_data
        return df


# compute mg 1 / 2 ###################################################################


def compute_all_mg_1_2(df):
    """
    Computes the mg 1/2 values for each construct in the mttr6 ttr mutant set.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data.

    Returns:
        pandas.DataFrame: A DataFrame with the computed mg 1/2 values for each construct.

    Raises:
        None

    Example:
        >>> df = pd.DataFrame({'name': ['A', 'A', 'B', 'B'], 'mg_conc': [1.0, 2.0, 1.0, 2.0], 'gaaa_avg': [0.5, 0.8, 0.3, 0.6]})
        >>> compute_all_mg_1_2(df)
           name  num_points  mg_1_2  mg_1_2_err  n  n_err  a_0  a_0_err
        0     A           2     1.5         0.1  3    0.2  0.4      0.05
        1     B           2     1.2         0.2  2    0.1  0.3      0.03
    """
    data = []
    df = df[df["mg_conc"] != 5.0]
    for name, g in df.groupby(["name"]):
        r = compute_mg_1_2(g["mg_conc"], g["gaaa_avg"])
        data.append(
            {
                "name": name,
                "num_points": len(g),
                "mg_1_2": r[0][0],
                "mg_1_2_err": r[1][0],
                "n": r[0][1],
                "n_err": r[1][1],
                "a_0": r[0][2],
                "a_0_err": r[1][2],
            }
        )
    df_results = pd.DataFrame(data)
    return df_results
