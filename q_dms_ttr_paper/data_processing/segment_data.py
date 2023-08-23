from seq_tools import SequenceStructure
from seq_tools.structure import find as seq_ss_find
from rna_secstruct import SecStruct, MotifSearchParams


def get_motif_data_in_dataframe(df, params):
    params = MotifSearchParams(**params)
    all_data = []
    for i, row in df.iterrows():
        ss = SecStruct(row["sequence"], row["structure"])
        motifs = ss.get_motifs(params)
        if len(motifs) != 1:
            raise ValueError("More than one motif found")
        data = []
        for s in motifs[0].strands:
            for e in s:
                data.append(row["data"][e])
        all_data.append(data)
    return all_data


def get_sec_struct_data_in_dataframe(df, sub_seq_struct, start=None, end=None):
    all_data = []
    for i, row in df.iterrows():
        ss = SequenceStructure(row["sequence"], row["structure"])
        r = seq_ss_find(ss, sub_seq_struct, start, end)
        if len(r) != 1:
            raise ValueError("More than one segment found")
        pos = []
        bounds = r[0]
        for r in bounds:
            pos.extend(list(range(r[0], r[1])))
        all_data.append([row["data"][p] for p in pos])
    return all_data


def get_wt_tlr_data_in_dataframe(df, start=None, end=None):
    # TODO okay this is a problem have to search in a specfic direction.
    # m_ss = structure.SequenceStructure("CCUAAG&UAUGG", "((...(&)..))")
    m_ss = SequenceStructure("UAUGG&CCUAAG", "(..((&))...)")
    all_data = []
    for i, row in df.iterrows():
        ss = SequenceStructure(row["sequence"], row["structure"])
        r = seq_ss_find(ss, m_ss, start, end)
        if len(r) != 1:
            raise ValueError("More than one segment found")
        bounds = r[0]
        data = (
            row["data"][bounds[1][0] : bounds[1][1]]
            + row["data"][bounds[0][0] : bounds[0][1]]
        )
        all_data.append(data)
    return all_data
