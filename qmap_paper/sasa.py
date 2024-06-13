import freesasa
from biopandas.pdb import PandasPdb
import pandas as pd
from typing import List


def compute_solvent_accessability(pdb_path: str) -> pd.DataFrame:
    """
    Computes the solvent accessibility of atoms in a protein structure.

    Args:
        pdb_path: The path to the PDB file.

    Returns:
        A pandas DataFrame containing the solvent accessibility information for each atom.

    Raises:
        FileNotFoundError: If the PDB file specified by pdb_path does not exist.
    """
    Nt: List[str] = []
    Nt_num: List[int] = []
    SASA: List[float] = []
    ppdb = PandasPdb()
    ppdb.read_pdb(f"{pdb_path}")
    ATOM = ppdb.df["ATOM"]
    resname = ATOM["residue_name"]
    atom = ATOM["atom_name"]
    resi_number = ATOM["residue_number"]
    length = len(ATOM.index)
    params = freesasa.Parameters({"algorithm": freesasa.LeeRichards, "probe-radius": 1})
    structure = freesasa.Structure(f"{pdb_path}")
    result = freesasa.calc(structure, params)
    for i in range(length):
        if atom[i] == "N1":
            if resname[i] == "A":
                selection = freesasa.selectArea(
                    (
                        f"n1,(name N1) and (resn A) and (resi {resi_number[i]}) ",
                        f"n3, (name N3) and (resn C) and (resi {resi_number[i]})",
                    ),
                    structure,
                    result,
                )
                sel1 = selection["n1"]
                Nt.append(resname[i])
                Nt_num.append(resi_number[i])
                SASA.append(sel1)

            elif resname[i] == "C":
                selection = freesasa.selectArea(
                    (
                        f"n1,(name N1) and (resn A) and (resi {resi_number[i]}) ",
                        f"n3, (name N3) and (resn C) and (resi {resi_number[i]})",
                    ),
                    structure,
                    result,
                )
                sel2 = selection["n3"]

                Nt.append(resname[i])
                Nt_num.append(resi_number[i])
                SASA.append(sel2)

            else:
                Nt.append(resname[i])
                Nt_num.append(resi_number[i])
                SASA.append(0)

    d = {"Nt": Nt, "Nt_num": Nt_num, "SASA": SASA}
    df = pd.DataFrame(d)
    return df
