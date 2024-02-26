import freesasa
from biopandas.pdb import PandasPdb
import sys
import pandas as pd


def SASA(name):
    Nt = []
    Nt_num = []
    SASA = []
    ppdb = PandasPdb()
    ppdb.read_pdb(f"{name}")
    ATOM = ppdb.df["ATOM"]
    resname = ATOM["residue_name"]
    atom = ATOM["atom_name"]
    resi_number = ATOM["residue_number"]
    length = len(ATOM.index)
    params = freesasa.Parameters(
        {"algorithm": freesasa.LeeRichards, "probe-radius": 0.5}
    )
    structure = freesasa.Structure(f"{name}")
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


if __name__ == "__main__":
    print(SASA(sys.argv[1]))
