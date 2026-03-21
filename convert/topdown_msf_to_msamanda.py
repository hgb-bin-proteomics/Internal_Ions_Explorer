#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas"
# ]
# ///

import sys
import sqlite3
import pandas as pd


def read_msf(msf_file: str) -> dict[str, pd.DataFrame]:
    conn = sqlite3.connect(msf_file)
    proteoforms = pd.read_sql_query(
        "SELECT * FROM TargetProteoformSpectrumMatchs", conn
    )
    modifications = pd.read_sql_query("SELECT * FROM FoundModifications", conn)
    conn.close()
    return {"proteoforms": proteoforms, "modifications": modifications}


def to_msamanda(msf_file: str) -> int:
    msf = read_msf(msf_file)
    ##
    mods = dict()
    for i, row in msf["modifications"].iterrows():
        name = str(row["Name"]).strip()
        abbreviation = str(row["Abbreviation"]).strip()
        mass = float(row["DeltaMonoisotopicMass"])
        if name not in mods:
            mods[name] = mass
        if abbreviation not in mods:
            mods[abbreviation] = mass
    ##
    title = list()
    sequence = list()
    modifications = list()
    protein_accessions = list()
    score = list()
    mz = list()
    charge = list()
    rt = list()
    filename = list()
    id = list()
    for i, row in msf["proteoforms"].iterrows():
        title.append(
            f"controllerType=0 controllerNumber=1 scan={row['FragmentationScans']}"
        )
        sequence.append(row["ModifiedSequence"])
        modifications_str = ""
        if not pd.isna(row["Modifications"]):
            for mod in str(row["Modifications"]).split(";"):
                loc = mod.split("(")[0].strip()
                name = ")".join("(".join(mod.split("(")[1:]).split(")")[:-1]).strip()
                mass = mods[name]
                modifications_str += f"{loc}({name}|{mass}|variable);"
        modifications.append(modifications_str.rstrip(";"))
        if pd.isna(row["ParentProteinAccessions"]):
            protein_accessions.append("sp|UNKNOWN")
        else:
            protein_accessions.append(row["ParentProteinAccessions"])
        score.append(row["CScore"])
        mz.append(row["MassOverCharge"])
        charge.append(row["Charge"])
        rt.append(row["RetentionTime"])
        filename.append(row["SpectrumFileName"])
        id.append(i)
    amanda = pd.DataFrame(
        {
            "Title": title,
            "Sequence": sequence,
            "Modifications": modifications,
            "Protein Accessions": protein_accessions,
            "Amanda Score": score,
            "m/z": mz,
            "Charge": charge,
            "RT": rt,
            "Filename": filename,
            "Id": id,
        }
    )
    amanda.to_csv(f"{msf_file}.csv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    exit(to_msamanda(sys.argv[1]))
