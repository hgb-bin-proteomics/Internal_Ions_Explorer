#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas",
#   "psm-utils"
# ]
# ///

import sqlite3
import argparse
import pandas as pd
from psm_utils.io import read_file

__version = "1.0.0"
__data = "2026-03-22"


def read_msf(msf_file: str) -> dict[str, pd.DataFrame]:
    conn = sqlite3.connect(msf_file)
    proteoforms = pd.read_sql_query(
        "SELECT * FROM TargetProteoformSpectrumMatchs", conn
    )
    modifications = pd.read_sql_query("SELECT * FROM FoundModifications", conn)
    conn.close()
    return {"proteoforms": proteoforms, "modifications": modifications}


def to_msamanda(msf_file: str, verify: bool = True) -> int:
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
    if verify:
        psms = read_file(f"{msf_file}.csv", filetype="msamanda")
        print(f"Successfully read {len(psms)} PSMs from file!")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="topdown_msf_to_msamanda.py",
        description="Converts a MSF file with top-down proteomics results to MS Amanda format.",
        epilog="(c) Bioinformatics Research Group, FH OÖ Campus Hagenberg, 2026",
    )
    parser.add_argument(
        dest="msf",
        help="MSF file to convert to MS Amanda format.",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--check",
        dest="verify",
        action="store_true",
        help="Check PSMs with psm_utils.",
    )
    parser.add_argument("--version", action="version", version=__version)
    args = parser.parse_args(argv)

    return to_msamanda(args.msf, args.verify)


if __name__ == "__main__":
    exit(main())
