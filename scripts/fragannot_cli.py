#!/usr/bin/env -S uv run python
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter
from io import BytesIO
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/internal_ions_mpl_cache")

from psm_utils.io import FILETYPES
from pyteomics import mass

from internal_ions.util import constants


class LocalUpload(BytesIO):
    def __init__(self, path: Path):
        super().__init__(path.read_bytes())
        self.name = str(path)


class SpectrumAliases:
    def __init__(self, spectra: SpectrumFile, aliases: dict[str, str]):
        self.spectra = spectra
        self.aliases = aliases
        self.name = spectra.name

    def get_by_id(self, spectrum_id: str):
        return self.spectra.get_by_id(self.aliases.get(spectrum_id, spectrum_id))


def positive_float(value: str) -> float:
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return number


def comma_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_charges(value: str):
    if value == "auto":
        return "auto"
    charges = comma_list(value)
    if not charges:
        raise argparse.ArgumentTypeError("must contain at least one charge")
    for charge in charges:
        if not re.fullmatch(r"[+-]?\d+", charge) or int(charge) == 0:
            raise argparse.ArgumentTypeError(f"invalid charge: {charge}")
    return charges


def parse_losses(value: str) -> list[str]:
    if value.lower() in {"none", "no", "false"}:
        return [""]
    losses = [item.strip() for item in value.split(",")]
    for loss in losses:
        if loss:
            try:
                mass.calculate_mass(loss, absolute=True)
            except Exception as exc:
                raise argparse.ArgumentTypeError(f"invalid neutral loss formula: {loss}") from exc
    if "" not in losses:
        losses.append("")
    return losses


def ions(value: str, direction: str) -> list[str]:
    selected = comma_list(value)
    choices = sorted(k for k, v in constants.ion_direction.items() if v == direction)
    bad = [ion for ion in selected if ion not in choices]
    if bad:
        raise argparse.ArgumentTypeError(f"invalid {direction} ion(s): {', '.join(bad)}")
    if not selected:
        raise argparse.ArgumentTypeError(f"must contain at least one {direction} ion")
    return selected


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="fragannot_cli.py", description="Run Fragannot from the command line.")
    p.add_argument("--spectrum", nargs="+", required=True, type=Path, help="MGF spectrum file(s).")
    id_group = p.add_mutually_exclusive_group(required=True)
    id_group.add_argument("--ident", type=Path, help="Identification file. Only valid with one --spectrum.")
    id_group.add_argument(
        "--ident-template",
        help="Identification template for multiple spectra. Supports {stem}, {name}, {dir}, {spectrum}.",
    )
    p.add_argument("--filetype", default="infer", help="Identification file type: infer or a psm_utils type.")
    p.add_argument("--out-dir", type=Path, help="Output directory. Defaults to each spectrum file directory.")
    p.add_argument("--tolerance", type=positive_float, default=0.02, help="Fragment mass tolerance in Da.")
    p.add_argument("--nterm", type=lambda x: ions(x, "n-term"), default="b", help="Comma-delimited n-term ions.")
    p.add_argument("--cterm", type=lambda x: ions(x, "c-term"), default="y", help="Comma-delimited c-term ions.")
    p.add_argument("--charges", type=parse_charges, default="+1", help="Comma-delimited charges, e.g. +1,+2.")
    p.add_argument("--losses", type=parse_losses, default="H2O,", help="Comma-delimited neutral losses. Use none for no losses.")
    p.add_argument("--deisotope", action=argparse.BooleanOptionalAction, default=True, help="Enable/disable deisotoping.")
    p.add_argument(
        "--all-identifications",
        "--all-identifications-per-spectrum",
        action="store_true",
        help="Process all identifications assigned to the same spectrum_id.",
    )
    p.add_argument("--no-json", action="store_true", help="Do not write result_<stem>.json.")
    p.add_argument("--check-only", action="store_true", help="Validate inputs and parameters without running Fragannot.")
    return p


def ident_for(args, spectrum: Path) -> Path:
    if args.ident is not None:
        return args.ident
    return Path(args.ident_template.format(
        stem=spectrum.stem,
        name=spectrum.name,
        dir=str(spectrum.parent),
        spectrum=str(spectrum),
    ))


def validate(args) -> list[tuple[Path, Path, Path]]:
    if args.filetype != "infer" and args.filetype not in FILETYPES:
        raise SystemExit(f"Unknown --filetype {args.filetype!r}. Valid: infer, {', '.join(sorted(FILETYPES))}")
    if args.ident is not None and len(args.spectrum) != 1:
        raise SystemExit("--ident can only be used with one --spectrum; use --ident-template for multiple spectra.")

    jobs = []
    for spectrum in args.spectrum:
        if not spectrum.is_file():
            raise SystemExit(f"Spectrum file not found: {spectrum}")
        if spectrum.suffix.lower() != ".mgf":
            raise SystemExit(f"Only .mgf spectrum files are supported: {spectrum}")
        ident = ident_for(args, spectrum)
        if not ident.is_file():
            raise SystemExit(f"Identification file not found for {spectrum}: {ident}")
        out_dir = args.out_dir or spectrum.parent
        jobs.append((spectrum, ident, out_dir))
    return jobs


def uniquify_duplicate_spectrum_ids(psms):
    counts = Counter(psm.spectrum_id for psm in psms)
    seen = Counter()
    aliases = {}
    output_ids = {}

    for psm in psms:
        original = psm.spectrum_id
        if counts[original] == 1:
            continue
        seen[original] += 1
        unique = f"{original}__psm_{seen[original]}"
        psm.spectrum_id = unique
        aliases[unique] = original
        output_ids[unique] = f"{original}_{seen[original]}"

    return aliases, output_ids, counts


def run_job(args, spectrum_file: Path, ident_file: Path, out_dir: Path) -> None:
    from psm_utils.io import read_file

    from internal_ions.fragannot.fragannot_call import fragannot_call
    from internal_ions.util.converter import JSONConverter
    from internal_ions.util.spectrumio import SpectrumFile

    out_dir.mkdir(parents=True, exist_ok=True)
    spectra = SpectrumFile(LocalUpload(spectrum_file))
    psms = read_file(str(ident_file), filetype=args.filetype)
    aliases, output_ids, original_counts = {}, {}, Counter()

    spectra_arg = spectra
    if args.all_identifications:
        aliases, output_ids, original_counts = uniquify_duplicate_spectrum_ids(psms)
        spectra_arg = SpectrumAliases(spectra, aliases)

    result = fragannot_call(
        spectra_arg,
        psms,
        tolerance=args.tolerance,
        nterm_fragment_types=args.nterm,
        cterm_fragment_types=args.cterm,
        charges=args.charges,
        losses=args.losses,
        deisotope=args.deisotope,
    )

    if args.all_identifications:
        for entry in result.values():
            original = aliases.get(entry["spectrum_id"], entry["spectrum_id"])
            entry["spectrum_id"] = output_ids.get(entry["spectrum_id"], entry["spectrum_id"])
            entry["nr_idents_with_same_rank"] = original_counts[original]

    fragment_df, spectrum_df = JSONConverter().to_dataframes(result)
    stem = spectrum_file.stem
    fragment_path = out_dir / f"fragment_centric_{stem}.csv"
    spectrum_path = out_dir / f"spectrum_centric_{stem}.csv"
    fragment_df.to_csv(fragment_path, index=False)
    spectrum_df.to_csv(spectrum_path, index=False)
    if not args.no_json:
        (out_dir / f"result_{stem}.json").write_text(json.dumps(result), encoding="utf-8")
    print(f"Wrote {fragment_path}")
    print(f"Wrote {spectrum_path}")


def main(argv: list[str]) -> int:
    args = parser().parse_args(argv)
    jobs = validate(args)
    if args.check_only:
        for spectrum, ident, out_dir in jobs:
            print(f"OK: {spectrum} + {ident} -> {out_dir}")
        return 0
    for spectrum, ident, out_dir in jobs:
        run_job(args, spectrum, ident, out_dir)
    return 0


raise SystemExit(main(sys.argv[1:]))
