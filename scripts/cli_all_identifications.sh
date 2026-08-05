#!/bin/bash
set -euo pipefail
shopt -s nullglob
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export UV_CACHE_DIR=.uv-cache
uv sync

files=("$@")
if [ "$#" -eq 0 ]; then
  files=(/data/tmp/InternalPhospho/OXP*.mgf)
fi

for file in "${files[@]}"
do
  base="$(basename "${file%.*}")"
  ident_file="/home/veit/devel/Bioinformatics/InternalIonsExploration/phospho_${base}_decoy.tsv"
  out_dir="/home/veit/devel/Bioinformatics/InternalIonsExploration/"

  uv run python - "$file" "$ident_file" "$out_dir" "$base" <<'PY'
from collections import Counter
from io import BytesIO
from pathlib import Path
import json
import sys

from psm_utils.io import read_file

from internal_ions.fragannot.fragannot_call import fragannot_call
from internal_ions.util.converter import JSONConverter
from internal_ions.util.spectrumio import SpectrumFile


class LocalUpload(BytesIO):
    def __init__(self, path):
        path = Path(path)
        super().__init__(path.read_bytes())
        self.name = str(path)


class SpectrumAliases:
    def __init__(self, spectra, aliases):
        self.spectra = spectra
        self.aliases = aliases
        self.name = spectra.name

    def get_by_id(self, spectrum_id):
        return self.spectra.get_by_id(self.aliases.get(spectrum_id, spectrum_id))


def uniquify_duplicate_spectrum_ids(psms):
    counts = Counter(psm.spectrum_id for psm in psms)
    seen = Counter()
    aliases = {}
    originals = {}

    for psm in psms:
        original = psm.spectrum_id
        if counts[original] == 1:
            continue
        seen[original] += 1
        unique = f"{original}__psm_{seen[original]}"
        psm.spectrum_id = unique
        aliases[unique] = original
        originals[unique] = original

    return aliases, originals, counts


spectrum_file, ident_file, out_dir, base = sys.argv[1:]
out_dir = Path(out_dir)
out_dir.mkdir(exist_ok=True)

spectra = SpectrumFile(LocalUpload(spectrum_file))
psms = read_file(ident_file, filetype="msms")
aliases, originals, original_counts = uniquify_duplicate_spectrum_ids(psms)

result = fragannot_call(
    SpectrumAliases(spectra, aliases),
    psms,
    tolerance=0.02,
    nterm_fragment_types=["b"],
    cterm_fragment_types=["y"],
    charges=["+1"],
    losses=["H2O", ""],
    deisotope=True,
)

for key, entry in result.items():
    original = originals.get(entry["spectrum_id"], entry["spectrum_id"])
    entry["spectrum_id"] = original
    entry["nr_idents_with_same_rank"] = original_counts[original]

fragment_df, spectrum_df = JSONConverter().to_dataframes(result)
fragment_df.to_csv(out_dir / f"fragment_centric_{base}.csv", index=False)
spectrum_df.to_csv(out_dir / f"spectrum_centric_{base}.csv", index=False)
(out_dir / f"result_{base}.json").write_text(json.dumps(result), encoding="utf-8")

print(f"Wrote {out_dir / f'fragment_centric_{base}.csv'}")
print(f"Wrote {out_dir / f'spectrum_centric_{base}.csv'}")
PY
done
