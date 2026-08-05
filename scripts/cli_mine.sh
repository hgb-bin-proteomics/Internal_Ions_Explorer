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
  ident_file="$(dirname "$file")/Spike_${base}.tsv"
  out_dir="$(dirname "$file")"

  uv run python - "$file" "$ident_file" "$out_dir" "$base" <<'PY'
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


spectrum_file, ident_file, out_dir, base = sys.argv[1:]
out_dir = Path(out_dir)
out_dir.mkdir(exist_ok=True)

spectra = SpectrumFile(LocalUpload(spectrum_file))
psms = read_file(ident_file, filetype="msms")

result = fragannot_call(
    spectra,
    psms,
    tolerance=0.02,
    nterm_fragment_types=["b"],
    cterm_fragment_types=["y"],
    charges=["+1"],
    losses=["H2O", ""],
    deisotope=True,
)

fragment_df, spectrum_df = JSONConverter().to_dataframes(result)
fragment_df.to_csv(out_dir / f"fragment_centric_{base}.csv", index=False)
spectrum_df.to_csv(out_dir / f"spectrum_centric_{base}.csv", index=False)
(out_dir / f"result_{base}.json").write_text(json.dumps(result), encoding="utf-8")

print(f"Wrote {out_dir / f'fragment_centric_{base}.csv'}")
print(f"Wrote {out_dir / f'spectrum_centric_{base}.csv'}")
PY
done
