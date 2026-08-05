#!/bin/bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export UV_CACHE_DIR=.uv-cache
uv sync

uv run python - <<'PY'
from io import BytesIO
from pathlib import Path
import json

from psm_utils.io import read_file

from internal_ions.fragannot.fragannot_call import fragannot_call
from internal_ions.util.converter import JSONConverter
from internal_ions.util.spectrumio import SpectrumFile


class LocalUpload(BytesIO):
    def __init__(self, path):
        path = Path(path)
        super().__init__(path.read_bytes())
        self.name = str(path)


spectrum_file = "SPECTRUMFILE"
ident_file = "IDENTFILE"
out_dir = Path("out")
out_dir.mkdir(exist_ok=True)

spectra = SpectrumFile(LocalUpload(spectrum_file))
psms = read_file(ident_file, filetype="infer")

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
fragment_df.to_csv(out_dir / "fragment_centric.csv", index=False)
spectrum_df.to_csv(out_dir / "spectrum_centric.csv", index=False)
(out_dir / "result.json").write_text(json.dumps(result), encoding="utf-8")

print(f"Wrote {out_dir / 'fragment_centric.csv'}")
print(f"Wrote {out_dir / 'spectrum_centric.csv'}")
PY
