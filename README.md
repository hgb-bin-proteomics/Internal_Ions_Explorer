# internal_ions

## Setup

- Requirements: **python 3.12+!**
- Packages: `pip install -r requirements.txt`
- Packages (alternative, exact version numbers): `pip install -r env.txt`

## Development Notes

- Implement plot functions in `util/tab*` either in the `example_plot.py` script in there or preferably write your
  own plotting function in a seperate script and throw it in the folder (to avoid merge conflicts).
- Frontend functions should be implemented directly in `tab*.py`.
- Frontend functions that are used accross several tabs go into `utils`.
- Getting data from other tabs: Access via `st.session_state` (examples given) in the `tab*.py` scripts.

## Example Data

- Can be found in `data`
- For plotting we use `fragment_centric.csv`, `spectrum_centric.csv` (and optionally `result.json`).

## Running the App

- Install requirements!
- Run `streamlit run streamlit_app.py`

## Running Fragannot from the Command Line

Use `scripts/fragannot_cli.sh` to create `fragment_centric_<spectrum>.csv`,
`spectrum_centric_<spectrum>.csv`, and `result_<spectrum>.json` without opening
the Streamlit app. The script runs from the repository root internally and uses
`uv sync` before starting the analysis.

Single spectrum/identification pair:

```bash
scripts/fragannot_cli.sh \
  --spectrum data/2022_mix2_rep1.mgf \
  --ident data/2022_mix2_rep1.mzid \
  --filetype mzid \
  --out-dir out \
  --tolerance 0.02 \
  --nterm b \
  --cterm y \
  --charges +1 \
  --losses H2O
```

Batch mode with one identification file per spectrum:

```bash
scripts/fragannot_cli.sh \
  --spectrum /data/tmp/InternalPhospho/OXP*.mgf \
  --ident-template '/path/to/identifications/phospho_{stem}_decoy.tsv' \
  --filetype msms \
  --out-dir /path/to/output
```

Useful flags:

- `--all-identifications`: process all identifications assigned to the same `spectrum_id`;
  duplicate spectrum IDs are exported as `<spectrum_id>_1`, `<spectrum_id>_2`, etc.
- `--check-only`: validate files and parameters without running Fragannot.
- `--deisotope` / `--no-deisotope`: enable or disable deisotoping.
- `--no-json`: write only the two CSV files.
- `--nterm b,c`, `--cterm y,z`, `--charges +1,+2`, `--losses H2O,NH3`: set annotation parameters.

The CLI validates file existence, `.mgf` spectrum input, known `psm_utils`
filetypes, positive tolerance, valid ion names, non-zero integer charges, and
neutral-loss formulas.

## Running the App via Docker

- Running this app via Docker is possible with:
  - `docker build . -f Dockerfile -t internalions`
  - `docker run -p 8501:8501 internalions`

## Server

- [Online](https://computproteomics.bmb.sdu.dk/app_direct/internal_ions/)
