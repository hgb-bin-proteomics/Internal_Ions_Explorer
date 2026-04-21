from pyteomics import mass
import logging

from psm_utils.io import FILETYPES
SUPPORTED_FILETYPES = list(FILETYPES)

REPO_OWNER = "hgb-bin-proteomics"
REPO_NAME = "Internal_Ions_Explorer"
DIV_COLOR = "rainbow"
FRAGANNOT_ION_NAMES = ["a", "b", "c", "cdot", "c-1", "c+1", "x", "y", "z", "zdot", "z+1", "z+2", "z+3"]

logger = logging.getLogger(__name__)


class HashableComp(mass.Composition):
    def __hash__(self):
        return hash(tuple(sorted(self.items())))


ion_comp = {key: HashableComp(val) for (key, val) in mass.std_ion_comp.items()}

for k in list(ion_comp):
    if k.endswith('dot'):
        ion_comp[k.replace('-', '')] = ion_comp.pop(k)
ion_comp['t'] = HashableComp({})

ion_cap_delta_mass = {
    name: mass.calculate_mass(composition=comp, absolute=False) for (name, comp) in ion_comp.items()
}

ion_direction = {
    "a": "n-term",
    "b": "n-term",
    "x": "c-term",
    "y": "c-term",
    "cdot": "n-term",
    "c": "n-term",
    "c-1": "n-term",
    "c+1": "n-term",
    "zdot": "c-term",
    "z": "c-term",
    "z+1": "c-term",
    "z+2": "c-term",
    "z+3": "c-term",
}


def _identical_internal_ions() -> dict[tuple[str, str], tuple[str, str]]:
    nterm = [key for key, val in ion_direction.items() if val == "n-term"]
    cterm = [key for key, val in ion_direction.items() if val == "c-term"]
    internal_ion_comps = {}
    ion_mapping = {}
    for n in nterm:
        for c in cterm:
            internal_ion_comps.setdefault((ion_comp[n] + ion_comp[c]), []).append((n, c))
    for values in internal_ion_comps.values():
        values.sort(key=lambda x: len(x[0]) + len(x[1]))
        for v in values:
            ion_mapping[v] = values[0]
    return ion_mapping

IDENTICAL_INTERNAL_IONS = _identical_internal_ions()

for key, val in IDENTICAL_INTERNAL_IONS.items():
    if key != val:
        logger.debug(f"Identical internal ions: {key} -> {val}")



# ---------------------------------------------------------------------------- #
#                               For visualization                              #
# ---------------------------------------------------------------------------- #


colors = [
    "#9b2226",
    "#005f73",
    "#ee9b00",
    "#0a9396",
    "#94d2bd",
    "#ca6702",
    "#e9d8a6",
    "#bb3e03",
    "#001219",
    "#006BA6",
    "#35A7FF",
    "#EFA8B8",
    "#BFACC8",
    "#476A6F",
    "#7067CF",
    "#364156",
    "#98FB98",
    "#8A2BE2",
    "#35682D",
    "#252850",
    "#7E7B52",
]

colors = colors + (["#808080"] * 1000)
