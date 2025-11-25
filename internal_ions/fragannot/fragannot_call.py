from .fragannot_numba import FragannotNumba

from ..util.spectrumio import SpectrumFile
from psm_utils.psm_list import PSMList


def fragannot_call(spectrum_file: SpectrumFile,
                   psms: PSMList,
                   tolerance: float,
                   nterm_fragment_types: list[str],
                   cterm_fragment_types: list[str],
                   charges: list[str],
                   losses: list[str],
                   deisotope: bool,
                   verbose: bool = False) -> dict:
    frag = FragannotNumba()
    fragannot_dict = frag.fragment_annotation(psms,
                                              spectrum_file,
                                              tolerance,
                                              nterm_fragment_types,
                                              cterm_fragment_types,
                                              charges,
                                              losses,
                                              deisotope,
                                              write_file=False)

    return fragannot_dict
