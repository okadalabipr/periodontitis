from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    'ext_up',
    'ext_sp',
    'lyz2',
    'defb3',
    'ctsk',
    'alpl',
    'k1',
    'ppp_1',
    'ppp_2',
    'dppp1r11',
    'tnf',
    'cxcl3',
    'pdcd1',
    'pdcd2',
    'arg_1',
    'arg_2',
    'dpdcd',
    'darg1',
    'il10_a',
    'inh_a',
    'il10_b',
    'inh_b',
    'dil10',
    'tlrsele',
    'tek_1',
    'tek_2',
    'dtek',
    'ccl2',
    'csf1',
    'gja',
    'ostdeg',
    'ccl9',
    'rankl',
    'mmp3',
    'gpx',
    'Ligand',
    'EGF',
    'HRG',
    'no_ligand',
]
for idx, name in enumerate(NAMES):
    exec("{} = {:d}".format(name, idx))

NUM: int = len(NAMES)

Parameters = make_dataclass(
    cls_name="Parameters",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

C = Parameters(**name2idx)

del name2idx
