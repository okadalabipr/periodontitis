from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    'ext_up',
    'ext_sp',
    'k1',
    'ppp_1',
    'ppp_2',
    'defb3',
    'tnf',
    'tlrsele',
    'tek_1',
    'tek_2',
    'lyz2',
    'cxcl3',
    'pdcd1',
    'pdcd2',
    'il10_a',
    'il10_b',
    'inh_a',
    'arg_2',
    'arg_1',
    'ccl2',
    'inh_b',
    'csf1',
    'ccl9',
    'ostdeg',
    'gja',
    'rankl',
    'ctsk',
    'alpl',
    'dpdcd',
    'dil10',
    'dtek',
    'dppp1r11',
    'darg1',
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
