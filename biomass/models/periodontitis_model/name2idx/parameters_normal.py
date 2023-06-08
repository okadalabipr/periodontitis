from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    'k1',
    'cxcl3',
    'tnf',
    'defb3',
    'lyz2',
    'ccl2',
    'csf1',
    'ccl9',
    'ctsk',
    'ext_sp',
    'ext_up',
    'il10_a',
    'il10_b',
    'inh_a',
    'inh_b',
    'mmp3',
    'gpx',
    'Ligand',
    'EGF',
    'HRG',
    'no_ligand',
    'arg_1',
    'arg_2',
    'pdcd1',
    'pdcd2',
    'ppp_1',
    'ppp_2',
    'tlrsele',
    'tek_1',
    'tek_2',
    'ostdeg',
    'gja',
    'rankl',
    'alpl',
    'dpdcd',
    'dil10',
    'dtek',
    'dppp1r11',
    'darg1',
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
