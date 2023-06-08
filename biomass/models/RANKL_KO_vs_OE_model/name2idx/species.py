from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    'bact',
    'tlr2',
    'infiltrate',
    'migrate',
    'epi',
    'ost',
    'boneres',
    'ext',
    'il10',
    'arg',
    'pdcd',
    'ppp1r11',
    'sele',
    'tek',
    'ob'
]
for idx, name in enumerate(NAMES):
    exec(f"{name} = {idx:d}")

NUM: int = len(NAMES)


Species = make_dataclass(
    cls_name="Species",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

V = Species(**name2idx)

del name2idx
