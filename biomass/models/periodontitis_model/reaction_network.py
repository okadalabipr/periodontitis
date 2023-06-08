from typing import Dict, List

from biomass.analysis.reaction import is_duplicate


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (target='reaction').
    """

    def __init__(self) -> None:
        self.reactions: Dict[str, List[int]] = {
            "Bact_up": [11],
            "Bact_down": [4,5],
            "Tlr_up": [1],
            "Tlr_down": [21],
            "infiltrate_up": [2],
            "infiltrate_down": [13,17,19],
            "migrate_up": [6,8],
            "migrate_down": [15],
            "epi_down": [3],
            "ost_up": [7,26],
            "ost_down": [25],
            "il10_up": [12,14],
            "il10_down": [30],
            "arg_up": [16],
            "boneres_up": [28],
            "boneres_down": [9],
            "pdcd_up": [18],
            "pdcd_down": [29],
            "ppp1r11_up": [20],
            "ppp1r11_down": [32],
            "sele_up": [22],
            "sele_down": [24],
            "tek_up": [23],
            "tek_down": [31],
            "ob_down": [27],
        }

    def group(self):
        """
        Group reactions according to biological processes
        """
        biological_processes = []
        for process, indices in self.reactions.items():
            if not isinstance(indices, list):
                raise TypeError("Use list for reaction indices in {}".format(process))
            biological_processes.append(indices)

        if not is_duplicate(self.reactions, biological_processes):
            return biological_processes
