import numpy as np

from biomass.dynamics.solver import get_steady_state, solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of self.obs_names.

    t : range
        Simulation time span.

    conditions : list of strings
        Experimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If :obj:`None`, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in ``sim.conditions`` will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    
    
    
    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names = [
            'Bacteria',
            'Toll-like_receptor_2',
            'Infiltrating_immune_cells',
            'Migrating_immune_cells',
            'Epithelial_barrier',
            'Osteoclasts',
            'Bone_integrity',
            'il10',
            'arg',
            'pdcd',
            'ppp1r11',
            'E-selectin',
            'VEGFA',
            'Osteoblasts'
        ]
        self.t: range = range(49)

        conlist=['RANKL_OE', 'RANKL_KO']
        self.conditions: list = conlist
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t), len(self.conditions))
        )
        self.normalization: dict = {}
        #for observable in self.obs_names:
        #    self.normalization[observable] = {"timepoint": None, "condition": []}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)
        
    def simulate(self, x, y0, _perturbation=None):
        if _perturbation is not None:
            self.perturbation = _perturbation
        # get steady state
        #x[C.Ligand] = x[C.no_ligand]  # No ligand
        #y0 = get_steady_state(self.diffeq, y0, tuple(x))
        if not y0:
            return False
        # add ligand
        x[C.ext_up] *= 10
        for i, condition in enumerate(self.conditions):
            y0 = [0] * V.NUM
            if condition == 'RANKL_OE':
                y0[V.ext] = 1
                y0[V.bact]=0.00231
                y0[V.tlr2]=0.203
                y0[V.infiltrate]=0.152
                y0[V.migrate]=0.154
                y0[V.epi]=0.856
                y0[V.ost]=0
                y0[V.boneres]=0.917
                y0[V.il10]=0.0486
                y0[V.arg]=0.405
                y0[V.pdcd]=0.0966
                y0[V.ppp1r11]=0.440
                y0[V.sele] = 0.0653
                y0[V.tek] = 0.213
                y0[V.ob] = 0.654

                x[C.rankl] *= 10
            elif condition == 'RANKL_KO':
                
                y0[V.ext] = 1
                y0[V.bact]=0
                y0[V.tlr2]=0.362
                y0[V.infiltrate]=0.0613
                y0[V.migrate]=0
                y0[V.epi]=1.0
                y0[V.ost]=0.420
                y0[V.boneres]=0.840
                y0[V.il10]=0.106
                y0[V.arg]=0.593
                y0[V.pdcd]=0.0470
                y0[V.ppp1r11]=0.140
                y0[V.sele] = 0.236
                y0[V.tek] = 0.309
                y0[V.ob] = 0.896

                x[C.rankl] = 0

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index('Bacteria'), :, i] = (
                    sol.y[V.bact,:]
                )
                self.simulations[self.obs_names.index('Toll-like_receptor_2'), :, i] = (
                    sol.y[V.tlr2,:]
                )
                self.simulations[self.obs_names.index('Infiltrating_immune_cells'), :, i] = (
                    sol.y[V.infiltrate,:]
                )
                self.simulations[self.obs_names.index('Migrating_immune_cells'), :, i] = (
                    sol.y[V.migrate,:]
                )

                self.simulations[self.obs_names.index('Epithelial_barrier'), :, i] = (
                    sol.y[V.epi,:]
                )
                self.simulations[self.obs_names.index('Osteoclasts'), :, i] = (
                    sol.y[V.ost,:]
                )

                self.simulations[self.obs_names.index('Bone_integrity'), :, i] = (
                    sol.y[V.boneres,:]
                )

                self.simulations[self.obs_names.index('il10'), :, i] = (
                    sol.y[V.il10,:]
                )
       
                self.simulations[self.obs_names.index('arg'), :, i] = (
                    sol.y[V.arg,:]
                )

                self.simulations[self.obs_names.index('pdcd'), :, i] = (
                    sol.y[V.pdcd,:]
                )

                self.simulations[self.obs_names.index('ppp1r11'), :, i] = (
                    sol.y[V.ppp1r11,:]
                )

                self.simulations[self.obs_names.index('E-selectin'), :, i] = (
                    sol.y[V.sele,:]
                )
   
                self.simulations[self.obs_names.index('VEGFA'), :, i] = (
                    sol.y[V.tek,:]
                )
                
                self.simulations[self.obs_names.index('Osteoblasts'), :, i] = (
                    sol.y[V.ob,:]
                )
    def set_data(self):
        return
    @staticmethod
    def get_timepoint(obs_name):
        return [0, 12, 24, 36, 48]
