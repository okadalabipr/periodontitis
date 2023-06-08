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

        conlist=['Healthy', 'Perio']
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
        
        for i, condition in enumerate(self.conditions):
            y0 = [0] * V.NUM
            if condition == 'Healthy':
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
            elif condition == 'Perio':
                x[C.ext_up] *= 10
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
        # Gram-neg only
        self.experiments[self.obs_names.index('Bacteria')] = {
            'Healthy': [0.00231, 0.0538, 0.000, 0.00294, 0.00223], 
            'Perio': [0, 1.000, 0.844, 0.767, 0.364], 
        }
        self.error_bars[self.obs_names.index('Bacteria')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.004, 0.0505, 0.0, 0.00508, 0.00386]], 
            'Perio': [sd/np.sqrt(3) for sd in [0, 0.0527, 0.0933, 0.0939, 0.480]], 
        }

        #TLR2
        self.experiments[self.obs_names.index('Toll-like_receptor_2')] = {
            'Healthy': [0.203, 0.230, 0.211, 0, 0.0828], 
            'Perio': [0.362, 0.851, 0.588, 1.000, 0.688], 
        }
        self.error_bars[self.obs_names.index('Toll-like_receptor_2')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.259, 0.0445, 0.170, 0.118, 0.107]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.163, 0.256, 0.519, 0.179, 0.212]], 
        }

        
        self.experiments[self.obs_names.index('il10')] = {
            'Healthy': [0.0486, 0.0521, 0.0698, 0.0, 0.0666], 
            'Perio': [0.106, 0.583, 1.000, 0.396, 0.338], 
        }
        self.error_bars[self.obs_names.index('il10')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
            'Perio': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
        }

        self.experiments[self.obs_names.index('Infiltrating_immune_cells')] = {
            'Healthy': [0.152, 0.502, 0.302, 0.0878, 0.0], 
            'Perio': [0.0613, 0.595, 0.498, 1.000, 0.875], 
        }
        self.error_bars[self.obs_names.index('Infiltrating_immune_cells')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.171, 0.242, 0.369, 0.123, 0.0925]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.0312, 0.330, 0.339, 0.136, 0.210]], 
        }

        self.experiments[self.obs_names.index('Migrating_immune_cells')] = {
            'Healthy': [0.154, 0.574, 0.0697, 0.107, 0.0677], 
            'Perio': [0, 0.589, 0.712, 1.000, 0.858], 
        }
        self.error_bars[self.obs_names.index('Migrating_immune_cells')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.133, 0.223, 0.215, 0.452, 0.159]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.0262, 0.158, 0.224, 0.561, 0.244]], 
        }

        self.experiments[self.obs_names.index('Epithelial_barrier')] = {
            'Healthy': [0.856, 0.713, 0.774, 0.699, 0.835], 
            'Perio': [1.0, 0.454, 0.316, 0.0, 0.631], 
        }
        self.error_bars[self.obs_names.index('Epithelial_barrier')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.321, 0.213, 0.257, 0.292, 0.180]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.0906, 0.0885, 0.270, 0.230, 0.131]], 
        }

        #Itgb3
        self.experiments[self.obs_names.index('Osteoclasts')] = {
            'Healthy': [0, 0.963, 0.340, 0.143, 0.00166], 
            'Perio': [0.420, 1.000, 0.297, 0.740, 0.516], 
        }
        self.error_bars[self.obs_names.index('Osteoclasts')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.403, 0.347, 0.412, 0.515, 0.441]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.816, 0.492, 0.353, 0.386, 0.193]], 
        }

        #CTSCAN
        self.experiments[self.obs_names.index('Bone_integrity')] = {
            'Healthy': [0.917, 0.919, 0.816, 1.000, 1.000], 
            'Perio': [0.840, 0.481, 0.189, 0.000, 0.158], 
        }
        self.error_bars[self.obs_names.index('Bone_integrity')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.00975, 0.00716, 0.00121, 0.00151, 0.00397]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.0259, 0.0560, 0.0879, 0.257, 0.0495]], 
        }
    
        #THIS IS VEGFA
        self.experiments[self.obs_names.index('VEGFA')] = {
            'Healthy': [0.213, 0.260, 0.0574, 0.0153, 0.0], 
            'Perio': [0.309, 0.688, 1.000, 0.413, 0.401], 
        }
        self.error_bars[self.obs_names.index('VEGFA')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
            'Perio': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
        }
        
        self.experiments[self.obs_names.index('arg')] = {
            'Healthy': [0.405, 0.0, 0.317, 0.360, 0.317], 
            'Perio': [0.593, 0.829, 0.674, 0.910, 1.000], 
        }
        self.error_bars[self.obs_names.index('arg')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
            'Perio': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
        }

        self.experiments[self.obs_names.index('pdcd')] = {
            'Healthy': [0.0966, 0.316, 0.0654, 0.0514, 0.0], 
            'Perio': [0.0470, 0.211, 0.957, 1.000, 0.471], 
        }
        self.error_bars[self.obs_names.index('pdcd')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
            'Perio': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
        }

        self.experiments[self.obs_names.index('ppp1r11')] = {
            'Healthy': [0.440, 1.000, 0.0293, 0.0, 0.107], 
            'Perio': [0.140, 0.265, 0.426, 0.363, 0.219], 
        }
        self.error_bars[self.obs_names.index('ppp1r11')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
            'Perio': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
        }
    
        self.experiments[self.obs_names.index('E-selectin')] = {
            'Healthy': [0.0653, 0.0, 0.00884, 0.0224, 0.0297], 
            'Perio': [0.236, 0.740, 0.537, 1.000, 0.283], 
        }
        self.error_bars[self.obs_names.index('E-selectin')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.0309, 0.0268, 0.00948, 0.0107, 0.04]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.125, 0.161, 0.181, 0.859, 0.05]], 
        }

        self.experiments[self.obs_names.index('Osteoblasts')] = {
            'Healthy': [0.654, 0.941, 0.722, 0.694, 0.496], 
            'Perio': [0.896, 0.590, 1.000, 0.0, 0.391], 
        }
        self.error_bars[self.obs_names.index('Osteoblasts')] = {
            'Healthy': [sd/np.sqrt(3) for sd in [0.384, 0.265, 0.0539, 0.0697, 0.171]], 
            'Perio': [sd/np.sqrt(3) for sd in [0.338, 0.244, 0.109, 0.236, 0.135]], 
        }

    @staticmethod
    def get_timepoint(obs_name):
        return [0, 12, 24, 36, 48]
