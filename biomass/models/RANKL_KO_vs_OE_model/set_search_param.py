import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .set_model import initial_values, param_values


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize."""

    def __init__(self):
        # parameters
        self.idx_params = [
            C.k1,
            C.lyz2,
            C.defb3,
            C.cxcl3,
            C.tnf,
            C.ccl2,
            C.ccl9,
            C.csf1,
            C.ctsk,
            C.ext_sp,
            C.ext_up,
            C.il10_a,
            C.il10_b,
            C.mmp3,
            C.inh_a,
            C.inh_b,
            C.arg_1,
            C.arg_2,
            C.pdcd1,
            C.pdcd2,
            C.ppp_1,
            C.ppp_2,
            C.tlrsele,
            C.tek_1,
            C.tek_2,
            C.ostdeg,
            C.gja,
            C.rankl,
            C.alpl,
            C.dpdcd,
            C.dil10,
            C.dtek,
            C.dppp1r11,
            C.darg1
        ]

        # initial values
        self.idx_initials = [
            # V.(specie)
        ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = initialize_search_param(
            parameters=C.NAMES,
            species=V.NAMES,
            param_values=x,
            initial_values=y0,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound
            search_rgn[1, j] = search_param[i] * 10.0  # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 0.5  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 2.0  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound, upper_bound]

        

        search_rgn = convert_scale(
            region=search_rgn,
            parameters=C.NAMES,
            species=V.NAMES,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        # constraints ----------------------------------------------------------
        
        # ----------------------------------------------------------------------

        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10 ** (indiv_gene * (search_rgn[1, :] - search_rgn[0, :]) + search_rgn[0, :])

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (np.log10(indiv) - search_rgn[0, :]) / (search_rgn[1, :] - search_rgn[0, :])

        return indiv_gene
