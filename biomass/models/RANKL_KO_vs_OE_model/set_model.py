from .name2idx import C, V


class DifferentialEquation(object):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Kinetic equations"""
        # v : flux vector
        v = {}
        v[1] = x[C.k1] * y[V.bact]
        v[2] = x[C.cxcl3] * y[V.tlr2]
        v[3] = x[C.tnf] * (y[V.tlr2]/(1+y[V.tlr2])) * (y[V.epi]/((1+y[V.epi])**(2)))
        v[4] = x[C.lyz2] * (y[V.infiltrate]/(1+y[V.infiltrate])) * (y[V.bact]/((1+y[V.bact])**(2)))
        v[5] = x[C.defb3] * (y[V.epi]/(1+y[V.epi])) * (y[V.bact]/((1+y[V.bact])**(2)))
        v[6] = x[C.ccl2] * y[V.tlr2]
        v[7] = x[C.ccl9] * y[V.migrate]
        v[8] = x[C.csf1] * y[V.sele]
        v[9] = x[C.ctsk] * (y[V.ost]/(1+y[V.ost])) * (y[V.boneres]/((1+y[V.boneres])**(2)))
        v[10] = x[C.ext_sp] * y[V.ext]
        v[11] = x[C.ext_up] * y[V.ext]
        v[12] = x[C.il10_a] * y[V.infiltrate]
        v[13] = x[C.il10_b] * (y[V.il10]/(1+y[V.il10])) * (y[V.infiltrate]/((1+y[V.infiltrate])**(2)))
        v[14] = x[C.inh_a] * y[V.migrate]
        v[15] = x[C.inh_b] * (y[V.il10]/(1+y[V.il10]))* (y[V.migrate]/((1+y[V.migrate])**(2)))
        v[16] = x[C.arg_1] * y[V.migrate]
        v[17] = x[C.arg_2] * y[V.arg] * (y[V.infiltrate]/((1+y[V.infiltrate])**(2)))
        v[18] = x[C.pdcd1] * y[V.infiltrate]
        v[19] = x[C.pdcd2] * (y[V.pdcd]/(1+y[V.pdcd])) * (y[V.infiltrate]/((1+y[V.infiltrate])**(2)))
        v[20] = x[C.ppp_1] * y[V.tlr2]
        v[21] = x[C.ppp_2] * (y[V.ppp1r11]/(1+y[V.ppp1r11])) * (y[V.tlr2]/((1+y[V.tlr2])**(2)))
        v[22] = x[C.tlrsele] * y[V.tlr2]
        v[23] = x[C.tek_1] * y[V.sele]
        v[24] = x[C.tek_2] * (y[V.tek]/(1+y[V.tek])) * (y[V.sele]/((1+y[V.sele])**(2)))
        v[25] = x[C.ostdeg] * y[V.ost]
        v[26] = x[C.rankl] * y[V.ob]
        v[27] = x[C.gja] * (y[V.ost]/(1+y[V.ost])) * (y[V.ob]/((1+y[V.ob])**(2)))
        v[28] = x[C.alpl] * y[V.ob]
        v[29] = x[C.dpdcd] * y[V.pdcd]
        v[30] = x[C.dil10] * y[V.il10]
        v[31] = x[C.dtek] * y[V.tek]
        v[32] = x[C.dppp1r11] * y[V.ppp1r11]
        v[33] = x[C.darg1] * y[V.arg]

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        dydt[V.ext] = - v[10]
        dydt[V.bact] = v[11] - v[4] - v[5] 
        dydt[V.tlr2] = v[1] - v[21]
        dydt[V.infiltrate] = v[2] - v[13] - v[17] - v[19]
        dydt[V.migrate] = v[6] + v[8] - v[7] - v[15] 
        dydt[V.epi] = - v[3]
        dydt[V.ost] =  v[7]  - v[25] + v[26]
        dydt[V.il10] = v[12] + v[14] - v[30]
        dydt[V.arg] = v[16] - v[33]
        dydt[V.boneres] = v[28] - v[9]
        dydt[V.pdcd] = v[18] - v[29]
        dydt[V.ppp1r11] = v[20] - v[32]
        dydt[V.sele] = v[22] - v[24]
        dydt[V.tek] = v[23] - v[31]
        dydt[V.ob] = - v[27]

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM
    
    x[C.k1] = 2.954608484347956976e-01
    x[C.lyz2] = 2.784841990991109739e+00
    x[C.defb3] = 6.294138599364786091e-02
    x[C.cxcl3] = 2.094971408219670483e-01
    x[C.tnf] = 3.664167778561581157e-01
    x[C.ccl2] = 7.521701236445179628e-01
    x[C.ccl9] = 9.879375011256332773e-01
    x[C.csf1] = 4.447155277497639309e-02
    x[C.ctsk] = 2.348736049176628760e-01
    x[C.ext_sp] = 7.112247021035033256e+00
    x[C.ext_up] = 2.670571211752287510e+00
    x[C.il10_a] = 6.625345950371550607e-01
    x[C.il10_b] = 1.370108129670456580e-02
    x[C.mmp3] = 4.184156467059738094e-01
    x[C.inh_a] = 5.343098222515851176e-02
    x[C.inh_b] = 4.695234826085960440e-03
    x[C.arg_1] = 1.001320977133113022e-01
    x[C.arg_2] = 3.159359021463684081e-02
    x[C.pdcd1] = 1.340304969070014618e+00
    x[C.pdcd2] = 1.925396862979747548e+00
    x[C.ppp_1] = 9.481744416231211448e-01
    x[C.ppp_2] = 1.325196415323223464e-01
    x[C.tlrsele] = 5.325670953792007722e-01
    x[C.tek_1] = 1.000652014467619644e+00
    x[C.tek_2] = 2.748402524407094649e+00
    x[C.ostdeg] = 7.206624553959842250e-02
    x[C.gja] = 1.961995579353989147e-02
    x[C.rankl] = 4.225959855214848210e-03
    x[C.alpl] = 9.975530828801275340e-02
    x[C.dpdcd] = 0.3
    x[C.dil10] = 0.3
    x[C.dtek] = 0.3
    x[C.dppp1r11] = 0.3

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM
    y0[V.ext] = 1
    y0[V.epi] = 2
    y0[V.boneres] = 5
    y0[V.ob] = 1

    return y0
