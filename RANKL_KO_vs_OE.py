from biomass.models import RANKL_KO_vs_OE_model
from biomass import run_analysis
from tqdm import tqdm
from biomass import Model, run_simulation,optimize
model = Model(RANKL_KO_vs_OE_model.__package__).create()

##simulate model
run_simulation(model, viz_type='best', show_all=False, stdev=True)