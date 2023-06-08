from biomass.models import CCL2_KO_model
from biomass import run_analysis
from tqdm import tqdm
from biomass import Model, run_simulation,optimize
model = Model(CCL2_KO_model.__package__).create()

##simulate model
run_simulation(model, viz_type='best', show_all=False, stdev=True)