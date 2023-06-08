from biomass.models import periodontitis_model
from biomass import run_analysis
from tqdm import tqdm
from biomass import Model, run_simulation,optimize
model = Model(periodontitis_model.__package__).create()

for x_id in tqdm(range(1, 21)):
    optimize(model, x_id=x_id, options={
            "popsize": 3,
            "max_generation": 100,
            "allowable_error": 0.5,
            "local_search_method": "DE",
            "maxiter": 50,
            "overwrite":True,
        })