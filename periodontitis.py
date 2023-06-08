from biomass.models import periodontitis_model
from biomass import run_analysis
from tqdm import tqdm
from biomass import Model, run_simulation,optimize
model = Model(periodontitis_model.__package__).create()

##simulate model
run_simulation(model, viz_type='best', show_all=False, stdev=True)

##sensitivity analysis
run_analysis(model, target='parameter', metric='integral', style='barplot',options = {
        'excluded_params': [
            'EGF', 'HRG', 'no_ligand', 'Ligand', 'gpx','mmp3'
        ]
    })

##output parameter range
from biomass.result import OptimizationResults
res = OptimizationResults(model)
# Export estimated parameters in CSV format
res.to_csv()
# Visualize estimated parameter sets
res.savefig(figsize=(25,5), boxplot_kws={"orient": "v"})