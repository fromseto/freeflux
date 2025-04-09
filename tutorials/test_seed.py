'''Example of Flux estimation at steady state with a E. coli model using synthetic data.
'''


from os import makedirs
import pandas as pd
from freeflux import Model

# solver = 'slsqp'
# solver = 'trust-constr'
solver = 'ralg'


MODEL_FILE = '../models/ecoli/experimental_data/reactions.xlsx' 
MEASURED_MDVS = '../models/ecoli/experimental_data/measured_MDVs.xlsx'
MEASURED_FLUXES = '../models/ecoli/experimental_data/measured_fluxes.xlsx'
OUT_DIR = f'/home/user/tmp/flux_seed/{solver}'

DILUTION_FROM = [
    'CO2u', 
    'Alau', 
    'Glyu', 
    'Valu', 
    'Leuu', 
    'Ileu', 
    'Seru', 
    'Pheu', 
    'Aspu', 
    'Gluu', 
    'Tyru'
]


# estimate fluxes at steady state
def ecoli_steady_state_fitting():

    ecoli = Model('ecoli')
    ecoli.read_from_file(MODEL_FILE)
    
    with ecoli.fitter('ss') as fit:
        # specify the lableing strategy, 
        # use this method for every labeled substrate
        fit.set_labeling_strategy(
            'Glc.ex', 
            labeling_pattern = ['100000', '111111'], 
            percentage = [0.77, 0.205], 
            purity = [0.99, 0.985]
        )
        
        # read measurements
        fit.set_measured_MDVs_from_file(MEASURED_MDVS)
        fit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set upper and lower bounds for fluxes
        fit.set_flux_bounds('all', bounds = [-100, 100]) 

        # solve the fluxes
        fit.prepare(
            dilution_from = DILUTION_FROM, 
            n_jobs = 30
        )

        initial_seed = 42
        while True:
            # Create an initial random number generator
            # initial_rng = np.random.default_rng(initial_seed)
            res = fit.solve(solver = solver, seed=42, max_iters=1000)
            print(res.optimization_successful)
            if res.optimization_successful:
                break
        # res = fit.solve(solver = solver, seed=100, max_iters=3000)
        # print(res.optimization_successful)
    
    # save the results
    pd.Series(res.opt_net_fluxes).to_csv(
        OUT_DIR+'/estimated_net_fluxes.csv'
    )
    pd.Series(res.opt_total_fluxes).to_csv(
        OUT_DIR+'/estimated_total_fluxes.csv'
    )

    net_cis = res.estimate_confidence_intervals(
        which = 'net', 
        confidence_level = 0.95
    )
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_csv(
        OUT_DIR+'/netflux_le_CIs.csv'
    )
    
    total_cis = res.estimate_confidence_intervals(
        which = 'total', 
        confidence_level = 0.95
    )
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_csv(
        OUT_DIR+'/totalflux_le_CIs.csv'
    )

    # normal probability plot of residuals
    res.plot_normal_probability(show_fig = False, output_dir = OUT_DIR)
    
    # compare simulations and measurements
    res.plot_simulated_vs_measured_MDVs(show_fig = False, output_dir = OUT_DIR)
    res.plot_simulated_vs_measured_fluxes(show_fig = False, output_dir = OUT_DIR)
    
    # export the contribution matrix
    res.estimate_contribution_matrix(which = 'net').to_csv(
        OUT_DIR+'/netflux_contribMat.csv'
    )
    res.estimate_contribution_matrix(which = 'total').to_csv(
        OUT_DIR+'/totalflux_contribMat.csv'
    )
    
    # export the sensitivity matrix
    res.estimate_sensitivity(which = 'net').to_csv(
        OUT_DIR+'/netflux_senMat.csv'
    )
    res.estimate_sensitivity(which = 'total').to_csv(
        OUT_DIR+'/totalflux_senMat.csv'
    )

if __name__ == '__main__':
    ecoli_steady_state_fitting()
