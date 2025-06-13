import unittest
import os
from collections import OrderedDict
import numpy as np
import pandas as pd # Though not directly used in this test logic, often good to have for debugging
import logging

# Assuming freeflux is installed or PYTHONPATH is set
try:
    from freeflux.core.model import Model
    from freeflux.analysis.fit import Fitter
except ImportError:
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from freeflux.core.model import Model
    from freeflux.analysis.fit import Fitter

# Suppress warnings and info logs during tests for cleaner output
# logging.basicConfig(level=logging.ERROR) # Suppress info/warnings from freeflux
# warnings.filterwarnings('ignore')

class TestFittingMultiExperiment(unittest.TestCase):

    def setUp(self):
        """Set up paths to data files."""
        self.base_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(self.base_dir, 'data')

        self.reactions_file = os.path.join(self.data_dir, 'reactions_simple.tsv')
        self.mdvs_exp1_file = os.path.join(self.data_dir, 'mdvs_exp1.tsv')
        # mdvs_exp2_multi.tsv is not directly used here, data is set manually for better control

        # Ensure reactions_simple.tsv exists (as it's loaded in tests)
        if not os.path.exists(self.reactions_file):
            if not os.path.exists(self.data_dir):
                os.makedirs(self.data_dir)
            with open(self.reactions_file, 'w') as f:
                f.write("#reaction_ID\treactant_IDs(atom)\tproduct_IDs(atom)\treversibility\n")
                f.write("v1\tA(a)\tB(a)\t0\n")
                f.write("v2\tB(a)\tC(a)\t0\n")

        # Ensure mdvs_exp1.tsv for backward compatibility test
        if not os.path.exists(self.mdvs_exp1_file):
            with open(self.mdvs_exp1_file, 'w') as f:
                f.write("fragment_ID\tmean\tsd\n")
                f.write("B_a\t0.5,0.5\t0.01,0.01\n") # Corresponds to 50% labeled A
                f.write("C_a\t0.5,0.5\t0.01,0.01\n") # Corresponds to 50% labeled A (via B)


    def test_single_experiment_fitting_backward_compatibility(self):
        model = Model(name="TestFitSingleExp")
        model.read_from_file(self.reactions_file) # A(a)->B(a), B(a)->C(a)

        fit = Fitter(model)

        # Labeling: A is 50% 13C1 (pattern '1'), 50% natural (pattern '0')
        # Purity of labeled part is 100%.
        fit.set_labeling_strategy(
            labeled_substrate='A',
            labeling_pattern=['1'], # Only specify the labeled part
            percentage=[0.5],      # 50% is this labeled form
            purity=[1.0]           # Labeled part is 100% pure 13C
            # experiment_id defaults to 'exp0'
        )
        # The other 50% will be natural 'A' by default.
        # So, B and C should be roughly 0.5 m0, 0.5 m1 if A is single carbon.
        # mdvs_exp1.tsv has B_a [0.5,0.5] and C_a [0.5,0.5] (changed from original for consistency)

        fit.set_measured_MDVs_from_file(self.mdvs_exp1_file)
        fit.set_flux_bounds('all', bounds=[-100, 100])
        fit.set_measured_flux('v1', mean=1.0, sd=0.1)

        # Suppress solver output for tests if possible, or check FreeFlux options
        fit.prepare(n_jobs=1)
        res = fit.solve(show_progress=False, tol=1e-4) # Added tol for potentially faster convergence

        self.assertTrue(res.success, f"Fitting failed for single experiment. Solver message: {res.message}")
        self.assertIsNotNone(res.opt_net_fluxes, "Net fluxes should be calculated.")
        if res.opt_net_fluxes is not None:
            self.assertIn('v1', res.opt_net_fluxes)
            self.assertAlmostEqual(res.opt_net_fluxes['v1'], 1.0, delta=0.2, msg="v1 flux not close to target.")
            # For A->B->C, v2 should be equal to v1 at steady state.
            if 'v2' in res.opt_net_fluxes:
                 self.assertAlmostEqual(res.opt_net_fluxes['v2'], res.opt_net_fluxes['v1'], delta=1e-3, msg="v2 should equal v1.")
        self.assertTrue(res.ssr > 0, "SSR should be positive.")


    def test_multi_experiment_fitting_simple_case(self):
        model = Model(name="TestFitMultiExp")
        model.read_from_file(self.reactions_file)

        fit = Fitter(model)

        # Experiment A: 'A' is 80% 13C labeled (m1=0.8), purity 1.0
        fit.set_labeling_strategy(
            labeled_substrate='A',
            labeling_pattern=['1'],
            percentage=[0.8],
            purity=[1.0],
            experiment_id='expA'
        )
        # Experiment B: 'A' is natural (effectively 100% pattern '0')
        fit.set_labeling_strategy(
            labeled_substrate='A',
            labeling_pattern=['0'],
            percentage=[1.0],
            purity=[1.0], # Purity of '0' pattern is nominal
            experiment_id='expB'
        )

        # Set Measured MDVs consistent with labeling and a target flux (e.g., v1=0.7)
        # For expA (A is 80% m1, 20% m0 (natural approx [0.99,0.01]))
        # Effective input A for expA: 0.8*[0,1] + 0.2*[0.99,0.01] = [0.198, 0.802]
        fit.set_measured_MDV('B_a', np.array([0.198, 0.802]), np.array([0.01, 0.01]), experiment_id='expA')
        fit.set_measured_MDV('C_a', np.array([0.198, 0.802]), np.array([0.01, 0.01]), experiment_id='expA')

        # For expB (A is natural approx [0.99,0.01])
        fit.set_measured_MDV('B_a', np.array([0.989, 0.011]), np.array([0.01, 0.01]), experiment_id='expB')
        # No C_a for expB to make it slightly different

        fit.set_flux_bounds('all', bounds=[-100, 100])
        fit.set_measured_flux('v1', mean=0.7, sd=0.05)

        fit.prepare(n_jobs=1)
        res = fit.solve(show_progress=False, tol=1e-4)

        self.assertTrue(res.success, f"Fitting failed for multi-experiment. Solver message: {res.message}")
        self.assertIsNotNone(res.opt_net_fluxes)
        if res.opt_net_fluxes is not None:
            self.assertIn('v1', res.opt_net_fluxes)
            self.assertAlmostEqual(res.opt_net_fluxes['v1'], 0.7, delta=0.15, msg="v1 flux not close to target for multi-exp.")
            if 'v2' in res.opt_net_fluxes:
                 self.assertAlmostEqual(res.opt_net_fluxes['v2'], res.opt_net_fluxes['v1'], delta=1e-3, msg="v2 should equal v1 for multi-exp.")

        self.assertTrue(res.ssr > 0, "SSR should be positive for multi-exp.")
        # A very small SSR is expected if data is perfectly consistent.
        # Given the direct mapping A->B->C, and consistent MDVs, SSR should be very low.
        self.assertTrue(res.ssr < 1.0, f"SSR {res.ssr} seems too high for this simple consistent case.")


if __name__ == '__main__':
    unittest.main()
