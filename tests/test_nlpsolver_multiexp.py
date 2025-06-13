import unittest
import os
from collections import OrderedDict
import numpy as np
import pandas as pd
import logging
from unittest.mock import MagicMock # Added for mocking

# Assuming freeflux is installed or PYTHONPATH is set
try:
    from freeflux.core.model import Model
    from freeflux.core.emu import EMU
    from freeflux.core.metabolite import Metabolite
    from freeflux.utils.utils import Calculator
    from freeflux.solver.nlpsolver import MFAModel
    # For setup convenience, though Fitter/Simulator aren't directly tested here
    from freeflux.analysis.fit import Fitter
except ImportError:
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from freeflux.core.model import Model
    from freeflux.core.emu import EMU
    from freeflux.core.metabolite import Metabolite
    from freeflux.utils.utils import Calculator
    from freeflux.solver.nlpsolver import MFAModel
    from freeflux.analysis.fit import Fitter

# Suppress warnings during tests for cleaner output
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

class TestNLPSolverMultiExperiment(unittest.TestCase):

    def setUp(self):
        """Set up a model and calculator for each test."""
        self.base_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(self.base_dir, 'data')
        self.reactions_file = os.path.join(self.data_dir, 'reactions_simple.tsv')

        self.model = Model(name="TestNLPSolverModel")

        if not os.path.exists(self.reactions_file):
            if not os.path.exists(self.data_dir):
                os.makedirs(self.data_dir)
            with open(self.reactions_file, 'w') as f:
                f.write("#reaction_ID\treactant_IDs(atom)\tproduct_IDs(atom)\treversibility\n")
                f.write("v1\tA(a)\tB(a)\t0\n")
                f.write("v2\tB(a)\tC(a)\t0\n")
        self.model.read_from_file(self.reactions_file)

        self.model.labeling_strategy = OrderedDict([
            ('exp1', OrderedDict([('A', [['1'], [1.0], [1.0]]))])),
            ('exp2', OrderedDict([('A', [['0'], [1.0], [1.0]]))]))
        ])

        self.model.measured_MDVs = OrderedDict([
            ('exp1', OrderedDict([
                ('B_a', [np.array([0.0, 1.0]), np.array([0.01, 0.01])]),
                ('C_a', [np.array([0.0, 1.0]), np.array([0.01, 0.01])]),
            ])),
            ('exp2', OrderedDict([
                ('B_a', [np.array([0.989, 0.011]), np.array([0.01, 0.01])]),
            ]))
        ])

        all_frags = set()
        for exp_id_key in self.model.measured_MDVs:
            all_frags.update(self.model.measured_MDVs[exp_id_key].keys())
        self.model.target_EMUs = list(all_frags)

        self.model.totalfluxids = ['v1', 'v2']
        self.model.netfluxids = ['v1', 'v2']
        # For A->B->C, if v1 is free flux u[0], then total fluxes v = [u[0], u[0]]
        self.model.null_space = np.array([[1.0], [1.0]])
        self.model.transform_matrix = np.eye(2) # Assuming net fluxes = total fluxes (irreversible)

        self.calculator = Calculator(self.model)
        self.calculator._calculate_measured_MDVs_inversed_covariance_matrix()

        self.mfamodel = MFAModel(self.model, fit_measured_fluxes=False, solver='slsqp')

    def test_mfamodel_objective_multi_experiment(self):
        u_known = np.array([1.0])

        mock_sim_B_exp1 = np.array([0.0, 1.0])
        mock_sim_C_exp1 = np.array([0.0, 1.0])
        mock_sim_B_exp2 = np.array([0.989, 0.011])

        def mock_mdv_calculator_func(params_u, experiment_id):
            if experiment_id == 'exp1':
                return {'B_a': mock_sim_B_exp1, 'C_a': mock_sim_C_exp1}
            elif experiment_id == 'exp2':
                return {'B_a': mock_sim_B_exp2}
            return {}

        self.mfamodel.calculator._calculate_MDVs_for_experiment = MagicMock(
            side_effect=mock_mdv_calculator_func
        )

        residuals_exp1_Ba = mock_sim_B_exp1 - self.model.measured_MDVs['exp1']['B_a'][0]
        residuals_exp1_Ca = mock_sim_C_exp1 - self.model.measured_MDVs['exp1']['C_a'][0]
        residuals_exp2_Ba = mock_sim_B_exp2 - self.model.measured_MDVs['exp2']['B_a'][0]

        all_residuals = np.concatenate([residuals_exp1_Ba, residuals_exp1_Ca, residuals_exp2_Ba])
        inv_cov_matrix = self.model.measured_MDVs_inv_cov
        obj_val_expected = 0.5 * all_residuals.T @ inv_cov_matrix @ all_residuals

        self.mfamodel.build_objective()
        obj_val_calculated = self.mfamodel.f(u_known)

        self.assertAlmostEqual(obj_val_calculated, obj_val_expected, places=6)

        self.mfamodel.calculator._calculate_MDVs_for_experiment.assert_any_call(u_known, experiment_id='exp1')
        self.mfamodel.calculator._calculate_MDVs_for_experiment.assert_any_call(u_known, experiment_id='exp2')
        self.assertEqual(self.mfamodel.calculator._calculate_MDVs_for_experiment.call_count, 2)

    def test_mfamodel_gradient_multi_experiment(self):
        u_known = np.array([1.0])

        # Mock simulated MDVs for obj_func to calculate residuals
        mock_sim_B_exp1 = np.array([0.0, 1.0]) # Perfect match for exp1 B_a
        mock_sim_C_exp1 = np.array([0.0, 0.9]) # Slight mismatch for exp1 C_a
        mock_sim_B_exp2 = np.array([0.98, 0.02])# Slight mismatch for exp2 B_a

        def mock_mdv_calculator_func_for_grad(params_u, experiment_id):
            if experiment_id == 'exp1':
                return {'B_a': mock_sim_B_exp1, 'C_a': mock_sim_C_exp1}
            elif experiment_id == 'exp2':
                return {'B_a': mock_sim_B_exp2}
            return {}

        self.mfamodel.calculator._calculate_MDVs_for_experiment = MagicMock(
            side_effect=mock_mdv_calculator_func_for_grad
        )

        # Mock the stacked derivative matrix from calculator
        # Order: B_a (exp1)[m0,m1], C_a (exp1)[m0,m1], B_a (exp2)[m0,m1] -> 6 rows
        # Assuming 1 free flux -> 1 col. Shape (6,1)
        der_B_exp1_m1 = 0.1  # d(B_a_m1)/du1 for exp1
        der_C_exp1_m1 = 0.2  # d(C_a_m1)/du1 for exp1
        der_B_exp2_m1 = 0.05 # d(B_a_m1)/du1 for exp2

        # Assuming d(m0)/du = -d(m1)/du due to sum(m_i)=1
        predefined_stacked_der_matrix = np.array([
            [-der_B_exp1_m1], [der_B_exp1_m1],
            [-der_C_exp1_m1], [der_C_exp1_m1],
            [-der_B_exp2_m1], [der_B_exp2_m1]
        ])

        # _calculate_MDVs_and_derivatives_p is called by obj_grad
        self.mfamodel.calculator._calculate_MDVs_and_derivatives_p = MagicMock(
            return_value=(None, predefined_stacked_der_matrix) # First val (simMDVs) is ignored by obj_grad
        )

        self.mfamodel.build_objective()
        self.mfamodel.build_gradient()

        # Call obj_func to populate self.mfamodel.residuals_mdv
        self.mfamodel.f(u_known)
        self.assertIsNotNone(getattr(self.mfamodel, 'residuals_mdv', None), "residuals_mdv not set by obj_func")
        self.assertEqual(self.mfamodel.residuals_mdv.shape[0], 6) # 2+2+2 elements

        # Manual gradient calculation: grad = residuals.T @ inv_cov @ derivatives
        inv_cov_matrix = self.model.measured_MDVs_inv_cov
        grad_expected_mdv_part = self.mfamodel.residuals_mdv.T @ inv_cov_matrix @ predefined_stacked_der_matrix
        grad_expected = grad_expected_mdv_part.flatten()

        grad_calculated = self.mfamodel.df(u_known) # Call obj_grad

        self.assertTrue(np.allclose(grad_calculated, grad_expected, atol=1e-6),
                        f"Gradient mismatch: \nExpected: {grad_expected}\nCalculated: {grad_calculated}")

        # Verify mock calls
        self.mfamodel.calculator._calculate_MDVs_for_experiment.assert_any_call(u_known, experiment_id='exp1')
        self.mfamodel.calculator._calculate_MDVs_for_experiment.assert_any_call(u_known, experiment_id='exp2')
        self.mfamodel.calculator._calculate_MDVs_and_derivatives_p.assert_called_once()


if __name__ == '__main__':
    unittest.main()
