import unittest
import os
from collections import OrderedDict
import numpy as np
import pandas as pd
import logging

# Assuming freeflux is installed or PYTHONPATH is set
try:
    from freeflux.core.model import Model
    from freeflux.core.emu import EMU
    from freeflux.core.metabolite import Metabolite
    from freeflux.core.mdv import MDV, get_natural_MDV, get_substrate_MDV
    from freeflux.analysis.fit import Fitter # To help setup model state for some tests
    from freeflux.analysis.simulate import Simulator # To help setup model state for some tests
    from freeflux.utils.utils import Calculator
except ImportError:
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from freeflux.core.model import Model
    from freeflux.core.emu import EMU
    from freeflux.core.metabolite import Metabolite
    from freeflux.core.mdv import MDV, get_natural_MDV, get_substrate_MDV
    from freeflux.analysis.fit import Fitter
    from freeflux.analysis.simulate import Simulator
    from freeflux.utils.utils import Calculator

# Suppress warnings during tests for cleaner output, specifically RuntimeWarning from numpy/scipy linalg
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning) # For pandas future warnings if any

class TestCalculatorMultiExperiment(unittest.TestCase):

    def setUp(self):
        """Set up a model and calculator for each test."""
        self.base_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(self.base_dir, 'data')
        self.reactions_file = os.path.join(self.data_dir, 'reactions_simple.tsv')

        self.model = Model(name="TestCalculatorModel")
        # It's good practice to load reactions if some calculator methods might depend on network structure indirectly
        # For some specific unit tests below, reactions might not be strictly needed if methods are self-contained enough.
        if os.path.exists(self.reactions_file):
            self.model.read_from_file(self.reactions_file)
        else:
            # Create a minimal reactions_simple.tsv if it doesn't exist for some reason
            # This ensures tests can run even if previous file creation steps were interrupted.
            if not os.path.exists(self.data_dir):
                os.makedirs(self.data_dir)
            with open(self.reactions_file, 'w') as f:
                f.write("#reaction_ID\treactant_IDs(atom)\tproduct_IDs(atom)\treversibility\n")
                f.write("v1\tA(a)\tB(a)\t0\n")
                f.write("v2\tB(a)\tC(a)\t0\n")
            self.model.read_from_file(self.reactions_file)


        self.calculator = Calculator(self.model)

        # Define EMUs
        self.emu_A_a = EMU('A_a', Metabolite('A', 'a'), 'a')
        self.emu_B_a = EMU('B_a', Metabolite('B', 'a'), 'a')
        self.emu_C_a = EMU('C_a', Metabolite('C', 'a'), 'a')

        # Define some default labeling strategies and measured MDVs for testing
        self.model.labeling_strategy = OrderedDict([
            ('exp1', OrderedDict([
                ('A', [['1'], [1.0], [0.99]]), # Fully labeled A in exp1
            ])),
            ('exp2', OrderedDict([
                ('A', [['0'], [1.0], [1.0]]),  # Fully unlabeled A in exp2 (natural abundance equivalent)
            ]))
        ])

        self.model.measured_MDVs = OrderedDict([
            ('exp1', OrderedDict([
                ('B_a', [np.array([0.01, 0.99]), np.array([0.001, 0.001])]),
                ('C_a', [np.array([0.01, 0.99]), np.array([0.001, 0.001])])
            ])),
            ('exp2', OrderedDict([
                ('B_a', [np.array([0.98, 0.02]), np.array([0.001, 0.001])]),
            ]))
        ])

        # Populate target_EMUs - essential for some calculator methods
        all_frags = set()
        for exp_id_key in self.model.measured_MDVs:
            all_frags.update(self.model.measured_MDVs[exp_id_key].keys())
        self.model.target_EMUs = list(all_frags)


    def test_get_experiment_labeling_strategy(self):
        strategy_exp1 = self.calculator._get_experiment_labeling_strategy('exp1')
        self.assertIn('A', strategy_exp1)
        self.assertEqual(strategy_exp1['A'], [['1'], [1.0], [0.99]])

        strategy_exp2 = self.calculator._get_experiment_labeling_strategy('exp2')
        self.assertIn('A', strategy_exp2)
        self.assertEqual(strategy_exp2['A'], [['0'], [1.0], [1.0]])

        # Test non-existent experiment
        # Suppress logging for this specific test to avoid clutter, or check log output
        with self.assertLogs(level='WARNING') as log_cm:
            strategy_non_existent = self.calculator._get_experiment_labeling_strategy('non_existent_exp')
            self.assertEqual(strategy_non_existent, {})
        self.assertTrue(any("Labeling strategy for experiment_id 'non_existent_exp' not found" in msg for msg in log_cm.output))


        # Test backward compatibility
        old_style_model = Model("OldStyle")
        old_style_model.labeling_strategy = {'A': [['1'], [1.0], [0.99]]} # Not an OrderedDict initially
        old_style_calculator = Calculator(old_style_model)

        # The Model __init__ now forces labeling_strategy to be OrderedDict.
        # To test the Calculator's backward compatibility, we'd have to bypass Model's init property.
        # Instead, let's simulate the state Calculator would see if Model somehow had an old dict:
        self.model.labeling_strategy = {'A': [['1'], [1.0], [0.99]]} # old style
        strategy_exp0_old_style = self.calculator._get_experiment_labeling_strategy('exp0')
        self.assertIn('A', strategy_exp0_old_style)
        self.assertEqual(strategy_exp0_old_style['A'], [['1'], [1.0], [0.99]])

        # Test with a new style dict that doesn't match 'exp0' directly
        self.model.labeling_strategy = OrderedDict([('some_other_exp', {'A': [['1'], [1.0], [0.99]]})])
        with self.assertLogs(level='WARNING'): # Expect warning for 'exp0' not found
             strategy_exp0_new_style_miss = self.calculator._get_experiment_labeling_strategy('exp0')
        self.assertEqual(strategy_exp0_new_style_miss, {})


    def test_get_substrate_MDVs_for_experiment(self):
        # Define end_substrates on the model for this test
        self.model.end_substrates = ['A']

        # Test exp1 (A is labeled)
        substrate_mdvs_exp1 = self.calculator._get_substrate_MDVs_for_experiment([self.emu_A_a, self.emu_B_a], 'exp1')
        self.assertIn(self.emu_A_a, substrate_mdvs_exp1)
        # Expected for A in exp1: 99% 13C1, so m1 should be high. m0 low.
        # get_substrate_MDV('a', ['1'], [1.0], [0.99]) -> MDV for 'a' (1 atom)
        # Purity 0.99 for '1': Atom is 13C with 0.99 prob, 12C with 0.01 prob.
        # Pattern '1' means this atom is considered labeled.
        # For a single atom 'a': m0 = (1-purity), m1 = purity if atom is labeled.
        # Here, pattern is '1', percentage 100%, purity 0.99. Atom 'a' is the first atom.
        # So for emu_A_a (single atom 'a'): m0 = 0.01, m1 = 0.99
        np.testing.assert_array_almost_equal(substrate_mdvs_exp1[self.emu_A_a].value, np.array([0.01, 0.99]), decimal=4)
        self.assertNotIn(self.emu_B_a, substrate_mdvs_exp1, "B_a is not an end_substrate, should not be in result")

        # Test exp2 (A is unlabeled -> natural abundance)
        substrate_mdvs_exp2 = self.calculator._get_substrate_MDVs_for_experiment([self.emu_A_a], 'exp2')
        self.assertIn(self.emu_A_a, substrate_mdvs_exp2)
        natural_A_a_mdv = get_natural_MDV(self.emu_A_a.size) # size is 1 for 'a'
        np.testing.assert_array_almost_equal(substrate_mdvs_exp2[self.emu_A_a].value, natural_A_a_mdv.value, decimal=4)

    def test_calculate_measured_MDVs_inversed_covariance_matrix(self):
        self.calculator._calculate_measured_MDVs_inversed_covariance_matrix()
        inv_cov = self.model.measured_MDVs_inv_cov
        self.assertIsNotNone(inv_cov)

        # exp1: B_a (2 isotopomers), C_a (2 isotopomers) -> 4 rows
        # exp2: B_a (2 isotopomers) -> 2 rows
        # Total rows/cols = 2+2+2 = 6
        self.assertEqual(inv_cov.shape, (6, 6))

        expected_variances = []
        expected_variances.extend(np.array([0.001, 0.001])**2) # exp1, B_a
        expected_variances.extend(np.array([0.001, 0.001])**2) # exp1, C_a
        expected_variances.extend(np.array([0.001, 0.001])**2) # exp2, B_a

        expected_diag = 1.0 / np.array(expected_variances)
        np.testing.assert_array_almost_equal(np.diag(inv_cov), expected_diag)
        # Check if it's diagonal
        self.assertEqual(np.count_nonzero(inv_cov - np.diag(np.diag(inv_cov))), 0)

    def test_generate_random_MDVs_multi_experiment_structure(self):
        original_mdvs = deepcopy(self.model.measured_MDVs)
        self.calculator._generate_random_MDVs() # This will create its own backup

        self.assertIn('exp1', self.model.measured_MDVs)
        self.assertIn('B_a', self.model.measured_MDVs['exp1'])

        # Check if means changed but SDs are the same
        self.assertFalse(np.array_equal(self.model.measured_MDVs['exp1']['B_a'][0], original_mdvs['exp1']['B_a'][0]))
        np.testing.assert_array_equal(self.model.measured_MDVs['exp1']['B_a'][1], original_mdvs['exp1']['B_a'][1])

        self.assertFalse(np.array_equal(self.model.measured_MDVs['exp2']['B_a'][0], original_mdvs['exp2']['B_a'][0]))
        np.testing.assert_array_equal(self.model.measured_MDVs['exp2']['B_a'][1], original_mdvs['exp2']['B_a'][1])

        self.calculator._reset_measured_MDVs()
        np.testing.assert_array_almost_equal(self.model.measured_MDVs['exp1']['B_a'][0], original_mdvs['exp1']['B_a'][0])
        np.testing.assert_array_almost_equal(self.model.measured_MDVs['exp1']['B_a'][1], original_mdvs['exp1']['B_a'][1])
        np.testing.assert_array_almost_equal(self.model.measured_MDVs['exp2']['B_a'][0], original_mdvs['exp2']['B_a'][0])


if __name__ == '__main__':
    unittest.main()
