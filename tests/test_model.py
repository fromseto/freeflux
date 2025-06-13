import unittest
import os
from collections import OrderedDict
import numpy as np

# Adjust import path if tests are run from root or within tests/
# Assuming freeflux is installed or PYTHONPATH is set to find freeflux package from project root
try:
    from freeflux.core.model import Model
    from freeflux.analysis.fit import Fitter
    from freeflux.analysis.simulate import Simulator # Added Simulator import
except ImportError:
    # Fallback for running script directly from tests directory
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from freeflux.core.model import Model
    from freeflux.analysis.fit import Fitter
    from freeflux.analysis.simulate import Simulator # Added Simulator import


class TestModelMultiExperiment(unittest.TestCase):

    def setUp(self):
        """Set up paths to data files."""
        self.base_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(self.base_dir, 'data')

        self.reactions_file = os.path.join(self.data_dir, 'reactions_simple.tsv')
        self.mdvs_exp1_file = os.path.join(self.data_dir, 'mdvs_exp1.tsv')
        self.mdvs_exp2_multi_file = os.path.join(self.data_dir, 'mdvs_exp2_multi.tsv')

        # Mock files are assumed to be created by the agent/setup.
        # If running standalone, ensure these files exist with the specified content.

    def test_load_mdvs_default_experiment(self):
        model = Model(name="TestModelDefaultMDV")
        model.read_from_file(self.reactions_file) # Load reactions for target_EMU testing

        fitter = Fitter(model)
        fitter.set_measured_MDVs_from_file(self.mdvs_exp1_file)

        self.assertIn('exp0', model.measured_MDVs, "Default experiment 'exp0' not found.")
        self.assertEqual(len(model.measured_MDVs), 1, "There should be only one experiment ('exp0').")

        exp0_data = model.measured_MDVs['exp0']
        self.assertIn('B_a', exp0_data)
        self.assertIn('C_a', exp0_data)
        np.testing.assert_array_almost_equal(exp0_data['B_a'][0], np.array([0.5, 0.5]))
        np.testing.assert_array_almost_equal(exp0_data['B_a'][1], np.array([0.01, 0.01]))
        np.testing.assert_array_almost_equal(exp0_data['C_a'][0], np.array([0.4, 0.6]))
        np.testing.assert_array_almost_equal(exp0_data['C_a'][1], np.array([0.01, 0.01]))

        fitter._decompose_network(n_jobs=1)
        self.assertIsInstance(model.target_EMUs, list)
        self.assertCountEqual(model.target_EMUs, ['B_a', 'C_a'], "target_EMUs for default experiment.")

    def test_load_mdvs_multiple_experiments(self):
        model = Model(name="TestModelMultiMDV")
        model.read_from_file(self.reactions_file)

        fitter = Fitter(model)
        fitter.set_measured_MDVs_from_file(self.mdvs_exp2_multi_file)

        self.assertEqual(len(model.measured_MDVs), 2, "Should contain data for 2 experiments.")
        self.assertIn('expA', model.measured_MDVs)
        self.assertIn('expB', model.measured_MDVs)

        expA_data = model.measured_MDVs['expA']
        self.assertIn('B_a', expA_data)
        self.assertIn('C_a', expA_data)
        np.testing.assert_array_almost_equal(expA_data['B_a'][0], np.array([0.5, 0.5]))
        np.testing.assert_array_almost_equal(expA_data['C_a'][1], np.array([0.01, 0.01]))

        expB_data = model.measured_MDVs['expB']
        self.assertIn('B_a', expB_data)
        self.assertIn('D_a', expB_data)
        np.testing.assert_array_almost_equal(expB_data['B_a'][0], np.array([0.6, 0.4]))
        np.testing.assert_array_almost_equal(expB_data['D_a'][0], np.array([0.2, 0.8]))

        fitter._decompose_network(n_jobs=1)
        self.assertIsInstance(model.target_EMUs, list)
        self.assertCountEqual(model.target_EMUs, ['B_a', 'C_a', 'D_a'], "target_EMUs aggregated from multiple experiments.")

    def test_initial_labeling_strategy_type(self):
        model = Model(name="TestInitialLabeling")
        self.assertIsInstance(model.labeling_strategy, OrderedDict, "labeling_strategy should be OrderedDict on init.")

    def test_set_labeling_strategy_single_experiment_default(self):
        model = Model(name="TestLabelingDefault")
        sim = Simulator(model) # Simulator is where set_labeling_strategy is defined

        sim.set_labeling_strategy(
            labeled_substrate='Sub1',
            labeling_pattern=['1'],
            percentage=[1.0],
            purity=[0.99]
        )
        # experiment_id defaults to 'exp0'
        self.assertIn('exp0', model.labeling_strategy)
        self.assertIn('Sub1', model.labeling_strategy['exp0'])
        self.assertEqual(model.labeling_strategy['exp0']['Sub1'], [['1'], [1.0], [0.99]])
        self.assertIsInstance(model.labeling_strategy['exp0'], OrderedDict, "Inner dict for exp0 should be OrderedDict.")

    def test_set_labeling_strategy_single_experiment_explicit_id(self):
        model = Model(name="TestLabelingExplicit")
        sim = Simulator(model)

        sim.set_labeling_strategy(
            labeled_substrate='Sub1',
            labeling_pattern=['1'],
            percentage=[1.0],
            purity=[0.99],
            experiment_id='expTest'
        )
        self.assertIn('expTest', model.labeling_strategy)
        self.assertNotIn('exp0', model.labeling_strategy) # Ensure it doesn't also create exp0
        self.assertIn('Sub1', model.labeling_strategy['expTest'])
        self.assertEqual(model.labeling_strategy['expTest']['Sub1'], [['1'], [1.0], [0.99]])
        self.assertIsInstance(model.labeling_strategy['expTest'], OrderedDict)

    def test_set_labeling_strategy_multiple_experiments(self):
        model = Model(name="TestLabelingMulti")
        sim = Simulator(model)

        sim.set_labeling_strategy('SubA', ['1'], [1.0], [0.99], experiment_id='expA')
        sim.set_labeling_strategy('SubB', ['1'], [1.0], [0.99], experiment_id='expB')
        sim.set_labeling_strategy('SubCommon', ['0'], [1.0], [0.99], experiment_id='expA')
        sim.set_labeling_strategy('SubCommon', ['1'], [0.5], [0.98], experiment_id='expB')

        self.assertIn('expA', model.labeling_strategy)
        self.assertIn('expB', model.labeling_strategy)

        self.assertIn('SubA', model.labeling_strategy['expA'])
        self.assertIn('SubCommon', model.labeling_strategy['expA'])
        self.assertEqual(model.labeling_strategy['expA']['SubCommon'], [['0'], [1.0], [0.99]])

        self.assertIn('SubB', model.labeling_strategy['expB'])
        self.assertIn('SubCommon', model.labeling_strategy['expB'])
        self.assertEqual(model.labeling_strategy['expB']['SubCommon'], [['1'], [0.5], [0.98]])

    def test_unset_labeling_strategy(self):
        model = Model(name="TestUnsetLabeling")
        sim = Simulator(model)

        sim.set_labeling_strategy('Sub1', ['1'], [1.0], [0.99], experiment_id='expA')
        sim.set_labeling_strategy('Sub2', ['0'], [1.0], [0.99], experiment_id='expA')
        sim.set_labeling_strategy('Sub3', ['1'], [1.0], [0.99], experiment_id='expB')

        # Unset Sub1 from expA
        sim._unset_labeling_strategy('expA', 'Sub1') # Calling private method for direct testing
        self.assertNotIn('Sub1', model.labeling_strategy['expA'])
        self.assertIn('Sub2', model.labeling_strategy['expA'])
        self.assertIn('expB', model.labeling_strategy) # expB should be unaffected

        # Unset Sub2 from expA (making expA empty)
        sim._unset_labeling_strategy('expA', 'Sub2')
        self.assertNotIn('expA', model.labeling_strategy, "Experiment 'expA' should be removed if empty.")
        self.assertIn('expB', model.labeling_strategy) # expB should still be there
        self.assertIn('Sub3', model.labeling_strategy['expB'])


if __name__ == '__main__':
    unittest.main()
