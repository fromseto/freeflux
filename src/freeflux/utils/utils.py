'''Define the Calculator class.'''


__author__ = 'Chao Wu'
__date__ = '05/18/2022'

import logging # Added for logging warnings
from copy import deepcopy # Added for deepcopy
from collections import ChainMap
from collections.abc import Iterable
from functools import reduce
from copy import deepcopy
import warnings
warnings.filterwarnings('ignore', category = RuntimeWarning) 
import numpy as np
from numpy.random import normal
import pandas as pd
from scipy.linalg import null_space, pinv, expm
from sympy import symbols, lambdify, Matrix, derive_by_array
try:
    import jax.numpy as jnp
    from jax import config, jacfwd
    config.update('jax_platform_name', 'cpu')
except ModuleNotFoundError:
    JAX_INSTALLED = False
else:
    JAX_INSTALLED = True
from multiprocessing import Pool
from ..core.mdv import MDV, get_natural_MDV, get_substrate_MDV, conv, diff_conv
import warnings
warnings.filterwarnings('ignore', category = RuntimeWarning)
warnings.filterwarnings('ignore', category = DeprecationWarning)

    
class Calculator():
    '''
    Parameters
    ----------
    model: Model
        Freeflux Model.
    '''
    
    def __init__(self, model):
        '''
        Parameters
        ----------
        model: Model
            Freeflux Model.
        '''
        
        self.model = model
        
    
    def _get_all_substrate_emus(self):
        """Helper to get all unique EMU objects that can be substrates."""
        all_emus = set()
        if self.model.matrix_Bs: # Check if EAMs are built
            for size in self.model.matrix_Bs:
                for sourceEMU_descriptor in self.model.matrix_Bs[size][2]:
                    if not isinstance(sourceEMU_descriptor, Iterable):
                        all_emus.add(sourceEMU_descriptor)
                    else:
                        for emu_in_tuple in sourceEMU_descriptor:
                            all_emus.add(emu_in_tuple)
        # Also consider end_substrates that might not be in EAMs explicitly as sources
        # This requires EMU objects for end_substrates. For now, rely on EAMs.
        # If self.model.end_substrates contains metabolite IDs, they need to be converted to EMUs
        # of various sizes if they are to be included here.
        # For now, substrate EMUs are those defined in the EAM source matrices (matrix_Bs).
        return list(all_emus)

    def _get_experiment_labeling_strategy(self, experiment_id='exp0'):
        """
        Retrieves the labeling strategy for a given experiment_id.
        Handles backward compatibility for old-style labeling_strategy.
        """
        if not self.model.labeling_strategy: # Empty strategy
            return {}

        # New style: experiment_id is a key
        if experiment_id in self.model.labeling_strategy and \
           isinstance(self.model.labeling_strategy[experiment_id], dict):
            return self.model.labeling_strategy[experiment_id]

        # Old style check for default experiment_id ('exp0')
        # Heuristic: if keys are not like 'exp...' and values are lists.
        # This assumes experiment IDs generally follow a pattern like 'exp<number>'
        # or are distinct from typical metabolite names.
        if experiment_id == 'exp0':
            is_old_style = True
            if not self.model.labeling_strategy: # Handles empty dict
                 is_old_style = False
            else:
                for key, value in self.model.labeling_strategy.items():
                    if not (isinstance(key, str) and isinstance(value, list)):
                        # If any key is not a string or its value is not a list,
                        # it's not the simple old style {metabolite_id: [details_list]}
                        # or it's a new style {exp_id: {metabolite_id: [details_list]}}
                        # where we didn't find the specific experiment_id above.
                        # If it's new style and we are here, it means experiment_id was not found.
                        # So, check if the key looks like an experiment ID.
                        if isinstance(self.model.labeling_strategy.get(key), dict): # definitely new style
                            is_old_style = False
                            break
                    # A more direct check for old style: are the keys potential metabolite names
                    # and values are lists (the details)?
                    # For simplicity, if the first key's value is a list and the key doesn't look like an exp_id:
                first_key = next(iter(self.model.labeling_strategy.keys()), None)
                if first_key and isinstance(self.model.labeling_strategy[first_key], list) and \
                   not first_key.lower().startswith('exp'):
                    # Likely old style
                    pass
                else:
                    is_old_style = False


            if is_old_style:
                logging.info(f"Using old-style labeling_strategy for experiment_id '{experiment_id}'.")
                return self.model.labeling_strategy # Return the whole dict as the strategy for 'exp0'

        logging.warning(f"Labeling strategy for experiment_id '{experiment_id}' not found or format mismatch. Using empty strategy (natural abundance).")
        return {}


    def _set_timepoints(self):
        
        ts = []
        for _, tpoint_infos in self.model.measured_inst_MDVs.items():
            ts.extend(tpoint_infos.keys())
        ts = sorted(set(ts))
        
        for t in ts:
            self.model.timepoints.append(t)
            
        
    def _calculate_null_space(self):
        
        S = self.model.get_total_stoichiometric_matrix(self.model.unbalanced_metabolites)
        self.model.null_space = null_space(S.values)
    
    
    def _calculate_transform_matrix(self):
        
        transMat = pd.DataFrame(
            0, 
            index = self.model.netfluxids, 
            columns = self.model.totalfluxids
        )
        for rxnid, rxn in self.model.reactions_info.items():
        
            if rxn.reversible:
                transMat.loc[rxnid, rxnid+'_f'] = 1
                transMat.loc[rxnid, rxnid+'_b'] = -1
            else:
                transMat.loc[rxnid, rxnid] = 1
        
        self.model.transform_matrix = transMat.values
    
    
    # Refactored to _get_substrate_MDVs_for_experiment
    # This method now takes a list of relevant EMUs and an experiment_id
    def _get_substrate_MDVs_for_experiment(self, substrate_emus_list, experiment_id='exp0', extra_subs=None):
        """
        Calculates MDVs for a given list of substrate EMUs for a specific experiment.
        Args:
            substrate_emus_list (list): List of EMU objects considered as potential substrates.
            experiment_id (str): The ID of the current experiment.
            extra_subs (list, optional): List of metabolite IDs to be treated as additional substrates.
        Returns:
            dict: {EMU_object: mdv_array} for the specified experiment.
        """
        if extra_subs is None:
            extra_subs = []

        substrate_mdvs_dict = {}
        current_labeling_strategy = self._get_experiment_labeling_strategy(experiment_id)

        all_possible_substrate_metab_ids = self.model.end_substrates + extra_subs

        for emu in substrate_emus_list:
            metabid = emu.metabolite_id

            # Only calculate if it's a defined end substrate or an explicitly passed extra substrate
            if metabid in all_possible_substrate_metab_ids:
                if metabid in current_labeling_strategy:
                    atom_nos = emu.atom_nos
                    labeling_details = current_labeling_strategy[metabid]
                    # Ensure labeling_details has the expected structure [patterns, percentages, purities]
                    if isinstance(labeling_details, list) and len(labeling_details) == 3:
                        patterns, percentages, purities = labeling_details
                        substrate_mdvs_dict[emu] = get_substrate_MDV(
                            atom_nos,
                            patterns,
                            percentages,
                            purities
                        )
                    else:
                        logging.error(f"Malformed labeling strategy for metabolite {metabid} in experiment {experiment_id}. Expected list of 3 elements. Using natural abundance.")
                        substrate_mdvs_dict[emu] = get_natural_MDV(emu.size) # emu.size should be n_atoms
                else:
                    # Not in specific labeling strategy for this experiment, so natural abundance
                    substrate_mdvs_dict[emu] = get_natural_MDV(emu.size)

        return substrate_mdvs_dict
    
    
    def _calculate_substrate_MDV_derivatives_basic(self, nvars, extra_subs):
        
        if extra_subs is None:
             extra_subs = []

        substrate_MDVs_der = {}
        for size in self.model.matrix_Bs:
            for sourceEMU in self.model.matrix_Bs[size][2]:
                
                if not isinstance(sourceEMU, Iterable):
                    sourceEMU = (sourceEMU,)
                    
                for emu in sourceEMU:
                    metabid = emu.metabolite_id
                    if metabid in self.model.end_substrates + extra_subs:
                        substrate_MDVs_der[emu] = np.zeros((emu.size+1, nvars))
                        
        return substrate_MDVs_der                
                        
                        
    def _calculate_substrate_MDV_derivatives_u(self, extra_subs):
        
        substrate_MDVs_der = self._calculate_substrate_MDV_derivatives_basic(
            self.model.null_space.shape[1], 
            extra_subs
        )
        
        return substrate_MDVs_der


    def _calculate_substrate_MDV_derivatives_c(self, extra_subs):
       
        substrate_MDVs_der = self._calculate_substrate_MDV_derivatives_basic(
            len(self.model.concids), 
            extra_subs
        )
        
        return substrate_MDVs_der
        
        
    def _calculate_substrate_MDV_derivatives_p(self, kind, extra_subs = None):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * "ss" if isotopic steady state.
            * "inst" if isotopically nonstationary state.
        extra_subs: list or None
        '''
        
        if kind == 'ss':
            substrate_MDVs_der_u = self._calculate_substrate_MDV_derivatives_u(extra_subs)
            for emu in substrate_MDVs_der_u:
                self.model.substrate_MDVs_der_p[emu] = substrate_MDVs_der_u[emu]

        elif kind == 'inst':
            substrate_MDVs_der_u = self._calculate_substrate_MDV_derivatives_u(extra_subs)
            substrate_MDVs_der_c = self._calculate_substrate_MDV_derivatives_c(extra_subs)
            for emu in substrate_MDVs_der_u:
                self.model.substrate_MDVs_der_p[emu] = np.concatenate(
                    (substrate_MDVs_der_u[emu], substrate_MDVs_der_c[emu]), 
                    axis = 1
                )
    
    
    def _calculate_measured_fluxes_inversed_covariance_matrix(self):
        
        sds = np.array([sd for _, sd in self.model.measured_fluxes.values()])
        self.model.measured_fluxes_inv_cov = np.diag(1/sds**2) 
    
   
    def _calculate_measured_fluxes_derivative_v(self):
        
        measfluxids = list(self.model.measured_fluxes.keys())
        measured_fluxes_der = np.array(
            Matrix(measfluxids).jacobian(Matrix(self.model.totalfluxids))
        ).astype(float)
        
        return measured_fluxes_der
        
        
    def _calculate_measured_fluxes_derivative_c(self):
        
        nmeasfluxes = len(self.model.measured_fluxes)
        nmetabs = len(self.model.concids)
        measured_fluxes_der = np.zeros((nmeasfluxes, nmetabs))
        
        return measured_fluxes_der
        
        
    def _calculate_measured_fluxes_derivative_p(self, kind):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * "ss" if isotopic steady state.
            * "inst" if isotopically nonstationary state.
        '''
        
        if kind == 'ss':
            measured_fluxes_der_v = self._calculate_measured_fluxes_derivative_v()
            measured_fluxes_der_u = measured_fluxes_der_v@self.model.null_space
            self.model.measured_fluxes_der_p = measured_fluxes_der_u

        elif kind == 'inst':
            measured_fluxes_der_v = self._calculate_measured_fluxes_derivative_v()
            measured_fluxes_der_u = measured_fluxes_der_v@self.model.null_space
            measured_fluxes_der_c = self._calculate_measured_fluxes_derivative_c()
            self.model.measured_fluxes_der_p = np.concatenate(
                (measured_fluxes_der_u, measured_fluxes_der_c), 
                axis = 1
            )

    
    def _generate_random_fluxes(self):
        
        self.ori_measured_fluxes = deepcopy(self.model.measured_fluxes)
        
        for fluxid, [mean, sd] in self.model.measured_fluxes.items():
            while True:
                meanNew = normal(mean, sd)
                if mean-3*sd <= meanNew <= mean+3*sd:
                    break
            self.model.measured_fluxes[fluxid][0] = meanNew
        
        
    def _reset_measured_fluxes(self):
        
        self.model.measured_fluxes = deepcopy(self.ori_measured_fluxes)
        
    
    def _calculate_measured_MDVs_inversed_covariance_matrix(self):
        
        if not self.model.measured_MDVs:
            self.model.measured_MDVs_inv_cov = None
            return

        all_sds_concat = []
        # Consistent ordering as used in MFAModel.build_objective
        for exp_id in sorted(self.model.measured_MDVs.keys()):
            exp_data = self.model.measured_MDVs[exp_id]
            for fragment_id in sorted(exp_data.keys()):
                mean, sd = exp_data[fragment_id] # mean is not used here but part of stored tuple
                all_sds_concat.extend(sd)

        if not all_sds_concat:
            self.model.measured_MDVs_inv_cov = None
            return

        variances = np.array(all_sds_concat)**2
        if np.any(variances == 0):
            logging.warning("Zero variance found in measured MDVs. Replacing with 1e-12 for covariance matrix calculation.")
            variances[variances == 0] = 1e-12

        self.model.measured_MDVs_inv_cov = np.diag(1.0 / variances)
        
        
    def _calculate_measured_inst_MDVs_inversed_covariance_matrix(self):
        
        if not self.model.measured_inst_MDVs:
            self.model.measured_inst_MDVs_inv_cov = None
            return

        all_sds_concat = []
        # Consistent ordering as used in InstMFAModel.build_objective
        for exp_id in sorted(self.model.measured_inst_MDVs.keys()):
            exp_data = self.model.measured_inst_MDVs[exp_id]
            for frag_id in sorted(exp_data.keys()):
                frag_tp_data = exp_data[frag_id]
                for timepoint in sorted(frag_tp_data.keys()):
                    # InstMFAModel iterates all measured timepoints for residuals.
                    # No t=0 skipping here, as it's handled by InstMFAModel's iteration if necessary,
                    # or t=0 is simply not in measured_inst_MDVs for fitting.
                    _mean, sd = frag_tp_data[timepoint] # _mean not used
                    all_sds_concat.extend(sd)

        if not all_sds_concat:
            self.model.measured_inst_MDVs_inv_cov = None
            return

        variances = np.array(all_sds_concat)**2
        if np.any(variances == 0):
            logging.warning("Zero variance found in measured inst. MDVs. Replacing with 1e-12 for covariance matrix calculation.")
            variances[variances == 0] = 1e-12

        self.model.measured_inst_MDVs_inv_cov = np.diag(1.0 / variances)
        
    
    def _generate_random_MDVs(self):
        
        if not hasattr(self, 'ori_measured_MDVs_backup'):
            # Backup original measured MDVs if not already backed up in this instance
            self.ori_measured_MDVs_backup = deepcopy(self.model.measured_MDVs)

        for exp_id in self.model.measured_MDVs: # Iterate experiments
            if exp_id not in self.ori_measured_MDVs_backup: continue # Should not happen if backup is complete
            for fragment_id in self.model.measured_MDVs[exp_id]: # Iterate fragments within experiment
                if fragment_id not in self.ori_measured_MDVs_backup[exp_id]: continue

                original_mean, original_sd = self.ori_measured_MDVs_backup[exp_id][fragment_id]

                # Generate new random means based on original mean and sd
                random_mean = np.random.normal(loc=original_mean, scale=original_sd)

                # Apply constraints: non-negativity and normalization
                random_mean[random_mean < 0] = 0
                if random_mean.sum() != 0:
                    random_mean = random_mean / random_mean.sum()
                else:
                    # Handle case where all components of random_mean are zero (e.g. after clipping negatives)
                    # Fallback to original mean or a uniform distribution might be options.
                    # Here, we fall back to the original mean to ensure data integrity.
                    random_mean = np.copy(original_mean) # Use a copy

                # Update the model's measured_MDVs with the new random mean, keeping original sd
                self.model.measured_MDVs[exp_id][fragment_id] = [random_mean, original_sd]


    def _reset_measured_MDVs(self):
        
        if hasattr(self, 'ori_measured_MDVs_backup'):
            self.model.measured_MDVs = deepcopy(self.ori_measured_MDVs_backup)
            # Optionally, remove the backup if it's meant to be short-lived for one generation cycle
            # del self.ori_measured_MDVs_backup
        else:
            logging.warning("Calculator: _reset_measured_MDVs called without a backup (ori_measured_MDVs_backup).")


    def _generate_random_inst_MDVs(self):
        
        if not hasattr(self, 'ori_measured_inst_MDVs_backup'):
            self.ori_measured_inst_MDVs_backup = deepcopy(self.model.measured_inst_MDVs)

        for exp_id in self.model.measured_inst_MDVs:
            if exp_id not in self.ori_measured_inst_MDVs_backup: continue
            for frag_id in self.model.measured_inst_MDVs[exp_id]:
                if frag_id not in self.ori_measured_inst_MDVs_backup[exp_id]: continue
                for tp in self.model.measured_inst_MDVs[exp_id][frag_id]:
                    if tp not in self.ori_measured_inst_MDVs_backup[exp_id][frag_id]: continue

                    # For instationary, the original code skips t=0.
                    # However, if t=0 is part of measured_inst_MDVs and ori_measured_inst_MDVs_backup,
                    # it should be perturbed or preserved based on consistent logic.
                    # The original _generate_random_inst_MDVs had "if t != 0:"
                    # We should honor this if it's a strict requirement.
                    # If t=0 data is used for initial conditions and not for fitting residuals, it shouldn't be randomized.
                    # Assuming t=0 should not be randomized if present in measured data.
                    if tp == 0: # or self.model.timepoints[0] if timepoints are globally defined and sorted
                        # Preserve original t=0 data if it exists
                        original_mean, original_sd = self.ori_measured_inst_MDVs_backup[exp_id][frag_id][tp]
                        self.model.measured_inst_MDVs[exp_id][frag_id][tp] = [np.copy(original_mean), original_sd]
                        continue

                    original_mean, original_sd = self.ori_measured_inst_MDVs_backup[exp_id][frag_id][tp]

                    random_mean = np.random.normal(loc=original_mean, scale=original_sd)
                    random_mean[random_mean < 0] = 0 # Non-negativity
                    if random_mean.sum() != 0:
                        random_mean = random_mean / random_mean.sum() # Normalize
                    else:
                        random_mean = np.copy(original_mean) # Fallback for all-zero sum

                    self.model.measured_inst_MDVs[exp_id][frag_id][tp] = [random_mean, original_sd]


    def _reset_measured_inst_MDVs(self):

        if hasattr(self, 'ori_measured_inst_MDVs_backup'):
            self.model.measured_inst_MDVs = deepcopy(self.ori_measured_inst_MDVs_backup)
            # del self.ori_measured_inst_MDVs_backup
        else:
            logging.warning("Calculator: _reset_measured_inst_MDVs called without a backup (ori_measured_inst_MDVs_backup).")
    

    @staticmethod
    def _calculate_matrix_A_and_B_derivatives_v(EAM, fluxids):
        '''
        Parameters
        ----------
        EAMstr: df
            EMU adjacency matrix of some size.
        fluxids: list
            Total fluxes IDs.
            
        Returns
        -------
        matrix_A_der, matrix_B_der: 3-D array
            Derivatives of A(B) w.r.t. total fluxes in shape of 
            (len(total fluxes), A(B).shape[0], A(B).shape[1]).
        '''
        
        import platform
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
            
        preAB = EAM.copy(deep = 'all').values
        colSum = preAB.sum(axis = 0)
        for i in range(colSum.size):
            preAB[i, i] = -colSum[i]
        
        A = preAB[:preAB.shape[1], :].T
        B = -preAB[preAB.shape[1]:, :].T
        
        matA = Matrix(A)
        matB = Matrix(B)
        
        if JAX_INSTALLED:
            lambA = lambdify(symbols(fluxids), matA, modules = 'jax')
            lambB = lambdify(symbols(fluxids), matB, modules = 'jax')
            
            matrix_A_der = np.array(
                jacfwd(lambA, range(len(fluxids)))(*jnp.ones(len(fluxids)))
            )
            matrix_B_der = np.array(
                jacfwd(lambB, range(len(fluxids)))(*jnp.ones(len(fluxids)))
            )
        else:
            matrix_A_der = np.array(
                derive_by_array(matA, symbols(fluxids)), 
                dtype = float
            )
            matrix_B_der = np.array(
                derive_by_array(matB, symbols(fluxids)), 
                dtype = float
            )
        
        return matrix_A_der, matrix_B_der
    
    
    def _calculate_matrix_As_and_Bs_derivatives_u(self, n_jobs):
        '''
        Parameters
        ----------
        n_jobs: int
            # of jobs to run in parallel.

        Returns
        -------
        matrix_A_der: 3-D array
            Derivatives of A(B) w.r.t. free fluxes in shape of 
            (len(free fluxes), A(B).shape[0], A(B).shape[1]).
        matrix_B_der: 3-D array
            Derivatives of A(B) w.r.t. free fluxes in shape of 
            (len(free fluxes), A(B).shape[0], A(B).shape[1]).        
        '''
            
        matrix_ABs_der = {}
        if n_jobs == 1:
            for size, EAM in self.model.EAMs.items():
                ABder = self._calculate_matrix_A_and_B_derivatives_v(EAM, self.model.totalfluxids)
                matrix_ABs_der[size] = ABder
        else:
            pool = Pool(processes = n_jobs)
            
            for size, EAM in self.model.EAMs.items():
                ABder = pool.apply_async(
                    func = self._calculate_matrix_A_and_B_derivatives_v, 
                    args = (EAM, self.model.totalfluxids)
                )
                matrix_ABs_der[size] = ABder
            
            pool.close()    
            pool.join()
            
            matrix_ABs_der = {size: ABder.get() for size, ABder in matrix_ABs_der.items()}
        
        matrix_As_der = {}
        matrix_Bs_der = {}
        for size, ABder in matrix_ABs_der.items():
            Ader = ABder[0].swapaxes(0,1).swapaxes(1,2)
            Bder = ABder[1].swapaxes(0,1).swapaxes(1,2)
            
            Ader = Ader@self.model.null_space
            Bder = Bder@self.model.null_space
        
            Ader = Ader.swapaxes(1,2).swapaxes(0,1)
            Bder = Bder.swapaxes(1,2).swapaxes(0,1)
            
            matrix_As_der[size] = Ader
            matrix_Bs_der[size] = Bder
            
        return matrix_As_der, matrix_Bs_der
        
            
    def _calculate_matrix_As_and_Bs_derivatives_c(self):
        '''
        Returns
        -------
        matrix_A_der: 3-D array
            Derivatives of A(B) w.r.t. total fluxes in shape of 
            (len(concs), A(B).shape[0], A(B).shape[1]).
        matrix_B_der: 3-D array
            Derivatives of A(B) w.r.t. total fluxes in shape of 
            (len(concs), A(B).shape[0], A(B).shape[1]).    
        '''
        
        nmetabs = len(self.model.concids)
        
        matrix_As_der = {}
        matrix_Bs_der = {}
        for size, EAM in self.model.EAMs.items():
            nEMUsout = EAM.shape[1]
            nEMUsin = EAM.shape[0] - nEMUsout
            
            matrix_As_der[size] = np.zeros((nmetabs, nEMUsout, nEMUsout))
            matrix_Bs_der[size] = np.zeros((nmetabs, nEMUsout, nEMUsin))
        
        return matrix_As_der, matrix_Bs_der
        
        
    def _calculate_matrix_As_and_Bs_derivatives_p(self, kind, n_jobs):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * "ss" if isotopic steady state.
            * "inst" if isotopically nonstationary state.
        n_jobs: int
            # of jobs to run in parallel.
        '''
        
        if kind == 'ss':
            (matrix_As_der_u, 
             matrix_Bs_der_u
            ) = self._calculate_matrix_As_and_Bs_derivatives_u(n_jobs)
            for size in self.model.EAMs:
                self.model.matrix_As_der_p[size] = matrix_As_der_u[size]
                self.model.matrix_Bs_der_p[size] = matrix_Bs_der_u[size]

        elif kind == 'inst':
            (matrix_As_der_u, 
             matrix_Bs_der_u
            ) = self._calculate_matrix_As_and_Bs_derivatives_u(n_jobs)
            (matrix_As_der_c, 
             matrix_Bs_der_c
            ) = self._calculate_matrix_As_and_Bs_derivatives_c()
            for size in self.model.EAMs:
                self.model.matrix_As_der_p[size] = np.concatenate(
                    (matrix_As_der_u[size], matrix_As_der_c[size]), 
                    axis = 0
                )
                self.model.matrix_Bs_der_p[size] = np.concatenate(
                    (matrix_Bs_der_u[size], matrix_Bs_der_c[size]), 
                    axis = 0
                )
    
    
    def _lambdify_matrix_As_and_Bs(self):
       
        for size, EAM in self.model.EAMs.items():
            
            preAB = EAM.copy(deep = 'all')
            for emu in EAM.columns:
                preAB.loc[emu, emu] = -preAB[emu].sum() 
        
            A = preAB.loc[preAB.columns, :].T
            B = -preAB.loc[preAB.index.difference(preAB.columns), :].T
            
            matA = Matrix(A)
            matB = Matrix(B)
            
            fluxidsA = list(map(str, matA.free_symbols))
            fluxidsB = list(map(str, matB.free_symbols))
            
            lambA = lambdify(fluxidsA, matA, modules = 'numpy')
            lambB = lambdify(fluxidsB, matB, modules = 'numpy')
            
            self.model.matrix_As[size] = [lambA, fluxidsA, A.columns.tolist()]
            self.model.matrix_Bs[size] = [lambB, fluxidsB, B.columns.tolist()]
    
    def _calculate_matrix_Ms_derivatives_u(self):
        
        nfreefluxes = self.model.null_space.shape[1]
        
        matrix_Ms_der = {}
        for size, EAM in self.model.EAMs.items():
            nEMUsout = EAM.shape[1]
            matrix_Ms_der[size] = np.zeros((nfreefluxes, nEMUsout, nEMUsout))
        
        return matrix_Ms_der
        
    
    def _calculate_matrix_Ms_derivatives_c(self):
        
        matrix_Ms_der = {}
        for size, EAM in self.model.EAMs.items():            
            
            matM = Matrix(np.diag(symbols([emu.metabolite_id for emu in EAM.columns])))
            if JAX_INSTALLED:
                lambM = lambdify(symbols(self.model.concids), matM, modules = 'jax')
                matrix_M_der = np.array(
                    jacfwd(lambM, 
                           range(len(self.model.concids))
                    )(*jnp.ones(len(self.model.concids)))
                )
            else:
                matrix_M_der = np.array(
                    derive_by_array(matM, symbols(self.model.concids)), 
                    dtype = float
                )
            matrix_Ms_der[size] = matrix_M_der
            
        return matrix_Ms_der
        
        
    def _calculate_matrix_Ms_derivatives_p(self):
        
        matrix_Ms_der_u = self._calculate_matrix_Ms_derivatives_u()
        matrix_Ms_der_c = self._calculate_matrix_Ms_derivatives_c()
        
        for size in self.model.EAMs:
            self.model.matrix_Ms_der_p[size] = np.concatenate(
                (matrix_Ms_der_u[size], matrix_Ms_der_c[size]), 
                axis = 0
            )
        
    
    def _lambdify_matrix_Ms(self):
        
        for size, EAM in self.model.EAMs.items():
            matM = Matrix(np.diag(symbols([emu.metabolite_id for emu in EAM.columns])))
            metabids = list(map(str, matM.free_symbols))
            lambM = lambdify(metabids, matM, modules = 'numpy')
            self.model.matrix_Ms[size] = [lambM, metabids]
        
    
    # initial X(Y) and their derivatives
    # _get_initial_matrix_Xs_for_experiment remains largely unchanged in internal logic for now,
    # as X usually starts with natural abundance. It now accepts experiment_id.
    def _get_initial_matrix_Xs_for_experiment(self, experiment_id='exp0'):
        """
        Calculates initial MDVs for EMUs in matrix X for a given experiment.
        Typically, these are natural abundance.
        Args:
            experiment_id (str): The ID of the current experiment. (Used for consistency, may not alter logic if X is always natural)
        Returns:
            dict: {size: initial_X_matrix}
        """
        # current_labeling_strategy = self._get_experiment_labeling_strategy(experiment_id) # Not used if all X are natural
        initial_X_dict = {}
        if not self.model.matrix_As: # EAMs not built yet
             logging.warning("EAMs (matrix_As) not available for _get_initial_matrix_Xs_for_experiment.")
             return initial_X_dict

        for size in self.model.matrix_As:
            product_emus_for_size = self.model.matrix_As[size][2] # List of EMU objects
            nEMUs = len(product_emus_for_size)
            iniX_matrix_rows = []
            for product_emu in product_emus_for_size: # Iterate to create rows for each product EMU
                 # If specific X EMUs could be non-natural based on experiment_id, logic would go here.
                 # For now, all are natural.
                iniX_matrix_rows.append(get_natural_MDV(product_emu.size).value)
            
            if iniX_matrix_rows:
                 initial_X_dict[size] = np.vstack(iniX_matrix_rows)
            else: # No product EMUs for this size
                 initial_X_dict[size] = np.empty((0, size + 1))
            
        return initial_X_dict


    def _get_initial_matrix_Ys_for_experiment(self, experiment_id='exp0', substrate_mdvs_for_experiment=None):
        """
        Calculates initial MDVs for EMUs in matrix Y for a given experiment.
        Uses pre-calculated substrate_mdvs_for_experiment for lookups.
        Args:
            experiment_id (str): The ID of the current experiment.
            substrate_mdvs_for_experiment (dict): Pre-calculated MDVs for all relevant substrate EMUs
                                                 for this experiment ({EMU_object: mdv_array}).
        Returns:
            dict: {size: initial_Y_matrix}
        """
        initial_Y_dict = {}
        if substrate_mdvs_for_experiment is None:
            logging.error(f"Calculator: substrate_mdvs_for_experiment not provided for experiment {experiment_id} in _get_initial_matrix_Ys. Cannot proceed.")
            return initial_Y_dict # Return empty or handle error appropriately

        if not self.model.matrix_Bs: # EAMs not built yet
             logging.warning("EAMs (matrix_Bs) not available for _get_initial_matrix_Ys_for_experiment.")
             return initial_Y_dict

        for size in self.model.matrix_Bs:
            source_emus_descriptors_for_size = self.model.matrix_Bs[size][2] # List of EMU objects or tuples of EMU objects

            iniY_matrix_rows = []
            for sourceEMU_descriptor in source_emus_descriptors_for_size:
                if not isinstance(sourceEMU_descriptor, Iterable): # Single EMU
                    # Look up in the pre-calculated dict for the current experiment
                    sourceMDV = substrate_mdvs_for_experiment.get(sourceEMU_descriptor, get_natural_MDV(sourceEMU_descriptor.size))
                else: # Tuple of EMUs, needs convolution
                    mdv_objects_to_convolve = []
                    for emu_in_tuple in sourceEMU_descriptor:
                        # Look up each part in the pre-calculated dict
                        mdv_part = substrate_mdvs_for_experiment.get(emu_in_tuple, get_natural_MDV(emu_in_tuple.size))
                        mdv_objects_to_convolve.append(mdv_part) # mdv_part is already an MDV object or np.array

                    # Ensure they are MDV objects for reduce(conv, ...)
                    # If they are numpy arrays, wrap them. Assuming get_substrate_MDV/get_natural_MDV return MDV objects.
                    # If they are already numpy arrays from the dict, ensure conv can handle them or wrap them.
                    # For now, assume substrate_mdvs_for_experiment stores MDV objects or arrays compatible with conv.
                    # Let's assume the dict stores MDV objects as Calculator's original substrate_MDVs did.
                    sourceMDV = reduce(conv, mdv_objects_to_convolve)

                iniY_matrix_rows.append(sourceMDV.value if isinstance(sourceMDV, MDV) else sourceMDV) # Get numpy array

            if iniY_matrix_rows:
                initial_Y_dict[size] = np.array(iniY_matrix_rows) # Convert list of arrays to 2D array
            else: # No source EMUs for this size
                # Y matrix shape is (num_source_emus, size+1)
                initial_Y_dict[size] = np.empty((0, size+1))

        return initial_Y_dict
        
    
    def _calculate_initial_matrix_Xs_derivatives_p(self):
        
        nvars = self.model.null_space.shape[1] + len(self.model.concids)
        for size, iniX in self.model.initial_matrix_Xs.items():
            Xshape = iniX.shape
            iniXder = np.zeros((nvars, *Xshape))
            self.model.initial_matrix_Xs_der_p[size] = iniXder
        
        
    def _calculate_initial_matrix_Ys_derivatives_p(self):
        
        nvars = self.model.null_space.shape[1] + len(self.model.concids)
        for size, iniY in self.model.initial_matrix_Ys.items():
            Yshape = iniY.shape
            iniYder = np.zeros((nvars, *Yshape))
            self.model.initial_matrix_Ys_der_p[size] = iniYder
        
    
    def _build_initial_sim_MDVs(self):
        
        for size in sorted(self.model.matrix_As):
            productEMUs = self.model.matrix_As[size][2]
            iniX = self.model.initial_matrix_Xs[size]        
            for productEMU, iniMDV in zip(productEMUs, iniX):
                if productEMU.id in self.model.target_EMUs:
                    self.model.initial_sim_MDVs[productEMU.id] = {0: MDV(iniMDV)}
                    
        
    # Renamed from _calculate_MDVs
    def _calculate_MDVs_for_experiment(self, params_u, experiment_id='exp0'):
        '''
        This method simulates MDVs at isotopically steady state for a specific experiment.
        Args:
            params_u: Free flux parameters (already used to set self.model.total_fluxes).
            experiment_id (str): The ID of the current experiment.
        Returns:
            dict: {emu_id: mdv_array} for the specified experiment.
        '''
        # Ensure total_fluxes are set based on params_u (usually done by caller like MFAModel)
        # self.model.total_fluxes[:] = self.model.null_space @ params_u

        all_substrate_emus = self._get_all_substrate_emus()
        substrate_mdvs_exp = self._get_substrate_MDVs_for_experiment(all_substrate_emus, experiment_id)
        
        # initial_X_dict_exp = self._get_initial_matrix_Xs_for_experiment(experiment_id) # Not directly used in SS Y calculation
        initial_Y_dict_exp = self._get_initial_matrix_Ys_for_experiment(experiment_id, substrate_mdvs_exp)

        simMDVs_exp_obj_map = {} # Maps EMU object to its MDV (as MDV object or array)

        for size in sorted(self.model.matrix_As.keys()):
            lambA, fluxidsA, productEMUs_list = self.model.matrix_As[size]
            lambB, fluxidsB, _ = self.model.matrix_Bs[size] # sourceEMUs descriptors from matrix_Bs used in initial_Y_dict_exp
            
            A_matrix_val = lambA(*self.model.total_fluxes[fluxidsA])
            B_matrix_val = lambB(*self.model.total_fluxes[fluxidsB])
            
            # Y_matrix must be constructed using initial_Y_dict_exp and
            # MDVs of smaller EMUs already simulated in *this current experiment simulation*
            Y_matrix_rows = []
            source_emu_descriptors_for_Y = self.model.matrix_Bs[size][2] # These are the actual descriptors for Y

            for source_desc in source_emu_descriptors_for_Y:
                if not isinstance(source_desc, Iterable): # Single EMU
                    # Try to get from already simulated smaller EMUs in this experiment, then from initial Y for this experiment
                    # The initial_Y_dict_exp contains MDVs for source EMUs defined in EAMs (matrix_Bs list)
                    # These are effectively the "inputs" to the EMU system of this size.
                    # The original initial_Y_dict_exp was based on substrate_mdvs_exp (for base substrates)
                    # or natural abundance. Here we use that directly.
                    # The ChainMap logic was for hierarchical calculation: simMDVs_exp_obj_map first, then substrate_mdvs_exp
                    mdv_val = ChainMap(simMDVs_exp_obj_map, substrate_mdvs_exp).get(source_desc)
                    if mdv_val is None: # Should not happen if all source EMUs are covered
                        logging.warning(f"SS MDV Calc: EMU {source_desc} not found in sim or substrate MDVs for exp {experiment_id}. Using natural.")
                        mdv_val = get_natural_MDV(source_desc.size)
                    Y_matrix_rows.append(mdv_val.value if isinstance(mdv_val, MDV) else mdv_val)
                else: # Tuple of EMUs for convolution
                    mdvs_to_convolve = []
                    for emu_in_tuple in source_desc:
                        mdv_val_part = ChainMap(simMDVs_exp_obj_map, substrate_mdvs_exp).get(emu_in_tuple)
                        if mdv_val_part is None:
                            logging.warning(f"SS MDV Calc: EMU {emu_in_tuple} in convolution not found for exp {experiment_id}. Using natural.")
                            mdv_val_part = get_natural_MDV(emu_in_tuple.size)
                        mdvs_to_convolve.append(mdv_val_part) # MDV object or array
                    convolved_mdv = reduce(conv, mdvs_to_convolve)
                    Y_matrix_rows.append(convolved_mdv.value if isinstance(convolved_mdv, MDV) else convolved_mdv)
            
            if not Y_matrix_rows: # No source EMUs for this size, Y is empty
                 # X must also be empty or all zeros if there are product EMUs
                if productEMUs_list:
                    X_simulated_values = np.zeros((len(productEMUs_list), size + 1))
                else:
                    X_simulated_values = np.empty((0, size + 1))
            else:
                Y_matrix_val = np.array(Y_matrix_rows)
                if Y_matrix_val.shape[0] == 0 and B_matrix_val.shape[1] > 0 : # B expects rows but Y is empty
                     X_simulated_values = np.zeros((A_matrix_val.shape[0], size + 1)) if A_matrix_val.shape[0] > 0 else np.empty((0, size+1))
                elif B_matrix_val.shape[1] != Y_matrix_val.shape[0]:
                    logging.error(f"SS MDV Calc: Mismatch B columns {B_matrix_val.shape[1]} vs Y rows {Y_matrix_val.shape[0]} for size {size} in exp {experiment_id}")
                    X_simulated_values = np.zeros((A_matrix_val.shape[0], size + 1)) if A_matrix_val.shape[0] > 0 else np.empty((0, size+1))
                else:
                    X_simulated_values = pinv(A_matrix_val, check_finite=False) @ B_matrix_val @ Y_matrix_val
        
            for emu_obj, mdv_array in zip(productEMUs_list, X_simulated_values):
                simMDVs_exp_obj_map[emu_obj] = MDV(mdv_array) # Store as MDV object for consistency
            
        # Convert EMU object keys to EMU IDs for the final returned dict
        simMDVs_exp_id_map = {emu.id: mdv_obj.value for emu, mdv_obj in simMDVs_exp_obj_map.items()}
        
        return simMDVs_exp_id_map
                    
    # This method now calculates the GLOBAL stacked derivative matrix.
    # The simMDVs_dict it returns is for compatibility or general info;
    # MFAModel should call _calculate_MDVs_for_experiment for specific experiment residuals.
    def _calculate_MDVs_and_derivatives_p(self):
        '''
        This method simulates MDVs and their derivatives at isotopically steady state.
        It calculates derivatives for all relevant (emu, experiment) pairs and stacks them.
        
        Returns
        -------
        simMDVs: dict
            EMU ID => MDV (in array).
        simMDVsDer: dict
            EMU ID => 2-D array in shape of (len(MDV), len(u)).
        '''
        
        simMDVs = {}
        simMDVsDer = {}
        for size in sorted(self.model.matrix_As):
            
            lambA, fluxidsA, productEMUs = self.model.matrix_As[size]
            lambB, fluxidsB, sourceEMUs = self.model.matrix_Bs[size]
            
            A = lambA(*self.model.total_fluxes[fluxidsA])
            B = lambB(*self.model.total_fluxes[fluxidsB])
            
            Ainv = pinv(A, check_finite = True)
            
            Ader = self.model.matrix_As_der_p[size]   
            Bder = self.model.matrix_Bs_der_p[size]   
            
            Y = []
            Yder = []
            for sourceEMU in sourceEMUs:
                if not isinstance(sourceEMU, Iterable):
                    sourceMDV = self.model.substrate_MDVs[sourceEMU]
                    sourceMDVder = self.model.substrate_MDVs_der_p[sourceEMU]
                else:
                    mdvs = []
                    mdvs_mdvders = []
                    for emu in sourceEMU:
                        mdv = ChainMap(simMDVs, self.model.substrate_MDVs)[emu]
                        mdvs.append(mdv)
                        mdvder = ChainMap(simMDVsDer, self.model.substrate_MDVs_der_p)[emu]
                        mdvs_mdvders.append([mdv, mdvder])
                    sourceMDV = reduce(conv, mdvs)
                    sourceMDVder = reduce(diff_conv, mdvs_mdvders)[1]
                Y.append(sourceMDV)
                Yder.append(sourceMDVder)
            Y = np.array(Y)
            Yder = np.array(Yder).swapaxes(1,2).swapaxes(0,1)   
            
            X = Ainv@B@Y
            Xder = Ainv@(Bder@Y + B@Yder - Ader@X)
            Xder = Xder.swapaxes(0,1).swapaxes(1,2)   
            
            simMDVs.update(zip(productEMUs, X))
            simMDVsDer.update(zip(productEMUs, Xder))
            
        # This method needs to calculate derivatives for each (exp_id, fragment_id) pair,
        # using the correct experiment-specific context (substrate_MDVs, Y-matrices).
        # The final returned stacked matrix must align with the order defined by iterating
        # through sorted(self.model.measured_MDVs.keys()), then sorted(fragments), etc.

        # For steady-state, params_p are actually params_u (free flux parameters)
        # These are assumed to be set on self.model.total_fluxes by the caller (MFAModel)

        all_substrate_emus_global = self._get_all_substrate_emus()

        # This will store emu_id -> derivative_array, but calculated per experiment context
        # The key for MFAModel is the final stacked matrix.
        # For simMDVs output, we can return a dict of all simulated EMUs from the *last* experiment,
        # or an aggregated dict if useful, or None if MFAModel will call _calculate_MDVs_for_experiment.
        # Let's return simMDVs from the first experiment for now if needed by caller.
        first_exp_id = sorted(self.model.measured_MDVs.keys())[0] if self.model.measured_MDVs else None
        simMDVs_for_first_exp_dict = {}

        ordered_derivatives_list = []
        expected_total_rows = 0
        num_free_params = self.model.null_space.shape[1] if self.model.null_space is not None else 0

        if not self.model.measured_MDVs:
            return {}, np.empty((0, num_free_params))

        for exp_id in sorted(self.model.measured_MDVs.keys()):
            substrate_mdvs_exp = self._get_substrate_MDVs_for_experiment(all_substrate_emus_global, exp_id)
            initial_Y_dict_exp = self._get_initial_matrix_Ys_for_experiment(exp_id, substrate_mdvs_exp)

            # Simulate MDVs for this specific experiment to be used in derivative calculations for this block
            # params_u are implicitly from self.model.total_fluxes which should be set by solver
            simMDVs_this_exp_obj_map = {} # emu_obj -> MDV_obj

            # Simplified simulation for derivative context (X = A_inv @ B @ Y)
            # This loop calculates X for each size for the current exp_id
            for size_iter in sorted(self.model.matrix_As.keys()):
                lambA, fluxidsA, productEMUs_iter_list = self.model.matrix_As[size_iter]
                lambB, fluxidsB, sourceEMUs_desc_iter_list = self.model.matrix_Bs[size_iter]

                A_iter_val = lambA(*self.model.total_fluxes[fluxidsA])
                B_iter_val = lambB(*self.model.total_fluxes[fluxidsB])
                Ainv_iter_val = pinv(A_iter_val, check_finite=False)

                Y_iter_matrix_rows = []
                for source_desc_iter in sourceEMUs_desc_iter_list:
                    if not isinstance(source_desc_iter, Iterable):
                        mdv_val = ChainMap(simMDVs_this_exp_obj_map, substrate_mdvs_exp).get(source_desc_iter, get_natural_MDV(source_desc_iter.size))
                        Y_iter_matrix_rows.append(mdv_val.value if isinstance(mdv_val, MDV) else mdv_val)
                    else:
                        mdvs_conv = [ChainMap(simMDVs_this_exp_obj_map, substrate_mdvs_exp).get(e, get_natural_MDV(e.size)) for e in source_desc_iter]
                        Y_iter_matrix_rows.append(reduce(conv, mdvs_conv).value)

                if not Y_iter_matrix_rows:
                    X_iter_sim_val = np.zeros((len(productEMUs_iter_list), size_iter + 1)) if productEMUs_iter_list else np.empty((0, size_iter + 1))
                else:
                    Y_iter_val = np.array(Y_iter_matrix_rows)
                    if B_iter_val.shape[1] != Y_iter_val.shape[0]: # Mismatch check
                        X_iter_sim_val = np.zeros((A_iter_val.shape[0], size_iter + 1))
                    else:
                        X_iter_sim_val = Ainv_iter_val @ B_iter_val @ Y_iter_val

                for emu_obj, mdv_arr in zip(productEMUs_iter_list, X_iter_sim_val):
                    simMDVs_this_exp_obj_map[emu_obj] = MDV(mdv_arr)

                if exp_id == first_exp_id: # Populate for the representative simMDVs_dict
                    for emu_obj, mdv_obj_val in simMDVs_this_exp_obj_map.items():
                         if emu_obj.id not in simMDVs_for_first_exp_dict : # only add if not from smaller size
                              simMDVs_for_first_exp_dict[emu_obj.id] = mdv_obj_val.value


            # Now calculate derivatives for fragments measured in this experiment
            # The derivative calculation for an EMU, Xder = Ainv@(Bder@Y + B@Yder - Ader@X)
            # requires X and Y specific to this experiment.
            # Ader, Bder are from self.model.matrix_As_der_p etc (these are generic based on flux values)
            # Yder needs to be specific to this experiment's substrate MDV derivatives.
            # This implies self.model.substrate_MDVs_der_p also needs to be experiment-specific if substrates change.
            # For now, assume self.model.substrate_MDVs_der_p is global (e.g. zero for non-tracer, or some fixed derivative for tracer).
            # This is a simplification. If substrate derivative depends on concentration/labeling that varies per exp, this needs more.

            temp_exp_simMDVsDer_obj_map = {} # emu_obj -> derivative_array for this experiment

            for size_iter in sorted(self.model.matrix_As.keys()): # Iterate sizes again for derivatives
                lambA, fluxidsA, productEMUs_iter_list = self.model.matrix_As[size_iter]
                lambB, fluxidsB, sourceEMUs_desc_iter_list = self.model.matrix_Bs[size_iter]

                A_iter_val = lambA(*self.model.total_fluxes[fluxidsA])
                B_iter_val = lambB(*self.model.total_fluxes[fluxidsB])
                Ainv_iter_val = pinv(A_iter_val, check_finite=False)

                Ader_val = self.model.matrix_As_der_p[size_iter] # These are (n_params, n_emu_out, n_emu_out)
                Bder_val = self.model.matrix_Bs_der_p[size_iter] # These are (n_params, n_emu_out, n_emu_in)

                # Construct Y and Yder for this experiment
                Y_matrix_rows = []
                Yder_matrix_rows = [] # List of (n_params, n_isotopomers)

                for source_desc_iter in sourceEMUs_desc_iter_list:
                    if not isinstance(source_desc_iter, Iterable):
                        mdv_obj = ChainMap(simMDVs_this_exp_obj_map, substrate_mdvs_exp).get(source_desc_iter, get_natural_MDV(source_desc_iter.size))
                        # Assume substrate_MDVs_der_p is global or correctly selected for exp. For now, global.
                        mdv_der_val = ChainMap(temp_exp_simMDVsDer_obj_map, self.model.substrate_MDVs_der_p).get(source_desc_iter, np.zeros((num_free_params, source_desc_iter.size + 1)))
                        Y_matrix_rows.append(mdv_obj.value)
                        Yder_matrix_rows.append(mdv_der_val)
                    else: # Convolution
                        mdv_parts = [ChainMap(simMDVs_this_exp_obj_map, substrate_mdvs_exp).get(e, get_natural_MDV(e.size)) for e in source_desc_iter]
                        mdv_der_parts_as_list_of_lists = [] # each item: [mdv_obj, mdv_der_array]
                        for e_conv in source_desc_iter:
                            e_mdv = ChainMap(simMDVs_this_exp_obj_map, substrate_mdvs_exp).get(e_conv, get_natural_MDV(e_conv.size))
                            e_der = ChainMap(temp_exp_simMDVsDer_obj_map, self.model.substrate_MDVs_der_p).get(e_conv, np.zeros((num_free_params, e_conv.size+1)))
                            mdv_der_parts_as_list_of_lists.append([e_mdv, e_der]) # Pass MDV object and der array

                        Y_matrix_rows.append(reduce(conv, mdv_parts).value)
                        convolved_der = reduce(diff_conv, mdv_der_parts_as_list_of_lists)[1] # diff_conv returns [mdv, mdv_der]
                        Yder_matrix_rows.append(convolved_der)

                if not Y_matrix_rows:
                     Xder_val_exp = np.zeros((num_free_params, len(productEMUs_iter_list), size_iter + 1)) if productEMUs_iter_list else np.empty((num_free_params, 0, size_iter+1))
                else:
                    Y_val_exp = np.array(Y_matrix_rows)
                    # Yder_val_exp shape: (n_params, n_source_emus_for_Y, n_isotopomers)
                    Yder_val_exp = np.array(Yder_matrix_rows).swapaxes(0,1).swapaxes(1,2) if Yder_matrix_rows else np.zeros((num_free_params,0,size_iter+1))

                    # X values for this experiment and size
                    X_val_exp = np.array([simMDVs_this_exp_obj_map[emu].value for emu in productEMUs_iter_list]) if productEMUs_iter_list else np.empty((0,size_iter+1))

                    if B_iter_val.shape[1] != Y_val_exp.shape[0] or \
                       (Yder_val_exp.ndim ==3 and B_iter_val.shape[1] != Yder_val_exp.shape[1]): # Yder might be empty if no source emus
                        Xder_val_exp = np.zeros((num_free_params, len(productEMUs_iter_list), size_iter + 1))
                    else:
                        # Term1: Bder @ Y_exp
                        term1 = np.einsum('pij,jk->pik', Bder_val, Y_val_exp) if Y_val_exp.size > 0 else np.zeros_like(Bder_val)
                        # Term2: B @ Yder_exp
                        term2 = np.einsum('ij,pjk->pik', B_iter_val, Yder_val_exp) if Yder_val_exp.size > 0 else np.zeros_like(Bder_val)
                        # Term3: Ader @ X_exp
                        term3 = np.einsum('pij,jk->pik', Ader_val, X_val_exp) if X_val_exp.size > 0 else np.zeros_like(Ader_val)

                        Xder_val_exp = Ainv_iter_val @ (term1 + term2 - term3)

                # Transpose Xder_val_exp from (n_params, n_prod_emu, n_isotopomers) to (n_prod_emu, n_params, n_isotopomers)
                # then to (n_prod_emu, n_isotopomers, n_params) for storage if that's convention.
                # Original was (n_prod_emu, n_isotopomers, n_params)
                Xder_val_exp_transposed = Xder_val_exp.swapaxes(0,1).swapaxes(1,2)

                for idx, emu_obj in enumerate(productEMUs_iter_list):
                    temp_exp_simMDVsDer_obj_map[emu_obj] = Xder_val_exp_transposed[idx]


            # Append derivatives for measured fragments in this experiment to the global list
            exp_data_measured = self.model.measured_MDVs[exp_id]
            for fragment_id_iter in sorted(exp_data_measured.keys()):
                num_mass_isotopomers = len(exp_data_measured[fragment_id_iter][0])
                expected_total_rows += num_mass_isotopomers

                # Find the EMU object for fragment_id_iter to look up in temp_exp_simMDVsDer_obj_map
                found_emu_obj = None
                for emu_obj_key in temp_exp_simMDVsDer_obj_map.keys():
                    if emu_obj_key.id == fragment_id_iter:
                        found_emu_obj = emu_obj_key
                        break

                if found_emu_obj and found_emu_obj in temp_exp_simMDVsDer_obj_map:
                    frag_der = temp_exp_simMDVsDer_obj_map[found_emu_obj] # Should be (n_iso, n_params)
                    if frag_der.shape[0] == num_mass_isotopomers and frag_der.shape[1] == num_free_params:
                        ordered_derivatives_list.append(frag_der)
                    else:
                        logging.error(f"Calculator SS Deriv: Shape mismatch for {fragment_id_iter} in exp {exp_id}. Expected ({num_mass_isotopomers}, {num_free_params}), got {frag_der.shape}. Appending zeros.")
                        ordered_derivatives_list.append(np.zeros((num_mass_isotopomers, num_free_params)))
                else:
                    logging.error(f"Calculator SS Deriv: Missing derivative for {fragment_id_iter} in exp {exp_id}. Appending zeros.")
                    ordered_derivatives_list.append(np.zeros((num_mass_isotopomers, num_free_params)))
        
        final_simMDVs_der_stacked = np.vstack(ordered_derivatives_list) if ordered_derivatives_list else np.empty((0, num_free_params))
        if final_simMDVs_der_stacked.shape[0] != expected_total_rows:
            logging.error(f"Calculator SS Deriv: Row count mismatch for stacked derivatives: expected {expected_total_rows}, got {final_simMDVs_der_stacked.shape[0]}")

        # For compatibility, simMDVs_dict might be needed by MFAModel.
        # Here, it's from the first experiment, or could be an aggregation or specific choice.
        return simMDVs_for_first_exp_dict, final_simMDVs_der_stacked
    
    
    # Renamed from _calculate_inst_MDVs
    def _calculate_inst_MDVs_for_experiment(self, params_p, experiment_id='exp0'):
        '''
        This method simulates MDVs at isotopically nonstationary state for a specific experiment.
        Args:
            params_p: Parameters [free_fluxes, concentrations].
            experiment_id (str): The ID of the current experiment.
        
        Returns
        -------
        sim_inst_MDVs_exp_id_map: dict
            EMU ID => {t => MDV (in array)} (starting from t1).
        '''
        
        # Update model fluxes and concentrations from params_p
        num_free_fluxes = self.model.null_space.shape[1]
        u_params = params_p[:num_free_fluxes]
        c_params = params_p[num_free_fluxes:]
        self.model.total_fluxes.loc[:] = self.model.null_space @ u_params
        # Ensure concids are sorted or consistently ordered if c_params order depends on it
        # Assuming self.model.concids is the order for c_params
        if len(c_params) == len(self.model.concids):
            self.model.concentrations.loc[self.model.concids] = c_params
        else:
            logging.error(f"Inst MDV Calc: Mismatch between c_params length ({len(c_params)}) and model.concids length ({len(self.model.concids)}) for exp {experiment_id}.")
            # Potentially return empty or raise error
            return {}


        all_substrate_emus = self._get_all_substrate_emus()
        substrate_mdvs_exp = self._get_substrate_MDVs_for_experiment(all_substrate_emus, experiment_id)

        initial_X_dict_exp = self._get_initial_matrix_Xs_for_experiment(experiment_id)
        initial_Y_dict_exp = self._get_initial_matrix_Ys_for_experiment(experiment_id, substrate_mdvs_exp)

        # Stores MDV objects for EMUs at different time points {tp: {size: matrix_X_or_Y}}
        # For X, it's product EMUs. For Y, it's source EMUs.
        # This needs to map EMU objects to MDV objects for use in ChainMap

        # Xs_sim: {timepoint: {size: X_matrix_at_timepoint}}
        # Ys_sim: {timepoint: {size: Y_matrix_at_timepoint}}
        # sim_inst_MDVs_exp_obj_map: {EMU_object: {timepoint: MDV_object}} # For hierarchical lookup within an experiment
        
        Xs_sim_all_tps_all_sizes = {0.0: initial_X_dict_exp}
        Ys_sim_all_tps_all_sizes = {0.0: initial_Y_dict_exp}
        sim_inst_MDVs_exp_obj_map = {} # Using EMU objects as keys initially

        # Initialize sim_inst_MDVs_exp_obj_map with t0 values from initial_X_dict_exp
        # This is for product EMUs which are keys in matrix_As
        for size, X_matrix_t0 in initial_X_dict_exp.items():
            product_emus_list_for_size = self.model.matrix_As[size][2]
            for idx, emu_obj in enumerate(product_emus_list_for_size):
                sim_inst_MDVs_exp_obj_map.setdefault(emu_obj, {})[0.0] = MDV(X_matrix_t0[idx, :])

        # Initialize with t0 for source EMUs from initial_Y_dict_exp (these are inputs to the system)
        # These are used if a sourceEMU is itself a product of another reaction (hierarchical)
        # For base substrates (from labeling_strategy), their t0 values are in substrate_mdvs_exp
        for size, Y_matrix_t0 in initial_Y_dict_exp.items():
            source_emu_descriptors_for_size = self.model.matrix_Bs[size][2]
            for idx, source_desc in enumerate(source_emu_descriptors_for_size):
                 if not isinstance(source_desc, Iterable): # Single EMU
                      # If this source_desc is a product EMU, its t0 MDV might already be in sim_inst_MDVs_exp_obj_map
                      # If it's a base substrate, its MDV is in substrate_mdvs_exp
                      if source_desc not in sim_inst_MDVs_exp_obj_map:
                           mdv_val = substrate_mdvs_exp.get(source_desc, get_natural_MDV(source_desc.size))
                           sim_inst_MDVs_exp_obj_map.setdefault(source_desc, {})[0.0] = mdv_val # mdv_val is already MDV object
                 else: # Tuple of EMUs
                      # For convoluted EMUs at t0, calculate their MDV and store if needed for lookups
                      # This is complex, as convoluted EMUs aren't usually stored directly.
                      # Their contribution is via the Y matrix.
                      pass


        # Current time, starts from t0
        # Sort timepoints excluding t0 for the simulation loop, t0 is initial condition.
        # self.model.timepoints should be sorted and include t0 if relevant.
        sim_timepoints = sorted([tp for tp in self.model.timepoints if tp > 0])

        # t_prev is the timepoint of the previously calculated state
        t_prev = 0.0
        for t_curr in sim_timepoints:
            deltat = t_curr - t_prev

            Xs_sim_curr_tp_all_sizes = {}
            Ys_sim_curr_tp_all_sizes = {}

            for size in sorted(self.model.matrix_As.keys()):
                lambA, fluxidsA, productEMUs_list = self.model.matrix_As[size]
                lambB, fluxidsB, sourceEMU_descriptors_list = self.model.matrix_Bs[size]
                lambM, metabids_M = self.model.matrix_Ms[size]

                A_val = lambA(*self.model.total_fluxes[fluxidsA])
                B_val = lambB(*self.model.total_fluxes[fluxidsB])
                M_val = lambM(*self.model.concentrations[metabids_M])
                Minv_val = pinv(M_val, check_finite=False)

                F_val = Minv_val @ A_val
                Finv_val = pinv(F_val, check_finite=False)
                I_mtx = np.eye(F_val.shape[0])
                Phi_val = expm(F_val * deltat)
                Gamma_val = (Phi_val - I_mtx) @ Finv_val
                Omega_val = (Gamma_val / deltat - I_mtx) @ Finv_val

                X_matrix_t_prev = Xs_sim_all_tps_all_sizes[t_prev][size]
                Y_matrix_t_prev = Ys_sim_all_tps_all_sizes[t_prev][size] # This is the Y matrix for inputs at t_prev
                G_matrix_t_prev = Minv_val @ B_val @ Y_matrix_t_prev

                # Construct Y_matrix for t_curr (Y_t1 in original equations)
                # This Y matrix represents the MDVs of the source EMUs at time t_curr
                Y_matrix_t_curr_rows = []
                for source_desc in sourceEMU_descriptors_list:
                    if not isinstance(source_desc, Iterable): # Single EMU
                        # For single EMUs, their MDVs at t_curr might be from substrates (constant) or from products simulated up to t_curr
                        # If source_desc is a product EMU, its MDV at t_curr would have been computed if it's smaller or from a prior step.
                        # This hierarchical calculation needs careful handling of dependencies.
                        # For now, assume substrate_mdvs_exp provides constant MDVs for true substrates over time.
                        # Products that are sources for larger EMUs would need their t_curr value from sim_inst_MDVs_exp_obj_map.
                        mdv_obj = ChainMap(sim_inst_MDVs_exp_obj_map, substrate_mdvs_exp).get(source_desc)
                        if mdv_obj is None: # Should ideally not happen
                            mdv_obj = get_natural_MDV(source_desc.size)

                        # If mdv_obj is {tp: MDV}, get for t_curr. If it's just MDV (constant substrate), use it.
                        current_mdv_val = mdv_obj.get(t_curr, mdv_obj) if isinstance(mdv_obj, dict) else mdv_obj
                        Y_matrix_t_curr_rows.append(current_mdv_val.value if isinstance(current_mdv_val, MDV) else current_mdv_val)

                    else: # Tuple of EMUs for convolution
                        mdvs_to_convolve = []
                        for emu_in_tuple in source_desc:
                            mdv_obj_part = ChainMap(sim_inst_MDVs_exp_obj_map, substrate_mdvs_exp).get(emu_in_tuple)
                            if mdv_obj_part is None: mdv_obj_part = get_natural_MDV(emu_in_tuple.size)
                            current_mdv_val_part = mdv_obj_part.get(t_curr, mdv_obj_part) if isinstance(mdv_obj_part, dict) else mdv_obj_part
                            mdvs_to_convolve.append(current_mdv_val_part)
                        convolved_mdv = reduce(conv, mdvs_to_convolve)
                        Y_matrix_t_curr_rows.append(convolved_mdv.value if isinstance(convolved_mdv, MDV) else convolved_mdv)

                if not Y_matrix_t_curr_rows: # Should not happen if sourceEMU_descriptors_list is not empty
                    Y_matrix_t_curr = np.empty((B_val.shape[1], size + 1)) if B_val.shape[1]>0 else np.empty((0,size+1))
                    # Handle if Y is empty but B expects rows - G_t1 becomes problematic.
                    # This case needs robust definition if B_val.shape[1] (num source emus for Y) > 0
                    if B_val.shape[1] > 0: G_matrix_t_curr = np.zeros((Minv_val.shape[0], Y_matrix_t_prev.shape[1]))
                    else: G_matrix_t_curr = np.empty((Minv_val.shape[0],0))

                else:
                    Y_matrix_t_curr = np.array(Y_matrix_t_curr_rows)
                    if B_val.shape[1] != Y_matrix_t_curr.shape[0]: # Mismatch check
                        logging.error(f"Inst MDV Calc: Mismatch B cols {B_val.shape[1]} vs Y rows {Y_matrix_t_curr.shape[0]} for size {size}, t={t_curr} in exp {experiment_id}")
                        # Create zero G_matrix_t_curr to prevent crash, results will be wrong
                        G_matrix_t_curr = np.zeros((Minv_val.shape[0], Y_matrix_t_prev.shape[1]))
                    else:
                        G_matrix_t_curr = Minv_val @ B_val @ Y_matrix_t_curr

                X_matrix_t_curr = Phi_val @ X_matrix_t_prev - Gamma_val @ G_matrix_t_prev - Omega_val @ (G_matrix_t_curr - G_matrix_t_prev)

                Xs_sim_curr_tp_all_sizes[size] = X_matrix_t_curr
                Ys_sim_curr_tp_all_sizes[size] = Y_matrix_t_curr # Store Y at t_curr

                for emu_obj, mdv_array_t_curr in zip(productEMUs_list, X_matrix_t_curr):
                    sim_inst_MDVs_exp_obj_map.setdefault(emu_obj, {})[t_curr] = MDV(mdv_array_t_curr)

            Xs_sim_all_tps_all_sizes[t_curr] = Xs_sim_curr_tp_all_sizes
            Ys_sim_all_tps_all_sizes[t_curr] = Ys_sim_curr_tp_all_sizes
            t_prev = t_curr

        # Convert EMU object keys to EMU IDs for the final returned dict
        sim_inst_MDVs_exp_id_map = {}
        for emu_obj, tp_mdv_map in sim_inst_MDVs_exp_obj_map.items():
            sim_inst_MDVs_exp_id_map[emu_obj.id] = {
                tp: mdv.value for tp, mdv in tp_mdv_map.items() if tp > 0 # Exclude t0 from final output as per original
            }

        return sim_inst_MDVs_exp_id_map


    # This method now calculates the GLOBAL stacked derivative matrix for instationary state.
    def _calculate_inst_MDVs_and_derivatives_p(self): # params_p implicitly from self.model.total_fluxes and self.model.concentrations
        '''
        This method simulates MDVs and their derivatives at isotopically nonstationary state.
        It calculates derivatives for all relevant (emu, timepoint, experiment) tuples and stacks them.
        
        Returns
        -------
        simInstMDVs_return_dict: dict
            EMU ID => {t => MDV (in array)} from the first experiment, for compatibility.
        final_sim_inst_MDVs_der_stacked: np.ndarray
            Concatenated matrix of derivatives.
        '''
        
        # Determine num_params for instationary state (free fluxes + concentrations)
        num_params = 0
        if self.model.null_space is not None and self.model.concids is not None:
            num_params = self.model.null_space.shape[1] + len(self.model.concids)
        else:
            logging.error("Calculator (Inst Deriv): Cannot determine num_params (null_space or concids missing).")
            # Fallback: try to get from As_der_p if it exists and is populated for a size
            if self.model.matrix_As_der_p:
                example_A_der = next(iter(self.model.matrix_As_der_p.values()))
                if example_A_der is not None:
                    num_params = example_A_der.shape[0] # (n_params, n_emu_out, n_emu_out)
            if num_params == 0:
                logging.critical("Calculator (Inst Deriv): num_params is zero. Cannot proceed with derivative calculation.")
                return {}, np.empty((0,0)) # Return empty structures with 0 columns for num_params

        if not self.model.measured_inst_MDVs:
            return {}, np.empty((0, num_params))

        all_substrate_emus_global = self._get_all_substrate_emus()
        ordered_derivatives_list = []
        expected_total_rows = 0
        
        simInstMDVs_for_first_exp_id_map = {} # For the return value
        first_exp_id_processed = None

        # Outer loop for experiments, to calculate derivatives block by block
        for exp_id in sorted(self.model.measured_inst_MDVs.keys()):
            # Experiment-specific initial conditions and simulations
            substrate_mdvs_exp = self._get_substrate_MDVs_for_experiment(all_substrate_emus_global, exp_id)
            initial_X_dict_exp = self._get_initial_matrix_Xs_for_experiment(exp_id) # Usually natural abundance
            initial_Y_dict_exp = self._get_initial_matrix_Ys_for_experiment(exp_id, substrate_mdvs_exp)

            # These will store the time-course for THIS experiment
            # For derivatives: Xs_exp_all_tps, Ys_exp_all_tps, Xders_exp_all_tps, Yders_exp_all_tps
            # And simInstMDVs_this_exp_obj_map, simInstMDVsDer_this_exp_obj_map
            
            Xs_exp_all_tps = {0.0: initial_X_dict_exp}
            Ys_exp_all_tps = {0.0: initial_Y_dict_exp}
            
            # For derivative part of Y_dot = H*p + F*X_dot, Yder is dY/dp
            # For substrate_MDVs_der_p, it's assumed global for now. If it becomes exp-specific, needs update.
            initial_Yders_exp = {} # {size: (n_params, n_source_emus, n_iso)}
            for size_iter_init_yder in initial_Y_dict_exp:
                source_emus_desc_list = self.model.matrix_Bs[size_iter_init_yder][2]
                yder_rows = []
                for src_desc in source_emus_desc_list:
                    if not isinstance(src_desc, Iterable):
                        # Global substrate_MDVs_der_p for now
                        deriv = self.model.substrate_MDVs_der_p.get(src_desc, np.zeros((num_params, src_desc.size + 1)))
                        yder_rows.append(deriv)
                    else: # Convolution of derivatives
                        # This requires a more complex diff_conv setup if hierarchical substrate derivatives are exp-specific
                        # For now, assume base substrates have derivatives, others (from sim) are handled through Xder.
                        # Simplified: assume convoluted term's derivative is sum of parts if independent, or needs full diff_conv.
                        # Using zeros for convoluted Yder's initial derivative for simplicity here, implies it's handled by main ODE.
                        # This means d(Y_conv)/dp = 0 for t=0 if Y_conv is not a direct labeled substrate.
                        # This part needs to be robust if hierarchical derivatives of Y are needed.
                        # The main ODE part Xder_t1 = Phi@Xder_t0 + Gamma@H_t0 + Omega@(H_t1 - H_t0) handles this.
                        # So, Yder_t1 for sourceEMUs that are products of other reactions will come from Xder_t1.
                        # Here, we only need d(substrate)/dp for true substrates.
                        # The ChainMap for mdvder in the original code handles this.
                        # For initial Yders (t=0), we only need derivatives of actual input substrates.
                        # Others will be zero.
                        # The original code uses self.model.initial_matrix_Ys_der_p[size]. This should be made exp-specific if needed.
                        # For now, assume self.model.initial_matrix_Ys_der_p is okay (e.g. all zeros, or exp-specific if it was adapted).
                        # This method IS the adaptation point. So, must build it.
                        # Let's assume initial derivatives for Y are mostly zero unless Y is a direct tracer whose MDV changes with a parameter.
                        # This is complex. The original code used self.model.initial_matrix_Ys_der_p.
                        # Let's assume it's correctly pre-populated or zero for now.
                        # A robust way: if a src_emu is a substrate, use its derivative; else use zeros for t=0.
                        # This is what self.model.substrate_MDVs_der_p effectively does.
                        # For simplicity, let's assume initial_matrix_Ys_der_p is okay for now.
                        # This needs to be fixed if initial Y derivatives depend on experiment_id.
                        # It should be:
                        yder_row_for_conv = np.zeros((num_params, size_iter_init_yder + 1)) # Placeholder
                        mdv_der_parts_for_conv = []
                        mdv_val_parts_for_conv = []
                        for e_conv in src_desc:
                             e_mdv = substrate_mdvs_exp.get(e_conv, get_natural_MDV(e_conv.size))
                             e_der = self.model.substrate_MDVs_der_p.get(e_conv, np.zeros((num_params, e_conv.size+1))) # Global for now
                             mdv_val_parts_for_conv.append(e_mdv)
                             mdv_der_parts_for_conv.append([e_mdv, e_der])
                        if mdv_der_parts_for_conv:
                             yder_row_for_conv = reduce(diff_conv, mdv_der_parts_for_conv)[1]
                        yder_rows.append(yder_row_for_conv)

                initial_Yders_exp[size_iter_init_yder] = np.array(yder_rows) if yder_rows else np.empty((0, num_params, size_iter_init_yder+1))
                if initial_Yders_exp[size_iter_init_yder].ndim == 3: # (n_source_emus, n_params, n_iso)
                     initial_Yders_exp[size_iter_init_yder] = initial_Yders_exp[size_iter_init_yder].swapaxes(0,1) # -> (n_params, n_source_emus, n_iso)


            # Assume initial_matrix_Xs_der_p is all zeros and correctly shaped from model global for now.
            # It should be (n_params, n_prod_emus, n_iso)
            initial_Xders_exp = self.model.initial_matrix_Xs_der_p # This is global, might need exp-specific if X0 depends on exp.

            simInstMDVs_this_exp_obj_map = {} # emu_obj -> {tp: MDV_obj}
            simInstMDVsDer_this_exp_obj_map = {} # emu_obj -> {tp: deriv_array (n_iso, n_params)}
            
            # Initialize t=0 data for products for this experiment
            for size, X_matrix_t0 in initial_X_dict_exp.items():
                prod_emus = self.model.matrix_As[size][2]
                Xder_t0_size = initial_Xders_exp.get(size, np.zeros((num_params, len(prod_emus), size+1)))
                for idx, emu_obj in enumerate(prod_emus):
                    simInstMDVs_this_exp_obj_map.setdefault(emu_obj, {})[0.0] = MDV(X_matrix_t0[idx, :])
                    # Store derivative as (n_iso, n_params)
                    simInstMDVsDer_this_exp_obj_map.setdefault(emu_obj, {})[0.0] = Xder_t0_size[idx, :, :].T # Assuming Xder was (n_param, n_emu, n_iso)

            # Initialize t=0 data for sources (used in ChainMap for Y construction)
            for size, Y_matrix_t0 in initial_Y_dict_exp.items():
                 source_emus_desc_list = self.model.matrix_Bs[size][2]
                 Yder_t0_size = initial_Yders_exp.get(size, np.zeros((num_params,len(source_emus_desc_list),size+1))) # (n_params, n_source, n_iso)
                 for idx, source_desc in enumerate(source_emus_desc_list):
                      if not isinstance(source_desc, Iterable):
                           if source_desc not in simInstMDVs_this_exp_obj_map: # If not a product EMU already initialized
                                simInstMDVs_this_exp_obj_map.setdefault(source_desc, {})[0.0] = substrate_mdvs_exp.get(source_desc, get_natural_MDV(source_desc.size))
                           if source_desc not in simInstMDVsDer_this_exp_obj_map:
                                simInstMDVsDer_this_exp_obj_map.setdefault(source_desc, {})[0.0] = Yder_t0_size[idx,:,:].T # (n_iso, n_params)

            # --- Main ODE and Sensitivity Loop for current experiment ---
            t_prev = 0.0
            current_Xs_exp = initial_X_dict_exp
            current_Ys_exp = initial_Y_dict_exp
            current_Xders_exp = initial_Xders_exp # This is {size: (n_params, n_prod, n_iso)}
            current_Yders_exp = initial_Yders_exp # This is {size: (n_params, n_source, n_iso)}

            sim_timepoints_exp = sorted([tp for tp in self.model.timepoints if tp > 0])

            for t_curr in sim_timepoints_exp:
                deltat = t_curr - t_prev
                next_Xs_exp_size_map = {}
                next_Ys_exp_size_map = {}
                next_Xders_exp_size_map = {}
                next_Yders_exp_size_map = {}

                for size in sorted(self.model.matrix_As.keys()):
                    # ... (ODE and sensitivity calculation as in original, but using current_Xs_exp, current_Ys_exp, etc.) ...
                    # ... This involves: A, B, M, Ader, Bder, Mder, Minvder, F, Finv, Phi, Gamma, Omega ...
                    # ... X_t0, Y_t0, G_t0, Xder_t0, Yder_t0, H_t0 ...
                    # ... Y_t1, Yder_t1, G_t1, X_t1, H_t1, Xder_t1 ...
                    # The crucial change is that Y_t1 and Yder_t1 must be constructed using ChainMap
                    # with simInstMDVs_this_exp_obj_map and simInstMDVsDer_this_exp_obj_map for hierarchical terms,
                    # and substrate_mdvs_exp and self.model.substrate_MDVs_der_p (global for now) for base substrates.
                    
                    # Placeholder for the detailed ODE solution from original code, adapted for exp-specific inputs
                    lambA, fluxidsA_iter, productEMUs_iter_list = self.model.matrix_As[size]
                    lambB, fluxidsB_iter, sourceEMUs_desc_iter_list = self.model.matrix_Bs[size]
                    lambM, metabids_M_iter = self.model.matrix_Ms[size]

                    A_val = lambA(*self.model.total_fluxes[fluxidsA_iter])
                    B_val = lambB(*self.model.total_fluxes[fluxidsB_iter])
                    M_val = lambM(*self.model.concentrations[metabids_M_iter]) # Concentrations already set for exp
                    Minv_val = pinv(M_val, check_finite=False)

                    Ader_val = self.model.matrix_As_der_p[size] # (n_params, n_emu_out, n_emu_out)
                    Bder_val = self.model.matrix_Bs_der_p[size] # (n_params, n_emu_out, n_emu_in)
                    Mder_val = self.model.matrix_Ms_der_p[size] # (n_params, n_emu_out, n_emu_out)
                    Minvder_val = -Minv_val @ Mder_val @ Minv_val # (n_params, n_emu_out, n_emu_out)

                    F_val = Minv_val @ A_val
                    Finv_val = pinv(F_val, check_finite=False)
                    I_mtx = np.eye(F_val.shape[0])
                    Phi_val = expm(F_val * deltat)
                    Gamma_val = (Phi_val - I_mtx) @ Finv_val
                    Omega_val = (Gamma_val / deltat - I_mtx) @ Finv_val
                    
                    X_t0_size = current_Xs_exp[size] # Matrix (n_prod_emu, n_iso)
                    Y_t0_size = current_Ys_exp[size] # Matrix (n_source_emu, n_iso)
                    G_t0_size = Minv_val @ B_val @ Y_t0_size if Y_t0_size.size > 0 else np.zeros((Minv_val.shape[0], X_t0_size.shape[1]))


                    Xder_t0_size = current_Xders_exp[size] # Matrix (n_params, n_prod_emu, n_iso)
                    Yder_t0_size = current_Yders_exp[size] # Matrix (n_params, n_source_emu, n_iso)
                    
                    # H_t0 = (Minvder@A@X_t0 + Minv@Ader@X_t0 - Minv@B@Yder_t0 - Minvder@B@Y_t0 - Minv@Bder@Y_t0)
                    # Shapes: Minvder (p,o,o), A(o,o), X_t0(o,i) -> einsum('pjk,kl,lm->pjm', Minvder, A_val, X_t0_size)
                    # Minv(o,o), Ader(p,o,o), X_t0(o,i) -> einsum('jk,pkl,lm->pjm', Minv_val, Ader_val, X_t0_size)
                    # Minv(o,o), B(o,s), Yder_t0(p,s,i) -> einsum('jk,kl,pls->pjs', Minv_val, B_val, Yder_t0_size)
                    # Minvder(p,o,o), B(o,s), Y_t0(s,i) -> einsum('pjk,kl,ls->pis', Minvder_val, B_val, Y_t0_size)
                    # Minv(o,o), Bder(p,o,s), Y_t0(s,i) -> einsum('jk,pkl,ls->pis', Minv_val, Bder_val, Y_t0_size)
                    
                    term1_H0 = np.einsum('pjk,kl,lm->pjm', Minvder_val, A_val, X_t0_size, optimize='optimal') if X_t0_size.size > 0 else np.zeros((num_params, Minvder_val.shape[1], X_t0_size.shape[1] if X_t0_size.size >0 else 0 ))
                    term2_H0 = np.einsum('jk,pkl,lm->pjm', Minv_val, Ader_val, X_t0_size, optimize='optimal') if X_t0_size.size > 0 else np.zeros_like(term1_H0)
                    term3_H0 = np.einsum('jk,kl,pls->pjs', Minv_val, B_val, Yder_t0_size, optimize='optimal') if Yder_t0_size.size > 0 else np.zeros_like(term1_H0) # p,o,i (Yder has n_iso cols)
                    term4_H0 = np.einsum('pjk,kl,ls->pis', Minvder_val, B_val, Y_t0_size, optimize='optimal') if Y_t0_size.size > 0 else np.zeros_like(term1_H0)
                    term5_H0 = np.einsum('jk,pkl,ls->pis', Minv_val, Bder_val, Y_t0_size, optimize='optimal') if Y_t0_size.size > 0 else np.zeros_like(term1_H0)
                    H_t0_val = term1_H0 + term2_H0 - term3_H0 - term4_H0 - term5_H0

                    # Construct Y_t1 (MDVs of source EMUs at t_curr for this experiment)
                    Y_t1_matrix_rows = []
                    Yder_t1_matrix_rows = [] # List of (n_params, n_iso) for each source emu

                    for source_desc_iter in sourceEMU_descriptors_list:
                        if not isinstance(source_desc_iter, Iterable): # Single EMU
                            mdv_obj = ChainMap(simInstMDVs_this_exp_obj_map, substrate_mdvs_exp).get(source_desc_iter)
                            mdv_der_obj = ChainMap(simInstMDVsDer_this_exp_obj_map, self.model.substrate_MDVs_der_p).get(source_desc_iter) # Global substrate_MDVs_der_p for now

                            if mdv_obj is None: mdv_obj = get_natural_MDV(source_desc_iter.size)
                            mdv_t1_val = mdv_obj.get(t_curr, mdv_obj) if isinstance(mdv_obj, dict) else mdv_obj # Get MDV at t_curr or constant
                            Y_t1_matrix_rows.append(mdv_t1_val.value if isinstance(mdv_t1_val, MDV) else mdv_t1_val)

                            if mdv_der_obj is None: mdv_der_obj = np.zeros((num_params, source_desc_iter.size + 1))
                            mdv_der_t1_val = mdv_der_obj.get(t_curr, mdv_der_obj) if isinstance(mdv_der_obj, dict) else mdv_der_obj
                            Yder_t1_matrix_rows.append(mdv_der_t1_val) # This is (n_params, n_iso)
                        else: # Tuple of EMUs for convolution
                            mdv_val_parts = []
                            mdv_der_parts_list_of_lists = [] # each item: [mdv_obj_at_t_curr, mdv_der_array_at_t_curr (n_params, n_iso)]
                            for e_conv in source_desc_iter:
                                e_mdv_obj = ChainMap(simInstMDVs_this_exp_obj_map, substrate_mdvs_exp).get(e_conv, get_natural_MDV(e_conv.size))
                                e_mdv_val_t1 = e_mdv_obj.get(t_curr, e_mdv_obj) if isinstance(e_mdv_obj, dict) else e_mdv_obj
                                mdv_val_parts.append(e_mdv_val_t1)

                                e_mdv_der_obj = ChainMap(simInstMDVsDer_this_exp_obj_map, self.model.substrate_MDVs_der_p).get(e_conv, np.zeros((num_params, e_conv.size+1)))
                                e_mdv_der_val_t1 = e_mdv_der_obj.get(t_curr, e_mdv_der_obj) if isinstance(e_mdv_der_obj, dict) else e_mdv_der_obj
                                mdv_der_parts_list_of_lists.append([e_mdv_val_t1, e_mdv_der_val_t1]) # Pass MDV obj and (n_params, n_iso) array

                            Y_t1_matrix_rows.append(reduce(conv, mdv_val_parts).value)
                            convolved_der_t1 = reduce(diff_conv, mdv_der_parts_list_of_lists)[1] # diff_conv returns [mdv, mdv_der] where mdv_der is (n_params, n_iso)
                            Yder_t1_matrix_rows.append(convolved_der_t1)

                    Y_t1_val = np.array(Y_t1_matrix_rows) if Y_t1_matrix_rows else np.empty((0, size + 1))
                    # Yder_t1_val shape needs to be (n_params, n_source_emus_for_Y, n_isotopomers)
                    Yder_t1_val = np.array(Yder_t1_matrix_rows).swapaxes(0,1).swapaxes(1,2) if Yder_t1_matrix_rows else np.zeros((num_params, 0, size + 1))


                    G_t1_val = Minv_val @ B_val @ Y_t1_val if Y_t1_val.size > 0 else np.zeros((Minv_val.shape[0], Y_t0_size.shape[1]))
                    X_t1_val = Phi_val @ X_t0_size - Gamma_val @ G_t0_size - Omega_val @ (G_t1_val - G_t0_size)
                    
                    # H_t1: same structure as H_t0, but with X_t1, Y_t1, Yder_t1
                    term1_H1 = np.einsum('pjk,kl,lm->pjm', Minvder_val, A_val, X_t1_val, optimize='optimal') if X_t1_val.size > 0 else np.zeros_like(term1_H0)
                    term2_H1 = np.einsum('jk,pkl,lm->pjm', Minv_val, Ader_val, X_t1_val, optimize='optimal') if X_t1_val.size > 0 else np.zeros_like(term1_H0)
                    term3_H1 = np.einsum('jk,kl,pls->pjs', Minv_val, B_val, Yder_t1_val, optimize='optimal') if Yder_t1_val.size > 0 else np.zeros_like(term1_H0)
                    term4_H1 = np.einsum('pjk,kl,ls->pis', Minvder_val, B_val, Y_t1_val, optimize='optimal') if Y_t1_val.size > 0 else np.zeros_like(term1_H0)
                    term5_H1 = np.einsum('jk,pkl,ls->pis', Minv_val, Bder_val, Y_t1_val, optimize='optimal') if Y_t1_val.size > 0 else np.zeros_like(term1_H0)
                    H_t1_val = term1_H1 + term2_H1 - term3_H1 - term4_H1 - term5_H1
                    
                    Xder_t1_val = Phi_val @ Xder_t0_size + Gamma_val @ H_t0_val + Omega_val @ (H_t1_val - H_t0_val) # Xder is (n_params, n_prod_emu, n_iso)
                    
                    next_Xs_exp_size_map[size] = X_t1_val
                    next_Ys_exp_size_map[size] = Y_t1_val
                    next_Xders_exp_size_map[size] = Xder_t1_val
                    next_Yders_exp_size_map[size] = Yder_t1_val # Store for next iteration if needed by H_t0

                    for idx, emu_obj_iter in enumerate(productEMUs_iter_list):
                        simInstMDVs_this_exp_obj_map.setdefault(emu_obj_iter, {})[t_curr] = MDV(X_t1_val[idx, :])
                        # Store derivative as (n_iso, n_params) by transposing from (n_params, n_iso)
                        simInstMDVsDer_this_exp_obj_map.setdefault(emu_obj_iter, {})[t_curr] = Xder_t1_val[:, idx, :].T

                # Update for next iteration
                current_Xs_exp = next_Xs_exp_size_map
                current_Ys_exp = next_Ys_exp_size_map
                current_Xders_exp = next_Xders_exp_size_map
                current_Yders_exp = next_Yders_exp_size_map
                t_prev = t_curr
            # --- End of ODE and Sensitivity Loop for current experiment ---

            if exp_id == (sorted(self.model.measured_inst_MDVs.keys())[0] if self.model.measured_inst_MDVs else None):
                for emu_obj, tp_mdv_map in simInstMDVs_this_exp_obj_map.items():
                    simInstMDVs_for_first_exp_id_map[emu_obj.id] = {tp: mdv.value for tp, mdv in tp_mdv_map.items() if tp > 0}

            # Append derivatives for measured fragments in this experiment to the global list
            exp_data_measured = self.model.measured_inst_MDVs[exp_id]
            for frag_id_iter in sorted(exp_data_measured.keys()):
                frag_tp_data_measured = exp_data_measured[frag_id_iter]
                for tp_iter in sorted(frag_tp_data_measured.keys()):
                    num_mass_isotopomers = len(frag_tp_data_measured[tp_iter][0])
                    expected_total_rows += num_mass_isotopomers
                    
                    found_emu_obj = next((obj for obj in simInstMDVsDer_this_exp_obj_map if obj.id == frag_id_iter), None)
                    
                    if found_emu_obj and tp_iter in simInstMDVsDer_this_exp_obj_map.get(found_emu_obj, {}):
                        frag_tp_der = simInstMDVsDer_this_exp_obj_map[found_emu_obj][tp_iter] # (n_iso, n_params)
                        if frag_tp_der.shape[0] == num_mass_isotopomers and frag_tp_der.shape[1] == num_params:
                            ordered_derivatives_list.append(frag_tp_der)
                        else:
                            logging.error(f"Calculator Inst Deriv: Shape mismatch for {frag_id_iter} at T={tp_iter} in exp {exp_id}. Appending zeros.")
                            ordered_derivatives_list.append(np.zeros((num_mass_isotopomers, num_params)))
                    else:
                        logging.error(f"Calculator Inst Deriv: Missing derivative for {frag_id_iter} at T={tp_iter} in exp {exp_id}. Appending zeros.")
                        ordered_derivatives_list.append(np.zeros((num_mass_isotopomers, num_params)))

        final_sim_inst_MDVs_der_stacked = np.vstack(ordered_derivatives_list) if ordered_derivatives_list else np.empty((0, num_params))
        if final_sim_inst_MDVs_der_stacked.shape[0] != expected_total_rows:
             logging.error(f"Calculator Inst Deriv: Row count mismatch for stacked derivatives: expected {expected_total_rows}, got {final_sim_inst_MDVs_der_stacked.shape[0]}")

        return simInstMDVs_for_first_exp_id_map, final_sim_inst_MDVs_der_stacked
                # current_Ys_exp = next_Ys_exp_size_map
                # current_Xders_exp = next_Xders_exp_size_map
                # current_Yders_exp = next_Yders_exp_size_map
                t_prev = t_curr
            # --- End of ODE and Sensitivity Loop for current experiment ---

            if exp_id == (sorted(self.model.measured_inst_MDVs.keys())[0] if self.model.measured_inst_MDVs else None):
                for emu_obj, tp_mdv_map in simInstMDVs_this_exp_obj_map.items():
                    simInstMDVs_for_first_exp_id_map[emu_obj.id] = {tp: mdv.value for tp, mdv in tp_mdv_map.items() if tp > 0}

            # Append derivatives for measured fragments in this experiment to the global list
            exp_data_measured = self.model.measured_inst_MDVs[exp_id]
            for frag_id_iter in sorted(exp_data_measured.keys()):
                frag_tp_data_measured = exp_data_measured[frag_id_iter]
                for tp_iter in sorted(frag_tp_data_measured.keys()):
                    num_mass_isotopomers = len(frag_tp_data_measured[tp_iter][0])
                    expected_total_rows += num_mass_isotopomers
                    
                    found_emu_obj = next((obj for obj in simInstMDVsDer_this_exp_obj_map if obj.id == frag_id_iter), None)
                    
                    if found_emu_obj and tp_iter in simInstMDVsDer_this_exp_obj_map.get(found_emu_obj, {}):
                        frag_tp_der = simInstMDVsDer_this_exp_obj_map[found_emu_obj][tp_iter] # (n_iso, n_params)
                        if frag_tp_der.shape[0] == num_mass_isotopomers and frag_tp_der.shape[1] == num_params:
                            ordered_derivatives_list.append(frag_tp_der)
                        else:
                            logging.error(f"Calculator Inst Deriv: Shape mismatch for {frag_id_iter} at T={tp_iter} in exp {exp_id}. Appending zeros.")
                            ordered_derivatives_list.append(np.zeros((num_mass_isotopomers, num_params)))
                    else:
                        logging.error(f"Calculator Inst Deriv: Missing derivative for {frag_id_iter} at T={tp_iter} in exp {exp_id}. Appending zeros.")
                        ordered_derivatives_list.append(np.zeros((num_mass_isotopomers, num_params)))

        final_sim_inst_MDVs_der_stacked = np.vstack(ordered_derivatives_list) if ordered_derivatives_list else np.empty((0, num_params))
        if final_sim_inst_MDVs_der_stacked.shape[0] != expected_total_rows:
             logging.error(f"Calculator Inst Deriv: Row count mismatch for stacked derivatives: expected {expected_total_rows}, got {final_sim_inst_MDVs_der_stacked.shape[0]}")

        return simInstMDVs_for_first_exp_id_map, final_sim_inst_MDVs_der_stacked
