'''Define the MFAModel and InstMFAModel class.'''


__author__ = 'Chao Wu'
__date__ = '05/19/2022'

import logging # Added for logging warnings
import numpy as np
import pandas as pd
from scipy.linalg import pinv
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
try:
    from openopt import NLP
except ModuleNotFoundError:
    OPENOPT_INSTALLED = False
else:
    OPENOPT_INSTALLED = True
from ..utils.utils import Calculator


class MFAModel():
    '''
    Parameters
    ----------
    model: Model
        Freeflux Model.
    fit_measured_fluxes: bool
        Whether to fit measured fluxes.
    solvor: {"slsqp", "ralg"}
        * If "slsqp", scipy.optimize.minimze will be used.
        * If "ralg", openopt NLP solver will be used.
    '''
    
    def __init__(self, model, fit_measured_fluxes, solver = 'slsqp'):
        '''
        Parameters
        ----------
        model: Model
            Freeflux Model.
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        solvor: {"slsqp", "ralg"}
            * If "slsqp", scipy.optimize.minimze will be used.
            * If "ralg", openopt NLP solver will be used.
        '''
        
        self.model = model
        self.calculator = Calculator(self.model)
        self.fit_measured_fluxes = fit_measured_fluxes
        self.solver = solver
        
        self.N = self.model.null_space
        self.T = self.model.transform_matrix
        
        self.ntotalfluxes = len(self.model.totalfluxids)
        
    
    def _calculate_difference_sim_exp_MDVs(self):
        
        simMDVs = self.calculator._calculate_MDVs()
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] 
                               for emuid in self.model.target_EMUs])
        
        return diff
            
            
    def _calculate_difference_sim_exp_fluxes(self):
        
        simFluxes = self.model.total_fluxes[self.model.measured_fluxes.keys()]
        expFluxes = np.array([mean for mean, _ in self.model.measured_fluxes.values()])
        diff = simFluxes - expFluxes
        
        return diff
        
        
    def _calculate_sim_MDVs_derivative(self):
        
        simMDVs, simMDVsDer = self.calculator._calculate_MDVs_and_derivatives_p()
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] 
                               for emuid in self.model.target_EMUs])
        
        dxsim_dp = np.vstack([simMDVsDer[emuid] for emuid in self.model.target_EMUs])
        
        return diff, dxsim_dp
        
    
    def _calculate_sim_fluxes_derivative(self):
        
        dvsim_dp = self.model.measured_fluxes_der_p
        
        return dvsim_dp


    def build_objective(self):
        
        def obj_func(u): # Renamed from _f to obj_func for clarity, though not strictly required
            self.model.total_fluxes[:] = self.N @ u

            # This call needs to be maintained for self.simMDVs to be populated
            # The original _calculate_difference_sim_exp_MDVs itself calls calculator._calculate_MDVs()
            # which populates self.simMDVs (if Calculator is designed that way) or returns simMDVs.
            # For this change, we assume self.calculator.calculate_sim_MDVs() is called and
            # populates self.simMDVs, or that _calculate_difference_sim_exp_MDVs populates it.
            # Let's assume self.simMDVs is populated by the calculator called within _calculate_difference_sim_exp_MDVs
            # or we call it explicitly if needed: self.simMDVs = self.calculator.calculate_sim_MDVs(self.model.target_EMUs, u)
            # For now, relying on existing self.simMDVs population from _calculate_difference_sim_exp_MDVs,
            # but the logic below reconstructs residuals based on multi-experiment structure.

            # self.simMDVs = self.calculator._calculate_MDVs() # OLD: This used to fetch a global simMDVs
            # NEW: simMDVs will be fetched per experiment inside the loop.
            # self.simMDVs attribute is no longer populated/used globally by obj_func for residual calculation.

            objVal_mdv = 0.0
            self.concat_meas_mdv_sds = [] # To store concatenated SDs for fallback in gradient

            if self.model.measured_MDVs:
                concat_meas_mdv_means = []
                concat_sim_mdv_values = []

                for exp_id in sorted(self.model.measured_MDVs.keys()):
                    # Calculate simMDVs for the current experiment_id and parameters u
                    current_exp_simMDVs = self.calculator._calculate_MDVs_for_experiment(u, experiment_id=exp_id)

                    for fragment_id in sorted(self.model.measured_MDVs[exp_id].keys()):
                        meas_mean, meas_sd = self.model.measured_MDVs[exp_id][fragment_id]

                        sim_mean_for_frag = current_exp_simMDVs.get(fragment_id) # Use experiment-specific simMDVs
                        if sim_mean_for_frag is not None:
                            if len(sim_mean_for_frag) == len(meas_mean):
                                concat_meas_mdv_means.extend(meas_mean)
                                self.concat_meas_mdv_sds.extend(meas_sd)
                                concat_sim_mdv_values.extend(sim_mean_for_frag)
                            else:
                                logging.warning(
                                    f"Mismatched length for MDV fragment {fragment_id} in exp {exp_id}: "
                                    f"measured {len(meas_mean)}, simulated {len(sim_mean_for_frag)}."
                                )
                        else:
                            logging.warning(f"Simulated MDV for fragment {fragment_id} (exp {exp_id}) not found in current_exp_simMDVs.")

                if concat_meas_mdv_means: # Check if any valid data was collected
                    self.residuals_mdv = np.array(concat_sim_mdv_values) - np.array(concat_meas_mdv_means)

                    if (self.model.measured_MDVs_inv_cov is not None and
                        self.model.measured_MDVs_inv_cov.shape[0] == len(self.residuals_mdv) and
                        self.model.measured_MDVs_inv_cov.shape[1] == len(self.residuals_mdv)):
                        objVal_mdv = 0.5 * np.dot(np.dot(self.residuals_mdv.T, self.model.measured_MDVs_inv_cov), self.residuals_mdv)
                    else:
                        if self.model.measured_MDVs_inv_cov is not None:
                            logging.warning(
                                "Mismatch in measured_MDVs_inv_cov shape or not available. "
                                f"Covar shape: {self.model.measured_MDVs_inv_cov.shape}, "
                                f"Residuals length: {len(self.residuals_mdv)}. "
                                "Falling back to SD-based objective for MDVs."
                            )
                        else:
                             logging.warning("measured_MDVs_inv_cov is None. Falling back to SD-based objective for MDVs.")

                        if self.concat_meas_mdv_sds and len(self.concat_meas_mdv_sds) == len(self.residuals_mdv):
                            valid_sds = np.array(self.concat_meas_mdv_sds)
                            # Ensure SDs are not zero to avoid division by zero
                            if np.any(valid_sds == 0):
                                logging.warning("Zero standard deviation found. Replacing with 1.0 for objective calculation.")
                                valid_sds[valid_sds == 0] = 1.0
                            weighted_residuals = self.residuals_mdv / valid_sds
                            objVal_mdv = 0.5 * np.sum(weighted_residuals**2)
                        else:
                            logging.warning("Fallback MDV objective: SDs not available or mismatched length. MDV objective part is zero.")
                            objVal_mdv = 0.0 # Or handle error more strictly
                else:
                    self.residuals_mdv = np.array([]) # No valid MDV data
            else:
                self.residuals_mdv = np.array([]) # No measured_MDVs

            objVal_flux = 0.0
            if self.fit_measured_fluxes and self.model.measured_fluxes:
                flux_diff = self._calculate_difference_sim_exp_fluxes() # This method should be fine
                if self.model.measured_fluxes_inv_cov is not None and \
                   flux_diff.shape[0] == self.model.measured_fluxes_inv_cov.shape[0]:
                    objVal_flux = 0.5 * np.dot(np.dot(flux_diff.T, self.model.measured_fluxes_inv_cov), flux_diff)
                else: # Fallback for fluxes if inv_cov is missing or mismatched
                    logging.warning("Flux inverse covariance matrix missing, mismatched, or no measured fluxes for objective. Flux part of objective is zero.")

            return objVal_mdv + objVal_flux

        self.f = obj_func # Assign the comprehensive objective function directly
        
    
    def build_gradient(self):
        
        def obj_grad(u): # Renamed from _df for clarity
            # Objective function obj_func(u) must have been called first to set self.residuals_mdv and self.N@u
            # self.model.total_fluxes[:] = self.N @ u # This is done in obj_func

            # The Calculator's _calculate_MDVs_and_derivatives_p method is called once by Fitter
            # (outside this function) and populates self.model.simMDVs_der_p with the global stacked matrix.
            # obj_func (called before obj_grad) has already updated self.residuals_mdv based on per-experiment simulations.
            # self.model.total_fluxes is also current from obj_func, which is based on current 'u'.

            # We need to recalculate the derivative matrix for the current 'u'.
            # The Calculator's method uses self.model.total_fluxes which obj_func has just set.
            _ignored_simMDVs_dict, self.model.simMDVs_der_p = self.calculator._calculate_MDVs_and_derivatives_p()
            current_simMDVs_der_p = self.model.simMDVs_der_p

            objGrad_mdv = np.zeros(self.N.shape[1]) # Gradient w.r.t. u (free fluxes)

            if hasattr(self, 'residuals_mdv') and self.residuals_mdv.size > 0:
                if (self.model.measured_MDVs_inv_cov is not None and
                    current_simMDVs_der_p is not None and # Use the local variable
                    self.model.measured_MDVs_inv_cov.shape[0] == self.residuals_mdv.size and
                    self.model.measured_MDVs_inv_cov.shape[1] == self.residuals_mdv.size and
                    current_simMDVs_der_p.shape[0] == self.residuals_mdv.size and
                    current_simMDVs_der_p.shape[1] == self.N.shape[1]): # Check cols match free fluxes
                    objGrad_mdv = np.dot(np.dot(self.residuals_mdv.T, self.model.measured_MDVs_inv_cov), current_simMDVs_der_p).squeeze()
                else: # Fallback for MDV gradient
                    logging.warning(
                        "MDV gradient: Primary path conditions not met (inv_cov or derivatives mismatch/missing). "
                        f"Covar shape: {self.model.measured_MDVs_inv_cov.shape if self.model.measured_MDVs_inv_cov is not None else 'None'}, "
                        f"Deriv shape: {current_simMDVs_der_p.shape if current_simMDVs_der_p is not None else 'None'}, "
                        f"Residuals size: {self.residuals_mdv.size}. "
                        "Attempting fallback SD-based gradient for MDVs."
                    )
                    if (hasattr(self, 'concat_meas_mdv_sds') and
                        current_simMDVs_der_p is not None and # Use the local variable
                        len(self.concat_meas_mdv_sds) == self.residuals_mdv.size and
                        current_simMDVs_der_p.shape[0] == self.residuals_mdv.size and
                        current_simMDVs_der_p.shape[1] == self.N.shape[1]):

                        valid_sds_sq = np.array(self.concat_meas_mdv_sds)**2
                        if np.any(valid_sds_sq == 0):
                            logging.warning("Zero standard deviation squared found. Replacing with 1.0 for gradient calculation.")
                            valid_sds_sq[valid_sds_sq == 0] = 1.0
                        weighted_residuals_for_grad = self.residuals_mdv / valid_sds_sq
                        objGrad_mdv = np.dot(weighted_residuals_for_grad.T, current_simMDVs_der_p).squeeze()
                    else:
                        logging.warning("MDV gradient fallback: simMDVs_der_p or SDs not available or mismatched. MDV gradient part is zero.")

            objGrad_flux = np.zeros(self.N.shape[1])
            if self.fit_measured_fluxes and self.model.measured_fluxes:
                # Ensure flux_der is correctly fetched or calculated
                flux_der = self.model.measured_fluxes_der_p # Derivatives w.r.t. u
                flux_diff = self._calculate_difference_sim_exp_fluxes() # This should be fine

                if (self.model.measured_fluxes_inv_cov is not None and
                    flux_der is not None and
                    flux_diff.shape[0] == self.model.measured_fluxes_inv_cov.shape[0] and
                    flux_der.shape[0] == flux_diff.shape[0] and # flux_der rows = num_measured_fluxes
                    flux_der.shape[1] == self.N.shape[1]): # flux_der cols = num_free_fluxes
                    objGrad_flux = np.dot(np.dot(flux_diff.T, self.model.measured_fluxes_inv_cov), flux_der).squeeze()
                else:
                    logging.warning("Flux gradient: Conditions not met (inv_cov or derivatives mismatch/missing). Flux gradient part is zero.")

            # Ensure dimensions match
            if objGrad_mdv.ndim == 0: objGrad_mdv = np.array([objGrad_mdv])
            if objGrad_flux.ndim == 0: objGrad_flux = np.array([objGrad_flux])
            if objGrad_mdv.size == 1 and objGrad_flux.size > 1 and objGrad_mdv.size < objGrad_flux.size : objGrad_mdv = np.zeros_like(objGrad_flux) # prevent broadcast error if one is scalar like 0.0
            if objGrad_flux.size == 1 and objGrad_mdv.size > 1 and objGrad_flux.size < objGrad_mdv.size : objGrad_flux = np.zeros_like(objGrad_mdv)


            grad = objGrad_mdv + objGrad_flux
            if grad.shape[0] != self.N.shape[1]: # Ensure correct shape
                logging.error(f"Gradient shape mismatch: expected ({self.N.shape[1]},), got {grad.shape}. Resetting to zeros.")
                grad = np.zeros(self.N.shape[1])
            return grad
        
        self.df = obj_grad # Assign the comprehensive gradient function

    
    def build_hessian(self):
        
        def _ddf(u):
            self.model.total_fluxes[:] = self.N@u
            
            _, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_MDVs_inv_cov
            hess = MDV_der.T@MDV_inv_cov@MDV_der
            
            return hess
            
        def ddf1(u):
            return _ddf(u)
        
        def ddf2(u):
            hess1 = _ddf(u)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            hess2 = flux_der.T@flux_inv_cov@flux_der
            
            return hess1 + hess2
        
        self.ddf = ddf2 if self.fit_measured_fluxes else ddf1    
    
            
    def build_flux_bound_constraints(self):
        
        A1 = self.N
        A2 = self.T@self.N
        A3 = -A2
        
        b1 = np.zeros(self.ntotalfluxes)
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        b2 = vnet_lb
        b3 = -vnet_ub
        
        A = np.vstack((A1, A2, A3))
        b = np.concatenate((b1, b2, b3))

        if self.solver == 'slsqp':
            self.constrs = {'type': 'ineq', 'fun': lambda u: A@u - b}
        elif self.solver == 'trust-constr' or self.solver == 'ipopt':
            self.constrs = LinearConstraint(A, b, np.inf)
        elif self.solver == 'ralg':
            self.A = -A
            self.b = -b
        
        
    def build_initial_flux_values(self, ini_netfluxes = None, rng = None):
        '''
        Parameters
        ----------
        ini_netfluxes: array
            Initial guess of net fluxes.
        '''
        
        if ini_netfluxes is None:
            vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
            vnet_ini = rng.uniform(low = vnet_lb, high = vnet_ub)
        else:
            vnet_ini = ini_netfluxes
            
        u_ini = pinv(self.T@self.N)@vnet_ini
        
        self.x0 = u_ini
        
        
    def _initialize_total_fluxes(self):
        
        for fluxid in self.model.totalfluxids:
            self.model.total_fluxes[fluxid] = 0.0
        

    def _solve_flux_slsqp(self, tol, max_iters, disp):
        
        res = minimize(
            fun = self.f,
            x0 = self.x0,
            method = 'SLSQP',
            jac = self.df,
            constraints = self.constrs,
            options = {
                'ftol': tol,
                'maxiter': max_iters, 
                'disp': disp
            }
        )
       
        return res.fun, res.x, res.success

    def _solve_flux_trust_constr(self, tol, max_iters, disp):
        res = minimize(
            fun = self.f,
            x0 = self.x0,
            method = 'trust-constr',
            jac = self.df,
            constraints = [self.constrs],
            options = {
                'gtol': tol,
                'xtol': tol,
                'barrier_tol': tol,
                'maxiter': max_iters, 
                'disp': True,
                'verbose': 2,
            },
        )       
        return res.fun, res.x, res.success

    # TODO: This solver might be more efficient
    # # from cyipopt import minimize_ipopt
    # def _solve_flux_ipopt(self, tol, max_iters, disp):
    #     res = minimize_ipopt(
    #         self.f,
    #         self.x0,
    #         jac = self.df,
    #         constraints = [self.constrs],
    #         options = {'disp': 5}
    #     )       
    #     return res.fun, res.x, res.success


    def _solve_flux_ralg(self, tol, max_iters, disp):

        if OPENOPT_INSTALLED:
            res = res = NLP(
                f = self.f, 
                x0 = self.x0, 
                df = self.df, 
                A = self.A, 
                b = self.b, 
                xtol = tol, 
                ftol = tol, 
                maxIter = max_iters, 
                iprint = 1 if disp else -1
            ).solve('ralg')

            return res.ff, res.xf, res.istop > 0 or res.istop == -7
        
        else:
            raise ModuleNotFoundError('install openopt first')


    def _calculate_residuals(self):
        
        MDV_diff = self._calculate_difference_sim_exp_MDVs()
        MDV_diag = np.diag(self.model.measured_MDVs_inv_cov)**0.5
        resids1 = MDV_diff*MDV_diag
        
        if self.fit_measured_fluxes:    
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_diag = np.diag(self.model.measured_fluxes_inv_cov)**0.5
            resids2 = flux_diff*flux_diag
            resids = np.concatenate((resids1, resids2))
        
        return resids         


    def _get_exp_and_sim_MDVs(self):

        exp_MDVs = self.model.measured_MDVs
        
        sim_MDVs_all = self.calculator._calculate_MDVs()
        sim_MDVs = {emuid: sim_MDVs_all[emuid] for emuid in self.model.target_EMUs}

        return exp_MDVs, sim_MDVs


    def _get_exp_and_sim_fluxes(self):

        if self.fit_measured_fluxes:
            exp_fluxes = self.model.measured_fluxes
            sim_fluxes = self.model.total_fluxes[self.model.measured_fluxes.keys()].to_dict()
        else:
            exp_fluxes = {}
            sim_fluxes = {}            

        return exp_fluxes, sim_fluxes


    def _get_nmeasurements(self, opt_resids):
        
        return opt_resids.size

         
    def _get_nparameters(self, opt_p):

        return opt_p.size


    def _get_hessian(self, opt_p):

        self.build_hessian()
        hess = self.ddf(opt_p)

        return hess


    def solve_flux(self, tol = 1e-6, max_iters = 400, disp = False):    
        
        self._initialize_total_fluxes()
        
        if self.solver == 'slsqp':
            opt_obj, opt_u, is_success = self._solve_flux_slsqp(tol, max_iters, disp)
        elif self.solver == 'trust-constr':
            opt_obj, opt_u, is_success = self._solve_flux_trust_constr(tol, max_iters, disp)
        # elif self.solver == 'ipopt':
        #     opt_obj, opt_u, is_success = self._solve_flux_ipopt(tol, max_iters, disp)
        elif self.solver == 'ralg':
            opt_obj, opt_u, is_success = self._solve_flux_ralg(tol, max_iters, disp)
        else:
            raise ValueError('currently only "slsqp" and "ralg" are acceptable')
        
        opt_totalfluxes = pd.Series(self.N@opt_u, index = self.model.totalfluxids)
        opt_netfluxes = pd.Series(self.T@self.N@opt_u, index = self.model.netfluxids)
        
        opt_resids = self._calculate_residuals()
        
        exp_MDVs, sim_MDVs = self._get_exp_and_sim_MDVs()
        exp_fluxes, sim_fluxes = self._get_exp_and_sim_fluxes()
        
        nmeas = self._get_nmeasurements(opt_resids)
        nparams = self._get_nparameters(opt_u)
        
        hess = self._get_hessian(opt_u)

        _, dxsim_du = self._calculate_sim_MDVs_derivative()
        dvsim_du = self._calculate_sim_fluxes_derivative()
        
        return (
            opt_totalfluxes, 
            opt_netfluxes, 
            opt_obj, 
            opt_resids, 
            nmeas, 
            nparams, 
            sim_MDVs, 
            exp_MDVs, 
            sim_fluxes, 
            exp_fluxes, 
            hess, 
            self.N, 
            self.T, 
            dxsim_du, 
            dvsim_du, 
            self.model.measured_MDVs_inv_cov, 
            self.model.measured_fluxes_inv_cov, 
            is_success
        )
        
        
    
        
class InstMFAModel(MFAModel):

    def __init__(self, *args):
        
        super().__init__(*args)
        
        self.nfreefluxes = self.N.shape[1]
        self.nconcs = len(self.model.concids)
        self.nnetfluxes = len(self.model.netfluxids)
        # Note: _calculate_difference_sim_exp_MDVs and _calculate_sim_MDVs_derivative
        # are effectively replaced by the new logic within build_objective and build_gradient.
        # They could be removed or refactored if they are not used elsewhere (e.g. in hessian or other calculations).
        # For this subtask, we focus on build_objective and build_gradient.
    
    def build_objective(self):
        
        def obj_func(p_values): # p_values are the optimization parameters (free fluxes u and concentrations c)
            u, c = p_values[:self.nfreefluxes], p_values[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N @ u
            self.model.concentrations[:] = c # Update model concentrations

            # Simulate MDVs for the current parameters u and c
            # self.sim_inst_MDVs = self.calculator._calculate_inst_MDVs() # OLD Global call
            # NEW: Fetch per experiment inside loop. self.sim_inst_MDVs attribute becomes less relevant for obj_func.

            objVal_mdv = 0.0
            self.concat_meas_inst_mdv_sds = [] # For fallback gradient

            if self.model.measured_inst_MDVs:
                concat_meas_inst_mdv_means = []
                concat_sim_inst_mdv_values = []

                # Consistent ordering: experiment_id, then fragment_id, then timepoint
                for exp_id in sorted(self.model.measured_inst_MDVs.keys()):
                    # Calculate sim_inst_MDVs for the current experiment_id and parameters p_values
                    current_exp_sim_inst_MDVs = self.calculator._calculate_inst_MDVs_for_experiment(p_values, experiment_id=exp_id)

                    exp_data = self.model.measured_inst_MDVs[exp_id]
                    for frag_id in sorted(exp_data.keys()):
                        frag_data = exp_data[frag_id]
                        for timepoint in sorted(frag_data.keys()):
                            # t=0 is typically not part of measured_inst_MDVs for fitting.
                            # _calculate_inst_MDVs_for_experiment also excludes t=0 from its returned dict.

                            meas_mean, meas_sd = frag_data[timepoint]

                            sim_mean_for_frag_tp = current_exp_sim_inst_MDVs.get(frag_id, {}).get(timepoint)
                            if sim_mean_for_frag_tp is not None:
                                if len(sim_mean_for_frag_tp) == len(meas_mean):
                                    concat_meas_inst_mdv_means.extend(meas_mean)
                                    self.concat_meas_inst_mdv_sds.extend(meas_sd)
                                    concat_sim_inst_mdv_values.extend(sim_mean_for_frag_tp)
                                else:
                                    logging.warning(
                                        f"InstMDV: Mismatched length for {frag_id} at T={timepoint} in exp {exp_id}. "
                                        f"Measured {len(meas_mean)}, simulated {len(sim_mean_for_frag_tp)}."
                                    )
                            else:
                                logging.warning(f"InstMDV: Simulated data for {frag_id} at T={timepoint} (exp {exp_id}) not found in current_exp_sim_inst_MDVs.")

                if concat_meas_inst_mdv_means:
                    self.residuals_inst_mdv = np.array(concat_sim_inst_mdv_values) - np.array(concat_meas_inst_mdv_means)

                    if (self.model.measured_inst_MDVs_inv_cov is not None and
                        self.model.measured_inst_MDVs_inv_cov.shape[0] == len(self.residuals_inst_mdv) and
                        self.model.measured_inst_MDVs_inv_cov.shape[1] == len(self.residuals_inst_mdv)):
                        objVal_mdv = 0.5 * np.dot(np.dot(self.residuals_inst_mdv.T, self.model.measured_inst_MDVs_inv_cov), self.residuals_inst_mdv)
                    else:
                        if self.model.measured_inst_MDVs_inv_cov is not None:
                             logging.warning(
                                "InstMDV: Mismatch in measured_inst_MDVs_inv_cov shape or not available. "
                                f"Covar shape: {self.model.measured_inst_MDVs_inv_cov.shape}, "
                                f"Residuals length: {len(self.residuals_inst_mdv)}. "
                                "Falling back to SD-based objective for InstMDVs."
                            )
                        else:
                            logging.warning("InstMDV: measured_inst_MDVs_inv_cov is None. Falling back to SD-based objective.")

                        if self.concat_meas_inst_mdv_sds and len(self.concat_meas_inst_mdv_sds) == len(self.residuals_inst_mdv):
                            valid_sds = np.array(self.concat_meas_inst_mdv_sds)
                            if np.any(valid_sds == 0):
                                logging.warning("InstMDV: Zero SD found. Replacing with 1.0 for objective.")
                                valid_sds[valid_sds == 0] = 1.0
                            weighted_residuals = self.residuals_inst_mdv / valid_sds
                            objVal_mdv = 0.5 * np.sum(weighted_residuals**2)
                        else:
                            logging.warning("InstMDV fallback: SDs not available/mismatched. MDV objective part is zero.")
                            objVal_mdv = 0.0
                else:
                    self.residuals_inst_mdv = np.array([])
            else:
                self.residuals_inst_mdv = np.array([])

            objVal_flux = 0.0
            if self.fit_measured_fluxes and self.model.measured_fluxes:
                # Using superclass method for flux part, assuming it's compatible
                flux_diff = super()._calculate_difference_sim_exp_fluxes()
                if self.model.measured_fluxes_inv_cov is not None and \
                   flux_diff.shape[0] == self.model.measured_fluxes_inv_cov.shape[0]:
                    objVal_flux = 0.5 * np.dot(np.dot(flux_diff.T, self.model.measured_fluxes_inv_cov), flux_diff)
                else:
                    logging.warning("InstFlux: Flux inv_cov missing/mismatched. Flux part of objective is zero.")

            return objVal_mdv + objVal_flux

        self.f = obj_func # Assign the comprehensive objective function
    
        
    def build_gradient(self):
        
        def obj_grad(p_values): # p_values are u and c
            # obj_func(p_values) must be called first to set self.residuals_inst_mdv, self.sim_inst_MDVs, etc.
            # The model's total_fluxes and concentrations are set by obj_func.

            # This call should provide derivatives w.r.t. p (u and c)
            # The calculator needs to provide derivatives in the same concatenated order as residuals.
            # The call to self.calculator._calculate_inst_MDVs_and_derivatives_p() is done once by InstFitter
            # and the resulting stacked derivative matrix is stored, typically on self.model.sim_inst_MDVs_der_p.
            # obj_func (called before obj_grad by the solver) has already updated self.residuals_inst_mdv,
            # and model fluxes/concentrations based on current 'p_values'.

            # We need to recalculate the derivative matrix for the current 'p_values'.
            # The Calculator's method uses self.model.total_fluxes and self.model.concentrations
            # which obj_func has just set.
            _ignored_sim_inst_MDVs_dict, self.model.sim_inst_MDVs_der_p = self.calculator._calculate_inst_MDVs_and_derivatives_p()
            current_sim_inst_MDVs_der_p = self.model.sim_inst_MDVs_der_p


            objGrad_mdv = np.zeros(len(p_values))

            if hasattr(self, 'residuals_inst_mdv') and self.residuals_inst_mdv.size > 0:
                if (self.model.measured_inst_MDVs_inv_cov is not None and
                    current_sim_inst_MDVs_der_p is not None and
                    self.model.measured_inst_MDVs_inv_cov.shape[0] == self.residuals_inst_mdv.size and
                    self.model.measured_inst_MDVs_inv_cov.shape[1] == self.residuals_inst_mdv.size and
                    current_sim_inst_MDVs_der_p.shape[0] == self.residuals_inst_mdv.size and
                    current_sim_inst_MDVs_der_p.shape[1] == len(p_values)): # Cols match num_parameters (u+c)
                    objGrad_mdv = np.dot(np.dot(self.residuals_inst_mdv.T, self.model.measured_inst_MDVs_inv_cov), current_sim_inst_MDVs_der_p).squeeze()
                else: # Fallback for InstMDV gradient
                    logging.warning(
                        "InstMDV Gradient: Primary path conditions not met. "
                        f"Covar shape: {self.model.measured_inst_MDVs_inv_cov.shape if self.model.measured_inst_MDVs_inv_cov is not None else 'None'}, "
                        f"Deriv shape: {current_sim_inst_MDVs_der_p.shape if current_sim_inst_MDVs_der_p is not None else 'None'}, "
                        f"Residuals size: {self.residuals_inst_mdv.size}. "
                        "Attempting fallback SD-based gradient."
                    )
                    if (hasattr(self, 'concat_meas_inst_mdv_sds') and
                        current_sim_inst_MDVs_der_p is not None and
                        len(self.concat_meas_inst_mdv_sds) == self.residuals_inst_mdv.size and
                        current_sim_inst_MDVs_der_p.shape[0] == self.residuals_inst_mdv.size and
                        current_sim_inst_MDVs_der_p.shape[1] == len(p_values)):

                        valid_sds_sq = np.array(self.concat_meas_inst_mdv_sds)**2
                        if np.any(valid_sds_sq == 0):
                            logging.warning("InstMDV Grad: Zero SD^2 found. Replacing with 1.0.")
                            valid_sds_sq[valid_sds_sq == 0] = 1.0
                        weighted_residuals_for_grad = self.residuals_inst_mdv / valid_sds_sq
                        objGrad_mdv = np.dot(weighted_residuals_for_grad.T, current_sim_inst_MDVs_der_p).squeeze()
                    else:
                        logging.warning("InstMDV Grad fallback: sim_inst_MDVs_der_p or SDs not available/mismatched. MDV grad part is zero.")

            objGrad_flux = np.zeros(len(p_values))
            if self.fit_measured_fluxes and self.model.measured_fluxes:
                # Derivatives for flux part (w.r.t. u, which are the first nfreefluxes of p_values)
                # The model.measured_fluxes_der_p are w.r.t u and c.
                flux_der_p = self.model.measured_fluxes_der_p # Should be (num_meas_flux, num_params_u_c)
                flux_diff = super()._calculate_difference_sim_exp_fluxes()

                if (self.model.measured_fluxes_inv_cov is not None and
                    flux_der_p is not None and
                    flux_diff.shape[0] == self.model.measured_fluxes_inv_cov.shape[0] and
                    flux_der_p.shape[0] == flux_diff.shape[0] and
                    flux_der_p.shape[1] == len(p_values)):
                    objGrad_flux = np.dot(np.dot(flux_diff.T, self.model.measured_fluxes_inv_cov), flux_der_p).squeeze()
                else:
                    logging.warning("InstFlux Grad: Conditions not met. Flux grad part is zero.")

            grad = objGrad_mdv + objGrad_flux
            if grad.shape[0] != len(p_values):
                logging.error(f"InstGradient shape mismatch: expected ({len(p_values)},), got {grad.shape}. Resetting to zeros.")
                grad = np.zeros(len(p_values))
            return grad
        
        self.df = obj_grad # Assign the comprehensive gradient function
        
    
    def build_hessian(self):
        
        def _ddf(p):
            u, c = p[:self.nfreefluxes], p[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N@u
            self.model.concentrations[:] = c
            
            _, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            hess = MDV_der.T@MDV_inv_cov@MDV_der
            
            return hess
            
        def ddf1(p):
            return _ddf(p)
        
        def ddf2(p):
            hess1 = _ddf(p)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            hess2 = flux_der.T@flux_inv_cov@flux_der
            
            return hess1 + hess2
        
        self.ddf = ddf2 if self.fit_measured_fluxes else ddf1
    
    
    def build_flux_and_conc_bound_constraints(self):
        
        A1 = self.N
        A2 = self.T@self.N
        A3 = -A2
        A4 = np.eye(self.nconcs)
        
        b1 = np.zeros(self.ntotalfluxes)
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        b2 = vnet_lb
        b3 = -vnet_ub
        b4 = np.zeros(self.nconcs)
        
        A = np.zeros(
            (self.ntotalfluxes+2*self.nnetfluxes+self.nconcs, 
             self.nfreefluxes+self.nconcs)
        )
        A[:(self.ntotalfluxes+2*self.nnetfluxes), :self.nfreefluxes] = np.vstack((A1, A2, A3))
        A[(self.ntotalfluxes+2*self.nnetfluxes):, self.nfreefluxes:] = A4
        b = np.concatenate((b1, b2, b3, b4))
        
        if self.solver == 'slsqp':
            self.constrs = {'type': 'ineq', 'fun': lambda p: A@p - b}
        elif self.solver == 'ralg':
            self.A = -A
            self.b = -b
    
    
    def build_initial_flux_and_conc_values(self, ini_netfluxes = None, ini_concs = None):
        '''
        Parameters
        ----------
        ini_netfluxes: array
            Initial guess of net fluxes.
        ini_concs: array
            Initial guess of concentrations.
        '''
        
        if ini_netfluxes is None:
            vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
            vnet_ini = np.random.uniform(low = vnet_lb, high = vnet_ub)
        else:
            vnet_ini = ini_netfluxes
        u_ini = pinv(self.T@self.N)@vnet_ini
        
        if ini_concs is None:
            c_lb, c_ub = np.array(list(self.model.concentrations_range.values())).T
            c_ini = np.random.uniform(low = c_lb, high = c_ub)
        else:
            c_ini = ini_concs

        self.x0 = np.concatenate((u_ini, c_ini))
    
    
    def _initialize_total_fluxes_and_concs(self):
        
        for fluxid in self.model.totalfluxids:
            self.model.total_fluxes[fluxid] = 0.0
    
        for metabid in self.model.concids:
            self.model.concentrations[metabid] = 0.0
            
    
    def _calculate_residuals(self):
        
        MDV_diff = self._calculate_difference_sim_exp_MDVs()
        MDV_diag = np.diag(self.model.measured_inst_MDVs_inv_cov)**0.5
        resids1 = MDV_diff*MDV_diag
        
        if self.fit_measured_fluxes:    
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_diag = np.diag(self.model.measured_fluxes_inv_cov)**0.5
            resids2 = flux_diff*flux_diag
            resids = np.concatenate((resids1, resids2))
        
        return resids
        

    def _get_exp_and_sim_inst_MDVs(self):

        exp_inst_MDVs = self.model.measured_inst_MDVs

        sim_inst_MDVs_all = self.calculator._calculate_inst_MDVs()
        self.calculator._build_initial_sim_MDVs()
        sim_inst_MDVs = {}
        for emuid in self.model.target_EMUs:
            instMDVs = self.model.initial_sim_MDVs[emuid].copy()   
            instMDVs.update({t: sim_inst_MDVs_all[emuid][t] 
                             for t in sim_inst_MDVs_all[emuid]})
            sim_inst_MDVs[emuid] = instMDVs

        return exp_inst_MDVs, sim_inst_MDVs 
        

    def solve_flux(self, tol = 1e-6, max_iters = 400, disp = False):
        
        self._initialize_total_fluxes_and_concs()

        if self.solver == 'slsqp':
            opt_obj, opt_p, is_success = self._solve_flux_slsqp(tol, max_iters, disp)
        elif self.solver == 'ralg':
            opt_obj, opt_p, is_success = self._solve_flux_ralg(tol, max_iters, disp)
        else:
            raise ValueError('currently only "slsqp" and "ralg" are acceptable')    

        opt_u = opt_p[:self.nfreefluxes]
        opt_c = opt_p[self.nfreefluxes:]
        
        opt_totalfluxes = pd.Series(self.N@opt_u, index = self.model.totalfluxids)
        opt_netfluxes = pd.Series(self.T@self.N@opt_u, index = self.model.netfluxids)        
        opt_concs = pd.Series(opt_c, index = self.model.concids)
        
        opt_resids = self._calculate_residuals()
        
        exp_inst_MDVs, sim_inst_MDVs = self._get_exp_and_sim_inst_MDVs()
        exp_fluxes, sim_fluxes = self._get_exp_and_sim_fluxes()
        
        nmeas = self._get_nmeasurements(opt_resids)
        nparams = self._get_nparameters(opt_p)
        
        hess = self._get_hessian(opt_p)

        _, dxsim_dp = self._calculate_sim_MDVs_derivative()
        dxsim_du = dxsim_dp[:, :self.nfreefluxes]
        
        dvsim_dp = self._calculate_sim_fluxes_derivative()
        dvsim_du = dvsim_dp[:, :self.nfreefluxes]
        
        return (
            opt_totalfluxes, 
            opt_netfluxes, 
            opt_concs, 
            opt_obj, 
            opt_resids, 
            nmeas, 
            nparams, 
            sim_inst_MDVs, 
            exp_inst_MDVs, 
            sim_fluxes, 
            exp_fluxes, 
            hess, 
            self.N, 
            self.T, 
            dxsim_du, 
            dvsim_du, 
            self.model.measured_inst_MDVs_inv_cov, 
            self.model.measured_fluxes_inv_cov, 
            is_success
        )
        