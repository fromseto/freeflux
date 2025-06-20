'''Define the MFAModel and InstMFAModel class.'''


__author__ = 'Chao Wu'
__date__ = '05/19/2022'

import numpy as np
import pandas as pd
from scipy.linalg import pinv as scipy_pinv # Renamed to avoid conflict if jax.scipy.linalg.pinv is used locally
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
try:
    from openopt import NLP
except ModuleNotFoundError:
    OPENOPT_INSTALLED = False
else:
    OPENOPT_INSTALLED = True

try:
    import jax
    import jax.numpy as jnp
    from ..utils.jax_utils import core_calculate_mdvs_jax, core_calculate_mdvs_and_derivatives_jax
    JAX_AVAILABLE = True
except ImportError:
    JAX_AVAILABLE = False

from ..utils.utils import Calculator
from functools import partial # For jax.jit static_argnames


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
    
    def __init__(self, model, fit_measured_fluxes, solver = 'slsqp', use_jax=False):
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
        use_jax: bool
            Whether to use JAX for objective and gradient calculation.
        '''
        
        self.model = model
        self.calculator = Calculator(self.model) # Still needed for non-JAX parts or data prep
        self.fit_measured_fluxes = fit_measured_fluxes
        self.solver = solver
        self.use_jax = use_jax and JAX_AVAILABLE
        
        self.N = self.model.null_space
        self.T = self.model.transform_matrix
        
        self.ntotalfluxes = len(self.model.totalfluxids)

        if self.use_jax:
            self._prepare_jax_data()
            self._jit_jax_functions()
        
    def _prepare_jax_data(self):
        """Prepares JAX-compatible static data structures."""
        if not self.use_jax: return

        self.N_jax = jnp.array(self.N)

        # These are expected to be prepared by Fitter.prepare() and stored on model:
        # self.model.matrix_As_jax_static_data, self.model.matrix_Bs_jax_static_data
        # self.model.substrate_MDVs_jax_static_data
        # self.model.matrix_As_der_p_jax_static_data, self.model.matrix_Bs_der_p_jax_static_data
        # self.model.substrate_MDVs_der_p_jax_static_data
        # For now, assume they exist and are pytrees of JAX arrays / callables
        self.matrix_As_jax_static = self.model.matrix_As_jax_static_data
        self.matrix_Bs_jax_static = self.model.matrix_Bs_jax_static_data
        self.substrate_MDVs_jax_static = self.model.substrate_MDVs_jax_static_data

        # Derivatives (for gradient calculation via core_calculate_mdvs_and_derivatives_jax)
        self.matrix_As_der_p_jax_static = self.model.matrix_As_der_p_jax_static_data
        self.matrix_Bs_der_p_jax_static = self.model.matrix_Bs_der_p_jax_static_data
        self.substrate_MDVs_der_p_jax_static = self.model.substrate_MDVs_der_p_jax_static_data

        # Measured MDVs (means)
        self.measured_MDVs_means_jax = {
            k: jnp.array(v[0]) for k, v in self.model.measured_MDVs.items()
        }
        self.measured_MDVs_inv_cov_jax = jnp.array(self.model.measured_MDVs_inv_cov)

        if self.fit_measured_fluxes:
            self.measured_fluxes_means_jax = jnp.array([
                mean for mean, sd in self.model.measured_fluxes.values()
            ])
            # Need a way to map measured flux keys to indices if order matters for residual vector
            self.measured_flux_ids_ordered = list(self.model.measured_fluxes.keys())
            # Assuming model.totalfluxids_map_jax: dict str_id -> int_idx exists
            self.measured_flux_indices_in_total_jax = jnp.array([
                self.model.totalfluxids_map_jax[fid] for fid in self.measured_flux_ids_ordered
            ])
            self.measured_fluxes_inv_cov_jax = jnp.array(self.model.measured_fluxes_inv_cov)
            if hasattr(self.model, 'measured_fluxes_der_p_jax'): # Expected from Fitter.prepare
                 self.measured_fluxes_der_p_jax = self.model.measured_fluxes_der_p_jax
            else: # Fallback or raise error
                 self.measured_fluxes_der_p_jax = jnp.array(self.model.measured_fluxes_der_p)


        self.target_EMU_ids_tuple_jax = tuple(self.model.target_EMUs)
        # EAMs keys are sizes. Ensure model.EAMs_jax_sorted_keys exists from Fitter.prepare()
        self.sorted_emu_sizes_tuple_jax = self.model.EAMs_jax_sorted_keys
        self.num_free_fluxes_jax = self.N.shape[1]

        # Define static argument names for JIT compilation
        # For _core_objective_fn_jax
        self.static_obj_argnames = (
            "N_jax", "matrix_As_static", "matrix_Bs_static",
            "substrate_MDVs_jax_static", "measured_MDVs_means_jax",
            "measured_MDVs_inv_cov_jax", "target_EMU_ids_tuple_jax",
            "sorted_emu_sizes_tuple_jax", "fit_measured_fluxes_static"
        )
        if self.fit_measured_fluxes:
            self.static_obj_argnames += (
                "measured_flux_ids_ordered_static", # For consistent ordering of flux residuals
                "measured_flux_indices_in_total_jax_static",
                "measured_fluxes_means_jax_static",
                "measured_fluxes_inv_cov_jax_static"
            )

        # For _core_gradient_fn_jax (which uses core_calculate_mdvs_and_derivatives_jax)
        # This will be jax.grad of _core_objective_fn_jax
        # Alternatively, if we define a _core_residual_vector_fn_jax, then objective is sum of squares
        # and grad can be derived using that. Let's stick to grad of objective for now.

    def _jit_jax_functions(self):
        if not self.use_jax: return

        # The function to be differentiated by JAX
        # It takes u_jax and all other static data (as JAX arrays/pytrees)
        # and returns a scalar objective value.

        # --- Define the core JAX objective function ---
        def _core_objective_fn_jax(
            u_jax,
            N_jax, matrix_As_static, matrix_Bs_static,
            substrate_MDVs_jax_static, measured_MDVs_means_jax,
            measured_MDVs_inv_cov_jax, target_EMU_ids_tuple_jax,
            sorted_emu_sizes_tuple_jax, fit_measured_fluxes_static,
            # Optional flux args, only used if fit_measured_fluxes_static is True
            measured_flux_ids_ordered_static=None,
            measured_flux_indices_in_total_jax_static=None,
            measured_fluxes_means_jax_static=None,
            measured_fluxes_inv_cov_jax_static=None
            ):

            total_fluxes_jax = N_jax @ u_jax

            sim_MDVs_dict_jax = core_calculate_mdvs_jax(
                total_fluxes_jax,
                matrix_As_static,
                matrix_Bs_static,
                substrate_MDVs_jax_static,
                target_EMU_ids_tuple_jax, # Pass along, though core_calculate_mdvs_jax doesn't currently filter
                sorted_emu_sizes_tuple_jax
            )

            # MDV residuals
            mdv_residuals_list = []
            for emu_id in target_EMU_ids_tuple_jax:
                sim_mdv = sim_MDVs_dict_jax[emu_id]
                exp_mdv = measured_MDVs_means_jax[emu_id]
                mdv_residuals_list.append(sim_mdv - exp_mdv)

            mdv_residuals_vector = jnp.concatenate(mdv_residuals_list)
            obj_mdv = mdv_residuals_vector.T @ measured_MDVs_inv_cov_jax @ mdv_residuals_vector

            obj_total = obj_mdv

            if fit_measured_fluxes_static:
                # Flux residuals
                # Need to select the *simulated* fluxes corresponding to *measured* fluxes
                sim_measured_fluxes_jax = total_fluxes_jax.take(measured_flux_indices_in_total_jax_static)
                flux_residuals_vector = sim_measured_fluxes_jax - measured_fluxes_means_jax_static
                obj_flux = flux_residuals_vector.T @ measured_fluxes_inv_cov_jax_static @ flux_residuals_vector
                obj_total += obj_flux

            return obj_total

        self._jitted_core_objective_fn = jax.jit(
            _core_objective_fn_jax, static_argnames=self.static_obj_argnames
        )

        # Gradient of the core objective function
        self._jitted_core_gradient_fn = jax.jit(
            jax.grad(_core_objective_fn_jax, argnums=0), # Differentiate w.r.t. u_jax (arg 0)
            static_argnames=self.static_obj_argnames
        )

        # Hessian of the core objective function
        self._jitted_core_hessian_fn = jax.jit(
            jax.hessian(_core_objective_fn_jax, argnums=0), # Hessian w.r.t. u_jax (arg 0)
            static_argnames=self.static_obj_argnames
        )
    
    # Make the core objective function a static method or a regular method callable by JAX
    # For simplicity here, it's defined locally in _jit_jax_functions. If it were a class method:
    # @staticmethod (or regular method if it needs to be part of the class for other reasons)
    # def _actual_core_objective_fn_jax(u_jax, N_jax, ...): ...
    # Then in _jit_jax_functions:
    # self._jitted_core_objective_fn = jax.jit(MFAModel._actual_core_objective_fn_jax, ...)
    # self._jitted_core_gradient_fn = jax.jit(jax.grad(MFAModel._actual_core_objective_fn_jax, ...), ...)
    # self._jitted_core_hessian_fn = jax.jit(jax.hessian(MFAModel._actual_core_objective_fn_jax, ...), ...)
    # For now, the local definition within _jit_jax_functions is used.


    # --- Original methods for reference/fallback ---
    def _calculate_difference_sim_exp_MDVs(self):
        # This is the original numpy-based calculation
        simMDVs = self.calculator._calculate_MDVs() # Uses self.model.total_fluxes (numpy)
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] 
                               for emuid in self.model.target_EMUs])
        return diff
            
    def _calculate_difference_sim_exp_fluxes(self):
        # Original numpy-based
        simFluxes = self.model.total_fluxes[list(self.model.measured_fluxes.keys())] # Pandas series indexing
        expFluxes = np.array([mean for mean, _ in self.model.measured_fluxes.values()])
        diff = simFluxes.to_numpy() - expFluxes # Ensure numpy array for operation
        return diff
        
    def _calculate_sim_MDVs_derivative(self):
        # Original numpy-based
        # This sets self.model.total_fluxes, then calls calculator
        simMDVs, simMDVsDer = self.calculator._calculate_MDVs_and_derivatives_p()
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] 
                               for emuid in self.model.target_EMUs])
        dxsim_dp = np.vstack([simMDVsDer[emuid] for emuid in self.model.target_EMUs])
        return diff, dxsim_dp
        
    def _calculate_sim_fluxes_derivative(self):
        # Original numpy-based
        dvsim_dp = self.model.measured_fluxes_der_p # This is d(total_fluxes)/du * N
        return dvsim_dp

    # --- Build objective and gradient using JAX if enabled ---
    def build_objective(self):
        if self.use_jax:
            # Prepare static args for the jitted function call
            static_args = {
                "N_jax": self.N_jax,
                "matrix_As_static": self.matrix_As_jax_static,
                "matrix_Bs_static": self.matrix_Bs_jax_static,
                "substrate_MDVs_jax_static": self.substrate_MDVs_jax_static,
                "measured_MDVs_means_jax": self.measured_MDVs_means_jax,
                "measured_MDVs_inv_cov_jax": self.measured_MDVs_inv_cov_jax,
                "target_EMU_ids_tuple_jax": self.target_EMU_ids_tuple_jax,
                "sorted_emu_sizes_tuple_jax": self.sorted_emu_sizes_tuple_jax,
                "fit_measured_fluxes_static": self.fit_measured_fluxes
            }
            if self.fit_measured_fluxes:
                static_args.update({
                    "measured_flux_ids_ordered_static": self.measured_flux_ids_ordered,
                    "measured_flux_indices_in_total_jax_static": self.measured_flux_indices_in_total_jax,
                    "measured_fluxes_means_jax_static": self.measured_fluxes_means_jax,
                    "measured_fluxes_inv_cov_jax_static": self.measured_fluxes_inv_cov_jax
                })

            def f_jax_wrapper(u_numpy):
                u_jax = jnp.array(u_numpy)
                # JAX operations happen inside _jitted_core_objective_fn
                obj_val_jax = self._jitted_core_objective_fn(u_jax, **static_args)
                return float(obj_val_jax)
            self.f = f_jax_wrapper
        else: # Original numpy-based objective
            def _f_numpy(u_np):
                self.model.total_fluxes.iloc[:] = self.N @ u_np # Update pandas Series

                MDV_diff = self._calculate_difference_sim_exp_MDVs() # Uses self.model.total_fluxes
                MDV_inv_cov = self.model.measured_MDVs_inv_cov
                obj = MDV_diff @ MDV_inv_cov @ MDV_diff
                return obj
            
            def f1_numpy(u_np):
                return _f_numpy(u_np)

            def f2_numpy(u_np):
                obj1 = _f_numpy(u_np)
                # flux_diff uses self.model.total_fluxes which was updated by _f_numpy
                flux_diff = self._calculate_difference_sim_exp_fluxes()
                flux_inv_cov = self.model.measured_fluxes_inv_cov
                obj2 = flux_diff @ flux_inv_cov @ flux_diff
                return obj1 + obj2

            self.f = f2_numpy if self.fit_measured_fluxes else f1_numpy
        
    def build_gradient(self):
        if self.use_jax:
            # Static args are the same as for the objective
            static_args = {
                "N_jax": self.N_jax,
                "matrix_As_static": self.matrix_As_jax_static,
                "matrix_Bs_static": self.matrix_Bs_jax_static,
                "substrate_MDVs_jax_static": self.substrate_MDVs_jax_static,
                "measured_MDVs_means_jax": self.measured_MDVs_means_jax,
                "measured_MDVs_inv_cov_jax": self.measured_MDVs_inv_cov_jax,
                "target_EMU_ids_tuple_jax": self.target_EMU_ids_tuple_jax,
                "sorted_emu_sizes_tuple_jax": self.sorted_emu_sizes_tuple_jax,
                "fit_measured_fluxes_static": self.fit_measured_fluxes
            }
            if self.fit_measured_fluxes:
                 static_args.update({
                    "measured_flux_ids_ordered_static": self.measured_flux_ids_ordered,
                    "measured_flux_indices_in_total_jax_static": self.measured_flux_indices_in_total_jax,
                    "measured_fluxes_means_jax_static": self.measured_fluxes_means_jax,
                    "measured_fluxes_inv_cov_jax_static": self.measured_fluxes_inv_cov_jax
                })

            def df_jax_wrapper(u_numpy):
                u_jax = jnp.array(u_numpy)
                grad_val_jax = self._jitted_core_gradient_fn(u_jax, **static_args)
                return np.array(grad_val_jax)
            self.df = df_jax_wrapper
        else: # Original numpy-based gradient
            # Important: self.model.total_fluxes must be set before calling derivative funcs
            def _df_numpy(u_np):
                # This state modification is problematic for purity if we were to JIT this part directly.
                # Here, it's part of the non-JAX path.
                self.model.total_fluxes.iloc[:] = self.N @ u_np

                MDV_diff, MDV_der = self._calculate_sim_MDVs_derivative() # Uses self.model.total_fluxes
                MDV_inv_cov = self.model.measured_MDVs_inv_cov
                grad = MDV_der.T @ MDV_inv_cov @ MDV_diff
                return grad

            def df1_numpy(u_np):
                return _df_numpy(u_np)
            
            def df2_numpy(u_np):
                # _df_numpy updates self.model.total_fluxes
                grad1 = _df_numpy(u_np)

                # These use the updated self.model.total_fluxes
                flux_der = self._calculate_sim_fluxes_derivative()
                flux_diff = self._calculate_difference_sim_exp_fluxes()
                flux_inv_cov = self.model.measured_fluxes_inv_cov
                grad2 = flux_der.T @ flux_inv_cov @ flux_diff
                return grad1 + grad2
            
            self.df = df2_numpy if self.fit_measured_fluxes else df1_numpy

    
    def build_hessian(self):
        if self.use_jax:
            # Prepare static args for the jitted Hessian function call
            static_args = {
                "N_jax": self.N_jax,
                "matrix_As_static": self.matrix_As_jax_static,
                "matrix_Bs_static": self.matrix_Bs_jax_static,
                "substrate_MDVs_jax_static": self.substrate_MDVs_jax_static,
                "measured_MDVs_means_jax": self.measured_MDVs_means_jax,
                "measured_MDVs_inv_cov_jax": self.measured_MDVs_inv_cov_jax,
                "target_EMU_ids_tuple_jax": self.target_EMU_ids_tuple_jax,
                "sorted_emu_sizes_tuple_jax": self.sorted_emu_sizes_tuple_jax,
                "fit_measured_fluxes_static": self.fit_measured_fluxes
            }
            if self.fit_measured_fluxes:
                 static_args.update({
                    "measured_flux_ids_ordered_static": self.measured_flux_ids_ordered,
                    "measured_flux_indices_in_total_jax_static": self.measured_flux_indices_in_total_jax,
                    "measured_fluxes_means_jax_static": self.measured_fluxes_means_jax,
                    "measured_fluxes_inv_cov_jax_static": self.measured_fluxes_inv_cov_jax
                })

            def ddf_jax_wrapper(u_numpy):
                u_jax = jnp.array(u_numpy)
                # Ensure _jitted_core_hessian_fn is available (defined in _jit_jax_functions)
                hess_val_jax = self._jitted_core_hessian_fn(u_jax, **static_args)
                return np.array(hess_val_jax)
            
            # Provide Hessian if solver is not SLSQP or ralg (which don't typically use it directly from user)
            # or if explicitly configured to use it.
            if self.solver not in ['slsqp', 'ralg']: # e.g. for 'trust-constr'
                self.ddf = ddf_jax_wrapper
            else:
                self.ddf = None # SLSQP and ralg can work without it or use approximations.
        else:
            self._build_original_hessian()

    def _build_original_hessian(self):
        # This is the original logic extracted
        def _ddf_numpy(u_np):
            # Critical: self.model.total_fluxes must be set based on u_np
            # This happens if called after obj/grad in scipy, but direct call needs care.
            # For safety, recalculate here if not sure about state.
            # However, original code implies it's called in a context where total_fluxes is current.
            # self.model.total_fluxes.iloc[:] = self.N @ u_np # Ensure state for original calc

            _, MDV_der = self._calculate_sim_MDVs_derivative() # Uses current self.model.total_fluxes
            MDV_inv_cov = self.model.measured_MDVs_inv_cov
            hess = MDV_der.T @ MDV_inv_cov @ MDV_der
            return hess
            
        def ddf1_numpy(u_np):
            return _ddf_numpy(u_np)
        
        def ddf2_numpy(u_np):
            hess1 = _ddf_numpy(u_np)
            flux_der = self._calculate_sim_fluxes_derivative() # Uses current self.model.total_fluxes
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            hess2 = flux_der.T @ flux_inv_cov @ flux_der
            return hess1 + hess2
        
        self.ddf = ddf2_numpy if self.fit_measured_fluxes else ddf1_numpy
    
            
    def build_flux_bound_constraints(self):
        
        A1 = self.N
        A2 = self.T @ self.N
        A3 = -A2
        
        b1 = np.zeros(self.ntotalfluxes)
        # Ensure net_fluxes_range values are numpy arrays for consistency
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        b2 = vnet_lb
        b3 = -vnet_ub
        
        A = np.vstack((A1, A2, A3))
        b = np.concatenate((b1, b2, b3))

        if self.solver == 'slsqp':
            # SLSQP constraints fun must return a 1D array
            self.constrs = {'type': 'ineq', 'fun': lambda u: (A @ u - b).ravel()}
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

    def __init__(self, model, fit_measured_fluxes, solver='slsqp', use_jax=False):
        # Call MFAModel's __init__ but without use_jax, as InstMFAModel handles its own JAX setup
        # Or, pass use_jax=False explicitly if MFAModel.__init__ uses it.
        # MFAModel.__init__ was: def __init__(self, model, fit_measured_fluxes, solver = 'slsqp', use_jax=False):
        super().__init__(model, fit_measured_fluxes, solver, use_jax=False) # Initialize base non-JAX parts
        
        self.use_jax_inst = use_jax and JAX_AVAILABLE # Specific flag for instationary JAX
        
        self.nfreefluxes = self.N.shape[1]
        self.nconcs = len(self.model.concids if hasattr(self.model, 'concids') else [])
        self.nnetfluxes = len(self.model.netfluxids)

        if self.use_jax_inst:
            if not (hasattr(self.model, 'jax_prepared_inst') and self.model.jax_prepared_inst):
                # This might indicate that InstFitter.prepare with JAX options was not called.
                warnings.warn("InstMFAModel initialized with use_jax=True, but JAX data for instationary model seems unprepared. JAX path may fail.")
            self._prepare_jax_data_inst()
            self._jit_jax_functions_inst()

    def _prepare_jax_data_inst(self):
        """Prepares JAX-compatible static data for InstMFAModel."""
        if not self.use_jax_inst: return

        # Data common with MFAModel, ensure it's JAXified if not already by superclass or here
        if not hasattr(self, 'N_jax'): # Could be set by MFAModel if its _prepare_jax_data was called
            self.N_jax = jnp.array(self.N)

        # Instationary specific JAX data structures from model (prepared by Calculator/InstFitter)
        self.matrix_As_jax_static_inst = self.model.matrix_As_jax_static_data
        self.matrix_Bs_jax_static_inst = self.model.matrix_Bs_jax_static_data
        self.matrix_Ms_jax_static_inst = self.model.matrix_Ms_jax_static_data

        self.substrate_MDVs_jax_static_inst = self.model.substrate_MDVs_jax_static_data
        
        self.matrix_As_der_p_jax_static_inst = self.model.matrix_As_der_p_jax_static_data # d/dp where p=(u,c)
        self.matrix_Bs_der_p_jax_static_inst = self.model.matrix_Bs_der_p_jax_static_data
        self.matrix_Ms_der_p_jax_static_inst = self.model.matrix_Ms_der_p_jax_static_data
        self.substrate_MDVs_der_p_jax_static_inst = self.model.substrate_MDVs_der_p_jax_static_data

        self.initial_Xs_jax_inst = self.model.initial_matrix_Xs_jax
        self.initial_Ys_jax_inst = self.model.initial_matrix_Ys_jax
        self.initial_Xs_der_p_jax_inst = self.model.initial_matrix_Xs_der_p_jax
        self.initial_Ys_der_p_jax_inst = self.model.initial_matrix_Ys_der_p_jax

        # Measured instMDVs (means) - complex structure {emu_id: {time: array}}
        # For JAX, this needs to be a pytree of JAX arrays.
        self.measured_inst_MDVs_means_jax = {
            emu_id: {t: jnp.array(mdv_data[0]) for t, mdv_data in time_data.items()}
            for emu_id, time_data in self.model.measured_inst_MDVs.items()
        }
        self.measured_inst_MDVs_inv_cov_jax = jnp.array(self.model.measured_inst_MDVs_inv_cov)
        self.timepoints_jax = jnp.array(self.model.timepoints)

        if self.fit_measured_fluxes:
            self.measured_fluxes_means_jax_inst = jnp.array([
                mean for mean, sd in self.model.measured_fluxes.values()
            ])
            self.measured_flux_ids_ordered_inst = list(self.model.measured_fluxes.keys())
            self.measured_flux_indices_in_total_jax_inst = jnp.array([
                self.model.totalfluxids_map_jax[fid] for fid in self.measured_flux_ids_ordered_inst
            ])
            self.measured_fluxes_inv_cov_jax_inst = jnp.array(self.model.measured_fluxes_inv_cov)
            # d(flux)/dp where p=(u,c)
            self.measured_fluxes_der_p_jax_inst = self.model.measured_fluxes_der_p_jax # Should be (n_params_inst, n_meas_fluxes)

        self.target_EMU_ids_tuple_jax_inst = tuple(self.model.target_EMUs)
        self.sorted_emu_sizes_tuple_jax_inst = self.model.EAMs_jax_sorted_keys
        self.num_total_params_inst = self.nfreefluxes + self.nconcs


        # Define static argument names for JIT compilation of instationary core objective
        self.static_obj_argnames_inst = (
            "N_jax", "matrix_As_static_inst", "matrix_Bs_static_inst", "matrix_Ms_static_inst",
            "substrate_MDVs_jax_static_inst",
            "initial_Xs_jax_inst", "initial_Ys_jax_inst",
            "measured_inst_MDVs_means_jax", "measured_inst_MDVs_inv_cov_jax",
            "timepoints_jax", "target_EMU_ids_tuple_jax_inst", "sorted_emu_sizes_tuple_jax_inst",
            "fit_measured_fluxes_static", "nfreefluxes_static" # nfreefluxes needed to split p into u,c
        )
        if self.fit_measured_fluxes:
            self.static_obj_argnames_inst += (
                "measured_flux_ids_ordered_static_inst",
                "measured_flux_indices_in_total_jax_static_inst",
                "measured_fluxes_means_jax_static_inst",
                "measured_fluxes_inv_cov_jax_static_inst"
            )

    def _jit_jax_functions_inst(self):
        if not self.use_jax_inst: return

        def _core_objective_fn_inst_jax(
            p_jax, # Concatenated [u_jax, c_jax]
            # Static args start here
            N_jax, matrix_As_static_inst, matrix_Bs_static_inst, matrix_Ms_static_inst,
            substrate_MDVs_jax_static_inst,
            initial_Xs_jax_inst, initial_Ys_jax_inst,
            measured_inst_MDVs_means_jax, measured_inst_MDVs_inv_cov_jax,
            timepoints_jax, target_EMU_ids_tuple_jax_inst, sorted_emu_sizes_tuple_jax_inst,
            fit_measured_fluxes_static, nfreefluxes_static,
            # Optional static flux args
            measured_flux_ids_ordered_static_inst=None,
            measured_flux_indices_in_total_jax_static_inst=None,
            measured_fluxes_means_jax_static_inst=None,
            measured_fluxes_inv_cov_jax_static_inst=None
        ):
            u_jax = p_jax[:nfreefluxes_static]
            c_jax = p_jax[nfreefluxes_static:]
            total_fluxes_jax = N_jax @ u_jax

            # This function is currently a placeholder in jax_utils.py
            sim_inst_MDVs_dict_jax = core_calculate_inst_mdvs_jax(
                initial_Xs_jax_inst, initial_Ys_jax_inst, timepoints_jax,
                total_fluxes_jax, c_jax,
                matrix_As_static_inst, matrix_Bs_static_inst, matrix_Ms_static_inst,
                substrate_MDVs_jax_static_inst,
                target_EMU_ids_tuple_jax_inst, sorted_emu_sizes_tuple_jax_inst
            )

            mdv_residuals_list = []
            for emu_id in target_EMU_ids_tuple_jax_inst:
                if emu_id in sim_inst_MDVs_dict_jax:
                    for t in timepoints_jax: # Iterate over all timepoints
                        if t == 0: continue # Skip t=0 for residuals usually
                        if t in sim_inst_MDVs_dict_jax[emu_id] and \
                           emu_id in measured_inst_MDVs_means_jax and \
                           t in measured_inst_MDVs_means_jax[emu_id]:

                            sim_mdv_t = sim_inst_MDVs_dict_jax[emu_id][t]
                            exp_mdv_t = measured_inst_MDVs_means_jax[emu_id][t]
                            mdv_residuals_list.append(sim_mdv_t - exp_mdv_t)

            if not mdv_residuals_list: # Should not happen if there's measured data
                obj_mdv = 0.0
            else:
                mdv_residuals_vector = jnp.concatenate(mdv_residuals_list)
                obj_mdv = mdv_residuals_vector.T @ measured_inst_MDVs_inv_cov_jax @ mdv_residuals_vector

            obj_total = obj_mdv

            if fit_measured_fluxes_static:
                sim_measured_fluxes_jax = total_fluxes_jax.take(measured_flux_indices_in_total_jax_static_inst)
                flux_residuals_vector = sim_measured_fluxes_jax - measured_fluxes_means_jax_static_inst
                obj_flux = flux_residuals_vector.T @ measured_fluxes_inv_cov_jax_static_inst @ flux_residuals_vector
                obj_total += obj_flux

            return obj_total

        self._jitted_core_objective_fn_inst = jax.jit(
            _core_objective_fn_inst_jax, static_argnames=self.static_obj_argnames_inst
        )

        self._jitted_core_gradient_fn_inst = jax.jit(
            jax.grad(_core_objective_fn_inst_jax, argnums=0), # Differentiate w.r.t. p_jax
            static_argnames=self.static_obj_argnames_inst
        )

        self._jitted_core_hessian_fn_inst = jax.jit(
            jax.hessian(_core_objective_fn_inst_jax, argnums=0), # Hessian w.r.t. p_jax
            static_argnames=self.static_obj_argnames_inst
        )

    # --- Original methods for InstMFAModel ---
    def _calculate_difference_sim_exp_MDVs(self):
        # This is the original numpy-based calculation for instationary
        # It updates self.model.total_fluxes and self.model.concentrations first
        # then calls self.calculator._calculate_inst_MDVs()
        # For this method to be called from original build_objective, p has to be split.
        # This method is specific to the non-JAX path.

        # The p argument is not passed here, assumes self.model attributes are set.
        # This is how original InstMFAModel._f calls it.
        simMDVs_all_emus_all_times = self.calculator._calculate_inst_MDVs() # numpy based
        expMDVs_all_emus_all_times = self.model.measured_inst_MDVs

        diff_list = []
        for emuid in self.model.target_EMUs:
            if emuid in simMDVs_all_emus_all_times and emuid in expMDVs_all_emus_all_times:
                sim_data_for_emu = simMDVs_all_emus_all_times[emuid]
                exp_data_for_emu = expMDVs_all_emus_all_times[emuid]
                for t_point in exp_data_for_emu: # Iterate over measured timepoints
                    if t_point != 0 and t_point in sim_data_for_emu:
                        diff_list.append(sim_data_for_emu[t_point] - exp_data_for_emu[t_point][0])

        if not diff_list: return np.array([]) # Handle case with no valid residuals
        return np.concatenate(diff_list)


    def _calculate_sim_MDVs_derivative(self):
        # Original numpy-based for instationary. Assumes model state is set.
        simMDVs_all_times, simMDVsDer_all_times = self.calculator._calculate_inst_MDVs_and_derivatives_p()
        expMDVs_all_times = self.model.measured_inst_MDVs

        diff_list = []
        dxsim_dp_list = []

        for emuid in self.model.target_EMUs:
            if emuid in simMDVs_all_times and emuid in expMDVs_all_times:
                sim_data_for_emu = simMDVs_all_times[emuid]
                exp_data_for_emu = expMDVs_all_times[emuid]
                sim_der_for_emu = simMDVsDer_all_times[emuid] # {time: deriv_array (n_coeffs, n_params)}

                for t_point in exp_data_for_emu:
                    if t_point != 0 and t_point in sim_data_for_emu and t_point in sim_der_for_emu:
                        diff_list.append(sim_data_for_emu[t_point] - exp_data_for_emu[t_point][0])
                        # Derivative from calculator is (n_coeffs, n_params). For vstack, it's fine.
                        dxsim_dp_list.append(sim_der_for_emu[t_point])

        if not diff_list: return np.array([]), np.array([]).reshape(0, self.nfreefluxes + self.nconcs) # Adjust shape for empty

        diff_vector = np.concatenate(diff_list)
        # dxsim_dp_stacked will be (total_coeffs_all_times, n_params)
        dxsim_dp_stacked = np.vstack(dxsim_dp_list) if dxsim_dp_list else np.array([]).reshape(0, self.nfreefluxes + self.nconcs)

        return diff_vector, dxsim_dp_stacked


    def build_objective(self):
        if self.use_jax_inst:
            static_args_inst = {
                "N_jax": self.N_jax, # Assuming this is prepared correctly
                "matrix_As_static_inst": self.matrix_As_jax_static_inst,
                "matrix_Bs_static_inst": self.matrix_Bs_jax_static_inst,
                "matrix_Ms_static_inst": self.matrix_Ms_jax_static_inst,
                "substrate_MDVs_jax_static_inst": self.substrate_MDVs_jax_static_inst,
                "initial_Xs_jax_inst": self.initial_Xs_jax_inst,
                "initial_Ys_jax_inst": self.initial_Ys_jax_inst,
                "measured_inst_MDVs_means_jax": self.measured_inst_MDVs_means_jax,
                "measured_inst_MDVs_inv_cov_jax": self.measured_inst_MDVs_inv_cov_jax,
                "timepoints_jax": self.timepoints_jax,
                "target_EMU_ids_tuple_jax_inst": self.target_EMU_ids_tuple_jax_inst,
                "sorted_emu_sizes_tuple_jax_inst": self.sorted_emu_sizes_tuple_jax_inst,
                "fit_measured_fluxes_static": self.fit_measured_fluxes,
                "nfreefluxes_static": self.nfreefluxes
            }
            if self.fit_measured_fluxes:
                static_args_inst.update({
                    "measured_flux_ids_ordered_static_inst": self.measured_flux_ids_ordered_inst,
                    "measured_flux_indices_in_total_jax_static_inst": self.measured_flux_indices_in_total_jax_inst,
                    "measured_fluxes_means_jax_static_inst": self.measured_fluxes_means_jax_inst,
                    "measured_fluxes_inv_cov_jax_static_inst": self.measured_fluxes_inv_cov_jax_inst
                })

            def f_jax_inst_wrapper(p_numpy):
                p_jax = jnp.array(p_numpy)
                obj_val_jax = self._jitted_core_objective_fn_inst(p_jax, **static_args_inst)
                return float(obj_val_jax)
            self.f = f_jax_inst_wrapper
        else: # Original numpy-based objective for InstMFAModel
            def _f_inst_numpy(p_np): # p_np is [u_np, c_np]
                u_np, c_np = p_np[:self.nfreefluxes], p_np[self.nfreefluxes:]
                # Update model state for calculator methods
                self.model.total_fluxes.iloc[:] = self.N @ u_np
                self.model.concentrations.iloc[:] = c_np # Assuming concentrations is a pandas Series

                MDV_diff = self._calculate_difference_sim_exp_MDVs() # Uses updated model state
                if MDV_diff.size == 0: return 0.0 # No residuals to compute objective from
                MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov # NumPy array
                obj = MDV_diff @ MDV_inv_cov @ MDV_diff
                return obj

            def f1_inst_numpy(p_np):
                return _f_inst_numpy(p_np)

            def f2_inst_numpy(p_np):
                obj1 = _f_inst_numpy(p_np) # This also updates model state

                # _calculate_difference_sim_exp_fluxes uses self.model.total_fluxes
                flux_diff = super()._calculate_difference_sim_exp_fluxes() # Use MFAModel's method
                if flux_diff.size == 0: return obj1 # No flux data to fit
                flux_inv_cov = self.model.measured_fluxes_inv_cov # NumPy array
                obj2 = flux_diff @ flux_inv_cov @ flux_diff
                return obj1 + obj2

            self.f = f2_inst_numpy if self.fit_measured_fluxes else f1_inst_numpy


    def build_gradient(self):
        if self.use_jax_inst:
            static_args_inst = { # Same as for objective
                "N_jax": self.N_jax,
                "matrix_As_static_inst": self.matrix_As_jax_static_inst,
                "matrix_Bs_static_inst": self.matrix_Bs_jax_static_inst,
                "matrix_Ms_static_inst": self.matrix_Ms_jax_static_inst,
                "substrate_MDVs_jax_static_inst": self.substrate_MDVs_jax_static_inst,
                "initial_Xs_jax_inst": self.initial_Xs_jax_inst,
                "initial_Ys_jax_inst": self.initial_Ys_jax_inst,
                "measured_inst_MDVs_means_jax": self.measured_inst_MDVs_means_jax,
                "measured_inst_MDVs_inv_cov_jax": self.measured_inst_MDVs_inv_cov_jax,
                "timepoints_jax": self.timepoints_jax,
                "target_EMU_ids_tuple_jax_inst": self.target_EMU_ids_tuple_jax_inst,
                "sorted_emu_sizes_tuple_jax_inst": self.sorted_emu_sizes_tuple_jax_inst,
                "fit_measured_fluxes_static": self.fit_measured_fluxes,
                "nfreefluxes_static": self.nfreefluxes
            }
            if self.fit_measured_fluxes:
                static_args_inst.update({
                    "measured_flux_ids_ordered_static_inst": self.measured_flux_ids_ordered_inst,
                    "measured_flux_indices_in_total_jax_static_inst": self.measured_flux_indices_in_total_jax_inst,
                    "measured_fluxes_means_jax_static_inst": self.measured_fluxes_means_jax_inst,
                    "measured_fluxes_inv_cov_jax_static_inst": self.measured_fluxes_inv_cov_jax_inst
                })

            def df_jax_inst_wrapper(p_numpy):
                p_jax = jnp.array(p_numpy)
                grad_val_jax = self._jitted_core_gradient_fn_inst(p_jax, **static_args_inst)
                return np.array(grad_val_jax)
            self.df = df_jax_inst_wrapper
        else: # Original numpy-based gradient for InstMFAModel
            def _df_inst_numpy(p_np): # p_np is [u_np, c_np]
                u_np, c_np = p_np[:self.nfreefluxes], p_np[self.nfreefluxes:]
                # Update model state
                self.model.total_fluxes.iloc[:] = self.N @ u_np
                self.model.concentrations.iloc[:] = c_np

                MDV_diff, MDV_der = self._calculate_sim_MDVs_derivative() # Uses updated model state
                                                                        # MDV_der shape (total_coeffs, n_params_inst)
                if MDV_diff.size == 0: return np.zeros_like(p_np)

                MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
                grad = MDV_der.T @ MDV_inv_cov @ MDV_diff # (n_params_inst, total_coeffs) @ (...) @ (total_coeffs,) -> (n_params_inst,)
                return grad

            def df1_inst_numpy(p_np):
                return _df_inst_numpy(p_np)

            def df2_inst_numpy(p_np):
                grad1 = _df_inst_numpy(p_np) # Also updates model state

                # _calculate_sim_fluxes_derivative is from MFAModel, uses self.model.measured_fluxes_der_p
                # This measured_fluxes_der_p should be d(flux)/dp where p=(u,c) for instationary.
                # It was set by Calculator._calculate_measured_fluxes_derivative_p('inst')
                flux_der = super()._calculate_sim_fluxes_derivative() # (n_meas_fluxes, n_params_inst)
                flux_diff = super()._calculate_difference_sim_exp_fluxes() # (n_meas_fluxes,)

                if flux_diff.size == 0: return grad1

                flux_inv_cov = self.model.measured_fluxes_inv_cov
                grad2 = flux_der.T @ flux_inv_cov @ flux_diff # (n_params_inst, n_meas_fluxes) @ (...) @ (n_meas_fluxes,) -> (n_params_inst,)
                return grad1 + grad2

            self.df = df2_inst_numpy if self.fit_measured_fluxes else df1_inst_numpy

    # build_hessian for InstMFAModel would follow similar logic if JAX path is taken
    # For now, it will inherit MFAModel's build_hessian.
    # If JAX is used for InstMFAModel, MFAModel's build_hessian might try to use
    # steady-state JAX data if not careful.
    # Override build_hessian for InstMFAModel:
    def build_hessian(self):
        if self.use_jax_inst:
            static_args_inst = { # Same static args as objective/gradient for instationary
                "N_jax": self.N_jax,
                "matrix_As_static_inst": self.matrix_As_jax_static_inst,
                "matrix_Bs_static_inst": self.matrix_Bs_jax_static_inst,
                "matrix_Ms_static_inst": self.matrix_Ms_jax_static_inst,
                "substrate_MDVs_jax_static_inst": self.substrate_MDVs_jax_static_inst,
                "initial_Xs_jax_inst": self.initial_Xs_jax_inst,
                "initial_Ys_jax_inst": self.initial_Ys_jax_inst,
                "measured_inst_MDVs_means_jax": self.measured_inst_MDVs_means_jax,
                "measured_inst_MDVs_inv_cov_jax": self.measured_inst_MDVs_inv_cov_jax,
                "timepoints_jax": self.timepoints_jax,
                "target_EMU_ids_tuple_jax_inst": self.target_EMU_ids_tuple_jax_inst,
                "sorted_emu_sizes_tuple_jax_inst": self.sorted_emu_sizes_tuple_jax_inst,
                "fit_measured_fluxes_static": self.fit_measured_fluxes,
                "nfreefluxes_static": self.nfreefluxes
            }
            if self.fit_measured_fluxes:
                static_args_inst.update({
                    "measured_flux_ids_ordered_static_inst": self.measured_flux_ids_ordered_inst,
                    "measured_flux_indices_in_total_jax_static_inst": self.measured_flux_indices_in_total_jax_inst,
                    "measured_fluxes_means_jax_static_inst": self.measured_fluxes_means_jax_inst,
                    "measured_fluxes_inv_cov_jax_static_inst": self.measured_fluxes_inv_cov_jax_inst
                })

            def ddf_jax_inst_wrapper(p_numpy):
                p_jax = jnp.array(p_numpy)
                hess_val_jax = self._jitted_core_hessian_fn_inst(p_jax, **static_args_inst)
                return np.array(hess_val_jax)

            if self.solver not in ['slsqp', 'ralg']:
                self.ddf = ddf_jax_inst_wrapper
            else:
                self.ddf = None
        else:
            self._build_original_hessian_inst()

    def _build_original_hessian_inst(self):
        # Original Hessian logic for InstMFAModel
        def _ddf_inst_numpy(p_np):
            u_np, c_np = p_np[:self.nfreefluxes], p_np[self.nfreefluxes:]
            self.model.total_fluxes.iloc[:] = self.N @ u_np
            self.model.concentrations.iloc[:] = c_np

            _, MDV_der = self._calculate_sim_MDVs_derivative() # (total_coeffs, n_params_inst)
            if MDV_der.size == 0: return np.zeros((len(p_np), len(p_np)))

            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            hess = MDV_der.T @ MDV_inv_cov @ MDV_der # Gauss-Newton part
            return hess

        def ddf1_inst_numpy(p_np):
            return _ddf_inst_numpy(p_np)

        def ddf2_inst_numpy(p_np):
            hess1 = _ddf_inst_numpy(p_np) # Also updates model state

            flux_der = super()._calculate_sim_fluxes_derivative() # (n_meas_fluxes, n_params_inst)
            if flux_der.size == 0: return hess1

            flux_inv_cov = self.model.measured_fluxes_inv_cov
            hess2 = flux_der.T @ flux_inv_cov @ flux_der
            return hess1 + hess2

        self.ddf = ddf2_inst_numpy if self.fit_measured_fluxes else ddf1_inst_numpy
    
    def _calculate_difference_sim_exp_MDVs(self):
        
        simMDVs = self.calculator._calculate_inst_MDVs()
        expMDVs = self.model.measured_inst_MDVs
        diff = np.concatenate([simMDVs[emuid][t] - expMDVs[emuid][t][0] 
                               for emuid in self.model.target_EMUs
                               for t in expMDVs[emuid] if t != 0])
        
        return diff
        
        
    def _calculate_sim_MDVs_derivative(self):
        
        simMDVs, simMDVsDer = self.calculator._calculate_inst_MDVs_and_derivatives_p()
        expMDVs = self.model.measured_inst_MDVs
        diff = np.concatenate([simMDVs[emuid][t] - expMDVs[emuid][t][0] 
                               for emuid in self.model.target_EMUs
                               for t in expMDVs[emuid] if t != 0])
        
        dxsim_dp = np.vstack([simMDVsDer[emuid][t] 
                              for emuid in self.model.target_EMUs
                              for t in expMDVs[emuid] if t != 0])
        
        return diff, dxsim_dp
        
    
    def build_objective(self):
        
        def _f(p):
            u, c = p[:self.nfreefluxes], p[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N@u
            self.model.concentrations[:] = c
            
            MDV_diff = self._calculate_difference_sim_exp_MDVs()
            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            
            obj = MDV_diff@MDV_inv_cov@MDV_diff
            
            return obj
        
        def f1(p):
            return _f(p)
            
        def f2(p):
            obj1 = _f(p)
            
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            obj2 = flux_diff@flux_inv_cov@flux_diff
            
            return obj1 + obj2
            
        self.f = f2 if self.fit_measured_fluxes else f1
    
        
    def build_gradient(self):
        
        def _df(p):
            u, c = p[:self.nfreefluxes], p[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N@u
            self.model.concentrations[:] = c
            
            MDV_diff, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            grad = MDV_der.T@MDV_inv_cov@MDV_diff
            
            return grad
            
        def df1(p):
            return _df(p)
            
        def df2(p):
            grad1 = _df(p)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            grad2 = flux_der.T@flux_inv_cov@flux_diff
            
            return grad1 + grad2
        
        self.df = df2 if self.fit_measured_fluxes else df1
        
    
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
        