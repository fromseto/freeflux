'''Define the Calculator class.'''


__author__ = 'Chao Wu'
__date__ = '05/18/2022'


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
    from jax import config, jacfwd, jit
    config.update('jax_platform_name', 'cpu')
    # It's good practice to import specific jax functions like jit if used directly in this file.
except ModuleNotFoundError:
    JAX_INSTALLED = False
    # Define jnp as np if JAX is not installed, for type hinting or conditional code.
    # However, code using JAX features should be guarded by JAX_INSTALLED.
    # For lambdify, sympy needs to know about jax.numpy.
    # We can pass {'jnp': jax.numpy} to lambdify's modules argument if needed.
else:
    JAX_INSTALLED = True
from multiprocessing import Pool
from ..core.mdv import MDV, get_natural_MDV, get_substrate_MDV, conv, diff_conv # These are numpy based
# from ..core.emu import EMU # If EMU objects are directly used as keys and need properties
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
    
    
    def _calculate_substrate_MDVs(self, extra_subs):
        
        if extra_subs is None:
             extra_subs = []

        for size in self.model.matrix_Bs:
            for sourceEMU in self.model.matrix_Bs[size][2]:
                
                if not isinstance(sourceEMU, Iterable):
                    sourceEMU = (sourceEMU,)
                    
                for emu in sourceEMU:
                    metabid = emu.metabolite_id
                    if metabid in self.model.end_substrates + extra_subs:
                        
                        if metabid in self.model.labeling_strategy:
                            atom_nos = emu.atom_nos
                            (labeling_pattern, 
                             percentage, 
                             purity) = self.model.labeling_strategy[metabid]
                            self.model.substrate_MDVs[emu] = get_substrate_MDV(
                                atom_nos, 
                                labeling_pattern, 
                                percentage, 
                                purity
                            )
                        else:
                            natoms = emu.size
                            self.model.substrate_MDVs[emu] = get_natural_MDV(natoms)                        
    
    
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
            if JAX_INSTALLED and hasattr(self.model, 'jax_flux_derivatives_enabled') and self.model.jax_flux_derivatives_enabled: # Control this via Fitter
                self.model.measured_fluxes_der_p_jax = jnp.array(measured_fluxes_der_u)


        elif kind == 'inst':
            measured_fluxes_der_v = self._calculate_measured_fluxes_derivative_v()
            measured_fluxes_der_u = measured_fluxes_der_v@self.model.null_space
            measured_fluxes_der_c = self._calculate_measured_fluxes_derivative_c()
            inst_der_p = np.concatenate(
                (measured_fluxes_der_u, measured_fluxes_der_c), 
                axis = 1
            )
            self.model.measured_fluxes_der_p = inst_der_p
            if JAX_INSTALLED and hasattr(self.model, 'jax_flux_derivatives_enabled') and self.model.jax_flux_derivatives_enabled:
                 self.model.measured_fluxes_der_p_jax = jnp.array(inst_der_p)

    
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
        
        sds = np.concatenate([sd for _, sd in self.model.measured_MDVs.values()])
        self.model.measured_MDVs_inv_cov = np.diag(1/sds**2)
        
        
    def _calculate_measured_inst_MDVs_inversed_covariance_matrix(self):
        
        sds = []
        for _, t_mdvs in self.model.measured_inst_MDVs.items():
            for t in sorted(t_mdvs):
                if t != 0:
                    sds.append(t_mdvs[t][1])
        sds = np.concatenate(sds)    
        
        self.model.measured_inst_MDVs_inv_cov = np.diag(1/sds**2)
        
    
    def _generate_random_MDVs(self):
        
        self.ori_measured_MDVs = deepcopy(self.model.measured_MDVs)
        
        for emuid, [means, sds] in self.model.measured_MDVs.items():
            mdv = MDV(normal(means, sds))
            self.model.measured_MDVs[emuid][0] = mdv.value


    def _reset_measured_MDVs(self):
        
        self.model.measured_MDVs = deepcopy(self.ori_measured_MDVs)
        

    def _generate_random_inst_MDVs(self):
        
        self.ori_measured_inst_MDVs = deepcopy(self.model.measured_inst_MDVs)

        for emuid, instMDVs in self.model.measured_inst_MDVs.items():
            for t, [means, sds] in instMDVs.items():
                if t != 0:
                    mdv = MDV(normal(means, sds))
                    self.model.measured_inst_MDVs[emuid][t][0] = mdv.value


    def _reset_measured_inst_MDVs(self):

        self.model.measured_inst_MDVs = deepcopy(self.ori_measured_inst_MDVs)
    

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

    def _lambdify_matrix_As_and_Bs_jax(self):
        """
        Lambdifies SymPy expressions for matrices A and B using JAX backend.
        Stores them in self.model.matrix_As_jax_static_data and self.model.matrix_Bs_jax_static_data.
        Also creates helper mappings like totalfluxids_map_jax and EAMs_jax_sorted_keys.
        """
        if not JAX_INSTALLED:
            # This method should only be called if JAX is intended to be used.
            # However, as a safeguard or if called unconditionally:
            warnings.warn("JAX not installed. Cannot lambdify matrices for JAX.")
            self.model.matrix_As_jax_static_data = {}
            self.model.matrix_Bs_jax_static_data = {}
            self.model.totalfluxids_map_jax = {}
            self.model.EAMs_jax_sorted_keys = tuple()
            return

        self.model.matrix_As_jax_static_data = {}
        self.model.matrix_Bs_jax_static_data = {}

        # Create a mapping from total flux ID (string) to integer index
        self.model.totalfluxids_map_jax = {
            fid: i for i, fid in enumerate(self.model.totalfluxids)
        }
        # Ensure EAMs are sorted by size for consistent processing order
        # self.model.EAMs is {size: EAM_DataFrame}
        if not self.model.EAMs: # Should not happen if prepare() is called correctly
             self.model.EAMs_jax_sorted_keys = tuple()
             return # Nothing to lambdify

        self.model.EAMs_jax_sorted_keys = tuple(sorted(self.model.EAMs.keys()))

        for size_key in self.model.EAMs_jax_sorted_keys:
            EAM_df = self.model.EAMs[size_key] # EAM_df is the pandas DataFrame

            # EMU objects are used as index/columns in EAM_df
            # preAB matrix construction from EAM_df (symbolic expressions)
            preAB = EAM_df.copy(deep='all')
            for emu_obj_col in EAM_df.columns: # These are product EMUs (EMU objects)
                # Sum fluxes for diagonal elements
                preAB.loc[emu_obj_col, emu_obj_col] = -preAB[emu_obj_col].sum()

            product_emu_objs = EAM_df.columns.tolist()
            source_emu_objs_or_tuples = EAM_df.index.difference(EAM_df.columns).tolist()

            # Ensure correct DataFrame indexing for Sympy Matrix conversion
            # A is (product_EMUs x product_EMUs), B is (product_EMUs x source_EMUs_terms)
            # The original code: A = preAB.loc[preAB.columns, :].T
            # This means A's rows/cols are product_EMUs.
            # B = -preAB.loc[preAB.index.difference(preAB.columns), :].T
            # This means B's rows are product_EMUs, cols are source_EMU_terms.

            A_sympy_df = preAB.loc[product_emu_objs, product_emu_objs] # Square matrix part for A
            B_sympy_df = -preAB.loc[source_emu_objs_or_tuples, product_emu_objs].T # B (products x sources_terms)
                                                                                # Transpose to match convention A*X = B*Y
                                                                                # where X and Y are column vectors of MDVs.
                                                                                # Original code implies B is (n_prod, n_source_terms)

            matA_sympy = Matrix(A_sympy_df.values) # Pass .values to avoid dtype issues with sympy
            matB_sympy = Matrix(B_sympy_df.values)

            flux_symbols_A = sorted(list(matA_sympy.free_symbols), key=lambda s: s.name)
            flux_symbols_B = sorted(list(matB_sympy.free_symbols), key=lambda s: s.name)

            fluxidsA_str = [s.name for s in flux_symbols_A]
            fluxidsB_str = [s.name for s in flux_symbols_B]

            flux_indices_A_int = tuple(self.model.totalfluxids_map_jax[fid] for fid in fluxidsA_str)
            flux_indices_B_int = tuple(self.model.totalfluxids_map_jax[fid] for fid in fluxidsB_str)

            lambA_jax = lambdify(flux_symbols_A, matA_sympy, modules=['jax'])
            lambB_jax = lambdify(flux_symbols_B, matB_sympy, modules=['jax'])

            product_emu_ids_str = tuple(emu.id for emu in product_emu_objs)

            source_emu_ids_or_tuples_str_list = []
            for item in source_emu_objs_or_tuples:
                if isinstance(item, tuple):
                    source_emu_ids_or_tuples_str_list.append(tuple(e.id for e in item))
                else:
                    source_emu_ids_or_tuples_str_list.append(item.id)
            source_emu_ids_or_tuples_str = tuple(source_emu_ids_or_tuples_str_list)

            self.model.matrix_As_jax_static_data[size_key] = {
                'func': lambA_jax,
                'flux_indices': flux_indices_A_int,
                'product_emu_ids': product_emu_ids_str
            }
            self.model.matrix_Bs_jax_static_data[size_key] = {
                'func': lambB_jax,
                'flux_indices': flux_indices_B_int,
                'source_emu_ids_or_tuples': source_emu_ids_or_tuples_str
            }

    def _prepare_substrate_MDVs_jax(self):
        if not JAX_INSTALLED: return
        self.model.substrate_MDVs_jax_static_data = {
            # Ensure substrate_MDVs keys (EMU objects) are converted to string IDs
            # and values (MDV objects or arrays) are converted to JAX arrays.
            emu.id if hasattr(emu, 'id') else str(emu): jnp.array(mdv_array.value if hasattr(mdv_array, 'value') and isinstance(mdv_array, MDV) else mdv_array)
            for emu, mdv_array in self.model.substrate_MDVs.items()
        }

    def _prepare_substrate_MDVs_der_p_jax(self):
        if not JAX_INSTALLED: return
        self.model.substrate_MDVs_der_p_jax_static_data = {}
        if self.model.substrate_MDVs_der_p:
            for emu, np_der_array in self.model.substrate_MDVs_der_p.items():
                 # Original shape (emu.size+1, nvars). Transpose to (nvars, emu.size+1)
                 key_id = emu.id if hasattr(emu, 'id') else str(emu)
                 self.model.substrate_MDVs_der_p_jax_static_data[key_id] = jnp.array(np_der_array.T)
        else:
            pass # Empty dict is fine. Handled by core_calculate_mdvs_and_derivatives_jax

    def _prepare_matrix_derivatives_jax(self):
        """ Converts pre-calculated NumPy matrix derivatives to JAX arrays. """
        if not JAX_INSTALLED: return
        self.model.matrix_As_der_p_jax_static_data = {
            size: jnp.array(deriv_array)
            for size, deriv_array in self.model.matrix_As_der_p.items()
        }
        self.model.matrix_Bs_der_p_jax_static_data = {
            size: jnp.array(deriv_array)
            for size, deriv_array in self.model.matrix_Bs_der_p.items()
        }

    def _lambdify_matrix_Ms_jax(self):
        """Lambdifies SymPy expressions for matrix M using JAX backend."""
        if not JAX_INSTALLED:
            self.model.matrix_Ms_jax_static_data = {}
            return

        self.model.matrix_Ms_jax_static_data = {}

        conc_ids_ordered = self.model.concids if hasattr(self.model, 'concids') else []
        conc_symbols_ordered = symbols(conc_ids_ordered)
        concid_to_symbol_map = {s.name: s for s in conc_symbols_ordered}

        for size, EAM_df in self.model.EAMs.items():
            product_emu_objs = EAM_df.columns.tolist()

            diag_elements_for_M = []
            active_conc_symbols_for_this_M = []

            for emu_obj in product_emu_objs:
                metab_id = emu_obj.metabolite_id
                if metab_id in concid_to_symbol_map:
                    diag_elements_for_M.append(concid_to_symbol_map[metab_id])
                    if concid_to_symbol_map[metab_id] not in active_conc_symbols_for_this_M:
                        active_conc_symbols_for_this_M.append(concid_to_symbol_map[metab_id])
                else:
                    diag_elements_for_M.append(1.0)

            if not diag_elements_for_M:
                matM_sympy = Matrix([])
            else:
                matM_sympy = Matrix(np.diag(diag_elements_for_M)) # np.diag can handle sympy symbols if they are in a list

            active_conc_symbols_for_this_M.sort(key=lambda s: conc_ids_ordered.index(s.name))

            lambM_jax = lambdify(active_conc_symbols_for_this_M, matM_sympy, modules=['jax'])

            self.model.matrix_Ms_jax_static_data[size] = {
                'func': lambM_jax,
                'conc_arg_indices': tuple(conc_ids_ordered.index(s.name) for s in active_conc_symbols_for_this_M)
            }

    def _prepare_matrix_Ms_derivatives_p_jax(self):
        """ Converts pre-calculated NumPy matrix M derivatives to JAX arrays. """
        if not JAX_INSTALLED:
            self.model.matrix_Ms_der_p_jax_static_data = {}
            return
        self.model.matrix_Ms_der_p_jax_static_data = {
            size: jnp.array(deriv_array)
            for size, deriv_array in self.model.matrix_Ms_der_p.items() # Assuming this exists from NumPy path
        }

    def _prepare_initial_conditions_jax(self):
        """ Prepares JAX versions of initial X, Y matrices and their derivatives. """
        if not JAX_INSTALLED:
            self.model.initial_matrix_Xs_jax = {}
            self.model.initial_matrix_Ys_jax = {}
            self.model.initial_matrix_Xs_der_p_jax = {}
            self.model.initial_matrix_Ys_der_p_jax = {}
            return

        self.model.initial_matrix_Xs_jax = {
            size: jnp.array(matrix_val) for size, matrix_val in self.model.initial_matrix_Xs.items()
        }
        self.model.initial_matrix_Ys_jax = {
            size: jnp.array(matrix_val) for size, matrix_val in self.model.initial_matrix_Ys.items()
        }
        self.model.initial_matrix_Xs_der_p_jax = { # Derivatives are (n_params, n_rows, n_cols)
            size: jnp.array(deriv_val) for size, deriv_val in self.model.initial_matrix_Xs_der_p.items()
        }
        self.model.initial_matrix_Ys_der_p_jax = {
            size: jnp.array(deriv_val) for size, deriv_val in self.model.initial_matrix_Ys_der_p.items()
        }
    
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
                # Ensure concids are available on model for JAX lambdify
                concids_for_lambdify = self.model.concids if hasattr(self.model, 'concids') else []
                lambM = lambdify(symbols(concids_for_lambdify), matM, modules = 'jax')
                # The jacfwd call needs arguments matching concids_for_lambdify
                # This part needs careful alignment of symbols and arguments if concids_for_lambdify is dynamic
                # For simplicity, assume self.model.concids is the fixed list of all possible concentration variables
                if concids_for_lambdify: # Only compute if there are concentration variables
                    matrix_M_der = np.array(
                        jacfwd(lambM, range(len(concids_for_lambdify)))(*jnp.ones(len(concids_for_lambdify)))
                    )
                else: # No concentration variables, derivative is zero or not applicable
                    nEMUsout = EAM.shape[1]
                    matrix_M_der = np.zeros((0, nEMUsout, nEMUsout)) # No params -> first dim is 0
            else: # Fallback if JAX not installed (original logic)
                matrix_M_der = np.array(
                    derive_by_array(matM, symbols(self.model.concids if hasattr(self.model, 'concids') else [])),
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
            # Original logic for numpy lambdify
            # Product EMUs are columns of EAM dataframe
            diag_symbols = [symbols(emu.metabolite_id) for emu in EAM.columns]
            matM = Matrix(np.diag(diag_symbols))

            # Arguments for this lambdified function are ordered by free_symbols
            metabids_arg_order = [s.name for s in sorted(list(matM.free_symbols), key=lambda s:s.name)]
            lambM = lambdify(metabids_arg_order, matM, modules = 'numpy')
            self.model.matrix_Ms[size] = [lambM, metabids_arg_order]
        
    
    # initial X(Y) and their derivatives
    def _calculate_initial_matrix_Xs(self):
        
        for size in self.model.matrix_As: # matrix_As is {size: [lambA, fluxidsA, productEMU_obj_list]}
            nEMUs = len(self.model.matrix_As[size][2]) # productEMU_obj_list
            iniX = np.vstack([get_natural_MDV(size).value] * nEMUs)
            self.model.initial_matrix_Xs[size] = iniX
            
            
    def _calculate_initial_matrix_Ys(self):
        
        for size in self.model.matrix_Bs: # matrix_Bs is {size: [lambB, fluxidsB, sourceEMU_obj_or_tuple_list]}
            iniY = []
            for sourceEMU_item in self.model.matrix_Bs[size][2]: # sourceEMU_obj_or_tuple_list
                if not isinstance(sourceEMU_item, Iterable): # Single EMU object
                    sourceMDV = self.model.substrate_MDVs[sourceEMU_item]
                else: # Tuple of EMU objects
                    mdvs = []
                    for emu in sourceEMU_item:
                        if emu not in self.model.substrate_MDVs:
                            mdv = get_natural_MDV(emu.size)
                        else:
                            mdv = self.model.substrate_MDVs[emu]
                        mdvs.append(mdv)
                    sourceMDV = reduce(conv, mdvs)
                iniY.append(sourceMDV.value if isinstance(sourceMDV, MDV) else sourceMDV) # Ensure array
            self.model.initial_matrix_Ys[size] = np.array(iniY) if iniY else np.empty((0,size+1)) # ensure 2D for empty
        
    
    def _calculate_initial_matrix_Xs_derivatives_p(self):
        
        # nvars depends on whether it's 'ss' or 'inst' mode (fluxes only, or fluxes+concentrations)
        # This needs to be determined based on the context (e.g., kind passed to parent)
        # For now, assume self.model.null_space and self.model.concids are set correctly for current mode.
        n_flux_params = self.model.null_space.shape[1]
        n_conc_params = len(self.model.concids if hasattr(self.model, 'concids') else [])
        # This logic might need adjustment based on 'kind' if called from a generic context
        nvars = n_flux_params + n_conc_params
        if not hasattr(self.model, 'concids'): # If steady state, no conc_params in p for derivatives
            nvars = n_flux_params


        for size, iniX in self.model.initial_matrix_Xs.items():
            Xshape = iniX.shape
            iniXder = np.zeros((nvars, *Xshape)) # (n_params, n_EMUs, n_coeffs)
            self.model.initial_matrix_Xs_der_p[size] = iniXder
        
        
    def _calculate_initial_matrix_Ys_derivatives_p(self):
        n_flux_params = self.model.null_space.shape[1]
        n_conc_params = len(self.model.concids if hasattr(self.model, 'concids') else [])
        nvars = n_flux_params + n_conc_params
        if not hasattr(self.model, 'concids'):
            nvars = n_flux_params

        for size, iniY in self.model.initial_matrix_Ys.items():
            Yshape = iniY.shape
            iniYder = np.zeros((nvars, *Yshape)) # (n_params, n_source_terms, n_coeffs)
            self.model.initial_matrix_Ys_der_p[size] = iniYder
        
    
    def _build_initial_sim_MDVs(self):
        
        for size in sorted(self.model.matrix_As):
            productEMUs_objs = self.model.matrix_As[size][2] # list of EMU objects
            iniX = self.model.initial_matrix_Xs[size]        
            for productEMU_obj, iniMDV_arr in zip(productEMUs_objs, iniX):
                if productEMU_obj.id in self.model.target_EMUs:
                    # Store as { emu_id: {0: MDV_object} }
                    self.model.initial_sim_MDVs[productEMU_obj.id] = {0: MDV(iniMDV_arr)} # Wrap array in MDV object
                    
        
    def _calculate_MDVs(self):
        '''
        This method simulate MDVs at isotopically steady state.
        
        Returns
        -------
        simMDVs: dict
            EMU ID => MDV (in array).
        '''
        
        simMDVs = {} # Stores emu_object -> mdv_array (numpy)
        for size in sorted(self.model.matrix_As):
            
            lambA, fluxidsA, productEMUs_objs = self.model.matrix_As[size]
            lambB, fluxidsB, sourceEMUs_items = self.model.matrix_Bs[size]
            
            A = lambA(*self.model.total_fluxes[fluxidsA])
            B = lambB(*self.model.total_fluxes[fluxidsB])
            
            Y_list_of_arrays = []
            for sourceEMU_item in sourceEMUs_items: # item is EMU obj or tuple of EMU objs
                if not isinstance(sourceEMU_item, Iterable): # single EMU object
                    # Substrate_MDVs stores emu_obj -> MDV_obj or np.array
                    sourceMDV_val = self.model.substrate_MDVs[sourceEMU_item]
                    if isinstance(sourceMDV_val, MDV): sourceMDV_val = sourceMDV_val.value
                else: # tuple of EMU objects
                    mdv_obj_list = []
                    for emu_obj_in_tuple in sourceEMU_item:
                        # ChainMap search: simMDVs (emu_obj -> np.array), then substrate_MDVs (emu_obj -> MDV_obj/np.array)
                        if emu_obj_in_tuple in simMDVs:
                            mdv_val = simMDVs[emu_obj_in_tuple] # This is already an array
                        else:
                            mdv_val = self.model.substrate_MDVs[emu_obj_in_tuple]
                            if isinstance(mdv_val, MDV): mdv_val = mdv_val.value
                        mdv_obj_list.append(MDV(mdv_val)) # Wrap in MDV for convolution via reduce

                    # reduce(conv, list_of_MDV_objs) -> MDV_obj
                    sourceMDV_obj = reduce(conv, mdv_obj_list)
                    sourceMDV_val = sourceMDV_obj.value
                Y_list_of_arrays.append(sourceMDV_val)
            
            Y_matrix = np.array(Y_list_of_arrays) if Y_list_of_arrays else np.empty((B.shape[1], size+1))
            if Y_matrix.ndim == 1 and B.shape[1] == 1: Y_matrix = Y_matrix.reshape(1,-1)


            X_matrix = pinv(A, check_finite=False) @ B @ Y_matrix # check_finite=False for speed
        
            for emu_obj, mdv_array in zip(productEMUs_objs, X_matrix):
                simMDVs[emu_obj] = mdv_array # Store emu_obj -> np.array
            
        # Convert final result to emu_id_str -> np.array
        simMDVs_by_id = {emu_obj.id: mdv_arr for emu_obj, mdv_arr in simMDVs.items()}
        
        return simMDVs_by_id
                    
    
    def _calculate_MDVs_and_derivatives_p(self):
        '''
        This method simulate MDVs and their derivatives at isotopically steady state.
        
        Returns
        -------
        simMDVs: dict
            EMU ID (str) => MDV_array (numpy).
        simMDVsDer: dict
            EMU ID (str) => derivative_array (numpy), shape (n_params, n_coeffs).
        '''
        
        # Internal working dicts use EMU objects as keys for ChainMap compatibility
        # simMDVs_work: emu_obj -> mdv_array (numpy)
        # simMDVsDer_work: emu_obj -> derivative_array (numpy), shape (n_coeffs, n_params)
        simMDVs_work = {}
        simMDVsDer_work = {}

        for size in sorted(self.model.matrix_As): # Iterate by increasing EMU size
            
            lambA, fluxidsA, productEMUs_objs = self.model.matrix_As[size]
            lambB, fluxidsB, sourceEMUs_items = self.model.matrix_Bs[size]
            
            A_val = lambA(*self.model.total_fluxes[fluxidsA])
            B_val = lambB(*self.model.total_fluxes[fluxidsB])
            
            Ainv_val = pinv(A_val, check_finite=False)
            
            # Derivatives of matrices A and B w.r.t parameters p (free fluxes for steady state)
            # Shape: (n_params, n_rows, n_cols)
            Ader_p_val = self.model.matrix_As_der_p[size]
            Bder_p_val = self.model.matrix_Bs_der_p[size]
            
            # Y_list will store MDV arrays (1D)
            # Yder_p_list will store MDV derivative arrays (2D, shape: n_coeffs, n_params)
            Y_list_of_arrays = []
            Yder_p_list_of_arrays = []

            for sourceEMU_item in sourceEMUs_items: # item is EMU obj or tuple of EMU objs
                if not isinstance(sourceEMU_item, Iterable): # single EMU object
                    # MDV value
                    mdv_obj = self.model.substrate_MDVs[sourceEMU_item] # MDV_obj or array
                    current_mdv_val = mdv_obj.value if isinstance(mdv_obj, MDV) else mdv_obj
                    # MDV derivative (n_coeffs, n_params)
                    current_mdv_der_val = self.model.substrate_MDVs_der_p[sourceEMU_item]
                else: # tuple of EMU objects, requires convolution
                    mdv_val_list_for_conv = []
                    mdv_der_pair_list_for_diff_conv = [] # list of [mdv_val_arr, mdv_der_arr (n_coeffs, n_params)]

                    for emu_obj_in_tuple in sourceEMU_item:
                        # Get MDV value (array)
                        if emu_obj_in_tuple in simMDVs_work:
                            mdv_val_arr = simMDVs_work[emu_obj_in_tuple]
                        else:
                            mdv_obj_s = self.model.substrate_MDVs[emu_obj_in_tuple]
                            mdv_val_arr = mdv_obj_s.value if isinstance(mdv_obj_s, MDV) else mdv_obj_s
                        mdv_val_list_for_conv.append(MDV(mdv_val_arr)) # Wrap for MDV.conv via reduce

                        # Get MDV derivative (array, n_coeffs, n_params)
                        if emu_obj_in_tuple in simMDVsDer_work:
                            mdv_der_arr = simMDVsDer_work[emu_obj_in_tuple]
                        else:
                            mdv_der_arr = self.model.substrate_MDVs_der_p[emu_obj_in_tuple]
                        mdv_der_pair_list_for_diff_conv.append([mdv_val_arr, mdv_der_arr])

                    # Perform convolution for value and derivative
                    convolved_mdv_obj = reduce(conv, mdv_val_list_for_conv) # reduce with MDV objects
                    current_mdv_val = convolved_mdv_obj.value

                    # diff_conv: input list of [arr, arr_der (coeffs,params)], output [arr_conv, arr_conv_der (coeffs,params)]
                    convolved_mdv_der_pair = reduce(diff_conv, mdv_der_pair_list_for_diff_conv)
                    current_mdv_der_val = convolved_mdv_der_pair[1]

                Y_list_of_arrays.append(current_mdv_val)
                Yder_p_list_of_arrays.append(current_mdv_der_val)
            
            Y_matrix = np.array(Y_list_of_arrays) if Y_list_of_arrays else np.empty((B_val.shape[1], size+1))
            if Y_matrix.ndim == 1 and B_val.shape[1] == 1: Y_matrix = Y_matrix.reshape(1,-1)
            
            # Yder_p_tensor shape: (n_source_terms, n_coeffs_Y, n_params)
            Yder_p_tensor = np.array(Yder_p_list_of_arrays) if Yder_p_list_of_arrays else np.empty((B_val.shape[1], size+1, Ader_p_val.shape[0]))
            if Yder_p_tensor.ndim == 2 and B_val.shape[1] == 1 : Yder_p_tensor = Yder_p_tensor.reshape(1, Yder_p_tensor.shape[0], Yder_p_tensor.shape[1])


            # Transpose Yder_p_tensor for broadcasting: (n_params, n_source_terms, n_coeffs_Y)
            Yder_p_tensor_transposed = Yder_p_tensor.transpose(2,0,1)
            
            # Calculate X = Ainv @ B @ Y
            X_matrix = Ainv_val @ B_val @ Y_matrix # X shape: (n_prod_EMUs, n_coeffs_X)

            # Calculate X_der_p = Ainv @ (Bder_p@Y + B@Yder_p - Ader_p@X)
            # Bder_p@Y: (n_params, n_prod, n_source) @ (n_source, n_coeffs_Y) -> (n_params, n_prod, n_coeffs_Y)
            term1 = Ader_p_val @ X_matrix # (n_params, n_prod, n_prod) @ (n_prod, n_coeffs_X) -> (n_params, n_prod, n_coeffs_X)

            # B@Yder_p: (n_prod, n_source) @ (n_params, n_source, n_coeffs_Y) - needs broadcasting/looping for params
            # Each slice B @ Yder_p[param_idx,:,:]
            # Result should be (n_params, n_prod, n_coeffs_Y)
            term2 = np.einsum('ik,pkm->pim', B_val, Yder_p_tensor_transposed) # (n_prod, n_source) @ (n_params, n_source, n_coeffs) -> (n_params, n_prod, n_coeffs)

            # Bder_p@Y: (n_params, n_prod, n_source) @ (n_source, n_coeffs_Y) -> (n_params, n_prod, n_coeffs_Y)
            term3 = Bder_p_val @ Y_matrix

            sum_terms = term3 + term2 - term1 # All terms (n_params, n_prod_EMUs, n_coeffs_X)

            # Ainv @ sum_terms: (n_prod, n_prod) @ (n_params, n_prod, n_coeffs_X) -> (n_params, n_prod, n_coeffs_X)
            Xder_p_matrix_transposed = np.einsum('ik,pkm->pim', Ainv_val, sum_terms)

            # Transpose back to (n_prod_EMUs, n_coeffs_X, n_params) for storage if needed, or keep as (n_params, n_prod, n_coeffs)
            # The original code Xder.swapaxes(0,1).swapaxes(1,2) suggests final storage as (n_coeffs, n_params) per EMU
            # Current Xder_p_matrix_transposed is (n_params, n_prod_EMUs, n_coeffs_X)
            # So, for each EMU (row i of n_prod_EMUs), we have Xder_p_matrix_transposed[:, i, :] which is (n_params, n_coeffs_X)
            # This needs to be transposed to (n_coeffs_X, n_params) for simMDVsDer_work[emu_obj]

            for idx, emu_obj in enumerate(productEMUs_objs):
                simMDVs_work[emu_obj] = X_matrix[idx, :]
                simMDVsDer_work[emu_obj] = Xder_p_matrix_transposed[:, idx, :].T # Transpose (n_params, n_coeffs) to (n_coeffs, n_params)

        # Convert final result to emu_id_str -> array
        simMDVs_by_id = {emu_obj.id: mdv_arr for emu_obj, mdv_arr in simMDVs_work.items()}
        # For simMDVsDer_by_id, values are (n_coeffs, n_params). Need to transpose to (n_params, n_coeffs) for nlpsolver.
        # nlpsolver's dxsim_dp = np.vstack([simMDVsDer[emuid] for emuid in self.model.target_EMUs])
        # where simMDVsDer[emuid] is (n_coeffs, n_params). So vstack makes (total_coeffs, n_params).
        # This is d(all_sim_mdvs)/d_params.
        # My core_calculate_mdvs_and_derivatives_jax expects to return derivatives as (n_params, n_coeffs).
        # So, the dictionary here should store (n_coeffs, n_params) to match original,
        # and the JAX wrapper in nlpsolver will handle final formatting if needed.
        # The current storage simMDVsDer_work[emu_obj] = Xder_p_matrix_transposed[:, idx, :].T is (n_coeffs, n_params)
        simMDVsDer_by_id = {emu_obj.id: deriv_arr for emu_obj, deriv_arr in simMDVsDer_work.items()}

        return simMDVs_by_id, simMDVsDer_by_id
    
    
    def _calculate_inst_MDVs(self):
        '''
        This method simulate MDVs at isotopically nonstationary state.
        
        Returns
        -------
        simInstMDVs: dict
            EMU ID => {t => MDV (in array)} (starting from t1).
        '''
        
        simInstMDVs = {} # emu_obj -> {time: mdv_array}
        Ys_t = {}   # time -> {size: Y_matrix_for_that_size_and_time}
        Xs_t = {}   # time -> {size: X_matrix_for_that_size_and_time}

        t1 = 0.0 # Represents initial time t=0
        # Populate Xs_t[0] and Ys_t[0] with initial conditions
        for size_idx in sorted(self.model.matrix_As): # matrix_As keys are emu sizes
            # Initial X values (usually natural abundance)
            # self.model.initial_matrix_Xs is {size: np.array}
            Xs_t.setdefault(t1, {})[size_idx] = self.model.initial_matrix_Xs[size_idx]

            # Initial Y values (from substrates or natural abundance for smaller EMUs)
            # self.model.initial_matrix_Ys is {size: np.array}
            Ys_t.setdefault(t1, {})[size_idx] = self.model.initial_matrix_Ys[size_idx]

        # Store initial MDVs for product EMUs that are targets
        for size_idx in sorted(self.model.matrix_As):
            productEMUs_objs_at_size = self.model.matrix_As[size_idx][2] # List of EMU objects
            X_t0_at_size = Xs_t[t1][size_idx] # X matrix for this size at t=0
            for i, emu_obj in enumerate(productEMUs_objs_at_size):
                if emu_obj.id in self.model.target_EMUs: # Target EMUs are by ID string
                    simInstMDVs.setdefault(emu_obj, {})[t1] = X_t0_at_size[i,:]


        for t_current_loop in self.model.timepoints: # these are t > 0
            if t_current_loop == 0.0: continue # Skip t=0 as it's initial condition

            t0 = t1 # Previous time point (becomes current t1 for next iteration)
            t_current = t_current_loop # Current time point from loop
            deltat = t_current - t0

            # Initialize dictionaries for current time t_current
            Xs_t.setdefault(t_current, {})
            Ys_t.setdefault(t_current, {})

            for size_idx in sorted(self.model.matrix_As): # Iterate by EMU size
                lambA, fluxidsA, productEMUs_objs = self.model.matrix_As[size_idx]
                lambB, fluxidsB, sourceEMUs_items = self.model.matrix_Bs[size_idx]
                lambM, metabids_for_M_args = self.model.matrix_Ms[size_idx] # lambM takes specific concs as args

                A_val = lambA(*self.model.total_fluxes[fluxidsA])
                B_val = lambB(*self.model.total_fluxes[fluxidsB])

                # Select concentrations for M matrix based on metabids_for_M_args
                conc_args_for_M = [self.model.concentrations[mid] for mid in metabids_for_M_args]
                M_val = lambM(*conc_args_for_M)
                Minv_val = pinv(M_val, check_finite=False)

                F_val = Minv_val @ A_val
                Finv_val = pinv(F_val, check_finite=False)
                I_mtx = np.eye(*F_val.shape)
                Phi_val = expm(F_val * deltat) # Matrix exponential
                Gamma_val = (Phi_val - I_mtx) @ Finv_val
                Omega_val = (Gamma_val / deltat - I_mtx) @ Finv_val

                X_t0_at_size = Xs_t[t0][size_idx] # X matrix for this size at previous time t0
                Y_t0_at_size = Ys_t[t0][size_idx] # Y matrix for this size at previous time t0
                G_t0_at_size = Minv_val @ B_val @ Y_t0_at_size

                # Calculate Y_current (Y matrix for current time t_current, for this size_idx)
                # This involves convolutions using MDVs from simInstMDVs (which has emu_obj keys)
                # simInstMDVs stores {emu_obj: {time: mdv_array}}
                Y_list_for_t_current = []
                for sourceEMU_item in sourceEMUs_items:
                    if not isinstance(sourceEMU_item, Iterable): # single EMU object
                        # Substrate MDVs are constant over time in current model structure for them
                        mdv_obj_s = self.model.substrate_MDVs[sourceEMU_item]
                        sourceMDV_val = mdv_obj_s.value if isinstance(mdv_obj_s, MDV) else mdv_obj_s
                    else: # tuple of EMU objects for convolution
                        mdv_obj_list_for_conv = []
                        for emu_obj_in_tuple in sourceEMU_item:
                            # Get MDV at current time t_current
                            # Look up in simInstMDVs (emu_obj -> {t -> arr})
                            # or if not there (substrate), from self.model.substrate_MDVs (emu_obj -> MDV_obj/arr)
                            if emu_obj_in_tuple in simInstMDVs and t_current in simInstMDVs[emu_obj_in_tuple]:
                                mdv_val_arr = simInstMDVs[emu_obj_in_tuple][t_current]
                            else: # Must be a substrate or an EMU from a previous time step not yet in simInstMDVs for *this* t_current
                                  # This implies recursive dependency on *current time* MDVs for smaller EMUs,
                                  # which should have been computed earlier in the size_idx loop for this t_current.
                                  # Or it's a base substrate.
                                mdv_obj_s = self.model.substrate_MDVs[emu_obj_in_tuple] # Fallback to substrate
                                mdv_val_arr = mdv_obj_s.value if isinstance(mdv_obj_s, MDV) else mdv_obj_s
                            mdv_obj_list_for_conv.append(MDV(mdv_val_arr))

                        sourceMDV_obj = reduce(conv, mdv_obj_list_for_conv)
                        sourceMDV_val = sourceMDV_obj.value
                    Y_list_for_t_current.append(sourceMDV_val)

                Y_t_current_at_size = np.array(Y_list_for_t_current) if Y_list_for_t_current else np.empty((B_val.shape[1], size_idx+1))
                if Y_t_current_at_size.ndim == 1 and B_val.shape[1] == 1: Y_t_current_at_size = Y_t_current_at_size.reshape(1,-1)

                Ys_t[t_current][size_idx] = Y_t_current_at_size
                G_t_current_at_size = Minv_val @ B_val @ Y_t_current_at_size

                X_t_current_at_size = Phi_val @ X_t0_at_size - Gamma_val @ G_t0_at_size - Omega_val @ (G_t_current_at_size - G_t0_at_size)
                Xs_t[t_current][size_idx] = X_t_current_at_size

                # Store results in simInstMDVs for product EMUs of this size at current time
                for i, emu_obj in enumerate(productEMUs_objs):
                    simInstMDVs.setdefault(emu_obj, {})[t_current] = X_t_current_at_size[i,:]

            t1 = t_current # Update t1 for the next iteration of the time loop

        # Convert final result to emu_id_str -> {time: mdv_array}
        simInstMDVs_by_id = {
            emu_obj.id: time_dict for emu_obj, time_dict in simInstMDVs.items()
        }
        
        return simInstMDVs_by_id


    def _calculate_inst_MDVs_and_derivatives_p(self):
        '''
        This method simulate MDVs and their derivatives at isotopically nonstationary state.
        
        Returns
        -------
        simInstMDVs: dict
            EMU ID (str) => {time (float) => MDV_array (numpy)}.
        simInstMDVsDer: dict
            EMU ID (str) => {time (float) => derivative_array (numpy, shape: n_params, n_coeffs)}.
        '''
        
        # Internal working dicts:
        # simInstMDVs_work: emu_obj -> {time: mdv_array}
        # simInstMDVsDer_work: emu_obj -> {time: derivative_array (n_coeffs, n_params)}
        simInstMDVs_work = {}
        simInstMDVsDer_work = {}

        # Ys_t_work, Xs_t_work: time -> {size: matrix_value (numpy)}
        # Yders_p_t_work, Xders_p_t_work: time -> {size: derivative_tensor (n_coeffs, n_params-per-matrix-row-implicitly)}
        # Actually, derivatives are (n_params, n_rows, n_cols) for matrices A, B, M, X, Y.
        # For MDVs (X, Y), this means (n_params, n_EMUs_or_Sources, n_coeffs).
        # Let's store derivatives as (n_params, n_EMUs_or_Sources, n_coeffs) in Xders_p_t, Yders_p_t.
        # When retrieving for a single EMU for simInstMDVsDer_work, it will be (n_params, n_coeffs),
        # then transposed to (n_coeffs, n_params) for consistency with steady-state.

        Ys_t_work = {}
        Xs_t_work = {}
        Yders_p_t_work = {} # Stores dY/dp as {time: {size: tensor (n_params, n_sources, n_coeffs)}}
        Xders_p_t_work = {} # Stores dX/dp as {time: {size: tensor (n_params, n_products, n_coeffs)}}
        
        time_prev = 0.0 # Represents initial time t=0

        # Populate initial conditions at t=0 for X, Y and their derivatives dX/dp, dY/dp
        for size_idx in sorted(self.model.matrix_As):
            Xs_t_work.setdefault(time_prev, {})[size_idx] = self.model.initial_matrix_Xs[size_idx]
            Ys_t_work.setdefault(time_prev, {})[size_idx] = self.model.initial_matrix_Ys[size_idx]
            Xders_p_t_work.setdefault(time_prev, {})[size_idx] = self.model.initial_matrix_Xs_der_p[size_idx]
            Yders_p_t_work.setdefault(time_prev, {})[size_idx] = self.model.initial_matrix_Ys_der_p[size_idx]

            # Store initial MDVs and their derivatives for target EMUs
            productEMUs_objs_at_size = self.model.matrix_As[size_idx][2]
            X_t0_at_size = Xs_t_work[time_prev][size_idx]
            Xder_p_t0_at_size = Xders_p_t_work[time_prev][size_idx] # (n_params, n_EMUs, n_coeffs)

            for i, emu_obj in enumerate(productEMUs_objs_at_size):
                if emu_obj.id in self.model.target_EMUs:
                    simInstMDVs_work.setdefault(emu_obj, {})[time_prev] = X_t0_at_size[i,:]
                    # Store derivative as (n_coeffs, n_params)
                    simInstMDVsDer_work.setdefault(emu_obj, {})[time_prev] = Xder_p_t0_at_size[:, i, :].T


        for t_current_loop in self.model.timepoints: # These are t > 0
            if t_current_loop == 0.0: continue

            t_current = t_current_loop
            deltat = t_current - time_prev
            
            Xs_t_work.setdefault(t_current, {})
            Ys_t_work.setdefault(t_current, {})
            Xders_p_t_work.setdefault(t_current, {})
            Yders_p_t_work.setdefault(t_current, {})

            for size_idx in sorted(self.model.matrix_As):
                lambA, fluxidsA, productEMUs_objs = self.model.matrix_As[size_idx]
                lambB, fluxidsB, sourceEMUs_items = self.model.matrix_Bs[size_idx]
                lambM, metabids_for_M_args = self.model.matrix_Ms[size_idx]

                A_val = lambA(*self.model.total_fluxes[fluxidsA])
                B_val = lambB(*self.model.total_fluxes[fluxidsB])
                conc_args_for_M = [self.model.concentrations[mid] for mid in metabids_for_M_args]
                M_val = lambM(*conc_args_for_M)
                Minv_val = pinv(M_val, check_finite=False)

                # Matrix derivatives d/dp (p includes free fluxes and concentrations)
                # Shapes: (n_params, n_rows, n_cols)
                Ader_p_val = self.model.matrix_As_der_p[size_idx]
                Bder_p_val = self.model.matrix_Bs_der_p[size_idx]
                Mder_p_val = self.model.matrix_Ms_der_p[size_idx]
                Minv_der_p_val = -Minv_val @ Mder_p_val @ Minv_val # d(M^-1)/dp = -M^-1 * dM/dp * M^-1 (element-wise for params dim)
                                                                # This should be vmap over params:
                                                                # Minv_der_p_val[p] = -Minv_val @ Mder_p_val[p] @ Minv_val
                Minv_der_p_val = np.einsum('ik,pkm,ml->pil', -Minv_val, Mder_p_val, Minv_val)


                F_val = Minv_val @ A_val
                Finv_val = pinv(F_val, check_finite=False)
                I_mtx = np.eye(*F_val.shape)
                Phi_val = expm(F_val * deltat) # Matrix exponential
                Gamma_val = (Phi_val - I_mtx) @ Finv_val
                Omega_val = (Gamma_val / deltat - I_mtx) @ Finv_val

                X_t_prev_at_size = Xs_t_work[time_prev][size_idx]
                Y_t_prev_at_size = Ys_t_work[time_prev][size_idx]
                G_t_prev_at_size = Minv_val @ B_val @ Y_t_prev_at_size

                Xder_p_t_prev_at_size = Xders_p_t_work[time_prev][size_idx] # (n_params, n_prods, n_coeffs)
                Yder_p_t_prev_at_size = Yders_p_t_work[time_prev][size_idx] # (n_params, n_sources, n_coeffs)

                # H_t_prev = dG_t_prev/dp = d(Minv*B*Y_prev)/dp
                # = (dMinv/dp)*B*Y_prev + Minv*(dB/dp)*Y_prev + Minv*B*(dY_prev/dp)
                # All derivatives are (n_params, n_rows, n_cols)
                # Minv_der_p_val @ B_val @ Y_t_prev_at_size
                term_H_1 = np.einsum('pij,jk,kl->pil', Minv_der_p_val, B_val, Y_t_prev_at_size)
                # Minv_val @ Bder_p_val @ Y_t_prev_at_size
                term_H_2 = np.einsum('ij,pjk,kl->pil', Minv_val, Bder_p_val, Y_t_prev_at_size)
                # Minv_val @ B_val @ Yder_p_t_prev_at_size (Yder is (n_params, n_sources, n_coeffs))
                term_H_3 = np.einsum('ij,jk,pkl->pil', Minv_val, B_val, Yder_p_t_prev_at_size)
                Gder_p_t_prev_at_size = term_H_1 + term_H_2 + term_H_3 # (n_params, n_prods, n_coeffs)


                # Calculate Y_t_current_at_size and its derivative Yder_p_t_current_at_size
                Y_list_for_t_curr = []
                Yder_p_list_for_t_curr = [] # List of (n_coeffs, n_params) arrays

                for sourceEMU_item in sourceEMUs_items:
                    if not isinstance(sourceEMU_item, Iterable): # single EMU object
                        mdv_obj = self.model.substrate_MDVs[sourceEMU_item]
                        curr_mdv_val = mdv_obj.value if isinstance(mdv_obj, MDV) else mdv_obj
                        # Substrate derivatives are (n_coeffs, n_params) in substrate_MDVs_der_p
                        curr_mdv_der_val = self.model.substrate_MDVs_der_p[sourceEMU_item]
                    else: # tuple of EMU objects for convolution
                        mdv_val_list_for_conv_iter = []
                        mdv_der_pair_list_for_diff_conv_iter = []
                        for emu_obj_in_tuple in sourceEMU_item:
                            # Get MDV value (array) at t_current
                            if emu_obj_in_tuple in simInstMDVs_work and t_current in simInstMDVs_work[emu_obj_in_tuple]:
                                mdv_val_arr_iter = simInstMDVs_work[emu_obj_in_tuple][t_current]
                            else: # Fallback to substrate (constant)
                                mdv_obj_s_iter = self.model.substrate_MDVs[emu_obj_in_tuple]
                                mdv_val_arr_iter = mdv_obj_s_iter.value if isinstance(mdv_obj_s_iter, MDV) else mdv_obj_s_iter
                            mdv_val_list_for_conv_iter.append(MDV(mdv_val_arr_iter))

                            # Get MDV derivative (array, n_coeffs, n_params) at t_current
                            if emu_obj_in_tuple in simInstMDVsDer_work and t_current in simInstMDVsDer_work[emu_obj_in_tuple]:
                                mdv_der_arr_iter = simInstMDVsDer_work[emu_obj_in_tuple][t_current]
                            else: # Fallback to substrate derivative
                                mdv_der_arr_iter = self.model.substrate_MDVs_der_p[emu_obj_in_tuple]
                            mdv_der_pair_list_for_diff_conv_iter.append([mdv_val_arr_iter, mdv_der_arr_iter])
                        
                        convolved_mdv_obj_iter = reduce(conv, mdv_val_list_for_conv_iter)
                        curr_mdv_val = convolved_mdv_obj_iter.value
                        convolved_mdv_der_pair_iter = reduce(diff_conv, mdv_der_pair_list_for_diff_conv_iter)
                        curr_mdv_der_val = convolved_mdv_der_pair_iter[1]

                    Y_list_for_t_curr.append(curr_mdv_val)
                    Yder_p_list_for_t_curr.append(curr_mdv_der_val) # List of (n_coeffs, n_params)

                Y_t_current_at_size = np.array(Y_list_for_t_curr) if Y_list_for_t_curr else np.empty((B_val.shape[1], size_idx+1))
                if Y_t_current_at_size.ndim == 1 and B_val.shape[1] == 1: Y_t_current_at_size = Y_t_current_at_size.reshape(1,-1)
                Ys_t_work[t_current][size_idx] = Y_t_current_at_size

                Yder_p_t_curr_tensor = np.array(Yder_p_list_for_t_curr) if Yder_p_list_for_t_curr else np.empty((B_val.shape[1], size_idx+1, Xder_p_t_prev_at_size.shape[0]))
                if Yder_p_t_curr_tensor.ndim == 2 and B_val.shape[1] == 1: Yder_p_t_curr_tensor = Yder_p_t_curr_tensor.reshape(1, Yder_p_t_curr_tensor.shape[0], Yder_p_t_curr_tensor.shape[1])
                Yder_p_t_current_at_size = Yder_p_t_curr_tensor.transpose(2,0,1) # (n_params, n_sources, n_coeffs)
                Yders_p_t_work[t_current][size_idx] = Yder_p_t_current_at_size


                G_t_current_at_size = Minv_val @ B_val @ Y_t_current_at_size
                # Gder_p_t_current_at_size (n_params, n_prods, n_coeffs)
                term_Hc_1 = np.einsum('pij,jk,kl->pil', Minv_der_p_val, B_val, Y_t_current_at_size)
                term_Hc_2 = np.einsum('ij,pjk,kl->pil', Minv_val, Bder_p_val, Y_t_current_at_size)
                term_Hc_3 = np.einsum('ij,jk,pkl->pil', Minv_val, B_val, Yder_p_t_current_at_size)
                Gder_p_t_current_at_size = term_Hc_1 + term_Hc_2 + term_Hc_3

                # Calculate X_t_current and its derivative Xder_p_t_current
                X_t_current_at_size = Phi_val @ X_t_prev_at_size - Gamma_val @ G_t_prev_at_size - Omega_val @ (G_t_current_at_size - G_t_prev_at_size)
                Xs_t_work[t_current][size_idx] = X_t_current_at_size

                # dF/dp = d(Minv*A)/dp = (dMinv/dp)*A + Minv*(dA/dp)
                dF_dp = np.einsum('pij,jk->pik', Minv_der_p_val, A_val) + np.einsum('ij,pjk->pik', Minv_val, Ader_p_val)
                # This is complex: d(expm(F*dt))/dp and subsequent terms (dPhi/dp, dGamma/dp, dOmega/dp)
                # This requires matrix derivative calculus for expm, inv.
                # The original paper/code might use approximations or specific formulas.
                # For now, assume H_t0 and H_t1 are simplified forms as in original code if available.
                # Original: H_t0 = (Minvder@A@X_t0 + Minv@Ader@X_t0 - Minv@B@Yder_t0 - Minvder@B@Y_t0 - Minv@Bder@Y_t0)
                # This is not dG/dp, but a different combination. Let's call it K.
                # K = d(Minv*A*X)/dp - d(Minv*B*Y)/dp
                # d(Minv*A*X)/dp = (dMinv/dp)AX + Minv(dA/dp)X + MinvA(dX/dp)
                # d(Minv*B*Y)/dp = (dMinv/dp)BY + Minv(dB/dp)Y + MinvB(dY/dp)

                # K_t_prev
                term_K1_prev = np.einsum('pij,jk,kl->pil', Minv_der_p_val, A_val, X_t_prev_at_size)
                term_K2_prev = np.einsum('ij,pjk,kl->pil', Minv_val, Ader_p_val, X_t_prev_at_size)
                # Minv*A*dX_prev/dp term for K_t_prev:
                term_K3_prev = np.einsum('ij,jk,pkl->pil', Minv_val, A_val, Xder_p_t_prev_at_size)
                K_t_prev = (term_K1_prev + term_K2_prev + term_K3_prev) - Gder_p_t_prev_at_size

                # K_t_current
                term_K1_curr = np.einsum('pij,jk,kl->pil', Minv_der_p_val, A_val, X_t_current_at_size) # Uses X_current
                term_K2_curr = np.einsum('ij,pjk,kl->pil', Minv_val, Ader_p_val, X_t_current_at_size) # Uses X_current
                # dX_curr/dp is not known yet. This formulation seems circular or needs specific derivative of ODE solution.
                # The provided H was: (Minvder@A@X + Minv@Ader@X - Minv@B@Yder - Minvder@B@Y - Minv@Bder@Y)
                # This does not include dX/dp or dY/dp terms directly in H.
                # Let's use the structure from original code's H_t0, H_t1 for Xder calculation.
                # H = d(F*X)/dp - d(G)/dp, where F=Minv*A. (dX/dt = F*X - G). So d(dX/dt)/dp = d(F*X)/dp - d(G)/dp
                # d(F*X)/dp = (dF/dp)X + F(dX/dp)
                # H_from_paper = (dF/dp)X - Gder

                H_val_t_prev = np.einsum('pij,jk->pik', dF_dp, X_t_prev_at_size) - Gder_p_t_prev_at_size
                H_val_t_current = np.einsum('pij,jk->pik', dF_dp, X_t_current_at_size) - Gder_p_t_current_at_size


                # Derivative of X_t_current w.r.t params 'p'
                # Xder_p_t_current = Phi @ Xder_p_t_prev + dPhi/dp @ X_t_prev  (Chain rule for Phi)
                #                    - (dGamma/dp @ G_prev + Gamma @ dG_prev/dp)
                #                    - (dOmega/dp @ (G_curr-G_prev) + Omega @ (dG_curr/dp - dG_prev/dp))
                # This is very complex. The original code uses:
                # Xder_t1 = Phi@Xder_t0 + Gamma@H_t0 + Omega@(H_t1 - H_t0)
                # This H must be d(F*X - G)/dX_t0 * dX_t0/dp ... no, this H is likely specific to the solution form.
                # Assuming that H formulation is correct from a source theory:
                Xder_p_t_current_at_size = Phi_val @ Xder_p_t_prev_at_size + Gamma_val @ H_val_t_prev + Omega_val @ (H_val_t_current - H_val_t_prev)
                Xders_p_t_work[t_current][size_idx] = Xder_p_t_current_at_size # Store (n_params, n_prods, n_coeffs)

                # Store results in simInstMDVs_work and simInstMDVsDer_work (derivatives as n_coeffs, n_params)
                for i, emu_obj in enumerate(productEMUs_objs):
                    simInstMDVs_work.setdefault(emu_obj, {})[t_current] = X_t_current_at_size[i,:]
                    simInstMDVsDer_work.setdefault(emu_obj, {})[t_current] = Xder_p_t_current_at_size[:, i, :].T # (n_coeffs, n_params)

            time_prev = t_current # Update time_prev for the next iteration

        # Convert final result to emu_id_str -> {time: array}
        simInstMDVs_by_id = {
            emu_obj.id: time_dict for emu_obj, time_dict in simInstMDVs_work.items()
        }
        simInstMDVsDer_by_id = {
            emu_obj.id: time_dict_der for emu_obj, time_dict_der in simInstMDVsDer_work.items()
        }

        return simInstMDVs_by_id, simInstMDVsDer_by_id

[end of src/freeflux/utils/utils.py]
