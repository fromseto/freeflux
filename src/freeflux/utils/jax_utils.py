"""JAX utility functions for freeflux."""

import jax
import jax.numpy as jnp
from jax.scipy.linalg import pinv


def jax_conv(arr1, arr2):
    """JAX equivalent of polynomial convolution for 1D arrays (like numpy.convolve)."""
    if arr1 is None or arr2 is None:
        # This case should ideally be handled by ensuring valid (non-None) arrays are passed,
        # or by defining specific behavior (e.g., if one is identity MDV [1.0]).
        raise ValueError("jax_conv received None input. This should be handled before calling.")
    return jnp.convolve(arr1, arr2)

def jax_diff_conv(mdv1_mdv1der_pair, mdv2_mdv2der_pair, num_params_for_zeros=None):
    """
    Calculates (conv(mdv1, mdv2), d(conv(mdv1, mdv2))/dp) using JAX.
    mdv1, mdv2 are 1D JAX arrays.
    mdv1der, mdv2der are 2D JAX arrays (n_params, n_coeffs_for_mdv), or None.
    num_params_for_zeros: integer, required if both mdv1der and mdv2der are None,
                          to correctly shape the zero derivative.

    Returns: (convolved_mdv [1D], convolved_mdv_derivative [2D: (n_params, n_coeffs_conv)])
             The derivative part can be None if both input derivatives are None and num_params_for_zeros is not given.
    """
    mdv1, mdv1der = mdv1_mdv1der_pair
    mdv2, mdv2der = mdv2_mdv2der_pair

    # Ensure mdv1 and mdv2 are not None
    if mdv1 is None or mdv2 is None:
        raise ValueError("MDV inputs to jax_diff_conv cannot be None.")

    convolved_mdv = jnp.convolve(mdv1, mdv2)

    term1_der = None
    if mdv1der is not None:
        # mdv1der shape: (n_params, len(mdv1))
        # mdv2 shape: (len(mdv2),)
        term1_der = jax.vmap(lambda d1_row: jnp.convolve(d1_row, mdv2))(mdv1der)
        # term1_der shape: (n_params, len(convolved_mdv))

    term2_der = None
    if mdv2der is not None:
        # mdv1 shape: (len(mdv1),)
        # mdv2der shape: (n_params, len(mdv2))
        term2_der = jax.vmap(lambda d2_row: jnp.convolve(mdv1, d2_row))(mdv2der)
        # term2_der shape: (n_params, len(convolved_mdv))

    if term1_der is not None and term2_der is not None:
        convolved_mdv_der = term1_der + term2_der
    elif term1_der is not None:
        convolved_mdv_der = term1_der
    elif term2_der is not None:
        convolved_mdv_der = term2_der
    else:
        # Both derivatives are None (conceptually zero).
        if num_params_for_zeros is not None:
            convolved_mdv_der = jnp.zeros((num_params_for_zeros, len(convolved_mdv)))
        else:
            # This case should ideally be avoided by providing num_params_for_zeros
            # if there's a possibility of all derivatives being None.
            # For JIT, shapes must be consistent.
            raise ValueError("num_params_for_zeros must be provided if all input derivatives can be None.")


    return convolved_mdv, convolved_mdv_der


def core_calculate_mdvs_jax(
    total_fluxes_jax,
    matrix_As_static_data,
    matrix_Bs_static_data,
    substrate_MDVs_jax,
    target_EMU_ids_tuple, # Currently unused in core calc, filtering is done outside
    sorted_emu_sizes_tuple
    ):
    """
    Calculates simulated Mass Distribution Vectors (MDVs) using JAX.
    This function is designed to be JIT-compatible.
    Static data (matrices, EMU lists) are passed in pytrees.
    """
    sim_MDVs_dict = {}

    for size in sorted_emu_sizes_tuple:
        A_data = matrix_As_static_data[size]
        B_data = matrix_Bs_static_data[size]

        lambA_jax = A_data['func']
        flux_indices_A = A_data['flux_indices']
        product_EMU_ids = A_data['product_emu_ids']

        lambB_jax = B_data['func']
        flux_indices_B = B_data['flux_indices']
        source_EMU_ids_or_tuples = B_data['source_emu_ids_or_tuples']

        fluxes_for_A = total_fluxes_jax.take(jnp.array(flux_indices_A))
        fluxes_for_B = total_fluxes_jax.take(jnp.array(flux_indices_B))

        A = lambA_jax(*fluxes_for_A)
        B = lambB_jax(*fluxes_for_B)

        Y_parts = []
        for source_item in source_EMU_ids_or_tuples:
            if isinstance(source_item, str):
                source_emu_id = source_item
                mdv = sim_MDVs_dict.get(source_emu_id, substrate_MDVs_jax.get(source_emu_id))
                if mdv is None: raise ValueError(f"MDV not found for {source_emu_id} in size {size}")
                Y_parts.append(mdv)
            else:
                mdvs_to_convolve = []
                for emu_id_in_tuple in source_item:
                    mdv = sim_MDVs_dict.get(emu_id_in_tuple, substrate_MDVs_jax.get(emu_id_in_tuple))
                    if mdv is None: raise ValueError(f"MDV not found for {emu_id_in_tuple} in convolution for size {size}")
                    mdvs_to_convolve.append(mdv)

                if mdvs_to_convolve:
                    current_conv = mdvs_to_convolve[0]
                    for i in range(1, len(mdvs_to_convolve)):
                        current_conv = jax_conv(current_conv, mdvs_to_convolve[i])
                    Y_parts.append(current_conv)
                else: # Should not happen with valid model structure
                    raise ValueError(f"Empty convolution list for size {size}")


        num_source_terms = B.shape[1]
        # mdv_len_this_size = size + 1 # This was an assumption about Y matrix content.
                                     # Each row in Y is an MDV. These MDVs must have a length compatible with B's columns.
                                     # The EMU formulation implies that B projects/combines these source MDVs
                                     # into contributions for product EMUs of current `size`.
                                     # The X = pinv(A)@B@Y implies Y's rows are MDVs that B can operate on.
                                     # The resulting X will have rows of length `size+1`.

        if Y_parts:
            # Y must be a 2D array (matrix) for B @ Y.
            # Each element of Y_parts is a 1D MDV array.
            # Their lengths can vary if they come from EMUs of different sizes.
            # This is a CRITICAL POINT: The original code `Y = np.array(Y)` implicitly assumes all MDVs in Y_parts
            # have the same length to form a 2D array. This is true if all source EMUs (or convolutions thereof)
            # for a given product size `s` also result in MDVs corresponding to size `s`.
            # This seems to be an implicit assumption of the (A*X = B*Y) formulation per size.
            # Let's assume all mdvs in Y_parts for a given `size` effectively have length `size+1`.
            mdv_len_check = size + 1
            Y = jnp.array(Y_parts)
            if Y.ndim == 1 and num_source_terms == 1: # Single source term, Y_parts contained one 1D array
                 Y = Y.reshape(1, -1)

            if Y.shape[0] != num_source_terms:
                raise ValueError(f"Mismatch in Y parts ({Y.shape[0]}) and B matrix columns ({num_source_terms}) for size {size}")
            if Y.shape[1] != mdv_len_check:
                 # This is a deviation from the simple assumption.
                 # If this happens, B must be structured to handle it (e.g. padded EMUs).
                 # For now, stick to the assumption that Y rows are all length `size+1`.
                 raise ValueError(f"Y row length {Y.shape[1]} does not match expected {mdv_len_check} for size {size}")

        elif num_source_terms > 0:
            raise ValueError(f"Y_parts is empty but B expects {num_source_terms} inputs for size {size}")
        else: # num_source_terms == 0, Y_parts is empty
            # If B has 0 columns, Y should be (0, mdv_len_for_products_of_this_size)
            Y = jnp.zeros((0, size + 1))

        X = pinv(A) @ B @ Y

        for i, emu_id in enumerate(product_EMU_ids):
            sim_MDVs_dict[emu_id] = X[i, :]

    return sim_MDVs_dict


def core_calculate_mdvs_and_derivatives_jax(
    total_fluxes_jax,
    matrix_As_static_data,
    matrix_Bs_static_data,
    substrate_MDVs_jax,
    matrix_As_der_p_static_data,
    matrix_Bs_der_p_static_data,
    substrate_MDVs_der_p_jax,
    target_EMU_ids_tuple, # Currently unused
    sorted_emu_sizes_tuple,
    num_free_fluxes
    ):
    """
    Calculates simulated MDVs and their derivatives w.r.t. parameters (free fluxes `u`) using JAX.
    """
    sim_MDVs_dict = {}
    sim_MDVs_der_dict = {}

    for size in sorted_emu_sizes_tuple:
        A_data = matrix_As_static_data[size]
        B_data = matrix_Bs_static_data[size]

        lambA_jax = A_data['func']
        flux_indices_A = A_data['flux_indices']
        product_EMU_ids = A_data['product_emu_ids']

        lambB_jax = B_data['func']
        flux_indices_B = B_data['flux_indices']
        source_EMU_ids_or_tuples = B_data['source_emu_ids_or_tuples']

        fluxes_for_A = total_fluxes_jax.take(jnp.array(flux_indices_A))
        fluxes_for_B = total_fluxes_jax.take(jnp.array(flux_indices_B))

        A = lambA_jax(*fluxes_for_A)
        B = lambB_jax(*fluxes_for_B)
        Ainv = pinv(A)

        Ader_p = matrix_As_der_p_static_data[size]
        Bder_p = matrix_Bs_der_p_static_data[size]

        Y_parts = []
        Yder_p_parts = []

        mdv_len_for_X = size + 1 # Product EMUs of this size have this MDV length

        for source_item in source_EMU_ids_or_tuples:
            if isinstance(source_item, str):
                source_emu_id = source_item
                mdv = sim_MDVs_dict.get(source_emu_id, substrate_MDVs_jax.get(source_emu_id))
                if mdv is None: raise ValueError(f"MDV not found for {source_emu_id} in size {size} (derivative calc)")

                mdv_der_p = sim_MDVs_der_dict.get(source_emu_id, substrate_MDVs_der_p_jax.get(source_emu_id))
                if mdv_der_p is None:
                    mdv_der_p = jnp.zeros((num_free_fluxes, len(mdv)))

                Y_parts.append(mdv)
                Yder_p_parts.append(mdv_der_p)
            else:
                mdv_val_list_for_conv = []
                mdv_der_pair_list_for_diff_conv = []

                for emu_id_in_tuple in source_item:
                    mdv_val = sim_MDVs_dict.get(emu_id_in_tuple, substrate_MDVs_jax.get(emu_id_in_tuple))
                    if mdv_val is None: raise ValueError(f"MDV not found for {emu_id_in_tuple} in convolution for size {size} (derivative calc)")

                    mdv_der_val = sim_MDVs_der_dict.get(emu_id_in_tuple, substrate_MDVs_der_p_jax.get(emu_id_in_tuple))
                    if mdv_der_val is None:
                        mdv_der_val = jnp.zeros((num_free_fluxes, len(mdv_val)))

                    mdv_val_list_for_conv.append(mdv_val)
                    mdv_der_pair_list_for_diff_conv.append([mdv_val, mdv_der_val])

                if mdv_val_list_for_conv:
                    current_conv_mdv = mdv_val_list_for_conv[0]
                    for i in range(1, len(mdv_val_list_for_conv)):
                        current_conv_mdv = jax_conv(current_conv_mdv, mdv_val_list_for_conv[i])
                    Y_parts.append(current_conv_mdv)

                    current_conv_mdv_and_der_pair = mdv_der_pair_list_for_diff_conv[0]
                    for i in range(1, len(mdv_der_pair_list_for_diff_conv)):
                        next_pair = mdv_der_pair_list_for_diff_conv[i]
                        current_conv_mdv_and_der_pair = jax_diff_conv(current_conv_mdv_and_der_pair, next_pair, num_params_for_zeros=num_free_fluxes)
                    Yder_p_parts.append(current_conv_mdv_and_der_pair[1])
                else:
                     raise ValueError(f"Empty convolution list for size {size} (derivative calc)")

        num_source_terms = B.shape[1]
        # Assumption: All MDVs in Y_parts (after convolution) correspond to the length expected by B for this size iteration.
        # This length is `size+1`.
        expected_Y_row_len = size + 1

        if Y_parts:
            Y = jnp.array(Y_parts)
            if Y.ndim == 1 and num_source_terms == 1 : Y = Y.reshape(1,-1) # Ensure Y is 2D

            Yder_p_temp = jnp.array(Yder_p_parts) # (n_source_terms, n_params, current_mdv_len)
            if Yder_p_temp.ndim == 2 and num_source_terms == 1: Yder_p_temp = Yder_p_temp.reshape(1, num_free_fluxes, -1) # Ensure 3D

            Yder_p = jnp.transpose(Yder_p_temp, (1,0,2)) # (n_params, n_source_terms, current_mdv_len)

            if Y.shape[0] != num_source_terms:
                raise ValueError(f"Y shape {Y.shape} mismatch B columns {num_source_terms} for size {size}")
            if Y.shape[1] != expected_Y_row_len:
                raise ValueError(f"Y row length {Y.shape[1]} does not match expected {expected_Y_row_len} for size {size}")
            if Yder_p.shape[1] != num_source_terms or Yder_p.shape[2] != expected_Y_row_len:
                 raise ValueError(f"Yder_p shape {Yder_p.shape} mismatch B columns {num_source_terms} or MDV length {expected_Y_row_len} for size {size}")

        elif num_source_terms > 0:
             raise ValueError(f"Y_parts is empty but B expects {num_source_terms} inputs for size {size} (derivative calc)")
        else: # num_source_terms == 0
            Y = jnp.zeros((0, expected_Y_row_len))
            Yder_p = jnp.zeros((num_free_fluxes, 0, expected_Y_row_len))


        X = Ainv @ B @ Y  # X shape: (n_prod_emus, mdv_len_for_X)

        term1 = jax.vmap(lambda Bder_slice: Bder_slice @ Y)(Bder_p)
        term2 = jax.vmap(lambda Yder_slice: B @ Yder_slice)(Yder_p)
        term3 = jax.vmap(lambda Ader_slice: Ader_slice @ X)(Ader_p)
        sum_terms = term1 + term2 - term3
        Xder_p = jax.vmap(lambda sum_slice: Ainv @ sum_slice)(sum_terms) # Xder_p shape: (n_params, n_prod_emus, mdv_len_for_X)

        for i, emu_id in enumerate(product_EMU_ids):
            sim_MDVs_dict[emu_id] = X[i, :]
            sim_MDVs_der_dict[emu_id] = Xder_p[:, i, :] # Store as (n_params, mdv_len_for_X)

    return sim_MDVs_dict, sim_MDVs_der_dict
