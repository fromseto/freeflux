'''Define data input functions.'''


__author__ = 'Chao Wu'
__date__ = '05/23/2022' # Updated date can be considered later if project policy dictates


from os.path import splitext, isfile
import re
import pandas as pd


def read_model_from_file(file):
    '''
    Parameters
    ----------
    file: file path
        tsv or excel file with reactions.
    '''
    
    ext = splitext(file)[1]
    if re.search(r'tsv', ext):
        data = pd.read_csv(
            file, 
            sep = '\t', 
            comment = '#', 
            header = None, 
            names = ['subs', 'pros', 'rev'], 
            index_col = 0
        )
    elif re.search(r'xls', ext):
        data = pd.read_excel(
            file, 
            comment = '#', 
            header = None, 
            names = ['subs', 'pros', 'rev'], 
            index_col = 0
        )
    else:
        raise TypeError('can only read from .tsv or excel file')
    
    data = data.dropna() # Original had data.dropna(), assuming it's for rows with all NaN
    
    return data


def read_preset_values_from_file(file):
    '''
    Parameters
    ----------
    file: file path
        tsv or excel file.
    '''
    
    ext = splitext(file)[1]
    if re.search(r'tsv', ext):
        data = pd.read_csv(
            file, 
            sep = '\t', 
            comment = '#', 
            header = None, 
            names = ['value'], 
            index_col = 0, 
        ).squeeze("columns") # pandas > 1.0, squeeze() needs axis/columns
    elif re.search(r'xls', ext):
        data = pd.read_excel(
            file, 
            comment = '#', 
            header = None, 
            names = ['value'], 
            index_col = 0, 
        ).squeeze("columns") # pandas > 1.0, squeeze() needs axis/columns
    else:
        raise TypeError('can only read from .tsv or excel file')
    
    data = data.dropna()

    return data
        

def read_initial_values(ini, ids):
    '''
    Parameters
    ----------
    ini: ser of file in .tsv or .xlsx
        Initial values.
    ids: list
        IDs of fluxes or concentrations in correct order.
    '''

    if isinstance(ini, str) and isfile(ini): # Check if string and is a file
        ini_values = read_preset_values_from_file(ini)
        # Reindex to match the order and content of ids, fill missing with NaN
        ini_values = ini_values.reindex(ids)
    elif isinstance(ini, pd.Series):
        # Reindex to match the order and content of ids, fill missing with NaN
        ini_values = ini.reindex(ids)
    else:
        raise ValueError('initial values should be in pd.Series or a valid file path')

    return ini_values


def read_measurements_from_file(file, inst_data=False):
    '''
    Reads measurement data from .tsv or Excel files.
    The file is expected to have a header row (first non-commented line).
    Columns are identified by attempting to match common names for fragment IDs,
    time (if applicable), mean, standard deviation, and experiment ID.

    Parameters
    ----------
    file: str
        Path to the .tsv or Excel file.
    inst_data: bool, optional
        If True, expects a 'time' column and sets a multi-index
        ('fragment_ID', 'time'). Otherwise, sets 'fragment_ID' as index.
        Default is False.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with 'fragment_ID' (and 'time' if `inst_data=True`) as index.
        Data columns will be 'mean', 'sd', and 'experiment_id' (if found in the file).

    Raises
    ------
    TypeError
        If the file is not a .tsv or Excel file.
    ValueError
        If required columns (fragment ID, mean, sd, time (if `inst_data`)) are not found.
    '''
    ext = splitext(file)[1]
    # Use pandas to read, assuming header is the first non-commented line
    if re.search(r'tsv', ext, re.IGNORECASE): # Added IGNORECASE
        df = pd.read_csv(file, sep='\t', comment='#', header=0)
    elif re.search(r'xls(?:x|m)?', ext, re.IGNORECASE): # Added IGNORECASE and support for xlsx/xlsm
        df = pd.read_excel(file, comment='#', header=0)
    else:
        raise TypeError(f"Unsupported file type: {ext}. Can only read from .tsv or Excel files.")

    df = df.dropna(how='all') # Drop rows where all elements are NaN

    # Normalize column names
    rename_map = {}
    # Using sets for faster lookups
    potential_fragment_id_names = {'fragment_id', 'fragmentid', 'id', 'emu_id', 'emuid'}
    potential_time_names = {'time', 'timepoint', 'time_point'}
    potential_mean_names = {'mean', 'average'}
    potential_sd_names = {'sd', 'stdev', 'stddev', 'standard_deviation', 'std'}
    potential_experiment_id_names = {'experiment_id', 'experimentid', 'exp_id', 'expid', 'experiment'}

    current_columns = list(df.columns) # Keep original for error messages if needed

    for col in current_columns:
        col_lower_stripped = str(col).lower().replace('_', '').replace(' ', '') # Ensure col is string
        if col_lower_stripped in potential_fragment_id_names:
            rename_map[col] = 'fragment_ID'
        elif col_lower_stripped in potential_time_names:
            rename_map[col] = 'time'
        elif col_lower_stripped in potential_mean_names:
            rename_map[col] = 'mean'
        elif col_lower_stripped in potential_sd_names:
            rename_map[col] = 'sd'
        elif col_lower_stripped in potential_experiment_id_names:
            rename_map[col] = 'experiment_id'

    df.rename(columns=rename_map, inplace=True)

    # Validate required columns
    required_cols = ['fragment_ID', 'mean', 'sd']
    if inst_data:
        required_cols.append('time')

    missing_req_cols = [req_col for req_col in required_cols if req_col not in df.columns]
    if missing_req_cols:
        # Try to find the original names for a better error message
        original_column_names_info = []
        for req_col in missing_req_cols:
            found_original = False
            for original_name, mapped_name in rename_map.items():
                if mapped_name == req_col and original_name in current_columns : # check if original_name was actually in the df
                    original_column_names_info.append(f"'{req_col}' (tried to map from '{original_name}')")
                    found_original = True
                    break
            if not found_original:
                 original_column_names_info.append(f"'{req_col}' (standard name, or mapping source not found/unclear from original: {current_columns})")
        
        raise ValueError(
            f"Required column(s) {', '.join(original_column_names_info)} not found in file {file}. "
            f"Detected columns after normalization: {list(df.columns)}. "
            f"Original columns read from file: {current_columns}."
        )

    # Set index
    index_cols_to_set = ['fragment_ID']
    if inst_data:
        if 'time' not in df.columns: # Should be caught by missing_req_cols above, but defensive check
             raise ValueError(f"'time' column is required for inst_data=True but not found after normalization in {file}.")
        index_cols_to_set.append('time')
        # Convert time column to float for consistent indexing if it's not already
        try:
            df['time'] = df['time'].astype(float)
        except ValueError as e:
            raise ValueError(f"Could not convert 'time' column to float in {file}. Error: {e}")


    try:
        df.set_index(index_cols_to_set, inplace=True)
    except KeyError as e:
        raise ValueError(f"Failed to set index with columns {index_cols_to_set}. One or more not found. Error: {e}")


    # Select and order data columns
    final_data_columns = []
    if 'mean' in df.columns: # Should always be true due to required_cols check
        final_data_columns.append('mean')
    if 'sd' in df.columns: # Should always be true
        final_data_columns.append('sd')
    if 'experiment_id' in df.columns: # Optional
        final_data_columns.append('experiment_id')

    # Ensure no other columns are present by re-assigning df
    # This also handles cases where a column was in rename_map but not in final_data_columns
    try:
        df = df[final_data_columns]
    except KeyError as e:
        # This might happen if 'mean' or 'sd' were somehow dropped or misnamed after rename and before this selection
        raise ValueError(f"Error selecting final data columns ('mean', 'sd', 'experiment_id' if present). Missing: {e}. Current df columns: {list(df.columns)}")

    return df
