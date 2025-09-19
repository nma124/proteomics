#!/usr/bin/env python3
"""
Proteomics PRM Data Processing Script

Processes Skyline PRM export data with heavy peptide dilution scheme to produce 
quantification results with regression analysis and quality control metrics.

Based on the analysis pipeline from notebook 005-ms-peggy-2-update.ipynb
"""

import pandas as pd
import numpy as np
from sklearn import linear_model, metrics
import pathlib
import warnings

# Suppress sklearn warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

def fix_dilution(row: str):
    """Fix dilution numbering"""
    str_to_check = row.split('_')[0]
    if len(str_to_check) == 1:
        dilution = 'D' + row
    else:
        dilution = row
    return dilution

def grp_area_ratio_n_area_cols(grouped_df: pd.DataFrame):
    """Compute the ratio of areas with higher precursor mz value as numerator."""
    df = grouped_df[['precursor mz', 'area']].copy()
    df = df.sort_values(by='precursor mz', ascending=True)
    
    if len(df) < 2:
        return pd.Series({'area_ratio': np.nan, 'area_min': np.nan, 'area_max': np.nan})
    
    area_ratio = df['area'].iloc[1] / df['area'].iloc[0] if df['area'].iloc[0] != 0 else np.inf
    area_min, area_max = df['area'].iloc[0], df['area'].iloc[1]

    col_vals = {'area_ratio': area_ratio, 'area_min': area_min, 'area_max': area_max}
    return pd.Series(col_vals)

def get_dilution_from_replicate(rep_string: str):
    """Get the dilution symbol from the replicate name."""
    str_to_check = rep_string.split('_')[0]
    if str_to_check.startswith('col'):
        dilution = 'D0'
    else:
        dilution = str_to_check
    return dilution

def get_peptide_dilution_comb(row):
    """Function to create a new column from peptide and dilution name."""
    peptide = row['peptide']
    replicate = row['replicate']
    dilution = get_dilution_from_replicate(rep_string=replicate)
    return peptide + '_' + dilution

def get_rep_plot_cat(row):
    """Get the plot category."""
    replicate = row['replicate']
    peptide = row['peptide']
    fragment_ion = row['fragment ion']

    rep_elements = replicate.split('_')[0:3]
    cat_name = [peptide] + [fragment_ion]  # + rep_elements
    cat_name = '_'.join(cat_name)
    return cat_name

def get_linear_fit(data: pd.DataFrame, plot_cat: str):
    """Get linear regression fit parameters for a plot category."""
    cat_data = data[data['plot_cat'] == plot_cat]
    
    if len(cat_data) < 2:
        return (np.nan, np.nan, np.nan)
    
    x = cat_data['area_ratio'].values
    y = cat_data['heavy_conc'].values
    
    # Remove infinite and NaN values
    valid_mask = np.isfinite(x) & np.isfinite(y)
    x = x[valid_mask]
    y = y[valid_mask]
    
    if len(x) < 2:
        return (np.nan, np.nan, np.nan)

    x = x[:, np.newaxis]
    y = y[:, np.newaxis]

    model = linear_model.LinearRegression()
    model.fit(x, y)
    y_fit = model.predict(x)

    r2 = metrics.r2_score(y, y_fit)
    intercept = model.intercept_.squeeze()
    grad = model.coef_[0].squeeze()

    return (r2, intercept, grad)

def get_fit_param_cat(row):
    """Get the fit param category."""
    replicate = row['replicate']
    peptide = row['peptide']

    rep_elements = replicate.split('_')[1:3]
    cat_name = [peptide] + rep_elements
    cat_name = '_'.join(cat_name)
    return cat_name

def qtest(data, right=True):
    """Simple Q-test implementation."""
    sorted_data = sorted(data)
    
    if len(sorted_data) < 2:
        return 0.0

    if right:
        gap = sorted_data[-1] - sorted_data[-2]
    else:
        gap = sorted_data[1] - sorted_data[0]

    try:
        range_val = sorted_data[-1] - sorted_data[0]
        q_val = gap / range_val if range_val != 0 else gap
    except ZeroDivisionError:
        q_val = gap

    return q_val

def get_grp_fit_agg(grouped_df: pd.DataFrame):
    """Get the aggregate values for each group."""
    df = grouped_df[['intercept', 'gradient', 'R2']].copy()
    
    # Convert to numeric and handle mixed types
    df['gradient'] = pd.to_numeric(df['gradient'], errors='coerce')
    df['intercept'] = pd.to_numeric(df['intercept'], errors='coerce')
    df['R2'] = pd.to_numeric(df['R2'], errors='coerce')
    
    gradient = df['gradient'].dropna().values
    intercept = df['intercept'].dropna().values
    r2 = df['R2'].dropna().values

    # Additional check for finite values
    try:
        gradient = gradient[np.isfinite(gradient)]
    except TypeError:
        gradient = gradient[~pd.isna(gradient)]  # fallback for mixed types
    
    try:
        intercept = intercept[np.isfinite(intercept)]
    except TypeError:
        intercept = intercept[~pd.isna(intercept)]
    
    try:
        r2 = r2[np.isfinite(r2)]
    except TypeError:
        r2 = r2[~pd.isna(r2)]

    mean_grad = gradient.mean() if len(gradient) > 0 else np.nan
    stdv_grad = gradient.std() if len(gradient) > 1 else np.nan
    cov_grad = stdv_grad / mean_grad if mean_grad != 0 else np.nan
    qtest_grad_right = qtest(gradient, right=True) if len(gradient) > 1 else 0.0

    mean_intercept = intercept.mean() if len(intercept) > 0 else np.nan
    stdv_intercept = intercept.std() if len(intercept) > 1 else np.nan
    cov_intercept = stdv_intercept / mean_intercept if mean_intercept != 0 else np.nan
    qtest_intercept_right = qtest(intercept, right=True) if len(intercept) > 1 else 0.0

    mean_r2 = r2.mean() if len(r2) > 0 else np.nan
    stdv_r2 = r2.std() if len(r2) > 1 else np.nan
    cov_r2 = stdv_r2 / mean_r2 if mean_r2 != 0 else np.nan
    qtest_r2_right = qtest(r2, right=True) if len(r2) > 1 else 0.0

    cols_dict = {
        'mean_grad': mean_grad, 'stdv_grad': stdv_grad,
        'cov_grad': cov_grad, 'qtest_grad': qtest_grad_right,
        'mean_intercept': mean_intercept, 'stdv_intercept': stdv_intercept,
        'cov_intercept': cov_intercept, 'qtest_intercept': qtest_intercept_right,
        'mean_r2': mean_r2, 'stdv_r2': stdv_r2,
        'cov_r2': cov_r2, 'qtest_r2': qtest_r2_right
    }

    return pd.Series(cols_dict)

def process_prm_data(skyline_file: str, dilution_file: str, output_file: str):
    """
    Main processing function to analyze PRM data.
    
    Args:
        skyline_file: Path to Skyline PRM export CSV
        dilution_file: Path to peptide dilution concentration CSV  
        output_file: Path for output CSV
    """
    
    print("Loading Skyline PRM data...")
    # Load Skyline data
    trans_res_df_raw = pd.read_csv(skyline_file)
    trans_res_df_raw.columns = trans_res_df_raw.columns.map(str.lower)
    
    print(f"Loaded {trans_res_df_raw.shape[0]} rows of raw data")
    print(f"Peptides found: {trans_res_df_raw.peptide.unique()}")
    
    # Fix replicate naming
    temp_df = trans_res_df_raw.copy()
    temp_df['replicate'] = trans_res_df_raw['replicate'].apply(fix_dilution)
    
    # Group by peptide, replicate, product charge, fragment ion and filter for pairs
    cols_to_group_by = ['peptide', 'replicate', 'product charge', 'fragment ion']
    print(f"\nGroup counts: {temp_df.groupby(cols_to_group_by)['area'].count().value_counts()}")
    
    # Filter groups that contain exactly two elements (for area ratio calculation)
    temp_df = temp_df.groupby(cols_to_group_by).filter(lambda g: g['area'].count() == 2)
    print(f"After filtering for pairs: {temp_df.shape[0]} rows")
    
    # Calculate area ratios
    print("Calculating area ratios...")
    df_area_ratio = temp_df.groupby(cols_to_group_by).apply(grp_area_ratio_n_area_cols).reset_index()
    print(f"Area ratios calculated: {df_area_ratio.shape[0]} rows")
    
    # Load peptide dilution concentrations
    print("Loading peptide dilution concentrations...")
    df_peptide_dilution_conc = pd.read_csv(dilution_file)
    
    # Reshape the dilution data
    df_peptide_dilution_conc = df_peptide_dilution_conc.melt(
        id_vars='Peptides',
        value_vars=['D1 (ng/mL)', 'D2 (ng/mL)', 'D3 (ng/mL)', 'D4 (ng/mL)', 
                   'D5 (ng/mL)', 'D6 (ng/mL)', 'D7 (ng/mL)', 'D0  (ng/mL)'],
        var_name='dilution',
        value_name='heavy_conc'
    )
    
    # Clean up dilution column names
    df_peptide_dilution_conc['dilution'] = df_peptide_dilution_conc['dilution'].str.extract(r'(D\d+)')
    
    # Create peptide_dilution_name column
    peptide_dilution_name_col = 'peptide_dilution_name'
    df_peptide_dilution_conc['peptide_dilution_name'] = df_peptide_dilution_conc['Peptides'] + '_' + df_peptide_dilution_conc['dilution']
    
    print(f"Dilution data shape: {df_peptide_dilution_conc.shape}")
    
    # Add peptide_dilution_name to area ratio data
    df_area_ratio[peptide_dilution_name_col] = df_area_ratio.apply(get_peptide_dilution_comb, axis=1)
    
    # Merge area ratios with concentrations
    print("Merging area ratios with concentrations...")
    df_area_ratio_conc = pd.merge(df_area_ratio, right=df_peptide_dilution_conc, how='left',
                                 on=peptide_dilution_name_col, suffixes=('', '_y'))
    df_area_ratio_conc.drop(df_area_ratio_conc.filter(regex='_y$').columns, axis=1, inplace=True)
    df_area_ratio_conc = df_area_ratio_conc.dropna(axis=0)  # drop nan

    # Drop infinite values
    df_area_ratio_conc.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_area_ratio_conc.dropna(subset=["area_ratio"], how="all", inplace=True)
    
    print(f"After merging and cleaning: {df_area_ratio_conc.shape[0]} rows")
    
    # Create plot categories
    df_area_ratio_conc['plot_cat'] = df_area_ratio_conc.apply(get_rep_plot_cat, axis=1)
    df_area_ratio_conc['plot_cat_grp'] = df_area_ratio_conc['plot_cat'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    df_area_ratio_conc['order_comp'] = df_area_ratio_conc['replicate'].apply(lambda x: x.split('_')[-1])
    df_area_ratio_conc['rep_partial'] = df_area_ratio_conc['replicate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    
    # Get linear fit parameters
    print("Calculating linear regression parameters...")
    plot_cats = df_area_ratio_conc['plot_cat'].unique()
    
    cats = []
    r2s = []
    intercepts = []
    grads = []

    for cat in plot_cats:
        r2, intercept, grad = get_linear_fit(df_area_ratio_conc, cat)
        cats.append(cat)
        r2s.append(r2)
        intercepts.append(intercept)
        grads.append(grad)

    df_plot_cat_fit_params = pd.DataFrame({
        'plot_cat': cats,
        'R2': r2s,
        'intercept': intercepts,
        'gradient': grads
    })

    # Merge fit parameters
    df_area_ratio_conc_n_fit_params = pd.merge(df_area_ratio_conc,
                                              right=df_plot_cat_fit_params,
                                              on='plot_cat',
                                              how='left')
    
    # Get fit parameter aggregates
    print("Calculating fit parameter aggregates...")
    temp = df_area_ratio_conc_n_fit_params.copy()
    temp['fit_param_grp'] = temp.apply(get_fit_param_cat, axis=1)

    cols_to_group_by = ['fit_param_grp']
    df_fit_param_agg = temp.groupby(cols_to_group_by).apply(get_grp_fit_agg).reset_index()

    df_w_fit_param_agg = pd.merge(left=temp, right=df_fit_param_agg, how='left', on=cols_to_group_by)
    
    print(f"Final output shape: {df_w_fit_param_agg.shape}")
    
    # Save output
    print(f"Saving results to {output_file}...")
    df_w_fit_param_agg.to_csv(output_file, index=True)
    
    print("Processing complete!")
    return df_w_fit_param_agg

def main():
    """Main function to run the processing pipeline."""
    # Define file paths
    data_dir = pathlib.Path("data/raw/peggy_2")
    skyline_file = data_dir / "skyline raw data_heavy_di_expanded_cleaned_AHS_3+only_rmHL0.csv"
    dilution_file = data_dir / "peptide_dilution_conc_peggy.csv"
    output_file = "Akin_PRM_heavy_di_expanded_AHS_3+only_rmHL0_output.csv"
    
    # Check if input files exist
    if not skyline_file.exists():
        print(f"Error: Skyline data file not found: {skyline_file}")
        return
    
    if not dilution_file.exists():
        print(f"Error: Dilution data file not found: {dilution_file}")
        return
    
    # Process the data
    result_df = process_prm_data(str(skyline_file), str(dilution_file), output_file)
    
    # Print summary statistics
    print(f"\n=== SUMMARY ===")
    print(f"Output file: {output_file}")
    print(f"Shape: {result_df.shape[0]} rows × {result_df.shape[1]} columns")
    print(f"Columns: {list(result_df.columns)}")
    print(f"Peptides processed: {result_df['peptide'].nunique()}")
    print(f"Fragment ions: {result_df['fragment ion'].nunique()}")
    
    # Check for key regression columns
    regression_cols = ['R2', 'intercept', 'gradient', 'mean_grad', 'stdv_grad', 'cov_grad', 
                      'mean_intercept', 'stdv_intercept', 'cov_intercept', 'mean_r2', 'stdv_r2', 'cov_r2']
    present_cols = [col for col in regression_cols if col in result_df.columns]
    print(f"Regression/QC columns present: {len(present_cols)}/{len(regression_cols)}")
    print(f"Present: {present_cols}")

if __name__ == "__main__":
    main()