#!/usr/bin/env python3
"""
Script to run the heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv file 
through the process_prm_data.py pipeline.
"""

import sys
import pathlib
from .process_prm_data import process_prm_data

def main():
    """Run the processing pipeline with the heavy_1st input file."""
    
    # Define file paths relative to project root (data/input and data/output)
    project_root = pathlib.Path(__file__).resolve().parent.parent
    input_dir = project_root / "data" / "input"
    output_dir = project_root / "data" / "output"

    skyline_file = input_dir / "heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv"
    dilution_file = input_dir / "peptide_dilution_conc_peggy.csv"
    output_file = output_dir / "heavy_1st_processed_output.csv"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if input files exist
    if not skyline_file.exists():
        print(f"Error: Skyline data file not found: {skyline_file}")
        print("Hint: Place the file in data/input or update the path accordingly.")
        return 1
    
    if not dilution_file.exists():
        print(f"Error: Dilution data file not found: {dilution_file}")
        print("Hint: Place the file in data/input or update the path accordingly.")
        return 1
    
    print(f"Processing {skyline_file} with dilution data from {dilution_file}")
    print(f"Output will be saved to: {output_file}")
    print("-" * 60)
    
    # Process the data
    try:
        result_df = process_prm_data(str(skyline_file), str(dilution_file), str(output_file))
        
        # Print summary statistics
        print(f"\n=== PROCESSING COMPLETE ===")
        print(f"Input file: {skyline_file}")
        print(f"Output file: {output_file}")
        print(f"Shape: {result_df.shape[0]} rows × {result_df.shape[1]} columns")
        print(f"Peptides processed: {result_df['peptide'].nunique()}")
        print(f"Fragment ions: {result_df['fragment ion'].nunique()}")
        print(f"Dilutions: {result_df['dilution'].nunique()}")
        print(f"Replicates: {result_df['replicate'].nunique()}")
        
        # Check regression quality
        if 'mean_r2' in result_df.columns:
            mean_r2 = result_df['mean_r2'].mean()
            print(f"Average R² across regressions: {mean_r2:.3f}")
        
        # Show sample of results
        print(f"\nFirst few rows of key columns:")
        key_cols = ['peptide', 'replicate', 'fragment ion', 'area_ratio', 'heavy_conc', 'R2']
        available_cols = [col for col in key_cols if col in result_df.columns]
        print(result_df[available_cols].head())
        
        return 0
        
    except Exception as e:
        print(f"Error during processing: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())