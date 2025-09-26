#!/usr/bin/env python3
"""
Main entry point for Proteomics PRM Data Processing

This script provides a command-line interface for processing Skyline PRM export data.
"""

import sys
import argparse
from pathlib import Path

# Add scripts directory to path for imports
sys.path.append(str(Path(__file__).parent / "scripts"))

from scripts.process_prm_data import process_prm_data


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description="Process Skyline PRM export data with heavy peptide dilution analysis",
        epilog="Example: python main.py data/input/skyline_data.csv data/input/dilutions.csv -o results.csv"
    )
    
    parser.add_argument("skyline_file", help="Path to Skyline PRM export CSV file")
    parser.add_argument("dilution_file", help="Path to peptide dilution concentration CSV file")
    parser.add_argument("-o", "--output", default="prm_analysis_output.csv", 
                       help="Output CSV file path (default: prm_analysis_output.csv)")
    parser.add_argument("--version", action="version", version="PRM Processing v1.0.0")
    
    args = parser.parse_args()
    
    # Validate input files
    skyline_path = Path(args.skyline_file)
    dilution_path = Path(args.dilution_file)
    
    if not skyline_path.exists():
        print(f"Error: Skyline data file not found: {skyline_path}")
        sys.exit(1)
    
    if not dilution_path.exists():
        print(f"Error: Dilution data file not found: {dilution_path}")
        sys.exit(1)
    
    # Process the data
    print(f"Processing {skyline_path} with dilution data from {dilution_path}")
    print(f"Output will be saved to: {args.output}")
    print("-" * 60)
    
    try:
        result_df = process_prm_data(str(skyline_path), str(dilution_path), args.output)
        
        # Print summary
        print(f"\n=== PROCESSING COMPLETE ===")
        print(f"Input file: {skyline_path}")
        print(f"Output file: {args.output}")
        print(f"Shape: {result_df.shape[0]} rows × {result_df.shape[1]} columns")
        print(f"Peptides processed: {result_df['peptide'].nunique()}")
        print(f"Fragment ions: {result_df['fragment ion'].nunique()}")
        if 'dilution' in result_df.columns:
            print(f"Dilutions: {result_df['dilution'].nunique()}")
        if 'replicate' in result_df.columns:
            print(f"Replicates: {result_df['replicate'].nunique()}")
        
        # Check regression quality
        if 'mean_r2' in result_df.columns:
            mean_r2 = result_df['mean_r2'].mean()
            print(f"Average R² across regressions: {mean_r2:.3f}")
        
        print("\nProcessing completed successfully!")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()