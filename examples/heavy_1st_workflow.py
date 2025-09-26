#!/usr/bin/env python3
"""
Example Workflow: Processing Heavy_1st Dataset

This script demonstrates how to process the heavy_1st_expanded dataset
using the proteomics PRM processing pipeline.

Usage:
    python examples/heavy_1st_workflow.py

Or run from the repository root:
    python -m examples.heavy_1st_workflow
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from scripts.process_prm_data import process_prm_data


def run_heavy_1st_workflow():
    """
    Demonstrates complete workflow for processing heavy_1st dataset.
    
    This workflow:
    1. Loads the heavy_1st Skyline export data
    2. Applies dilution concentration information
    3. Calculates area ratios and performs regression analysis
    4. Generates quality control metrics
    5. Saves processed results
    """
    
    print("="*70)
    print("HEAVY_1ST DATASET PROCESSING WORKFLOW")
    print("="*70)
    
    # Define file paths
    data_dir = project_root / "data"
    input_dir = data_dir / "input"
    output_dir = data_dir / "output"
    
    skyline_file = input_dir / "heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv"
    dilution_file = input_dir / "peptide_dilution_conc_peggy.csv"
    output_file = output_dir / "heavy_1st_workflow_results.csv"
    
    # Check if input files exist
    if not skyline_file.exists():
        print(f"❌ Error: Skyline data file not found: {skyline_file}")
        print(f"Please ensure the file exists in the data/input directory.")
        return False
    
    if not dilution_file.exists():
        print(f"❌ Error: Dilution data file not found: {dilution_file}")
        print(f"Please ensure the file exists in the data/input directory.")
        return False
    
    print(f"📁 Input files:")
    print(f"   Skyline data: {skyline_file.name}")
    print(f"   Dilution data: {dilution_file.name}")
    print(f"📊 Output file: {output_file.name}")
    print("-" * 70)
    
    try:
        # Process the data
        print("🔬 Starting PRM data processing...")
        result_df = process_prm_data(str(skyline_file), str(dilution_file), str(output_file))
        
        # Analysis summary
        print("\\n" + "="*70)
        print("✅ PROCESSING COMPLETE - ANALYSIS SUMMARY")
        print("="*70)
        
        # Basic statistics
        print(f"📊 Dataset Statistics:")
        print(f"   • Total rows processed: {result_df.shape[0]}")
        print(f"   • Total columns: {result_df.shape[1]}")
        print(f"   • Peptides analyzed: {result_df['peptide'].nunique()}")
        print(f"   • Fragment ions: {result_df['fragment ion'].nunique()}")
        
        if 'dilution' in result_df.columns:
            print(f"   • Dilution levels: {result_df['dilution'].nunique()}")
        
        if 'replicate' in result_df.columns:
            print(f"   • Replicates: {result_df['replicate'].nunique()}")
        
        # Regression quality assessment
        if 'mean_r2' in result_df.columns:
            mean_r2 = result_df['mean_r2'].mean()
            print(f"\\n🎯 Regression Quality:")
            print(f"   • Average R² across regressions: {mean_r2:.3f}")
            
            # R² quality interpretation
            if mean_r2 >= 0.8:
                print("   • Quality: Excellent (R² ≥ 0.8) 🌟")
            elif mean_r2 >= 0.6:
                print("   • Quality: Good (R² ≥ 0.6) ✅")
            elif mean_r2 >= 0.4:
                print("   • Quality: Fair (R² ≥ 0.4) ⚠️")
            else:
                print("   • Quality: Poor (R² < 0.4) ❌")
        
        # Key output columns summary
        key_columns = ['area_ratio', 'heavy_conc', 'R2', 'intercept', 'gradient']
        present_key_cols = [col for col in key_columns if col in result_df.columns]
        
        print(f"\\n📋 Key Analysis Columns Generated: {len(present_key_cols)}/{len(key_columns)}")
        for col in present_key_cols:
            print(f"   ✓ {col}")
        
        # Sample data preview
        print(f"\\n🔍 Sample Results (first 3 rows):")
        display_cols = ['peptide', 'replicate', 'fragment ion', 'area_ratio', 'heavy_conc', 'R2']
        available_display_cols = [col for col in display_cols if col in result_df.columns]
        
        sample_data = result_df[available_display_cols].head(3)
        for idx, row in sample_data.iterrows():
            print(f"   Row {idx + 1}: {dict(row)}")
        
        print(f"\\n💾 Results saved to: {output_file}")
        print(f"📈 Ready for downstream analysis and visualization!")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during processing: {e}")
        print("\\nTroubleshooting tips:")
        print("1. Check that input files are valid CSV format")
        print("2. Ensure peptide names match between skyline and dilution files")
        print("3. Verify column names in input files match expected format")
        return False


def main():
    """Main function for running the workflow."""
    success = run_heavy_1st_workflow()
    
    if success:
        print("\\n🎉 Workflow completed successfully!")
        print("\\nNext steps:")
        print("• Examine the output file for detailed results")
        print("• Use the processed data for quantification analysis")
        print("• Create visualizations of calibration curves")
        print("• Apply quality filters based on R² values")
    else:
        print("\\n💥 Workflow failed. Please check the error messages above.")
        sys.exit(1)


if __name__ == "__main__":
    main()