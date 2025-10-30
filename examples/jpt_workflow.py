#!/usr/bin/env python3
"""
Example Workflow: Processing JPT Dataset

Processes JPT mass spectrometry and peptide concentration data using the
reusable scripts/process_jpt_data module. Follows the heavy_1st workflow pattern:
1. Loads Skyline PRM export and dilution metadata
2. Aggregates MS intensities by experimental condition
3. Computes regression metrics (R², slope, intercept)
4. Generates QC summary with per-fragment analysis
5. Saves results for plotting

Usage:
    python examples/jpt_workflow.py

Or run from the repository root:
    python -m examples.jpt_workflow
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from scripts.process_jpt_data import process_jpt_data


def run_jpt_workflow():
    """
    Demonstrates complete workflow for processing JPT dataset.
    
    This workflow:
    1. Loads the JPT Skyline export data
    2. Loads JPT peptide concentration metadata
    3. Aggregates MS intensities by experimental condition
    4. Merges with concentration data for calibration
    5. Computes linear regression metrics per fragment/condition
    6. Generates comprehensive QC report
    7. Saves processed results
    """
    
    print("="*70)
    print("JPT DATASET PROCESSING WORKFLOW")
    print("="*70)
    
    # Define file paths
    data_dir = project_root / "data"
    input_dir = data_dir / "input"
    output_dir = data_dir / "output"
    output_dir.mkdir(exist_ok=True)
    
    ms_file = input_dir / "JPT_1_3_-_HHT_1.2_-_Cytosolic.csv"
    conc_file = input_dir / "JPT1-3_peptide_conc.csv"
    output_file = output_dir / "jpt_workflow_results.csv"
    
    # Check if input files exist
    if not ms_file.exists():
        print(f"❌ Error: MS data file not found: {ms_file}")
        print(f"   Please ensure the file exists in the data/input directory.")
        return False
    
    if not conc_file.exists():
        print(f"❌ Error: Concentration data file not found: {conc_file}")
        print(f"   Please ensure the file exists in the data/input directory.")
        return False
    
    print(f"\n📁 Input files:")
    print(f"   MS data: {ms_file.name}")
    print(f"   Concentration data: {conc_file.name}")
    print(f"\n📊 Output file: {output_file.name}")
    print("-" * 70)
    
    try:
        # Process the data
        print("\n🔬 Starting JPT data processing...")
        result_df = process_jpt_data(str(ms_file), str(conc_file), str(output_file))
        
        # Analysis summary
        print("\n" + "="*70)
        print("✅ PROCESSING COMPLETE - ANALYSIS SUMMARY")
        print("="*70)
        
        # Basic statistics
        print(f"\n📊 Dataset Statistics:")
        print(f"   • Total rows processed: {result_df.shape[0]}")
        print(f"   • Total columns: {result_df.shape[1]}")
        print(f"   • Unique peptides: {result_df['peptide'].nunique()}")
        print(f"   • Unique conditions: {result_df['condition'].nunique()}")
        print(f"   • Unique fragments: {result_df['fragment_ion_charged'].nunique()}")
        print(f"   • Days analyzed: {sorted(result_df['day'].unique())}")
        
        # Regression quality assessment
        if 'R2' in result_df.columns:
            # Get unique regression categories and their R² values
            unique_regs = result_df.groupby('regression_category')['R2'].first()
            mean_r2 = unique_regs.mean()
            
            print(f"\n🎯 Regression Quality ({len(unique_regs)} fragment/condition combinations):")
            print(f"   • Average R² across regressions: {mean_r2:.3f}")
            
            # R² quality interpretation
            if mean_r2 >= 0.8:
                quality = "Excellent (R² ≥ 0.8) 🌟"
            elif mean_r2 >= 0.6:
                quality = "Good (R² ≥ 0.6) ✅"
            elif mean_r2 >= 0.4:
                quality = "Fair (R² ≥ 0.4) ⚠️"
            else:
                quality = "Poor (R² < 0.4) ❌"
            print(f"   • Quality: {quality}")
            
            # Top performers
            top_r2 = unique_regs.nlargest(5)
            print(f"\n   Top 5 regressions by R²:")
            for cat, r2 in top_r2.items():
                print(f"      • {cat}: R² = {r2:.3f}")
        
        # Key output columns summary
        key_columns = ['area_mean', 'concentration_ng_mL', 'R2', 'intercept', 'slope', 
                      'mean_r2', 'mean_slope', 'mean_intercept']
        present_key_cols = [col for col in key_columns if col in result_df.columns]
        
        print(f"\n📋 Key Analysis Columns Generated: {len(present_key_cols)}/{len(key_columns)}")
        for col in present_key_cols:
            print(f"   ✓ {col}")
        
        # Concentration range
        if 'concentration_ng_mL' in result_df.columns:
            conc_range = result_df['concentration_ng_mL'].describe()
            print(f"\n📈 Concentration Range (ng/mL):")
            print(f"   • Min: {conc_range['min']:.3e}")
            print(f"   • Mean: {conc_range['mean']:.3e}")
            print(f"   • Max: {conc_range['max']:.3e}")
        
        # Sample data preview
        print(f"\n🔍 Sample Results (first 3 rows):")
        display_cols = ['peptide', 'condition', 'fragment_ion_charged', 'concentration_ng_mL', 'area_mean', 'R2', 'slope']
        available_display_cols = [col for col in display_cols if col in result_df.columns]
        
        sample_data = result_df[available_display_cols].drop_duplicates(subset=['peptide', 'condition', 'fragment_ion_charged']).head(3)
        for idx, row in sample_data.iterrows():
            print(f"   Row {idx + 1}:")
            for col in available_display_cols:
                val = row[col]
                if isinstance(val, float):
                    print(f"      {col}: {val:.3e}" if abs(val) < 0.001 else f"      {col}: {val:.3f}")
                else:
                    print(f"      {col}: {val}")
        
        print(f"\n💾 Results saved to: {output_file}")
        print(f"📈 Ready for plotting and visualization!")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during processing: {e}")
        print("\n🔧 Troubleshooting tips:")
        print("   1. Check that input files are valid CSV format")
        print("   2. Ensure peptide names match between MS and concentration files")
        print("   3. Verify replicate names follow expected format (D#_...)")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main function for running the workflow."""
    success = run_jpt_workflow()
    
    if success:
        print("\n🎉 Workflow completed successfully!")
        print("\n📚 Next steps:")
        print("   • Run: python scripts/plot_jpt.py")
        print("   • This will generate calibration curves under reports/jpt/plots/")
        print("   • Examine regression metrics to identify high-quality peptides")
        print("   • Filter results based on R² thresholds for downstream analysis")
    else:
        print("\n💥 Workflow failed. Please check the error messages above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
