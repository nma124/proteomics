#!/usr/bin/env python3
"""
JPT Quick Export Script
=======================
Minimal script to create plots and export CSV files from JPT data

Usage:
    python jpt_quick_export.py
"""

import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
from typing import Tuple

# Configuration
DATA_INPUT_DIR = pathlib.Path("data/input")
DATA_OUTPUT_DIR = pathlib.Path("data/output")
DATA_OUTPUT_DIR.mkdir(exist_ok=True)

MS_FILE = "JPT_1_3_-_HHT_1.2_-_Cytosolic.csv"
CONC_FILE = "JPT1-3_peptide_conc.csv"


def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load MS and concentration data"""
    ms_df = pd.read_csv(DATA_INPUT_DIR / MS_FILE)
    ms_df.columns = ms_df.columns.str.lower().str.strip()
    
    conc_df = pd.read_csv(DATA_INPUT_DIR / CONC_FILE)
    
    return ms_df, conc_df


def extract_conditions(replicate_col: pd.Series) -> pd.DataFrame:
    """Extract experimental conditions from replicate names"""
    df = pd.DataFrame({'replicate': replicate_col})
    df['day'] = replicate_col.str.extract(r'D(\d+)')[0]
    df['condition'] = replicate_col.str.extract(r'(cyto|SDC_Sup)')
    df['wash_type'] = replicate_col.str.extract(r'_(W{1,2})_')
    return df


def process_ms_data(ms_df: pd.DataFrame) -> pd.DataFrame:
    """Process and aggregate MS data by condition"""
    conditions = extract_conditions(ms_df['replicate'])
    ms_combined = pd.concat([ms_df, conditions[['day', 'condition', 'wash_type']]], axis=1)
    
    # Filter precursor ions
    precursor = ms_combined[ms_combined['fragment ion'] == 'precursor'].copy()
    
    # Aggregate
    agg = precursor.groupby(['peptide', 'protein', 'day', 'condition', 'wash_type'])['area'].agg(
        ['mean', 'std', 'count', 'sum']
    ).reset_index()
    
    agg.columns = ['peptide', 'protein', 'day', 'condition', 'wash_type', 'area_mean', 'area_std', 'count', 'area_sum']
    return agg


def process_concentration_data(conc_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process concentration data to long format and calculate fold changes"""
    # Long format
    conc_long = pd.melt(
        conc_df,
        id_vars=['Peptides'],
        var_name='timepoint',
        value_name='concentration_ng_mL'
    )
    conc_long['day'] = pd.to_numeric(conc_long['timepoint'].str.extract(r'D(\d+)')[0], errors='coerce')
    
    # Fold changes
    pivot = conc_long.pivot(index='Peptides', columns='day', values='concentration_ng_mL')
    fold_change = pivot.div(pivot.iloc[:, 0], axis=0)
    fold_change = fold_change.replace([np.inf, -np.inf], np.nan)
    
    return conc_long, fold_change


def create_plots(ms_df: pd.DataFrame, conc_df: pd.DataFrame, conc_long: pd.DataFrame):
    """Create summary visualization"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('JPT Data Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Fragment ion types
    fragment_counts = ms_df['fragment ion'].value_counts()
    axes[0, 0].bar(range(len(fragment_counts)), fragment_counts.values, color='steelblue')
    axes[0, 0].set_xticks(range(len(fragment_counts)))
    axes[0, 0].set_xticklabels(fragment_counts.index, rotation=45, ha='right')
    axes[0, 0].set_title('Fragment Ion Distribution')
    axes[0, 0].set_ylabel('Count')
    
    # Plot 2: Area distribution (log scale)
    log_areas = np.log10(ms_df['area'] + 1)
    axes[0, 1].hist(log_areas, bins=50, alpha=0.7, color='teal', edgecolor='black')
    axes[0, 1].set_title('Peak Area Distribution (Log10)')
    axes[0, 1].set_xlabel('Log10(Area + 1)')
    axes[0, 1].set_ylabel('Frequency')
    
    # Plot 3: Top peptides
    top_pep = ms_df['peptide'].value_counts().head(10)
    axes[0, 2].barh(range(len(top_pep)), top_pep.values, color='coral')
    axes[0, 2].set_yticks(range(len(top_pep)))
    axes[0, 2].set_yticklabels(top_pep.index, fontsize=9)
    axes[0, 2].set_title('Top 10 Peptides')
    axes[0, 2].set_xlabel('Count')
    
    # Plot 4: Peptide concentration heatmap
    conc_matrix = conc_df.set_index('Peptides').iloc[:10]
    im = axes[1, 0].imshow(conc_matrix.values, aspect='auto', cmap='RdYlGn')
    axes[1, 0].set_xticks(range(len(conc_matrix.columns)))
    axes[1, 0].set_xticklabels(conc_matrix.columns, rotation=45, ha='right', fontsize=9)
    axes[1, 0].set_yticks(range(len(conc_matrix.index)))
    axes[1, 0].set_yticklabels(conc_matrix.index, fontsize=9)
    axes[1, 0].set_title('Peptide Concentrations (ng/mL)')
    plt.colorbar(im, ax=axes[1, 0])
    
    # Plot 5: Time course
    top_pep_conc = conc_long['Peptides'].unique()[:5]
    for peptide in top_pep_conc:
        pep_data = conc_long[conc_long['Peptides'] == peptide].sort_values('day')
        label = peptide[:15] + '...' if len(peptide) > 15 else peptide
        axes[1, 1].plot(pep_data['day'], pep_data['concentration_ng_mL'], marker='o', label=label, linewidth=2)
    
    axes[1, 1].set_title('Concentration Time Course (Top 5 Peptides)')
    axes[1, 1].set_xlabel('Day')
    axes[1, 1].set_ylabel('Concentration (ng/mL)')
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].set_yscale('log')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Plot 6: Summary stats
    stats_text = f"""Summary Statistics

MS Data:
  Total measurements: {len(ms_df):,}
  Unique peptides: {ms_df['peptide'].nunique()}
  Unique proteins: {ms_df['protein'].nunique()}
  Unique replicates: {ms_df['replicate'].nunique()}

Area (Peak Intensity):
  Mean: {ms_df['area'].mean():.2e}
  Median: {ms_df['area'].median():.2e}
  Std Dev: {ms_df['area'].std():.2e}
  Min: {ms_df['area'].min():.2e}
  Max: {ms_df['area'].max():.2e}

Concentration Data:
  Peptides: {len(conc_df)}
  Time points: {conc_df.shape[1] - 1}
"""
    axes[1, 2].text(0.05, 0.95, stats_text, transform=axes[1, 2].transAxes,
                    fontsize=9, verticalalignment='top', fontfamily='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    
    # Save figure
    output_path = DATA_OUTPUT_DIR / "jpt_data_summary.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved: {output_path}")
    plt.close()


def export_csv_files(ms_agg: pd.DataFrame, conc_long: pd.DataFrame, fold_change: pd.DataFrame):
    """Export processed data to CSV"""
    
    # 1. Aggregated MS data by condition
    output_path = DATA_OUTPUT_DIR / "jpt_peptide_areas_by_condition.csv"
    ms_agg.to_csv(output_path, index=False)
    print(f"✓ CSV saved: {output_path}")
    
    # 2. Concentration timecourse
    output_path = DATA_OUTPUT_DIR / "jpt_concentration_timecourse.csv"
    conc_long.to_csv(output_path, index=False)
    print(f"✓ CSV saved: {output_path}")
    
    # 3. Fold changes
    output_path = DATA_OUTPUT_DIR / "jpt_concentration_fold_changes.csv"
    fold_change.to_csv(output_path)
    print(f"✓ CSV saved: {output_path}")


def main():
    """Main execution"""
    print("="*60)
    print("JPT QUICK EXPORT")
    print("="*60)
    
    print("\n1. Loading data...")
    ms_df, conc_df = load_data()
    print(f"   MS data: {len(ms_df)} measurements")
    print(f"   Concentration data: {len(conc_df)} peptides")
    
    print("\n2. Processing data...")
    ms_agg = process_ms_data(ms_df)
    conc_long, fold_change = process_concentration_data(conc_df)
    print(f"   Aggregated to {len(ms_agg)} peptide-condition combinations")
    
    print("\n3. Creating plots...")
    create_plots(ms_df, conc_df, conc_long)
    
    print("\n4. Exporting CSV files...")
    export_csv_files(ms_agg, conc_long, fold_change)
    
    print("\n" + "="*60)
    print("COMPLETE!")
    print(f"Output directory: {DATA_OUTPUT_DIR.absolute()}")
    print("="*60)


if __name__ == "__main__":
    main()
