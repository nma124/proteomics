#!/usr/bin/env python3
"""
JPT Files Processing Script
==========================

This script processes JPT (peptide/protein) data files including:
1. Mass spectrometry data with peptide measurements
2. Peptide concentration data across time points

Author: WARP AI Assistant
Date: 2025-01-21
"""

import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class JPTDataProcessor:
    """Class to process JPT peptide and mass spectrometry data"""
    
    def __init__(self, data_dir: str = "data/input"):
        self.data_dir = pathlib.Path(data_dir)
        self.output_dir = pathlib.Path("data/output")
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize data containers
        self.ms_data = None
        self.peptide_conc_data = None
        self.processed_data = {}
    
    def load_ms_data(self, filename: str = "JPT_1_3_-_HHT_1.2_-_Cytosolic.csv") -> pd.DataFrame:
        """Load mass spectrometry data from JPT CSV file"""
        file_path = self.data_dir / filename
        print(f"Loading MS data from: {file_path}")
        
        try:
            self.ms_data = pd.read_csv(file_path)
            # Standardize column names
            self.ms_data.columns = self.ms_data.columns.str.lower().str.strip()
            print(f"Loaded MS data: {self.ms_data.shape[0]} rows, {self.ms_data.shape[1]} columns")
            print(f"Columns: {list(self.ms_data.columns)}")
            return self.ms_data
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            return None
    
    def load_peptide_concentrations(self, filename: str = "JPT1-3_peptide_conc.csv") -> pd.DataFrame:
        """Load peptide concentration data"""
        file_path = self.data_dir / filename
        print(f"Loading peptide concentration data from: {file_path}")
        
        try:
            self.peptide_conc_data = pd.read_csv(file_path)
            print(f"Loaded peptide concentration data: {self.peptide_conc_data.shape[0]} rows, {self.peptide_conc_data.shape[1]} columns")
            print(f"Columns: {list(self.peptide_conc_data.columns)}")
            return self.peptide_conc_data
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            return None
    
    def parse_replicate_info(self, replicate_col: pd.Series) -> pd.DataFrame:
        """Parse replicate column to extract experimental conditions"""
        replicate_df = pd.DataFrame()
        
        # Extract day, condition, and other metadata from replicate names
        replicate_df['replicate'] = replicate_col
        replicate_df['day'] = replicate_col.str.extract(r'D(\d+)')[0]
        replicate_df['condition'] = replicate_col.str.extract(r'(cyto|SDC_Sup)')
        replicate_df['wash_type'] = replicate_col.str.extract(r'_(W{1,2})_')
        replicate_df['volume'] = replicate_col.str.extract(r'(\d+)uL')
        replicate_df['sample_id'] = replicate_col.str.extract(r'_(\d+)$')
        
        return replicate_df
    
    def summarize_ms_data(self) -> Dict:
        """Generate summary statistics for mass spectrometry data"""
        if self.ms_data is None:
            print("No MS data loaded. Call load_ms_data() first.")
            return {}
        
        summary = {}
        
        # Basic statistics
        summary['total_measurements'] = len(self.ms_data)
        summary['unique_peptides'] = self.ms_data['peptide'].nunique()
        summary['unique_proteins'] = self.ms_data['protein'].nunique()
        summary['unique_replicates'] = self.ms_data['replicate'].nunique()
        
        # Fragment ion types
        summary['fragment_types'] = self.ms_data['fragment ion'].value_counts().to_dict()
        
        # Area statistics
        summary['area_stats'] = {
            'mean': self.ms_data['area'].mean(),
            'median': self.ms_data['area'].median(),
            'std': self.ms_data['area'].std(),
            'min': self.ms_data['area'].min(),
            'max': self.ms_data['area'].max()
        }
        
        # Peptide frequency
        summary['top_peptides'] = self.ms_data['peptide'].value_counts().head(10).to_dict()
        
        return summary
    
    def process_peptide_areas_by_condition(self) -> pd.DataFrame:
        """Process and aggregate peptide areas by experimental conditions"""
        if self.ms_data is None:
            print("No MS data loaded.")
            return pd.DataFrame()
        
        # Parse replicate information
        replicate_info = self.parse_replicate_info(self.ms_data['replicate'])
        
        # Combine with main data
        ms_with_conditions = pd.concat([self.ms_data, replicate_info[['day', 'condition', 'wash_type']]], axis=1)
        
        # Filter for precursor ions (main peptide signal)
        precursor_data = ms_with_conditions[ms_with_conditions['fragment ion'] == 'precursor'].copy()
        
        # Aggregate by peptide, day, and condition
        aggregated = precursor_data.groupby(['peptide', 'protein', 'day', 'condition', 'wash_type']).agg({
            'area': ['mean', 'std', 'count', 'sum']
        }).reset_index()
        
        # Flatten column names
        aggregated.columns = ['_'.join(col).strip('_') for col in aggregated.columns.values]
        aggregated.columns = [col.replace('_', '') if col.endswith('_') else col for col in aggregated.columns]
        
        self.processed_data['peptide_areas_by_condition'] = aggregated
        return aggregated
    
    def analyze_time_course(self) -> pd.DataFrame:
        """Analyze peptide concentration time course data"""
        if self.peptide_conc_data is None:
            print("No peptide concentration data loaded.")
            return pd.DataFrame()
        
        # Melt the dataframe to long format
        conc_long = pd.melt(
            self.peptide_conc_data, 
            id_vars=['Peptides'], 
            var_name='timepoint', 
            value_name='concentration_ng_mL'
        )
        
        # Extract day information
        conc_long['day'] = conc_long['timepoint'].str.extract(r'D(\d+)')[0]
        conc_long['day'] = pd.to_numeric(conc_long['day'], errors='coerce')
        
        # Calculate fold changes relative to D0 (if D0 has non-zero values)
        pivot_conc = conc_long.pivot(index='Peptides', columns='day', values='concentration_ng_mL')
        
        # Calculate fold changes (relative to D1 since D0 is zero)
        fold_changes = pivot_conc.div(pivot_conc.iloc[:, 0], axis=0)  # Relative to first non-zero timepoint
        fold_changes = fold_changes.replace([np.inf, -np.inf], np.nan)
        
        self.processed_data['concentration_timecourse'] = conc_long
        self.processed_data['fold_changes'] = fold_changes
        
        return conc_long
    
    def create_summary_plots(self, save_plots: bool = True):
        """Create summary plots of the JPT data"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('JPT Data Analysis Summary', fontsize=16, fontweight='bold')
        
        # Plot 1: Fragment ion distribution (MS data)
        if self.ms_data is not None:
            fragment_counts = self.ms_data['fragment ion'].value_counts()
            axes[0, 0].bar(range(len(fragment_counts)), fragment_counts.values)
            axes[0, 0].set_xticks(range(len(fragment_counts)))
            axes[0, 0].set_xticklabels(fragment_counts.index, rotation=45, ha='right')
            axes[0, 0].set_title('Fragment Ion Distribution')
            axes[0, 0].set_ylabel('Count')
        
        # Plot 2: Area distribution (MS data)
        if self.ms_data is not None:
            # Log transform for better visualization
            log_areas = np.log10(self.ms_data['area'] + 1)
            axes[0, 1].hist(log_areas, bins=50, alpha=0.7, edgecolor='black')
            axes[0, 1].set_title('Area Distribution (Log10)')
            axes[0, 1].set_xlabel('Log10(Area + 1)')
            axes[0, 1].set_ylabel('Frequency')
        
        # Plot 3: Top peptides by measurement count
        if self.ms_data is not None:
            top_peptides = self.ms_data['peptide'].value_counts().head(10)
            axes[0, 2].barh(range(len(top_peptides)), top_peptides.values)
            axes[0, 2].set_yticks(range(len(top_peptides)))
            axes[0, 2].set_yticklabels(top_peptides.index)
            axes[0, 2].set_title('Top 10 Peptides by Measurement Count')
            axes[0, 2].set_xlabel('Count')
        
        # Plot 4: Peptide concentration heatmap
        if self.peptide_conc_data is not None:
            conc_matrix = self.peptide_conc_data.set_index('Peptides').iloc[:10]  # Top 10 peptides
            im = axes[1, 0].imshow(conc_matrix.values, aspect='auto', cmap='viridis')
            axes[1, 0].set_xticks(range(len(conc_matrix.columns)))
            axes[1, 0].set_xticklabels(conc_matrix.columns, rotation=45)
            axes[1, 0].set_yticks(range(len(conc_matrix.index)))
            axes[1, 0].set_yticklabels(conc_matrix.index)
            axes[1, 0].set_title('Peptide Concentrations (Top 10)')
            plt.colorbar(im, ax=axes[1, 0], label='ng/mL')
        
        # Plot 5: Time course for selected peptides
        if 'concentration_timecourse' in self.processed_data:
            conc_data = self.processed_data['concentration_timecourse']
            top_peptides_conc = conc_data['Peptides'].value_counts().head(5).index
            
            for i, peptide in enumerate(top_peptides_conc):
                peptide_data = conc_data[conc_data['Peptides'] == peptide]
                axes[1, 1].plot(peptide_data['day'], peptide_data['concentration_ng_mL'], 
                               marker='o', label=peptide[:10] + '...' if len(peptide) > 10 else peptide)
            
            axes[1, 1].set_title('Concentration Time Course (Top 5 Peptides)')
            axes[1, 1].set_xlabel('Day')
            axes[1, 1].set_ylabel('Concentration (ng/mL)')
            axes[1, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            axes[1, 1].set_yscale('log')
        
        # Plot 6: Summary statistics
        if self.ms_data is not None:
            stats_text = f"""
            Total Measurements: {len(self.ms_data):,}
            Unique Peptides: {self.ms_data['peptide'].nunique():,}
            Unique Proteins: {self.ms_data['protein'].nunique():,}
            Unique Replicates: {self.ms_data['replicate'].nunique():,}
            
            Area Statistics:
            Mean: {self.ms_data['area'].mean():.2e}
            Median: {self.ms_data['area'].median():.2e}
            Std: {self.ms_data['area'].std():.2e}
            """
            axes[1, 2].text(0.1, 0.9, stats_text, transform=axes[1, 2].transAxes, 
                            fontsize=10, verticalalignment='top', fontfamily='monospace')
            axes[1, 2].set_title('Data Summary Statistics')
            axes[1, 2].axis('off')
        
        plt.tight_layout()
        
        if save_plots:
            plot_path = self.output_dir / "jpt_data_summary.png"
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            print(f"Summary plots saved to: {plot_path}")
        
        plt.show()
    
    def export_processed_data(self):
        """Export all processed data to CSV files"""
        print("Exporting processed data...")
        
        # Export processed peptide areas
        if 'peptide_areas_by_condition' in self.processed_data:
            output_path = self.output_dir / "jpt_peptide_areas_processed.csv"
            self.processed_data['peptide_areas_by_condition'].to_csv(output_path, index=False)
            print(f"Peptide areas by condition saved to: {output_path}")
        
        # Export concentration timecourse
        if 'concentration_timecourse' in self.processed_data:
            output_path = self.output_dir / "jpt_concentration_timecourse.csv"
            self.processed_data['concentration_timecourse'].to_csv(output_path, index=False)
            print(f"Concentration timecourse saved to: {output_path}")
        
        # Export fold changes
        if 'fold_changes' in self.processed_data:
            output_path = self.output_dir / "jpt_concentration_fold_changes.csv"
            self.processed_data['fold_changes'].to_csv(output_path)
            print(f"Concentration fold changes saved to: {output_path}")
        
        # Export summary statistics
        if self.ms_data is not None:
            summary = self.summarize_ms_data()
            
            # Create a summary dataframe
            summary_data = []
            for key, value in summary.items():
                if isinstance(value, dict):
                    for subkey, subvalue in value.items():
                        summary_data.append({'category': key, 'metric': subkey, 'value': subvalue})
                else:
                    summary_data.append({'category': 'basic_stats', 'metric': key, 'value': value})
            
            summary_df = pd.DataFrame(summary_data)
            output_path = self.output_dir / "jpt_summary_statistics.csv"
            summary_df.to_csv(output_path, index=False)
            print(f"Summary statistics saved to: {output_path}")
    
    def run_full_analysis(self):
        """Run complete JPT data analysis pipeline"""
        print("="*60)
        print("JPT DATA ANALYSIS PIPELINE")
        print("="*60)
        
        # Load data
        print("\n1. Loading data files...")
        self.load_ms_data()
        self.load_peptide_concentrations()
        
        # Process data
        print("\n2. Processing mass spectrometry data...")
        if self.ms_data is not None:
            peptide_areas = self.process_peptide_areas_by_condition()
            print(f"Processed {len(peptide_areas)} peptide-condition combinations")
        
        print("\n3. Analyzing concentration time course...")
        if self.peptide_conc_data is not None:
            conc_timecourse = self.analyze_time_course()
            print(f"Analyzed {len(conc_timecourse)} concentration measurements")
        
        # Generate summaries
        print("\n4. Generating summary statistics...")
        if self.ms_data is not None:
            summary = self.summarize_ms_data()
            print("Summary statistics generated")
        
        # Create plots
        print("\n5. Creating summary plots...")
        self.create_summary_plots()
        
        # Export results
        print("\n6. Exporting processed data...")
        self.export_processed_data()
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print("="*60)
        print(f"Results saved to: {self.output_dir.absolute()}")


def main():
    """Main function to run JPT data processing"""
    print("JPT Files Processing Script")
    print("Processing mass spectrometry and peptide concentration data...")
    
    # Initialize processor
    processor = JPTDataProcessor()
    
    # Run full analysis
    processor.run_full_analysis()
    
    return processor


if __name__ == "__main__":
    # Run the analysis
    processor = main()
    
    # Optional: Interactive analysis
    print("\nProcessor object available as 'processor' for interactive analysis")
    print("Available methods:")
    print("- processor.ms_data (DataFrame)")
    print("- processor.peptide_conc_data (DataFrame)")
    print("- processor.processed_data (Dict of processed DataFrames)")
    print("- processor.summarize_ms_data()")
    print("- processor.create_summary_plots()")