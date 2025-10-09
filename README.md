# Proteomics PRM Data Processing

[![Python](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A comprehensive toolkit for processing Skyline PRM (Parallel Reaction Monitoring) export data with heavy peptide dilution schemes. Performs quantification analysis with regression modeling and quality control metrics.

## 🚀 Quick Start

```bash
# Clone the repository
git clone <repository-url>
cd proteomics

# Install dependencies
pip install -r requirements.txt

# Run the example workflow
python examples/heavy_1st_workflow.py

# Or use the command-line interface
python main.py data/input/skyline_data.csv data/input/dilutions.csv -o results.csv
```

## ✅ Recent Success: Heavy_1st Dataset Analysis

**Input:** `heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv`  
**Output:** `heavy_1st_processed_output.csv` (105 KB)

### Key Results:
- **📊 Data processed:** 624 rows → 220 rows (after filtering and analysis)
- **🧬 Peptide analyzed:** AHSLNPDASGSSCSLAR (single peptide)
- **📈 Fragment ions:** 11 different fragment ions
- **🔬 Dilution series:** 8 concentration levels (D0-D7)
- **🔄 Replicates:** 24 different replicate conditions
- **📏 Average R²:** 0.692 (good linear correlation)

### What the Analysis Did:
1. **Area Ratio Calculation:** Computed ratios between light and heavy peptide peak areas
2. **Dilution Series Analysis:** Merged data with concentration information
3. **Linear Regression:** Fitted regression models for each fragment ion
4. **Quality Metrics:** Generated R², gradients, intercepts, and statistical measures
5. **Aggregation:** Calculated mean, standard deviation, and coefficient of variation for regression parameters

### Output File Contains:
- Area ratios for each dilution/fragment ion combination
- Linear regression parameters (R², slope, intercept)
- Quality control metrics (Q-tests, coefficients of variation)
- Aggregated statistics across replicates

The processed data is ready for downstream analysis, visualization, or quantification workflows. The average R² of 0.692 indicates reasonably good linear relationships between area ratios and heavy peptide concentrations.

## 📁 Project Structure

```
proteomics/
├── 📜 README.md                 # This file
├── 📜 requirements.txt          # Python dependencies
├── 📜 setup.py                 # Installation script
├── 📜 main.py                  # Command-line interface
├── 📂 scripts/                 # Core processing modules
│   ├── 📜 __init__.py
│   ├── 📜 process_prm_data.py   # Main processing functions
│   └── 📜 run_heavy_1st.py     # Heavy_1st specific runner
├── 📂 data/
│   ├── 📂 input/               # Input data files
│   │   ├── 📄 heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv
│   │   └── 📄 peptide_dilution_conc_peggy.csv
│   ├── 📂 output/              # Processed results
│   │   ├── 📄 heavy_1st_processed_output.csv
│   │   └── 📄 Akin_PRM_heavy_di_expanded_AHS_3+only_rmHL0_output.csv
│   └── 📂 raw/                 # Original raw data
├── 📂 examples/                # Example workflows and tutorials
│   ├── 📜 __init__.py
│   └── 📜 heavy_1st_workflow.py # Complete workflow demonstration
├── 📂 docs/                    # Documentation
├── 📂 tests/                   # Unit tests
├── 📂 config/                  # Configuration files
└── 📂 notebooks/               # Jupyter analysis notebooks
```

## 🔧 Installation

### Method 1: Direct Installation
```bash
pip install -r requirements.txt
```

### Method 2: Development Installation
```bash
pip install -e .
```

### Method 3: With Optional Dependencies
```bash
# Install with development tools
pip install -e ".[dev]"

# Install with documentation tools
pip install -e ".[docs]"
```

## 💻 Usage

### Command Line Interface

```bash
# Basic usage
python main.py skyline_data.csv dilution_data.csv

# Specify output file
python main.py skyline_data.csv dilution_data.csv -o my_results.csv

# Get help
python main.py --help
```

### Python API

```python
from scripts.process_prm_data import process_prm_data

# Process your data
result_df = process_prm_data(
    skyline_file="data/input/skyline_export.csv",
    dilution_file="data/input/dilution_concentrations.csv",
    output_file="results.csv"
)

# Access the processed DataFrame
print(f"Processed {result_df.shape[0]} rows")
print(f"Average R²: {result_df['mean_r2'].mean():.3f}")
```

### Example Workflows

```bash
# Run the complete heavy_1st workflow
python examples/heavy_1st_workflow.py

# Or as a module
python -m examples.heavy_1st_workflow
```

## 📊 Input Data Format

### Skyline Export File
Expected columns:
- `Peptide`: Peptide sequence
- `Protein`: Protein identifier
- `Replicate`: Sample replicate name (with dilution info)
- `Precursor Mz`: Precursor mass-to-charge ratio
- `Precursor Charge`: Precursor charge state
- `Product Mz`: Fragment ion mass-to-charge ratio
- `Product Charge`: Fragment ion charge state
- `Fragment Ion`: Fragment ion identifier
- `Area`: Peak area

### Dilution Concentration File
Expected format:
- `Peptides`: Peptide sequences (matching Skyline export)
- `D1 (ng/mL)`, `D2 (ng/mL)`, ..., `D7 (ng/mL)`: Dilution concentrations
- `D0 (ng/mL)`: Control/blank concentration

## 📈 Output Data

The processed output contains:

- **Original data**: Peptide, replicate, fragment ion information
- **Area ratios**: Calculated ratios between light and heavy peptides
- **Regression parameters**: R², slope, intercept for each fragment ion
- **Quality metrics**: Statistical measures for regression quality
- **Aggregated statistics**: Mean, std dev, coefficient of variation across replicates

## 🧪 Analysis Pipeline

1. **Data Loading**: Read Skyline export and dilution concentration files
2. **Data Filtering**: Keep only groups with exactly 2 elements (light/heavy pairs)
3. **Area Ratio Calculation**: Compute heavy/light area ratios
4. **Concentration Mapping**: Merge with dilution concentrations
5. **Linear Regression**: Fit models for area_ratio vs concentration
6. **Quality Assessment**: Calculate R², Q-tests, coefficients of variation
7. **Aggregation**: Summarize statistics across technical replicates

## 🔍 Quality Control

The pipeline includes several QC metrics:

- **R² values**: Measure of linear fit quality
- **Q-tests**: Outlier detection for regression parameters
- **Coefficients of variation**: Measure of replicate consistency
- **Data completeness**: Tracking of filtered vs retained data points

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Note**: This tool was developed for processing Skyline PRM export data in proteomics workflows. It has been tested with the heavy_1st dataset and similar experimental designs involving heavy peptide dilution series.

# proteomics
Analysis of body-wide organ network transcriptome for different mouse strains
