# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

This is a bioinformatics research project focused on analyzing body-wide organ network transcriptome data across different mouse strains. The project implements likelihood-based statistical methods to compare organ-organ and strain-strain relationships in biological data.

## Development Environment Setup

### Dependencies
```bash
pip install -r requirements.txt
```

The project uses specific versions of core scientific computing libraries:
- matplotlib==3.3.4
- numpy==1.20.2
- pandas==1.2.4

### Running Jupyter Notebooks
```bash
jupyter lab notebooks/
```

### Python Package Structure
The main code is in `src/utils.py` and can be imported as:
```python
from src import utils
```

## Key Architecture Components

### Data Processing Pipeline
The codebase implements a specialized likelihood computation framework for proteomics data:

1. **Data Label Standardization**: Functions like `complete_label()` and `fix_time_label()` standardize column/row naming conventions following the pattern `organ_strain_timepoint_condition`

2. **Variable Extraction System**: The code uses positional indexing to extract variables from underscore-separated labels:
   - Index 0: Organ names (eye, kidney, liver, lung, muscle, pancreas, small.intestine, spleen)
   - Index 1: Strain names (AJ, BL, NOD, SJL)  
   - Index 2: Time points (with '0h' as default)

3. **Likelihood Analysis Engine**: Two main analytical functions:
   - `get_strain_var_likelihood()`: Computes strain-strain likelihood matrices
   - `get_organ_var_likelihood()`: Computes organ-organ likelihood matrices

### Core Statistical Methods

The likelihood computation uses Gaussian probability distributions:
- Reference distribution parameters (μ, σ) calculated from one variable pair
- Test data likelihood computed against reference distribution
- Geometric normalization applied: `L^(1/n)` where n is sample size
- Results aggregated across all combinations

### Data Structure Conventions

Input data format: Correlation matrices with labels following `organ_strain_timepoint` pattern
- Symmetric matrices with matching row/column labels
- NaN values in upper triangles (lower triangular analysis)
- Data filtering supports common variable selection across groups

## Common Development Tasks

### Running Analysis Notebooks
```bash
# Navigate to notebooks directory and start analysis
cd notebooks/
jupyter lab 006-organ-strain-likelihood.ipynb
```

### Working with the Utils Module
The main analysis functions require these parameters:
- `var_str_idx`: Index position for the primary variable (0 for organs, 1 for strains)
- `var_2_idx`: Index position for the conditioning variable  
- `data`: Input correlation matrix
- `geom_normalize=True`: Apply geometric normalization
- `use_common_var_2=True`: Filter to common variables across groups

### Data File Locations
- Raw data expected in `data/raw/` (gitignored)
- Notebook outputs in `notebooks/`
- Generated reports in `reports/`

## Debugging and Development Notes

- The codebase includes extensive debug modes in likelihood functions
- Set `debug=True` in likelihood functions to see intermediate calculations
- Specific test cases are hardcoded for validation (e.g., AJ-BL strain pairs, bone-brain organ pairs)
- The ranking system uses Excel output with conditional formatting for result visualization
