# Unified CLI Guide

## Overview

The unified CLI automatically detects your input format and processes accordingly. **No need to specify format or use separate workflows!**

## Basic Usage

```bash
python main.py <ms_data.csv> <concentrations.csv> -o <output.csv>
```

That's it! The system will:
1. 🔍 Inspect your input files
2. 🎯 Detect the format automatically
3. 🔧 Route to the appropriate processing pipeline
4. 💾 Save results with proper structure

## Examples

### Heavy/Light Paired Peptides
```bash
python main.py \
  data/input/heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv \
  data/input/peptide_dilution_conc_peggy.csv \
  -o results_heavy.csv
```

**Detected as:** `PAIRED RATIO` mode
- Calculates area ratios (heavy/light)
- Fits calibration curves per fragment ion

### JPT Format (also has pairs)
```bash
python main.py \
  data/input/JPT_1_3_-_HHT_1.2_-_Cytosolic.csv \
  data/input/JPT1-3_peptide_conc.csv \
  -o results_jpt.csv
```

**Detected as:** `PAIRED RATIO` mode (JPT also uses isotope pairs!)
- Same processing as heavy_1st
- Handles multiple experimental conditions

### Future Formats
Any new format with the same column structure will work automatically!

## What Gets Detected

The system inspects:

1. **Light/Heavy Pairs**: Multiple precursor m/z values per peptide/replicate/fragment
   - ✓ Found → Use paired ratio analysis
   - ✗ Not found → Use single intensity analysis

2. **Experimental Conditions**: Replicate naming patterns
   - `cyto`, `SDC_Sup`, `W_`, `WW_` → Multiple conditions detected
   - Affects aggregation strategy

3. **Concentration Columns**: D0, D1, D2, etc.
   - Automatically maps dilution points
   - Works with any number of dilution levels

## Input File Requirements

### MS Data File (Skyline Export)
**Required columns:**
- `Peptide`
- `Protein`
- `Replicate`
- `Precursor Mz`
- `Precursor Charge`
- `Product Mz`
- `Product Charge`
- `Fragment Ion`
- `Area`

### Concentration File
**Required columns:**
- `Peptides`: Peptide sequences (must match MS file)
- `D0 (ng/mL)`, `D1 (ng/mL)`, ...: Concentration values

**Note:** Column order and exact number of dilution points don't matter!

## Processing Modes

### Mode 1: PAIRED RATIO
**Triggered when:** Light/heavy pairs detected

**Processing:**
1. Groups data by peptide/replicate/fragment
2. Calculates area ratios (heavy/light)
3. Merges with concentration data
4. Fits linear regression: `area_ratio ~ concentration`
5. Computes QC metrics (R², CV, Q-tests)

**Output columns:**
- `area_ratio`, `area_min`, `area_max`
- `R2`, `intercept`, `gradient`
- `mean_r2`, `std_r2`, `cov_r2`
- `mean_grad`, `std_grad`, `cov_grad`
- `mean_intercept`, `std_intercept`, `cov_intercept`
- Q-test values

### Mode 2: SINGLE INTENSITY (with conditions)
**Triggered when:** No pairs, but conditions detected

**Processing:**
1. Extracts conditions from replicate names
2. Aggregates intensities by peptide/condition/fragment
3. Fits regression: `intensity ~ concentration` per condition
4. Computes aggregated QC metrics

**Output columns:**
- `area_mean`, `area_std`, `area_count`
- `concentration_ng_mL`
- `R2`, `slope`, `intercept`
- `mean_r2`, `std_r2`, `cv_r2`
- Similar aggregations for slope/intercept

### Mode 3: SINGLE BASIC
**Triggered when:** No pairs, no conditions

**Processing:**
- Simplified single-peptide analysis
- Direct intensity vs concentration regression

## Command Line Options

```
python main.py -h
```

**Arguments:**
- `ms_file`: Path to mass spectrometry CSV (required)
- `concentration_file`: Path to concentration CSV (required)
- `-o, --output`: Output file path (default: `prm_analysis_output.csv`)
- `--version`: Show version

## Adding New Formats

To support a new format variation:

1. **Same column structure?** → It already works! Just run it.

2. **New detection pattern needed?**
   - Edit `scripts/unified_processor.py`
   - Update `_detect_format()` method
   - Add new condition checks

3. **New processing logic needed?**
   - Create new processor in `scripts/process_<format>_data.py`
   - Add new mode in `PRMDataProcessor.process()`
   - Update `get_processing_mode()`

## Troubleshooting

### "Format detected incorrectly"
Check your data structure:
```bash
head -5 your_ms_file.csv
head -5 your_concentration_file.csv
```

Ensure:
- Column names match expected format (case-insensitive)
- Peptide names are consistent between files
- Replicate names follow patterns if using conditions

### "No results generated"
Common causes:
1. Peptide names don't match between MS and concentration files
2. Missing concentration values (NaN in concentration file)
3. Insufficient data points for regression (< 2 points)

### "Processing mode unexpected"
The detection logic prioritizes:
1. Light/heavy pairs (highest priority)
2. Multiple conditions (if no pairs)
3. Basic single-peptide (fallback)

If you want to force a specific behavior, you can modify the concentration file or replicate naming to trigger the desired mode.

## Performance

- **Heavy_1st**: ~624 rows → 220 results in ~2 seconds
- **JPT**: ~9,398 rows → 2,869 results in ~5 seconds
- Scales linearly with input size

## Future Enhancements

Potential additions:
- [ ] Config file support for custom detection rules
- [ ] Batch processing multiple file pairs
- [ ] Integration with plotting tools
- [ ] Web interface for parameter selection
- [ ] Support for non-Skyline export formats
