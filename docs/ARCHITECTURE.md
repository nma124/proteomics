# Unified Processor Architecture

## Overview

The unified processor uses a **detect-and-route** pattern to handle multiple input formats with a single CLI.

## Architecture Flow

```
┌─────────────────────────────────────────────────────────────┐
│                       USER INPUT                            │
│   python main.py ms_data.csv concentrations.csv -o out.csv  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                   UNIFIED PROCESSOR                          │
│              (scripts/unified_processor.py)                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
           ┌─────────────────────┐
           │  FORMAT DETECTION   │
           │                     │
           │  1. Load files      │
           │  2. Inspect columns │
           │  3. Check patterns  │
           └──────────┬──────────┘
                      │
         ┌────────────┼────────────┐
         │            │            │
         ▼            ▼            ▼
    ┌────────┐  ┌─────────┐  ┌──────────┐
    │ Pairs? │  │ Condi-  │  │ Dilution │
    │ Check  │  │ tions?  │  │ Columns? │
    └────┬───┘  └────┬────┘  └────┬─────┘
         │           │            │
         └───────────┼────────────┘
                     │
         ┌───────────┴───────────┐
         │                       │
         ▼                       ▼
    HAS PAIRS              NO PAIRS
         │                       │
         │                  ┌────┴─────┐
         │                  │          │
         │                  ▼          ▼
         │            HAS COND    NO COND
         │                  │          │
         ▼                  ▼          ▼
┌──────────────┐   ┌──────────────┐  ┌──────────────┐
│ PAIRED_RATIO │   │SINGLE_INTENS.│  │ SINGLE_BASIC │
│   MODE       │   │    MODE      │  │    MODE      │
└──────┬───────┘   └──────┬───────┘  └──────┬───────┘
       │                  │                 │
       ▼                  ▼                 ▼
┌─────────────┐   ┌──────────────┐  ┌──────────────┐
│ Heavy/Light │   │ JPT Format   │  │ Basic Single │
│ Processor   │   │ Processor    │  │ Processor    │
└──────┬──────┘   └──────┬───────┘  └──────┬───────┘
       │                  │                 │
       └──────────────────┼─────────────────┘
                          │
                          ▼
               ┌────────────────────┐
               │  REGRESSION ENGINE │
               │                    │
               │ • Area ratio calc  │
               │ • Linear fit       │
               │ • QC metrics       │
               │ • Aggregation      │
               └─────────┬──────────┘
                         │
                         ▼
                 ┌──────────────┐
                 │ OUTPUT FILE  │
                 │   (.csv)     │
                 └──────────────┘
```

## Detection Logic

### Format Features Detected

| Feature | Detection Method | Impact |
|---------|-----------------|---------|
| **Light/Heavy Pairs** | Multiple precursor m/z per peptide/replicate/fragment | → PAIRED_RATIO mode |
| **Experimental Conditions** | Keywords in replicate names (`cyto`, `SDC_Sup`, `W_`, `WW_`) | → Aggregation by condition |
| **Dilution Columns** | Columns starting with `D` in concentration file | → Auto-map concentrations |

### Decision Tree

```
IF has_light_heavy_pairs:
    mode = "paired_ratio"
    → Use Heavy/Light processor
    
ELIF has_multiple_conditions:
    mode = "single_intensity"
    → Use JPT processor (with conditions)
    
ELSE:
    mode = "single_basic"
    → Use simplified single-peptide processor
```

## Module Structure

```
proteomics/
├── main.py                          # CLI entry point
├── scripts/
│   ├── unified_processor.py         # 🎯 Format detection & routing
│   ├── process_prm_data.py          # Heavy/Light processing
│   ├── process_jpt_data.py          # JPT/condition processing
│   └── [future processors]          # Easy to add!
└── docs/
    ├── UNIFIED_CLI_GUIDE.md         # User guide
    └── ARCHITECTURE.md              # This file
```

## Key Classes

### `PRMDataProcessor`
**Location:** `scripts/unified_processor.py`

**Responsibilities:**
- Load and inspect input files
- Detect format characteristics
- Route to appropriate processor
- Provide format summary

**Methods:**
```python
__init__(ms_file, concentration_file)  # Load data
_detect_format() → dict                # Inspect structure
get_processing_mode() → str            # Determine mode
process(output_file) → DataFrame       # Execute pipeline
print_format_summary()                 # Display detection results
```

## Adding New Formats

### Step 1: Add Detection Pattern
Edit `scripts/unified_processor.py`:

```python
def _detect_format(self) -> dict:
    # ... existing checks ...
    
    # Add new pattern check
    if self._check_for_new_pattern():
        format_info['has_new_feature'] = True
```

### Step 2: Update Mode Selection
```python
def get_processing_mode(self) -> str:
    if self.has_pairs:
        return 'paired_ratio'
    elif self.has_new_feature:  # NEW
        return 'new_format_mode'
    # ... rest of logic ...
```

### Step 3: Add Processor
Create `scripts/process_new_format_data.py`:

```python
def process_new_format_data(ms_file, conc_file, output_file):
    # Your processing logic
    return result_df
```

### Step 4: Route to Processor
```python
def _process_new_format(self, output_file: str) -> pd.DataFrame:
    from scripts.process_new_format_data import process_new_format_data
    return process_new_format_data(self.ms_file, self.concentration_file, output_file)
```

## Benefits of This Architecture

### ✅ Extensibility
- Add new formats without changing CLI
- Processors are independent modules
- Detection logic is centralized

### ✅ User-Friendly
- Single command for all formats
- No manual format specification
- Clear feedback about detection

### ✅ Maintainability
- Each processor handles one format
- Easy to test individual components
- Clear separation of concerns

### ✅ Robustness
- Fallback modes for edge cases
- Detailed format inspection
- Informative error messages

## Processing Modes Explained

### PAIRED_RATIO Mode
**Input:** Light and heavy isotope pairs  
**Output:** Area ratios with calibration curves

**Pipeline:**
1. Group by peptide/replicate/fragment
2. Filter for pairs (n=2 per group)
3. Calculate area ratios (heavy/light)
4. Merge with concentrations
5. Fit linear regression per fragment
6. Compute QC metrics
7. Aggregate across replicates

**Used by:** Heavy_1st, JPT (both have isotope pairs!)

### SINGLE_INTENSITY Mode
**Input:** Single peptide measurements with conditions  
**Output:** Intensity-based calibration per condition

**Pipeline:**
1. Extract conditions from replicate names
2. Aggregate areas by peptide/condition/fragment
3. Merge with concentrations
4. Fit regression per condition
5. Compute aggregated QC metrics

**Used by:** (Future formats without pairs but with conditions)

### SINGLE_BASIC Mode
**Input:** Single peptide measurements, no conditions  
**Output:** Simple intensity calibration

**Pipeline:**
1. Aggregate intensities by peptide/fragment
2. Merge with concentrations
3. Fit simple regression
4. Compute basic QC metrics

**Used by:** (Fallback for minimal formats)

## Future Enhancements

### Planned Features
- [ ] Config file for custom detection rules
- [ ] Plugin system for external processors
- [ ] Format validation and suggestions
- [ ] Batch processing mode
- [ ] Format conversion utilities

### Potential New Formats
- Label-free quantification (LFQ)
- SILAC (stable isotope labeling)
- TMT/iTRAQ (isobaric tags)
- Data-independent acquisition (DIA)

Each can be added without changing the core architecture!
