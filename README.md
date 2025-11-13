# Proteomics PRM Data Processing

[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Automated processing toolkit for Skyline PRM (Parallel Reaction Monitoring) export data with intelligent format detection and comprehensive statistical analysis.

## 🚀 Quick Start

```bash
# Clone and install
git clone https://github.com/nma124/proteomics.git
cd proteomics
pip install -r requirements.txt

# Process your data (auto-detects format!)
python main.py <ms_data.csv> <concentrations.csv> -o results.csv
```

### What it does:
- ✅ Automatically detects input format (Heavy/Light paired, JPT, single peptide)
- ✅ Calculates area ratios and performs linear regression analysis
- ✅ Generates quality control metrics (R², CV, Q-tests)
- ✅ Outputs processed data ready for downstream analysis

## 📊 Example Results

Processing a heavy peptide dilution series dataset:
- **Input:** 624 MS measurements across 8 dilution levels (D0-D7)
- **Output:** 220 rows with regression analysis and QC metrics
- **Analysis:** Area ratio calculations, linear regression (R² = 0.692), statistical validation
- **QC Metrics:** Coefficient of variation, Q-tests, replicate consistency checks

## 📁 Project Structure

```
proteomics/
├── main.py              # ➡️ CLI entry point
├── requirements.txt     # Python dependencies
├── setup.py             # Package installation
├── scripts/             # Core processing modules
│   ├── unified_processor.py   # Format detection & routing
│   ├── process_prm_data.py    # Heavy/Light processing
│   └── process_jpt_data.py    # JPT format processing
├── src/                 # Statistical analysis utilities
├── data/                # Data directory (gitignored)
├── docs/                # Documentation
├── examples/            # Example workflows
└── tests/               # Unit tests
```

## 🔧 Installation

```bash
# Clone repository
git clone https://github.com/nma124/proteomics.git
cd proteomics

# Install dependencies
pip install -r requirements.txt

# Optional: Install as package
pip install -e .
```

**Requirements:** Python 3.9+, pandas, numpy, scikit-learn

## 💻 Usage

### Basic Usage

```bash
python main.py <ms_data.csv> <concentrations.csv> -o <output.csv>
```

**Arguments:**
- `ms_data.csv` - Skyline PRM export file
- `concentrations.csv` - Peptide dilution concentrations
- `-o, --output` - Output file path (optional, default: `prm_analysis_output.csv`)

### Examples

```bash
# Process heavy/light paired peptides
python main.py heavy_1st_data.csv peptide_conc.csv -o results.csv

# Process JPT format
python main.py jpt_data.csv jpt_concentrations.csv -o jpt_results.csv

# Get help
python main.py --help
```

### Python API

```python
from scripts.unified_processor import process_prm_unified

# Process data programmatically
result_df = process_prm_unified(
    ms_file="ms_data.csv",
    concentration_file="concentrations.csv",
    output_file="results.csv"
)

print(f"Processed {result_df.shape[0]} rows")
```

## 📄 Input Data Format

### MS Data File (Skyline PRM Export)
Required columns: `Peptide`, `Protein`, `Replicate`, `Precursor Mz`, `Precursor Charge`, `Product Mz`, `Product Charge`, `Fragment Ion`, `Area`

### Concentration File
Required format: `Peptides` column + dilution columns (`D0 (ng/mL)`, `D1 (ng/mL)`, ..., `D7 (ng/mL)`)

## 📈 Output

Processed CSV file containing:
- Area ratios (light/heavy peptides)
- Linear regression metrics (R², slope, intercept)
- Quality control statistics (CV, Q-tests)
- Aggregated statistics across replicates

## 🔬 Analysis Pipeline

1. **Format Detection** - Automatically identifies input format
2. **Data Processing** - Calculates area ratios, filters outliers
3. **Regression Analysis** - Fits linear models for concentration curves
4. **QC Metrics** - Generates quality control statistics
5. **Output Generation** - Saves processed data with all metrics

## 🔗 Additional Resources

- **Detailed CLI Guide**: See [docs/UNIFIED_CLI_GUIDE.md](docs/UNIFIED_CLI_GUIDE.md)
- **Architecture Details**: See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)  
- **Web Interface**: A web-based interface is currently under development for users who prefer a graphical interface

## 📝 License

MIT License - See [LICENSE](LICENSE) for details.

## 👥 Contributing

This tool is designed for processing Skyline PRM export data in proteomics workflows. Contributions and issues are welcome via GitHub.

---

**Repository**: https://github.com/nma124/proteomics
