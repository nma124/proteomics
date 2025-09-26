# Documentation

This directory contains detailed documentation for the Proteomics PRM Data Processing toolkit.

## Contents

- **API Reference**: Detailed function and class documentation
- **User Guide**: Step-by-step tutorials and examples
- **Developer Guide**: Information for contributors and developers
- **FAQ**: Frequently asked questions and troubleshooting

## Quick Links

- [Installation Guide](installation.md)
- [Getting Started](getting_started.md)
- [API Documentation](api.md)
- [Examples and Tutorials](examples.md)

## Building Documentation

To build the documentation locally:

```bash
# Install documentation dependencies
pip install -e ".[docs]"

# Build HTML documentation
cd docs
sphinx-build -b html . _build
```

The documentation will be available in `_build/index.html`.