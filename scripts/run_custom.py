#!/usr/bin/env python
import argparse
import sys
from pathlib import Path
import pandas as pd

try:
    from scripts.process_prm_data import process_prm_data
except Exception as e:
    sys.stderr.write(
        "ERROR: Could not import scripts.process_prm_data.\n"
        "Run this from your repo root so 'scripts' is importable, "
        "and ensure scripts/process_prm_data.py exists.\n"
        f"Underlying import error: {e}\n"
    )
    sys.exit(1)

def to_csv_if_needed(path_str: str, sheet, tmp_dir: Path) -> str:
    p = Path(path_str)
    if not p.exists():
        raise FileNotFoundError(f"Input not found: {p}")
    if p.suffix.lower() == ".csv":
        return str(p)
    if p.suffix.lower() in (".xlsx", ".xls"):
        # Convert Excel -> CSV (first sheet by default)
        if sheet is None:
            sheet = 0
        try:
            df = pd.read_excel(p, sheet_name=sheet)
        except Exception as e:
            raise RuntimeError(
                f"Failed reading Excel file: {p}. "
                f"If it has multiple sheets, pass --skyline-sheet/--dilutions-sheet.\nOriginal error: {e}"
            )
        out_name = p.stem.replace(" ", "_") + ".csv"
        out_path = tmp_dir / out_name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path, index=False)
        print(f"[convert] Wrote temporary CSV: {out_path}")
        return str(out_path)
    raise ValueError(f"Unsupported file extension for {p.name}. Use .csv or .xlsx/.xls")

def parse_sheet(val):
    if val is None:
        return None
    try:
        return int(val)
    except ValueError:
        return val  # treat as sheet name

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Run PRM pipeline on Skyline + dilution files (CSV or Excel)."
    )
    parser.add_argument("--skyline", required=True, help="Skyline export path (.csv or .xlsx/.xls)")
    parser.add_argument("--dilutions", required=True, help="Dilution-concentration table (.csv or .xlsx/.xls)")
    parser.add_argument("-o", "--out", default=None, help="Output CSV path (default: <skyline_basename>__processed.csv)")
    parser.add_argument("--skyline-sheet", default=None, help="Excel sheet (name or 0-based index) for Skyline file")
    parser.add_argument("--dilutions-sheet", default=None, help="Excel sheet (name or 0-based index) for dilutions file")
    parser.add_argument("--tmp-dir", default="data/tmp", help="Folder for temporary CSVs when converting Excel")

    args = parser.parse_args(argv)

    skyline_sheet = parse_sheet(args.skyline_sheet)
    dilutions_sheet = parse_sheet(args.dilutions_sheet)

    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    skyline_csv = to_csv_if_needed(args.skyline, skyline_sheet, tmp_dir)
    dilutions_csv = to_csv_if_needed(args.dilutions, dilutions_sheet, tmp_dir)

    # Default output name
    out_path = args.out
    if out_path is None:
        base = Path(args.skyline).stem.replace(" ", "_")
        out_path = f"{base}__processed.csv"

    print(f"[run] Skyline:   {skyline_csv}")
    print(f"[run] Dilutions: {dilutions_csv}")
    print(f"[run] Output:    {out_path}")

    try:
        df = process_prm_data(skyline_csv, dilutions_csv, out_path)
    except Exception as e:
        sys.stderr.write("Pipeline failed.\n" + str(e) + "\n")
        return 1

    print(f"[done] Wrote {out_path} with {len(df)} rows.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
