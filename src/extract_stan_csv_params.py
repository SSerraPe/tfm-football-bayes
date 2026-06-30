#!/usr/bin/env python3
"""
Column-pass extractor for CmdStan CSV output files.

Reads a single chain CSV file and extracts specific parameter columns
without loading the full file into memory. CmdStan CSVs have comment lines
starting with '#', then a header row, then data rows.

Usage:
  python extract_stan_csv_params.py <csv_path> <chain_id> <out_csv> \
    param1:col_idx1 param2:col_idx2 ...

Arguments:
  csv_path   Path to the CmdStan chain CSV file
  chain_id   Integer chain identifier (appended as .chain column)
  out_csv    Output CSV path
  param:idx  One or more name:col_idx pairs (col_idx is 0-based from the
             CmdStan CSV header, NOT the Python column index — they are the
             same: 0 = first column = lp__)

Example:
  python extract_stan_csv_params.py chain1.csv 1 out.csv nu:7 sigma_e.1:91605
"""

import sys
import csv


def extract(csv_path, chain_id, out_csv, param_col_map):
    """
    param_col_map: list of (name, 0-based-col-index) tuples
    """
    col_indices = [idx for _, idx in param_col_map]
    col_names = [name for name, _ in param_col_map]

    with open(csv_path, "r", newline="") as f_in, \
         open(out_csv, "w", newline="") as f_out:

        writer = csv.writer(f_out)
        writer.writerow([".chain", ".iteration"] + col_names)

        reader = csv.reader(f_in)
        header_seen = False
        iteration = 0

        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if not header_seen:
                # First non-comment row is the header — skip it
                header_seen = True
                continue
            # Data row
            iteration += 1
            try:
                values = [row[i] for i in col_indices]
            except IndexError:
                # Trailing partial rows (rare); skip
                continue
            writer.writerow([chain_id, iteration] + values)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    csv_path = sys.argv[1]
    chain_id = int(sys.argv[2])
    out_csv = sys.argv[3]
    raw_pairs = sys.argv[4:]

    param_col_map = []
    for pair in raw_pairs:
        name, idx_str = pair.rsplit(":", 1)
        param_col_map.append((name, int(idx_str)))

    extract(csv_path, chain_id, out_csv, param_col_map)
