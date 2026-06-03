#!/bin/bash

echo "Downloading PDB..."
python scripts/download_pdb.py

echo "Cleaning structure..."
python scripts/clean_structure.py

echo "Extracting sequence..."
python scripts/extract_sequence.py

echo "Running ESM mutation generation..."
python esm/generate_mutations.py

echo "Running Rosetta ddG scan..."
python rosetta/ddg_scan.py

echo "Selecting top 5..."
python analysis/select_top5.py

echo "DONE."