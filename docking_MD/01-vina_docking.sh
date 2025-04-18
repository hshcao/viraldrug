#!/bin/bash

## Script to perform molecular docking using AutoDock Vina v1.2.7
## Website to download Vina: https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7
## After downloading, rename the binary to `vina` and make it executable:
## chmod +x vina

# --------------------------------------
# Docking for Bioactive Ligands
# --------------------------------------

for file in ../pdbqtLigands/*.pdbqt; do
  baseName="${file##*/}"                   # Extract filename (e.g., ligand1.pdbqt)
  prefix="${baseName%.pdbqt}"              # Remove extension (e.g., ligand1)
  logFile="${prefix}.log"                  # Define the log filename

  # Run docking only if the log file does not already exist
  if [[ ! -f "$logFile" ]]; then
    ./vina --config config.txt \           # Configuration for docking (receptor, grid box, etc.)
           --ligand "$file" \              # Input ligand file
           --out ${prefix}.pdbqt \       # Output docked pose
           > "$logFile"                    # Save docking log (binding affinity, etc.)
  fi
done

# Create a directory to store bioactive docking results
mkdir -p bioactiveResults                 # -p prevents error if directory already exists
mv EBC*.log bioactiveResults/            # Move all log files starting with "EBC" to results folder

# --------------------------------------
# Docking for Control Ligands
# --------------------------------------

for file in ../controls/*.pdbqt; do
  baseName="${file##*/}"                   # Extract filename
  prefix="${baseName%.pdbqt}"              # Remove extension
  logFile="${prefix}.log"                  # Define the log filename

  # Run docking only if the log file does not already exist
  if [[ ! -f "$logFile" ]]; then
    ./vina --config config.txt \           # Configuration remains the same
           --ligand "$file" \              # Input control ligand
           --out ${prefix}.pdbqt \            # Output is overwritten for each control (consider changing if needed)
           > "$logFile"                    # Save docking log
  fi
done

# --------------------------------------
# Extract Energy Scores
# --------------------------------------

./extract_energy.sh > ranked_energy.log   # Run custom script to extract and rank binding affinities