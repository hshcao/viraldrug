#!/usr/bin/env bash

set -u
set -o pipefail

# ============================================================
# Run PLIP for top 100 CIF files for ONE strain at a time.
#
# Usage:
#   ./run_plip_one_strain_top100.sh congo
#   ./run_plip_one_strain_top100.sh portugal
#   ./run_plip_one_strain_top100.sh sudan
#   ./run_plip_one_strain_top100.sh ref
#
# Input folder structure:
#   congoAF3cif/OPG027cif_top200/001_EBC-03080_model.cif
#   congoAF3cif/OPG092cif_top200/001_EBC-157193_model.cif
#
# Output folder structure:
#   PLIP_top100_by_strain/congo/plip_congo_OPG027_R001_EBC-03080/
#   PLIP_top100_by_strain/congo/plip_congo_OPG092_R001_EBC-157193/
#
# Each PLIP output folder contains:
#   *_report.txt
#   *_report.xml
#   *.pse
#   *.png
# ----------------
# How to run it:
# screen -S plip_congo
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate plip
# cd /NFS/hcao/af3output/4strain8protCIFs
# ./run_plip_one_strain_top100.sh congo 2>&1 | tee run_plip_congo_top100_screen.log
# Check results:
# find PLIP_top100_by_strain/congo \
#   -mindepth 1 \
#   -maxdepth 1 \
#   -type d \
#   -name 'plip_congo_OPG*_R*_EBC-*' \
#   | wc -l
# find PLIP_top100_by_strain/congo -type f -name '*.pse' | wc -l
# find PLIP_top100_by_strain/congo -type f -name '*.png' | wc -l
# find PLIP_top100_by_strain/congo -type f -name '*_report.txt' | wc -l
# find PLIP_top100_by_strain/congo -type f -name '*_report.xml' | wc -l
# ============================================================

TOP_N=100

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 {congo|portugal|sudan|ref}"
    exit 1
fi

STRAIN="$1"

case "$STRAIN" in
    congo)
        STRAIN_CIF_DIR="congoAF3cif"
        ;;
    portugal)
        STRAIN_CIF_DIR="portugalAF3cif"
        ;;
    sudan)
        STRAIN_CIF_DIR="sudanAF3cif"
        ;;
    ref)
        STRAIN_CIF_DIR="refAF3cif"
        ;;
    *)
        echo "ERROR: unknown strain: $STRAIN"
        echo "Allowed values: congo, portugal, sudan, ref"
        exit 1
        ;;
esac

BASE_DIR="$(pwd)"
INPUT_DIR="${BASE_DIR}/${STRAIN_CIF_DIR}"

OUT_ROOT="${BASE_DIR}/PLIP_top100_by_strain"
OUT_BASE="${OUT_ROOT}/${STRAIN}"
PDB_DIR="${OUT_ROOT}/_converted_pdbs/${STRAIN}"
LOG_FILE="${OUT_BASE}/plip_${STRAIN}_top${TOP_N}_run.log"
STATUS_FILE="${OUT_BASE}/plip_${STRAIN}_top${TOP_N}_status.tsv"

# 1 = skip output folders that already have report.txt
# 0 = rerun everything
SKIP_DONE=1

mkdir -p "$OUT_BASE"
mkdir -p "$PDB_DIR"

echo -e "rank\tstrain\tprotein\tligand\tcif_file\tpdb_file\tplip_output\tstatus\tmessage" > "$STATUS_FILE"

echo "PLIP strain top-${TOP_N} run started: $(date)" | tee "$LOG_FILE"
echo "Strain: $STRAIN" | tee -a "$LOG_FILE"
echo "Input folder: $INPUT_DIR" | tee -a "$LOG_FILE"
echo "Output folder: $OUT_BASE" | tee -a "$LOG_FILE"
echo "PDB folder: $PDB_DIR" | tee -a "$LOG_FILE"
echo "Top N: $TOP_N" | tee -a "$LOG_FILE"
echo "Output naming: plip_${STRAIN}_OPGxxx_R001_EBC-xxxxx" | tee -a "$LOG_FILE"
echo | tee -a "$LOG_FILE"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input folder not found: $INPUT_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

for cmd in obabel plip pymol; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: command not found: $cmd" | tee -a "$LOG_FILE"
        exit 1
    fi
done

protein_count=0
job_count=0

while IFS= read -r -d '' protein_dir; do

    protein_folder="$(basename -- "$protein_dir")"
    protein="${protein_folder%%cif_top*}"

    protein_count=$((protein_count + 1))

    echo "============================================================" | tee -a "$LOG_FILE"
    echo "Protein: $protein" | tee -a "$LOG_FILE"
    echo "Protein CIF folder: $protein_dir" | tee -a "$LOG_FILE"
    echo "============================================================" | tee -a "$LOG_FILE"

    while IFS= read -r -d '' cif_file; do

        cif_base="$(basename -- "$cif_file")"

        rank="${cif_base%%_*}"
        ligand="$(echo "$cif_base" | sed -E 's/^[0-9]{3}_//; s/_model\.cif$//')"

        if [[ ! "$rank" =~ ^[0-9]{3}$ ]]; then
            echo "SKIP bad rank: $cif_file" | tee -a "$LOG_FILE"
            echo -e "NA\t${STRAIN}\t${protein}\tNA\t${cif_file}\tNA\tNA\tSKIP_BAD_RANK\tBad rank in filename" >> "$STATUS_FILE"
            continue
        fi

        rank_num=$((10#$rank))

        if (( rank_num > TOP_N )); then
            continue
        fi

        if [[ ! "$ligand" =~ ^EBC-[0-9]+$ ]]; then
            echo "SKIP bad ligand name: $cif_file" | tee -a "$LOG_FILE"
            echo -e "${rank}\t${STRAIN}\t${protein}\t${ligand}\t${cif_file}\tNA\tNA\tSKIP_BAD_LIGAND\tBad ligand name in filename" >> "$STATUS_FILE"
            continue
        fi

        job_count=$((job_count + 1))

        plip_name="plip_${STRAIN}_${protein}_R${rank}_${ligand}"
        pdb_file="${PDB_DIR}/${plip_name}.pdb"
        plip_out="${OUT_BASE}/${plip_name}"

        if [[ "$SKIP_DONE" -eq 1 ]] && [[ -d "$plip_out" ]] && find "$plip_out" -maxdepth 1 -type f -name '*_report.txt' | grep -q .; then
            echo "SKIP existing: $plip_out" | tee -a "$LOG_FILE"
            echo -e "${rank}\t${STRAIN}\t${protein}\t${ligand}\t${cif_file}\t${pdb_file}\t${plip_out}\tSKIP_DONE\tExisting report found" >> "$STATUS_FILE"
            continue
        fi

        echo "------------------------------------------------------------" | tee -a "$LOG_FILE"
        echo "Rank: R${rank}" | tee -a "$LOG_FILE"
        echo "Strain: $STRAIN" | tee -a "$LOG_FILE"
        echo "Protein: $protein" | tee -a "$LOG_FILE"
        echo "Ligand: $ligand" | tee -a "$LOG_FILE"
        echo "CIF: $cif_file" | tee -a "$LOG_FILE"
        echo "PDB: $pdb_file" | tee -a "$LOG_FILE"
        echo "PLIP output: $plip_out" | tee -a "$LOG_FILE"

        mkdir -p "$plip_out"

        if ! obabel "$cif_file" -O "$pdb_file" </dev/null >> "$LOG_FILE" 2>&1; then
            echo "FAILED obabel: $cif_file" | tee -a "$LOG_FILE"
            echo -e "${rank}\t${STRAIN}\t${protein}\t${ligand}\t${cif_file}\t${pdb_file}\t${plip_out}\tOBABEL_FAILED\tCIF to PDB conversion failed" >> "$STATUS_FILE"
            continue
        fi

        if [[ ! -s "$pdb_file" ]]; then
            echo "FAILED empty PDB: $pdb_file" | tee -a "$LOG_FILE"
            echo -e "${rank}\t${STRAIN}\t${protein}\t${ligand}\t${cif_file}\t${pdb_file}\t${plip_out}\tEMPTY_PDB\tPDB file missing or empty after conversion" >> "$STATUS_FILE"
            continue
        fi

        if ! plip \
            -f "$pdb_file" \
            -o "$plip_out" \
            -x \
            -t \
            -y \
            -p </dev/null >> "$LOG_FILE" 2>&1; then

            echo "FAILED PLIP: $pdb_file" | tee -a "$LOG_FILE"
            echo -e "${rank}\t${STRAIN}\t${protein}\t${ligand}\t${cif_file}\t${pdb_file}\t${plip_out}\tPLIP_FAILED\tPLIP analysis failed" >> "$STATUS_FILE"
            continue
        fi

        report_count="$(find "$plip_out" -maxdepth 1 -type f -name '*_report.txt' | wc -l)"
        xml_count="$(find "$plip_out" -maxdepth 1 -type f -name '*_report.xml' | wc -l)"
        pse_count="$(find "$plip_out" -maxdepth 1 -type f -name '*.pse' | wc -l)"
        png_count="$(find "$plip_out" -maxdepth 1 -type f -name '*.png' | wc -l)"

        echo "DONE: $plip_out" | tee -a "$LOG_FILE"
        echo "reports=${report_count}, xml=${xml_count}, pse=${pse_count}, png=${png_count}" | tee -a "$LOG_FILE"

        echo -e "${rank}\t${STRAIN}\t${protein}\t${ligand}\t${cif_file}\t${pdb_file}\t${plip_out}\tDONE\treport=${report_count};xml=${xml_count};pse=${pse_count};png=${png_count}" >> "$STATUS_FILE"

    done < <(
        find "$protein_dir" \
            -maxdepth 1 \
            -type f \
            -name '[0-9][0-9][0-9]_EBC-*_model.cif' \
            -print0 \
        | sort -z
    )

done < <(
    find "$INPUT_DIR" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name 'OPG*cif_top*' \
        -print0 \
    | sort -z
)

echo | tee -a "$LOG_FILE"
echo "Proteins processed: $protein_count" | tee -a "$LOG_FILE"
echo "Candidate top-${TOP_N} PLIP jobs: $job_count" | tee -a "$LOG_FILE"
echo "PLIP strain top-${TOP_N} run finished: $(date)" | tee -a "$LOG_FILE"
echo "Status file: $STATUS_FILE" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"

