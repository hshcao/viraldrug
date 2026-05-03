#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import pandas as pd


# Ref strain target proteins
PROTEINS = [
    "OPG198",
    "OPG027",
    "OPG199",
    "OPG120",
    "OPG154",
    "OPG174",
    "OPG092",
    "OPG112",
]


def ligand_to_jebc(ligand_id: str) -> str:
    """
    Convert Excel ligandID to Ref AF3 folder/file format.

    Examples:
        EBC-13852 -> jebc-13852
        EBC-00009 -> jebc-00009
        ebc-00009 -> jebc-00009
    """
    ligand_id = str(ligand_id).strip()

    if ligand_id.startswith("EBC-"):
        return ligand_id.replace("EBC-", "jebc-", 1)

    if ligand_id.startswith("ebc-"):
        return ligand_id.replace("ebc-", "jebc-", 1)

    if ligand_id.startswith("jebc-"):
        return ligand_id

    raise ValueError(f"Unexpected ligandID format: {ligand_id}")


def normalize_output_ligand_id(ligand_id: str) -> str:
    """
    Standardize output file names to EBC-*.

    Examples:
        jebc-13852 -> EBC-13852
        ebc-13852  -> EBC-13852
        EBC-13852  -> EBC-13852
    """
    ligand_id = str(ligand_id).strip()

    if ligand_id.startswith("EBC-"):
        return ligand_id

    if ligand_id.startswith("ebc-"):
        return ligand_id.replace("ebc-", "EBC-", 1)

    if ligand_id.startswith("jebc-"):
        return ligand_id.replace("jebc-", "EBC-", 1)

    return ligand_id


def read_protein_sheet(excel_file: Path, protein: str) -> pd.DataFrame:
    """
    Read one protein sheet from the Excel file.
    Required columns: ligandID, iptmRank, ptmRank.
    """
    try:
        df = pd.read_excel(excel_file, sheet_name=protein)
    except ValueError:
        raise ValueError(f"Could not find sheet named {protein} in {excel_file}")

    required_columns = {"ligandID", "iptmRank", "ptmRank"}
    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(
            f"Sheet {protein} is missing columns: {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    return df


def find_cif_file(protein_af3_dir: Path, jebc_id: str) -> Path | None:
    """
    Find the model CIF file for one ligand.

    Expected Ref structure:
        OPG198_af3output/jebc-13852/jebc-13852_model.cif

    Also handles possible suffixed folders:
        OPG112_af3output/jebc-00002_20250416_105957/jebc-00002_model.cif
    """
    candidate_dirs = []

    # Exact folder match
    exact_dir = protein_af3_dir / jebc_id
    if exact_dir.is_dir():
        candidate_dirs.append(exact_dir)

    # Folders with possible timestamp suffixes
    candidate_dirs.extend(sorted(protein_af3_dir.glob(f"{jebc_id}*")))

    # Remove duplicate directories while preserving order
    unique_dirs = []
    seen = set()

    for d in candidate_dirs:
        if d.is_dir() and d not in seen:
            unique_dirs.append(d)
            seen.add(d)

    for d in unique_dirs:
        exact_cif = d / f"{jebc_id}_model.cif"

        if exact_cif.is_file():
            return exact_cif

        # Fallback: any *_model.cif in that ligand folder
        matches = sorted(d.glob("*_model.cif"))
        if matches:
            return matches[0]

    return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Copy top Ref AF3 model CIF files using iptmRank then ptmRank. "
            "Outputs are organized into refAF3cif/ with rank-numbered CIF names."
        )
    )

    parser.add_argument(
        "--excel",
        default="ref_af3drugs_all10prot_ranked.xlsx",
        help="Excel file containing sheets named OPG198, OPG027, etc.",
    )

    parser.add_argument(
        "--base",
        default=".",
        help="Base folder containing OPG198_af3output, OPG027_af3output, etc.",
    )

    parser.add_argument(
        "--top",
        type=int,
        default=100,
        help="Number of top ligands to copy per protein.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview only; do not copy files.",
    )

    args = parser.parse_args()

    base_dir = Path(args.base).resolve()
    excel_file = Path(args.excel).resolve()
    ref_output_dir = base_dir / "refAF3cif"

    if not excel_file.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_file}")

    print("========================================")
    print("Strain:             Ref")
    print(f"Excel file:         {excel_file}")
    print(f"Base folder:        {base_dir}")
    print(f"Main output folder: {ref_output_dir}")
    print(f"Proteins:           {', '.join(PROTEINS)}")
    print(f"Top N:              {args.top}")
    print(f"Dry run:            {args.dry_run}")
    print("Ranking rule:       iptmRank ascending, then ptmRank ascending")
    print("Output naming:      001_EBC-xxxxx_model.cif")
    print("========================================")

    report_rows = []

    if not args.dry_run:
        ref_output_dir.mkdir(exist_ok=True)

    for protein in PROTEINS:
        print(f"\n=== Processing {protein} ===")

        protein_af3_dir = base_dir / f"{protein}_af3output"

        if not protein_af3_dir.is_dir():
            print(f"WARNING: Missing AF3 folder: {protein_af3_dir}")
            report_rows.append(
                ["ref", protein, "", "", "", "", "MISSING_AF3_FOLDER"]
            )
            continue

        try:
            df = read_protein_sheet(excel_file, protein)
        except ValueError as e:
            print(f"WARNING: {e}")
            report_rows.append(
                ["ref", protein, "", "", "", "", "MISSING_EXCEL_SHEET"]
            )
            continue

        # Ranking logic:
        # lowest iptmRank is best; lowest ptmRank breaks ties
        df_top = (
            df.dropna(subset=["ligandID", "iptmRank", "ptmRank"])
            .sort_values(["iptmRank", "ptmRank"], ascending=[True, True])
            .head(args.top)
        )

        output_dir = ref_output_dir / f"{protein}cif_top{args.top}"

        if not args.dry_run:
            output_dir.mkdir(exist_ok=True)

        for rank_number, (_, row) in enumerate(df_top.iterrows(), start=1):
            ligand_id = str(row["ligandID"]).strip()
            jebc_id = ligand_to_jebc(ligand_id)

            cif_file = find_cif_file(protein_af3_dir, jebc_id)

            output_ligand_id = normalize_output_ligand_id(ligand_id)
            dest_name = f"{rank_number:03d}_{output_ligand_id}_model.cif"
            dest_file = output_dir / dest_name

            if cif_file is None:
                print(f"MISSING: {protein} {ligand_id} expected {jebc_id}_model.cif")
                report_rows.append(
                    [
                        "ref",
                        protein,
                        rank_number,
                        ligand_id,
                        jebc_id,
                        "",
                        "MISSING",
                    ]
                )
                continue

            if args.dry_run:
                print(f"DRY-RUN: {cif_file} -> {dest_file}")
                report_rows.append(
                    [
                        "ref",
                        protein,
                        rank_number,
                        ligand_id,
                        jebc_id,
                        str(cif_file),
                        "DRY_RUN",
                    ]
                )
            else:
                shutil.copy2(cif_file, dest_file)
                print(f"COPIED: {cif_file} -> {dest_file}")
                report_rows.append(
                    [
                        "ref",
                        protein,
                        rank_number,
                        ligand_id,
                        jebc_id,
                        str(cif_file),
                        "COPIED",
                    ]
                )

    report = pd.DataFrame(
        report_rows,
        columns=[
            "strain",
            "protein",
            "rank_number",
            "ligandID",
            "jebcID",
            "source_cif",
            "status",
        ],
    )

    if not args.dry_run:
        ref_output_dir.mkdir(exist_ok=True)

    report_file = ref_output_dir / f"ref_copy_top{args.top}_af3_cif_report.tsv"
    report.to_csv(report_file, sep="\t", index=False)

    print("\n========================================")
    print(f"Report written to: {report_file}")
    print("Done.")
    print("========================================")


if __name__ == "__main__":
    main()
