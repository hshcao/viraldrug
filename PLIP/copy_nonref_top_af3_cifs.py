#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import pandas as pd


STRAIN_CONFIG = {
    "congo": {
        "json_dir": "jsonCongo_output",
        "excel": "congo_af3drugs_all10prot_ranked.xlsx",
        "output_dir": "congoAF3cif",
    },
    "portugal": {
        "json_dir": "jsonPortugal_output",
        "excel": "portugal_af3drugs_all10prot_ranked.xlsx",
        "output_dir": "portugalAF3cif",
    },
    "sudan": {
        "json_dir": "jsonSudan_output",
        "excel": "sudan_af3drugs_all10prot_ranked.xlsx",
        "output_dir": "sudanAF3cif",
    },
}


PROTEINS = [
    "OPG027",
    "OPG112",
    "OPG174",
    "OPG198",
]


def ligand_to_ebc(ligand_id: str) -> str:
    ligand_id = str(ligand_id).strip()

    if ligand_id.startswith("EBC-"):
        return ligand_id.replace("EBC-", "ebc-", 1)

    if ligand_id.startswith("ebc-"):
        return ligand_id

    if ligand_id.startswith("jebc-"):
        return ligand_id.replace("jebc-", "ebc-", 1)

    raise ValueError(f"Unexpected ligandID format: {ligand_id}")


def normalize_output_ligand_id(ligand_id: str) -> str:
    ligand_id = str(ligand_id).strip()

    if ligand_id.startswith("EBC-"):
        return ligand_id

    if ligand_id.startswith("ebc-"):
        return ligand_id.replace("ebc-", "EBC-", 1)

    if ligand_id.startswith("jebc-"):
        return ligand_id.replace("jebc-", "EBC-", 1)

    return ligand_id


def read_protein_sheet(excel_file: Path, protein: str) -> pd.DataFrame:
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


def find_cif_file(protein_dir: Path, ebc_id: str) -> Path | None:
    candidate_dirs = []

    exact_dir = protein_dir / ebc_id
    if exact_dir.is_dir():
        candidate_dirs.append(exact_dir)

    candidate_dirs.extend(sorted(protein_dir.glob(f"{ebc_id}*")))

    unique_dirs = []
    seen = set()

    for d in candidate_dirs:
        if d.is_dir() and d not in seen:
            unique_dirs.append(d)
            seen.add(d)

    for d in unique_dirs:
        exact_cif = d / f"{ebc_id}_model.cif"

        if exact_cif.is_file():
            return exact_cif

        matches = sorted(d.glob("*_model.cif"))
        if matches:
            return matches[0]

    return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Copy top AF3 model CIF files for Congo, Portugal, or Sudan "
            "using iptmRank then ptmRank, for OPG027, OPG112, OPG174, and OPG198."
        )
    )

    parser.add_argument(
        "--strain",
        required=True,
        choices=["congo", "portugal", "sudan"],
        help="Strain to process: congo, portugal, or sudan.",
    )

    parser.add_argument(
        "--base",
        default=".",
        help="Base folder containing json strain output folders, Excel files, and AF3cif output folders.",
    )

    parser.add_argument(
        "--excel",
        default=None,
        help="Optional Excel file. If omitted, uses the default Excel file for the strain.",
    )

    parser.add_argument(
        "--json-dir",
        default=None,
        help="Optional AF3 output folder. If omitted, uses the default output folder for the strain.",
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
    config = STRAIN_CONFIG[args.strain]

    excel_file = Path(args.excel) if args.excel else base_dir / config["excel"]
    json_dir = Path(args.json_dir) if args.json_dir else base_dir / config["json_dir"]
    strain_output_dir = base_dir / config["output_dir"]

    excel_file = excel_file.resolve()
    json_dir = json_dir.resolve()
    strain_output_dir = strain_output_dir.resolve()

    if not excel_file.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_file}")

    if not json_dir.is_dir():
        raise NotADirectoryError(f"AF3 output folder not found: {json_dir}")

    print("========================================")
    print(f"Strain:            {args.strain}")
    print(f"Excel file:        {excel_file}")
    print(f"AF3 output folder: {json_dir}")
    print(f"Main output folder:{strain_output_dir}")
    print(f"Proteins:          {', '.join(PROTEINS)}")
    print(f"Top N:             {args.top}")
    print(f"Dry run:           {args.dry_run}")
    print("Ranking rule:      iptmRank ascending, then ptmRank ascending")
    print("========================================")

    report_rows = []

    if not args.dry_run:
        strain_output_dir.mkdir(exist_ok=True)

    for protein in PROTEINS:
        print(f"\n=== Processing {protein} ===")

        protein_dir = json_dir / protein

        if not protein_dir.is_dir():
            print(f"WARNING: Missing protein folder: {protein_dir}")
            report_rows.append(
                [args.strain, protein, "", "", "", "", "MISSING_PROTEIN_FOLDER"]
            )
            continue

        try:
            df = read_protein_sheet(excel_file, protein)
        except ValueError as e:
            print(f"WARNING: {e}")
            report_rows.append(
                [args.strain, protein, "", "", "", "", "MISSING_EXCEL_SHEET"]
            )
            continue

        df_top = (
            df.dropna(subset=["ligandID", "iptmRank", "ptmRank"])
            .sort_values(["iptmRank", "ptmRank"], ascending=[True, True])
            .head(args.top)
        )

        output_dir = strain_output_dir / f"{protein}cif_top{args.top}"

        if not args.dry_run:
            output_dir.mkdir(exist_ok=True)

        for rank_number, (_, row) in enumerate(df_top.iterrows(), start=1):
            ligand_id = str(row["ligandID"]).strip()
            ebc_id = ligand_to_ebc(ligand_id)

            cif_file = find_cif_file(protein_dir, ebc_id)

            output_ligand_id = normalize_output_ligand_id(ligand_id)
            dest_name = f"{rank_number:03d}_{output_ligand_id}_model.cif"
            dest_file = output_dir / dest_name

            if cif_file is None:
                print(f"MISSING: {protein} {ligand_id} expected {ebc_id}_model.cif")
                report_rows.append(
                    [
                        args.strain,
                        protein,
                        rank_number,
                        ligand_id,
                        ebc_id,
                        "",
                        "MISSING",
                    ]
                )
                continue

            if args.dry_run:
                print(f"DRY-RUN: {cif_file} -> {dest_file}")
                report_rows.append(
                    [
                        args.strain,
                        protein,
                        rank_number,
                        ligand_id,
                        ebc_id,
                        str(cif_file),
                        "DRY_RUN",
                    ]
                )
            else:
                shutil.copy2(cif_file, dest_file)
                print(f"COPIED: {cif_file} -> {dest_file}")
                report_rows.append(
                    [
                        args.strain,
                        protein,
                        rank_number,
                        ligand_id,
                        ebc_id,
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
            "ebcID",
            "source_cif",
            "status",
        ],
    )

    report_file = strain_output_dir / f"{args.strain}_copy_top{args.top}_af3_cif_report.tsv"

    if not args.dry_run:
        strain_output_dir.mkdir(exist_ok=True)

    report.to_csv(report_file, sep="\t", index=False)

    print("\n========================================")
    print(f"Report written to: {report_file}")
    print("Done.")
    print("========================================")


if __name__ == "__main__":
    main()
