#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 University of Utah

import csv
import os
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

GENE_COL_DEFAULTS = ["Element symbol", "Gene symbol"]
SUBTYPE_COL_DEFAULTS = ["Subtype", "Element subtype"]
CLASS_COL_DEFAULTS = ["Class"]
SUBCLASS_COL_DEFAULTS = ["Subclass"]

#--------
# Helpers
#--------



def simplify_gene_names(gene_name):
    """
    Remove/replace special characters for Bayesian network structure learning
    in ReGAIN Module 2.
    """
    if gene_name is None:
        return None
    gene_name = str(gene_name)
    gene_name = gene_name.replace("'", "p")
    gene_name = gene_name.replace('"', "pp")
    gene_name = gene_name.replace(".", "_")
    gene_name = gene_name.replace("(", "")
    gene_name = gene_name.replace(")", "")
    gene_name = gene_name.replace("-", "_")
    gene_name = gene_name.replace("/", "_")
    return gene_name

def detect_delimiter(file_path: Path) -> str:
    """
    Infer whether a file is comma- or tab-delimited.
    Defaults to comma if detection is ambiguous.
    """
    with file_path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)

    if not sample.strip():
        return ","

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        comma_count = first_line.count(",")
        tab_count = first_line.count("\t")
        return "\t" if tab_count > comma_count else ","
    
def find_amrfinder_files(directory: Path, expected_extension_string: str) -> List[Path]:
    """
    Finder AMRFinder output files in a directory using a required filename ending

    Default expected ending should remain '.amrfinder.csv' to preserve the regain AMR -> regain matrix workflow
    """
    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")

    if not directory.is_dir():
        raise NotADirectoryError(f"Not a directory: {directory}")

    files = [
        p for p in sorted(directory.iterdir())
        if p.is_file() and p.name.endswith(expected_extension_string)
    ]

    if not files:
        raise FileNotFoundError(
            f"No files ending with '{expected_extension_string}' were found in {directory}\n"
        )

    return files

def read_amrfinder_table(file_path: Path) -> pd.DataFrame:
    """
    Read one AMRFinder output table using detected delimiter
    """
    delimiter = detect_delimiter(file_path)
    df = pd.read_csv(file_path, sep=delimiter, na_filter=False, dtype=str)
    df.columns = [str(c).strip() for c in df.columns]
    return df

def resolve_column_name(
    df: pd.DataFrame,
    user_value: str,
    fallback_options: List[str],
    arg_name: str,
    file_path: Path,
) -> str:
    """
    Resolve a column name from either:
        1. an explicit user-provided value or
        2. the first matching fallback column in fallback_options

    If user_value is provided and exists, use it exactly
    If user_value is provided and does not exist, raise an error
    If not fallback matches, raise an error
    """
    columns = [str(c).strip() for c in df.columns]

    # If a user provides the values
    if user_value is not None:
        if user_value in columns:
            return user_value
        raise ValueError(
            f"Requested {arg_name} column '{user_value}' was not found in {file_path.name}\n"
            f"Available columns include: {', '.join(columns)}\n"
        )
    
    # No user override; try fallback options in order
    for candidate in fallback_options:
        if candidate in columns:
            return candidate
    
    raise ValueError(
        f"Count not resolve {arg_name} column in {file_path.name}\n"
        f"Tried: {', '.join(fallback_options)}\n"
        f"Available columns include: {', '.join(columns)}\n"

    )

def validate_required_columns(
    df: pd.DataFrame,
    file_path: Path,
    subtype_col: str,
    gene_col: str,
    class_col: str,
    subclass_col: str,
) -> None:
    """
    Confirm required columns exist in one AMRFinder table
    """
    required = [subtype_col, gene_col, class_col, subclass_col]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValueError(
            f"Required column(s) missing in {file_path.name}: {', '.join(missing)}"
        )
    
def get_sample_id(file_path: Path, expected_extension_string: str) -> str:
    """
    Generate a sample ID from filename.

    Default logic preserves the regain AMR convention:
        Sample1.fna.amrfinder.csv -> Sample1
        Sample1.fa.amrfinder.csv  -> Sample1
        Sample1.fasta.amrfinder.csv -> Sample1
    """
    filename = file_path.name

    if filename.endswith(expected_extension_string):
        sample_id = filename[: -len(expected_extension_string)]
    else:
        sample_id = file_path.stem

    for fasta_ext in [".fna", ".fa", ".fasta", ".fas"]:
        if sample_id.endswith(fasta_ext):
            sample_id = sample_id[: -len(fasta_ext)]
            break

    return sample_id.rstrip(".")

def normalize_string(value) -> str:
    """
    Normalize cell values to stripped strings
    """
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() == "nan" else value

def filter_by_gene_type(
    df: pd.DataFrame,
    gene_type: str,
    subtype_col: str,
) -> pd.DataFrame:
    """
    Filter AMRFinder rows by regain matrix gene-type behavior
    """
    subtype_series = df[subtype_col].astype(str).str.strip()

    subtype_map = {
        "resistance": ["AMR", "METAL", "BIOCIDE", "POINT"],
        "virulence": ["VIRULENCE", "HEAT", "ACID"],
        "all": ["AMR", "METAL", "BIOCIDE", "POINT", "VIRULENCE", "HEAT", "ACID"],
    }

    allowed = subtype_map[gene_type]
    return df.loc[subtype_series.isin(allowed)].copy()

def combine_amrfinder_output_files(
    dataframes_by_file: Dict[Path, pd.DataFrame],
    output_file_path: Path,
    gene_col: str,
    expected_extension_string: str,
    verbose_gene_report: bool = False,
    source_col_name: str = "SourceFile",
) -> pd.DataFrame:
    """
    Combine parsed AMRFinder tables and write an unfiltered combined file

    Default behavior:
        Drop duplicate rows based on gene_col to create a gene dictionary

    Verbose behavior:
        Retain duplicate gene rows and append a source-file column    
    """
    if not dataframes_by_file:
        raise ValueError("No valid AMRFinder tables available to combine\n")

    dfs = []

    for file_path, df in dataframes_by_file.items():
        working = df.copy()

        if verbose_gene_report:
            working[source_col_name] = get_sample_id(file_path, expected_extension_string)
        dfs.append(working)

    combined_df = pd.concat(dfs, ignore_index=True)

    if not verbose_gene_report:
        combined_df = combined_df.drop_duplicates(subset=[gene_col]).copy()

    combined_df.to_csv(output_file_path, index=False, sep="\t")
    print(f"Combined results saved to {output_file_path}")

    return combined_df

def build_metadata(
    combined_df: pd.DataFrame,
    gene_col: str,
    class_col: str,
    subclass_col: str,
) -> pd.DataFrame:
    """
    Build metadata for unique genes

    Uses the first non-empty class/subclass encountered per gene
    """
    working = combined_df[[gene_col, class_col, subclass_col]].copy()
    working[gene_col] = working[gene_col].map(normalize_string)
    working[class_col] = working[class_col].map(normalize_string)
    working[subclass_col] = working[subclass_col].map(normalize_string)

    working = working.loc[working[gene_col] != ""].copy()

    metadata_records = []

    for gene, group in working.groupby(gene_col, sort=True):
        class_values = [v for v in group[class_col].tolist() if v != ""]
        subclass_values = [v for v in group[subclass_col].tolist() if v != ""]

        gene_class = class_values[0] if class_values else "virulence"
        gene_subclass = subclass_values[0] if subclass_values else "virulence"

        metadata_records.append({
            "Gene": gene,
            "GeneClass": str(gene_class).title(),
            "GeneSubClass": str(gene_subclass).title(),
        })

    metadata_df = pd.DataFrame(metadata_records)

    if metadata_df.empty:
        raise ValueError(
            "No genes remained after filtering. Check --gene-type and input AMRFinder tables\n"
        )

    return metadata_df

def build_presence_absence_matrix(
   file_tables: Dict[Path, pd.DataFrame],
   metadata_df: pd.DataFrame,
   gene_col: str,
   expected_extension_string: str, 
) -> pd.DataFrame:
    """
    Build a presence/absence matrix using exact gene-column membership
    """
    genes = metadata_df["Gene"].tolist()
    rows = []

    for file_path, df in file_tables.items():
        gene_values = (
            df[gene_col]
            .astype(str)
            .map(normalize_string)
        )
        gene_set = set(v for v in gene_values.tolist() if v != "")

        row = {"file": get_sample_id(file_path, expected_extension_string)}
        for gene in genes:
            row[gene] = 1 if gene in gene_set else 0

        rows.append(row)

    matrix_df = pd.DataFrame(rows)
    matrix_df = matrix_df[["file"] + genes]
    return matrix_df

def build_count_matrix(
    file_tables: Dict[Path, pd.DataFrame],
    metadata_df: pd.DataFrame,
    gene_col: str,
    expected_extension_string: str,
) -> pd.DataFrame:
    """
    Build and exact-count matrix using gene-column value counts
    """
    genes = metadata_df["Gene"].tolist()
    rows = []

    for file_path, df in file_tables.items():
        gene_values = (
            df[gene_col]
            .astype(str)
            .map(normalize_string)
        )
        gene_counts = gene_values[gene_values != ""].value_counts().to_dict()

        row = {"file": get_sample_id(file_path, expected_extension_string)}
        for gene in genes:
            row[gene] = int(gene_counts.get(gene, 0))

        rows.append(row)

    count_df = pd.DataFrame(rows)
    count_df = count_df[["file"] + genes]
    return count_df

def apply_min_max_filter(
    matrix_df: pd.DataFrame,
    required_min: int,
    required_max: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter matrix columns by min/max occurrence thresholds
    """
    feature_cols = [c for c in matrix_df.columns if c != "file"]

    sums_df = pd.DataFrame({
        "variable": feature_cols,
        "sum": [int(matrix_df[col].sum()) for col in feature_cols],
    })

    if required_max is None:
        cols_to_keep = sums_df.loc[sums_df["sum"] >= required_min, "variable"].tolist()
    else:
        cols_to_keep = sums_df.loc[
            (sums_df["sum"] >= required_min) & (sums_df["sum"] <= required_max),
            "variable"
        ].tolist()

    filtered_df = matrix_df[["file"] + cols_to_keep].copy()

    return filtered_df, sums_df

def synchronize_metadata_to_matrix(
    metadata_df: pd.DataFrame,
    matrix_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Keep only metadata rows for genes retained in the matrix
    """
    kept_genes = set(matrix_df.columns[1:])
    synced = metadata_df.loc[metadata_df["Gene"].isin(kept_genes)].copy()

    # Preserve matrix column order.
    order_map = {gene: i for i, gene in enumerate(matrix_df.columns[1:])}
    synced["__order"] = synced["Gene"].map(order_map)
    synced = synced.sort_values("__order").drop(columns="__order").reset_index(drop=True)

    return synced

#----------------------------
# Main entry point for cli.py
#----------------------------


def run(args):
    path = Path(args.directory)
    output_directory = Path(os.getcwd()) / args.output_dir
    
    if output_directory.exists():
        if args.force_overwrite:
            shutil.rmtree(output_directory)
        else:
            raise ValueError(
                f"\n{output_directory.name} already exists in the current working directory\n"
                "Rename or move the existing output directory to avoid overwriting, \n"
                "or re-run with --force-overwrite if you are sure you want to replace it\n"
            )
        
    output_directory.mkdir(parents=True, exist_ok=False)

    filtered_matrix_path = output_directory / "filtered_matrix.csv"
    metadata_path = output_directory / "metadata.csv"
    unfiltered_combined_path = output_directory / "combined_AMR_results_unfiltered.csv"

    # Optional new arguments with backwards-compatible defaults.
    expected_extension_string = getattr(args, "expected_extension_string", ".amrfinder.csv")

    user_subtype_col = getattr(args, "subtype_col", None)
    user_gene_col = getattr(args, "gene_col", None)
    user_class_col = getattr(args, "class_col", "Class")
    user_subclass_col = getattr(args, "subclass_col", "Subclass")

    gene_type = args.gene_type.lower()

    if gene_type not in {"resistance", "virulence", "all"}:
        raise ValueError(f"Unsupported gene_type: {args.gene_type}\n")

    files = find_amrfinder_files(path, expected_extension_string)

    raw_tables_by_file: Dict[Path, pd.DataFrame] = {}
    filtered_tables_by_file: Dict[Path, pd.DataFrame] = {}

    for file_path in files:
        try:
            df = read_amrfinder_table(file_path)

            resolved_gene_col = resolve_column_name(
                df=df,
                user_value=user_gene_col,
                fallback_options=GENE_COL_DEFAULTS,
                arg_name="--gene-col",
                file_path=file_path,
            )

            resolved_subtype_col = resolve_column_name(
                df=df,
                user_value=user_subtype_col,
                fallback_options=SUBTYPE_COL_DEFAULTS,
                arg_name="--subtype-col",
                file_path=file_path,
            )

            resolved_class_col = resolve_column_name(
                df=df,
                user_value=user_class_col,
                fallback_options=CLASS_COL_DEFAULTS,
                arg_name="--class-col",
                file_path=file_path,
            )

            resolved_subclass_col = resolve_column_name(
                df=df,
                user_value=user_subclass_col,
                fallback_options=SUBCLASS_COL_DEFAULTS,
                arg_name="--subclass-col",
                file_path=file_path,
            )

            # Rename per-file columns to stable internal names so the rest of the code
            # does not care which AMRFinder header variant was present
            df = df.rename(columns={
                resolved_gene_col: "GeneColumn",
                resolved_subtype_col: "SubtypeColumn",
                resolved_class_col: "ClassColumn",
                resolved_subclass_col: "SubclassColumn",
            })

            validate_required_columns(
                df=df,
                file_path=file_path,
                subtype_col="SubtypeColumn",
                gene_col="GeneColumn",
                class_col="ClassColumn",
                subclass_col="SubclassColumn",
            )

            raw_tables_by_file[file_path] = df
            filtered_tables_by_file[file_path] = filter_by_gene_type(
                df=df,
                gene_type=gene_type,
                subtype_col="SubtypeColumn",
            )

        except Exception as e:
            print(f"Skipping {file_path.name} due to an error: {e}\n")

    if not raw_tables_by_file:
        raise ValueError("No valid AMRFinder tables were available after parsing\n")

    # Combine filtered tables so downstream metadata and matrices reflect gene_type.
    combined_df = combine_amrfinder_output_files(
        dataframes_by_file=filtered_tables_by_file,
        output_file_path=unfiltered_combined_path,
        gene_col="GeneColumn",
        expected_extension_string=expected_extension_string,
        verbose_gene_report=args.verbose_gene_report,
    )

    metadata_df = build_metadata(
        combined_df=combined_df,
        gene_col="GeneColumn",
        class_col="ClassColumn",
        subclass_col="SubclassColumn",
    )

    matrix_df = build_presence_absence_matrix(
        file_tables=filtered_tables_by_file,
        metadata_df=metadata_df,
        gene_col="GeneColumn",
        expected_extension_string=expected_extension_string,
    )

    filtered_matrix_df, sums_df = apply_min_max_filter(
        matrix_df=matrix_df,
        required_min=args.min,
        required_max=args.max,
    )

    metadata_df = synchronize_metadata_to_matrix(
        metadata_df=metadata_df,
        matrix_df=filtered_matrix_df,
    )

    if not args.keep_gene_names:
        metadata_df["Gene"] = metadata_df["Gene"].apply(simplify_gene_names)
        filtered_matrix_df.columns = [
            simplify_gene_names(col) if col != "file" else col
            for col in filtered_matrix_df.columns
        ]

    metadata_df.to_csv(metadata_path, index=False)
    filtered_matrix_df.to_csv(filtered_matrix_path, index=False)

    if args.report_all:
        unfiltered_matrix_path = output_directory / "unfiltered_matrix.csv"
        count_matrix_df = build_count_matrix(
        file_tables=filtered_tables_by_file,
        metadata_df=build_metadata(
            combined_df=combined_df,
            gene_col="GeneColumn",
            class_col="ClassColumn",
            subclass_col="SubclassColumn",
        ),
        gene_col="GeneColumn",
        expected_extension_string=expected_extension_string,
    )

        if not args.keep_gene_names:
            count_matrix_df.columns = [
                simplify_gene_names(col) if col != "file" else col
                for col in count_matrix_df.columns
            ]

        count_matrix_df.to_csv(unfiltered_matrix_path, index=False)
        print(f"Unfiltered matrix saved to {unfiltered_matrix_path}\n")

    print(f"\nFiltered matrix saved to {filtered_matrix_path}\n")
    print(f"Metadata file saved to {metadata_path}\n")
