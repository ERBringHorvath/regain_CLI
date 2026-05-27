#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 University of Utah

"""
collapse_features.py

ReGAIN utility for collapsing presence/absence matrix features using a
user-provided Gene-to-bin mapping.

Expected CLI usage through cli.py:

    regain collapse-features \
        -i presence_absence.csv \
        -M metadata_with_bins.csv \
        -o collapsed_presence_absence.csv

Required metadata columns by default:
    Gene,bin

Example metadata:
    Gene,bin
    aac6p_Ib_cr5,aac6p_Ib_cr_group
    aac6p_Ib_cr6,aac6p_Ib_cr_group
    tetA,tet_group
    tetB,tet_group

Default behavior:
    - Features listed in metadata are collapsed according to the bin column.
    - Features not listed in metadata are kept unchanged unless --drop-unmapped is passed
    - Collapsed feature values are forced back to binary 0/1.
    - Missing matrix values cause the run to stop unless --missing-as-zero is used.

Outputs:
    - collapsed matrix
    - collapsed metadata file
    - matrix summary file
    - collapse summary file
    - collapse groups file
    - unmapped features report
    - missing values report, if missing values are detected
    - non-binary values report, if non-binary values are detected
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import csv
import sys

import pandas as pd

# Helpers
def sniff_delimiter(path: Path) -> str:
    """
    Determine whether a file is comma- or tab-delimited
    Defaults to comma if the delimiter cannot be confidently detected
    """
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)

    if not sample.strip():
        raise ValueError(f"Input file appears to be empty: {path}")

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        # Conservative fallback: whichever delimiter appears more often in the header.
        first_line = sample.splitlines()[0]
        if first_line.count("\t") > first_line.count(","):
            return "\t"
        return ","
    
def get_header(path: Path, delimiter: str) -> List[str]:
    """
    Read the header row from a delimited file
    """
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Input file appears to be empty: {path}")

    header = [h.strip() for h in header]
    if not header:
        raise ValueError(f"No header detected in file: {path}")

    return header

def find_duplicates(values: List[str]) -> List[str]:
    """
    Return duplicated values while preserving first duplicate encounter order
    """
    seen = set()
    duplicates = []
    duplicate_seen = set()

    for value in values:
        if value in seen and value not in duplicate_seen:
            duplicates.append(value)
            duplicate_seen.add(value)
        seen.add(value)

    return duplicates

def output_paths(output_file: Path) -> Dict[str, Path]:
    """
    Construct auxiliary output paths from the requested collapsed matrix path

    Example:
        collapsed_matrix.csv
        collapsed_matrix.metadata.csv
        collapsed_matrix.matrix_summary.csv
    """
    parent = output_file.parent
    stem = output_file.stem

    return {
        "matrix": output_file,
        "metadata": parent / f"{stem}.metadata.csv",
        "matrix_summary": parent / f"{stem}.matrix_summary.csv",
        "collapse_summary": parent / f"{stem}.collapse_summary.csv",
        "collapse_groups": parent / f"{stem}.collapse_groups.csv",
        "unmapped_features": parent / f"{stem}.unmapped_features.csv",
        "missing_values": parent / f"{stem}.missing_values.csv",
        "non_binary_values": parent / f"{stem}.non_binary_values.csv",
        "mapping_absent": parent / f"{stem}.mapped_features_absent_from_matrix.csv",
        "blank_bin_features": parent / f"{stem}.blank_bin_features.csv",
        "features_present_all": parent / f"{stem}.features_present_in_all_genomes.csv",
    }

def fail(message: str, exit_code: int = 1) -> None:
    """
    Print a clear error message and exit
    """
    print(f"\nERROR: {message}\n", file=sys.stderr)
    raise SystemExit(exit_code)

# Read inputs
def read_matrix(matrix_file: Path, id_col: Optional[str] = None) -> Tuple[pd.DataFrame, str]:
    """
    Read a presence/absence matrix

    If id_col is None, the first column is used as the genome/sample ID column
    Duplicate feature names are detected from the raw header before pandas can
    modify them.
    """
    delimiter = sniff_delimiter(matrix_file)
    header = get_header(matrix_file, delimiter)

    if len(header) < 2:
        raise ValueError(
            "Presence/absence matrix must contain at least one genome/sample ID "
            "column and one feature column."
        )

    duplicate_headers = find_duplicates(header)
    if duplicate_headers:
        dup_str = ", ".join(duplicate_headers)
        raise ValueError(
            f"Duplicate column names detected in matrix header: {dup_str}. "
            "Please resolve duplicate feature names before collapsing features"
        )

    if id_col is None:
        id_col = header[0]
    elif id_col not in header:
        raise ValueError(
            f"Requested ID column '{id_col}' was not found in matrix file. "
            f"Available columns include: {', '.join(header[:10])}"
        )

    matrix_df = pd.read_csv(matrix_file, delimiter=delimiter, dtype=str, keep_default_na=True)

    # Strip leading/trailing whitespace from column names and string values.
    matrix_df.columns = [str(c).strip() for c in matrix_df.columns]

    if matrix_df.empty:
        raise ValueError(f"Matrix file contains a header but no data rows: {matrix_file}")

    return matrix_df, id_col

def read_mapping_metadata(metadata_file: Path) -> pd.DataFrame:
    """
    Read a metadata/mapping file. Extra columns are allowed
    """
    delimiter = sniff_delimiter(metadata_file)
    header = get_header(metadata_file, delimiter)

    duplicate_headers = find_duplicates(header)
    if duplicate_headers:
        dup_str = ", ".join(duplicate_headers)
        raise ValueError(
            f"Duplicate column names detected in metadata header: {dup_str}"
        )

    mapping_df = pd.read_csv(metadata_file, delimiter=delimiter, dtype=str, keep_default_na=True)
    mapping_df.columns = [str(c).strip() for c in mapping_df.columns]

    if mapping_df.empty:
        raise ValueError(f"Metadata file contains a header but no data rows: {metadata_file}")

    return mapping_df

# Validation
def validate_matrix_ids(matrix_df: pd.DataFrame, id_col: str) -> None:
    """
    Validate genome/sample IDs
    """
    if id_col not in matrix_df.columns:
        raise ValueError(f"ID column '{id_col}' not found in matrix.")

    id_series = matrix_df[id_col]

    missing_ids = id_series.isna() | id_series.astype(str).str.strip().eq("")
    if missing_ids.any():
        n_missing = int(missing_ids.sum())
        raise ValueError(
            f"Matrix contains {n_missing} missing/blank genome IDs in column '{id_col}'"
        )

    duplicated_ids = id_series[id_series.duplicated()].dropna().astype(str).unique().tolist()
    if duplicated_ids:
        preview = ", ".join(duplicated_ids[:10])
        extra = "" if len(duplicated_ids) <= 10 else f" ... plus {len(duplicated_ids) - 10} more"
        raise ValueError(
            f"Duplicate genome/sample IDs detected in column '{id_col}': {preview}{extra}"
        )
    
def validate_mapping_metadata(
    mapping_df: pd.DataFrame,
    gene_col: str = "Gene",
    bin_col: str = "bin",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Validate and clean mapping metadata

    Rules:
        - gene_col and bin_col must exist
        - Gene values cannot be missing
        - Blank bin values are allowed and treated as "unmapped"
        - Only rows with nonblank bin values are used for collapsing
        - One Gene cannot map to multiple different nonblank bins
        - Duplicate identical Gene-bin rows are allowed and de-duplicated

    Returns:
        active_mapping:
            Rows with nonblank bin values. These are used for collapsing

        blank_bin_rows:
            Rows where Gene is present but bin is blank. These are reported
            as intentionally/unintentionally unmapped depending on user review
    """
    if gene_col not in mapping_df.columns:
        raise ValueError(
            f"Required mapping column '{gene_col}' not found in metadata file"
        )

    if bin_col not in mapping_df.columns:
        raise ValueError(
            f"Required mapping column '{bin_col}' not found in metadata file"
        )

    cleaned = mapping_df[[gene_col, bin_col]].copy()

    cleaned[gene_col] = cleaned[gene_col].astype(str).str.strip()
    cleaned[bin_col] = cleaned[bin_col].astype(str).str.strip()

    missing_like = {"", "nan", "None", "NA", "N/A", "na", "n/a"}

    gene_missing = cleaned[gene_col].isin(missing_like)
    if gene_missing.any():
        n_missing = int(gene_missing.sum())
        raise ValueError(
            f"Mapping metadata contains {n_missing} missing/blank values in '{gene_col}'"
        )

    blank_bin_mask = cleaned[bin_col].isin(missing_like)

    blank_bin_rows = cleaned.loc[blank_bin_mask].copy()
    blank_bin_rows = blank_bin_rows.rename(columns={
        gene_col: "Gene",
        bin_col: "bin"
    })
    blank_bin_rows["status"] = "blank_bin_treated_as_unmapped"

    active_mapping = cleaned.loc[~blank_bin_mask].copy()
    active_mapping = active_mapping.drop_duplicates()

    # Detect conflicting mappings: one Gene -> multiple nonblank bins.
    n_bins_per_gene = active_mapping.groupby(gene_col)[bin_col].nunique()
    conflicting = n_bins_per_gene[n_bins_per_gene > 1]

    if not conflicting.empty:
        conflict_genes = conflicting.index.tolist()
        preview = ", ".join(conflict_genes[:10])
        extra = "" if len(conflict_genes) <= 10 else f" ... plus {len(conflict_genes) - 10} more"
        raise ValueError(
            f"One or more genes map to multiple bins: {preview}{extra}. "
            "Each Gene must map to no more than one nonblank bin"
        )

    return active_mapping, blank_bin_rows


def collect_missing_values(
    feature_df: pd.DataFrame,
    id_series: pd.Series,
) -> pd.DataFrame:
    """
    Return a long-format report of missing feature values
    """
    records = []

    for row_idx, genome_id in enumerate(id_series.tolist()):
        for feature in feature_df.columns:
            value = feature_df.iat[row_idx, feature_df.columns.get_loc(feature)]

            is_missing = pd.isna(value) or str(value).strip() == ""
            if is_missing:
                records.append({
                    "data_row_number": row_idx + 2,  # +2 accounts for header and 1-based line numbers
                    "genome": genome_id,
                    "feature": feature,
                    "value": "" if pd.isna(value) else value,
                })

    return pd.DataFrame(records)

def coerce_binary_matrix(
    feature_df: pd.DataFrame,
    id_series: pd.Series,
    missing_as_zero: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Coerce feature columns to binary numeric 0/1

    Allowed values:
        0, 1, 0.0, 1.0, "0", "1", "0.0", "1.0"

    Missing values:
        - reported
        - filled with 0 only if missing_as_zero=True
        - otherwise trigger a stop in run()

    Returns:
        coerced_df
        missing_report_df
        non_binary_report_df
    """
    missing_report = collect_missing_values(feature_df, id_series)

    working = feature_df.copy()

    if not missing_report.empty and missing_as_zero:
        working = working.fillna("0")
        working = working.replace(r"^\s*$", "0", regex=True)

    allowed = {"0", "1", "0.0", "1.0"}
    non_binary_records = []

    coerced = pd.DataFrame(index=working.index)

    for feature in working.columns:
        values = working[feature]

        cleaned_values = []
        for row_idx, value in enumerate(values.tolist()):
            if pd.isna(value) or str(value).strip() == "":
                # Leave as missing for now. run() will decide whether this is fatal.
                cleaned_values.append(value)
                continue

            value_str = str(value).strip()

            if value_str not in allowed:
                non_binary_records.append({
                    "data_row_number": row_idx + 2,
                    "genome": id_series.iloc[row_idx],
                    "feature": feature,
                    "value": value,
                })
                cleaned_values.append(value)
            else:
                # Convert 0.0/1.0 to int-compatible values.
                cleaned_values.append(int(float(value_str)))

        coerced[feature] = cleaned_values

    non_binary_report = pd.DataFrame(non_binary_records)

    # Only convert to int if no missing values remain and no non-binary values exist.
    if missing_report.empty or missing_as_zero:
        if non_binary_report.empty:
            coerced = coerced.astype(int)

    return coerced, missing_report, non_binary_report

# Mapping and collapsing
def make_feature_mapping(
    matrix_features: List[str],
    mapping_df: pd.DataFrame,
    gene_col: str = "Gene",
    bin_col: str = "bin",
    keep_unmapped: bool = True,
) -> Tuple[Dict[str, str], pd.DataFrame, pd.DataFrame]:
    """
    Build a feature mapping dictionary from original matrix features to final
    output features.

    Returns:
        feature_map:
            Dict of original feature -> output feature

        unmapped_df:
            Features present in the matrix but absent from the mapping file

        absent_mapped_df:
            Features present in the mapping file but absent from the matrix
    """
    matrix_feature_set = set(matrix_features)

    raw_map = dict(zip(mapping_df[gene_col], mapping_df[bin_col]))

    feature_map = {}
    unmapped_records = []

    for feature in matrix_features:
        if feature in raw_map:
            feature_map[feature] = raw_map[feature]
        else:
            action = "kept_unchanged" if keep_unmapped else "dropped"
            unmapped_records.append({
                "Gene": feature,
                "action": action,
            })

            if keep_unmapped:
                feature_map[feature] = feature

    mapped_genes = set(raw_map.keys())
    absent_mapped = sorted(mapped_genes - matrix_feature_set)

    absent_mapped_df = pd.DataFrame({
        "Gene": absent_mapped,
        "bin": [raw_map[g] for g in absent_mapped],
        "status": ["present_in_metadata_absent_from_matrix"] * len(absent_mapped),
    })

    unmapped_df = pd.DataFrame(unmapped_records)

    return feature_map, unmapped_df, absent_mapped_df

def ordered_unique(values: List[str]) -> List[str]:
    """
    Return unique values preserving first-observed order
    """
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            out.append(value)
            seen.add(value)
    return out

def collapse_features_by_mapping(
    matrix_df: pd.DataFrame,
    binary_feature_df: pd.DataFrame,
    feature_map: Dict[str, str],
    id_col: str,
) -> pd.DataFrame:
    """
    Collapse matrix feature columns according to feature_map.

    If multiple original features map to the same output feature, the collapsed
    value is max(original values), preserving binary 0/1 behavior.
    """
    if not feature_map:
        raise ValueError(
            "No features remain after mapping. If you used --drop-unmapped, "
            "check that your metadata Gene values match matrix column names."
        )

    id_series = matrix_df[id_col].copy()

    selected_features = list(feature_map.keys())
    feature_df = binary_feature_df[selected_features].copy()

    renamed = feature_df.rename(columns=feature_map)

    # Collapse duplicate output feature columns using max, then restore original
    # first-observed output feature order.
    collapsed_features = renamed.T.groupby(level=0).max().T

    output_feature_order = ordered_unique([feature_map[f] for f in selected_features])
    collapsed_features = collapsed_features.reindex(columns=output_feature_order)

    collapsed_df = pd.concat([id_series, collapsed_features], axis=1)

    return collapsed_df

# Summarize and reports
def summarize_matrix(
    matrix_df: pd.DataFrame,
    binary_feature_df: Optional[pd.DataFrame],
    id_col: str,
    n_missing_values: int = 0,
    n_non_binary_values: int = 0,
) -> pd.DataFrame:
    """
    Create a matrix summary table
    """
    feature_cols = [c for c in matrix_df.columns if c != id_col]
    n_genomes = matrix_df.shape[0]
    n_features = len(feature_cols)

    duplicated_genome_ids = int(matrix_df[id_col].duplicated().sum())
    duplicated_feature_names = len(find_duplicates(feature_cols))

    summary = [
        {"metric": "genome_id_column", "value": id_col},
        {"metric": "n_genomes", "value": n_genomes},
        {"metric": "n_features", "value": n_features},
        {"metric": "duplicated_genome_ids", "value": duplicated_genome_ids},
        {"metric": "duplicated_feature_names", "value": duplicated_feature_names},
        {"metric": "missing_values", "value": n_missing_values},
        {"metric": "non_binary_values", "value": n_non_binary_values},
    ]

    if binary_feature_df is not None and n_non_binary_values == 0:
        # Only calculate feature-per-genome statistics if binary data are valid
        # If missing values remain, this calculation may fail; handle nicely 
        try:
            counts = binary_feature_df.astype(int).sum(axis=1)

            summary.extend([
                {"metric": "min_features_per_genome", "value": float(counts.min())},
                {"metric": "max_features_per_genome", "value": float(counts.max())},
                {"metric": "mean_features_per_genome", "value": float(round(counts.mean(), 2))},
                {"metric": "median_features_per_genome", "value": float(round(counts.median(), 2))},
                {"metric": "sd_features_per_genome", "value": float(round(counts.std(ddof=1), 2)) if len(counts) > 1 else 0.0},
            ])
        except Exception:
            summary.extend([
                {"metric": "min_features_per_genome", "value": "not_calculated"},
                {"metric": "max_features_per_genome", "value": "not_calculated"},
                {"metric": "mean_features_per_genome", "value": "not_calculated"},
                {"metric": "median_features_per_genome", "value": "not_calculated"},
                {"metric": "sd_features_per_genome", "value": "not_calculated"},
            ])

    return pd.DataFrame(summary)

def build_collapse_groups(feature_map: Dict[str, str]) -> pd.DataFrame:
    """
    Build one row per collapsed/final feature, including source features
    """
    grouped: Dict[str, List[str]] = {}

    for source_feature, final_feature in feature_map.items():
        grouped.setdefault(final_feature, []).append(source_feature)

    records = []

    for final_feature, source_features in grouped.items():
        n_source = len(source_features)

        if n_source > 1:
            collapse_status = "collapsed"
        elif final_feature == source_features[0]:
            collapse_status = "unchanged_unmapped_or_identity_mapped"
        else:
            collapse_status = "renamed"

        records.append({
            "Gene": final_feature,
            "n_source_features": n_source,
            "source_features": "; ".join(source_features),
            "collapse_status": collapse_status,
        })

    return pd.DataFrame(records)

def summarize_collapse(
    original_features: List[str],
    collapsed_features: List[str],
    feature_map: Dict[str, str],
    unmapped_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create a high-level collapse summary table
    """
    collapse_groups = build_collapse_groups(feature_map)

    if collapse_groups.empty:
        largest_group = "NA"
        largest_group_size = 0
    else:
        largest_idx = collapse_groups["n_source_features"].astype(int).idxmax()
        largest_group = collapse_groups.loc[largest_idx, "Gene"]
        largest_group_size = int(collapse_groups.loc[largest_idx, "n_source_features"])

    n_original_features = len(original_features)
    n_collapsed_features = len(collapsed_features)
    n_features_unmapped = 0 if unmapped_df.empty else unmapped_df.shape[0]
    n_features_mapped = len(feature_map) - n_features_unmapped

    # Number of final groups that represent more than one source feature.
    n_multi_feature_bins = int((collapse_groups["n_source_features"].astype(int) > 1).sum()) if not collapse_groups.empty else 0

    summary = [
        {"metric": "n_original_features", "value": n_original_features},
        {"metric": "n_collapsed_features", "value": n_collapsed_features},
        {"metric": "n_features_mapped", "value": n_features_mapped},
        {"metric": "n_features_unmapped", "value": n_features_unmapped},
        {"metric": "n_multi_feature_bins", "value": n_multi_feature_bins},
        {"metric": "largest_collapsed_group", "value": largest_group},
        {"metric": "largest_collapsed_group_size", "value": largest_group_size},
    ]

    return pd.DataFrame(summary)

def find_features_present_in_all_genomes(
    collapsed_df: pd.DataFrame,
    id_col: str,
) -> pd.DataFrame:
    """
    Identify collapsed features that are present in every genome/sample

    These features have no variance and may cause issues or provide no
    information for downstream Bayesian network structure learning
    """
    feature_cols = [c for c in collapsed_df.columns if c != id_col]
    n_genomes = collapsed_df.shape[0]

    records = []

    for feature in feature_cols:
        feature_sum = int(collapsed_df[feature].astype(int).sum())

        if feature_sum == n_genomes:
            records.append({
                "feature": feature,
                "n_genomes": n_genomes,
                "n_present": feature_sum,
                "percent_present": 100.0,
                "status": "present_in_all_genomes",
            })

    return pd.DataFrame(records)

def build_collapsed_metadata(feature_map: Dict[str, str]) -> pd.DataFrame:
    """
    Build metadata compatible with downstream ReGAIN modules

    The primary column is Gene, matching existing ReGAIN metadata conventions
    """
    return build_collapse_groups(feature_map)

def write_csv(df: pd.DataFrame, path: Path) -> None:
    """
    Write dataframe to CSV
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

#-------------
# Main entry point
#-------------

def run(args: argparse.Namespace) -> None:
    """
    Entry point used by ReGAIN cli.py.
    """
    input_file = Path(args.input)
    metadata_file = Path(args.metadata)
    output_file = Path(args.output_file)

    id_col = getattr(args, "id_col", None)
    gene_col = getattr(args, "gene_col", "Gene")
    bin_col = getattr(args, "bin_col", "bin")
    drop_unmapped = bool(getattr(args, "drop_unmapped", False))
    missing_as_zero = bool(getattr(args, "missing_as_zero", False))

    paths = output_paths(output_file)

    if not input_file.exists():
        fail(f"Input matrix file not found: {input_file}")

    if not metadata_file.exists():
        fail(f"Metadata/mapping file not found: {metadata_file}")

    try:
        matrix_df, id_col = read_matrix(input_file, id_col=id_col)
        validate_matrix_ids(matrix_df, id_col)

        mapping_df = read_mapping_metadata(metadata_file)
        mapping_df, blank_bin_rows = validate_mapping_metadata(
            mapping_df,
            gene_col=gene_col,
            bin_col=bin_col,
        )

        if not blank_bin_rows.empty:
            write_csv(blank_bin_rows, paths["blank_bin_features"])

        original_features = [c for c in matrix_df.columns if c != id_col]

        if not original_features:
            raise ValueError(
                "No feature columns found in input matrix after excluding the ID column."
            )

        feature_df = matrix_df[original_features].copy()
        id_series = matrix_df[id_col].copy()

        binary_feature_df, missing_report, non_binary_report = coerce_binary_matrix(
            feature_df=feature_df,
            id_series=id_series,
            missing_as_zero=missing_as_zero,
        )

        # Always write missing/non-binary reports when issues are detected.
        if not missing_report.empty:
            write_csv(missing_report, paths["missing_values"])

        if not non_binary_report.empty:
            write_csv(non_binary_report, paths["non_binary_values"])

        # Write matrix summary before potentially stopping.
        matrix_summary = summarize_matrix(
            matrix_df=matrix_df,
            binary_feature_df=binary_feature_df if non_binary_report.empty else None,
            id_col=id_col,
            n_missing_values=0 if missing_report.empty else missing_report.shape[0],
            n_non_binary_values=0 if non_binary_report.empty else non_binary_report.shape[0],
        )
        write_csv(matrix_summary, paths["matrix_summary"])

        if not missing_report.empty and not missing_as_zero:
            fail(
                f"Missing values were detected in the input matrix. "
                f"A missing-value report was written to: {paths['missing_values']}\n"
                "Re-run with --missing-as-zero if you want missing values treated as absence.\n"
                "WARNING: Only fill missing values with 0 if they are truly absent; "
                "replacing missing values with 0 will affect the network structure and "
                "result in inaccurate results"
            )

        if not non_binary_report.empty:
            fail(
                f"Non-binary values were detected in the input matrix. "
                f"A non-binary value report was written to: {paths['non_binary_values']}\n"
                "Only 0/1-style values are accepted for feature columns."
            )

        feature_map, unmapped_df, absent_mapped_df = make_feature_mapping(
            matrix_features=original_features,
            mapping_df=mapping_df,
            gene_col=gene_col,
            bin_col=bin_col,
            keep_unmapped=not drop_unmapped,
        )

        # Reports that are useful regardless of whether unmapped features exist.
        write_csv(unmapped_df, paths["unmapped_features"])
        write_csv(absent_mapped_df, paths["mapping_absent"])

        collapsed_df = collapse_features_by_mapping(
            matrix_df=matrix_df,
            binary_feature_df=binary_feature_df,
            feature_map=feature_map,
            id_col=id_col,
        )

        features_present_all_df = find_features_present_in_all_genomes(
            collapsed_df=collapsed_df,
            id_col=id_col
        )

        if not features_present_all_df.empty:
            write_csv(features_present_all_df, paths["features_present_all"])

        collapsed_features = [c for c in collapsed_df.columns if c != id_col]

        collapse_groups = build_collapse_groups(feature_map)
        collapse_summary = summarize_collapse(
            original_features=original_features,
            collapsed_features=collapsed_features,
            feature_map=feature_map,
            unmapped_df=unmapped_df,
        )

        if not features_present_all_df.empty:
            collapse_summary = pd.concat(
                [
                    collapse_summary,
                    pd.DataFrame([
                        {
                            "metric": "n_features_present_in_all_genomes",
                            "value": features_present_all_df.shape[0],
                        },
                        {
                            "metric": "features_present_in_all_genomes",
                            "value": "; ".join(features_present_all_df["feature"].tolist()),
                        },
                    ])
                ],
                ignore_index=True,
            )

        collapsed_metadata = build_collapsed_metadata(feature_map)

        write_csv(collapsed_df, paths["matrix"])
        write_csv(collapsed_metadata, paths["metadata"])
        write_csv(collapse_summary, paths["collapse_summary"])
        write_csv(collapse_groups, paths["collapse_groups"])

        print("\nReGAIN collapse-features complete.")
        print(f"Collapsed matrix:          {paths['matrix']}")
        print(f"Collapsed metadata:        {paths['metadata']}")
        print(f"Matrix summary:            {paths['matrix_summary']}")
        print(f"Collapse summary:          {paths['collapse_summary']}")
        print(f"Collapse groups:           {paths['collapse_groups']}")
        print(f"Unmapped features report:  {paths['unmapped_features']}")
        print(f"Absent mapped features:    {paths['mapping_absent']}")

        if not blank_bin_rows.empty:
            print(f"Blank bin features: {paths['blank_bin_features']}")

        if not missing_report.empty:
            print(f"Missing values report:     {paths['missing_values']}")

        if not features_present_all_df.empty:
            print(
                f"\n[WARNING] {features_present_all_df.shape[0]} collapsed feature(s) are present "
                "in all genomes and have no variance for Bayesian network structure learning"
            )
            print(f"Features present in all genomes: {paths['features_present_all']}")

        if drop_unmapped:
            print("\nNote: --drop-unmapped was used. Unmapped matrix features were excluded.")
        else:
            print("\nNote: Unmapped matrix features were kept unchanged by default.")

    except SystemExit:
        raise
    except Exception as exc:
        fail(str(exc))