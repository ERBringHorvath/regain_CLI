#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 University of Utah

"""
matrix_summary.py

ReGAIN utility for summarizing and validating presence/absence matrices.

Expected CLI usage through cli.py:

    regain matrix-summary \
        -i presence_absence.csv \
        -o presence_absence_summary.csv

Default behavior:
    - Uses the first column as the genome/sample ID column unless --id-col is supplied
    - Accepts CSV or TSV input
    - Validates binary feature columns
    - Reports missing values
    - Stops if missing values are detected unless --missing-as-zero is supplied
    - Stops if non-binary feature values are detected
    - Reports invariant features
    - Reports identical resistance/accessory profiles across genomes

Outputs:
    - matrix summary CSV
    - profile groups CSV
    - profile members CSV
    - missing values report, if missing values are detected
    - non-binary values report, if non-binary values are detected
    - features present in all genomes report, if applicable
    - features absent in all genomes report, if applicable
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import csv
import hashlib
import sys

import pandas as pd

# Helpers
def sniff_delimiter(path: Path) -> str:
    """
    Infer whether a file is comma- or tab-delimited
    Defaults to comma (expected ReGAIN matrix output) if the
        delimiter cannot be confidently called
    """
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)

    if not sample.strip():
        raise ValueError(f"Input file appears to be empty: {path}")

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
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
    Construct auxiliary output paths from the requested summary output path

    Example (-o matrix_summary.csv):
        matrix_summary.csv
        matrix_summary.profile_groups.csv
        matrix_summary.profile_members.csv
    """
    parent = output_file.parent
    stem = output_file.stem

    return {
        "summary": output_file,
        "profile_groups": parent / f"{stem}.profile_groups.csv",
        "profile_members": parent / f"{stem}.profile_members.csv",
        "feature_counts": parent / f"{stem}.feature_counts.csv",
        "missing_values": parent / f"{stem}.missing_values.csv",
        "non_binary_values": parent / f"{stem}.non_binary_values.csv",
        "features_present_all": parent / f"{stem}.features_present_in_all_genomes.csv",
        "features_absent_all": parent / f"{stem}.features_absent_in_all_genomes.csv",
    }

def write_csv(df: pd.DataFrame, path: Path) -> None:
    """
    Write output dataframe to CSV
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

def fail(message: str, exit_code: int = 1) -> None:
    """
    Print a clear error message and exit
    """
    print(f"ERROR: {message}\n", file=sys.stderr)
    raise SystemExit(exit_code)

# Read and validate
def read_matrix(matrix_file: Path, id_col: Optional[str] = None) -> Tuple[pd.DataFrame, str]:
    """
    Read a presence/absence matrix

    if id_col is None, the first column is assumed to be the genome/sample ID column
    Duplicate feature names are detected from the raw header before pandas can modify them
    """
    delimiter = sniff_delimiter(matrix_file)
    header = get_header(matrix_file, delimiter)

    if len(header) < 2:
        raise ValueError(
            "Presence/absence matrix must contain at least one genome/sample ID "
            "column and one feature column"
        )

    duplicate_headers = find_duplicates(header)
    if duplicate_headers:
        dup_str = ", ".join(duplicate_headers)
        raise ValueError(
            f"Duplicate column names detected in matrix header: {dup_str}. "
            "Please resolve duplicate feature names before running matrix-summary"
        )

    if id_col is None:
        id_col = header[0]
    elif id_col not in header:
        raise ValueError(
            f"Requested ID column '{id_col}' was not found in matrix file. "
            f"Available columns include: {', '.join(header[:10])}"
        )

    matrix_df = pd.read_csv(
        matrix_file,
        delimiter=delimiter,
        dtype=str,
        keep_default_na=True,
    )

    matrix_df.columns = [str(c).strip() for c in matrix_df.columns]

    if matrix_df.empty:
        raise ValueError(f"Matrix file contains a header but no data rows: {matrix_file}")

    return matrix_df, id_col

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
            f"Duplicate genome/sample IDs detected in column '{id_col}': {preview}{extra}\n"
        )
    
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
                    "data_row_number": row_idx + 2,
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
        0, 1, 1.0, 0.0, "1", "0", "1.0", "0.0"

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
                cleaned_values.append(int(float(value_str)))

        coerced[feature] = cleaned_values

    non_binary_report = pd.DataFrame(non_binary_records)

    if missing_report.empty or missing_as_zero:
        if non_binary_report.empty:
            coerced = coerced.astype(int)

    return coerced, missing_report, non_binary_report

#----- Matrix-level summaries -----

def summarize_matrix(
    matrix_df: pd.DataFrame,
    binary_feature_df: Optional[pd.DataFrame],
    id_col: str,
    n_missing_values: int = 0,
    n_non_binary_values: int = 0,
    profile_summary: Optional[Dict[str, object]] = None,
    n_features_present_all: int = 0,
    n_features_absent_all: int = 0,
) -> pd.DataFrame:
    """
    Create the matrix summary table
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
        try:
            counts = binary_feature_df.astype(int).sum(axis=1)

            summary.extend([
                {"metric": "min_features_per_genome", "value": float(counts.min())},
                {"metric": "max_features_per_genome", "value": float(counts.max())},
                {"metric": "mean_features_per_genome", "value": float(round(counts.mean(), 2))},
                {"metric": "median_features_per_genome", "value": float(counts.median())},
                {
                    "metric": "sd_features_per_genome",
                    "value": float(round(counts.std(ddof=1), 2)) if len(counts) > 1 else 0.0,
                },
            ])
        except Exception:
            summary.extend([
                {"metric": "min_features_per_genome", "value": "not_calculated"},
                {"metric": "max_features_per_genome", "value": "not_calculated"},
                {"metric": "mean_features_per_genome", "value": "not_calculated"},
                {"metric": "median_features_per_genome", "value": "not_calculated"},
                {"metric": "sd_features_per_genome", "value": "not_calculated"},
            ])

    summary.extend([
        {"metric": "n_features_present_in_all_genomes", "value": n_features_present_all},
        {"metric": "n_features_absent_in_all_genomes", "value": n_features_absent_all},
    ])

    if profile_summary is not None:
        for metric, value in profile_summary.items():
            summary.append({"metric": metric, "value": value})

    return pd.DataFrame(summary)

def find_features_present_in_all_genomes(
    binary_feature_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Identify ubiquitous features
    """
    n_genomes = binary_feature_df.shape[0]
    records = []

    for feature in binary_feature_df.columns:
        feature_sum = int(binary_feature_df[feature].astype(int).sum())

        if feature_sum == n_genomes:
            records.append({
                "feature": feature,
                "n_genomes": n_genomes,
                "n_present": feature_sum,
                "percent_present": 100.0,
                "status": "present_in_all_genomes",
            })

    return pd.DataFrame(records)

def find_features_absent_in_all_genomes(
    binary_feature_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Identify features that are absent in every genome/sample
    """
    n_genomes = binary_feature_df.shape[0]
    records = []

    for feature in binary_feature_df.columns:
        feature_sum = int(binary_feature_df[feature].astype(int).sum())

        if feature_sum == 0:
            records.append({
                "feature": feature,
                "n_genomes": n_genomes,
                "n_present": feature_sum,
                "percent_present": 0.0,
                "status": "absent_in_all_genomes",
            })

    return pd.DataFrame(records)

# Resistance/accessory profile summaries

def make_profile_hash(values: List[int]) -> str:
    """
    Create a stable short hash for a binary feature vector
    """
    vector_string = "".join(str(int(v)) for v in values)
    return hashlib.sha256(vector_string.encode("utf-8")).hexdigest()[:16]

def build_profile_members(
    binary_feature_df: pd.DataFrame,
    id_series: pd.Series,
) -> pd.DataFrame:
    """
    Assign each genome/sample to an identical profile group based on its complete
    binary feature vector

    Profile IDs are assigned in descending order of group size, then by first
    observed order for stability
    """
    records = []

    for idx, genome in enumerate(id_series.tolist()):
        values = binary_feature_df.iloc[idx].astype(int).tolist()
        profile_hash = make_profile_hash(values)
        n_features_present = int(sum(values))

        records.append({
            "genome": genome,
            "profile_hash": profile_hash,
            "n_features_present": n_features_present,
        })

    raw_members = pd.DataFrame(records)

    # Summarize hash group sizes.
    group_sizes = (
        raw_members
        .groupby("profile_hash", sort=False)
        .size()
        .reset_index(name="n_genomes")
    )

    # Track first observed order for stable tie-breaking.
    first_order = (
        raw_members
        .reset_index()
        .groupby("profile_hash", sort=False)["index"]
        .min()
        .reset_index(name="first_observed_index")
    )

    group_sizes = group_sizes.merge(first_order, on="profile_hash", how="left")
    group_sizes = group_sizes.sort_values(
        by=["n_genomes", "first_observed_index"],
        ascending=[False, True],
    ).reset_index(drop=True)

    profile_id_map = {
        row["profile_hash"]: f"profile_{i + 1:05d}"
        for i, row in group_sizes.iterrows()
    }

    raw_members["profile_id"] = raw_members["profile_hash"].map(profile_id_map)

    # Add profile size/status to each member.
    size_map = dict(zip(group_sizes["profile_hash"], group_sizes["n_genomes"]))
    raw_members["profile_group_size"] = raw_members["profile_hash"].map(size_map)
    raw_members["status"] = raw_members["profile_group_size"].apply(
        lambda x: "duplicate_profile" if int(x) > 1 else "unique_profile"
    )

    # Reorder columns.
    members = raw_members[
        [
            "genome",
            "profile_id",
            "profile_hash",
            "profile_group_size",
            "n_features_present",
            "status",
        ]
    ].copy()

    return members

def build_profile_groups(
    profile_members_df: pd.DataFrame,
    binary_feature_df: pd.DataFrame,
    id_series: pd.Series,
) -> pd.DataFrame:
    """
    Summarize identical profile groups

    Adds a semicolon-sep list of features present in each profile
    """
    records = []

    # Build lookup from genome ID to row index in the binary matrix.
    genome_to_index = {
        str(genome): idx
        for idx, genome in enumerate(id_series.astype(str).tolist())
    }

    grouped = profile_members_df.groupby(
        ["profile_id", "profile_hash"],
        sort=False,
    )

    for (profile_id, profile_hash), group in grouped:
        genomes = group["genome"].astype(str).tolist()
        n_genomes = len(genomes)
        n_features_present_values = group["n_features_present"].unique().tolist()

        # All members in the same binary profile should have the same count.
        n_features_present = int(n_features_present_values[0])

        # Use the first genome in the profile group as the representative vector.
        representative_genome = genomes[0]
        representative_idx = genome_to_index[representative_genome]

        present_features = binary_feature_df.columns[
            binary_feature_df.iloc[representative_idx].astype(int) == 1
        ].astype(str).tolist()

        status = "duplicate_profile" if n_genomes > 1 else "unique_profile"

        records.append({
            "profile_id": profile_id,
            "profile_hash": profile_hash,
            "features_present": "; ".join(present_features),
            "n_genomes": n_genomes,
            "genomes": "; ".join(genomes),
            "n_features_present": n_features_present,
            "status": status,
        })

    profile_groups = pd.DataFrame(records)

    if not profile_groups.empty:
        profile_groups = profile_groups.sort_values(
            by=["n_genomes", "profile_id"],
            ascending=[False, True],
        ).reset_index(drop=True)

    return profile_groups

def build_feature_counts(binary_feature_df: pd.DataFrame) -> pd.DataFrame:
    """
    Reports the occurrence of each variable/feature/gene

    Returns one row per feature
    """
    n_genomes = binary_feature_df.shape[0]

    records = []

    for feature in binary_feature_df.columns:
        count = int(binary_feature_df[feature].astype(int).sum())
        percent_present = (
            100.0 * count / n_genomes
            if n_genomes > 0
            else 0.0
        )

        records.append({
            "feature": feature,
            "count": count,
            "percent_present": round(percent_present, 2),
        })

    feature_counts_df = pd.DataFrame(records)

    if not feature_counts_df.empty:
        feature_counts_df = feature_counts_df.sort_values(
            by=["count", "feature"],
            ascending=[False, True],
        ).reset_index(drop=True)

    return feature_counts_df

def summarize_profiles(profile_groups_df: pd.DataFrame) -> Dict[str, object]:
    """
    Return profile-level summary metrics for the main summary table
    """
    if profile_groups_df.empty:
        return {
            "n_unique_profiles": 0,
            "n_duplicate_profile_groups": 0,
            "n_genomes_in_duplicate_profiles": 0,
            "percent_genomes_in_duplicate_profiles": 0.0,
            "largest_profile_group_size": 0,
        }

    n_unique_profiles = profile_groups_df.shape[0]

    duplicate_groups = profile_groups_df[
        profile_groups_df["status"] == "duplicate_profile"
    ]

    n_duplicate_profile_groups = duplicate_groups.shape[0]
    n_genomes_in_duplicate_profiles = int(duplicate_groups["n_genomes"].sum()) if not duplicate_groups.empty else 0
    n_total_genomes = int(profile_groups_df["n_genomes"].sum())
    percent_genomes_in_duplicate_profiles = (
        100.0 * n_genomes_in_duplicate_profiles / n_total_genomes
        if n_total_genomes > 0
        else 0.0
    )
    largest_profile_group_size = int(profile_groups_df["n_genomes"].max())

    return {
        "n_unique_profiles": n_unique_profiles,
        "n_duplicate_profile_groups": n_duplicate_profile_groups,
        "n_genomes_in_duplicate_profiles": n_genomes_in_duplicate_profiles,
        "percent_genomes_in_duplicate_profiles": round(percent_genomes_in_duplicate_profiles, 2),
        "largest_profile_group_size": largest_profile_group_size,
    }

#------------
# Entry point
#------------

def run(args: argparse.Namespace) -> None:
    """
    Entry point used by cli.py
    """
    input_file = Path(args.input)
    output_file = Path(args.output_file)

    id_col = getattr(args, "id_col", None)
    missing_as_zero = bool(getattr(args, "missing_as_zero", False))

    paths = output_paths(output_file)

    if not input_file.exists():
        fail(f"Input matrix file not found: {input_file}\n")

    try:
        matrix_df, id_col = read_matrix(input_file, id_col=id_col)
        validate_matrix_ids(matrix_df, id_col)

        feature_cols = [c for c in matrix_df.columns if c != id_col]

        if not feature_cols:
            raise ValueError(
                "No feature columns found in input matrix after excluding the ID column"
            )

        feature_df = matrix_df[feature_cols].copy()
        id_series = matrix_df[id_col].copy()

        binary_feature_df, missing_report, non_binary_report = coerce_binary_matrix(
            feature_df=feature_df,
            id_series=id_series,
            missing_as_zero=missing_as_zero,
        )

        if not missing_report.empty:
            write_csv(missing_report, paths["missing_values"])

        if not non_binary_report.empty:
            write_csv(non_binary_report, paths["non_binary_values"])

        # Stop early if missing values remain unresolved.
        if not missing_report.empty and not missing_as_zero:
            matrix_summary = summarize_matrix(
                matrix_df=matrix_df,
                binary_feature_df=None,
                id_col=id_col,
                n_missing_values=missing_report.shape[0],
                n_non_binary_values=non_binary_report.shape[0] if not non_binary_report.empty else 0,
            )
            write_csv(matrix_summary, paths["summary"])

            fail(
                f"Missing values were detected in the input matrix. "
                f"A missing-value report was written to: {paths['missing_values']} "
                "Re-run with --missing-as-zero if you want missing values treated as absence "
                "WARNING: Only fill missing values with 0 if they are truly absent; "
                "replacing missing values with 0 will affect downstream network structure"
            )

        # Stop early if non-binary values are present.
        if not non_binary_report.empty:
            matrix_summary = summarize_matrix(
                matrix_df=matrix_df,
                binary_feature_df=None,
                id_col=id_col,
                n_missing_values=missing_report.shape[0] if not missing_report.empty else 0,
                n_non_binary_values=non_binary_report.shape[0],
            )
            write_csv(matrix_summary, paths["summary"])

            fail(
                f"Non-binary values were detected in the input matrix. "
                f"A non-binary value report was written to: {paths['non_binary_values']} "
                "Only 0/1-style values are accepted for feature columns"
            )

        # At this point, binary_feature_df is safe to summarize.
        features_present_all_df = find_features_present_in_all_genomes(binary_feature_df)
        features_absent_all_df = find_features_absent_in_all_genomes(binary_feature_df)

        feature_counts_df = build_feature_counts(binary_feature_df)
        write_csv(feature_counts_df, paths["feature_counts"])

        if not features_present_all_df.empty:
            write_csv(features_present_all_df, paths["features_present_all"])

        if not features_absent_all_df.empty:
            write_csv(features_absent_all_df, paths["features_absent_all"])

        profile_members_df = build_profile_members(
            binary_feature_df=binary_feature_df,
            id_series=id_series,
        )

        profile_groups_df = build_profile_groups(
            profile_members_df=profile_members_df,
            binary_feature_df=binary_feature_df,
            id_series=id_series,
        )
        profile_summary = summarize_profiles(profile_groups_df)

        write_csv(profile_members_df, paths["profile_members"])
        write_csv(profile_groups_df, paths["profile_groups"])

        matrix_summary = summarize_matrix(
            matrix_df=matrix_df,
            binary_feature_df=binary_feature_df,
            id_col=id_col,
            n_missing_values=0,
            n_non_binary_values=0,
            profile_summary=profile_summary,
            n_features_present_all=features_present_all_df.shape[0],
            n_features_absent_all=features_absent_all_df.shape[0],
        )

        write_csv(matrix_summary, paths["summary"])

        print("\nReGAIN matrix-summary complete\n")
        print(f"Matrix summary:            {paths['summary']}")
        print(f"Profile groups:            {paths['profile_groups']}")
        print(f"Profile members:           {paths['profile_members']}")
        print(f"Feature counts:            {paths['feature_counts']}")

        if not missing_report.empty:
            print(f"Missing values report:     {paths['missing_values']}")

        if not non_binary_report.empty:
            print(f"Non-binary values report:  {paths['non_binary_values']}")

        if not features_present_all_df.empty:
            print(
                f"Warning: {features_present_all_df.shape[0]} feature(s) are present "
                "in all genomes and have no variance for Bayesian network structure learning\n"
            )
            print(f"Features present in all genomes: {paths['features_present_all']}\n")

        if not features_absent_all_df.empty:
            print(
                f"Warning: {features_absent_all_df.shape[0]} feature(s) are absent "
                "in all genomes and have no variance for Bayesian network structure learning\n"
            )
            print(f"Features absent in all genomes:  {paths['features_absent_all']}\n")

        if profile_summary["n_duplicate_profile_groups"] > 0:
            print(
                f"Profile redundancy: {profile_summary['n_genomes_in_duplicate_profiles']} \n"
                f"genome(s) occur in {profile_summary['n_duplicate_profile_groups']} \n"
                "duplicate resistance/accessory profile group(s)\n"
            )
            print(
                f"Largest duplicate/identical profile group size: \n"
                f"{profile_summary['largest_profile_group_size']}\n"
            )
        else:
            print("Profile redundancy: no duplicate resistance/accessory profiles detected\n")

        if missing_as_zero and not missing_report.empty:
            print(
                "Note: --missing-as-zero was used. Missing matrix values were treated as absence\n"
            )

    except SystemExit:
        raise
    except Exception as exc:
        fail(str(exc))