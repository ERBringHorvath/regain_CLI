#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 University of Utah

"""
network_analysis.py

ReGAIN utility for comparing two bnS/bnL results tables

Design:
    - Pair matching is undirected: Gene_A/Gene_B and Gene_B/Gene_A are treated
      as belonging to the same gene pair.
    - Directional statistics are preserved separately within each matched pair
    - Inputs must be ReGAIN bnS/bnL results tables

Expected CLI usage through cli.py:

    regain network-analysis \
        --network1 full_results.csv \
        --network2 reduced_results.csv \
        -o network_analysis_out

Default behavior:
    - Match unordered gene pairs across the two result tables
    - Preserve directional means separately:
        Gene_A -> Gene_B
        Gene_B -> Gene_A
    - Compute delta values as:
        network2_mean - network1_mean
    - Classify pair-level status using relative risk mean deltas by default

Outputs:
    network_analysis.comparison.csv
    network_analysis.summary.csv
    network_analysis.network1_only.csv
    network_analysis.network2_only.csv
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import csv
import sys
import statistics

import pandas as pd

#---------
# Defaults
#---------

DEFAULT_GENE_A_COL = "Gene_1"
DEFAULT_GENE_B_COL = "Gene_2"

DEFAULT_CPR_MEAN_COL = "Conditional_Probability_Mean"
DEFAULT_RR_MEAN_COL = "Relative_Risk_Mean"
DEFAULT_ARD_MEAN_COL = "Absolute_Risk_Mean"

DEFAULT_STATUS_METRIC = "ard"
DEFAULT_DELTA_THRESHOLD = 0.1

# Helpers

def fail(message: str, exit_code: int = 1) -> None:
    """
    Print a clear error message and quit
    """
    print(f"\nERROR: {message}\n", file=sys.stderr)
    raise SystemExit(exit_code)

def write_csv(df: pd.DataFrame, path: Path) -> None:
    """
    Write a dataframe to CSV
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

def sniff_delimiter(path: Path) -> str:
    """
    Infer whether a file is comma- or tab-delimited
    Defaults to comma (standard ReGAIN out) if delimiter cannot be 
    confidently detected
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

def output_paths(output_dir: Path) -> Dict[str, Path]:
    """
    Construct output paths for network-analysis
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    return {
        "comparison": output_dir / "network_analysis.comparison.csv",
        "summary": output_dir / "network_analysis.summary.csv",
        "network1_only": output_dir / "network_analysis.network1_only.csv",
        "network2_only": output_dir / "network_analysis.network2_only.csv",
    }

#--------------------
# Read in / normalize
#--------------------

def read_results_table(results_file: Path) -> pd.DataFrame:
    """
    Read a bnS/bnL results table
    """
    if not results_file.exists():
        raise ValueError(f"Results file not found: {results_file}")

    delimiter = sniff_delimiter(results_file)

    df = pd.read_csv(
        results_file,
        delimiter=delimiter,
        dtype=str,
        keep_default_na=True,
    )

    if df.empty:
        raise ValueError(f"Results file contains a header but no data rows: {results_file}")

    df.columns = [str(c).strip() for c in df.columns]

    return df

def validate_required_columns(
        df: pd.DataFrame,
        gene_a_col: str,
        gene_b_col: str,
        cpr_mean_col: str,
        rr_mean_col: str,
        ard_mean_col: str,
        label: str,
) -> None:
    """
    Confirm that the expected ReGAIN results columns are present
    """
    required = [gene_a_col, gene_b_col, cpr_mean_col, rr_mean_col, ard_mean_col]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValueError(
            f"{label} is missing required column(s): {', '.join(missing)}"
        )
    
def coerce_numeric_column(series: pd.Series, column_name: str, label: str) -> pd.Series:
    """
    Converts a results column to numeric and fail if values are non-numeric
    """
    numeric = pd.to_numeric(series, errors="coerce")

    bad_mask = numeric.isna() & ~series.isna() & series.astype(str).str.strip().ne("")
    if bad_mask.any():
        bad_values = series.loc[bad_mask].astype(str).unique().tolist()
        preview = ", ".join(bad_values[:10])
        extra = "" if len(bad_values) <= 10 else f" ... plus {len(bad_values) - 10} more"
        raise ValueError(
            f"{label} contains non-numeric values in '{column_name}': {preview}{extra}"
        )

    return numeric

def normalize_results_table(
    df: pd.DataFrame,
    label: str,
    gene_a_col: str = DEFAULT_GENE_A_COL,
    gene_b_col: str = DEFAULT_GENE_B_COL,
    cpr_mean_col: str = DEFAULT_CPR_MEAN_COL,
    rr_mean_col: str = DEFAULT_RR_MEAN_COL,
    ard_mean_col: str = DEFAULT_ARD_MEAN_COL,
) -> pd.DataFrame:
    """
    Normalize a ReGAIN results table to a standard schema

    The normalized table contains one row per directional relationship:
        gene_from -> gene_to

    It also stores an undirected pair key for matching
    """
    validate_required_columns(
        df=df,
        gene_a_col=gene_a_col,
        gene_b_col=gene_b_col,
        cpr_mean_col=cpr_mean_col,
        rr_mean_col=rr_mean_col,
        ard_mean_col=ard_mean_col,
        label=label,
    )

    working = df.copy()

    working[gene_a_col] = working[gene_a_col].astype(str).str.strip()
    working[gene_b_col] = working[gene_b_col].astype(str).str.strip()

    blank_gene_mask = (
        working[gene_a_col].eq("") |
        working[gene_b_col].eq("") |
        working[gene_a_col].isna() |
        working[gene_b_col].isna()
    )

    if blank_gene_mask.any():
        n_blank = int(blank_gene_mask.sum())
        raise ValueError(
            f"{label} contains {n_blank} row(s) with missing/blank gene names."
        )

    normalized = pd.DataFrame({
        "gene_from": working[gene_a_col],
        "gene_to": working[gene_b_col],
        "cpr_mean": coerce_numeric_column(working[cpr_mean_col], cpr_mean_col, label),
        "rr_mean": coerce_numeric_column(working[rr_mean_col], rr_mean_col, label),
        "ard_mean": coerce_numeric_column(working[ard_mean_col], ard_mean_col, label),
    })

    normalized["pair_gene_1"] = normalized.apply(
        lambda row: sorted([str(row["gene_from"]), str(row["gene_to"])])[0],
        axis=1,
    )
    normalized["pair_gene_2"] = normalized.apply(
        lambda row: sorted([str(row["gene_from"]), str(row["gene_to"])])[1],
        axis=1,
    )
    normalized["pair_key"] = normalized["pair_gene_1"] + "||" + normalized["pair_gene_2"]
    normalized["direction_key"] = normalized["gene_from"] + "->" + normalized["gene_to"]
    normalized["pair"] = normalized["pair_gene_1"] + "--" + normalized["pair_gene_2"]

    # Detect exact duplicate directional rows within the same results file.
    dup_mask = normalized.duplicated(subset=["pair_key", "direction_key"], keep=False)
    if dup_mask.any():
        dup_rows = normalized.loc[dup_mask, ["gene_from", "gene_to"]].drop_duplicates()
        preview_df = dup_rows.head(10)
        preview = "; ".join(
            f"{row.gene_from}->{row.gene_to}"
            for row in preview_df.itertuples(index=False)
        )
        extra = "" if len(dup_rows) <= 10 else f" ... plus {len(dup_rows) - 10} more"
        raise ValueError(
            f"{label} contains duplicate directional rows for the same gene pair: "
            f"{preview}{extra}"
        )

    return normalized

#----------------------
# Pair-level comparison
#----------------------

def directional_row_lookup(pair_df: pd.DataFrame) -> Dict[str, pd.Series]:
    """
    Build a lookup of direction_key -> row for one undirected pair
    """
    lookup = {}
    for _, row in pair_df.iterrows():
        lookup[str(row["direction_key"])] = row
    return lookup

def classify_pair_status(
    comparison_record: Dict[str, object],
    status_metric: str = DEFAULT_STATUS_METRIC,
    delta_threshold: float = DEFAULT_DELTA_THRESHOLD,
) -> str:
    """
    Classify a pair-level status

    Rules:
        - If present only in network1: lost_from_network2
        - If present only in network2: gained_in_network2
        - If shared in both:
            use the maximum  absolute directional delta from the chosen metric
            to classify retained / weakened / strengthened
    """
    n1_present = bool(comparison_record["present_in_network1"])
    n2_present = bool(comparison_record["present_in_network2"])

    if n1_present and not n2_present:
        return "lost_from_network2"

    if n2_present and not n1_present:
        return "gained_in_network2"

    if not n1_present and not n2_present:
        return "not_present"

    metric_prefix = {
        "cpr": "delta_cpr_mean",
        "rr": "delta_rr_mean",
        "ard": "delta_ard_mean",
    }.get(status_metric)

    if metric_prefix is None:
        raise ValueError(
            f"Unsupported status metric '{status_metric}'. Choose from: cpr, rr, ard."
        )

    directional_deltas = []

    for direction_suffix in ["_ab", "_ba"]:
        value = comparison_record.get(metric_prefix + direction_suffix, None)
        if value is not None and pd.notna(value):
            directional_deltas.append(float(value))

    if not directional_deltas:
        return "retained"

    max_abs_delta = max(directional_deltas, key=lambda x: abs(x))

    if max_abs_delta >= delta_threshold:
        return "strengthened"
    if max_abs_delta <= -delta_threshold:
        return "weakened"

    return "retained"

def compare_pair_group(
    pair_key: str,
    n1_pair_df: pd.DataFrame,
    n2_pair_df: pd.DataFrame,
    status_metric: str = DEFAULT_STATUS_METRIC,
    delta_threshold: str = DEFAULT_DELTA_THRESHOLD,
) -> Dict[str, object]:
    """
    Compare one undirected pair across two normalized results tables

    Matching is undirected at the pair level, but directional values are preserved
    separately for:
        pair_gene_1 -> pair_gene_2 (_ab)
        pair_gene_2 -> pair_gene_1 (_ba)
    """
    pair_gene_1, pair_gene_2 = pair_key.split("||")
    edge = f"{pair_gene_1}--{pair_gene_2}"

    n1_lookup = directional_row_lookup(n1_pair_df) if n1_pair_df is not None and not n1_pair_df.empty else {}
    n2_lookup = directional_row_lookup(n2_pair_df) if n2_pair_df is not None and not n2_pair_df.empty else {}

    direction_ab = f"{pair_gene_1}->{pair_gene_2}"
    direction_ba = f"{pair_gene_2}->{pair_gene_1}"

    def extract_metrics(lookup: Dict[str, pd.Series], direction_key: str) -> Dict[str, Optional[float]]:
        if direction_key not in lookup:
            return {
                "cpr_mean": None,
                "rr_mean": None,
                "ard_mean": None,
            }

        row = lookup[direction_key]
        return {
            "cpr_mean": float(row["cpr_mean"]),
            "rr_mean": float(row["rr_mean"]),
            "ard_mean": float(row["ard_mean"]),
        }

    n1_ab = extract_metrics(n1_lookup, direction_ab)
    n1_ba = extract_metrics(n1_lookup, direction_ba)
    n2_ab = extract_metrics(n2_lookup, direction_ab)
    n2_ba = extract_metrics(n2_lookup, direction_ba)

    def delta(v1: Optional[float], v2: Optional[float]) -> Optional[float]:
        if v1 is None or v2 is None:
            return None
        return v2 - v1

    record = {
        "Gene_1": pair_gene_1,
        "Gene_2": pair_gene_2,
        "pair": edge,

        "present_in_network1": bool(n1_lookup),
        "present_in_network2": bool(n2_lookup),
        "pair_presence": (
            "shared" if bool(n1_lookup) and bool(n2_lookup)
            else "network1_only" if bool(n1_lookup)
            else "network2_only" if bool(n2_lookup)
            else "absent"
        ),

        "direction_ab": direction_ab,
        "direction_ba": direction_ba,

        "network1_cpr_mean_ab": n1_ab["cpr_mean"],
        "network2_cpr_mean_ab": n2_ab["cpr_mean"],
        "delta_cpr_mean_ab": delta(n1_ab["cpr_mean"], n2_ab["cpr_mean"]),

        "network1_rr_mean_ab": n1_ab["rr_mean"],
        "network2_rr_mean_ab": n2_ab["rr_mean"],
        "delta_rr_mean_ab": delta(n1_ab["rr_mean"], n2_ab["rr_mean"]),

        "network1_ard_mean_ab": n1_ab["ard_mean"],
        "network2_ard_mean_ab": n2_ab["ard_mean"],
        "delta_ard_mean_ab": delta(n1_ab["ard_mean"], n2_ab["ard_mean"]),

        "network1_cpr_mean_ba": n1_ba["cpr_mean"],
        "network2_cpr_mean_ba": n2_ba["cpr_mean"],
        "delta_cpr_mean_ba": delta(n1_ba["cpr_mean"], n2_ba["cpr_mean"]),

        "network1_rr_mean_ba": n1_ba["rr_mean"],
        "network2_rr_mean_ba": n2_ba["rr_mean"],
        "delta_rr_mean_ba": delta(n1_ba["rr_mean"], n2_ba["rr_mean"]),

        "network1_ard_mean_ba": n1_ba["ard_mean"],
        "network2_ard_mean_ba": n2_ba["ard_mean"],
        "delta_ard_mean_ba": delta(n1_ba["ard_mean"], n2_ba["ard_mean"]),
    }

    record["status"] = classify_pair_status(
        comparison_record=record,
        status_metric=status_metric,
        delta_threshold=delta_threshold,
    )

    return record

def compare_network_tables(
    network1_df: pd.DataFrame,
    network2_df: pd.DataFrame,
    status_metric: str = DEFAULT_STATUS_METRIC,
    delta_threshold: str = DEFAULT_DELTA_THRESHOLD,
) -> pd.DataFrame:
    """
    Compare two normalized network results tables
    """
    pair_keys_1 = set(network1_df["pair_key"].tolist())
    pair_keys_2 = set(network2_df["pair_key"].tolist())
    all_pair_keys = sorted(pair_keys_1 | pair_keys_2)

    records = []

    for pair_key in all_pair_keys:
        n1_pair_df = network1_df.loc[network1_df["pair_key"] == pair_key].copy()
        n2_pair_df = network2_df.loc[network2_df["pair_key"] == pair_key].copy()

        record = compare_pair_group(
            pair_key=pair_key,
            n1_pair_df=n1_pair_df,
            n2_pair_df=n2_pair_df,
            status_metric=status_metric,
            delta_threshold=delta_threshold,
        )
        records.append(record)

    comparison_df = pd.DataFrame(records)

    if not comparison_df.empty:
        comparison_df = comparison_df.sort_values(
            by=["Gene_1", "Gene_2"],
            ascending=[True, True],
        ).reset_index(drop=True)

    return comparison_df

#----------
# Summaries
#----------

def safe_mean(values: List[Optional[float]]) -> Optional[float]:
    """
    Mean over non-nul numeric values
    """
    cleaned = [float(v) for v in values if v is not None and pd.notna(v)]
    if not cleaned:
        return None
    return float(statistics.mean(cleaned))

def build_summary_table(comparison_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a high-level summary of network comparison results
    """
    if comparison_df.empty:
        rows = [
            {"metric": "n_pairs_network1", "value": 0},
            {"metric": "n_pairs_network2", "value": 0},
            {"metric": "n_shared_pairs", "value": 0},
            {"metric": "n_network1_only_pairs", "value": 0},
            {"metric": "n_network2_only_pairs", "value": 0},
            {"metric": "n_retained_pairs", "value": 0},
            {"metric": "n_weakened_pairs", "value": 0},
            {"metric": "n_strengthened_pairs", "value": 0},
            {"metric": "mean_delta_cpr_ab", "value": None},
            {"metric": "mean_delta_cpr_ba", "value": None},
            {"metric": "mean_delta_rr_ab", "value": None},
            {"metric": "mean_delta_rr_ba", "value": None},
            {"metric": "mean_delta_ard_ab", "value": None},
            {"metric": "mean_delta_ard_ba", "value": None},
        ]
        return pd.DataFrame(rows)

    rows = [
        {
            "metric": "n_pairs_network1",
            "value": int(comparison_df["present_in_network1"].sum()),
        },
        {
            "metric": "n_pairs_network2",
            "value": int(comparison_df["present_in_network2"].sum()),
        },
        {
            "metric": "n_shared_pairs",
            "value": int((comparison_df["pair_presence"] == "shared").sum()),
        },
        {
            "metric": "n_network1_only_pairs",
            "value": int((comparison_df["pair_presence"] == "network1_only").sum()),
        },
        {
            "metric": "n_network2_only_pairs",
            "value": int((comparison_df["pair_presence"] == "network2_only").sum()),
        },
        {
            "metric": "n_retained_pairs",
            "value": int((comparison_df["status"] == "retained").sum()),
        },
        {
            "metric": "n_weakened_pairs",
            "value": int((comparison_df["status"] == "weakened").sum()),
        },
        {
            "metric": "n_strengthened_pairs",
            "value": int((comparison_df["status"] == "strengthened").sum()),
        },
        {
            "metric": "mean_delta_cpr_ab",
            "value": safe_mean(comparison_df["delta_cpr_mean_ab"].tolist()),
        },
        {
            "metric": "mean_delta_cpr_ba",
            "value": safe_mean(comparison_df["delta_cpr_mean_ba"].tolist()),
        },
        {
            "metric": "mean_delta_rr_ab",
            "value": safe_mean(comparison_df["delta_rr_mean_ab"].tolist()),
        },
        {
            "metric": "mean_delta_rr_ba",
            "value": safe_mean(comparison_df["delta_rr_mean_ba"].tolist()),
        },
        {
            "metric": "mean_delta_ard_ab",
            "value": safe_mean(comparison_df["delta_ard_mean_ab"].tolist()),
        },
        {
            "metric": "mean_delta_ard_ba",
            "value": safe_mean(comparison_df["delta_ard_mean_ba"].tolist()),
        },
    ]

    return pd.DataFrame(rows)

#-----------------
# Main entry point
#-----------------

def run(args: argparse.Namespace) -> None:
    """
    Entry point for cli.py
    """
    network1_file = Path(args.network1)
    network2_file = Path(args.network2)
    output_dir = Path(args.output_dir)

    gene_a_col = getattr(args, "gene_a_col", DEFAULT_GENE_A_COL)
    gene_b_col = getattr(args, "gene_b_col", DEFAULT_GENE_B_COL)
    cpr_mean_col = getattr(args, "cpr_mean_col", DEFAULT_CPR_MEAN_COL)
    rr_mean_col = getattr(args, "rr_mean_col", DEFAULT_RR_MEAN_COL)
    ard_mean_col = getattr(args, "ard_mean_col", DEFAULT_ARD_MEAN_COL)
    status_metric = getattr(args, "status_metric", DEFAULT_STATUS_METRIC)
    delta_threshold = float(getattr(args, "delta_threshold", DEFAULT_DELTA_THRESHOLD))

    paths = output_paths(output_dir)

    if delta_threshold < 0:
        fail("--delta-threshold must be greater than or equal to 0.")

    try:
        network1_raw = read_results_table(network1_file)
        network2_raw = read_results_table(network2_file)

        network1_df = normalize_results_table(
            df=network1_raw,
            label="network1",
            gene_a_col=gene_a_col,
            gene_b_col=gene_b_col,
            cpr_mean_col=cpr_mean_col,
            rr_mean_col=rr_mean_col,
            ard_mean_col=ard_mean_col,
        )

        network2_df = normalize_results_table(
            df=network2_raw,
            label="network2",
            gene_a_col=gene_a_col,
            gene_b_col=gene_b_col,
            cpr_mean_col=cpr_mean_col,
            rr_mean_col=rr_mean_col,
            ard_mean_col=ard_mean_col,
        )

        comparison_df = compare_network_tables(
            network1_df=network1_df,
            network2_df=network2_df,
            status_metric=status_metric,
            delta_threshold=delta_threshold,
        )

        summary_df = build_summary_table(comparison_df)

        network1_only_df = comparison_df.loc[
            comparison_df["pair_presence"] == "network1_only"
        ].copy()

        network2_only_df = comparison_df.loc[
            comparison_df["pair_presence"] == "network2_only"
        ].copy()

        write_csv(comparison_df, paths["comparison"])
        write_csv(summary_df, paths["summary"])
        write_csv(network1_only_df, paths["network1_only"])
        write_csv(network2_only_df, paths["network2_only"])

        print("\nReGAIN network-analysis complete.")
        print(f"Comparison table:          {paths['comparison']}")
        print(f"Summary table:             {paths['summary']}")
        print(f"Network1-only pairs:       {paths['network1_only']}")
        print(f"Network2-only pairs:       {paths['network2_only']}")

        n_shared = int((comparison_df["pair_presence"] == "shared").sum())
        n_lost = int((comparison_df["status"] == "lost_from_network2").sum())
        n_gained = int((comparison_df["status"] == "gained_in_network2").sum())
        n_retained = int((comparison_df["status"] == "retained").sum())
        n_weakened = int((comparison_df["status"] == "weakened").sum())
        n_strengthened = int((comparison_df["status"] == "strengthened").sum())

        print(f"\nShared undirected pairs:   {n_shared}")
        print(f"Lost from network2:        {n_lost}")
        print(f"Gained in network2:        {n_gained}")
        print(f"Retained:                  {n_retained}")
        print(f"Weakened:                  {n_weakened}")
        print(f"Strengthened:              {n_strengthened}")

        print(
            f"\nStatus classification used '{status_metric}' mean deltas "
            f"with threshold {delta_threshold}."
        )
        print(
            "Pair matching is undirected, but directional metrics are preserved "
            "separately for Gene_1->Gene_2 and Gene_2->Gene_1."
        )

    except SystemExit:
        raise
    except Exception as exc:
        fail(str(exc))