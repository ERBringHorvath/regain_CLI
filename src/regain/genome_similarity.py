#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 University of Utah

"""
genome_similarity.py

ReGAIN utility for genome similarity / clonal redundancy diagnostics

This module estimates pairwise genome similarity using either FastANI or Mash,
clusters genomes at user-defined or default thresholds, and writes user-facing 
reports summarizing potential clonal redundancy

Important:
    This module does NOT remove genomes, dereplicate matrices, or perform
    phylogenetic correction. It reports genome similarity structure so users can
    make informed decisions about downstream ReGAIN network interpretation

Expected CLI usage through cli.py:

    regain genome-similarity \
        -f genomes/ \
        -o genome_similarity_out \
        --method fastani \
        -T 8

    regain genome-similarity \
        -f genomes/ \
        -o genome_similarity_out \
        --method mash \
        --threshold 0.001 0.005 0.01 \
        -T 8

Default thresholds:
    fastani: 99.9 99.5 99.0
        Interpreted as ANI percent. Genomes cluster if ANI >= threshold

    mash: 0.001 0.005 0.01
        Interpreted as Mash distance. Genomes cluster if distance <= threshold

Outputs:
    genome_similarity.run_summary.csv
    genome_similarity.pairwise.csv
    genome_similarity.cluster_summary.csv
    genome_similarity.clusters_<threshold>.csv
    raw tool output files
    genome list files used for execution
"""

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import argparse
import csv
import shutil
import subprocess
import sys
import statistics

import pandas as pd

# Helpers
DEFAULT_FASTA_EXTENSIONS = [".fna", ".fa", ".fasta", ".ffn", ".fas"]
DEFAULT_FASTANI_THRESHOLDS = [99.9, 99.5, 99.0]
DEFAULT_MASH_THRESHOLDS = [0.001, 0.005, 0.01]

def fail(message: str, exit_code: int = 1) -> None:
    """
    Print a clear error message and exit
    """
    print(f"\nERROR: {message}\n", file=sys.stderr)
    raise SystemExit(exit_code)

def write_csv(df: pd.DataFrame, path: Path) -> None:
    """
    Write dataframe to CSV
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

def write_text_lines(lines: Iterable[str], path: Path) -> None:
    """
    Write one strng per line
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for line in lines:
            handle.write(f"{line}\n")

def check_external_tool(tool_name: str) -> str:
    """
    Confirm that an external tool exists in PATH
    """
    exe = shutil.which(tool_name)

    if exe is None:
        install_hint = {
            "fastANI": "mamba install -c conda-forge -c bioconda fastani",
            "mash": "mamba install -c conda-forge -c bioconda mash",
        }.get(tool_name, f"mamba install -c conda-forge -c bioconda {tool_name}")

        raise RuntimeError(
            f"Required tool '{tool_name}' was not found in PATH. "
            f"Suggested install command: {install_hint}"
        )

    return exe

def output_paths(output_dir: Path, method: str) -> Dict[str, Path]:
    """
    Build output paths
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "run_summary": output_dir / "genome_similarity.run_summary.csv",
        "pairwise": output_dir / "genome_similarity.pairwise.csv",
        "cluster_summary": output_dir / "genome_similarity.cluster_summary.csv",
        "selected_genomes": output_dir / "genome_similarity.selected_genomes.txt",
    }

    if method == "fastani":
        paths.update({
            "fastani_query_list": output_dir / "fastani_query_list.txt",
            "fastani_reference_list": output_dir / "fastani_reference_list.txt",
            "fastani_raw": output_dir / "fastani_raw.tsv",
        })
    elif method == "mash":
        paths.update({
            "mash_input_list": output_dir / "mash_input_genomes.txt",
            "mash_sketch": output_dir / "mash_genomes",
            "mash_sketch_file": output_dir / "mash_genomes.msh",
            "mash_raw": output_dir / "mash_raw.tsv",
        })

    return paths

def sanitize_threshold_for_filename(threshold: float) -> str:
    """
    Convert a threshold float into a filesystem-friendly string
    """
    return(str(threshold).replace(".", "p"))

def basename_without_fasta_suffix(path: Path, extensions: List[str]) -> str:
    """
    Generate a genome ID by stripping a recognized FASTA extension

    This handles names like:
        genome.fna -> genome
        genome.fasta -> genome
        ...
    """
    name = path.name

    for ext in sorted(extensions, key=len, reverse=True):
        if name.lower().endswith(ext.lower()):
            return name[: -len(ext)]

    return path.stem

#----- Genome selection -----
def normalize_extensions(extensions: List[str]) -> List[str]:
    """
    Ensure all extensions begin with a period and are lowercase
    """
    normalized = []

    for ext in extensions:
        ext = ext.strip()
        if not ext:
            continue
        if not ext.startswith("."):
            ext = f".{ext}"
        normalized.append(ext.lower())

    if not normalized:
        raise ValueError("At least one FASTA extension must be provided.")

    return normalized
    
def find_genome_files(
    fasta_dir: Path,
    extensions: List[str],
    recursive: bool = False,
) -> List[Path]:
    """
    Find FASTA files in the input directory
    """
    if not fasta_dir.exists():
        raise ValueError(f"FASTA directory not found: {fasta_dir}")

    if not fasta_dir.is_dir():
        raise ValueError(f"--fasta-dir must point to a directory: {fasta_dir}")

    pattern_iter = fasta_dir.rglob("*") if recursive else fasta_dir.glob("*")

    genome_files = []
    ext_set = set(e.lower() for e in extensions)

    for path in pattern_iter:
        if not path.is_file():
            continue

        if path.suffix.lower() in ext_set:
            genome_files.append(path.resolve())

    genome_files = sorted(genome_files, key=lambda p: str(p))

    if len(genome_files) < 2:
        raise ValueError(
            f"Found fewer than two genome FASTA files in {fasta_dir}. "
            f"Allowed extensions: {', '.join(extensions)}"
        )

    return genome_files

def read_genome_list(genome_list_file: Path) -> List[str]:
    """
    Read optional genome list

    Format:
        One filename or relativve path per row
        Blank lines are ignored
        Lines beginning with # are ignored
    """
    if not genome_list_file.exists():
        raise ValueError(f"Genome list file not found: {genome_list_file}")

    entries = []

    with genome_list_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            if value.startswith("#"):
                continue
            entries.append(value)

    if len(entries) < 2:
        raise ValueError(
            f"Genome list file must contain at least two genome entries: {genome_list_file}"
        )

    return entries

def resolve_genome_selection(
    fasta_dir: Path,
    genome_list_file: Optional[Path],
    extensions: List[str],
    recursive: bool = False,
) -> List[Path]:
    """
    Resolve the genome files to analyze

    If genome_list_file is not provided:
        Use all FASTA files discovered in fasta_dir

    if genome_list_file is provided:
        Treat entries as filenames or relative paths under fasta_dir
    """
    if genome_list_file is None:
        genome_files = find_genome_files(
            fasta_dir=fasta_dir,
            extensions=extensions,
            recursive=recursive,
        )
    else:
        entries = read_genome_list(genome_list_file)
        genome_files = []

        missing = []
        not_files = []

        for entry in entries:
            candidate = (fasta_dir / entry).resolve()

            if not candidate.exists():
                missing.append(entry)
                continue

            if not candidate.is_file():
                not_files.append(entry)
                continue

            genome_files.append(candidate)

        if missing:
            preview = ", ".join(missing[:10])
            extra = "" if len(missing) <= 10 else f" ... plus {len(missing) - 10} more"
            raise ValueError(
                f"The following genome-list entries were not found under {fasta_dir}: "
                f"{preview}{extra}"
            )

        if not_files:
            preview = ", ".join(not_files[:10])
            extra = "" if len(not_files) <= 10 else f" ... plus {len(not_files) - 10} more"
            raise ValueError(
                f"The following genome-list entries are not files: {preview}{extra}"
            )

        if len(genome_files) < 2:
            raise ValueError("At least two genome FASTA files are required.")

    # Remove duplicate exact paths while preserving order.
    seen = set()
    unique_files = []
    duplicates = []

    for path in genome_files:
        if path in seen:
            duplicates.append(path)
            continue
        seen.add(path)
        unique_files.append(path)

    if duplicates:
        dup_preview = ", ".join(str(p) for p in duplicates[:10])
        extra = "" if len(duplicates) <= 10 else f" ... plus {len(duplicates) - 10} more"
        raise ValueError(f"Duplicate genome paths detected: {dup_preview}{extra}")

    return unique_files

def build_genome_id_map(
    genome_files: List[Path],
    extensions: List[str],
) -> Dict[Path, str]:
    """
    Create genome IDs from FASTA filenames and de-duplicate if necessary
    """
    genome_id_map = {
        path: basename_without_fasta_suffix(path, extensions)
        for path in genome_files
    }

    id_to_paths: Dict[str, List[Path]] = {}

    for path, genome_id in genome_id_map.items():
        id_to_paths.setdefault(genome_id, []).append(path)

    duplicated_ids = {
        genome_id: paths
        for genome_id, paths in id_to_paths.items()
        if len(paths) > 1
    }

    if duplicated_ids:
        lines = []
        for genome_id, paths in list(duplicated_ids.items())[:10]:
            path_str = "; ".join(str(p) for p in paths)
            lines.append(f"{genome_id}: {path_str}")

        raise ValueError(
            "Duplicate genome IDs were generated from FASTA filenames. "
            "Please rename files or avoid recursive inputs that contain duplicate stems\n"
            + "\n".join(lines)
        )

    return genome_id_map

#------------------------
# External tool execution
#------------------------

def run_subprocess(command: List[str], stdout_path: Optional[Path] = None) -> None:
    """
    Run a subprocess command, optionally redirecting stdout to a file
    """
    try:
        if stdout_path is None:
            subprocess.run(command, check=True)
        else:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with stdout_path.open("w", encoding="utf-8") as out_handle:
                subprocess.run(command, check=True, stdout=out_handle)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"External command failed with exit code {exc.returncode}:\n"
            f"{' '.join(command)}"
        ) from exc
    
def run_fastani(
    genome_files: List[Path],
    paths: Dict[str, Path],
    threads: int,
) -> Path:
    """
    Run FastANI all-vs-all using the same genome list for query and reference
    """
    fastani_exe = check_external_tool("fastANI")

    genome_path_lines = [str(p) for p in genome_files]

    write_text_lines(genome_path_lines, paths["fastani_query_list"])
    write_text_lines(genome_path_lines, paths["fastani_reference_list"])

    command = [
        fastani_exe,
        "--ql", str(paths["fastani_query_list"]),
        "--rl", str(paths["fastani_reference_list"]),
        "-o", str(paths["fastani_raw"]),
        "-t", str(threads),
    ]

    run_subprocess(command)

    return paths["fastani_raw"]

def run_mash(
    genome_files: List[Path],
    paths: Dict[str, Path],
    threads: int,
    sketch_size: int = 10000,
    kmer_size: int = 21,
) -> Path:
    """
    Run Mash all-vs-all by sketching selected genomes and comparing the sketch
    file against itself
    """
    mash_exe = check_external_tool("mash")

    genome_path_lines = [str(p) for p in genome_files]
    write_text_lines(genome_path_lines, paths["mash_input_list"])

    # Mash takes file paths directly for sketching.
    # For large command lines, this could eventually be chunked or adapted,
    # but direct paths keep v1.8.0 simple.
    sketch_command = [
        mash_exe,
        "sketch",
        "-o", str(paths["mash_sketch"]),
        "-s", str(sketch_size),
        "-k", str(kmer_size),
        "-p", str(threads),
    ] + genome_path_lines

    run_subprocess(sketch_command)

    dist_command = [
        mash_exe,
        "dist",
        "-p", str(threads),
        str(paths["mash_sketch_file"]),
        str(paths["mash_sketch_file"]),
    ]

    run_subprocess(dist_command, stdout_path=paths["mash_raw"])

    return paths["mash_raw"]

#-------------------------
# Parsing pairwise results
#-------------------------

def path_to_genome_id(
    raw_path: str,
    genome_id_by_resolved_path: Dict[str, str],
    genome_id_by_name: Dict[str, str],
) -> str:
    """
    Resolve a raw path/name reported by an external tool back to a genome ID
    """
    resolved = str(Path(raw_path).resolve())

    if resolved in genome_id_by_resolved_path:
        return genome_id_by_resolved_path[resolved]

    name = Path(raw_path).name
    if name in genome_id_by_name:
        return genome_id_by_name[name]

    stem = Path(raw_path).stem
    if stem in genome_id_by_name:
        return genome_id_by_name[stem]

    # Fallback: use stem. This should rarely happen if genome maps are correct
    return stem

def parse_fastani_output(
    raw_file: Path,
    genome_id_map: Dict[Path, str],
) -> pd.DataFrame:
    """
    Parse FastANI output

    Expected columns:
        query reference ANI fragments_mapped total_fragments
    """
    genome_id_by_resolved_path = {
        str(path.resolve()): genome_id
        for path, genome_id in genome_id_map.items()
    }
    genome_id_by_name = {
        path.name: genome_id
        for path, genome_id in genome_id_map.items()
    }
    genome_id_by_name.update({
        path.stem: genome_id
        for path, genome_id in genome_id_map.items()
    })

    records = []

    if not raw_file.exists():
        raise ValueError(f"FastANI raw output not found: {raw_file}")

    with raw_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                parts = line.split()

            if len(parts) < 5:
                continue

            query_raw, reference_raw, ani_raw, fragments_mapped_raw, total_fragments_raw = parts[:5]

            genome_a = path_to_genome_id(
                query_raw,
                genome_id_by_resolved_path,
                genome_id_by_name,
            )
            genome_b = path_to_genome_id(
                reference_raw,
                genome_id_by_resolved_path,
                genome_id_by_name,
            )

            if genome_a == genome_b:
                continue

            ani = float(ani_raw)
            distance = 1.0 - (ani / 100.0)

            records.append({
                "genome_a": genome_a,
                "genome_b": genome_b,
                "similarity": ani,
                "distance": distance,
                "method": "fastani",
                "fragments_mapped": int(float(fragments_mapped_raw)),
                "total_fragments": int(float(total_fragments_raw)),
                "p_value": "",
                "matching_hashes": "",
                "status": "compared",
            })

    return pd.DataFrame(records)

def parse_mash_output(
    raw_file: Path,
    genome_id_map: Dict[Path, str],
) -> pd.DataFrame:
    """
    Parse Mash dist output

    Expected columns:
        reference query distance p_value matching_hashes
    """
    genome_id_by_resolved_path = {
        str(path.resolve()): genome_id
        for path, genome_id in genome_id_map.items()
    }
    genome_id_by_name = {
        path.name: genome_id
        for path, genome_id in genome_id_map.items()
    }
    genome_id_by_name.update({
        path.stem: genome_id
        for path, genome_id in genome_id_map.items()
    })

    records = []

    if not raw_file.exists():
        raise ValueError(f"Mash raw output not found: {raw_file}")

    with raw_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                parts = line.split()

            if len(parts) < 5:
                continue

            reference_raw, query_raw, distance_raw, p_value_raw, matching_hashes_raw = parts[:5]

            genome_a = path_to_genome_id(
                reference_raw,
                genome_id_by_resolved_path,
                genome_id_by_name,
            )
            genome_b = path_to_genome_id(
                query_raw,
                genome_id_by_resolved_path,
                genome_id_by_name,
            )

            if genome_a == genome_b:
                continue

            distance = float(distance_raw)
            similarity = 1.0 - distance

            records.append({
                "genome_a": genome_a,
                "genome_b": genome_b,
                "similarity": similarity,
                "distance": distance,
                "method": "mash",
                "fragments_mapped": "",
                "total_fragments": "",
                "p_value": p_value_raw,
                "matching_hashes": matching_hashes_raw,
                "status": "compared",
            })

    return pd.DataFrame(records)

def collapse_directional_pairs(pairwise_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert directional output into one undirected pair per genome pair

    If both A->B and B->A are present, numeric values are averaged
    """
    if pairwise_df.empty:
        return pairwise_df

    records = []

    working = pairwise_df.copy()
    working["pair_key"] = working.apply(
        lambda row: "||".join(sorted([str(row["genome_a"]), str(row["genome_b"])])),
        axis=1,
    )

    for pair_key, group in working.groupby("pair_key", sort=False):
        genomes = pair_key.split("||")
        genome_a, genome_b = genomes[0], genomes[1]

        method = group["method"].iloc[0]

        similarity_values = pd.to_numeric(group["similarity"], errors="coerce").dropna()
        distance_values = pd.to_numeric(group["distance"], errors="coerce").dropna()

        record = {
            "genome_a": genome_a,
            "genome_b": genome_b,
            "similarity": float(similarity_values.mean()) if not similarity_values.empty else "",
            "distance": float(distance_values.mean()) if not distance_values.empty else "",
            "method": method,
            "fragments_mapped": "",
            "total_fragments": "",
            "p_value": "",
            "matching_hashes": "",
            "status": "compared",
            "n_directional_results": group.shape[0],
        }

        if method == "fastani":
            fragments = pd.to_numeric(group["fragments_mapped"], errors="coerce").dropna()
            totals = pd.to_numeric(group["total_fragments"], errors="coerce").dropna()
            record["fragments_mapped"] = int(round(fragments.mean())) if not fragments.empty else ""
            record["total_fragments"] = int(round(totals.mean())) if not totals.empty else ""

        elif method == "mash":
            # If both directions are present, they should usually be identical.
            record["p_value"] = group["p_value"].iloc[0] if "p_value" in group else ""
            record["matching_hashes"] = group["matching_hashes"].iloc[0] if "matching_hashes" in group else ""

        records.append(record)

    collapsed = pd.DataFrame(records)

    if not collapsed.empty:
        collapsed = collapsed.sort_values(
            by=["genome_a", "genome_b"],
            ascending=[True, True],
        ).reset_index(drop=True)

    return collapsed

#-----------
# Clustering
#-----------

class UnionFind:
    """
    Small union-find function for connected-component clustering
    """

    def __init__(self, items: Iterable[str]):
        self.parent = {item: item for item in items}
        self.rank = {item: 0 for item in items}

    def find(self, item: str) -> str:
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, a: str, b: str) -> None:
        root_a = self.find(a)
        root_b = self.find(b)

        if root_a == root_b:
            return

        if self.rank[root_a] < self.rank[root_b]:
            self.parent[root_a] = root_b
        elif self.rank[root_a] > self.rank[root_b]:
            self.parent[root_b] = root_a
        else:
            self.parent[root_b] = root_a
            self.rank[root_a] += 1

def pair_passes_threshold(row: pd.Series, method: str, threshold: float) -> bool:
    """
    Decide whether a pair should be linked at a threshold
    """
    if method == "fastani":
        try:
            return float(row["similarity"]) >= threshold
        except Exception:
            return False

    if method == "mash":
        try:
            return float(row["distance"]) <= threshold
        except Exception:
            return False

    raise ValueError(f"Unsupported method for thresholding: {method}")

def cluster_pairs_by_threshold(
    pairwise_df: pd.DataFrame,
    genome_ids: List[str],
    method: str,
    threshold: float,
) -> pd.DataFrame:
    """
    Cluster genomes using connected components at similarity/distance cutoff
    """
    uf = UnionFind(genome_ids)

    for _, row in pairwise_df.iterrows():
        if pair_passes_threshold(row, method, threshold):
            uf.union(str(row["genome_a"]), str(row["genome_b"]))

    component_members: Dict[str, List[str]] = {}

    for genome in genome_ids:
        root = uf.find(genome)
        component_members.setdefault(root, []).append(genome)

    # Sort clusters by size descending, then first genome ID.
    clusters = sorted(
        component_members.values(),
        key=lambda members: (-len(members), sorted(members)[0]),
    )

    records = []

    for cluster_index, members in enumerate(clusters, start=1):
        cluster_id = f"cluster_{cluster_index:05d}"
        cluster_size = len(members)

        for genome in sorted(members):
            records.append({
                "threshold": threshold,
                "method": method,
                "genome": genome,
                "cluster_id": cluster_id,
                "cluster_size": cluster_size,
            })

    return pd.DataFrame(records)

def summarize_cluster_file(cluster_df: pd.DataFrame, method: str, threshold: float) -> Dict[str, object]:
    """
    Summarize one cluster assignment table
    """
    n_genomes = cluster_df.shape[0]

    cluster_sizes = (
        cluster_df[["cluster_id", "cluster_size"]]
        .drop_duplicates()
        ["cluster_size"]
        .astype(int)
        .tolist()
    )

    if not cluster_sizes:
        return {
            "method": method,
            "threshold": threshold,
            "n_genomes": n_genomes,
            "n_clusters": 0,
            "n_singleton_clusters": 0,
            "largest_cluster_size": 0,
            "percent_genomes_in_largest_cluster": 0.0,
            "percent_genomes_in_top_5_clusters": 0.0,
            "mean_cluster_size": 0.0,
            "median_cluster_size": 0.0,
        }

    cluster_sizes_sorted = sorted(cluster_sizes, reverse=True)
    n_clusters = len(cluster_sizes_sorted)
    n_singletons = sum(1 for size in cluster_sizes_sorted if size == 1)
    largest = cluster_sizes_sorted[0]
    top_5_sum = sum(cluster_sizes_sorted[:5])

    return {
        "method": method,
        "threshold": threshold,
        "n_genomes": n_genomes,
        "n_clusters": n_clusters,
        "n_singleton_clusters": n_singletons,
        "largest_cluster_size": largest,
        "percent_genomes_in_largest_cluster": round(100.0 * largest / n_genomes, 4) if n_genomes else 0.0,
        "percent_genomes_in_top_5_clusters": round(100.0 * top_5_sum / n_genomes, 4) if n_genomes else 0.0,
        "mean_cluster_size": round(float(statistics.mean(cluster_sizes_sorted)), 4),
        "median_cluster_size": round(float(statistics.median(cluster_sizes_sorted)), 4),
    }

#------------
# Run summary
#------------

def build_run_summary(
    method: str,
    fasta_dir: Path,
    genome_files: List[Path],
    genome_list_file: Optional[Path],
    thresholds: List[float],
    threads: int,
    recursive: bool,
    extensions: List[str],
    pairwise_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build a high-level run summary
    """
    possible_pairs = len(genome_files) * (len(genome_files) - 1) // 2
    observed_pairs = pairwise_df.shape[0]

    rows = [
        {"metric": "method", "value": method},
        {"metric": "fasta_dir", "value": str(fasta_dir)},
        {"metric": "genome_list_file", "value": str(genome_list_file) if genome_list_file else ""},
        {"metric": "n_genomes", "value": len(genome_files)},
        {"metric": "n_possible_undirected_pairs", "value": possible_pairs},
        {"metric": "n_observed_undirected_pairs", "value": observed_pairs},
        {"metric": "n_missing_undirected_pairs", "value": possible_pairs - observed_pairs},
        {"metric": "thresholds", "value": "; ".join(str(t) for t in thresholds)},
        {"metric": "threads", "value": threads},
        {"metric": "recursive", "value": recursive},
        {"metric": "extensions", "value": "; ".join(extensions)},
    ]

    if method == "fastani":
        rows.append({
            "metric": "threshold_interpretation",
            "value": "Genomes are linked if ANI percent is greater than or equal to threshold.",
        })
    elif method == "mash":
        rows.append({
            "metric": "threshold_interpretation",
            "value": "Genomes are linked if Mash distance is less than or equal to threshold. Similarity is reported as 1 - Mash distance and is not ANI.",
        })

    return pd.DataFrame(rows)

#-----------------
# Main entry point
#-----------------

def run(args: argparse.Namespace) -> None:
    """
    Entry point used by cli.py
    """
    fasta_dir = Path(args.fasta_dir)
    output_dir = Path(args.output_dir)
    method = getattr(args, "method", "fastani").lower()
    genome_list = getattr(args, "genome_list", None)
    genome_list_file = Path(genome_list) if genome_list else None
    recursive = bool(getattr(args, "recursive", False))
    threads = int(getattr(args, "threads", 2))
    extensions = normalize_extensions(getattr(args, "extensions", DEFAULT_FASTA_EXTENSIONS))

    thresholds = getattr(args, "threshold", None)

    if thresholds is None:
        if method == "fastani":
            thresholds = DEFAULT_FASTANI_THRESHOLDS
        elif method == "mash":
            thresholds = DEFAULT_MASH_THRESHOLDS
        else:
            fail(f"Unsupported method: {method}")
    else:
        thresholds = [float(t) for t in thresholds]

    if threads < 1:
        fail("--threads must be an integer greater than or equal to 1.")

    if method not in {"fastani", "mash"}:
        fail("--method must be one of: fastani, mash")

    paths = output_paths(output_dir, method)

    try:
        genome_files = resolve_genome_selection(
            fasta_dir=fasta_dir,
            genome_list_file=genome_list_file,
            extensions=extensions,
            recursive=recursive,
        )

        genome_id_map = build_genome_id_map(
            genome_files=genome_files,
            extensions=extensions,
        )

        genome_ids = [genome_id_map[path] for path in genome_files]

        selected_records = [
            {
                "genome_id": genome_id_map[path],
                "path": str(path),
            }
            for path in genome_files
        ]
        write_csv(pd.DataFrame(selected_records), paths["selected_genomes"])

        if method == "fastani":
            raw_output = run_fastani(
                genome_files=genome_files,
                paths=paths,
                threads=threads,
            )
            raw_pairwise_df = parse_fastani_output(
                raw_file=raw_output,
                genome_id_map=genome_id_map,
            )

        elif method == "mash":
            raw_output = run_mash(
                genome_files=genome_files,
                paths=paths,
                threads=threads,
                sketch_size=int(getattr(args, "sketch_size", 10000)),
                kmer_size=int(getattr(args, "kmer_size", 21)),
            )
            raw_pairwise_df = parse_mash_output(
                raw_file=raw_output,
                genome_id_map=genome_id_map,
            )

        else:
            raise ValueError(f"Unsupported method: {method}")

        pairwise_df = collapse_directional_pairs(raw_pairwise_df)
        write_csv(pairwise_df, paths["pairwise"])

        cluster_summary_records = []

        for threshold in thresholds:
            cluster_df = cluster_pairs_by_threshold(
                pairwise_df=pairwise_df,
                genome_ids=genome_ids,
                method=method,
                threshold=threshold,
            )

            threshold_label = sanitize_threshold_for_filename(threshold)
            cluster_path = output_dir / f"genome_similarity.clusters_{threshold_label}.csv"
            write_csv(cluster_df, cluster_path)

            cluster_summary_records.append(
                summarize_cluster_file(
                    cluster_df=cluster_df,
                    method=method,
                    threshold=threshold,
                )
            )

        cluster_summary_df = pd.DataFrame(cluster_summary_records)
        write_csv(cluster_summary_df, paths["cluster_summary"])

        run_summary_df = build_run_summary(
            method=method,
            fasta_dir=fasta_dir,
            genome_files=genome_files,
            genome_list_file=genome_list_file,
            thresholds=thresholds,
            threads=threads,
            recursive=recursive,
            extensions=extensions,
            pairwise_df=pairwise_df,
        )
        write_csv(run_summary_df, paths["run_summary"])

        print("\nReGAIN genome-similarity complete.")
        print(f"Method:                   {method}")
        print(f"Genomes analyzed:          {len(genome_files)}")
        print(f"Pairwise table:            {paths['pairwise']}")
        print(f"Cluster summary:           {paths['cluster_summary']}")
        print(f"Run summary:               {paths['run_summary']}")
        print(f"Selected genomes:          {paths['selected_genomes']}")

        if method == "fastani":
            print(f"FastANI raw output:        {paths['fastani_raw']}")
        elif method == "mash":
            print(f"Mash raw output:           {paths['mash_raw']}")
            print(
                "Note: Mash similarity is reported as 1 - Mash distance; "
                "it should not be interpreted as ANI."
            )

        possible_pairs = len(genome_files) * (len(genome_files) - 1) // 2
        observed_pairs = pairwise_df.shape[0]

        if observed_pairs < possible_pairs:
            print(
                f"\nWarning: observed {observed_pairs} of {possible_pairs} possible "
                "undirected genome pairs. Missing pairs may reflect divergent genomes, "
                "tool-specific filtering, or failed comparisons."
            )

        print(
            "\nNote: genome-similarity reports potential clonal/genome relatedness "
            "structure only. It does not remove genomes or correct downstream "
            "Bayesian network inference."
        )

    except SystemExit:
        raise
    except Exception as exc:
        fail(str(exc))


    