#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 University of Utah

import os
import sys
import shlex
import argparse
import subprocess
import re
import difflib
from pathlib import Path
from tqdm import tqdm

MAGENTA = "\033[95m"
CYAN = "\033[96m"
GREEN = "\033[92m"
RESET = "\033[0m"

#Include ARMfinder citation, green ANSI escape code
print("\033[92m" + "\nReGAIN utilizes AMRfinderPlus:\nFeldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728. Published 2021 Jun 16. doi:10.1038/s41598-021-91456-0 \n" + "\033[0m")

NUC_EXTS = {".fa", ".fna", ".ffn", ".fasta", ".fas", ".frn"}
PROT_EXTS = {".faa", ".fa", ".fasta"}

def _which(prog: str) -> str | None:
    from shutil import which
    return which(prog)

ORGANISM_TOKENS = {
    "Acinetobacter_baumannii",
    "Bordetella_pertussis",
    "Burkholderia_cepacia",
    "Burkholderia_mallei",
    "Burkholderia_pseudomallei",
    "Campylobacter",
    "Citrobacter_freundii",
    "Clostridioides_difficile",
    "Corynebacterium_diphtheriae",
    "Enterobacter_asburiae",
    "Enterobacter_cloacae",
    "Enterococcus_faecalis",
    "Enterococcus_faecium",
    "Escherichia",
    "Haemophilus_influenzae",
    "Klebsiella_oxytoca",
    "Klebsiella_pneumoniae",
    "Neisseria_gonorrhoeae",
    "Neisseria_meningitidis",
    "Pseudomonas_aeruginosa",
    "Salmonella",
    "Serratia_marcescens",
    "Staphylococcus_aureus",
    "Staphylococcus_pseudintermedius",
    "Streptococcus_agalactiae",
    "Streptococcus_pneumoniae",
    "Streptococcus_pyogenes",
    "Vibrio_cholerae",
    "Vibrio_parahaemolyticus",
    "Vibrio_vulnificus",
}

ORGANISM_SYNONYMS = {
    "e_coli": "Escherichia",
    "escherichia_coli": "Escherichia",
    "c_difficile": "Clostridioides_difficile",
    "s_aureus": "Staphylococcus_aureus",
    "staph_aureus": "Staphylococcus_aureus",
    "group_b_strep": "Streptococcus_agalactiae",
    "strep_pneumoniae": "Streptococcus_pneumoniae",
}

_LOWER_MAP = {t.lower(): t for t in ORGANISM_TOKENS}  # for canonicalization

def _canon_organism(user_text: str) -> tuple[str | None, list[str]]:
    """
    Return (canonical_token_or_None, suggestions).
    Normalizes whitespace, punctuation, and case for robust matching.
    """
    raw = user_text.strip()
    # normalize: punctuation/whitespace -> underscore, collapse repeats
    norm = re.sub(r"[^\w]+", "_", raw).strip("_").lower()
    # exact (case-insensitive) match to known tokens
    if norm in _LOWER_MAP:
        return _LOWER_MAP[norm], []
    # synonyms
    if norm in ORGANISM_SYNONYMS:
        return ORGANISM_SYNONYMS[norm], []
    # suggestions (closest 3)
    sug = difflib.get_close_matches(norm, list(_LOWER_MAP.keys()), n=3, cutoff=0.6)
    sug = [ _LOWER_MAP[s] for s in sug ]
    return None, sug
    
def _is_nucleotide_file(p: Path) -> bool:
    return p.suffix.lower() in NUC_EXTS

def _is_protein_file(p: Path) -> bool:
    return p.suffix.lower() in PROT_EXTS

def build_amrfinder_cmd(
    infile: Path,
    mode: str,
    gff: Path | None,
    output_file: Path,
    organism: str | None,
    threads: int | None,
    database: Path | None,
    plus: bool,
    name: str | None,
    quiet: bool,
    ident_min: float | None,
    mutation_all: Path | None,
    nucleotide_output: Path | None,
    nucleotide_flank5_output: Path | None,
    nucleotide_flank5_size: int | None,
    protein_output: Path | None,
    print_node: bool,
) -> list[str]:
    
    cmd = ["amrfinder"]

    if mode == "nucleotide":
        cmd += ["-n", str(infile)]
    elif mode == "protein":
        cmd += ["-p", str(infile)]
    elif mode == "combined":
        #Require GFF file in combined mode
        if gff is None:
            raise ValueError("Combined mode requires --gff to be provided")
        #In combined mode, AMRFinder expects -p proteins + -g gff, and optionally -n genomic
        #Here we assume infile is nucleotide or protein—choose based on extension
        if _is_protein_file(infile):
            cmd += ["-p", str(infile), "-g", str(gff)]
        elif _is_nucleotide_file(infile):
            #If user feeds nucleotide file, we still need proteins + gff to use combined mode properly
            #Raise a helpful error to avoid misleading runs
            raise ValueError("Combined mode expects a protein FASTA for -p plus a matching GFF. Provide a .faa file to use `--mode combined`.")
        else:
            raise ValueError("Unable to infer input for combined mode; please supply a .faa file.")
    else:
        raise ValueError(f"Unknown mode: {mode}")
    
    if plus:
        cmd += ["--plus"]
    if organism:
        cmd += ["-O", organism]
    if threads is not None:
        cmd += ["--threads", str(threads)]
    if database is not None:
        cmd += ["-d", str(database)]
    if name:
        cmd += ["--name", name]
    if quiet:
        cmd += ["-q"]
    if print_node:
        cmd += ["--print_node"]
    if ident_min is not None:
        # Pass through exactly as AMRFinder expects (0-1 float). See wiki.
        cmd += ["--ident_min", str(ident_min)]
    if mutation_all is not None:
        cmd += ["--mutation_all", str(mutation_all)]
    if nucleotide_output is not None:
        cmd += ["--nucleotide_output", str(nucleotide_output)]
    if nucleotide_flank5_output is not None:
        cmd += ["--nucleotide_flank5_output", str(nucleotide_flank5_output)]
    if nucleotide_flank5_size is not None:
        cmd += ["--nucleotide_flank5_size", str(nucleotide_flank5_size)]
    if protein_output is not None:
        cmd += ["--protein_output", str(protein_output)]

    cmd += ["-o", str(output_file)]
    return cmd

def run(args):
    if getattr(args, "organism_list", False):
        print("\nSupported --organism tokens:\n  " + ", ".join(sorted(ORGANISM_TOKENS)))
        return
    
    if not args.fasta_directory:
        print("\033[91mERROR:\033[0m -f/--fasta-directory is required unless --organism-list is used")
        sys.exit(2)

    if args.mutation_all and not args.organism:
        print("\033[91mERROR:\033[0m point mutation audit only available with -O flag")
        sys.exit(2)

    exe = _which("amrfinder")
    if not exe:
        print("\033[91mERROR:\033[0m 'amrfinder' not found on PATH\nn")
        sys.exit(1)

    print(CYAN + f"Using AMRFinder: {exe}" + RESET)

    directory = Path(args.fasta_directory).expanduser().resolve()
    if not directory.is_dir():
        print(f"\033[91mERROR:\033[0m Input direcotry not found: {directory}\n")
        sys.exit(1)

    organism = args.organism
    if organism:
        canonical, suggestions = _canon_organism(organism)
        if canonical:
            organism = canonical
            print(CYAN + f"Organism token: {organism}" + RESET)
        else:
            msg = f"Unrecognized organism '{args.organism}'."
            if suggestions:
                msg += " Did you mean: " + ", ".join(suggestions) + "?"
            print(CYAN + msg + " Proceeding without -O." + RESET)
            organism = None

    #If user supplied an organism, auto-enable --plus unless they explicitly disabled it
    if organism and args.no_plus:
        print("\033[93mWarning:\033[0m You provided -O/--organism and --no-plus. "
            "Proceeding without --plus as requested.")
    effective_plus = (organism is not None and not args.no_plus) or (organism is None and not args.no_plus)
                
    if args.threads:
        print(GREEN + f"Using {args.threads} cores\n" + RESET)
    else:
        print(MAGENTA + "Using AMRFinder default threads\n" + RESET)

    outdir = Path(args.output_dir).expanduser().resolve() if args.output_dir else directory / "AMRfinder_Results"
    outdir.mkdir(parents=True, exist_ok=True)
    
    all_files = [p for p in directory.iterdir() if p.is_file()]
    if args.mode == "protein":
        in_files = [p for p in all_files if _is_protein_file(p)]
    elif args.mode == "nucleotide":
        in_files = [p for p in all_files if _is_nucleotide_file(p)]
    elif args.mode == "combined":
        in_files = [p for p in all_files if _is_protein_file(p)]
        if args.gff is None:
            print(f"\033[91mERROR:\033[0m Combinned mode requires --gff <annotation.gff>\n")
            sys.exit(1)
    else:
        print("\033[91mERROR:\033[0m Unknown mode specified\n")
        sys.exit(1)

    if not in_files:
        print("\033[91mERROR:\033[0m No input FASTA files found matching the selected mode")
        return
    
    if args.mode == "combined":
        gff_path = Path(args.gff).expanduser().resolve()
        if not gff_path.exists():
            print(f"\033[91mERROR:\033[0m GFF file not found: {gff_path}")
            sys.exit(1)
    else:
        gff_path = None

    failures = []

    for infile in tqdm(in_files, desc="Processing files", unit="file"):
        stem = infile.stem
        outfile = outdir / f"{stem}.amrfinder.csv"

        nuc_out = (outdir / f"{stem}.regions.fna") if args.nucleotide_output else None
        nuc_flank_out = (outdir / f"{stem}.regions_flank5.fna") if args.nucleotide_flank5_output else None
        prot_out = (outdir / f"{stem}.proteins.faa") if args.protein_output else None
        mut_all = (outdir / f"{stem}.mutation_all.tsv") if args.mutation_all else None

        try:
            cmd = build_amrfinder_cmd(
                infile=infile,
                mode=args.mode,
                gff=gff_path,
                output_file=outfile,
                organism=organism,
                threads=args.threads,
                database=Path(args.database).expanduser().resolve() if args.database else None,
                plus=effective_plus,
                name=args.name,
                quiet=args.quiet,
                ident_min=args.ident_min,
                mutation_all=mut_all,
                nucleotide_output=nuc_out,
                nucleotide_flank5_output=nuc_flank_out,
                nucleotide_flank5_size=args.nucleotide_flank5_size,
                protein_output=prot_out,
                print_node=args.print_node,
            )

            proc = subprocess.run(cmd, text=True, capture_output=True)
            if proc.returncode != 0:
                failures.append((infile.name, proc.stderr.strip()))
        except Exception as e:
            failures.append((infile.name, str(e)))

    if failures:
        print("\033[91mSome runs failed:\033[0m\n")
        for fname, err in failures:
            print(f" - {fname}: {err[:300]}") #truncate long BLAST/HMMER log
    else:
        print(GREEN + "All AMRFinderPlus runs completed\n" + RESET)