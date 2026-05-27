#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 University of Utah

import argparse, logging, sys, subprocess, os, shutil, importlib
import contextlib, io, importlib, traceback
from importlib import import_module
from datetime import datetime
from . import __version__

HERE = os.path.dirname(os.path.realpath(__file__))

def _bash():
    b = shutil.which("bash")
    if not b:
        sys.exit("ERROR: bash not found in PATH")
    return b

def _rscript():
    r = shutil.which("Rscript")
    if not r:
        sys.exit("ERROR: Rscript not found in PATH")
    return r

def _path_in_pkg(*names):
    return os.path.join(HERE, *names)

# ---- Python module loader (safe) ----
available_modules = {}  # name -> callable

@contextlib.contextmanager
def _silence_stdio():
    _out, _err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = io.StringIO(), io.StringIO()
        yield
    finally:
        sys.stdout, sys.stderr = _out, _err

def safe_import(name, func_name):
    """
    Try to import .<name> and get <func_name>. On success, cache callable.
    Returns (callable_or_None, error_message_or_None).
    """
    pkg = __package__ or "regain"
    try:
        with _silence_stdio():
            mod = importlib.import_module(f".{name}", pkg)
        fn = getattr(mod, func_name)  # may raise AttributeError
        available_modules[name] = fn
        return fn, None
    except Exception as e:
        # short, helpful reason
        etype = type(e).__name__
        msg = f"{etype}: {e}"
        return None, msg

def set_py_handler(parser, module_name, func_name):
    """
    Attach a handler that imports and calls regain.<module_name>.<func_name> at runtime.
    """
    def _handler(args, _mod=module_name, _fn=func_name):
        try:
            mod = import_module(f"regain.{_mod}")
            fn = getattr(mod, _fn, None)
        except Exception as e:
            parser.error(
                f"Failed to import 'regain.{_mod}.{_fn}': {e}\n"
                "Run 'regain --module-health' for details."
            )
        if not callable(fn):
            parser.error(
                f"'regain.{_mod}.{_fn}' is not callable or was not found.\n"
                "Run 'regain --module-health' for details."
            )
        return fn(args)
    parser.set_defaults(func=_handler)

def run_bnL(args):
    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'bnCPQuery.sh')
    command = [
        path_to_shell_script,
        "--input", args.input,
        "--metadata", args.metadata,
        "--output-boot", args.output_boot,
        "--threads", str(args.threads) if args.threads is not None else "2",
        "--bootstraps", str(args.bootstraps),
        "--resamples", str(args.resamples),
        "--iss", str(args.iss if args.iss is not None else 10),
        "--cp-samples", str(args.cp_samples if args.cp_samples is not None else 10000),

    ]
    if getattr(args, "blacklist", None):
        command += ["--blacklist", args.blacklist]

    if getattr(args, "no_viz", False):
        command += ["--no-viz"]

    subprocess.run(command, check=True)

def run_bnS(args):
    dir_of_executable = os.path.dirname(os.path.realpath(__file__))
    path_to_shell_script = os.path.join(dir_of_executable, 'bnQueryGrain.sh')
    command = [
        path_to_shell_script,
        "--input", args.input,
        "--metadata", args.metadata,
        "--output-boot", args.output_boot,
        "--threads", str(args.threads) if args.threads is not None else "2",
        "--bootstraps", str(args.bootstraps),
        "--resamples", str(args.resamples),
        "--iss", str(args.iss if args.iss is not None else 10),
    ]
    if getattr(args, "blacklist", None):
        command += ["--blacklist", args.blacklist]
    
    if getattr(args, "no_viz", False):
        command += ["--no-viz"]

    subprocess.run(command, check=True)

def run_mva(args):
    import os, shutil, subprocess
    bash = shutil.which("bash") or "/bin/bash"
    script_dir = os.path.dirname(os.path.realpath(__file__))
    sh = os.path.join(script_dir, "mva.sh")

    cmd = [
        bash, sh,
        "--input", args.input,
        "--method", args.method,
        "--k", str(args.k if args.k is not None else 1),
        "--confidence", str(args.confidence),
        "--seed", str(getattr(args, "seed", 42)),
        "--label", getattr(args, "label", "auto"),
        "--point-size", str(getattr(args, "point_size", 3.5)),
        "--alpha", str(getattr(args, "alpha", 0.75)),
        "--pseudocount", str(getattr(args, "pseudocount", 1e-6)),
        "--png-out", getattr(args, "png_out", "MVA.png"),
        "--pdf-out", getattr(args, "pdf_out", "MVA.pdf"),
        "--coords-out", getattr(args, "coords_out", "MVA_coordinates.csv"),
        "--dist-out", getattr(args, "dist_out", "MVA_distance.csv"),
        "--pcoa-correction", getattr(args, "pcoa_correction", "auto"),
    ]
    if getattr(args, "no_ellipses", False):
        cmd.append("--no-ellipses")
    if getattr(args, "save_dist", False):
        cmd.append("--save-dist")

    subprocess.run(cmd, check=True)


def run_network(args):
    dir_of_executables = os.path.dirname(os.path.realpath(__file__))
    sh = os.path.join(dir_of_executables, 'bayesnetwork.sh')
    cmd = [
        sh,
        '--boot', args.network,
        '--input', args.input,
        '--metadata', args.metadata,
        '--stats', args.statistics_results,
        '--threshold', str(args.threshold),
        '--seed', str(args.seed),
        '--html-out', args.html_out,
        '--pdf-out', args.pdf_out,
        '--width-metric', args.width_metric,
        '--rr-threshold', str(args.rr_threshold),
    ]
    if getattr(args, 'blacklist', None):
        cmd += ['--blacklist', args.blacklist]
    subprocess.run(cmd, check=True)

def check_curate_range(value):
    ivalue = int(value)
    if ivalue < 0 or ivalue > 100:
        raise argparse.ArgumentTypeError(f" \033[91m{value} is out of allowed range (0-100)\033[0m \n")
    return ivalue

def check_extract_range(value):
    value = float(value)
    if value < 0 or value > 100:
        raise argparse.ArgumentTypeError(f" \033[91mValue is out of allowed range (0-100)\033[0m \n")
    return value

def _print_module_health():
    print("\n\033[96mReGAIN Module Status Report\033[0m")

    #Python modules to check: (module, func)
    py_targets = [
        ("AMR", "run"),
        ("matrix", "run"),
        ("curate", "run"),
        ("extract", "run"),
        ("combine", "run"),
        ("collapse_features", "run"),
        ("matrix_summary", "run"),
        ("genome_similarity", "run"),
        ("network_analysis", "run"),
    ]

    ok = []
    bad = []

    for mod, fn in py_targets:
        tryfn = safe_import(mod, fn)
        if tryfn is None:
            bad.append(f"{mod}.{fn}")
        else:
            ok.append(f"{mod}.{fn}")

    # Shell/R-backed commands: check presence + toolchain
    shell_targets = [
        ("bnS", "bnQueryGrain.sh", "bayesQueryGrain.r"),
        ("bnL", "bnCPQuery.sh",    "bayesCPQuery.r"),
        ("MVA", "mva.sh",          "MVA.r"),
        ("network", "bayesnetwork.sh", "bayesnetwork.r"),
    ]

    bash_ok = shutil.which("bash") is not None
    r_ok    = shutil.which("Rscript") is not None

    for label, sh, r in shell_targets:
        sh_path = _path_in_pkg(sh)
        r_path  = _path_in_pkg(r)
        parts = []
        parts.append("bash:OK" if bash_ok else "bash:MISSING")
        parts.append("Rscript:OK" if r_ok else "Rscript:MISSING")
        parts.append(f"{sh}:OK" if os.path.exists(sh_path) else f"{sh}:MISSING")
        parts.append(f"{r}:OK"  if os.path.exists(r_path)  else f"{r}:MISSING")
        status_line = ", ".join(parts)
        if all(s.endswith("OK") for s in parts):
            ok.append(label)
        else:
            bad.append(f"{label} ({status_line})")

    if ok:
        for name in ok:
            print(f" - {name}: \033[92mAvailable\033[0m")
    if bad:
        for name in bad:
            print(f" - {name}: \033[91mBroken or Missing\033[0m")

    print()  # blank line

def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog='regain',
        description='ReGAIN: Resistance Gene Association and Inference Network',
        epilog="For detailed instructions, visit: https://github.com/ERBringHorvath/regain_CLI"
    )
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('--module-health', action='store_true',
                        help='Check availability of submodules and external tooling, then exit.')
    parser.add_argument('--citation', action='store_true',
                        help="Print software and dependency citations and exit")

    subparsers = parser.add_subparsers(dest='command')

    #---------------- AMR (Python) ----------------
    amr_parser = subparsers.add_parser(
    'AMR',
    help='Module 1: Resistance/Virulence gene identification via AMRFinderPlus'
    )

    amr_parser.add_argument('-f', '--fasta-directory',
                            help='Directory containing FASTA files')
    amr_parser.add_argument('--mode', choices=['nucleotide', 'protein', 'combined'],
                            default='nucleotide',
                            help='Input mode: nucleotide, protein, or combined (protein + GFF annotation file)')
    amr_parser.add_argument('--gff', help='GFF for combined mode (required when --mode combined)')

    # AMRFinder+ pass-throughs
    amr_parser.add_argument('-O', '--organism',
                            help='Organism token for point-mutation screening and blacklisting')
    amr_parser.add_argument('-T', '--threads', type=int,
                            help='Threads for AMRFinderPlus (--threads)')
    amr_parser.add_argument('-D', '--database', metavar='DB_DIR',
                            help='Alternate AMRFinderPlus database directory (-d/--database)')
    amr_parser.add_argument('--no-plus', action='store_true',
                            help='Disable --plus (include virulence/stress genes by default)')
    amr_parser.add_argument('--name',
                            help='Prepend a sample/run identifier to AMRFinderPlus output (generates column "name")')
    amr_parser.add_argument('--quiet', action='store_true',
                            help='Suppress AMRFinderPlus status messages')
    amr_parser.add_argument('--ident-min', type=float,
                            help='Override minimum identity (0-1). Use only if you have a specific reason.')
    amr_parser.add_argument('--print-node', action='store_true',
                            help='Add hierarchy node column')
    amr_parser.add_argument('--mutation-all', action='store_true',
                            help='Emit a per-sample point-mutation audit file')
    amr_parser.add_argument('--nucleotide-output', action='store_true',
                            help='Write FASTA of detected nucleotide regions')
    amr_parser.add_argument('--nucleotide-flank5-output', action='store_true',
                            help='Write FASTA of detected nucleotide regions + 5\' flank')
    amr_parser.add_argument('--nucleotide-flank5-size', type=int,
                            help='Number of additional 5\' bases for flank output')
    amr_parser.add_argument('--protein-output', action='store_true',
                            help='Write FASTA of detected proteins from input protein FASTA')
    amr_parser.add_argument('-o', '--output-dir', help='Output directory path')
    amr_parser.add_argument('--organism-list', action='store_true',
                        help='Print the built-in list of supported -O tokens and exit')

    set_py_handler(amr_parser, 'AMR', 'run')

    #---------------- Matrix (Python) ----------------
    matrix_parser = subparsers.add_parser('matrix', help='Module 1.1: Dataset Creation...')

    matrix_parser.add_argument(
        '-d', '--directory', 
        required=True, 
        help='Input directory to CSV files to search'
    )

    matrix_parser.add_argument(
        '-o', '--output-dir',
        default='ReGAIN_Dataset',
        help='Optional output directory name. Default: ReGAIN_Dataset'
    )
    
    matrix_parser.add_argument(
        '--keep-gene-names', 
        action='store_true', 
        help='Replace special characters...'
    )

    matrix_parser.add_argument(
        '--min', 
        type=int, 
        default=5, 
        required=True, 
        help='Minimum required gene occurrences'
    )

    matrix_parser.add_argument(
        '--max', 
        type=int, 
        help='Maximum allowed gene occurrences'
    )

    matrix_parser.add_argument(
        '--gene-type', 
        required=True,
        choices=['resistance', 'virulence', 'all'],                       
        help='specify gene type: resistance, virulence, all'
    )

    matrix_parser.add_argument(
        '--report-all', 
        action='store_true', 
        help='Return actual gene count...'
    )

    matrix_parser.add_argument(
        '--expected-extension-string',
        default='amrfinder.csv',
        help='Expected filename ending for AMRFinder output files (default: .amrfinder.csv (regain AMR default out))'
    )

    matrix_parser.add_argument(
        '--subtype-col',
        default=None,
        help='Column name for AMRFinder subtype annotations. Default: "Element subtype" or "Subtype"'
    )

    matrix_parser.add_argument(
        '--gene-col',
        default=None,
        help='Column name for gene/element symbols. Default: "Element symbol" or "Gene symbol"'
    )

    matrix_parser.add_argument(
        '--class-col',
        default='Class',
        help='Column name for AMRFinder Class annotations. Default: Class'
    )

    matrix_parser.add_argument(
        '--subclass-col',
        default='Subclass',
        help='Column name for AMRFinder Subclass annotations. Default: Subclass'
    )

    matrix_parser.add_argument(
        '--verbose-gene-report',
        action='store_true',
        help=(
            'Write combined_AMR_results_unfiltered.csv as a verbose per-source report: '
            'Retain duplicate gene rows and append a source-file column. '
            'Default behavior drops duplicate gene rows to create a simple gene dictionary'
        )
    )

    matrix_parser.add_argument(
        '--force-overwrite',
        action='store_true',
        help='Delete an existing output directory before writing new results'
    )
    set_py_handler(matrix_parser, 'matrix', 'run')

    #---------------- bnS (shell/R) ----------------
    bnS_parser = subparsers.add_parser('bnS', help='Module 2: Bayesian Network (< 100 genes)')
    bnS_parser.add_argument('-i', '--input', required=True, help='Input data matrix in CSV format')
    bnS_parser.add_argument('-M', '--metadata', required=True, help='Input metadata file with genes to query')
    bnS_parser.add_argument('-o', '--output-boot', required=True, help='Output file name for Network')
    bnS_parser.add_argument('-T', '--threads', type=int, help='Dedicated threads')
    bnS_parser.add_argument('-n', '--bootstraps', required=True, type=int, help='Number of bootstraps to perform')
    bnS_parser.add_argument('-r', '--resamples', type=int, required=True, help='Number of resamples to perform')
    bnS_parser.add_argument('-b', '--blacklist', help='Optional blacklist CSV (no header): from,to')
    bnS_parser.add_argument('--iss', type=int, default=10, help='Imaginary sample size for BDe score (default: 10)')
    bnS_parser.add_argument('--no-viz', dest='no_viz', action='store_true', help='Skip HTML/PDF visualization')
    bnS_parser.set_defaults(func=run_bnS)

    #---------------- bnL (shell/R) ----------------
    bnL_parser = subparsers.add_parser('bnL', help='Module 2: Bayesian Network (≥ 100 genes)')
    bnL_parser.add_argument('-i', '--input', required=True, help='Input data matrix in CSV format')
    bnL_parser.add_argument('-M', '--metadata', required=True, help='Input metadata file with genes to query')
    bnL_parser.add_argument('-o', '--output-boot', required=True, help='Output file name for Network')
    bnL_parser.add_argument('-T', '--threads', type=int, help='Dedicated threads')
    bnL_parser.add_argument('-n', '--bootstraps', required=True, type=int, help='Number of bootstraps to perform')
    bnL_parser.add_argument('-r', '--resamples', required=True, type=int, help='Number of resamples to perform')
    bnL_parser.add_argument('-b', '--blacklist', help='Optional blacklist CSV (no header): from,to')
    bnL_parser.add_argument('--iss', type=int, default=10, help='Imaginary sample size for BDe score (default: 10)')
    bnL_parser.add_argument('--no-viz', dest='no_viz', action='store_true', help='Skip HTML/PDF visualization')
    bnL_parser.add_argument('--cp-samples', dest='cp_samples', type=int, default=10000,
                            help='Monte Carlo samples for cpquery (default: 10000)')
    bnL_parser.set_defaults(func=run_bnL)

    #---------------- Network (shell/R) ----------------
    network_parser = subparsers.add_parser('network', help='Visualize Bayesian network (HTML + PDF)')
    network_parser.add_argument('-N', '--network', required=True, help='Bootstrapped network .rds')
    network_parser.add_argument('-i', '--input', required=True, help='Input data matrix in CSV format')
    network_parser.add_argument('-M', '--metadata', required=True, help='Metadata CSV')
    network_parser.add_argument('-s', '--statistics-results', required=True, help='Stats CSV from bnS/bnL')
    network_parser.add_argument('--threshold', type=float, default=0.5, help='Averaged network threshold (default 0.5)')
    network_parser.add_argument('--seed', type=int, default=42, help='Layout seed (default 42)')
    network_parser.add_argument('--html-out', dest='html_out', default='Bayesian_Network.html', help='HTML output filename')
    network_parser.add_argument('--pdf-out', dest='pdf_out', default='Bayesian_Network.pdf', help='PDF output filename')
    network_parser.add_argument('-b', '--blacklist', help='Optional blacklist CSV (no header): from,to')
    network_parser.add_argument('--width-metric', choices=['auto','abs_mean','abs_ci','cp_ci','cp_mean'], default='auto',
                                help='Edge-width metric selection')
    network_parser.add_argument('--rr-threshold', type=float, default=1.0, help='RR color threshold (default 1.0)')
    network_parser.set_defaults(func=run_network)

    #---------------- MVA (shell/R) ----------------
    mva_parser = subparsers.add_parser('MVA', help='Multivariate analysis (PCoA + k-means + ellipses)')
    mva_parser.add_argument('-i','--input', required=True, help='Input data file in CSV format')
    mva_parser.add_argument('-m','--method', default='euclidean',
                            help='Distance: manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, chord, hellinger, aitchison, mahalanobis')
    mva_parser.add_argument('--k', type=int, default=1, help='k for k-means; 0 = auto (2..10)')
    mva_parser.add_argument('-C','--confidence', type=float, default=0.95, help='Ellipse confidence')
    mva_parser.add_argument('--seed', type=int, default=42, help="Set custom random seed (default = 42)")
    mva_parser.add_argument('--label', choices=['none','auto','all'], default='auto', help="Manage data point labels")
    mva_parser.add_argument('--point-size', type=float, default=3.5, help="Manage data point size")
    mva_parser.add_argument('--alpha', type=float, default=0.75, help="Manage data point opacity")
    mva_parser.add_argument('--no-ellipses', dest='no_ellipses', action='store_true', help="Do not show confidence ellipses")
    mva_parser.add_argument('--pseudocount', type=float, default=1e-6, help="Add a small pseudocount for negative eigenvalues")
    mva_parser.add_argument('--save-dist', action='store_true', help="Save distances to CSV")
    mva_parser.add_argument('--png-out', default='MVA.png', help="Manually save PNG (default = MVA.png)")
    mva_parser.add_argument('--pdf-out', default='MVA.pdf', help="Manually save PDF (default = MVA.pdf)")
    mva_parser.add_argument('--coords-out', default='MVA_coordinates.csv', help="Save coordinates to CSV (default = MVA_coordinates.csv)")
    mva_parser.add_argument('--dist-out', default='MVA_distance.csv', help="Manually save distance file (default = MVA_distance.csv)")
    mva_parser.add_argument('--pcoa-correction', choices=['auto','none','lingoes','cailliez'], default='auto', help="Apply PCoA correction (default = auto)")
    mva_parser.set_defaults(func=run_mva)

    #---------------- Curate (Python) ----------------
    curate_parser = subparsers.add_parser('curate', help='ReGAIN Curate: Create a curated dataset...')
    curate_parser.add_argument('-f', '--fasta-directory', required=True, help='Input directory to genomes in FASTA')
    curate_parser.add_argument('-q', '--query', required=True, help='Input directory to gene FASTA queries')
    curate_parser.add_argument('-T', '--threads', type=int, help='How many cores to dedicate')
    curate_parser.add_argument('--nucleotide-query', action='store_true', help='Use blastn for nucleotide queries')
    curate_parser.add_argument('--report-all', action='store_true', help='Report all BLAST hits...')
    curate_parser.add_argument('--evalue', type=float, help="Override e-value threshold (default = 1e-5)")
    curate_parser.add_argument('--perc', type=check_curate_range, help='Override percent identity threshold (default 90)')
    curate_parser.add_argument('--cov', type=check_curate_range, help='Override query coverage threshold (default 75)')
    curate_parser.add_argument('--min-seq-len', type=int, help='Minimum sequence length')
    curate_parser.add_argument('--keep-gene-names', action='store_true', help='Do not remove special characters')
    curate_parser.add_argument('--min', type=int, default=5, required=True, help='Minimum required gene occurrence')
    curate_parser.add_argument('--max', type=int, required=True, help='Maximum allowed gene occurrence')
    set_py_handler(curate_parser, 'curate', 'run')

    #---------------- Extract (Python) ----------------

    parser_extract = subparsers.add_parser('extract', help="ReGAIN Curate: extract sequences...")
    parser_extract.add_argument('-c', '--csv-path', required=True, help="Path to BLAST results files")
    parser_extract.add_argument('-f', '--fasta-directory', required=True, help="Path to reference FASTA assemblies")
    parser_extract.add_argument('-o', '--output-fasta', required=True, help="Output FASTA file name")
    parser_extract.add_argument('-T', '--threads', help="How many CPUs to dedicate")
    parser_extract.add_argument('--evalue', type=float, default=1e-5, help='Override e-value threshold (default = 1e-5)')
    parser_extract.add_argument('--perc', type=check_extract_range, default=90.0, help='Override percent identity threshold (default = 90.0)')
    parser_extract.add_argument('--cov', type=check_extract_range, default=75.0, help='Override query coverage threshold (default = 75.0)')
    parser_extract.add_argument('--translate', action='store_true', help='Translate using standard code.')
    set_py_handler(parser_extract, 'extract', 'run')

    #---------------- Combine (Python) ----------------

    parser_combine = subparsers.add_parser('combine', help='Combine datasets from ReGAIN AMR and ReGAIN Curate')
    parser_combine.add_argument('--matrix1', required=True, help='Path to ReGAIN AMR data matrix file')
    parser_combine.add_argument('--matrix2', required=True, help='Path to ReGAIN Curate data matrix file')
    parser_combine.add_argument('--metadata1', required=True, help='Path to ReGAIN AMR metadata file')
    parser_combine.add_argument('--metadata2', required=True, help='Path to ReGAIN Curate metadata file')
    parser_combine.add_argument('--delete-duplicates', action='store_true', help='Delete duplicate values')
    set_py_handler(parser_combine, 'combine', 'run')

    #---------------- Collapse (Python) ----------------

    collapse_parser = subparsers.add_parser(
        'collapse-features',
        help="Collapse presence/absence matrix features using a user-provided Gene-to-bin mapping file"
    )

    collapse_parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input presence/absence matrix in CSV or TSV format'
    )

    collapse_parser.add_argument(
        '-M', '--metadata',
        required=True,
        help='Metadata/mapping file containing Gene and bin columns'
    )

    collapse_parser.add_argument(
        '-o', '--output-file',
        required=True,
        help='Output collapsed presence/absence matrix'
    )

    collapse_parser.add_argument(
        '--id-col',
        default=None,
        help='Presence/absence matrix genome/sample ID column. Default: first column'
    )

    collapse_parser.add_argument(
        '--gene-col',
        default='Gene',
        help='Column in metadata containing original feature name (generated from `regain matrix`). Default: Gene'
    )

    collapse_parser.add_argument(
        '--bin-col',
        default='bin',
        help='Column in metadata containing collapsed feature names. Default: bin'
    )

    collapse_parser.add_argument(
        '--drop-unmapped',
        action='store_true',
        help='Drop matrix features not present in metadata. Default: keep unmapped features unchanged'
    )

    collapse_parser.add_argument(
        '--missing-as-zero',
        action='store_true',
        help='Fill missing matrix values with 0 after writing missing-value report'
    )
    set_py_handler(collapse_parser, "collapse_features", 'run')

    #---------------- Matrix Summary (Python) ----------------

    summary_parser = subparsers.add_parser(
        'matrix-summary',
        help='Summarize and validate a ReGAIN presence/absence matrix'
    )

    summary_parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input presence/absence matrix'
    )

    summary_parser.add_argument(
        '-o', '--output-file',
        required=True,
        help='Output matrix summary CSV file'
    )

    summary_parser.add_argument(
        '--id-col',
        default=None,
        help='Presence/absence matrix genome/sample ID column. Default: first column (header required)'
    )

    summary_parser.add_argument(
        '--missing-as-zero',
        action='store_true',
        help=(
            "Fill missing matrix values with 0 after writing missing-value report.\n"
            "NOTE: Only use if values truly are 'absent'. Treating empirically missing values as absent will "
            "affect Bayesian network structure learning results")
    )
    set_py_handler(summary_parser, 'matrix_summary', 'run')

    #---------------- Genome Similarity (Python) ----------------
    genome_similarity_parser = subparsers.add_parser(
        'genome-similarity',
        help='Estimate genome similarity and report potential clonal groups'
    )

    genome_similarity_parser.add_argument(
        '-f', '--fasta-dir',
        required=True,
        help='Directory containing genome FASTA files'
    )

    genome_similarity_parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for genome similarity reports'
    )

    genome_similarity_parser.add_argument(
        '--method',
        choices=['fastani', 'mash'],
        default='fastani',
        help='Genome similarity method. Default: fastani'
    )

    genome_similarity_parser.add_argument(
        '--threshold',
        nargs='+',
        type=float,
        default=None,
        help=(
            'One or more space-separated thresholds. For FastANI, thresholds are ANI percent values '
            'such as 99.9 99.5 99.0. For Mash, thresholds are Mash distance values '
            'such as 0.001 0.005 0.01. Defaults are method-specific (values displayed here)'
        )
    )

    genome_similarity_parser.add_argument(
        '--genome-list',
        default=None,
        help=(
            'Optional file containing one genome filename or relative path per row. '
            'Entries are resolved relative to --fasta-dir'
        )
    )

    genome_similarity_parser.add_argument(
        '-e', '--extensions',
        nargs='+',
        default=['.fna', '.fa', '.ffn', '.fas', '.fasta'],
        help='FASTA extensions (space-separated) to include. Default: .fna, .fa, .fas, .ffn, .fas, .fasta'
    )

    genome_similarity_parser.add_argument(
        '--recursive',
        action='store_true',
        help='Search --fasta-dir recursively for FASTA files (one or more levels deep; use with caution)'
    )

    genome_similarity_parser.add_argument(
        '-T', '--threads',
        type=int,
        default=4,
        help='Threads to pass to selected method. Default: 4'
    )

    genome_similarity_parser.add_argument(
        '--sketch-size',
        type=int,
        default=10000,
        help='Mash sketch size. Used only with --method mash. Default: 10000'
    )

    genome_similarity_parser.add_argument(
        '--kmer-size',
        type=int,
        default=21,
        help='Mash k-mer size. Used only with --method mash. Default: 21'
    )
    set_py_handler(genome_similarity_parser, 'genome_similarity', 'run')

    #---------------- Network Analysis (Python) ----------------

    network_analysis_parser = subparsers.add_parser(
        'network-analysis',
        help='Compare two ReGAIN bnS/bnL results tables'
    )

    network_analysis_parser.add_argument(
        '--network1',
        required=True,
        help='Baseline ReGAIN bnS/bnL probability results table (Query_Results.csv)'
    )

    network_analysis_parser.add_argument(
        '--network2',
        required=True,
        help='Comparison ReGAIN bnS/bnL results table (Query_Results.csv)'
    )

    network_analysis_parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for network comparison results'
    )

    network_analysis_parser.add_argument(
        '--gene-a-col',
        default='Gene_1',
        help='Gene_1 column name. Default: Gene_1'
    )

    network_analysis_parser.add_argument(
        '--gene-b-col',
        default='Gene_2',
        help='Gene_2 column name. Default: Gene_2'
    )

    network_analysis_parser.add_argument(
        '--cpr-mean-col',
        default='Conditional_Probability_Mean',
        help='Conditional probability mean column name. Default: Conditional_Probability_Mean'
    )

    network_analysis_parser.add_argument(
        '--rr-mean-col',
        default='Relative_Risk_Mean',
        help='Relative risk mean column name. Default: Relative_Risk_Mean'
    )

    network_analysis_parser.add_argument(
        '--ard-mean-col',
        default='Absolute_Risk_Mean',
        help='Absolute risk difference column name. Default: Absolute_Risk_Mean'
    )

    network_analysis_parser.add_argument(
        '--status-metric',
        choices=['cpr', 'rr', 'ard'],
        default='ard',
        help='Metric used to classify retained/weakened/strengthened. Default: ard'
    )

    network_analysis_parser.add_argument(
        '--delta-threshold',
        type=float,
        default=0.1,
        help='Minimum mean delta magnitude used to classify weakened/strengthened. Default: 0.1'
    )
    set_py_handler(network_analysis_parser, 'network_analysis', 'run')

    # ---- parse & early exit for health ----
    args = parser.parse_args(argv)
    if getattr(args, "module_health", False):
        _print_module_health()
        return 0
    
    if getattr(args, "citation", False):
        print("""
              
\033[92m\nResistance Gene Association and Inference Network (ReGAIN) — please cite:\033[0m
            Bring Horvath E, Stein M, Mulvey MA, Hernandez EJ, Winter JM.
            Resistance Gene Association and Inference Network (ReGAIN): A Bioinformatics Pipeline for Assessing Probabilistic Co-Occurrence Between Resistance Genes in Bacterial Pathogens.
            bioRxiv 2024.02.26.582197; doi: https://doi.org/10.1101/2024.02.26.582197
              
\033[92mAMRFinder — please cite:\033[0m
            Feldgarden M, Brover V, Haft DH. et al.
            Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates.
            Antimicrob Agents Chemother. 2019 Oct 22;63(11):e00483-19. doi: 10.1128/AAC.00483-19
              
\033[92mAMRFinderPlus — please cite:\033[0m
            Feldgarden M, Brover V, Gonzalez-Escalona N. et al. 
            AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. 
            Sci Rep 11, 12728 (2021). https://doi.org/10.1038/s41598-021-91456-0
          
""".strip())
        return 0

    # ---- logging ----
    if getattr(args, "command", None):
        logging.basicConfig(filename=f'{args.command}_log.txt', level=logging.INFO)
        logging.info(f'Program started at {datetime.now()}')
        logging.info(f'User command: regain {" ".join(argv)}')

    # ---- dispatch ----
    if hasattr(args, "func"):
        try:
            args.func(args)
        finally:
            logging.info(f'Program ended at {datetime.now()}')
    else:
        parser.print_help()
        return 2

if __name__ == "__main__":
    sys.exit(main())
