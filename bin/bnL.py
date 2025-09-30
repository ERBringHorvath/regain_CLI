#!/usr/bin/env python3
import argparse, os, shutil, subprocess, sys

def run(cmd):
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

def main():
    ap = argparse.ArgumentParser(description="Run Bayesian Network analysis (large datasets, CPQuery)")
    ap.add_argument("-i","--input", required=True, help="Input matrix CSV")
    ap.add_argument("-M","--metadata", required=True, help="Metadata CSV")
    ap.add_argument("-o","--output_boot", required=True, help="Output .rds (bootstrapped network)")
    ap.add_argument("-T","--threads", type=int, default=2, help="Threads")
    ap.add_argument("-n","--number_of_bootstraps", required=True, type=int, help="Bootstrap number (e.g., 300–500)")
    ap.add_argument("-r","--number_of_resamples", required=True, type=int, help="Number of resamples for querying")
    ap.add_argument("-b","--blacklist", help="Optional blacklist CSV (no header): from,to")
    ap.add_argument("--iss", type=int, default=10, help="Imaginary sample size for BDe score (default: 10)")
    args, unknown = ap.parse_known_args()

    bash = shutil.which("bash")
    if not bash:
        sys.exit("ERROR: bash not found in PATH")

    script_dir = os.path.dirname(os.path.realpath(__file__))
    sh_path = os.path.join(script_dir, "bnCPQuery.sh")
    if not os.path.exists(sh_path):
        sys.exit(f"ERROR: cannot find {sh_path}")

    cmd = [
        bash, sh_path,
        "--input", args.input,
        "--metadata", args.metadata,
        "--output_boot", args.output_boot,
        "--threads", str(args.threads),
        "--number_of_bootstraps", str(args.number_of_bootstraps),
        "--resamples", str(args.number_of_resamples),
        "--iss", str(args.iss),
    ]
    if args.blacklist:
        cmd += ["--blacklist", args.blacklist]

    # Forward any unknown extra flags to R
    cmd += unknown

    run(cmd)

if __name__ == "__main__":
    main()