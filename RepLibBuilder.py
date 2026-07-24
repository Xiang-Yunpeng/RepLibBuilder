#!/usr/bin/env python3
import argparse
import os
import sys
import ctypes

# Automatically determine the absolute path to the RepLibBuilder root directory
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MAPPING = os.path.join(BASE_DIR, 'modules', 'rating_module', 'te_mapping_rules.tsv')

# Single source of truth for the software version.
__version__ = "1.0.1"

# Import modules from the new rating_module directory
from modules.rating_module.repbase_fixer import RepbaseFixer
from modules.rating_module.dfam_fetcher import DfamFetcher
from modules.rating_module.evaluator import run_evaluate_pipeline
from modules.rating_module.lib_merger import RawLibraryMerger
from modules.rating_module.html_reporter import generate_html
from modules.rating_module.normalizer import TEDictionary
from utils.common import check_file, ensure_dir

# Import modules from the build directory
from modules.build.library_filter import LibraryFilter
from modules.build.cdhit_wrapper import CDHITWrapper

def setup_parsers():
    parser = argparse.ArgumentParser(
        description="""
==============================================================================
RepLibBuilder v1.0.1: Automated Construction of Species-Specific Repeat Libraries
==============================================================================
Two-Step Workflow:
  Step 1 (evaluate): Extract, standardize, diagnose, and merge public/denovo 
                     libraries. Generates comprehensive statistics for review.
  Step 2 (build):    Filter sequences based on diagnostic tags and perform 
                     CD-HIT clustering to finalize the custom library.
""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--version', action='version', version=f'RepLibBuilder v{__version__}')
    subparsers = parser.add_subparsers(dest='mode', required=True, help="Select operation mode")

    # ======================================================
    # Mode 1: Evaluate (Step 1)
    # ======================================================
    parser_eval = subparsers.add_parser(
        'evaluate', 
        help="Step 1: Evaluate, standardize, and merge all input libraries.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Public DBs
    pub_group = parser_eval.add_argument_group('Public Database Inputs')
    pub_group.add_argument('--repbase', metavar='FILE', help="Path to a Repbase FASTA export (raw GIRI format).")
    pub_group.add_argument('--dfam_db', metavar='DIR', help="Path to a Dfam FamDB directory (folder of partition .h5 files; root *.0.h5 required).")
    pub_group.add_argument('--tax_id', metavar='ID', help="NCBI Taxonomy ID (Required if --dfam_db is used).")

    # De novo DBs
    denovo_group = parser_eval.add_argument_group('De Novo Prediction Inputs')
    denovo_group.add_argument("--hite", metavar='FILE', help="Path to HiTE output FASTA file.")
    denovo_group.add_argument("--edta", metavar='FILE', help="Path to EDTA output FASTA file.")
    denovo_group.add_argument("--rm2", metavar='FILE', help="Path to RepeatModeler2 output FASTA file.")
    
    # Config
    config_group = parser_eval.add_argument_group('Configuration')
    config_group.add_argument("-o", "--out_dir", required=True, metavar='DIR', help="Base output directory.")
    config_group.add_argument("-m", "--mapping", default=DEFAULT_MAPPING, metavar='FILE', help="Do not enable this parameter unless you want to define your own dictionary.")
    config_group.add_argument("-t", "--threads", default=4, type=int, metavar='INT', help="Threads for TEsorter (Default: 4).")
    config_group.add_argument("-db", "--te_db", default='rexdb-metazoa', metavar='NAME', 
                        choices=['gydb', 'rexdb', 'rexdb-plant', 'rexdb-metazoa', 'rexdb-v3', 'rexdb-plantv3', 'rexdb-metazoav3', 'rexdb-pnas', 'rexdb-line', 'sine'], 
                        help="HMM database for TEsorter (Default: rexdb-metazoa).")

    # ======================================================
    # Mode 2: Build (Step 2)
    # ======================================================
    parser_build = subparsers.add_parser(
        'build', 
        help="Step 2: Filter and cluster the library based on Step 1 results.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    build_in_group = parser_build.add_argument_group('Input and Output')
    build_in_group.add_argument("-i", "--input", required=True, metavar='FILE', help="Input merged raw FASTA from Step 1 (e.g., step1_combined_raw.fa)")
    build_in_group.add_argument("-o", "--out_dir", required=True, metavar='DIR', help="Base output directory for Step 2.")
    build_in_group.add_argument("-m", "--mapping", default=DEFAULT_MAPPING, metavar='FILE', help="Path to TE mapping dictionary (Default is built-in).")
    build_in_group.add_argument("-t", "--threads", type=int, default=8, metavar='INT', help="Number of threads for parallel CD-HIT clustering (Default: 8).")
    
    filter_group = parser_build.add_argument_group('Filtering Options')
    filter_group.add_argument("--filter_db", help="Filter by DB and Class (e.g., RM2:Unknown,EDTA:Unknown)")
    filter_group.add_argument("--filter_tag", help="Filter by diagnostic tag (e.g., Unverified,Exempted)")
    filter_group.add_argument("--filter_length", help="Filter by length (e.g., EDTA:SINE:max:1000,RM2:DNA:min:50)")

    cluster_group = parser_build.add_argument_group('CD-HIT Clustering Options')
    cluster_group.add_argument("-c", type=float, default=0.8, help="Sequence identity threshold for CD-HIT (Default: 0.8)")
    cluster_group.add_argument("-aL", type=float, default=0.8, help="Alignment coverage for the longer sequence (Default: 0.8)")
    cluster_group.add_argument("-aS", type=float, default=0.8, help="Alignment coverage for the shorter sequence (Default: 0.8)")
    cluster_group.add_argument("--no_stratify", action="store_true",
                               help="Disable stratified (per-class) clustering; cluster all TE sequences in a\n"
                                    "single CD-HIT run instead of one run per TE class. Non-TE / structural\n"
                                    "sequences (satellites, simple/low-complexity repeats, snRNA/tRNA/rRNA/\n"
                                    "scRNA/etc.) are held out from clustering in either mode. -c/-aL/-aS still apply.")
    cluster_group.add_argument("--dedup_unknown", action="store_true",
                               help="Opt-in: after clustering, drop fully-Unknown consensus that a classified\n"
                                    "consensus already covers (cd-hit-est-2d, single pairwise). Raises the\n"
                                    "share of classified sequences at near-zero coverage cost. Off by default.")
    cluster_group.add_argument("--dedup_unknown_id", type=float, default=0.8,
                               help="Identity threshold (-c) for --dedup_unknown, applied to the shorter\n"
                                    "(Unknown) sequence. Must be within [0.75, 1.0]. Default: 0.8.")
    cluster_group.add_argument("--dedup_unknown_cov", type=float, default=0.8,
                               help="Coverage threshold (-aS) on the Unknown sequence for --dedup_unknown.\n"
                                    "Default: 0.8.")
    cluster_group.add_argument("--dedup_unknown_aL", type=float, default=0.0,
                               help="Coverage threshold (-aL) on the LONGER sequence for --dedup_unknown.\n"
                                    "Default: 0.0 (off). Raise (e.g. 0.8) for a stricter, more symmetric\n"
                                    "containment test that spares small Unknown fragments of large elements.")

    return parser

def set_process_name(name):
    """Set the process name in system monitor (max 15 chars)."""
    try:
        libc = ctypes.cdll.LoadLibrary('libc.so.6')
        b_name = name.encode('utf-8')[:15]
        libc.prctl(15, b_name, 0, 0, 0)
    except Exception:
        pass

def create_evaluate_directory_tree(base_dir):
    """Creates the standard directory structure for Step 1."""
    out_dir = os.path.abspath(base_dir)
    dirs = {
        'pub_repbase': os.path.join(out_dir, "01.public_db", "Repbase"),
        'pub_dfam':    os.path.join(out_dir, "01.public_db", "Dfam"),
        'denovo_hite': os.path.join(out_dir, "02.denovo_db", "HiTE"),
        'denovo_edta': os.path.join(out_dir, "02.denovo_db", "EDTA"),
        'denovo_rm2':  os.path.join(out_dir, "02.denovo_db", "RM2"),
        'merge_db':    os.path.join(out_dir, "03.merge_db"),
        'report':      os.path.join(out_dir, "04.report")
    }
    for path in dirs.values():
        os.makedirs(path, exist_ok=True)
    return dirs

def create_build_directory_tree(base_dir):
    """Creates the standard directory structure for Step 2."""
    out_dir = os.path.abspath(base_dir)
    dirs = {
        'filter':  os.path.join(out_dir, "01.filter"),
        'cluster': os.path.join(out_dir, "02.cluster")
    }
    for path in dirs.values():
        os.makedirs(path, exist_ok=True)
    return dirs

def run_evaluate(args):
    """Executes Step 1 pipeline."""
    dirs = create_evaluate_directory_tree(args.out_dir)
    merged_inputs = {}

    print(f"=== RepLibBuilder v{__version__} (Step 1: Evaluate) Started ===")

    # Pre-load the TE dictionary once; reused by RepbaseFixer (and conceptually
    # available for any other module in this pipeline that needs label
    # normalization). Loading here avoids re-parsing the mapping file later.
    te_dict = TEDictionary(args.mapping)

    # 1. Process Public Databases
    if args.repbase:
        print("\n>>> Processing Repbase...")
        fixer = RepbaseFixer(args.repbase, dirs['pub_repbase'], te_dict)
        merged_inputs['repbase'] = fixer.run()
    
    if args.dfam_db:
        if not args.tax_id:
            print("[Error] Dfam is enabled but missing --tax_id configuration.")
            sys.exit(1)
        print("\n>>> Processing Dfam...")
        fetcher = DfamFetcher(args.dfam_db, args.tax_id, None, dirs['pub_dfam'], te_dict)
        merged_inputs['dfam'] = fetcher.run()

    # 2. Process De Novo Databases
    if args.hite or args.edta or args.rm2:
        print("\n>>> Evaluating De Novo Predictions...")
        denovo_results = run_evaluate_pipeline(args, dirs)
        merged_inputs.update(denovo_results)
    else:
        print("\n>>> No De Novo libraries provided, skipping evaluation.")

    # 3. Merge raw data and 4. Generate Report
    if merged_inputs:
        print("\n>>> Pooling standard libraries...")
        merger = RawLibraryMerger(merged_inputs, dirs['merge_db'], dirs['report'])
        merger.run()

        print("\n>>> Generating Interactive HTML Report...")
        eval_json_path = os.path.join(dirs['report'], "evaluation_report.json")
        merge_json_path = os.path.join(dirs['report'], "merge_stats.json")
        out_html_path = os.path.join(dirs['report'], "evaluation_report.html")
        
        generate_html(eval_json_path, merge_json_path, out_html_path)

    else:
        print("[Error] No valid inputs provided to process.")

    print("\n" + "="*60)
    print(f"Step 1 Completed Successfully!")
    print(f"Directory tree created at: {os.path.abspath(args.out_dir)}")
    print("="*60)

def run_build(args):
    """Executes Step 2 pipeline."""
    dirs = create_build_directory_tree(args.out_dir)
    
    print(f"=== RepLibBuilder v{__version__} (Step 2: Build) Started ===")
    
    # Check input file
    input_fasta = check_file(args.input, "Merged raw FASTA from Step 1")
    
    # 1. Load TE Dictionary
    te_dict = TEDictionary(args.mapping)
    
    # 2. Execute Filtering
    print("\n>>> Phase 1: Filtering Library...")
    filtered_fasta_path = os.path.join(dirs['filter'], "step2_filtered.fa")
    
    filter_module = LibraryFilter(
        input_fasta=input_fasta,
        output_fasta=filtered_fasta_path,
        te_dict=te_dict,
        filter_db=args.filter_db,
        filter_tag=args.filter_tag,
        filter_length=args.filter_length
    )
    filter_module.run()

    # 3. Execute Clustering
    print("\n>>> Phase 2: Clustering Library (CD-HIT)...")
    cdhit_module = CDHITWrapper(
        input_fasta=filtered_fasta_path,
        output_dir=dirs['cluster'],
        te_dict=te_dict,
        c=args.c,
        aL=args.aL,
        aS=args.aS,
        threads=args.threads,
        stratify=not args.no_stratify,
        dedup_unknown=args.dedup_unknown,
        dedup_id=args.dedup_unknown_id,
        dedup_cov=args.dedup_unknown_cov,
        dedup_aL=args.dedup_unknown_aL
    )
    cdhit_module.run()

    print("\n" + "="*60)
    print(f"Step 2 Progress saved to: {os.path.abspath(args.out_dir)}")
    print("="*60)

def main():
    set_process_name("RepLibBuilder")
    parser = setup_parsers()
    args = parser.parse_args()

    if args.mode == 'evaluate':
        run_evaluate(args)
    elif args.mode == 'build':
        run_build(args)

if __name__ == "__main__":
    main()
