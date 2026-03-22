#!/usr/bin/env python3

"""
LR_FusionInspector_to_ORF.py

Pipeline:
  1. Run LRAA on FusionInspector junction read alignments to reconstruct fusion
     transcript isoforms.
  2. Extract reference peptides from FusionInspector annotations (finspector.gtf)
     for use as BLAST homology targets in TransDecoder.
  3. Run TransDecoder on the LRAA-assembled transcripts to predict ORFs.
"""

import argparse
import logging
import os
import subprocess
import sys


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

LRAA_DEFAULT = os.path.join(_SCRIPT_DIR, "LRAA", "LRAA")
GTF_TO_FEATURE_SEQS_DEFAULT = os.path.join(
    _SCRIPT_DIR, "LRAA", "util", "misc", "gtf_file_to_feature_seqs.pl"
)
TRANSDECODER_DEFAULT = os.path.join(_SCRIPT_DIR, "TransDecoder", "TransDecoder")


def run_cmd(cmd, description):
    logger.info("Running: %s", description)
    logger.info("CMD: %s", " ".join(cmd))
    result = subprocess.run(cmd, check=True)
    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required inputs
    parser.add_argument(
        "--genome",
        required=True,
        help="FusionInspector genome FASTA (finspector.fa)",
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="FusionInspector junction reads BAM (finspector.junction_reads.bam)",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="FusionInspector annotation GTF (finspector.gtf) used to extract reference peptides",
    )

    # Output
    parser.add_argument(
        "--output_dir",
        default="LR_FusionInspector_to_ORF_outdir",
        help="Output directory (default: %(default)s)",
    )

    # Tool paths
    parser.add_argument(
        "--lraa",
        default=LRAA_DEFAULT,
        help="Path to LRAA executable (default: %(default)s)",
    )
    parser.add_argument(
        "--gtf_to_feature_seqs",
        default=GTF_TO_FEATURE_SEQS_DEFAULT,
        help="Path to gtf_file_to_feature_seqs.pl (default: %(default)s)",
    )
    parser.add_argument(
        "--transdecoder",
        default=TRANSDECODER_DEFAULT,
        help="Path to TransDecoder executable (default: %(default)s)",
    )

    # LRAA options
    parser.add_argument(
        "--HiFi",
        action="store_true",
        default=False,
        help="Pass --HiFi flag to LRAA (high-fidelity mode) (default: False)",
    )
    parser.add_argument(
        "--num_parallel_contigs",
        type=int,
        default=None,
        help="Pass --num_parallel_contigs to LRAA",
    )

    # TransDecoder options
    parser.add_argument(
        "--blast_threads",
        type=int,
        default=1,
        help="Threads for TransDecoder homology search (default: %(default)s)",
    )

    return parser.parse_args()


def validate_paths(args):
    errors = []
    for label, path in [
        ("--genome", args.genome),
        ("--bam", args.bam),
        ("--gtf", args.gtf),
        ("--lraa", args.lraa),
        ("--gtf_to_feature_seqs", args.gtf_to_feature_seqs),
        ("--transdecoder", args.transdecoder),
    ]:
        if not os.path.exists(path):
            errors.append(f"{label}: path not found: {path}")
    if errors:
        for e in errors:
            logger.error(e)
        sys.exit(1)


def main():
    args = parse_args()

    # Resolve to absolute paths before any chdir
    args.genome = os.path.abspath(args.genome)
    args.bam = os.path.abspath(args.bam)
    args.gtf = os.path.abspath(args.gtf)
    args.lraa = os.path.abspath(os.path.expanduser(args.lraa))
    args.gtf_to_feature_seqs = os.path.abspath(
        os.path.expanduser(args.gtf_to_feature_seqs)
    )
    args.transdecoder = os.path.abspath(os.path.expanduser(args.transdecoder))

    validate_paths(args)

    os.makedirs(args.output_dir, exist_ok=True)
    output_dir = os.path.abspath(args.output_dir)

    os.chdir(output_dir)
    logger.info("Working directory: %s", output_dir)

    # -------------------------------------------------------------------------
    # Step 1: LRAA - reconstruct fusion transcript isoforms
    # -------------------------------------------------------------------------
    lraa_gtf = os.path.join(output_dir, "LRAA.gtf")

    lraa_cmd = [
        args.lraa,
        "--genome",
        args.genome,
        "--bam",
        args.bam,
        "--ME_only",
        "--min_reads_novel",
        "1",
    ]
    if args.HiFi:
        lraa_cmd.append("--HiFi")
    if args.num_parallel_contigs is not None:
        lraa_cmd += ["--num_parallel_contigs", str(args.num_parallel_contigs)]

    run_cmd(lraa_cmd, "LRAA: reconstruct fusion transcript isoforms")

    if not os.path.exists(lraa_gtf):
        logger.error("LRAA did not produce expected output: %s", lraa_gtf)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Step 2: Extract reference peptides from FusionInspector annotations
    # -------------------------------------------------------------------------
    ref_pep = os.path.join(output_dir, "finspector.gtf.pep")

    with open(ref_pep, "w") as pep_fh:
        logger.info("Running: extract reference peptides from finspector.gtf")
        gtf_pep_cmd = [
            "perl",
            args.gtf_to_feature_seqs,
            "--gtf_file",
            args.gtf,
            "--genome_fa",
            args.genome,
            "--seqType",
            "prot",
        ]
        logger.info("CMD: %s > %s", " ".join(gtf_pep_cmd), ref_pep)
        subprocess.run(gtf_pep_cmd, stdout=pep_fh, check=True)

    if os.path.getsize(ref_pep) == 0:
        logger.error("Reference peptide file is empty: %s", ref_pep)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Step 3: TransDecoder - predict ORFs in LRAA-assembled transcripts
    # -------------------------------------------------------------------------
    transdecoder_cmd = [
        args.transdecoder,
        "--genome",
        args.genome,
        "--gtf",
        lraa_gtf,
        "--blast_search_pep",
        ref_pep,
        "--blast_threads",
        str(args.blast_threads),
        "--no_refine_starts",
        "--single_best_only",
    ]

    run_cmd(transdecoder_cmd, "TransDecoder: predict ORFs in fusion transcripts")

    logger.info("Pipeline complete. Results in: %s", output_dir)


if __name__ == "__main__":
    main()
