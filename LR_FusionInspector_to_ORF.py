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
import csv
import logging
import os
import re
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
    parser.add_argument(
        "--fusions_tsv",
        default=None,
        help="Optional FusionInspector fusions abridged TSV used to annotate whether predicted ORFs span the fusion breakpoint",
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
    if args.fusions_tsv and not os.path.exists(args.fusions_tsv):
        errors.append(f"--fusions_tsv: path not found: {args.fusions_tsv}")
    if errors:
        for e in errors:
            logger.error(e)
        sys.exit(1)


def parse_gtf_attributes(attr_text):
    attrs = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', attr_text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def parse_gff3_attributes(attr_text):
    attrs = {}
    for entry in attr_text.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            attrs[key] = value
    return attrs


def load_fusion_breakpoints(fusions_tsv):
    fusion_to_breakpoints = {}
    with open(fusions_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fusion_name = row["#FusionName"]
            breakpoint_pair = (
                int(row["LeftLocalBreakpoint"]),
                int(row["RightLocalBreakpoint"]),
            )
            fusion_to_breakpoints.setdefault(fusion_name, [])
            if breakpoint_pair not in fusion_to_breakpoints[fusion_name]:
                fusion_to_breakpoints[fusion_name].append(breakpoint_pair)
    return fusion_to_breakpoints


def load_lraa_transcripts(lraa_gtf):
    transcripts = {}
    with open(lraa_gtf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            feature_type = fields[2]
            if feature_type not in {"transcript", "exon"}:
                continue

            contig = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if transcript_id is None:
                continue

            transcript_entry = transcripts.setdefault(
                transcript_id,
                {
                    "fusion_name": contig,
                    "strand": strand,
                    "exons": [],
                },
            )
            if feature_type == "exon":
                transcript_entry["exons"].append((start, end))

    for transcript_entry in transcripts.values():
        transcript_entry["exons"].sort()

    return transcripts


def assign_breakpoint_to_transcript(exons, breakpoint_candidates):
    if len(exons) < 2:
        return None

    observed_junctions = {
        (exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)
    }

    matches = [bp for bp in breakpoint_candidates if bp in observed_junctions]
    if not matches:
        return None

    return matches[0]


def transcript_position_of_left_breakpoint(exons, left_breakpoint):
    transcript_pos = 0
    for exon_start, exon_end in exons:
        if exon_start <= left_breakpoint <= exon_end:
            return transcript_pos + (left_breakpoint - exon_start + 1)
        transcript_pos += exon_end - exon_start + 1
    return None


def load_orf_predictions(transdecoder_genome_gff3):
    orfs = {}
    with open(transdecoder_genome_gff3) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_gff3_attributes(fields[8])

            if feature_type == "mRNA":
                orf_id = attrs["ID"]
                transcript_id = re.sub(r"\.p\d+$", "", orf_id)
                orfs[orf_id] = {
                    "transcript_id": transcript_id,
                    "fusion_name": fields[0],
                    "strand": fields[6],
                    "cds_segments": [],
                }
            elif feature_type == "CDS":
                parent_id = attrs["Parent"]
                if parent_id in orfs:
                    orfs[parent_id]["cds_segments"].append((start, end))

    for orf_entry in orfs.values():
        orf_entry["cds_segments"].sort()

    return orfs


def load_orf_pep_annotations(transdecoder_pep):
    annotations = {}
    with open(transdecoder_pep) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            orf_id = header.split()[0]
            type_match = re.search(r"ORF type:([^(,]+)", header)
            annotations[orf_id] = {
                "orf_type": type_match.group(1).strip() if type_match else "",
            }
    return annotations


def interval_overlap_len(a_start, a_end, b_start, b_end):
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    return max(0, overlap_end - overlap_start + 1)


def annotate_orf_breakpoint_spanning(
    lraa_gtf, fusions_tsv, transdecoder_genome_gff3, transdecoder_pep, output_tsv
):
    fusion_breakpoints = load_fusion_breakpoints(fusions_tsv)
    transcripts = load_lraa_transcripts(lraa_gtf)
    orfs = load_orf_predictions(transdecoder_genome_gff3)
    pep_annotations = load_orf_pep_annotations(transdecoder_pep)

    rows = []
    for orf_id, orf_entry in sorted(orfs.items()):
        transcript_id = orf_entry["transcript_id"]
        transcript_entry = transcripts.get(transcript_id)

        row = {
            "fusion_name": orf_entry["fusion_name"],
            "transcript_id": transcript_id,
            "orf_id": orf_id,
            "ORF_type": "",
            "breakpoint_match_found": "False",
            "left_local_breakpoint": "",
            "right_local_breakpoint": "",
            "breakpoint_transcript_pos": "",
            "cds_bases_5p_of_breakpoint": "0",
            "cds_bases_3p_of_breakpoint": "0",
            "cds_spans_fusion_breakpoint": "False",
            "breakpoint_codon_phase": "",
        }

        pep_annotation = pep_annotations.get(orf_id)
        if pep_annotation:
            row["ORF_type"] = pep_annotation["orf_type"]

        if transcript_entry is None:
            row["note"] = "transcript_not_found_in_lraa_gtf"
            rows.append(row)
            continue

        breakpoint_candidates = fusion_breakpoints.get(transcript_entry["fusion_name"], [])
        matched_breakpoint = assign_breakpoint_to_transcript(
            transcript_entry["exons"], breakpoint_candidates
        )

        if matched_breakpoint is None:
            row["note"] = "no_matching_fusion_junction_for_transcript"
            rows.append(row)
            continue

        left_breakpoint, right_breakpoint = matched_breakpoint
        cds_left = sum(
            interval_overlap_len(cds_start, cds_end, 1, left_breakpoint)
            for cds_start, cds_end in orf_entry["cds_segments"]
        )
        cds_right = sum(
            interval_overlap_len(cds_start, cds_end, right_breakpoint, sys.maxsize)
            for cds_start, cds_end in orf_entry["cds_segments"]
        )
        spans_breakpoint = cds_left > 0 and cds_right > 0

        row.update(
            {
                "breakpoint_match_found": "True",
                "left_local_breakpoint": str(left_breakpoint),
                "right_local_breakpoint": str(right_breakpoint),
                "breakpoint_transcript_pos": str(
                    transcript_position_of_left_breakpoint(
                        transcript_entry["exons"], left_breakpoint
                    )
                    or ""
                ),
                "cds_bases_5p_of_breakpoint": str(cds_left),
                "cds_bases_3p_of_breakpoint": str(cds_right),
                "cds_spans_fusion_breakpoint": str(spans_breakpoint),
                "breakpoint_codon_phase": str(cds_left % 3) if spans_breakpoint else "",
                "note": "",
            }
        )

        rows.append(row)

    fieldnames = [
        "fusion_name",
        "transcript_id",
        "orf_id",
        "ORF_type",
        "breakpoint_match_found",
        "left_local_breakpoint",
        "right_local_breakpoint",
        "breakpoint_transcript_pos",
        "cds_bases_5p_of_breakpoint",
        "cds_bases_3p_of_breakpoint",
        "cds_spans_fusion_breakpoint",
        "breakpoint_codon_phase",
        "note",
    ]
    with open(output_tsv, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    spanning_count = sum(
        row["cds_spans_fusion_breakpoint"] == "True" for row in rows
    )
    logger.info(
        "Annotated %d ORFs for fusion-breakpoint spanning; %d span the fusion breakpoint",
        len(rows),
        spanning_count,
    )


def main():
    args = parse_args()

    # Resolve to absolute paths before any chdir
    args.genome = os.path.abspath(args.genome)
    args.bam = os.path.abspath(args.bam)
    args.gtf = os.path.abspath(args.gtf)
    if args.fusions_tsv:
        args.fusions_tsv = os.path.abspath(args.fusions_tsv)
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

    if args.fusions_tsv:
        annotate_orf_breakpoint_spanning(
            lraa_gtf=lraa_gtf,
            fusions_tsv=args.fusions_tsv,
            transdecoder_genome_gff3=os.path.join(
                output_dir, "LRAA.cDNA.fasta.transdecoder.genome.gff3"
            ),
            transdecoder_pep=os.path.join(
                output_dir, "LRAA.cDNA.fasta.transdecoder.pep"
            ),
            output_tsv=os.path.join(output_dir, "fusion_orf_breakpoint_annotations.tsv"),
        )

    logger.info("Pipeline complete. Results in: %s", output_dir)


if __name__ == "__main__":
    main()
