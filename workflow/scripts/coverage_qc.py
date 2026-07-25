#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import os
import statistics
import subprocess
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate raw FASTQ-equivalent coverage, aligned genome-wide "
            "coverage, and CpG coverage metrics."
        )
    )

    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier."
    )

    parser.add_argument(
        "--cpg",
        required=True,
        help="Gzip-compressed methylKit CpG file."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Coordinate-sorted and indexed final BAM file."
    )

    parser.add_argument(
        "--fastp-json",
        required=True,
        help="fastp JSON report."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output coverage-QC TSV."
    )

    parser.add_argument(
        "--genome-size",
        required=True,
        type=int,
        help=(
            "Genome-size denominator used for vendor-style raw coverage."
        )
    )

    parser.add_argument(
        "--min-mapq",
        required=True,
        type=int,
        help="Minimum mapping quality used by samtools depth."
    )

    parser.add_argument(
        "--min-baseq",
        required=True,
        type=int,
        help="Minimum base quality used by samtools depth."
    )

    parser.add_argument(
        "--excluded-contigs",
        default="lambda,pUC19",
        help=(
            "Comma-separated contigs excluded from aligned human coverage."
        )
    )

    return parser.parse_args()


def format_float(
    value: Optional[float],
    decimals: int = 6
) -> str:
    if value is None:
        return "NA"

    return f"{value:.{decimals}f}"


def run_command(command: Sequence[str]) -> str:
    completed = subprocess.run(
        list(command),
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if completed.stderr:
        print(completed.stderr, end="")

    return completed.stdout


def calculate_raw_fastq_coverage(
    fastp_json: str,
    genome_size: int
) -> Dict[str, object]:
    """
    Calculate vendor-style raw coverage from fastp before-filtering bases.

    fastp summary.before_filtering.total_bases represents the bases present
    before fastp trimming and filtering.
    """
    if genome_size <= 0:
        raise ValueError(
            f"Genome-size denominator must be positive: {genome_size}"
        )

    with open(fastp_json, encoding="utf-8") as fh:
        fastp = json.load(fh)

    before_filtering = (
        fastp
        .get("summary", {})
        .get("before_filtering", {})
    )

    total_bases = before_filtering.get("total_bases")

    if total_bases is None:
        raise ValueError(
            "fastp JSON does not contain "
            "summary.before_filtering.total_bases"
        )

    try:
        total_bases = int(total_bases)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Invalid fastp before-filtering total_bases value: "
            f"{total_bases}"
        ) from exc

    raw_fastq_coverage = total_bases / genome_size

    return {
        "vendor_genome_size_denominator": genome_size,
        "raw_fastq_coverage": raw_fastq_coverage,
    }


def calculate_cpg_metrics(cpg_file: str) -> Dict[str, object]:
    """
    Calculate metrics for CpGs represented in a methylKit file.

    Expected methylKit columns:
        0: chrBase
        1: chr
        2: base
        3: strand
        4: coverage
        5: freqC
        6: freqT
    """
    coverages: List[int] = []

    methylated_sum = 0.0
    unmethylated_sum = 0.0

    skipped_lines = 0

    with gzip.open(cpg_file, "rt") as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 7:
                skipped_lines += 1
                print(
                    f"Skipping CpG line {line_number}: "
                    "fewer than seven columns"
                )
                continue

            try:
                coverage = int(fields[4])
                freq_c = float(fields[5])
                freq_t = float(fields[6])
            except ValueError:
                # This also skips a methylKit header row.
                continue

            if coverage < 0:
                skipped_lines += 1
                print(
                    f"Skipping CpG line {line_number}: "
                    f"negative coverage {coverage}"
                )
                continue

            coverages.append(coverage)

            # These counts are reconstructed from coverage and rounded
            # percentages, so they are approximate.
            methylated_sum += coverage * freq_c / 100.0
            unmethylated_sum += coverage * freq_t / 100.0

    mean_coverage: Optional[float] = None
    median_coverage: Optional[float] = None

    if coverages:
        mean_coverage = statistics.mean(coverages)
        median_coverage = statistics.median(coverages)

    total_observations = methylated_sum + unmethylated_sum

    methylation_fraction: Optional[float] = None

    if total_observations > 0:
        methylation_fraction = (
            methylated_sum / total_observations
        )

    print(f"CpG lines skipped: {skipped_lines}")

    return {
        "cpg_sites_called": len(coverages),
        "approx_total_methylated_counts": round(methylated_sum),
        "approx_total_unmethylated_counts": round(unmethylated_sum),
        "mean_coverage_called_cpgs": mean_coverage,
        "median_coverage_called_cpgs": median_coverage,
        "coverage_weighted_methylation_fraction":
            methylation_fraction,
    }


def get_included_reference_contigs(
    bam_file: str,
    excluded_contigs: Sequence[str]
) -> List[Tuple[str, int]]:
    """
    Get reference contigs and lengths from the BAM header, excluding
    specified spike-in contigs.
    """
    excluded = set(excluded_contigs)

    header = run_command(
        [
            "samtools",
            "view",
            "-H",
            bam_file
        ]
    )

    included_contigs: List[Tuple[str, int]] = []

    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue

        contig_name: Optional[str] = None
        contig_length: Optional[int] = None

        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                contig_name = field[3:]
            elif field.startswith("LN:"):
                contig_length = int(field[3:])

        if contig_name is None or contig_length is None:
            continue

        if contig_name in excluded:
            continue

        included_contigs.append(
            (contig_name, contig_length)
        )

    if not included_contigs:
        raise RuntimeError(
            "No reference contigs remained after excluding: "
            + ",".join(excluded_contigs)
        )

    return included_contigs


def write_contig_bed(
    contigs: Sequence[Tuple[str, int]],
    bed_file: str
) -> None:
    with open(bed_file, "w", encoding="utf-8") as out:
        for contig_name, contig_length in contigs:
            out.write(
                f"{contig_name}\t0\t{contig_length}\n"
            )


def calculate_aligned_base_coverage(
    bam_file: str,
    min_mapq: int,
    min_baseq: int,
    excluded_contigs: Sequence[str]
) -> Dict[str, object]:
    """
    Calculate mean base coverage across included reference contigs.

    samtools depth options:
        -aa: emit zero-coverage positions
        -s: count overlapping paired-end mates once
        -q: minimum base quality
        -Q: minimum mapping quality
        -b: restrict analysis to included reference intervals
    """
    contigs = get_included_reference_contigs(
        bam_file=bam_file,
        excluded_contigs=excluded_contigs
    )

    expected_reference_size = sum(
        length
        for _, length in contigs
    )

    if expected_reference_size <= 0:
        raise RuntimeError(
            "Aligned reference denominator is zero"
        )

    with tempfile.TemporaryDirectory(
        prefix="coverage_qc_"
    ) as temporary_directory:
        bed_file = os.path.join(
            temporary_directory,
            "included_reference_contigs.bed"
        )

        write_contig_bed(
            contigs=contigs,
            bed_file=bed_file
        )

        command = [
            "samtools",
            "depth",
            "-aa",
            "-s",
            "-q",
            str(min_baseq),
            "-Q",
            str(min_mapq),
            "-b",
            bed_file,
            bam_file
        ]

        print(
            "Running aligned-depth command: "
            + " ".join(command)
        )

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024
        )

        if process.stdout is None:
            raise RuntimeError(
                "Could not read samtools depth output"
            )

        observed_reference_positions = 0
        total_aligned_depth = 0

        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                continue

            try:
                depth = int(fields[2])
            except ValueError:
                continue

            observed_reference_positions += 1
            total_aligned_depth += depth

        stderr_text = ""

        if process.stderr is not None:
            stderr_text = process.stderr.read()

        return_code = process.wait()

        if stderr_text:
            print(stderr_text, end="")

        if return_code != 0:
            raise subprocess.CalledProcessError(
                return_code,
                command,
                stderr=stderr_text
            )

    if observed_reference_positions != expected_reference_size:
        raise RuntimeError(
            "Unexpected number of reference positions from samtools depth: "
            f"expected={expected_reference_size}, "
            f"observed={observed_reference_positions}"
        )

    mean_aligned_base_coverage = (
        total_aligned_depth / expected_reference_size
    )

    return {
        "mean_aligned_base_coverage":
            mean_aligned_base_coverage,
        "aligned_reference_size_denominator":
            expected_reference_size,
    }


def write_output(
    output_file: str,
    metrics: Dict[str, object]
) -> None:
    fieldnames = [
        "vendor_genome_size_denominator",
        "raw_fastq_coverage",
        "mean_aligned_base_coverage",
        "aligned_reference_size_denominator",
        "minimum_mapping_quality",
        "minimum_base_quality",
        "overlapping_mates_counted_once",
        "excluded_depth_contigs",
        "cpg_sites_called",
        "approx_total_methylated_counts",
        "approx_total_unmethylated_counts",
        "mean_coverage_called_cpgs",
        "median_coverage_called_cpgs",
        "coverage_weighted_methylation_fraction",
    ]

    output_directory = os.path.dirname(output_file)

    if output_directory:
        os.makedirs(
            output_directory,
            exist_ok=True
        )

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8"
    ) as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="raise"
        )

        writer.writeheader()
        writer.writerow(metrics)


def main() -> None:
    args = parse_arguments()

    excluded_contigs = [
        contig.strip()
        for contig in args.excluded_contigs.split(",")
        if contig.strip()
    ]

    print(f"Sample: {args.sample}")
    print(f"CpG input: {args.cpg}")
    print(f"BAM input: {args.bam}")
    print(f"fastp JSON: {args.fastp_json}")
    print(f"Output: {args.output}")
    print(
        "Excluded aligned-depth contigs: "
        + ",".join(excluded_contigs)
    )

    raw_metrics = calculate_raw_fastq_coverage(
        fastp_json=args.fastp_json,
        genome_size=args.genome_size
    )

    cpg_metrics = calculate_cpg_metrics(
        cpg_file=args.cpg
    )

    aligned_metrics = calculate_aligned_base_coverage(
        bam_file=args.bam,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        excluded_contigs=excluded_contigs
    )

    metrics: Dict[str, object] = {
        "vendor_genome_size_denominator":
            raw_metrics["vendor_genome_size_denominator"],

        "raw_fastq_coverage":
            format_float(
                raw_metrics["raw_fastq_coverage"]
            ),

        "mean_aligned_base_coverage":
            format_float(
                aligned_metrics[
                    "mean_aligned_base_coverage"
                ]
            ),

        "aligned_reference_size_denominator":
            aligned_metrics[
                "aligned_reference_size_denominator"
            ],

        "minimum_mapping_quality":
            args.min_mapq,

        "minimum_base_quality":
            args.min_baseq,

        "overlapping_mates_counted_once":
            "true",

        "excluded_depth_contigs":
            ",".join(excluded_contigs),

        "cpg_sites_called":
            cpg_metrics["cpg_sites_called"],

        "approx_total_methylated_counts":
            cpg_metrics[
                "approx_total_methylated_counts"
            ],

        "approx_total_unmethylated_counts":
            cpg_metrics[
                "approx_total_unmethylated_counts"
            ],

        "mean_coverage_called_cpgs":
            format_float(
                cpg_metrics[
                    "mean_coverage_called_cpgs"
                ]
            ),

        "median_coverage_called_cpgs":
            format_float(
                cpg_metrics[
                    "median_coverage_called_cpgs"
                ]
            ),

        "coverage_weighted_methylation_fraction":
            format_float(
                cpg_metrics[
                    "coverage_weighted_methylation_fraction"
                ]
            ),
    }

    write_output(
        output_file=args.output,
        metrics=metrics
    )

    print("")
    print("Coverage QC summary")
    print("-------------------")
    print(
        "Raw FASTQ coverage: "
        f"{metrics['raw_fastq_coverage']}x"
    )
    print(
        "Mean aligned base coverage: "
        f"{metrics['mean_aligned_base_coverage']}x"
    )
    print(
        "Called CpG sites: "
        f"{metrics['cpg_sites_called']}"
    )
    print(
        "Mean called-CpG coverage: "
        f"{metrics['mean_coverage_called_cpgs']}x"
    )
    print(
        "Median called-CpG coverage: "
        f"{metrics['median_coverage_called_cpgs']}x"
    )


if __name__ == "__main__":
    main()
