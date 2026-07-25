#!/usr/bin/env python3

import csv
import gzip
import json
import os
import statistics
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple


def format_float(
    value: Optional[float],
    decimals: int = 6,
) -> str:
    """Format a numeric QC value or return NA."""
    if value is None:
        return "NA"

    return f"{value:.{decimals}f}"


def read_first_tsv_row(filename: str) -> Dict[str, str]:
    """Read the first data row from a tab-delimited file."""
    with open(filename, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            return dict(row)

    return {}


def read_fastp_metrics(
    fastp_json: str,
    genome_size: int,
) -> Dict[str, object]:
    """
    Calculate vendor-style raw sequencing coverage from fastp.

    fastp summary.before_filtering.total_bases represents the bases supplied
    to fastp before trimming and quality filtering.
    """
    with open(fastp_json) as fh:
        fastp = json.load(fh)

    before = fastp.get("summary", {}).get("before_filtering", {})

    raw_total_bases = before.get("total_bases")
    raw_total_reads = before.get("total_reads")
    read1_mean_length = before.get("read1_mean_length")
    read2_mean_length = before.get("read2_mean_length")

    raw_fastq_coverage: Optional[float] = None

    try:
        raw_fastq_coverage = float(raw_total_bases) / genome_size
    except (TypeError, ValueError, ZeroDivisionError):
        pass

    return {
        "raw_fastq_coverage": raw_fastq_coverage,
        "raw_total_bases": raw_total_bases,
        "raw_total_reads": raw_total_reads,
        "read1_mean_length": read1_mean_length,
        "read2_mean_length": read2_mean_length,
    }


def calculate_cpg_metrics(cpg_file: str) -> Dict[str, object]:
    """
    Calculate CpG metrics from a methylKit-format file.

    Expected columns:
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

    with gzip.open(cpg_file, "rt") as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 7:
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
                # This also skips a header row when one is present.
                continue

            if coverage < 0:
                print(
                    f"Skipping CpG line {line_number}: "
                    f"negative coverage {coverage}"
                )
                continue

            coverages.append(coverage)

            # methylKit stores percentages. These reconstructed totals are
            # approximate if freqC and freqT have been rounded.
            methylated_sum += coverage * freq_c / 100.0
            unmethylated_sum += coverage * freq_t / 100.0

    cpg_sites_called = len(coverages)

    mean_coverage: Optional[float] = None
    median_coverage: Optional[float] = None

    if coverages:
        mean_coverage = statistics.mean(coverages)
        median_coverage = statistics.median(coverages)

    total_observations = methylated_sum + unmethylated_sum

    methylation_fraction: Optional[float] = None

    if total_observations > 0:
        methylation_fraction = methylated_sum / total_observations

    return {
        "cpg_sites_called": cpg_sites_called,
        "approx_total_methylated_counts": round(methylated_sum),
        "approx_total_unmethylated_counts": round(unmethylated_sum),
        "mean_coverage_called_cpgs": mean_coverage,
        "median_coverage_called_cpgs": median_coverage,
        "coverage_weighted_methylation_fraction": methylation_fraction,
    }


def run_command(
    command: Sequence[str],
) -> str:
    """Run a command safely and return standard output."""
    completed = subprocess.run(
        list(command),
        check=True,
        text=True,
        capture_output=True,
    )

    if completed.stderr:
        print(completed.stderr, end="")

    return completed.stdout


def get_reference_contigs(
    bam_file: str,
    excluded_contigs: Sequence[str],
) -> List[Tuple[str, int]]:
    """
    Read reference contig names and lengths from the BAM header.

    Spike-in contigs are excluded from the human aligned-depth denominator.
    """
    excluded = set(excluded_contigs)

    header = run_command(
        [
            "samtools",
            "view",
            "-H",
            bam_file,
        ]
    )

    contigs: List[Tuple[str, int]] = []

    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue

        name: Optional[str] = None
        length: Optional[int] = None

        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                length = int(field[3:])

        if name is None or length is None:
            continue

        if name in excluded:
            continue

        contigs.append((name, length))

    if not contigs:
        raise RuntimeError(
            "No reference contigs remained after exclusions"
        )

    return contigs


def write_reference_bed(
    contigs: Sequence[Tuple[str, int]],
    bed_file: str,
) -> None:
    """Write one whole-contig BED interval per included reference contig."""
    with open(bed_file, "w") as out:
        for contig, length in contigs:
            out.write(f"{contig}\t0\t{length}\n")


def calculate_aligned_base_coverage(
    bam_file: str,
    min_mapq: int,
    min_baseq: int,
    excluded_contigs: Sequence[str],
) -> Dict[str, object]:
    """
    Calculate mean aligned base depth across selected reference contigs.

    samtools depth:
        -aa includes zero-depth positions.
        -s suppresses overlapping mate double-counting.
        -q sets minimum base quality.
        -Q sets minimum mapping quality.
        -b restricts calculation to human reference contigs.
    """
    contigs = get_reference_contigs(
        bam_file=bam_file,
        excluded_contigs=excluded_contigs,
    )

    expected_reference_bases = sum(
        length
        for _, length in contigs
    )

    if expected_reference_bases <= 0:
        raise RuntimeError(
            "Reference denominator is zero"
        )

    with tempfile.TemporaryDirectory(
        prefix="coverage_qc_"
    ) as temp_dir:
        bed_file = os.path.join(
            temp_dir,
            "included_reference_contigs.bed",
        )

        write_reference_bed(
            contigs=contigs,
            bed_file=bed_file,
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
            bam_file,
        ]

        print(
            "Running:",
            " ".join(command),
        )

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024,
        )

        if process.stdout is None:
            raise RuntimeError(
                "Could not read samtools depth output"
            )

        observed_reference_bases = 0
        total_depth = 0

        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                continue

            try:
                depth = int(fields[2])
            except ValueError:
                continue

            observed_reference_bases += 1
            total_depth += depth

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
                stderr=stderr_text,
            )

    if observed_reference_bases != expected_reference_bases:
        raise RuntimeError(
            "samtools depth returned an unexpected number of positions: "
            f"expected={expected_reference_bases}, "
            f"observed={observed_reference_bases}"
        )

    mean_coverage = total_depth / expected_reference_bases

    return {
        "mean_aligned_base_coverage": mean_coverage,
        "aligned_reference_size_denominator": expected_reference_bases,
        "included_reference_contigs": len(contigs),
    }


def write_output(
    out_file: str,
    metrics: Dict[str, object],
) -> None:
    """Write the one-row coverage QC table."""
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

    output_directory = os.path.dirname(out_file)

    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    with open(out_file, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerow(metrics)


def main() -> None:
    sample = str(snakemake.wildcards.sample)

    cpg_file = str(snakemake.input.cpg)
    bam_file = str(snakemake.input.bam)
    fastp_json = str(snakemake.input.fastp_json)
    out_file = str(snakemake.output.tsv)

    genome_size = int(snakemake.params.genome_size)
    min_mapq = int(snakemake.params.min_mapq)
    min_baseq = int(snakemake.params.min_baseq)

    excluded_contigs = list(
        snakemake.params.excluded_contigs
    )

    print(f"Sample: {sample}")
    print(f"CpG file: {cpg_file}")
    print(f"BAM file: {bam_file}")
    print(f"fastp JSON: {fastp_json}")
    print(f"Output: {out_file}")

    fastp_metrics = read_fastp_metrics(
        fastp_json=fastp_json,
        genome_size=genome_size,
    )

    cpg_metrics = calculate_cpg_metrics(
        cpg_file=cpg_file,
    )

    aligned_metrics = calculate_aligned_base_coverage(
        bam_file=bam_file,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        excluded_contigs=excluded_contigs,
    )

    metrics: Dict[str, object] = {
        "vendor_genome_size_denominator": genome_size,
        "raw_fastq_coverage": format_float(
            fastp_metrics["raw_fastq_coverage"]
        ),
        "mean_aligned_base_coverage": format_float(
            aligned_metrics["mean_aligned_base_coverage"]
        ),
        "aligned_reference_size_denominator":
            aligned_metrics["aligned_reference_size_denominator"],
        "minimum_mapping_quality": min_mapq,
        "minimum_base_quality": min_baseq,
        "overlapping_mates_counted_once": "true",
        "excluded_depth_contigs": ",".join(excluded_contigs),
        "cpg_sites_called": cpg_metrics["cpg_sites_called"],
        "approx_total_methylated_counts":
            cpg_metrics["approx_total_methylated_counts"],
        "approx_total_unmethylated_counts":
            cpg_metrics["approx_total_unmethylated_counts"],
        "mean_coverage_called_cpgs": format_float(
            cpg_metrics["mean_coverage_called_cpgs"]
        ),
        "median_coverage_called_cpgs": format_float(
            cpg_metrics["median_coverage_called_cpgs"]
        ),
        "coverage_weighted_methylation_fraction": format_float(
            cpg_metrics[
                "coverage_weighted_methylation_fraction"
            ]
        ),
    }

    write_output(
        out_file=out_file,
        metrics=metrics,
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


main()
