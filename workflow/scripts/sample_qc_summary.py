#!/usr/bin/env python3

import csv
import json
import os
import subprocess
from typing import Dict, Sequence


sample = str(snakemake.wildcards.sample)

fastp_json = str(snakemake.input.fastp_json)
raw_bam = str(snakemake.input.raw_bam)
final_bam = str(snakemake.input.final_bam)
lambda_qc_file = str(snakemake.input.lambda_qc)
coverage_qc_file = str(snakemake.input.coverage_qc)

out_file = str(snakemake.output.tsv)


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

    return completed.stdout.strip()


def samtools_count(
    bam_file: str,
    extra_arguments: Sequence[str] = (),
) -> int:
    """Count BAM records with optional samtools view filters."""
    command = [
        "samtools",
        "view",
        "-c",
        *extra_arguments,
        bam_file,
    ]

    return int(run_command(command))


def read_first_tsv_row(filename: str) -> Dict[str, str]:
    """Read the first data row from a tab-delimited file."""
    with open(filename, newline="") as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t",
        )

        for row in reader:
            return dict(row)

    return {}


with open(fastp_json) as fh:
    fastp = json.load(fh)


summary_before = (
    fastp
    .get("summary", {})
    .get("before_filtering", {})
)

summary_after = (
    fastp
    .get("summary", {})
    .get("after_filtering", {})
)

duplication = fastp.get(
    "duplication",
    {},
)

insert_size = fastp.get(
    "insert_size",
    {},
)

filtering_result = fastp.get(
    "filtering_result",
    {},
)


raw_reads = samtools_count(
    raw_bam,
)

final_reads = samtools_count(
    final_bam,
)

raw_mapped_reads = samtools_count(
    raw_bam,
    extra_arguments=[
        "-F",
        "4",
    ],
)

final_mapped_reads = samtools_count(
    final_bam,
    extra_arguments=[
        "-F",
        "4",
    ],
)


percent_reads_retained = "NA"

if raw_reads > 0:
    percent_reads_retained = (
        f"{final_reads / raw_reads:.6f}"
    )


percent_mapped_reads_retained = "NA"

if raw_mapped_reads > 0:
    percent_mapped_reads_retained = (
        f"{final_mapped_reads / raw_mapped_reads:.6f}"
    )


raw_bam_size_bytes = os.path.getsize(
    raw_bam
)

final_bam_size_bytes = os.path.getsize(
    final_bam
)


idxstats_output = run_command(
    [
        "samtools",
        "idxstats",
        final_bam,
    ]
)

lambda_reads = 0
puc19_reads = 0
human_reads = 0

for line in idxstats_output.splitlines():
    fields = line.split("\t")

    if len(fields) < 4:
        continue

    chrom = fields[0]
    mapped = int(fields[2])

    if chrom == "lambda":
        lambda_reads += mapped
    elif chrom == "pUC19":
        puc19_reads += mapped
    elif chrom != "*":
        human_reads += mapped


lambda_metrics = read_first_tsv_row(
    lambda_qc_file
)

coverage_metrics = read_first_tsv_row(
    coverage_qc_file
)


fieldnames = [
    "sample",

    # fastp before filtering
    "fastp_before_total_reads",
    "fastp_before_total_bases",
    "fastp_before_q20_bases",
    "fastp_before_q30_bases",
    "fastp_before_q20_rate",
    "fastp_before_q30_rate",
    "fastp_before_read1_mean_length",
    "fastp_before_read2_mean_length",
    "fastp_before_gc_content",

    # fastp after filtering
    "fastp_after_total_reads",
    "fastp_after_total_bases",
    "fastp_after_q20_bases",
    "fastp_after_q30_bases",
    "fastp_after_q20_rate",
    "fastp_after_q30_rate",
    "fastp_after_read1_mean_length",
    "fastp_after_read2_mean_length",
    "fastp_after_gc_content",

    # fastp filtering and duplication
    "fastp_passed_filter_reads",
    "fastp_low_quality_reads",
    "fastp_too_many_N_reads",
    "fastp_too_short_reads",
    "fastp_too_long_reads",
    "fastp_duplication_rate",
    "fastp_peak_insert_size",

    # BAM retention and composition
    "raw_aligned_reads",
    "raw_mapped_reads",
    "final_blacklist_filtered_reads",
    "final_blacklist_filtered_mapped_reads",
    "percent_reads_retained",
    "percent_mapped_reads_retained",
    "raw_bam_size_bytes",
    "final_bam_size_bytes",
    "final_human_mapped_reads",
    "final_lambda_mapped_reads",
    "final_pUC19_mapped_reads",

    # Coverage QC
    "vendor_genome_size_denominator",
    "raw_fastq_coverage",
    "mean_aligned_base_coverage",
    "aligned_reference_size_denominator",
    "minimum_mapping_quality",
    "minimum_base_quality",
    "overlapping_mates_counted_once",
    "excluded_depth_contigs",

    # CpG QC
    "cpg_sites_called",
    "approx_total_methylated_counts",
    "approx_total_unmethylated_counts",
    "mean_coverage_called_cpgs",
    "median_coverage_called_cpgs",
    "coverage_weighted_methylation_fraction",

    # Lambda conversion-control QC
    "lambda_cpg_sites",
    "lambda_methylated_counts",
    "lambda_unmethylated_counts",
    "lambda_mean_methylation_fraction",
]


row = {
    "sample": sample,

    # fastp before filtering
    "fastp_before_total_reads":
        summary_before.get(
            "total_reads",
            "NA",
        ),
    "fastp_before_total_bases":
        summary_before.get(
            "total_bases",
            "NA",
        ),
    "fastp_before_q20_bases":
        summary_before.get(
            "q20_bases",
            "NA",
        ),
    "fastp_before_q30_bases":
        summary_before.get(
            "q30_bases",
            "NA",
        ),
    "fastp_before_q20_rate":
        summary_before.get(
            "q20_rate",
            "NA",
        ),
    "fastp_before_q30_rate":
        summary_before.get(
            "q30_rate",
            "NA",
        ),
    "fastp_before_read1_mean_length":
        summary_before.get(
            "read1_mean_length",
            "NA",
        ),
    "fastp_before_read2_mean_length":
        summary_before.get(
            "read2_mean_length",
            "NA",
        ),
    "fastp_before_gc_content":
        summary_before.get(
            "gc_content",
            "NA",
        ),

    # fastp after filtering
    "fastp_after_total_reads":
        summary_after.get(
            "total_reads",
            "NA",
        ),
    "fastp_after_total_bases":
        summary_after.get(
            "total_bases",
            "NA",
        ),
    "fastp_after_q20_bases":
        summary_after.get(
            "q20_bases",
            "NA",
        ),
    "fastp_after_q30_bases":
        summary_after.get(
            "q30_bases",
            "NA",
        ),
    "fastp_after_q20_rate":
        summary_after.get(
            "q20_rate",
            "NA",
        ),
    "fastp_after_q30_rate":
        summary_after.get(
            "q30_rate",
            "NA",
        ),
    "fastp_after_read1_mean_length":
        summary_after.get(
            "read1_mean_length",
            "NA",
        ),
    "fastp_after_read2_mean_length":
        summary_after.get(
            "read2_mean_length",
            "NA",
        ),
    "fastp_after_gc_content":
        summary_after.get(
            "gc_content",
            "NA",
        ),

    # fastp filtering and duplication
    "fastp_passed_filter_reads":
        filtering_result.get(
            "passed_filter_reads",
            "NA",
        ),
    "fastp_low_quality_reads":
        filtering_result.get(
            "low_quality_reads",
            "NA",
        ),
    "fastp_too_many_N_reads":
        filtering_result.get(
            "too_many_N_reads",
            "NA",
        ),
    "fastp_too_short_reads":
        filtering_result.get(
            "too_short_reads",
            "NA",
        ),
    "fastp_too_long_reads":
        filtering_result.get(
            "too_long_reads",
            "NA",
        ),
    "fastp_duplication_rate":
        duplication.get(
            "rate",
            "NA",
        ),
    "fastp_peak_insert_size":
        insert_size.get(
            "peak",
            "NA",
        ),

    # BAM retention and composition
    "raw_aligned_reads": raw_reads,
    "raw_mapped_reads": raw_mapped_reads,
    "final_blacklist_filtered_reads": final_reads,
    "final_blacklist_filtered_mapped_reads":
        final_mapped_reads,
    "percent_reads_retained":
        percent_reads_retained,
    "percent_mapped_reads_retained":
        percent_mapped_reads_retained,
    "raw_bam_size_bytes":
        raw_bam_size_bytes,
    "final_bam_size_bytes":
        final_bam_size_bytes,
    "final_human_mapped_reads":
        human_reads,
    "final_lambda_mapped_reads":
        lambda_reads,
    "final_pUC19_mapped_reads":
        puc19_reads,

    # Coverage QC
    "vendor_genome_size_denominator":
        coverage_metrics.get(
            "vendor_genome_size_denominator",
            "NA",
        ),
    "raw_fastq_coverage":
        coverage_metrics.get(
            "raw_fastq_coverage",
            "NA",
        ),
    "mean_aligned_base_coverage":
        coverage_metrics.get(
            "mean_aligned_base_coverage",
            "NA",
        ),
    "aligned_reference_size_denominator":
        coverage_metrics.get(
            "aligned_reference_size_denominator",
            "NA",
        ),
    "minimum_mapping_quality":
        coverage_metrics.get(
            "minimum_mapping_quality",
            "NA",
        ),
    "minimum_base_quality":
        coverage_metrics.get(
            "minimum_base_quality",
            "NA",
        ),
    "overlapping_mates_counted_once":
        coverage_metrics.get(
            "overlapping_mates_counted_once",
            "NA",
        ),
    "excluded_depth_contigs":
        coverage_metrics.get(
            "excluded_depth_contigs",
            "NA",
        ),

    # CpG QC
    "cpg_sites_called":
        coverage_metrics.get(
            "cpg_sites_called",
            "NA",
        ),
    "approx_total_methylated_counts":
        coverage_metrics.get(
            "approx_total_methylated_counts",
            "NA",
        ),
    "approx_total_unmethylated_counts":
        coverage_metrics.get(
            "approx_total_unmethylated_counts",
            "NA",
        ),
    "mean_coverage_called_cpgs":
        coverage_metrics.get(
            "mean_coverage_called_cpgs",
            "NA",
        ),
    "median_coverage_called_cpgs":
        coverage_metrics.get(
            "median_coverage_called_cpgs",
            "NA",
        ),
    "coverage_weighted_methylation_fraction":
        coverage_metrics.get(
            "coverage_weighted_methylation_fraction",
            "NA",
        ),

    # Lambda QC
    "lambda_cpg_sites":
        lambda_metrics.get(
            "lambda_cpg_sites",
            "NA",
        ),
    "lambda_methylated_counts":
        lambda_metrics.get(
            "lambda_methylated_counts",
            "NA",
        ),
    "lambda_unmethylated_counts":
        lambda_metrics.get(
            "lambda_unmethylated_counts",
            "NA",
        ),
    "lambda_mean_methylation_fraction":
        lambda_metrics.get(
            "lambda_mean_methylation_fraction",
            "NA",
        ),
}


output_directory = os.path.dirname(
    out_file
)

if output_directory:
    os.makedirs(
        output_directory,
        exist_ok=True,
    )


with open(
    out_file,
    "w",
    newline="",
) as out_fh:
    writer = csv.DictWriter(
        out_fh,
        fieldnames=fieldnames,
        delimiter="\t",
        extrasaction="raise",
    )

    writer.writeheader()
    writer.writerow(row)
