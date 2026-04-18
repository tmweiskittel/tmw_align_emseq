rule sample_qc_summary:
    input:
        fastp_json=str(QC_DIR / "fastp" / "{sample}.fastp.json"),
        raw_bam=str(BAM_DIR / "{sample}.aligned.sorted.bam"),
        raw_bai=str(BAM_DIR / "{sample}.aligned.sorted.bam.bai"),
        final_bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        final_bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz"),
        lambda_qc=str(SPIKEIN_DIR / "{sample}.lambda_qc.tsv"),
        coverage_qc=str(COVERAGE_DIR / "{sample}.coverage_qc.tsv")
    output:
        tsv=str(SUMMARY_DIR / "{sample}.qc_summary.tsv")
    log:
        str(LOCAL_PATH / "logs" / "sample_qc_summary" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SUMMARY_DIR} $(dirname {log})

        python - <<'PY' > {log} 2>&1
import os
import json
import csv
import subprocess

sample = "{wildcards.sample}"
fastp_json = "{input.fastp_json}"
raw_bam = "{input.raw_bam}"
final_bam = "{input.final_bam}"
lambda_qc_file = "{input.lambda_qc}"
coverage_qc_file = "{input.coverage_qc}"
out_file = "{output.tsv}"

def run_cmd(cmd):
    return subprocess.check_output(cmd, shell=True, text=True).strip()

# fastp
with open(fastp_json) as fh:
    fastp = json.load(fh)

summary_before = fastp.get("summary", {{}}).get("before_filtering", {{}})
summary_after = fastp.get("summary", {{}}).get("after_filtering", {{}})
duplication = fastp.get("duplication", {{}})
insert_size = fastp.get("insert_size", {{}})
filtering_result = fastp.get("filtering_result", {{}})

# BAM metrics
raw_reads = int(run_cmd(f"samtools view -c {raw_bam}"))
final_reads = int(run_cmd(f"samtools view -c {final_bam}"))
raw_mapped_reads = int(run_cmd(f"samtools view -c -F 4 {raw_bam}"))
final_mapped_reads = int(run_cmd(f"samtools view -c -F 4 {final_bam}"))

percent_reads_retained = "NA"
if raw_reads > 0:
    percent_reads_retained = f"{final_reads / raw_reads:.6f}"

percent_mapped_reads_retained = "NA"
if raw_mapped_reads > 0:
    percent_mapped_reads_retained = f"{final_mapped_reads / raw_mapped_reads:.6f}"

raw_bam_size_bytes = os.path.getsize(raw_bam)
final_bam_size_bytes = os.path.getsize(final_bam)

idxstats_output = run_cmd(f"samtools idxstats {final_bam}")
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
        lambda_reads = mapped
    elif chrom == "pUC19":
        puc19_reads = mapped
    elif chrom != "*":
        human_reads += mapped

# lambda QC
lambda_metrics = {{}}
with open(lambda_qc_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        lambda_metrics = row
        break

# coverage QC
coverage_metrics = {{}}
with open(coverage_qc_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        coverage_metrics = row
        break

fieldnames = [
    "sample",

    "fastp_before_total_reads",
    "fastp_before_total_bases",
    "fastp_before_q20_bases",
    "fastp_before_q30_bases",
    "fastp_before_q20_rate",
    "fastp_before_q30_rate",
    "fastp_before_read1_mean_length",
    "fastp_before_read2_mean_length",
    "fastp_before_gc_content",

    "fastp_after_total_reads",
    "fastp_after_total_bases",
    "fastp_after_q20_bases",
    "fastp_after_q30_bases",
    "fastp_after_q20_rate",
    "fastp_after_q30_rate",
    "fastp_after_read1_mean_length",
    "fastp_after_read2_mean_length",
    "fastp_after_gc_content",

    "fastp_passed_filter_reads",
    "fastp_low_quality_reads",
    "fastp_too_many_N_reads",
    "fastp_too_short_reads",
    "fastp_too_long_reads",
    "fastp_duplication_rate",
    "fastp_peak_insert_size",

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

    "cpg_sites_called",
    "total_methylated_counts",
    "total_unmethylated_counts",
    "mean_coverage",
    "median_coverage",
    "mean_methylation_fraction",
    "cpg_sites_ge_1x",
    "cpg_sites_ge_5x",
    "cpg_sites_ge_10x",
    "cpg_sites_ge_20x",

    "lambda_cpg_sites",
    "lambda_methylated_counts",
    "lambda_unmethylated_counts",
    "lambda_mean_methylation_fraction"
]

row = {{
    "sample": sample,

    "fastp_before_total_reads": summary_before.get("total_reads", "NA"),
    "fastp_before_total_bases": summary_before.get("total_bases", "NA"),
    "fastp_before_q20_bases": summary_before.get("q20_bases", "NA"),
    "fastp_before_q30_bases": summary_before.get("q30_bases", "NA"),
    "fastp_before_q20_rate": summary_before.get("q20_rate", "NA"),
    "fastp_before_q30_rate": summary_before.get("q30_rate", "NA"),
    "fastp_before_read1_mean_length": summary_before.get("read1_mean_length", "NA"),
    "fastp_before_read2_mean_length": summary_before.get("read2_mean_length", "NA"),
    "fastp_before_gc_content": summary_before.get("gc_content", "NA"),

    "fastp_after_total_reads": summary_after.get("total_reads", "NA"),
    "fastp_after_total_bases": summary_after.get("total_bases", "NA"),
    "fastp_after_q20_bases": summary_after.get("q20_bases", "NA"),
    "fastp_after_q30_bases": summary_after.get("q30_bases", "NA"),
    "fastp_after_q20_rate": summary_after.get("q20_rate", "NA"),
    "fastp_after_q30_rate": summary_after.get("q30_rate", "NA"),
    "fastp_after_read1_mean_length": summary_after.get("read1_mean_length", "NA"),
    "fastp_after_read2_mean_length": summary_after.get("read2_mean_length", "NA"),
    "fastp_after_gc_content": summary_after.get("gc_content", "NA"),

    "fastp_passed_filter_reads": filtering_result.get("passed_filter_reads", "NA"),
    "fastp_low_quality_reads": filtering_result.get("low_quality_reads", "NA"),
    "fastp_too_many_N_reads": filtering_result.get("too_many_N_reads", "NA"),
    "fastp_too_short_reads": filtering_result.get("too_short_reads", "NA"),
    "fastp_too_long_reads": filtering_result.get("too_long_reads", "NA"),
    "fastp_duplication_rate": duplication.get("rate", "NA"),
    "fastp_peak_insert_size": insert_size.get("peak", "NA"),

    "raw_aligned_reads": raw_reads,
    "raw_mapped_reads": raw_mapped_reads,
    "final_blacklist_filtered_reads": final_reads,
    "final_blacklist_filtered_mapped_reads": final_mapped_reads,
    "percent_reads_retained": percent_reads_retained,
    "percent_mapped_reads_retained": percent_mapped_reads_retained,
    "raw_bam_size_bytes": raw_bam_size_bytes,
    "final_bam_size_bytes": final_bam_size_bytes,

    "final_human_mapped_reads": human_reads,
    "final_lambda_mapped_reads": lambda_reads,
    "final_pUC19_mapped_reads": puc19_reads,

    "cpg_sites_called": coverage_metrics.get("cpg_sites_called", "NA"),
    "total_methylated_counts": coverage_metrics.get("total_methylated_counts", "NA"),
    "total_unmethylated_counts": coverage_metrics.get("total_unmethylated_counts", "NA"),
    "mean_coverage": coverage_metrics.get("mean_coverage", "NA"),
    "median_coverage": coverage_metrics.get("median_coverage", "NA"),
    "mean_methylation_fraction": coverage_metrics.get("mean_methylation_fraction", "NA"),
    "cpg_sites_ge_1x": coverage_metrics.get("cpg_sites_ge_1x", "NA"),
    "cpg_sites_ge_5x": coverage_metrics.get("cpg_sites_ge_5x", "NA"),
    "cpg_sites_ge_10x": coverage_metrics.get("cpg_sites_ge_10x", "NA"),
    "cpg_sites_ge_20x": coverage_metrics.get("cpg_sites_ge_20x", "NA"),

    "lambda_cpg_sites": lambda_metrics.get("lambda_cpg_sites", "NA"),
    "lambda_methylated_counts": lambda_metrics.get("lambda_methylated_counts", "NA"),
    "lambda_unmethylated_counts": lambda_metrics.get("lambda_unmethylated_counts", "NA"),
    "lambda_mean_methylation_fraction": lambda_metrics.get("lambda_mean_methylation_fraction", "NA"),
}}

with open(out_file, "w", newline="") as out_fh:
    writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
PY
        """
