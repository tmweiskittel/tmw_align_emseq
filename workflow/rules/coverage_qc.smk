rule coverage_qc:
    input:
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz")
    output:
        tsv=str(COVERAGE_DIR / "{sample}.coverage_qc.tsv")
    log:
        str(LOCAL_PATH / "logs" / "coverage_qc" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {COVERAGE_DIR} $(dirname {log})

        python3 - <<'PY' > {log} 2>&1
import gzip
import csv
import statistics

cpg_file = "{input.cpg}"
out_file = "{output.tsv}"

coverages = []
meth_sum = 0
unmeth_sum = 0
site_count = 0

ge_1 = 0
ge_5 = 0
ge_10 = 0
ge_20 = 0

with gzip.open(cpg_file, "rt") as fh:
    for line in fh:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue

        try:
            coverage = int(fields[4])
            meth = int(fields[5])
            unmeth = int(fields[6])
        except ValueError:
            continue

        site_count += 1
        coverages.append(coverage)
        meth_sum += meth
        unmeth_sum += unmeth

        if coverage >= 1:
            ge_1 += 1
        if coverage >= 5:
            ge_5 += 1
        if coverage >= 10:
            ge_10 += 1
        if coverage >= 20:
            ge_20 += 1

mean_coverage = "NA"
median_coverage = "NA"
mean_methylation_fraction = "NA"

if coverages:
    mean_coverage = f"{statistics.mean(coverages):.6f}"
    median_coverage = f"{statistics.median(coverages):.6f}"

if meth_sum + unmeth_sum > 0:
    mean_methylation_fraction = f"{meth_sum / (meth_sum + unmeth_sum):.6f}"

with open(out_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "cpg_sites_called",
        "total_methylated_counts",
        "total_unmethylated_counts",
        "mean_coverage",
        "median_coverage",
        "mean_methylation_fraction",
        "cpg_sites_ge_1x",
        "cpg_sites_ge_5x",
        "cpg_sites_ge_10x",
        "cpg_sites_ge_20x"
    ])
    writer.writerow([
        site_count,
        meth_sum,
        unmeth_sum,
        mean_coverage,
        median_coverage,
        mean_methylation_fraction,
        ge_1,
        ge_5,
        ge_10,
        ge_20
    ])
PY
        """
