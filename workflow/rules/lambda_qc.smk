rule lambda_spikein_qc:
    input:
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz")
    output:
        tsv=temp(str(SPIKEIN_DIR / "{sample}.lambda_qc.tsv"))
    log:
        str(LOCAL_PATH / "logs" / "lambda_spikein_qc" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SPIKEIN_DIR} $(dirname {log})

        python3 - <<'PY' > {log} 2>&1
import gzip
import csv

cpg_file = "{input.cpg}"
out_file = "{output.tsv}"

sites = 0
meth = 0.0
unmeth = 0.0

with gzip.open(cpg_file, "rt") as fh:
    for line in fh:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue

        if fields[0] == "chrBase":
            continue

        chrom = fields[1]

        if chrom != "lambda":
            continue

        try:
            coverage = int(fields[4])
            freqC = float(fields[5])
            freqT = float(fields[6])
        except ValueError:
            continue

        sites += 1
        meth += coverage * freqC / 100.0
        unmeth += coverage * freqT / 100.0

mean_methyl = "NA"
if meth + unmeth > 0:
    mean_methyl = "%.6f" % (meth / (meth + unmeth))

with open(out_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "lambda_cpg_sites",
        "lambda_methylated_counts",
        "lambda_unmethylated_counts",
        "lambda_mean_methylation_fraction"
    ])
    writer.writerow([
        sites,
        round(meth),
        round(unmeth),
        mean_methyl
    ])
PY
        """
