with gzip.open(cpg_file, "rt") as fh:
    for line in fh:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue

        if fields[0] == "chrBase":
            continue

        chrom = fields[1]  # methylKit chrom column, not chrBase

        if chrom != "lambda":
            continue

        try:
            coverage = int(fields[4])
            freqC = float(fields[5])
            freqT = float(fields[6])
        except ValueError:
            continue

        m = coverage * freqC / 100.0
        u = coverage * freqT / 100.0

        sites += 1
        meth += m
        unmeth += u

mean_methyl = "NA"
if meth + unmeth > 0:
    mean_methyl = meth / (meth + unmeth)

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
    "%.6f" % mean_methyl if mean_methyl != "NA" else "NA"
])
PY
        """
