#!/usr/bin/env python3

import sys
import csv
from pathlib import Path
from collections import defaultdict

pairs = defaultdict(dict)

for line in sys.stdin:
    gcs_path = line.strip()
    if not gcs_path:
        continue

    # remove gs:// prefix
    no_gs = gcs_path.replace("gs://", "")

    path = Path(no_gs)
    filename = path.name
    sample = path.parent.name
    folder = str(path.parent)

    if filename.endswith("_1.fq.gz") or filename.endswith("_R1.fq.gz"):
        pairs[sample]["R1_file"] = filename
        pairs[sample]["Path"] = folder
    elif filename.endswith("_2.fq.gz") or filename.endswith("_R2.fq.gz"):
        pairs[sample]["R2_file"] = filename
        pairs[sample]["Path"] = folder

writer = csv.DictWriter(
    sys.stdout,
    fieldnames=["Name", "R1_file", "R2_file", "Path"]
)

writer.writeheader()

for sample in sorted(pairs):
    row = pairs[sample]
    if "R1_file" in row and "R2_file" in row:
        writer.writerow({
            "Name": sample,
            "R1_file": row["R1_file"],
            "R2_file": row["R2_file"],
            "Path": row["Path"]
        })
