#!/usr/bin/env bash

set -euo pipefail

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

SAMPLE_SHEET_GS="gs://weiskittel-projects1/radnecrosis/raw_data/samples.csv"

# Search all existing project objects beneath this location.
CLOUD_SEARCH_ROOT="gs://weiskittel-projects1/radnecrosis"

# New QC files will be written here.
CLOUD_OUTPUT_ROOT="gs://weiskittel-projects1/radnecrosis/qc_backfill"

# Repository containing workflow/scripts.
REPO_PATH="${REPO_PATH:-$PWD}"

COVERAGE_SCRIPT="${REPO_PATH}/workflow/scripts/coverage_qc.py"
SUMMARY_SCRIPT="${REPO_PATH}/workflow/scripts/sample_qc_summary_once.py"

WORK_ROOT="${WORK_ROOT:-/tmp/radnecrosis_qc_backfill}"

GENOME_SIZE=3100000000
MIN_MAPQ=30
MIN_BASEQ=0
EXCLUDED_CONTIGS="lambda,pUC19"

SAMPLE_SHEET="${WORK_ROOT}/samples.csv"
CLOUD_INVENTORY="${WORK_ROOT}/cloud_inventory.txt"
FAILURE_LOG="${WORK_ROOT}/failed_samples.tsv"

mkdir -p "${WORK_ROOT}"

printf "sample\treason\n" > "${FAILURE_LOG}"

# ----------------------------------------------------------------------
# Validate programs and scripts
# ----------------------------------------------------------------------

for program in gcloud python3 samtools awk grep sort; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "Required program not found: ${program}" >&2
        exit 1
    fi
done

if [[ ! -f "${COVERAGE_SCRIPT}" ]]; then
    echo "Coverage script not found: ${COVERAGE_SCRIPT}" >&2
    exit 1
fi

if [[ ! -f "${SUMMARY_SCRIPT}" ]]; then
    echo "Summary script not found: ${SUMMARY_SCRIPT}" >&2
    exit 1
fi

# ----------------------------------------------------------------------
# Download sample sheet and create cloud inventory
# ----------------------------------------------------------------------

echo "Downloading sample sheet..."
gcloud storage cp \
    "${SAMPLE_SHEET_GS}" \
    "${SAMPLE_SHEET}"

echo "Inventorying cloud files..."
gcloud storage ls \
    --recursive \
    "${CLOUD_SEARCH_ROOT}/**" \
    > "${CLOUD_INVENTORY}"

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

find_unique_object() {
    local basename="$1"
    local matches
    local count

    matches=$(
        awk -v target="/${basename}" '
            substr($0, length($0) - length(target) + 1) == target {
                print
            }
        ' "${CLOUD_INVENTORY}"
    )

    count=$(
        printf "%s\n" "${matches}" \
            | awk 'NF > 0 {count++} END {print count + 0}'
    )

    if [[ "${count}" -eq 0 ]]; then
        echo "No cloud object found with basename: ${basename}" >&2
        return 1
    fi

    if [[ "${count}" -gt 1 ]]; then
        echo "Multiple cloud objects found with basename: ${basename}" >&2
        printf "%s\n" "${matches}" >&2
        return 2
    fi

    printf "%s\n" "${matches}"
}


download_required_file() {
    local basename="$1"
    local destination="$2"
    local cloud_object

    cloud_object=$(find_unique_object "${basename}")

    echo "Downloading ${cloud_object}"
    gcloud storage cp \
        "${cloud_object}" \
        "${destination}"
}


cleanup_sample() {
    local sample_dir="$1"

    rm -rf "${sample_dir}"
}


record_failure() {
    local sample="$1"
    local reason="$2"

    printf "%s\t%s\n" \
        "${sample}" \
        "${reason}" \
        >> "${FAILURE_LOG}"
}


# ----------------------------------------------------------------------
# Read the Name column from samples.csv
# ----------------------------------------------------------------------

mapfile -t SAMPLES < <(
    python3 - "${SAMPLE_SHEET}" <<'PY'
import csv
import sys

sample_sheet = sys.argv[1]

with open(
    sample_sheet,
    newline="",
    encoding="utf-8-sig",
) as fh:
    reader = csv.DictReader(fh)

    if reader.fieldnames is None:
        raise SystemExit("Sample sheet has no header")

    normalized = {
        name.strip(): name
        for name in reader.fieldnames
        if name is not None
    }

    if "Name" not in normalized:
        raise SystemExit(
            "Sample sheet does not contain a Name column. "
            f"Columns: {reader.fieldnames}"
        )

    source_column = normalized["Name"]

    for row in reader:
        sample = (row.get(source_column) or "").strip()

        if sample:
            print(sample)
PY
)

echo "Samples found: ${#SAMPLES[@]}"

# ----------------------------------------------------------------------
# Process samples
# ----------------------------------------------------------------------

for sample in "${SAMPLES[@]}"; do
    echo
    echo "============================================================"
    echo "Processing ${sample}"
    echo "============================================================"

    sample_dir="${WORK_ROOT}/${sample}"
    mkdir -p "${sample_dir}"

    raw_bam="${sample_dir}/${sample}.aligned.sorted.bam"
    final_bam="${sample_dir}/${sample}.aligned.sorted.filt.bl.bam"
    final_bai="${final_bam}.bai"

    cpg_file="${sample_dir}/${sample}.CpG.methylKit.gz"
    fastp_json="${sample_dir}/${sample}.fastp.json"
    lambda_qc="${sample_dir}/${sample}.lambda_qc.tsv"

    coverage_qc="${sample_dir}/${sample}.coverage_qc.tsv"
    summary_qc="${sample_dir}/${sample}.qc_summary.tsv"

    coverage_log="${sample_dir}/${sample}.coverage_qc.log"
    summary_log="${sample_dir}/${sample}.qc_summary.log"

    sample_failed=0

    # Adjust this lambda-QC basename here if the existing files use a
    # different naming convention.
    required_basenames=(
        "${sample}.aligned.sorted.bam"
        "${sample}.aligned.sorted.filt.bl.bam"
        "${sample}.aligned.sorted.filt.bl.bam.bai"
        "${sample}.CpG.methylKit.gz"
        "${sample}.fastp.json"
        "${sample}.lambda_qc.tsv"
    )

    destinations=(
        "${raw_bam}"
        "${final_bam}"
        "${final_bai}"
        "${cpg_file}"
        "${fastp_json}"
        "${lambda_qc}"
    )

    for index in "${!required_basenames[@]}"; do
        if ! download_required_file \
            "${required_basenames[$index]}" \
            "${destinations[$index]}"; then

            record_failure \
                "${sample}" \
                "Missing or ambiguous object: ${required_basenames[$index]}"

            sample_failed=1
            break
        fi
    done

    if [[ "${sample_failed}" -eq 1 ]]; then
        cleanup_sample "${sample_dir}"
        continue
    fi

    # Confirm the downloaded final BAM index is usable.
    if ! samtools quickcheck \
        "${raw_bam}" \
        "${final_bam}"; then

        record_failure \
            "${sample}" \
            "samtools quickcheck failed"

        cleanup_sample "${sample_dir}"
        continue
    fi

    if ! samtools idxstats \
        "${final_bam}" \
        >/dev/null; then

        record_failure \
            "${sample}" \
            "Final BAM index could not be used"

        cleanup_sample "${sample_dir}"
        continue
    fi

    echo "Running coverage QC..."

    if ! python3 "${COVERAGE_SCRIPT}" \
        --sample "${sample}" \
        --cpg "${cpg_file}" \
        --bam "${final_bam}" \
        --fastp-json "${fastp_json}" \
        --output "${coverage_qc}" \
        --genome-size "${GENOME_SIZE}" \
        --min-mapq "${MIN_MAPQ}" \
        --min-baseq "${MIN_BASEQ}" \
        --excluded-contigs "${EXCLUDED_CONTIGS}" \
        > "${coverage_log}" 2>&1; then

        record_failure \
            "${sample}" \
            "coverage_qc.py failed; see ${coverage_log}"

        # Preserve the failed log for inspection.
        gcloud storage cp \
            "${coverage_log}" \
            "${CLOUD_OUTPUT_ROOT}/logs/${sample}.coverage_qc.failed.log" \
            || true

        cleanup_sample "${sample_dir}"
        continue
    fi

    echo "Reconstructing sample QC summary..."

    if ! python3 "${SUMMARY_SCRIPT}" \
        --sample "${sample}" \
        --fastp-json "${fastp_json}" \
        --raw-bam "${raw_bam}" \
        --final-bam "${final_bam}" \
        --lambda-qc "${lambda_qc}" \
        --coverage-qc "${coverage_qc}" \
        --output "${summary_qc}" \
        > "${summary_log}" 2>&1; then

        record_failure \
            "${sample}" \
            "sample_qc_summary_once.py failed; see ${summary_log}"

        gcloud storage cp \
            "${summary_log}" \
            "${CLOUD_OUTPUT_ROOT}/logs/${sample}.qc_summary.failed.log" \
            || true

        cleanup_sample "${sample_dir}"
        continue
    fi

    echo "Uploading reconstructed QC..."

    gcloud storage cp \
        "${coverage_qc}" \
        "${CLOUD_OUTPUT_ROOT}/coverage/${sample}.coverage_qc.tsv"

    gcloud storage cp \
        "${summary_qc}" \
        "${CLOUD_OUTPUT_ROOT}/summary/${sample}.qc_summary.tsv"

    gcloud storage cp \
        "${coverage_log}" \
        "${CLOUD_OUTPUT_ROOT}/logs/${sample}.coverage_qc.log"

    gcloud storage cp \
        "${summary_log}" \
        "${CLOUD_OUTPUT_ROOT}/logs/${sample}.qc_summary.log"

    cleanup_sample "${sample_dir}"

    echo "Completed ${sample}"
done

# Upload the failure report even when it contains only the header.
gcloud storage cp \
    "${FAILURE_LOG}" \
    "${CLOUD_OUTPUT_ROOT}/qc_backfill_failures.tsv"

echo
echo "Backfill complete."
echo "Outputs: ${CLOUD_OUTPUT_ROOT}"
echo "Failure report: ${CLOUD_OUTPUT_ROOT}/qc_backfill_failures.tsv"
