#!/usr/bin/env bash

set -euo pipefail


# ============================================================================
# User configuration
# ============================================================================

# Sample sheet used by the original project.
SAMPLE_SHEET="${SAMPLE_SHEET:-/home/jupyter/data/samples.csv}"

# Root beneath which the existing per-sample pipeline outputs are stored.
# The script inventories this location and finds files by exact basename.
CLOUD_SEARCH_ROOT="gs://weiskittel-projects1/radnecrosis"

# Use a new destination initially so the original QC files are not overwritten.
CLOUD_OUTPUT_ROOT="gs://weiskittel-projects1/radnecrosis/qc_backfill"

# Local repository root. This should contain workflow/scripts/.
REPO_PATH="${REPO_PATH:-$PWD}"

# Existing local fastp JSON directory.
FASTP_JSON_DIR="${FASTP_JSON_DIR:-/home/jupyter/data/qc/fastp}"

# Local temporary working directory.
#
# Place this on the disk used by the original pipeline for large temporary
# files. Only one sample is retained here at a time.
WORK_ROOT="${WORK_ROOT:-/home/jupyter/data/qc_backfill_tmp}"

# Existing Python scripts.
COVERAGE_SCRIPT="${REPO_PATH}/workflow/scripts/coverage_qc.py"
SUMMARY_SCRIPT="${REPO_PATH}/workflow/scripts/sample_qc_summary_once.py"

# Coverage settings.
GENOME_SIZE=3100000000
MIN_MAPQ=30
MIN_BASEQ=0
EXCLUDED_CONTIGS="lambda,pUC19"


# ============================================================================
# Local paths
# ============================================================================

CLOUD_INVENTORY="${WORK_ROOT}/cloud_inventory.txt"
FAILURE_LOG="${WORK_ROOT}/qc_backfill_failures.tsv"
SUCCESS_LOG="${WORK_ROOT}/qc_backfill_successes.tsv"

mkdir -p "${WORK_ROOT}"

printf "sample\tstage\treason\n" > "${FAILURE_LOG}"
printf "sample\tcoverage_output\tsummary_output\n" > "${SUCCESS_LOG}"


# ============================================================================
# Validation
# ============================================================================
if [[ ! -s "${SAMPLE_SHEET}" ]]; then
     echo "ERROR: Sample sheet is missing or empty:" >&2
     echo "  ${SAMPLE_SHEET}" >&2
      exit 1
fi
    
for program in gcloud python3 samtools awk; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: ${program}" >&2
        exit 1
    fi
done

if [[ ! -f "${COVERAGE_SCRIPT}" ]]; then
    echo "ERROR: Coverage script not found:" >&2
    echo "  ${COVERAGE_SCRIPT}" >&2
    exit 1
fi

if [[ ! -f "${SUMMARY_SCRIPT}" ]]; then
    echo "ERROR: Summary script not found:" >&2
    echo "  ${SUMMARY_SCRIPT}" >&2
    exit 1
fi

if [[ ! -d "${FASTP_JSON_DIR}" ]]; then
    echo "ERROR: fastp JSON directory not found:" >&2
    echo "  ${FASTP_JSON_DIR}" >&2
    exit 1
fi


# ============================================================================
# Helper functions
# ============================================================================

record_failure() {
    local sample="$1"
    local stage="$2"
    local reason="$3"

    printf "%s\t%s\t%s\n" \
        "${sample}" \
        "${stage}" \
        "${reason}" \
        >> "${FAILURE_LOG}"
}


# Find exactly one cloud object whose basename matches the requested filename.
#
# This avoids assuming a specific uploaded directory structure. It also refuses
# to proceed if duplicate objects with the same basename exist.
find_unique_cloud_object() {
    local requested_basename="$1"
    local matches
    local count

    matches=$(
        awk -v target="/${requested_basename}" '
            length($0) >= length(target) &&
            substr($0, length($0) - length(target) + 1) == target {
                print $0
            }
        ' "${CLOUD_INVENTORY}"
    )

    count=$(
        printf "%s\n" "${matches}" |
            awk 'NF > 0 {n++} END {print n + 0}'
    )

    if [[ "${count}" -eq 0 ]]; then
        echo "No cloud object found for ${requested_basename}" >&2
        return 1
    fi

    if [[ "${count}" -gt 1 ]]; then
        echo "Multiple cloud objects found for ${requested_basename}:" >&2
        printf "%s\n" "${matches}" >&2
        return 2
    fi

    printf "%s\n" "${matches}"
}


download_required_object() {
    local requested_basename="$1"
    local local_destination="$2"
    local cloud_source

    cloud_source=$(
        find_unique_cloud_object "${requested_basename}"
    )

    echo "Downloading:"
    echo "  ${cloud_source}"
    echo "  -> ${local_destination}"

    gcloud storage cp \
        "${cloud_source}" \
        "${local_destination}"
}


upload_file() {
    local local_source="$1"
    local cloud_destination="$2"

    echo "Uploading:"
    echo "  ${local_source}"
    echo "  -> ${cloud_destination}"

    gcloud storage cp \
        "${local_source}" \
        "${cloud_destination}"
}


cleanup_sample_directory() {
    local sample_directory="$1"

    if [[ -n "${sample_directory}" && -d "${sample_directory}" ]]; then
        rm -rf "${sample_directory}"
    fi
}



# ============================================================================
# Inventory existing project outputs
# ============================================================================

echo "Creating cloud-object inventory beneath:"
echo "  ${CLOUD_SEARCH_ROOT}"

# The ** glob is intentional and is quoted so the local shell does not expand it.
gcloud storage ls \
    --recursive \
    "${CLOUD_SEARCH_ROOT}/**" \
    > "${CLOUD_INVENTORY}"

if [[ ! -s "${CLOUD_INVENTORY}" ]]; then
    echo "ERROR: Cloud inventory is empty." >&2
    exit 1
fi


# ============================================================================
# Read sample names
# ============================================================================

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

    normalized_names = {
        name.strip(): name
        for name in reader.fieldnames
        if name is not None
    }

    if "Name" not in normalized_names:
        raise SystemExit(
            "Sample sheet does not contain a Name column. "
            f"Columns found: {reader.fieldnames}"
        )

    name_column = normalized_names["Name"]
    seen = set()

    for row in reader:
        sample = (row.get(name_column) or "").strip()

        if not sample or sample in seen:
            continue

        seen.add(sample)
        print(sample)
PY
)

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: No samples were found in the sample sheet." >&2
    exit 1
fi

echo "Samples found: ${#SAMPLES[@]}"


# ============================================================================
# Process each sample sequentially
# ============================================================================

for sample in "${SAMPLES[@]}"; do
    echo
    echo "========================================================================"
    echo "Processing ${sample}"
    echo "========================================================================"

    sample_dir="${WORK_ROOT}/${sample}"

    # Remove a partial directory from an earlier failed attempt.
    cleanup_sample_directory "${sample_dir}"
    mkdir -p "${sample_dir}"

    fastp_json="${FASTP_JSON_DIR}/${sample}.fastp.json"

    raw_bam="${sample_dir}/${sample}.aligned.sorted.bam"

    final_bam="${sample_dir}/${sample}.aligned.sorted.filt.bl.bam"
    final_bai="${sample_dir}/${sample}.aligned.sorted.filt.bl.bam.bai"

    cpg_file="${sample_dir}/${sample}.CpG.methylKit.gz"
    lambda_qc="${sample_dir}/${sample}.lambda_qc.tsv"

    coverage_qc="${sample_dir}/${sample}.coverage_qc.tsv"
    summary_qc="${sample_dir}/${sample}.qc_summary.tsv"

    coverage_log="${sample_dir}/${sample}.coverage_qc.log"
    summary_log="${sample_dir}/${sample}.qc_summary.log"

    cloud_coverage_output=(
        "${CLOUD_OUTPUT_ROOT}/coverage/"
        "${sample}.coverage_qc.tsv"
    )
    cloud_coverage_output="${cloud_coverage_output[*]}"
    cloud_coverage_output="${cloud_coverage_output// /}"

    cloud_summary_output=(
        "${CLOUD_OUTPUT_ROOT}/summary/"
        "${sample}.qc_summary.tsv"
    )
    cloud_summary_output="${cloud_summary_output[*]}"
    cloud_summary_output="${cloud_summary_output// /}"

    cloud_coverage_log=(
        "${CLOUD_OUTPUT_ROOT}/logs/"
        "${sample}.coverage_qc.log"
    )
    cloud_coverage_log="${cloud_coverage_log[*]}"
    cloud_coverage_log="${cloud_coverage_log// /}"

    cloud_summary_log=(
        "${CLOUD_OUTPUT_ROOT}/logs/"
        "${sample}.qc_summary.log"
    )
    cloud_summary_log="${cloud_summary_log[*]}"
    cloud_summary_log="${cloud_summary_log// /}"


    # ------------------------------------------------------------------------
    # Validate the existing local fastp JSON before downloading large BAMs
    # ------------------------------------------------------------------------

    
    if [[ ! -s "${fastp_json}" ]]; then
        echo "ERROR: fastp JSON is missing or empty:"
        echo "  ${fastp_json}"

        record_failure \
            "${sample}" \
            "input_validation" \
            "Missing or empty fastp JSON"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if ! python3 - "${fastp_json}" <<'PY'
import json
import sys

filename = sys.argv[1]

with open(filename, encoding="utf-8") as fh:
    report = json.load(fh)

before = report.get("summary", {}).get("before_filtering", {})

if before.get("total_bases") is None:
    raise SystemExit(
        "fastp JSON is missing summary.before_filtering.total_bases"
    )
PY
    then
        record_failure \
            "${sample}" \
            "input_validation" \
            "fastp JSON failed validation"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi


    # ------------------------------------------------------------------------
    # Download required pipeline artifacts
    # ------------------------------------------------------------------------

    if ! download_required_object \
        "${sample}.aligned.sorted.bam" \
        "${raw_bam}"; then

        record_failure \
            "${sample}" \
            "download" \
            "Raw aligned BAM missing or ambiguous"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if ! download_required_object \
        "${sample}.aligned.sorted.filt.bl.bam" \
        "${final_bam}"; then

        record_failure \
            "${sample}" \
            "download" \
            "Final filtered BAM missing or ambiguous"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if ! download_required_object \
        "${sample}.aligned.sorted.filt.bl.bam.bai" \
        "${final_bai}"; then

        record_failure \
            "${sample}" \
            "download" \
            "Final filtered BAM index missing or ambiguous"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if ! download_required_object \
        "${sample}.CpG.methylKit.gz" \
        "${cpg_file}"; then

        record_failure \
            "${sample}" \
            "download" \
            "CpG methylKit file missing or ambiguous"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if ! download_required_object \
        "${sample}.lambda_qc.tsv" \
        "${lambda_qc}"; then

        record_failure \
            "${sample}" \
            "download" \
            "Lambda QC file missing or ambiguous"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi


    # ------------------------------------------------------------------------
    # Validate BAMs
    # ------------------------------------------------------------------------

    if ! samtools quickcheck \
        "${raw_bam}" \
        "${final_bam}"; then

        record_failure \
            "${sample}" \
            "bam_validation" \
            "samtools quickcheck failed"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    # This confirms that the final BAM index is found and usable.
    if ! samtools idxstats \
        "${final_bam}" \
        >/dev/null; then

        record_failure \
            "${sample}" \
            "bam_validation" \
            "Final BAM index is missing or unusable"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi


    # ------------------------------------------------------------------------
    # Generate corrected coverage QC
    # ------------------------------------------------------------------------

    echo "Running coverage QC for ${sample}..."

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

        echo "ERROR: coverage_qc.py failed for ${sample}."

        record_failure \
            "${sample}" \
            "coverage_qc" \
            "coverage_qc.py failed"

        if [[ -s "${coverage_log}" ]]; then
            upload_file \
                "${coverage_log}" \
                "${CLOUD_OUTPUT_ROOT}/logs/${sample}.coverage_qc.failed.log" \
                || true
        fi

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if [[ ! -s "${coverage_qc}" ]]; then
        record_failure \
            "${sample}" \
            "coverage_qc" \
            "Coverage script produced no output"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi


    # ------------------------------------------------------------------------
    # Reconstruct the complete QC summary
    # ------------------------------------------------------------------------

    echo "Running QC summary reconstruction for ${sample}..."

    if ! python3 "${SUMMARY_SCRIPT}" \
        --sample "${sample}" \
        --fastp-json "${fastp_json}" \
        --raw-bam "${raw_bam}" \
        --final-bam "${final_bam}" \
        --lambda-qc "${lambda_qc}" \
        --coverage-qc "${coverage_qc}" \
        --output "${summary_qc}" \
        > "${summary_log}" 2>&1; then

        echo "ERROR: sample_qc_summary_once.py failed for ${sample}."

        record_failure \
            "${sample}" \
            "qc_summary" \
            "sample_qc_summary_once.py failed"

        if [[ -s "${coverage_log}" ]]; then
            upload_file \
                "${coverage_log}" \
                "${CLOUD_OUTPUT_ROOT}/logs/${sample}.coverage_qc.log" \
                || true
        fi

        if [[ -s "${summary_log}" ]]; then
            upload_file \
                "${summary_log}" \
                "${CLOUD_OUTPUT_ROOT}/logs/${sample}.qc_summary.failed.log" \
                || true
        fi

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if [[ ! -s "${summary_qc}" ]]; then
        record_failure \
            "${sample}" \
            "qc_summary" \
            "Summary script produced no output"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi


    # ------------------------------------------------------------------------
    # Basic output validation
    # ------------------------------------------------------------------------

    if ! python3 - \
        "${coverage_qc}" \
        "${summary_qc}" \
        "${sample}" <<'PY'
import csv
import sys

coverage_file = sys.argv[1]
summary_file = sys.argv[2]
expected_sample = sys.argv[3]


def read_one_row(filename):
    with open(filename, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            raise RuntimeError(f"No header in {filename}")

        rows = list(reader)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one data row in {filename}; found {len(rows)}"
        )

    return reader.fieldnames, rows[0]


coverage_fields, coverage_row = read_one_row(coverage_file)
summary_fields, summary_row = read_one_row(summary_file)

required_coverage_fields = {
    "raw_fastq_coverage",
    "mean_aligned_base_coverage",
    "cpg_sites_called",
    "mean_coverage_called_cpgs",
    "median_coverage_called_cpgs",
}

missing_coverage = required_coverage_fields - set(coverage_fields)

if missing_coverage:
    raise RuntimeError(
        "Coverage output is missing fields: "
        + ", ".join(sorted(missing_coverage))
    )

required_summary_fields = {
    "sample",
    "fastp_before_total_bases",
    "raw_fastq_coverage",
    "mean_aligned_base_coverage",
    "mean_coverage_called_cpgs",
    "lambda_mean_methylation_fraction",
}

missing_summary = required_summary_fields - set(summary_fields)

if missing_summary:
    raise RuntimeError(
        "Summary output is missing fields: "
        + ", ".join(sorted(missing_summary))
    )

if summary_row.get("sample") != expected_sample:
    raise RuntimeError(
        f"Summary sample mismatch: expected {expected_sample!r}, "
        f"found {summary_row.get('sample')!r}"
    )
PY
    then
        record_failure \
            "${sample}" \
            "output_validation" \
            "Generated TSV failed validation"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi


    # ------------------------------------------------------------------------
    # Upload outputs and logs
    # ------------------------------------------------------------------------

    if ! upload_file \
        "${coverage_qc}" \
        "${cloud_coverage_output}"; then

        record_failure \
            "${sample}" \
            "upload" \
            "Coverage QC upload failed"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    if ! upload_file \
        "${summary_qc}" \
        "${cloud_summary_output}"; then

        record_failure \
            "${sample}" \
            "upload" \
            "QC summary upload failed"

        cleanup_sample_directory "${sample_dir}"
        continue
    fi

    upload_file \
        "${coverage_log}" \
        "${cloud_coverage_log}"

    upload_file \
        "${summary_log}" \
        "${cloud_summary_log}"

    printf "%s\t%s\t%s\n" \
        "${sample}" \
        "${cloud_coverage_output}" \
        "${cloud_summary_output}" \
        >> "${SUCCESS_LOG}"

    # Removes both large BAMs before the next sample starts.
    cleanup_sample_directory "${sample_dir}"

    echo "Completed ${sample}"
done


# ============================================================================
# Upload run-level reports
# ============================================================================

upload_file \
    "${FAILURE_LOG}" \
    "${CLOUD_OUTPUT_ROOT}/qc_backfill_failures.tsv"

upload_file \
    "${SUCCESS_LOG}" \
    "${CLOUD_OUTPUT_ROOT}/qc_backfill_successes.tsv"

echo
echo "QC backfill finished."
echo
echo "Output root:"
echo "  ${CLOUD_OUTPUT_ROOT}"
echo
echo "Local failure report:"
echo "  ${FAILURE_LOG}"
echo
echo "Local success report:"
echo "  ${SUCCESS_LOG}"
