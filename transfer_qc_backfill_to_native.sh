#!/usr/bin/env bash

set -euo pipefail

SAMPLE_SHEET="${SAMPLE_SHEET:-/home/jupyter/data/samples.csv}"
BACKFILL_ROOT="${BACKFILL_ROOT:-gs://weiskittel-projects1/radnecrosis/qc_backfill}"
CLOUD_SEARCH_ROOT="${CLOUD_SEARCH_ROOT:-gs://weiskittel-projects1/radnecrosis}"
WORK_DIR="${WORK_DIR:-/home/jupyter/data/promote_qc_backfill}"

INVENTORY="${WORK_DIR}/cloud_inventory.txt"
SUCCESS_LOG="${WORK_DIR}/transfer_successes.tsv"
FAILURE_LOG="${WORK_DIR}/transfer_failures.tsv"

mkdir -p "${WORK_DIR}"

printf "sample\tcoverage_destination\tsummary_destination\n" > "${SUCCESS_LOG}"
printf "sample\treason\n" > "${FAILURE_LOG}"

for program in gcloud python3 awk; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: ${program}" >&2
        exit 1
    fi
done

if [[ ! -s "${SAMPLE_SHEET}" ]]; then
    echo "ERROR: Sample sheet is missing or empty: ${SAMPLE_SHEET}" >&2
    exit 1
fi

echo "Inventorying native cloud outputs..."

gcloud storage ls \
    --recursive \
    "${CLOUD_SEARCH_ROOT}/**" \
    > "${INVENTORY}"

find_native_object() {
    local basename="$1"
    local matches
    local count

    matches=$(
        awk \
            -v target="/${basename}" \
            -v backfill="${BACKFILL_ROOT}/" '
            index($0, backfill) == 1 {
                next
            }

            length($0) >= length(target) &&
            substr($0, length($0) - length(target) + 1) == target {
                print $0
            }
        ' "${INVENTORY}"
    )

    count=$(
        printf "%s\n" "${matches}" |
            awk 'NF {n++} END {print n + 0}'
    )

    if [[ "${count}" -eq 0 ]]; then
        echo "No native object found for ${basename}" >&2
        return 1
    fi

    if [[ "${count}" -gt 1 ]]; then
        echo "Multiple native objects found for ${basename}:" >&2
        printf "%s\n" "${matches}" >&2
        return 2
    fi

    printf "%s\n" "${matches}"
}

mapfile -t SAMPLES < <(
    python3 - "${SAMPLE_SHEET}" <<'PY'
import csv
import sys

with open(
    sys.argv[1],
    newline="",
    encoding="utf-8-sig",
) as fh:
    reader = csv.DictReader(fh)

    if reader.fieldnames is None:
        raise SystemExit("Sample sheet has no header")

    normalized = {
        field.strip(): field
        for field in reader.fieldnames
        if field is not None
    }

    if "Name" not in normalized:
        raise SystemExit(
            "Sample sheet does not contain a Name column"
        )

    name_field = normalized["Name"]
    seen = set()

    for row in reader:
        sample = (row.get(name_field) or "").strip()

        if not sample or sample in seen:
            continue

        seen.add(sample)
        print(sample)
PY
)

for sample in "${SAMPLES[@]}"; do
    echo "Processing ${sample}"

    source_coverage="${BACKFILL_ROOT}/coverage/${sample}.coverage_qc.tsv"
    source_summary="${BACKFILL_ROOT}/summary/${sample}.qc_summary.tsv"

    if ! native_coverage=$(
        find_native_object "${sample}.coverage_qc.tsv"
    ); then
        printf "%s\t%s\n" \
            "${sample}" \
            "Native coverage destination missing or ambiguous" \
            >> "${FAILURE_LOG}"
        continue
    fi

    if ! native_summary=$(
        find_native_object "${sample}.qc_summary.tsv"
    ); then
        printf "%s\t%s\n" \
            "${sample}" \
            "Native summary destination missing or ambiguous" \
            >> "${FAILURE_LOG}"
        continue
    fi

    if ! gcloud storage cp \
        "${source_coverage}" \
        "${native_coverage}"; then

        printf "%s\t%s\n" \
            "${sample}" \
            "Coverage transfer failed" \
            >> "${FAILURE_LOG}"
        continue
    fi

    if ! gcloud storage cp \
        "${source_summary}" \
        "${native_summary}"; then

        printf "%s\t%s\n" \
            "${sample}" \
            "Summary transfer failed" \
            >> "${FAILURE_LOG}"
        continue
    fi

    printf "%s\t%s\t%s\n" \
        "${sample}" \
        "${native_coverage}" \
        "${native_summary}" \
        >> "${SUCCESS_LOG}"

    echo "Completed ${sample}"
done

echo
echo "Transfer complete."
echo "Successes: ${SUCCESS_LOG}"
echo "Failures: ${FAILURE_LOG}"
