from pathlib import Path

SAMPLESHEET = Path(DATA_PATH) / SAMPLE_SHEET
MD5_FILE = Path(DATA_PATH) / MD5_SHEET

rule download_metadata:
    output:
        samples=str(LOCAL_SAMPLES),
        md5=str(LOCAL_MD5)
    params:
        samples=str(SAMPLESHEET),
        md5=str(MD5_FILE)
    log:
        "logs/download_metadata.log"
    shell:
        r"""
        mkdir -p {LOCAL_PATH} logs
        gcloud storage cp gs://{params.samples} {output.samples} > {log} 2>&1
        gcloud storage cp gs://{params.md5} {output.md5} >> {log} 2>&1
        """
