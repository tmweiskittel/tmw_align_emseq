from pathlib import Path

SAMPLESHEET = Path(config['meta']['data_path']) / config['meta']['sample_sheet']
MD5_FILE = Path(config['meta']['data_path']) / config['meta']['md5_sheet']
LOCAL_PATH = Path(config['meta']['local_path'])

LOCAL_SAMPLES = LOCAL_PATH / "samples.csv"
LOCAL_MD5 = LOCAL_PATH / "md5.txt"

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
        """
        mkdir -p {LOCAL_PATH} logs
        gcloud storage cp {params.samples} {output.samples}
        gcloud storage cp {params.md5} {output.md5}
        """
