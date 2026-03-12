from pathlib import Path

SAMPLESHEET = Path(config['meta']['data_paths']) / config['meta']['sample_sheet']
MD5_FILE = Path(config['meta']['data_paths']) / config['meta']['md5_sheet']
LOCAL_PATH= Path(config['meta']['local_paths'])
rule download_metadata:
    output:
        samples=str(LOCAL_PATH / "samples.csv"),
        md5=str(LOCAL_PATH / "md5.txt")
    params:
        samples=SAMPLESHEET,
        md5=MD5_FILE
    log:
        "logs/download_metadata.log"
    shell:
        """
        mkdir -p {LOCAL_PATH} logs
        gcloud storage cp {params.samples} {output.samples}
        gcloud storage cp {params.md5} {output.md5}
        """
