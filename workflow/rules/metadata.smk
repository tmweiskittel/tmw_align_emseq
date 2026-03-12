from pathlib import Path

SAMPLESHEET = Path(config['meta']['data_paths']) / config['meta']['sample_sheet']
MD5_FILE = Path(config['meta']['data_paths']) / config['meta']['md5_sheet']
LOCAL_PATH= Path(config['meta']['local_paths'])
rule download_metadata:
    output:
        samples=LOCAL_PATH / samples.csv,
        md5=LOCAL_PATH / md5.txt
    params:
        samples=SAMPLESHEET,
        md5=MD5_FILE
    log:
        "logs/download_metadata.log"
    shell:
        r"""
        mkdir -p resources/metadata logs
        mkdir -p /home/jupyter/data
        gcloud storage cp {params.samples_uri} {output.samples}
        gcloud storage cp {params.md5_uri} {output.md5}
        """
