SAMPLESHEET = f"{config['paths']['metadata_dir']}/samples.csv"
MD5_FILE = f"{config['paths']['metadata_dir']}/fastq.md5"

rule download_metadata:
    output:
        samples=SAMPLESHEET,
        md5=MD5_FILE
    params:
        samples_uri=config["metadata"]["samples_csv"],
        md5_uri=config["metadata"]["md5_manifest"]
    log:
        "logs/download_metadata.log"
    shell:
        r"""
        mkdir -p resources/metadata logs

        gcloud storage cp {params.samples_uri} {output.samples}
        gcloud storage cp {params.md5_uri} {output.md5}
        """
