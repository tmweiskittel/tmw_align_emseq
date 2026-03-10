#!/usr/bin/env bash
set -euo pipefail

WORKDIR=/home/jupyter
REPO_URL=https://github.com/tmweiskittel/tmwalign.git
REPO_DIR=${WORKDIR}/pipeline

mkdir -p "${WORKDIR}"

# install basics
apt-get update
apt-get install -y git curl wget bzip2

# install micromamba if not already present
if [ ! -x /usr/local/bin/micromamba ]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
  mv bin/micromamba /usr/local/bin/micromamba
  rm -rf bin
fi

# clone repo if needed
if [ ! -d "${REPO_DIR}" ]; then
  git clone "${REPO_URL}" "${REPO_DIR}"
else
  cd "${REPO_DIR}"
  git pull
fi

# create base snakemake env if needed
export MAMBA_ROOT_PREFIX=/opt/micromamba
eval "$(micromamba shell hook -s bash)"

if ! micromamba env list | grep -q '^snakemake '; then
  micromamba create -y -n snakemake -c conda-forge -c bioconda snakemake
fi

echo "VM setup complete."
echo "To start working:"
echo "  cd ${REPO_DIR}"
echo "  micromamba activate snakemake"
echo "  snakemake --profile profiles/gcp"
