#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/home/jupyter"
REPO_URL="https://github.com/tmweiskittel/tmwalign.git"
REPO_DIR="${WORKDIR}/tmwalign"
MAMBA_ROOT_PREFIX="/opt/micromamba"
JUPYTER_BASHRC="/home/jupyter/.bashrc"

echo "Starting VM setup..."

mkdir -p "${WORKDIR}"

echo "Installing system packages..."
sudo apt-get update
sudo apt-get install -y git curl wget bzip2 ca-certificates

echo "Installing micromamba if needed..."
if [ ! -x /usr/local/bin/micromamba ]; then
  cd /tmp
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
  sudo mv bin/micromamba /usr/local/bin/micromamba
  rm -rf bin
fi

echo "Configuring micromamba for future shells..."
if ! grep -qxF 'export MAMBA_ROOT_PREFIX=/opt/micromamba' "${JUPYTER_BASHRC}" 2>/dev/null; then
  echo 'export MAMBA_ROOT_PREFIX=/opt/micromamba' >> "${JUPYTER_BASHRC}"
fi

if ! grep -qxF 'eval "$(micromamba shell hook -s bash)"' "${JUPYTER_BASHRC}" 2>/dev/null; then
  echo 'eval "$(micromamba shell hook -s bash)"' >> "${JUPYTER_BASHRC}"
fi

echo "Cloning or updating repo..."
if [ ! -d "${REPO_DIR}/.git" ]; then
  git clone "${REPO_URL}" "${REPO_DIR}"
else
  cd "${REPO_DIR}"
  git pull
fi

echo "Initializing micromamba in this shell..."
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX}"
eval "$(micromamba shell hook -s bash)"

echo "Creating snakemake environment if needed..."
if ! micromamba env list | awk '{print $1}' | grep -qx "snakemake"; then
  micromamba create -y -n snakemake -c conda-forge -c bioconda snakemake
fi

echo
echo "VM setup complete."
echo "Next steps:"
echo "  cd ${REPO_DIR}"
echo "  export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX}"
echo '  eval "$(micromamba shell hook -s bash)"'
echo "  micromamba activate snakemake"
echo "  snakemake --help"
