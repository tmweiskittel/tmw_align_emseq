#!/usr/bin/env bash
set -euo pipefail

WORKDIR="$HOME"
REPO_URL="https://github.com/tmweiskittel/tmwalign.git"
REPO_DIR="${WORKDIR}/tmwalign"
USER_BASHRC="${HOME}/.bashrc"

echo "Starting VM setup..."

mkdir -p "${WORKDIR}"

echo "Installing system packages..."
sudo apt-get update
sudo apt-get install -y git curl wget bzip2 ca-certificates

echo "Checking that conda is available..."
if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda was not found in PATH."
  echo "On this VM, conda may live at /opt/conda/bin/conda."
  exit 1
fi

echo "Configuring conda for future shells..."
if ! grep -qxF 'source /opt/conda/etc/profile.d/conda.sh' "${USER_BASHRC}" 2>/dev/null; then
  echo 'source /opt/conda/etc/profile.d/conda.sh' >> "${USER_BASHRC}"
fi

echo "Cloning or updating repo..."
if [ ! -d "${REPO_DIR}/.git" ]; then
  git clone "${REPO_URL}" "${REPO_DIR}"
else
  cd "${REPO_DIR}"
  git pull
fi

echo "Initializing conda in this shell..."
source /opt/conda/etc/profile.d/conda.sh

echo "Leaving active base environment if needed..."
conda deactivate >/dev/null 2>&1 || true

echo "Creating snakemake environment if needed..."
if ! conda env list | awk '{print $1}' | grep -qx "snakemake"; then
  conda create -y \
    -n snakemake \
    -c conda-forge \
    -c bioconda \
    --override-channels \
    snakemake-minimal
fi

echo
echo "VM setup complete."
echo "Next steps:"
echo "  cd ${REPO_DIR}"
echo "  source /opt/conda/etc/profile.d/conda.sh"
echo "  conda activate snakemake"
echo "  snakemake --help"
