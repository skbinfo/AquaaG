#!/usr/bin/env bash
set -euo pipefail

# ============================================
# AquaG Setup Script (portable / user-agnostic)
# ============================================

# --- Project root (directory where this script lives) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"

# --- Environment Names (must match environment_*.yml 'name' fields if set) ---
MAIN_ENV_NAME="assembly_tool_env"
BRAKER_ENV_NAME="braker_env"       # helper env for Docker CLI etc.

DOCKER_IMAGE="${DOCKER_IMAGE:-teambraker/braker3}"
DOCKER_TAG="${DOCKER_TAG:-latest}" # can be overridden by env var

echo "Project root detected as: ${PROJECT_ROOT}"

# --- 0. Check Conda/Mamba availability ---
echo "Checking for Conda..."
if ! command -v conda &> /dev/null; then
    echo "Error: Conda not found. Install Miniconda or Anaconda first."
    exit 1
fi

# Load conda functions
# shellcheck source=/dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"

if command -v mamba &> /dev/null; then
    echo "Mamba found! Using Mamba for faster environment creation. 🚀"
    SOLVER="mamba"
else
    echo "Mamba not found. Falling back to Conda. ⚠️"
    SOLVER="conda"
fi

# --- 1. Create Main Conda Environment ---
echo "Creating Main Conda environment '${MAIN_ENV_NAME}' (if not already present)..."

if conda env list | awk '{print $1}' | grep -qx "${MAIN_ENV_NAME}"; then
    echo "Environment '${MAIN_ENV_NAME}' already exists. Skipping creation."
else
    if [[ ! -f "${PROJECT_ROOT}/environment_main.yml" ]]; then
        echo "FATAL: environment_main.yml not found in ${PROJECT_ROOT}"
        exit 1
    fi
    # If the YAML already has a 'name:' field, Conda will use that.
    # If not, you can add '-n ${MAIN_ENV_NAME}' here.
    ${SOLVER} env create -f "${PROJECT_ROOT}/environment_main.yml"
    echo "✅ '${MAIN_ENV_NAME}' created successfully."
fi

# --- 2. Create helper env for Docker / CLI tools ---
echo "Creating helper env '${BRAKER_ENV_NAME}' (if not already present)..."

if conda env list | awk '{print $1}' | grep -qx "${BRAKER_ENV_NAME}"; then
    echo "Environment '${BRAKER_ENV_NAME}' already exists. Skipping creation."
else
    if [[ ! -f "${PROJECT_ROOT}/environment_braker.yml" ]]; then
        echo "FATAL: environment_braker.yml not found in ${PROJECT_ROOT}"
        exit 1
    fi
    ${SOLVER} env create -f "${PROJECT_ROOT}/environment_braker.yml"
    echo "✅ '${BRAKER_ENV_NAME}' created successfully."
fi

# Activate helper env
conda activate "${BRAKER_ENV_NAME}" || true

# --- 3. Ensure Docker CLI is available ---
if ! command -v docker &>/dev/null; then
    echo "Docker CLI not found in PATH."
    echo "Attempting to install docker-cli and docker-compose from conda-forge (client only)..."
    conda install -y -c conda-forge docker-cli docker-compose || true
fi

if ! command -v docker &>/dev/null; then
    echo "⚠️  Docker CLI still not found."
    echo "    Please install Docker engine (rootless or system-wide) yourself."
    echo "    See Docker docs: https://docs.docker.com/engine/install/"
    echo "    The pipeline will work once 'docker run hello-world' succeeds."
else
    echo "Docker CLI detected: $(docker --version || true)"
fi

# --- 4. Pull BRAKER image (optional but recommended) ---
if command -v docker &>/dev/null; then
    echo "Pulling BRAKER image ${DOCKER_IMAGE}:${DOCKER_TAG} ..."
    docker pull "${DOCKER_IMAGE}:${DOCKER_TAG}"
    echo "✅ Docker image pulled."
fi

# --- 5. AUGUSTUS config (force replace, always ours) ---

DEFAULT_AUG_DIR="${PROJECT_ROOT}/augustus_config"

# 1) Decide which directory to use on the host
if [[ -n "${AUGUSTUS_CONFIG_PATH:-}" ]]; then
    echo "Existing AUGUSTUS_CONFIG_PATH detected: ${AUGUSTUS_CONFIG_PATH}"
    read -r -p "Use this path as the target for AquaG AUGUSTUS config (it will be overwritten)? [Y/n] " ans
    ans=${ans:-Y}
    if [[ "${ans}" =~ ^[Yy]$ ]]; then
        AUGUSTUS_CONFIG_DIR="${AUGUSTUS_CONFIG_PATH}"
    else
        read -r -p "Enter directory for AUGUSTUS config [${DEFAULT_AUG_DIR}]: " user_aug
        AUGUSTUS_CONFIG_DIR="${user_aug:-$DEFAULT_AUG_DIR}"
    fi
else
    read -r -p "Enter directory for AUGUSTUS config [${DEFAULT_AUG_DIR}]: " user_aug
    AUGUSTUS_CONFIG_DIR="${user_aug:-$DEFAULT_AUG_DIR}"
fi

echo "AUGUSTUS config will be installed to: ${AUGUSTUS_CONFIG_DIR}"
echo "⚠️  Any existing contents in this directory will be removed."
read -r -p "Proceed and overwrite this directory? [Y/n] " overwrite_ans
overwrite_ans=${overwrite_ans:-Y}
if [[ ! "${overwrite_ans}" =~ ^[Yy]$ ]]; then
    echo "Aborting setup at user's request."
    exit 1
fi

# Remove any existing config dir and recreate
if [[ -d "${AUGUSTUS_CONFIG_DIR}" ]]; then
    rm -rf "${AUGUSTUS_CONFIG_DIR}"
fi
mkdir -p "${AUGUSTUS_CONFIG_DIR}"

# 2) Copy config from the BRAKER Docker image
if command -v docker &>/dev/null; then
    echo "Detecting AUGUSTUS config path inside BRAKER container..."

    # Ask container for its AUGUSTUS_CONFIG_PATH
    CONTAINER_AUG_CFG=$(docker run --rm "${DOCKER_IMAGE}:${DOCKER_TAG}" \
        bash -lc 'echo "${AUGUSTUS_CONFIG_PATH:-}"')

    if [[ -z "${CONTAINER_AUG_CFG}" ]]; then
        echo "⚠️  Container AUGUSTUS_CONFIG_PATH is not set."
        echo "    Trying common locations inside the container..."
        CONTAINER_AUG_CFG=$(docker run --rm "${DOCKER_IMAGE}:${DOCKER_TAG}" \
            bash -lc 'for d in /usr/share/augustus/config /usr/local/share/augustus/config /opt/augustus/config; do if [ -d "$d" ]; then echo "$d"; break; fi; done')
    fi

    if [[ -z "${CONTAINER_AUG_CFG}" ]]; then
        echo "FATAL: Could not locate AUGUSTUS config inside the container."
        echo "       Please check the image or copy AUGUSTUS config manually."
        exit 1
    fi

    echo "Using container AUGUSTUS config from: ${CONTAINER_AUG_CFG}"
    echo "Copying to host directory: ${AUGUSTUS_CONFIG_DIR}"

    docker run --rm "${DOCKER_IMAGE}:${DOCKER_TAG}" \
        bash -lc 'tar -C "'"${CONTAINER_AUG_CFG}"'" -czf - .' \
      | tar -C "${AUGUSTUS_CONFIG_DIR}" -xzf -

    echo "✅ AUGUSTUS config installed to ${AUGUSTUS_CONFIG_DIR}"
else
    echo "FATAL: Docker is not available, cannot copy AUGUSTUS config from BRAKER image."
    echo "       Install Docker or manually copy AUGUSTUS config to ${AUGUSTUS_CONFIG_DIR}."
    exit 1
fi

# 3) Export AUGUSTUS_CONFIG_PATH for future shells
if ! grep -q "AUGUSTUS_CONFIG_PATH" "${HOME}/.bashrc" 2>/dev/null; then
    {
      echo ""
      echo "# AUGUSTUS configuration path for BRAKER / AquaG"
      echo "export AUGUSTUS_CONFIG_PATH=\"${AUGUSTUS_CONFIG_DIR}\""
    } >> "${HOME}/.bashrc"
fi

echo "AUGUSTUS_CONFIG_PATH set to: ${AUGUSTUS_CONFIG_DIR}"

# --- 6. GeneMark-ETP license key (gm_key) ---
echo "--------------------------------------------------------------------"
echo "ACTION REQUIRED: GeneMark-ETP license key (gm_key)"
echo "Get gm_key from: https://exon.gatech.edu/GeneMark/license_download.cgi"
echo "We do NOT install GeneMark-ETP on the host; Docker image has /opt/ETP."
echo "--------------------------------------------------------------------"

# Default location: project-local .braker_keys (not in ~ by default)
DEFAULT_GM_KEY_TARGET="${PROJECT_ROOT}/.braker_keys/gm_key"
mkdir -p "$(dirname "${DEFAULT_GM_KEY_TARGET}")"

while true; do
    read -r -p "Enter FULL path to your 'gm_key' (uncompressed): " GM_KEY_PATH
    if [[ -f "${GM_KEY_PATH}" ]]; then
        cp "${GM_KEY_PATH}" "${DEFAULT_GM_KEY_TARGET}"
        chmod 600 "${DEFAULT_GM_KEY_TARGET}"
        echo "✅ License key copied to ${DEFAULT_GM_KEY_TARGET}"
        break
    else
        echo "ERROR: File not found: ${GM_KEY_PATH}"
    fi
done

# Export BRAKER_GM_KEY in .bashrc so Python can find it
if ! grep -q "BRAKER_GM_KEY" "${HOME}/.bashrc" 2>/dev/null; then
    {
      echo ""
      echo "# GeneMark-ETP license key for BRAKER (used by AquaG)"
      echo "export BRAKER_GM_KEY=\"${DEFAULT_GM_KEY_TARGET}\""
    } >> "${HOME}/.bashrc"
fi

conda deactivate || true

echo "--------------------------------------------------------------------"
echo "✅ Setup complete!"
echo "• Project root : ${PROJECT_ROOT}"
echo "• Main env     : ${MAIN_ENV_NAME}  (QUAST/Prokka/BUSCO/etc.)"
echo "• Helper env   : ${BRAKER_ENV_NAME} (Docker client, etc.)"
echo "• AUGUSTUS cfg : ${AUGUSTUS_CONFIG_DIR}"
echo "• gm_key path  : ${DEFAULT_GM_KEY_TARGET} (also in \$BRAKER_GM_KEY)"
echo ""
echo "IMPORTANT:"
echo "1) Open a new shell (or 'source ~/.bashrc') so AUGUSTUS_CONFIG_PATH and BRAKER_GM_KEY are available."
echo "2) Make sure 'docker run hello-world' succeeds (engine running, user in 'docker' group)."
echo "3) BRAKER container will use its own /opt/ETP plus your gm_key (mounted via BRAKER_GM_KEY)."
echo "--------------------------------------------------------------------"

