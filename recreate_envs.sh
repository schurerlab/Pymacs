#!/usr/bin/env bash
# ===================================================
# 🧬 Conda Environment Restoration Script
# Recreates cgenff and mdanalysis environments
# Installs required OS build tools for pip-compiled packages like DockQ
# ===================================================

set -euo pipefail

# -------------------------------
# Helpers
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CGENFF_YML="environment_cgenff.yml"
MDA_YML="environment_mdanalysis.yml"

require_file() {
  local file="$1"
  local label="$2"
  if [[ ! -f "$file" ]]; then
    echo "❌ Missing: $label ($file)"
    exit 1
  fi
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

env_exists() {
  local env_name="$1"
  conda env list | awk '{print $1}' | grep -Fxq "$env_name"
}

get_env_name() {
  local yml_file="$1"
  awk '/^name:/ {print $2; exit}' "$yml_file"
}

# -------------------------------
# Validate files
# -------------------------------
require_file "$CGENFF_YML" "cgenff environment file"
require_file "$MDA_YML" "mdanalysis environment file"

# -------------------------------
# Validate conda
# -------------------------------
if ! have_cmd conda; then
  echo "❌ Conda is not installed or not on PATH."
  echo "Install Miniconda first, then re-run this script."
  exit 1
fi

# Ensure conda shell functions are available in non-interactive shells
eval "$(conda shell.bash hook)"

# -------------------------------
# Install system prerequisites
# -------------------------------
echo "🔍 Checking required system tools..."

NEED_BUILD_TOOLS=0
for tool in gcc g++ make; do
  if ! have_cmd "$tool"; then
    NEED_BUILD_TOOLS=1
    break
  fi
done

if [[ "$NEED_BUILD_TOOLS" -eq 1 ]]; then
  echo "➡️  Missing compiler toolchain. Installing build-essential..."
  sudo apt update
  sudo apt install -y build-essential
else
  echo "✅ build-essential already available"
fi

# Optional but often useful for scientific Python builds
if ! have_cmd pkg-config; then
  echo "➡️  Installing pkg-config..."
  sudo apt install -y pkg-config
else
  echo "✅ pkg-config already available"
fi

# -------------------------------
# Read env names from YAML
# -------------------------------
CGENFF_ENV_NAME="$(get_env_name "$CGENFF_YML")"
MDA_ENV_NAME="$(get_env_name "$MDA_YML")"

if [[ -z "${CGENFF_ENV_NAME:-}" ]]; then
  echo "❌ Could not read environment name from $CGENFF_YML"
  exit 1
fi

if [[ -z "${MDA_ENV_NAME:-}" ]]; then
  echo "❌ Could not read environment name from $MDA_YML"
  exit 1
fi

echo "📦 Creating Conda environments..."

# -------------------------------
# Create cgenff environment
# -------------------------------
echo "➡️  Processing environment: $CGENFF_ENV_NAME"
if env_exists "$CGENFF_ENV_NAME"; then
  echo "⚠️  Environment '$CGENFF_ENV_NAME' already exists — skipping..."
else
  conda env create -f "$CGENFF_YML"
fi

# -------------------------------
# Create mdanalysis environment
# -------------------------------
echo "➡️  Processing environment: $MDA_ENV_NAME"
if env_exists "$MDA_ENV_NAME"; then
  echo "⚠️  Environment '$MDA_ENV_NAME' already exists — skipping..."
else
  conda env create -f "$MDA_YML"
fi

# -------------------------------
# Verify key package(s)
# -------------------------------
echo "🧪 Verifying key tools in '$MDA_ENV_NAME'..."

conda activate "$MDA_ENV_NAME"

# Check DockQ explicitly since it previously failed to build without gcc
if python -c "import DockQ" >/dev/null 2>&1; then
  echo "✅ DockQ Python module imports successfully"
else
  echo "⚠️  DockQ import failed. Attempting repair via pip..."
  pip install --no-cache-dir dockq==2.1.3
  python -c "import DockQ" >/dev/null 2>&1 && echo "✅ DockQ repaired successfully"
fi

# Optional CLI check
if command -v dockq >/dev/null 2>&1; then
  echo "✅ dockq CLI available"
else
  echo "⚠️  dockq CLI not found on PATH inside '$MDA_ENV_NAME'"
fi

conda deactivate

echo
echo "✅ All environments ready!"
echo "Use:"
echo "  conda activate $CGENFF_ENV_NAME"
echo "  conda activate $MDA_ENV_NAME"
echo
echo "💡 Recommendation: run this project from your Linux home, e.g. ~/projects/Pymacs"
echo "   instead of /mnt/c/... for better performance and fewer permission issues."