#!/bin/bash
# ===================================================
# 🧬 Conda Environment Restoration Script
# Recreates cgenff and mdanalysis environments
# ===================================================

set -euo pipefail

# Check for environment files
if [[ ! -f "environment_cgenff.yml" ]]; then
  echo "❌ Missing: cgenff_environment.yml"
  exit 1
fi

if [[ ! -f "environment_mdanalysis.yml" ]]; then
  echo "❌ Missing: mdanalysis_environment.yml"
  exit 1
fi

echo "📦 Creating Conda environments..."

# Create cgenff environment
echo "➡️  Creating environment: cgenff"
conda env create -f environment_cgenff.yml || {
  echo "⚠️  Environment 'cgenff' already exists — skipping..."
}

# Create mdanalysis environment
echo "➡️  Creating environment: mdanalysis"
conda env create -f environment_mdanalysis.yml || {
  echo "⚠️  Environment 'mdanalysis' already exists — skipping..."
}

echo
echo "✅ All environments ready!"
echo "Use 'conda activate cgenff' or 'conda activate mdanalysis' to start."
