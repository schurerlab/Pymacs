#!/bin/bash
# ===================================================
# 🧬 Conda Environment Restoration Script
# Recreates cgenff and mdanalysis environments
# ===================================================

set -euo pipefail

# Check for environment files
if [[ ! -f "cgenff_environment.yml" ]]; then
  echo "❌ Missing: cgenff_environment.yml"
  exit 1
fi

if [[ ! -f "mdanalysis_environment.yml" ]]; then
  echo "❌ Missing: mdanalysis_environment.yml"
  exit 1
fi

echo "📦 Creating Conda environments..."

# Create cgenff environment
echo "➡️  Creating environment: cgenff"
conda env create -f cgenff_environment.yml || {
  echo "⚠️  Environment 'cgenff' already exists — skipping..."
}

# Create mdanalysis environment
echo "➡️  Creating environment: mdanalysis"
conda env create -f mdanalysis_environment.yml || {
  echo "⚠️  Environment 'mdanalysis' already exists — skipping..."
}

echo
echo "✅ All environments ready!"
echo "Use 'conda activate cgenff' or 'conda activate mdanalysis' to start."
