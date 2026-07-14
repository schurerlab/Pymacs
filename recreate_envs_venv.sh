#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

RECREATE=false
DRY_RUN=false
SKIP_GROMACS_CHECK=false
PYTHON_BIN=""

VENV_ROOT=".venvs"
CGENFF_VENV="$VENV_ROOT/cgenff"
MDA_VENV="$VENV_ROOT/mdanalysis"
CGENFF_REQ="requirements_cgenff.txt"
MDA_REQ="requirements_mdanalysis.txt"

usage() {
  cat <<'EOF'
Usage: bash recreate_envs_venv.sh [options]

Create or refresh the two non-Conda PyMACS virtual environments:
  .venvs/cgenff
  .venvs/mdanalysis

Options:
  --recreate            Delete and rebuild both virtual environments
  --dry-run             Print planned actions without changing anything
  --python PATH         Use a specific Python interpreter
  --skip-gromacs-check  Skip the post-install GROMACS executable warning check
  -h, --help            Show this help message

Notes:
  - This script targets Linux, macOS, and WSL.
  - Conda/Mamba remains the recommended fully reproducible setup path.
  - GROMACS is not installed by this script and must be available separately.
EOF
}

log() {
  printf '[INFO] %s\n' "$*"
}

warn() {
  printf '[WARN] %s\n' "$*" >&2
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

run_cmd() {
  if [[ "$DRY_RUN" == true ]]; then
    printf '[DRY-RUN] '
    printf '%q ' "$@"
    printf '\n'
    return 0
  fi
  "$@"
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || die "Required file not found: $path"
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

python_version_tuple() {
  local python_bin="$1"
  "$python_bin" - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
}

python_version_supported() {
  local python_bin="$1"
  "$python_bin" - <<'PY' >/dev/null
import sys
major, minor = sys.version_info[:2]
raise SystemExit(0 if (major, minor) >= (3, 9) and (major, minor) <= (3, 12) else 1)
PY
}

detect_python() {
  if [[ -n "$PYTHON_BIN" ]]; then
    if [[ -x "$PYTHON_BIN" ]]; then
      python_version_supported "$PYTHON_BIN" || die "Specified Python must be between 3.9 and 3.12 for this package set: $PYTHON_BIN"
      return
    fi
    if have_cmd "$PYTHON_BIN"; then
      PYTHON_BIN="$(command -v "$PYTHON_BIN")"
      python_version_supported "$PYTHON_BIN" || die "Specified Python must be between 3.9 and 3.12 for this package set: $PYTHON_BIN"
      return
    fi
    die "Specified Python is not executable: $PYTHON_BIN"
    return
  fi

  local candidate
  local resolved
  local checked=()

  for candidate in \
    python3.12 python3.11 python3.10 python3.9 python3 python \
    /usr/bin/python3 /opt/homebrew/bin/python3 /usr/local/bin/python3; do
    if [[ "$candidate" == /* ]]; then
      [[ -x "$candidate" ]] || continue
      resolved="$candidate"
    else
      have_cmd "$candidate" || continue
      resolved="$(command -v "$candidate")"
    fi

    checked+=("$resolved")
    if python_version_supported "$resolved"; then
      PYTHON_BIN="$resolved"
      return
    fi
  done

  if [[ "${#checked[@]}" -eq 0 ]]; then
    die "No Python interpreter found. Install Python 3.9-3.12 and retry."
  fi

  die "No suitable Python interpreter found. Install Python 3.9-3.12 and retry. Checked: ${checked[*]}"
}

check_python_version() {
  local version_output
  version_output="$("$PYTHON_BIN" - <<'PY'
import sys
major, minor = sys.version_info[:2]
if (major, minor) < (3, 9) or (major, minor) > (3, 12):
    raise SystemExit(f"{major}.{minor}")
print(f"{sys.executable}|{major}.{minor}.{sys.version_info.micro}")
PY
)" || die "PyMACS venv setup requires Python 3.9-3.12 because the selected scientific stack does not install reliably on newer Python versions yet."
  log "Using Python: ${version_output#*|} (${version_output%%|*})"
}

create_venv() {
  local venv_path="$1"
  local selected_version
  selected_version="$(python_version_tuple "$PYTHON_BIN")"

  if [[ -d "$venv_path" && "$RECREATE" == true ]]; then
    log "Removing existing virtual environment: $venv_path"
    run_cmd rm -rf "$venv_path"
  fi

  if [[ -d "$venv_path" && -x "$venv_path/bin/python" ]]; then
    local existing_version
    existing_version="$(python_version_tuple "$venv_path/bin/python")"
    if [[ "$existing_version" != "$selected_version" ]]; then
      warn "Existing virtual environment at $venv_path uses Python $existing_version, but the selected interpreter is Python $selected_version. Rebuilding it."
      run_cmd rm -rf "$venv_path"
    fi
  fi

  if [[ -d "$venv_path" ]]; then
    log "Reusing existing virtual environment: $venv_path"
  else
    log "Creating virtual environment: $venv_path"
    run_cmd "$PYTHON_BIN" -m venv "$venv_path"
  fi
}

upgrade_packaging_tools() {
  local venv_python="$1"
  log "Upgrading pip/setuptools/wheel in ${venv_python%/bin/python}"
  run_cmd "$venv_python" -m pip install --upgrade pip setuptools wheel
}

install_requirements() {
  local venv_python="$1"
  local req_file="$2"
  local env_label="$3"

  log "Installing $env_label requirements from $req_file"
  if [[ "$DRY_RUN" == true ]]; then
    run_cmd "$venv_python" -m pip install -r "$req_file"
    return 0
  fi

  if ! "$venv_python" -m pip install -r "$req_file"; then
    cat >&2 <<EOF
[ERROR] pip install failed for $env_label using $req_file
Helpful next steps:
  1. Review the error above for missing compiler tools or system headers.
  2. Confirm your Python version is 3.9+ and supported by the selected wheels.
  3. If pip resolution or compiled wheels fail on this machine, use the Conda/Mamba path:
     conda env create -f environment_cgenff.yml
     conda env create -f environment_mdanalysis.yml
EOF
    exit 1
  fi
}

install_package() {
  local venv_python="$1"
  local env_label="$2"
  shift 2

  log "Installing extra package(s) for $env_label: $*"
  if [[ "$DRY_RUN" == true ]]; then
    run_cmd "$venv_python" -m pip install "$@"
    return 0
  fi

  if ! "$venv_python" -m pip install "$@"; then
    cat >&2 <<EOF
[ERROR] Extra pip install failed for $env_label: $*
Helpful next steps:
  1. Review the error above for missing compiler tools or system headers.
  2. Retry with a newer Python interpreter via --python if this platform has one.
  3. Fall back to the Conda/Mamba environment files for the most reproducible setup.
EOF
    exit 1
  fi
}

run_smoke_test() {
  local venv_python="$1"
  local env_label="$2"
  local smoke_script="$3"

  log "Running smoke test for $env_label"
  if [[ "$DRY_RUN" == true ]]; then
    printf '[DRY-RUN] %q - <<'\''PY'\''\n' "$venv_python"
    printf '%s\n' "$smoke_script"
    printf 'PY\n'
    return 0
  fi

  "$venv_python" - <<PY
$smoke_script
PY
}

print_known_limitations() {
  warn "Fallback venv note: PDBFixer/OpenMM are not installed in .venvs/cgenff because the currently available PyPI releases were not a safe cross-platform match during validation. Use the Conda/Mamba environments if you need the PDBFixer/OpenMM-assisted preparation path."
}

check_gromacs() {
  [[ "$SKIP_GROMACS_CHECK" == true ]] && return 0

  if [[ -n "${PYMACS_GMX_BIN:-}" ]]; then
    if command -v "$PYMACS_GMX_BIN" >/dev/null 2>&1; then
      log "Detected GROMACS via PYMACS_GMX_BIN=$PYMACS_GMX_BIN"
      return 0
    fi
    warn "PYMACS_GMX_BIN is set to '$PYMACS_GMX_BIN' but that executable was not found on PATH."
  fi

  local candidate
  for candidate in gmx gmx_mpi gmx-mpi; do
    if command -v "$candidate" >/dev/null 2>&1; then
      log "Detected GROMACS executable: $candidate"
      return 0
    fi
  done

  warn "No GROMACS executable detected. Install GROMACS separately and expose it as 'gmx', 'gmx_mpi', 'gmx-mpi', or set PYMACS_GMX_BIN."
}

print_instructions() {
  cat <<'EOF'

PyMACS venv setup complete.

Activate the setup / ligand parameterization environment:
  source .venvs/cgenff/bin/activate

Activate the simulation / analysis / reporting environment:
  source .venvs/mdanalysis/bin/activate

Example PyMACS usage with venvs:
  source .venvs/cgenff/bin/activate
  python 1_AutomateGromacs.py

  deactivate
  source .venvs/mdanalysis/bin/activate
  python 2_AutomateGromacs.py
  python 3A_AutomateGromacs.py
  python 4PDF4MD.py
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --recreate)
      RECREATE=true
      ;;
    --dry-run)
      DRY_RUN=true
      ;;
    --python)
      shift
      [[ $# -gt 0 ]] || die "--python requires a path argument"
      PYTHON_BIN="$1"
      ;;
    --skip-gromacs-check)
      SKIP_GROMACS_CHECK=true
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1. Use --help for usage."
      ;;
  esac
  shift
done

require_file "$CGENFF_REQ"
require_file "$MDA_REQ"
detect_python
check_python_version

log "Repository root: $SCRIPT_DIR"
run_cmd mkdir -p "$VENV_ROOT"

create_venv "$CGENFF_VENV"
create_venv "$MDA_VENV"

CGENFF_PY="$CGENFF_VENV/bin/python"
MDA_PY="$MDA_VENV/bin/python"

upgrade_packaging_tools "$CGENFF_PY"
upgrade_packaging_tools "$MDA_PY"

install_requirements "$CGENFF_PY" "$CGENFF_REQ" "cgenff"
install_requirements "$MDA_PY" "$MDA_REQ" "mdanalysis"
install_package "$MDA_PY" "mdanalysis" --no-deps "svglib==1.6.0"

CGENFF_SMOKE=$(cat <<'PY'
import sys
import Bio
import networkx
import numpy
from rdkit import Chem
print("cgenff venv Python:", sys.version)
print("cgenff smoke test OK")
PY
)

MDA_SMOKE=$(cat <<'PY'
import sys
print("mdanalysis venv Python:", sys.version, flush=True)
import DockQ
print("DockQ import OK", flush=True)
import MDAnalysis
print("MDAnalysis:", MDAnalysis.__version__, flush=True)
import matplotlib
print("matplotlib import OK", flush=True)
import mdtraj
print("mdtraj import OK", flush=True)
import networkx
print("networkx import OK", flush=True)
import numpy
print("numpy import OK", flush=True)
import pandas
print("pandas import OK", flush=True)
import reportlab
print("reportlab import OK", flush=True)
import ruptures
print("ruptures import OK", flush=True)
import seaborn
print("seaborn import OK", flush=True)
import svglib
print("svglib import OK", flush=True)
from rdkit import Chem
print("rdkit import OK", flush=True)
print("mdanalysis smoke test OK")
PY
)

run_smoke_test "$CGENFF_PY" "cgenff" "$CGENFF_SMOKE"
run_smoke_test "$MDA_PY" "mdanalysis" "$MDA_SMOKE"

print_known_limitations
check_gromacs
print_instructions
