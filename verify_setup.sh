#!/bin/bash
# Verification script - checks if all files are in place
# Run this after setup to verify everything worked

echo "════════════════════════════════════════════════════════════"
echo "WSL SETUP VERIFICATION"
echo "════════════════════════════════════════════════════════════"
echo ""

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

check_file() {
    if [ -f "$1" ]; then
        echo -e "${GREEN}✓${NC} $1"
        return 0
    else
        echo -e "${RED}✗${NC} $1"
        return 1
    fi
}

check_command() {
    if command -v "$1" &> /dev/null; then
        echo -e "${GREEN}✓${NC} $1 is installed"
        return 0
    else
        echo -e "${RED}✗${NC} $1 is NOT installed"
        return 1
    fi
}

echo "Checking documentation files..."
check_file "START_HERE.txt"
check_file "COPY_PASTE_COMMANDS.txt"
check_file "QUICK_REFERENCE.txt"
check_file "WSL_SETUP_GUIDE.txt"
check_file "README_WSL.md"
check_file "FILES_CREATED_SUMMARY.md"
check_file "QUICK_START.txt"

echo ""
echo "Checking setup scripts..."
check_file "setup_wsl_env.sh"

echo ""
echo "Checking analysis scripts..."
check_file "run_protonation_analysis.sh"
check_file "run_vina_analysis.sh"

echo ""
echo "Checking Python scripts..."
check_file "local_runner/wsl_vina_protonation_scorer.py"
check_file "local_runner/standalone_protonation_scorer.py"

echo ""
echo "Checking data files..."
check_file "data/1a0q/1a0q_protein_processed.pdbqt"
check_file "data/1a0q/1a0q_ligand.sdf"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "Checking environment..."
echo "════════════════════════════════════════════════════════════"
echo ""

check_command "python"
check_command "conda"
check_command "pip"

echo ""
if [ ! -z "$CONDA_DEFAULT_ENV" ]; then
    echo -e "${GREEN}✓${NC} Conda environment active: $CONDA_DEFAULT_ENV"
else
    echo -e "${YELLOW}⚠${NC}  No conda environment active (run: conda activate diffdock-vina)"
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "Checking Python packages..."
echo "════════════════════════════════════════════════════════════"
echo ""

python << 'PYEOF'
import sys

packages = [
    ('torch', 'PyTorch'),
    ('rdkit', 'RDKit'),
    ('meeko', 'Meeko'),
    ('vina', 'Vina'),
    ('torch_geometric', 'PyTorch Geometric')
]

all_ok = True
for module_name, display_name in packages:
    try:
        __import__(module_name)
        print(f"✓ {display_name} is installed")
    except ImportError:
        print(f"✗ {display_name} is NOT installed")
        all_ok = False

print()
if all_ok:
    print("✓ All packages installed!")
    sys.exit(0)
else:
    print("✗ Some packages missing - run setup_wsl_env.sh")
    sys.exit(1)
PYEOF

echo ""
echo "════════════════════════════════════════════════════════════"
echo "VERIFICATION SUMMARY"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "If all items show ✓, you're ready to go!"
echo ""
echo "Next steps:"
echo "1. cd ~/DiffDock"
echo "2. conda activate diffdock-vina"
echo "3. bash run_vina_analysis.sh"
echo ""
