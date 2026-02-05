#!/bin/bash
# WSL Environment Setup Script for DiffDock with Vina
# This script sets up a complete development environment in WSL

set -e  # Exit on any error

echo "=========================================="
echo "DiffDock + Vina WSL Setup"
echo "=========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_step() {
    echo -e "${BLUE}>>> $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# Step 1: Update system packages
print_step "Step 1: Updating system packages..."
sudo apt-get update
sudo apt-get upgrade -y
print_success "System packages updated"

# Step 2: Install build tools and dependencies
print_step "Step 2: Installing build tools and system dependencies..."
sudo apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    libboost-all-dev \
    python3-dev \
    python3-pip \
    python3-venv
print_success "Build tools installed"

# Step 3: Install Miniconda if not already installed
print_step "Step 3: Checking Miniconda installation..."
if ! command -v conda &> /dev/null; then
    print_step "Installing Miniconda..."
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh
    ~/miniconda3/bin/conda init bash
    source ~/.bashrc
    print_success "Miniconda installed"
else
    print_success "Miniconda already installed"
fi

# Step 4: Create conda environment
print_step "Step 4: Creating conda environment 'diffdock-vina'..."
conda create -n diffdock-vina python=3.9 -y
print_success "Conda environment created"

# Step 5: Activate environment and install PyTorch
print_step "Step 5: Installing PyTorch (CPU version)..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate diffdock-vina
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
print_success "PyTorch installed"

# Step 6: Install core dependencies
print_step "Step 6: Installing core molecular chemistry packages..."
pip install \
    rdkit \
    scipy \
    numpy \
    pandas \
    matplotlib \
    jupyter \
    biopython \
    prody
print_success "Core packages installed"

# Step 7: Install PyTorch Geometric
print_step "Step 7: Installing PyTorch Geometric..."
pip install torch-geometric
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
print_success "PyTorch Geometric installed"

# Step 8: Install Meeko
print_step "Step 8: Installing Meeko for ligand preparation..."
pip install meeko
print_success "Meeko installed"

# Step 9: Install Vina (THE KEY STEP - works on Linux!)
print_step "Step 9: Installing AutoDock Vina..."
pip install vina
print_success "Vina installed successfully!"

# Step 10: Verify installations
print_step "Step 10: Verifying installations..."
echo "Python version:"
python --version
echo ""
echo "Checking imports..."
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import rdkit; print('RDKit: OK')"
python -c "import meeko; print('Meeko: OK')"
python -c "from vina import Vina; print('Vina: OK')"
print_success "All verifications passed!"

echo ""
echo "=========================================="
echo "✓ WSL Environment Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Copy your project to WSL home: cp -r /mnt/c/Users/manne/DiffDock ~/DiffDock"
echo "2. Enter the project: cd ~/DiffDock"
echo "3. Activate environment: conda activate diffdock-vina"
echo "4. Run analysis: bash run_vina_analysis.sh"
echo ""
