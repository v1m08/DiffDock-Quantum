# DiffDock CPU Environment Setup Script for Windows
# Run this in PowerShell as Administrator if needed

Write-Host "Setting up DiffDock CPU environment..." -ForegroundColor Green

# Remove existing environment if it exists
Write-Host "Removing any existing diffdock-cpu environment..." -ForegroundColor Yellow
conda env remove -n diffdock-cpu -y 2>$null

# Create fresh conda environment
Write-Host "Creating fresh conda environment..." -ForegroundColor Yellow
conda create -n diffdock-cpu python=3.9.18 -y
if ($LASTEXITCODE -ne 0) { Write-Error "Failed to create conda environment"; exit 1 }

# Activate environment
Write-Host "Activating environment..." -ForegroundColor Yellow
conda activate diffdock-cpu
if ($LASTEXITCODE -ne 0) { Write-Error "Failed to activate environment"; exit 1 }

# Install basic conda packages
Write-Host "Installing basic conda packages..." -ForegroundColor Yellow
conda install setuptools=69.5.1 prody=2.2.0 scipy=1.12.0 -y
conda install -c conda-forge openbabel -y

# Install CPU PyTorch ecosystem
Write-Host "Installing PyTorch CPU..." -ForegroundColor Yellow
python -m pip install --upgrade pip wheel
python -m pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1+cpu --index-url https://download.pytorch.org/whl/cpu

# Install PyG for CPU
Write-Host "Installing PyTorch Geometric..." -ForegroundColor Yellow
python -m pip install torch-scatter==2.1.0+pt113cpu torch-sparse==0.6.16+pt113cpu torch-cluster==1.6.0+pt113cpu torch-spline-conv==1.2.1+pt113cpu --find-links https://data.pyg.org/whl/torch-1.13.1+cpu.html
python -m pip install torch-geometric==2.2.0 --find-links https://data.pyg.org/whl/torch-1.13.1+cpu.html

# Install remaining Python dependencies
Write-Host "Installing remaining Python dependencies..." -ForegroundColor Yellow
python -m pip install `
  dllogger@git+https://github.com/NVIDIA/dllogger.git `
  e3nn==0.5.1 `
  networkx==2.8.4 `
  pandas==1.5.1 `
  pybind11==2.11.1 `
  rdkit==2022.03.3 `
  scikit-learn==1.1.0 `
  pytorch-lightning==1.9.5 `
  torchmetrics==0.11.0 `
  gradio==3.50.* `
  meeko

# Try to install fair-esm (may be problematic)
Write-Host "Installing fair-esm..." -ForegroundColor Yellow
python -m pip install fair-esm[esmfold]==2.0.0

# Download and setup Vina if not already present
$vinaDir = "C:\tools\vina"
$vinaExe = "$vinaDir\vina.exe"

if (!(Test-Path $vinaExe)) {
    Write-Host "Downloading AutoDock Vina..." -ForegroundColor Yellow
    New-Item -ItemType Directory -Force -Path "C:\tools" | Out-Null
    
    # Download Vina binary
    $vinaUrl = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_win64.zip"
    $vinaZip = "$env:TEMP\vina.zip"
    
    try {
        Invoke-WebRequest -Uri $vinaUrl -OutFile $vinaZip -UseBasicParsing
        Expand-Archive -Path $vinaZip -DestinationPath "C:\tools\vina_temp" -Force
        
        # Find the actual vina.exe and copy it to our standard location
        $vinaFound = Get-ChildItem "C:\tools\vina_temp" -Filter "vina.exe" -Recurse | Select-Object -First 1
        if ($vinaFound) {
            Copy-Item $vinaFound.FullName -Destination $vinaExe
            Write-Host "Vina installed to $vinaExe" -ForegroundColor Green
        } else {
            Write-Warning "Could not find vina.exe in downloaded archive"
        }
        
        Remove-Item "C:\tools\vina_temp" -Recurse -Force -ErrorAction SilentlyContinue
        Remove-Item $vinaZip -Force -ErrorAction SilentlyContinue
    } catch {
        Write-Warning "Failed to download Vina automatically. Please download manually from GitHub releases."
    }
}

# Download and setup PDB2PQR if not present
$pdb2pqrDir = "C:\tools\pdb2pqr"
$pdb2pqrExe = "$pdb2pqrDir\pdb2pqr.exe"

if (!(Test-Path $pdb2pqrExe)) {
    Write-Host "You'll need to manually install PDB2PQR from:" -ForegroundColor Yellow
    Write-Host "https://github.com/Electrostatics/pdb2pqr/releases/latest" -ForegroundColor Cyan
    Write-Host "Install it to C:\tools\pdb2pqr or add it to your PATH" -ForegroundColor Cyan
}

# Update PATH for current session
Write-Host "Setting up PATH..." -ForegroundColor Yellow
$env:PATH = "$env:PATH;$vinaDir;$pdb2pqrDir"

# Test installations
Write-Host "Testing installations..." -ForegroundColor Yellow

# Test Python imports
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import torch_geometric; print(f'PyG: {torch_geometric.__version__}')"
python -c "import rdkit; print('RDKit: OK')"
python -c "import meeko; print('Meeko: OK')"

# Test binaries
Write-Host "Testing binaries..." -ForegroundColor Yellow
if (Test-Path $vinaExe) {
    Write-Host "Vina: OK ($vinaExe)" -ForegroundColor Green
} else {
    Write-Host "Vina: NOT FOUND" -ForegroundColor Red
}

try {
    obabel -h 2>$null | Out-Null
    Write-Host "OpenBabel: OK" -ForegroundColor Green
} catch {
    Write-Host "OpenBabel: NOT FOUND" -ForegroundColor Red
}

if (Test-Path $pdb2pqrExe) {
    Write-Host "PDB2PQR: OK ($pdb2pqrExe)" -ForegroundColor Green
} else {
    Write-Host "PDB2PQR: MANUAL INSTALL NEEDED" -ForegroundColor Yellow
}

Write-Host "`nSetup complete! To activate this environment in the future:" -ForegroundColor Green
Write-Host "conda activate diffdock-cpu" -ForegroundColor Cyan
Write-Host "`nIf you need to set PATH for Vina/PDB2PQR in new sessions:" -ForegroundColor Green
Write-Host "`$env:PATH = `"`$env:PATH;$vinaDir;$pdb2pqrDir`"" -ForegroundColor Cyan

Write-Host "`nNext steps:" -ForegroundColor Green
Write-Host "1. If PDB2PQR install is needed, download and install it" -ForegroundColor Yellow
Write-Host "2. Run: python local_runner/diffdock_run.py --help" -ForegroundColor Yellow
Write-Host "3. Then: python local_runner/protonation_scoring.py --help" -ForegroundColor Yellow