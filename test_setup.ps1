# Quick Test Script for DiffDock CPU Setup
Write-Host "Testing DiffDock CPU setup..." -ForegroundColor Green

# Activate environment
conda activate diffdock-cpu

# Set PATH for this session
$env:PATH = "$env:PATH;C:\tools\vina;C:\tools\pdb2pqr"

# Test basic functionality
Write-Host "Testing Python imports..." -ForegroundColor Yellow
python -c @"
try:
    import torch
    import torch_geometric
    import rdkit
    import meeko
    print('All Python imports successful')
    print(f'  - PyTorch: {torch.__version__}')
    print(f'  - PyG: {torch_geometric.__version__}')
    print('  - RDKit: OK')
    print('  - Meeko: OK')
except ImportError as e:
    print(f'Import error: {e}')
    exit(1)
"@

# Test binaries
Write-Host "Testing binaries..." -ForegroundColor Yellow
$vinaPath = "C:\tools\vina\vina.exe"
if (Test-Path $vinaPath) {
    Write-Host "Vina found at $vinaPath" -ForegroundColor Green
} else {
    Write-Host "Vina not found at $vinaPath" -ForegroundColor Red
    Write-Host "  Please run the setup script first" -ForegroundColor Yellow
}

try {
    $result = obabel -h 2>&1 | Out-String
    Write-Host "OpenBabel is working" -ForegroundColor Green
} catch {
    Write-Host "OpenBabel not found" -ForegroundColor Red
}

$pdb2pqrPath = "C:\tools\pdb2pqr\pdb2pqr.exe"
if (Test-Path $pdb2pqrPath) {
    Write-Host "PDB2PQR found at $pdb2pqrPath" -ForegroundColor Green
} else {
    Write-Host "PDB2PQR not found at $pdb2pqrPath" -ForegroundColor Red
}

# Test DiffDock scripts
Write-Host "Testing DiffDock scripts..." -ForegroundColor Yellow
try {
    python local_runner/diffdock_run.py --help > $null
    Write-Host "DiffDock runner script is working" -ForegroundColor Green
} catch {
    Write-Host "DiffDock runner script has issues" -ForegroundColor Red
}

try {
    python local_runner/protonation_scoring.py --help > $null
    Write-Host "Protonation scoring script is working" -ForegroundColor Green
} catch {
    Write-Host "Protonation scoring script has issues" -ForegroundColor Red
}

Write-Host "`nIf all tests passed, you can run:" -ForegroundColor Green
Write-Host "python local_runner/protonation_scoring.py --protein_pdb data/1a0q/1a0q_protein_processed.pdb --protein_pdbqt data/1a0q/1a0q_protein_processed.pdbqt --ligand data/1a0q/1a0q_ligand.sdf --center 0 0 0 --size 20 20 20 --out_dir results/test" -ForegroundColor Cyan