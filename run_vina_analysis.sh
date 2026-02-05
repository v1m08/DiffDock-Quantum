#!/bin/bash
# Accurate Vina-based scoring script
# Usage: bash run_vina_analysis.sh

echo "=========================================="
echo "PROTONATION ANALYSIS - Vina (Accurate)"
echo "=========================================="
echo ""

# Check if conda environment is active
if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
    echo "Activating conda environment..."
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate diffdock-vina
fi

echo "Environment: $CONDA_DEFAULT_ENV"
echo ""

# Verify Vina is installed
echo "Checking Vina installation..."
python -c "from vina import Vina; print('✓ Vina is available')" || {
    echo "✗ Vina not installed!"
    echo "Install with: pip install vina"
    exit 1
}

echo ""
echo "Running Vina docking analysis..."
echo "This may take 1-2 minutes..."
echo ""

# Run the Vina-based scorer
python local_runner/wsl_vina_protonation_scorer.py \
    --protein_pdbqt data/1a0q/1a0q_protein_processed.pdbqt \
    --ligand data/1a0q/1a0q_ligand.sdf \
    --center 0 0 0 \
    --size 20 20 20 \
    --out_dir results/vina_protonation_analysis

echo ""
echo "=========================================="
echo "✓ Vina Analysis Complete!"
echo "=========================================="
echo "Results saved to: results/vina_protonation_analysis/"
echo ""
echo "Compare results:"
echo "- Quick (simplified): results/quick_protonation_analysis/protonation_results.csv"
echo "- Vina (accurate):    results/vina_protonation_analysis/vina_protonation_results.csv"
echo ""
