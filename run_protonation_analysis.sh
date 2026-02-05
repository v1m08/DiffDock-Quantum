#!/bin/bash
# Quick run script for protonation-aware scoring (simplified version)
# Usage: bash run_protonation_analysis.sh

echo "=========================================="
echo "PROTONATION ANALYSIS - Quick Version"
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

# Run the standalone (simplified) scorer
python local_runner/standalone_protonation_scorer.py \
    --protein_pdbqt data/1a0q/1a0q_protein_processed.pdbqt \
    --ligand data/1a0q/1a0q_ligand.sdf \
    --center 0 0 0 \
    --size 20 20 20 \
    --out_dir results/quick_protonation_analysis

echo ""
echo "=========================================="
echo "✓ Analysis Complete!"
echo "=========================================="
echo "Results saved to: results/quick_protonation_analysis/"
echo ""
