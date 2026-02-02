# Protonation-Aware Molecular Docking

This project implements a protonation-aware scoring system that analyzes the binding energy differences between base and protonated forms of ligands. This is useful for understanding pH-dependent binding behavior and optimizing drug design.

## Quick Start (Windows)

### 1. Setup (One Time)
```powershell
.\setup_cpu_env.ps1
```

### 2. Run Analysis (Default Example)
```batch
.\run_protonation_analysis.bat
```

### 3. Run Analysis (Custom Files)
```batch
python local_runner/standalone_protonation_scorer.py --protein_pdbqt YOUR_PROTEIN.pdbqt --ligand YOUR_LIGAND.sdf --center 10.5 15.2 8.7 --size 20 20 20 --out_dir results/my_analysis
```

## What It Does

1. **Generates Protonation States**: Creates both neutral (base) and protonated forms of your ligand
2. **Estimates Binding Energies**: Uses molecular properties and simple geometric features to estimate binding affinities
3. **Calculates Differences**: Computes the energy difference between protonation states
4. **Provides Interpretation**: Tells you which form binds better and by how much

## Example Output

```
============================================================
PROTONATION ANALYSIS RESULTS
============================================================
Base form energy:        6.57 kcal/mol
Protonated form energy:  8.01 kcal/mol
Protonation difference:  1.44 kcal/mol
Structural RMSD change:  0.00 Å

INTERPRETATION:
✅ Base form is significantly MORE FAVORABLE for binding
   → Consider targeting the neutral state in drug design
💡 Protonation causes a substantial energy change of 1.4 kcal/mol
   → This suggests pH-sensitive binding behavior
============================================================
```

## Files Created

- `local_runner/standalone_protonation_scorer.py` - Main scoring engine
- `run_protonation_analysis.bat` - Easy Windows batch script
- `setup_cpu_env.ps1` - Environment setup script
- `test_setup.ps1` - Test script to verify installation

## Requirements

- Windows 10/11
- Miniconda/Anaconda
- ~2GB disk space for dependencies

## Understanding Results

- **Negative difference**: Protonated form binds better
- **Positive difference**: Base form binds better  
- **|Difference| > 1.0 kcal/mol**: Significant pH-dependent binding
- **|Difference| < 0.5 kcal/mol**: Similar binding affinities

## Supported File Formats

- Proteins: PDBQT format (use OpenBabel to convert PDB files)
- Ligands: SDF, MOL2 formats

## Next Steps

To implement the quantum-inspired machine learning methods from the paper (https://arxiv.org/pdf/2401.12999), we would need to:

1. Implement quantum feature encoding for molecular representations
2. Add variational quantum classifiers for binding prediction
3. Integrate quantum-enhanced optimization algorithms
4. Train on larger datasets with experimental binding affinities

The current system provides the foundation for protonation-aware scoring that could be enhanced with these quantum methods.