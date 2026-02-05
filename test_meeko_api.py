#!/usr/bin/env python3
"""Test meeko's proper PDBQT writing API for Vina"""
import sys
from pathlib import Path
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

# Load ligand
mol = Chem.SDMolSupplier('/mnt/c/Users/manne/DiffDock/data/1a0q/1a0q_ligand.sdf')[0]
mol = Chem.AddHs(mol)

if mol.GetNumConformers() == 0:
    AllChem.EmbedMolecule(mol, randomSeed=42)

AllChem.UFFOptimizeMolecule(mol)
AllChem.ComputeGasteigerCharges(mol)

# Use meeko's new API
mmp = MoleculePreparation()
setups = mmp.prepare(mol)

print(f"Got {len(setups)} setup(s)")
setup = setups[0]

# Write PDBQT using the proper API
writer = PDBQTWriterLegacy()
result = writer.write_string(setup)

# Result is a tuple (pdbqt_string, warnings)
if isinstance(result, tuple):
    pdbqt_string = result[0]
else:
    pdbqt_string = result

# Write to file
with open('/tmp/test_meeko_api.pdbqt', 'w') as f:
    f.write(pdbqt_string)

# Print first 30 lines
print("\n=== PDBQT Output (first 30 lines) ===")
for line in pdbqt_string.split('\n')[:30]:
    print(repr(line))
