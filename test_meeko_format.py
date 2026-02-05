#!/usr/bin/env python3
"""Test meeko PDBQT format generation"""
import sys
from pathlib import Path
from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem

# Load ligand from SDF
mol = Chem.SDMolSupplier('data/1a0q/1a0q_ligand.sdf')[0]

# Add explicit hydrogens
mol = Chem.AddHs(mol)

# Add 3D coords
if mol.GetNumConformers() == 0:
    AllChem.EmbedMolecule(mol, randomSeed=42)

# Generate charges
AllChem.ComputeGasteigerCharges(mol)

# Generate with meeko
mmp = MoleculePreparation()
mmp.prepare(mol)
mmp.write_pdbqt_file('/tmp/test_meeko.pdbqt')

# Print first 30 lines
with open('/tmp/test_meeko.pdbqt') as f:
    lines = f.readlines()[:30]
    for line in lines:
        print(line, end='')


