#!/usr/bin/env python3
"""Debug PDBQT cleaning"""
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.SDMolSupplier('/mnt/c/Users/manne/DiffDock/data/1a0q/1a0q_ligand.sdf')[0]
mol = Chem.AddHs(mol)
if mol.GetNumConformers() == 0:
    AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.UFFOptimizeMolecule(mol)
AllChem.ComputeGasteigerCharges(mol)

mmp = MoleculePreparation()
setups = mmp.prepare(mol)
writer = PDBQTWriterLegacy()
result = writer.write_string(setups[0])
pdbqt = result[0] if isinstance(result, tuple) else result

print("=== Before Cleaning ===")
for i, line in enumerate(pdbqt.split('\n')[:10], 1):
    print(f'{i}: {repr(line)}')

# Clean
cleaned = []
for line in pdbqt.split('\n'):
    if any(tag in line for tag in ['ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF']):
        continue
    cleaned.append(line)

print("\n=== After Cleaning ===")
for i, line in enumerate(cleaned[:10], 1):
    print(f'{i}: {repr(line)}')

# Write to file
with open('/tmp/test_cleaned.pdbqt', 'w') as f:
    f.write('\n'.join(cleaned))

print("\n=== File Contents ===")
with open('/tmp/test_cleaned.pdbqt') as f:
    for i, line in enumerate(f.readlines()[:10], 1):
        print(f'{i}: {repr(line)}')
