#!/usr/bin/env python3
"""Test writing PDBQT that Vina will accept"""
import sys
sys.path.insert(0, '/mnt/c/Users/manne/DiffDock')

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# Load ligand
mol = Chem.SDMolSupplier('/mnt/c/Users/manne/DiffDock/data/1a0q/1a0q_ligand.sdf')[0]
mol = Chem.AddHs(mol)

if mol.GetNumConformers() == 0:
    AllChem.EmbedMolecule(mol, randomSeed=42)

AllChem.UFFOptimizeMolecule(mol)
AllChem.ComputeGasteigerCharges(mol)

# Write manual PDBQT with exact format Vina expects
with open('/tmp/test_manual.pdbqt', 'w') as f:
    for i, atom in enumerate(mol.GetAtoms()):
        conf = mol.GetConformer()
        pos = conf.GetAtomPosition(i)
        
        try:
            charge = float(atom.GetDoubleProp('_GasteigerCharge'))
        except:
            charge = 0.0
        
        atom_type = atom.GetSymbol()
        if atom_type == 'O' and atom.GetTotalDegree() == 1:
            atom_type = 'OA'
        
        # Format: ATOM serial name resname chain resnum X Y Z occ Bfactor charge type
        # Example: ATOM      1  C   UNL     1      12.798  20.933  61.805  1.00  0.00     0.221 C
        f.write(f"ATOM  {i+1:5d}  {atom.GetSymbol():2s} {' ':1s}LIG A   1    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00    {charge:7.4f} {atom_type:2s}\n")

with open('/tmp/test_manual.pdbqt') as f:
    for line in f.readlines()[:20]:
        print(repr(line))
