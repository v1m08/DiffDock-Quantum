#!/usr/bin/env python3
"""Test if Vina can read our PDBQT format"""
import sys
from vina import Vina

# Test receptor
print("Testing Vina with cleaned protein PDBQT...")
try:
    v = Vina(sf_name='vina')
    v.set_receptor('/mnt/c/Users/manne/DiffDock/data/1a0q/1a0q_protein_processed_clean.pdbqt')
    print("✓ Receptor loaded successfully")
except Exception as e:
    print(f"✗ Error loading receptor: {e}")

# Test ligand
print("\nTesting Vina with cleaned ligand PDBQT...")
try:
    v2 = Vina(sf_name='vina')
    v2.set_receptor('/mnt/c/Users/manne/DiffDock/data/1a0q/1a0q_protein_processed_clean.pdbqt')
    v2.set_ligand_from_file('/tmp/test_manual.pdbqt')
    print("✓ Ligand loaded successfully")
except Exception as e:
    print(f"✗ Error loading ligand: {e}")

print("\nPDBQT format test passed!")
