#!/usr/bin/env python3
"""
Vina-based protonation-aware scoring system.
Uses AutoDock Vina for accurate binding energy calculations.
Designed for WSL/Linux with pip-installed vina package.
Uses pdb2pqr/custom PDBQT writing for ligand preparation.
"""

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import Tuple, Dict, Any
import subprocess
import re
import os
import shutil

from rdkit import Chem
from rdkit.Chem import MolStandardize, Descriptors, AllChem

try:
    from vina import Vina
    VINA_AVAILABLE = True
except ImportError:
    print("Warning: Vina module not found. Install with: pip install vina")
    VINA_AVAILABLE = False

# Check for vina command-line tool
try:
    result = subprocess.run(['vina', '--help'], capture_output=True, timeout=5)
    VINA_CLI_AVAILABLE = result.returncode == 0
except (FileNotFoundError, subprocess.TimeoutExpired):
    VINA_CLI_AVAILABLE = False


def embed_and_optimize(mol: Chem.Mol) -> Chem.Mol:
    """Add hydrogens, embed 3D coords, and optimize geometry."""
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def write_pdb_from_mol(mol: Chem.Mol, pdb_path: str) -> None:
    """Write minimal PDB file from RDKit molecule."""
    with open(pdb_path, 'w') as f:
        f.write("REMARK Written by RDKit\n")
        conf = mol.GetConformer()
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_name = atom.GetSymbol().rjust(2)
            f.write(
                f"ATOM  {i+1:5d} {atom_name:<4s}LIG A   1    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00           {atom.GetSymbol():>2s}\n"
            )
        f.write("END\n")


def clean_ligand_pdbqt(pdbqt_path: str) -> None:
    """Remove ROOT/BRANCH tags from ligand PDBQT in-place."""
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()
    
    # Remove ROOT/BRANCH/TORSDOF lines, keep ATOM records
    with open(pdbqt_path, 'w') as f:
        for line in lines:
            if any(tag in line for tag in ['ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF']):
                continue  # Skip these tags
            elif line.startswith('ATOM'):
                f.write(line)
            elif line.startswith('REMARK') or line.startswith('CONECT') or line.startswith('END'):
                f.write(line)


def write_pdbqt_from_mol(mol: Chem.Mol, pdbqt_path: str) -> None:
    """Write PDBQT file from RDKit molecule with Gasteiger charges - Vina 1.1.2 compatible."""
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except:
        pass
    
    with open(pdbqt_path, 'w') as f:
        for i, atom in enumerate(mol.GetAtoms()):
            conf = mol.GetConformer()
            pos = conf.GetAtomPosition(i)
            
            # Get partial charge
            try:
                charge = float(atom.GetDoubleProp('_GasteigerCharge'))
            except:
                charge = float(atom.GetFormalCharge())
            
            # AutoDock atom type
            atom_sym = atom.GetSymbol()
            atom_type = 'C'
            if atom_sym == 'N':
                atom_type = 'N'
            elif atom_sym == 'O':
                atom_type = 'OA'
            elif atom_sym == 'S':
                atom_type = 'S'
            elif atom_sym == 'P':
                atom_type = 'P'
            elif atom_sym == 'H':
                atom_type = 'HD'
            elif atom_sym == 'F':
                atom_type = 'F'
            elif atom_sym == 'Cl':
                atom_type = 'Cl'
            elif atom_sym == 'Br':
                atom_type = 'Br'
            
            # Standard PDB/PDBQT format
            # Column positions: 1-6(record), 7-11(serial), 13-16(name), 18(alt), 18-20(resname), 22(chain), 23-26(resnum), 31-38(x), 39-46(y), 47-54(z), 55-60(occ), 61-66(Bfactor), 77-78(element)
            serial = i + 1
            name = f" {atom_sym} "[:4]  # Pad to 4 chars
            
            f.write(f"ATOM  {serial:5d}{name:4s}  LIG A   1    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00 {charge:6.2f}          {atom_sym:>2s}\n")
        
        f.write("END\n")


def generate_base_and_protonated(mol: Chem.Mol) -> Tuple[Chem.Mol, Chem.Mol]:
    """Generate base and protonated forms of a molecule."""
    try:
        # Neutral (base) form: strip formal charges, then embed/optimize with hydrogens
        base = Chem.Mol(mol)
        for atom in base.GetAtoms():
            atom.SetFormalCharge(0)
        Chem.SanitizeMol(base)
        base = embed_and_optimize(base)

        # Protonated form: add a proton to best basic site (prefer aliphatic N, fallback O)
        rwmol = Chem.RWMol(mol)
        target_idx = None
        for atom in rwmol.GetAtoms():
            if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic() and atom.GetTotalValence() <= 3:
                target_idx = atom.GetIdx()
                break
        if target_idx is None:
            for atom in rwmol.GetAtoms():
                if atom.GetAtomicNum() == 8 and atom.GetTotalValence() <= 2:
                    target_idx = atom.GetIdx()
                    break
        if target_idx is not None:
            atom = rwmol.GetAtomWithIdx(target_idx)
            atom.SetFormalCharge(atom.GetFormalCharge() + 1)
            print(f"Protonated atom index {target_idx} (atomic num {atom.GetAtomicNum()})")
        else:
            print("No clear protonation site found; using original formal charges")

        protonated = rwmol.GetMol()
        Chem.SanitizeMol(protonated)
        protonated = embed_and_optimize(protonated)

        return base, protonated

    except Exception as e:
        print(f"Warning: Error generating protonation states: {e}")
        base = embed_and_optimize(mol)

        protonated = Chem.RWMol(mol)
        for atom in protonated.GetAtoms():
            if atom.GetAtomicNum() in [7, 8]:
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                break
        protonated = protonated.GetMol()
        protonated = embed_and_optimize(protonated)

        return base, protonated


def prepare_ligand_pdbqt(mol: Chem.Mol, output_path: str) -> bool:
    """Prepare ligand PDBQT using pdb2pqr (charges) then convert to Vina-compatible PDBQT."""
    try:
        # Ensure pdb2pqr exists
        pdb2pqr_exe = shutil.which("pdb2pqr")
        if not pdb2pqr_exe:
            print("Error: pdb2pqr not found in PATH. Install with: conda install -c conda-forge pdb2pqr")
            return False

        # Preserve existing geometry; if missing, embed and optimize
        mol = Chem.Mol(mol)
        if mol.GetNumConformers() == 0:
            mol = embed_and_optimize(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pdb_path = tmpdir / "ligand.pdb"
            pqr_path = tmpdir / "ligand.pqr"

            # Write minimal PDB from RDKit
            write_pdb_from_mol(mol, str(pdb_path))

            # Run pdb2pqr to compute charges/radii (keep geometry)
            cmd = [pdb2pqr_exe, "--ff=AMBER", "--keep-chain", "--noopt", str(pdb_path), str(pqr_path)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode != 0 or not pqr_path.exists():
                print("pdb2pqr failed:\n" + result.stderr)
                return False

            # Convert PQR to PDBQT (use charge column from PQR)
            pdbqt_lines = []
            with open(pqr_path, "r") as f:
                for line in f:
                    if not (line.startswith("ATOM") or line.startswith("HETATM")):
                        continue
                    # PQR columns: record, serial, name, resName, chain, resNum, x, y, z, charge, radius
                    parts = line.split()
                    if len(parts) < 10:
                        continue
                    record, serial, name, resn, chain, resi = parts[0:6]
                    x, y, z, charge = map(float, parts[6:10])

                    elem = ''.join([c for c in name if c.isalpha()])
                    elem = elem.capitalize() if elem else "C"
                    # AutoDock atom type heuristic
                    atom_type = {
                        'C': 'C', 'N': 'N', 'O': 'OA', 'S': 'S', 'P': 'P',
                        'H': 'HD', 'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I'
                    }.get(elem, elem)

                    atom_name = name[:4].ljust(4)
                    pdbqt_lines.append(
                        f"ATOM  {int(serial):5d} {atom_name:<4s}LIG A{int(resi):4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    {charge:7.4f} {atom_type:>2s}"
                    )

            if not pdbqt_lines:
                print("Error: no atoms parsed from pdb2pqr output")
                return False

            with open(output_path, "w") as f:
                f.write("\n".join(pdbqt_lines))
                f.write("\nEND\n")

            print(f"Ligand prepared with pdb2pqr: {output_path}")
            return True

    except Exception as e:
        print(f"Error preparing ligand: {e}")
        import traceback
        traceback.print_exc()
        return False


def clean_protein_pdbqt(pdbqt_path: str) -> str:
    """Remove PDBQT control tags that cause Vina issues."""
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()
    
    # Check if cleaning needed
    content = ''.join(lines)
    has_tags = any(tag in content for tag in ['ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH'])
    if not has_tags:
        return pdbqt_path
    
    # Write cleaned version - only keep valid PDBQT records
    clean_path = pdbqt_path.replace('.pdbqt', '_clean.pdbqt')
    with open(clean_path, 'w') as f:
        for line in lines:
            # Skip control tags
            if any(x in line for x in ['ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF']):
                continue
            # Keep valid records
            if any(line.startswith(tag) for tag in ['ATOM', 'HETATM', 'CONECT', 'END', 'REMARK']):
                f.write(line)
    
    print(f"Created cleaned PDBQT: {clean_path}")
    return clean_path


def run_vina_docking(protein_pdbqt: str, ligand_pdbqt: str, 
                     center: Tuple[float, float, float],
                     size: Tuple[float, float, float],
                     output_pdbqt: str) -> float:
    """Run AutoDock Vina and return binding energy."""
    
    if not VINA_AVAILABLE:
        print("Error: Vina not installed. Run: pip install vina")
        return None
    
    try:
        # Verify input files exist
        if not Path(protein_pdbqt).exists():
            print(f"Error: Protein PDBQT not found: {protein_pdbqt}")
            return None
        if not Path(ligand_pdbqt).exists():
            print(f"Error: Ligand PDBQT not found: {ligand_pdbqt}")
            return None
        
        # Clean protein PDBQT if needed (remove ROOT/ENDROOT tags)
        protein_pdbqt = clean_protein_pdbqt(protein_pdbqt)
        
        # Use vina command-line tool if available
        if VINA_CLI_AVAILABLE:
            print(f"Using AutoDock Vina command-line tool...")
            center_str = ' '.join(map(str, center))
            size_str = ' '.join(map(str, size))
            
            cmd = [
                'vina',
                '--receptor', protein_pdbqt,
                '--ligand', ligand_pdbqt,
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(size[0]),
                '--size_y', str(size[1]),
                '--size_z', str(size[2]),
                '--out', output_pdbqt,
                '--exhaustiveness', '8',
                '--num_modes', '20',
                '--energy_range', '3.0'
            ]
            
            print(f"Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                print(f"Vina error output:\n{result.stderr}")
                return None
            
            # Parse output to get best energy
            for line in result.stderr.split('\n'):
                if 'Affinity' in line and 'kcal/mol' in line:
                    try:
                        # Line format: "  1        -6.5      0.000      0.000"
                        parts = line.split()
                        if parts and '-' in parts[1]:
                            best_energy = float(parts[1])
                            print(f"Best binding energy: {best_energy:.2f} kcal/mol")
                            return best_energy
                    except (ValueError, IndexError):
                        pass
            
            return None
        
        # Fallback to Python API
        print(f"Using Vina Python API...")
        v = Vina(sf_name='vina')
        
        # Set up receptor
        print(f"Setting receptor: {protein_pdbqt}")
        v.set_receptor(protein_pdbqt)
        
        # Set up ligand
        print(f"Setting ligand: {ligand_pdbqt}")
        v.set_ligand_from_file(ligand_pdbqt)
        
        # Define search space
        print(f"Computing maps with center {center} and size {size}...")
        v.compute_vina_maps(center=center, box_size=size)
        
        # Run docking
        print(f"Running Vina docking...")
        v.dock(exhaustiveness=8, n_poses=20)
        
        # Write output
        v.write_poses(output_pdbqt, n_poses=20, energy_range=3.0)
        
        # Get best binding energy
        best_energy = v.energy_best()
        print(f"Best binding energy: {best_energy:.2f} kcal/mol")
        
        return best_energy
        
    except Exception as e:
        print(f"Error running Vina: {e}")
        import traceback
        traceback.print_exc()
        return None


def calculate_molecular_properties(mol: Chem.Mol) -> Dict[str, float]:
    """Calculate molecular properties."""
    try:
        return {
            'mw': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'tpsa': Descriptors.TPSA(mol),
            'charge': Chem.rdmolops.GetFormalCharge(mol),
            'num_rotatable': Descriptors.NumRotatableBonds(mol)
        }
    except:
        return {k: 0.0 for k in ['mw', 'logp', 'hbd', 'hba', 'tpsa', 'charge', 'num_rotatable']}


def score_protonation_states(protein_pdbqt: str, ligand_sdf: str, 
                            center: Tuple[float, float, float], 
                            size: Tuple[float, float, float],
                            out_dir: Path) -> Dict[str, Any]:
    """Score protonation states using Vina."""
    
    # Load ligand
    if ligand_sdf.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
        mol = supplier[0]
    elif ligand_sdf.endswith('.mol2'):
        mol = Chem.MolFromMol2File(ligand_sdf, removeHs=False)
    else:
        raise ValueError(f"Unsupported ligand format: {ligand_sdf}")
    
    if mol is None:
        raise ValueError(f"Could not read ligand from {ligand_sdf}")
    
    print(f"Original molecule: {Chem.MolToSmiles(mol)}")
    
    # Generate protonation states
    base_mol, prot_mol = generate_base_and_protonated(mol)
    
    print(f"Base form: {Chem.MolToSmiles(base_mol)} (charge: {Chem.rdmolops.GetFormalCharge(base_mol)})")
    print(f"Protonated form: {Chem.MolToSmiles(prot_mol)} (charge: {Chem.rdmolops.GetFormalCharge(prot_mol)})")
    
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Prepare ligands in PDBQT format
        base_pdbqt = tmpdir / "base.pdbqt"
        prot_pdbqt = tmpdir / "protonated.pdbqt"
        
        print("\nPreparing ligands...")
        if not prepare_ligand_pdbqt(base_mol, str(base_pdbqt)):
            raise ValueError("Failed to prepare base ligand")
        if not prepare_ligand_pdbqt(prot_mol, str(prot_pdbqt)):
            raise ValueError("Failed to prepare protonated ligand")
        
        # Run Vina docking for both forms
        print("\n" + "="*60)
        print("VINA DOCKING: Base Form")
        print("="*60)
        base_out = tmpdir / "base_out.pdbqt"
        base_energy = run_vina_docking(protein_pdbqt, str(base_pdbqt), center, size, str(base_out))
        
        print("\n" + "="*60)
        print("VINA DOCKING: Protonated Form")
        print("="*60)
        prot_out = tmpdir / "prot_out.pdbqt"
        prot_energy = run_vina_docking(protein_pdbqt, str(prot_pdbqt), center, size, str(prot_out))
        
        if base_energy is None or prot_energy is None:
            raise ValueError("Vina docking failed")
        
        # Calculate difference
        protonation_energy_diff = prot_energy - base_energy
        
        # Calculate properties
        base_props = calculate_molecular_properties(base_mol)
        prot_props = calculate_molecular_properties(prot_mol)
        
        return {
            'base_energy': base_energy,
            'protonated_energy': prot_energy,
            'protonation_difference': protonation_energy_diff,
            'base_properties': base_props,
            'protonated_properties': prot_props,
            'base_smiles': Chem.MolToSmiles(base_mol),
            'prot_smiles': Chem.MolToSmiles(prot_mol)
        }


def main():
    parser = argparse.ArgumentParser(description='Vina-based protonation scoring')
    parser.add_argument('--protein_pdbqt', required=True, help='Protein PDBQT file')
    parser.add_argument('--ligand', required=True, help='Ligand SDF/MOL2 file')
    parser.add_argument('--center', nargs=3, type=float, required=True, 
                       help='Binding site center (x y z)')
    parser.add_argument('--size', nargs=3, type=float, required=True,
                       help='Binding site size (x y z)')
    parser.add_argument('--out_dir', default='results/vina_protonation',
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Verify Vina is available
    if not VINA_AVAILABLE:
        print("Error: Vina module not available!")
        print("Install it with: pip install vina")
        return 1
    
    # Create output directory
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Analyzing protonation states with VINA for {args.ligand}...")
    print(f"Binding site center: {args.center}")
    print(f"Binding site size: {args.size}")
    print(f"Output directory: {out_dir}")
    
    try:
        results = score_protonation_states(
            args.protein_pdbqt, args.ligand,
            tuple(args.center), tuple(args.size),
            out_dir
        )
        
        print("\n" + "="*60)
        print("VINA BINDING AFFINITY ANALYSIS")
        print("="*60)
        print(f"Base form energy:        {results['base_energy']:.2f} kcal/mol")
        print(f"Protonated form energy:  {results['protonated_energy']:.2f} kcal/mol")
        print(f"Protonation difference:  {results['protonation_difference']:.2f} kcal/mol")
        
        print("\nMolecular Properties Comparison:")
        base_props = results['base_properties']
        prot_props = results['protonated_properties']
        
        print(f"{'Property':<15} {'Base':<12} {'Protonated':<14} {'Change':<10}")
        print("-" * 55)
        for prop in ['mw', 'logp', 'hbd', 'hba', 'charge']:
            base_val = base_props[prop]
            prot_val = prot_props[prop]
            change = prot_val - base_val
            print(f"{prop:<15} {base_val:<12.2f} {prot_val:<14.2f} {change:<10.2f}")
        
        # Save detailed results
        csv_path = out_dir / 'vina_protonation_results.csv'
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Metric', 'Base', 'Protonated', 'Difference', 'Units'])
            writer.writerow(['Binding Affinity (Vina)', f"{results['base_energy']:.2f}", 
                           f"{results['protonated_energy']:.2f}", 
                           f"{results['protonation_difference']:.2f}", 'kcal/mol'])
            writer.writerow(['Molecular Weight', f"{base_props['mw']:.2f}", 
                           f"{prot_props['mw']:.2f}", f"{prot_props['mw']-base_props['mw']:.2f}", 'Da'])
            writer.writerow(['LogP', f"{base_props['logp']:.2f}", 
                           f"{prot_props['logp']:.2f}", f"{prot_props['logp']-base_props['logp']:.2f}", ''])
            writer.writerow(['H-Bond Donors', f"{base_props['hbd']:.0f}", 
                           f"{prot_props['hbd']:.0f}", f"{prot_props['hbd']-base_props['hbd']:.0f}", ''])
            writer.writerow(['H-Bond Acceptors', f"{base_props['hba']:.0f}", 
                           f"{prot_props['hba']:.0f}", f"{prot_props['hba']-base_props['hba']:.0f}", ''])
            writer.writerow(['Formal Charge', f"{base_props['charge']:.0f}", 
                           f"{prot_props['charge']:.0f}", f"{prot_props['charge']-base_props['charge']:.0f}", ''])
        
        print(f"\nDetailed results saved to: {csv_path}")
        
        # Interpretation
        print("\n" + "="*60)
        print("INTERPRETATION:")
        if results['protonation_difference'] < -1.0:
            print("✅ Protonated form is SIGNIFICANTLY MORE FAVORABLE for binding")
            print("   → Protonation enhances binding affinity")
        elif results['protonation_difference'] > 1.0:
            print("✅ Base form is SIGNIFICANTLY MORE FAVORABLE for binding") 
            print("   → Deprotonation enhances binding affinity")
        else:
            print("⚠️  Both forms have SIMILAR binding affinities")
            print("   → Protonation state has minimal effect")
        
        if abs(results['protonation_difference']) > 2.0:
            print(f"\n💡 MAJOR DIFFERENCE: {abs(results['protonation_difference']):.1f} kcal/mol")
            print("   → This is a CRITICAL pH-dependent binding effect!")
        elif abs(results['protonation_difference']) > 1.0:
            print(f"\n💡 Moderate difference: {abs(results['protonation_difference']):.1f} kcal/mol")
            print("   → pH effects are relevant for drug development")
        else:
            print(f"\n💡 Small difference: {abs(results['protonation_difference']):.1f} kcal/mol")
            print("   → Protonation state has minor effect")
        
        print("="*60)
        
        # Save SMILES for reference
        smiles_file = out_dir / 'structures.txt'
        with open(smiles_file, 'w') as f:
            f.write(f"Base form SMILES:\n{results['base_smiles']}\n\n")
            f.write(f"Protonated form SMILES:\n{results['prot_smiles']}\n")
        
        print(f"Molecular structures saved to: {smiles_file}")
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
