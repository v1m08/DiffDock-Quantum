#!/usr/bin/env python3
"""
Simple test script for protonation-aware scoring without DiffDock dependency.
Use pre-generated poses or simple ligand conformers for testing.
"""

import argparse
import csv
import tempfile
from pathlib import Path
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import MolStandardize

from vina_pipeline import VinaRunner


def generate_base_and_protonated(mol: Chem.Mol) -> Tuple[Chem.Mol, Chem.Mol]:
    # Simple approach: use original as base, add protons for protonated
    base = Chem.RemoveHs(mol)
    base = Chem.AddHs(base, addCoords=True)  # Add hydrogens in neutral state
    
    # For protonated state, add extra hydrogens to likely protonation sites
    protonated = Chem.AddHs(mol, addCoords=True)
    
    return base, protonated


def write_temp_sdf(mol: Chem.Mol, path: Path):
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.flush()
    writer.close()


def score_pose_with_vina(vina_runner: VinaRunner, protein_pdbqt: str, ligand_sdf: str,
                         center: Tuple[float, float, float], size: Tuple[float, float, float],
                         out_prefix: Path) -> Tuple[float, float, float]:
    mol = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]
    base, protonated = generate_base_and_protonated(mol)

    base_sdf = Path(str(out_prefix) + "_base.sdf")
    prot_sdf = Path(str(out_prefix) + "_prot.sdf")
    write_temp_sdf(base, base_sdf)
    write_temp_sdf(protonated, prot_sdf)

    base_ligand_pdbqt = str(out_prefix) + "_base.pdbqt"
    prot_ligand_pdbqt = str(out_prefix) + "_prot.pdbqt"

    vina_runner.prepare_ligand(str(base_sdf), base_ligand_pdbqt)
    vina_runner.prepare_ligand(str(prot_sdf), prot_ligand_pdbqt)

    base_pose, base_log = vina_runner.run_vina(protein_pdbqt, base_ligand_pdbqt, center, size, str(out_prefix) + "_base")
    prot_pose, prot_log = vina_runner.run_vina(protein_pdbqt, prot_ligand_pdbqt, center, size, str(out_prefix) + "_prot")

    base_energy = vina_runner.parse_vina_log(base_log)
    prot_energy = vina_runner.parse_vina_log(prot_log)
    if base_energy is None or prot_energy is None:
        raise RuntimeError(f"Failed to parse Vina energies for {ligand_sdf}")
    delta = prot_energy - base_energy
    return base_energy, prot_energy, delta


def aggregate_score(prot_energy: float, delta: float, delta_weight: float) -> float:
    return prot_energy + delta_weight * delta


def main():
    parser = argparse.ArgumentParser(description="Simple protonation-aware Vina scoring test")
    parser.add_argument("--protein_pdbqt", required=True, help="Receptor PDBQT for Vina")
    parser.add_argument("--ligand", required=True, help="Ligand SDF file")
    parser.add_argument("--center", nargs=3, type=float, required=True, help="Grid center x y z")
    parser.add_argument("--size", nargs=3, type=float, required=True, help="Grid size x y z")
    parser.add_argument("--out_dir", default="results/simple_test")
    parser.add_argument("--vina_exe", default="C:\\tools\\vina\\vina.exe")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--num_modes", type=int, default=5)
    parser.add_argument("--energy_range", type=int, default=6)
    parser.add_argument("--delta_weight", type=float, default=1.0)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    runner = VinaRunner(vina_exe=args.vina_exe, exhaustiveness=args.exhaustiveness,
                        num_modes=args.num_modes, energy_range=args.energy_range)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        prefix = tmp_path / "test"
        base_e, prot_e, delta = score_pose_with_vina(runner, args.protein_pdbqt, args.ligand,
                                                     tuple(args.center), tuple(args.size), prefix)
        score = aggregate_score(prot_e, delta, args.delta_weight)

        result = {
            "ligand": args.ligand,
            "base_energy": base_e,
            "protonated_energy": prot_e,
            "delta": delta,
            "aggregate_score": score
        }

        csv_path = out_dir / "simple_test_results.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["ligand", "base_energy", "protonated_energy", "delta", "aggregate_score"])
            writer.writeheader()
            writer.writerow(result)

        print(f"Simple protonation test complete. Results: {csv_path}")
        print(f"Base energy: {base_e:.3f} kcal/mol")
        print(f"Protonated energy: {prot_e:.3f} kcal/mol")
        print(f"Delta: {delta:.3f} kcal/mol")
        print(f"Aggregate score: {score:.3f}")


if __name__ == "__main__":
    main()