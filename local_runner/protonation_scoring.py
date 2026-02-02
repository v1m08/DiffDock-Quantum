import argparse
import csv
import tempfile
from pathlib import Path
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import MolStandardize

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from diffdock_run import run_diffdock
from vina_pipeline import VinaRunner


def generate_base_and_protonated(mol: Chem.Mol) -> Tuple[Chem.Mol, Chem.Mol]:
    parent = MolStandardize.FragmentParent(mol)
    uncharger = MolStandardize.Uncharger()
    base = uncharger.uncharge(parent)
    protonated = MolStandardize.Reionize()(MolStandardize.Cleanup(parent))
    base = Chem.AddHs(base, addCoords=True)
    protonated = Chem.AddHs(protonated, addCoords=True)
    return base, protonated


def write_temp_sdf(mol: Chem.Mol, path: Path):
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.flush()
    writer.close()


def collect_pose_sdfs(diffdock_dir: Path) -> List[Path]:
    sdf_paths = sorted(diffdock_dir.glob("**/rank*.sdf"))
    if not sdf_paths:
        raise FileNotFoundError(f"No SDF poses found under {diffdock_dir}; run DiffDock first or check out_dir.")
    return sdf_paths


def score_pose_with_vina(vina_runner: VinaRunner, protein_pdbqt: str, pose_path: Path,
                         center: Tuple[float, float, float], size: Tuple[float, float, float],
                         out_prefix: Path) -> Tuple[float, float, float]:
    mol = Chem.SDMolSupplier(str(pose_path), removeHs=False)[0]
    base, protonated = generate_base_and_protonated(mol)

    base_sdf = out_prefix.with_suffix("_base.sdf")
    prot_sdf = out_prefix.with_suffix("_prot.sdf")
    write_temp_sdf(base, base_sdf)
    write_temp_sdf(protonated, prot_sdf)

    base_ligand_pdbqt = str(out_prefix.with_suffix("_base.pdbqt"))
    prot_ligand_pdbqt = str(out_prefix.with_suffix("_prot.pdbqt"))

    vina_runner.prepare_ligand(str(base_sdf), base_ligand_pdbqt)
    vina_runner.prepare_ligand(str(prot_sdf), prot_ligand_pdbqt)

    base_pose, base_log = vina_runner.run_vina(protein_pdbqt, base_ligand_pdbqt, center, size, str(out_prefix) + "_base")
    prot_pose, prot_log = vina_runner.run_vina(protein_pdbqt, prot_ligand_pdbqt, center, size, str(out_prefix) + "_prot")

    base_energy = vina_runner.parse_vina_log(base_log)
    prot_energy = vina_runner.parse_vina_log(prot_log)
    if base_energy is None or prot_energy is None:
        raise RuntimeError(f"Failed to parse Vina energies for {pose_path}")
    delta = prot_energy - base_energy
    return base_energy, prot_energy, delta


def aggregate_score(prot_energy: float, delta: float, delta_weight: float) -> float:
    return prot_energy + delta_weight * delta


def main():
    parser = argparse.ArgumentParser(description="DiffDock poses re-ranked by protonation-aware Vina energies")
    parser.add_argument("--protein_pdb", required=False, default="data/1a0q/1a0q_protein_processed.pdb")
    parser.add_argument("--protein_pdbqt", required=True, help="Receptor PDBQT for Vina")
    parser.add_argument("--ligand", required=False, default="data/1a0q/1a0q_ligand.sdf")
    parser.add_argument("--config", default="default_inference_args.yaml")
    parser.add_argument("--out_dir", default="results/protonation_scoring")
    parser.add_argument("--center", nargs=3, type=float, required=True, help="Grid center x y z")
    parser.add_argument("--size", nargs=3, type=float, required=True, help="Grid size x y z")
    parser.add_argument("--samples_per_complex", type=int, default=10)
    parser.add_argument("--batch_size", type=int, default=1)
    parser.add_argument("--vina_exe", default="C:\\tools\\vina\\vina.exe")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--num_modes", type=int, default=5)
    parser.add_argument("--energy_range", type=int, default=6)
    parser.add_argument("--delta_weight", type=float, default=1.0, help="Penalty weight for binding energy change after protonation")
    parser.add_argument("--skip_diffdock", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    diffdock_dir = out_dir / "diffdock_poses"

    if not args.skip_diffdock:
        run_diffdock(args.protein_pdb, args.ligand, str(diffdock_dir), args.config, args.samples_per_complex, args.batch_size)

    pose_paths = collect_pose_sdfs(diffdock_dir)

    runner = VinaRunner(vina_exe=args.vina_exe, exhaustiveness=args.exhaustiveness,
                        num_modes=args.num_modes, energy_range=args.energy_range)

    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for idx, pose_path in enumerate(pose_paths):
            prefix = tmp_path / f"pose{idx+1}"
            base_e, prot_e, delta = score_pose_with_vina(runner, args.protein_pdbqt, pose_path,
                                                         tuple(args.center), tuple(args.size), prefix)
            score = aggregate_score(prot_e, delta, args.delta_weight)
            results.append({
                "pose": str(pose_path),
                "base_energy": base_e,
                "protonated_energy": prot_e,
                "delta": delta,
                "aggregate_score": score
            })

    results = sorted(results, key=lambda r: r["aggregate_score"])
    csv_path = out_dir / "protonation_scores.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pose", "base_energy", "protonated_energy", "delta", "aggregate_score"])
        writer.writeheader()
        writer.writerows(results)

    print(f"Protonation-aware ranking complete. {len(results)} poses scored. Results: {csv_path}")


if __name__ == "__main__":
    main()
