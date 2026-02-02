import argparse
import subprocess
from pathlib import Path
from typing import Tuple, Optional

from meeko import MoleculePreparation
from meeko import PDBQTMolecule


class VinaRunner:
    def __init__(self, vina_exe: str = "C:\\tools\\vina\\vina.exe", exhaustiveness: int = 16, num_modes: int = 5, energy_range: int = 6):
        self.vina_exe = vina_exe
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.energy_range = energy_range

    def prepare_ligand(self, ligand_sdf: str, out_pdbqt: str):
        from rdkit import Chem
        mol = Chem.SDMolSupplier(ligand_sdf)[0]
        mol = Chem.AddHs(mol)  # Add explicit hydrogens
        prep = MoleculePreparation()
        mol_prep = prep.prepare(mol)
        # mol_prep is a list of prepared molecules, take the first one
        prepared_mol = mol_prep[0]
        # Write directly using the preparation result
        prep.write_pdbqt_file(out_pdbqt)
        return out_pdbqt

    def run_vina(self, protein_pdbqt: str, ligand_pdbqt: str, center: Tuple[float, float, float],
                 size: Tuple[float, float, float], out_prefix: str, extra_args: Optional[list] = None) -> Tuple[str, str]:
        cmd = [
            self.vina_exe,
            "--receptor", protein_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--energy_range", str(self.energy_range),
            "--out", f"{out_prefix}_poses.pdbqt",
            "--log", f"{out_prefix}_log.txt"
        ]
        if extra_args:
            cmd.extend(extra_args)
        subprocess.run(cmd, check=True)
        return f"{out_prefix}_poses.pdbqt", f"{out_prefix}_log.txt"

    @staticmethod
    def parse_vina_log(log_path: str) -> Optional[float]:
        """Return best-mode binding energy (kcal/mol) from a Vina log."""
        best_energy = None
        with open(log_path, "r") as f:
            for line in f:
                if line.strip().startswith("REMARK VINA RESULT"):
                    parts = line.strip().split()
                    if len(parts) >= 5:
                        try:
                            best_energy = float(parts[3])
                            break
                        except ValueError:
                            continue
        return best_energy


def main():
    parser = argparse.ArgumentParser(description="Run AutoDock Vina on prepared pdbqt files")
    parser.add_argument("--protein_pdbqt", required=True)
    parser.add_argument("--ligand_sdf", required=True)
    parser.add_argument("--center", nargs=3, type=float, required=True, help="grid center x y z")
    parser.add_argument("--size", nargs=3, type=float, required=True, help="grid size x y z")
    parser.add_argument("--out", default="results/vina/run1")
    parser.add_argument("--vina_exe", default="C:\\tools\\vina\\vina.exe")
    args = parser.parse_args()

    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    runner = VinaRunner(vina_exe=args.vina_exe)
    ligand_pdbqt = runner.prepare_ligand(args.ligand_sdf, f"{args.out}_ligand.pdbqt")
    runner.run_vina(args.protein_pdbqt, ligand_pdbqt, tuple(args.center), tuple(args.size), args.out)
    print(f"Vina finished. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
