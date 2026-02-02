import argparse
import subprocess
from pathlib import Path

from vina_pipeline import VinaRunner


def run_pdb2pqr(input_pdb: str, out_pqr: str, ph: float = 7.0, ff: str = "AMBER"):
    cmd = [
        "C:\\tools\\pdb2pqr\\pdb2pqr.exe",
        f"--ff={ff}",
        f"--with-ph={ph}",
        input_pdb,
        out_pqr
    ]
    subprocess.run(cmd, check=True)
    return out_pqr


def main():
    parser = argparse.ArgumentParser(description="Protonate with PDB2PQR then run Vina")
    parser.add_argument("--protein_pdb", required=True)
    parser.add_argument("--ligand_sdf", required=True)
    parser.add_argument("--center", nargs=3, type=float, required=True)
    parser.add_argument("--size", nargs=3, type=float, required=True)
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--forcefield", default="AMBER")
    parser.add_argument("--vina_exe", default="C:\\tools\\vina\\vina.exe")
    parser.add_argument("--out", default="results/pdb2pqr_vina/run1")
    args = parser.parse_args()

    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    pqr_path = f"{args.out}_protein.pqr"
    run_pdb2pqr(args.protein_pdb, pqr_path, ph=args.ph, ff=args.forcefield)

    # AutoDock Vina needs pdbqt; assume user provides protein pdbqt externally or converts separately.
    # Here we expect an existing protein pdbqt alongside the input pdb name with .pdbqt extension.
    protein_pdbqt_guess = Path(args.protein_pdb).with_suffix(".pdbqt")
    if not protein_pdbqt_guess.exists():
        raise FileNotFoundError(f"Expected protein pdbqt at {protein_pdbqt_guess}; please prepare it with AutoDockTools or OBabel.")

    vina_runner = VinaRunner(vina_exe=args.vina_exe)
    ligand_pdbqt = vina_runner.prepare_ligand(args.ligand_sdf, f"{args.out}_ligand.pdbqt")
    vina_runner.run_vina(str(protein_pdbqt_guess), ligand_pdbqt, tuple(args.center), tuple(args.size), args.out)

    print(f"PDB2PQR + Vina finished. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
