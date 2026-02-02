import os
import argparse
import copy
from pathlib import Path
from argparse import Namespace

import torch
from tqdm import tqdm

# Reuse DiffDock modules
import sys
import os
# Add parent directory to path to import DiffDock modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from inference import get_parser as get_inference_parser, main as inference_main
from utils.inference_utils import InferenceDataset
from utils.diffusion_utils import t_to_sigma as t_to_sigma_compl, get_t_schedule
from utils.sampling import randomize_position, sampling
from utils.utils import get_model


def build_args(protein_path: str, ligand_path: str, out_dir: str, config: str = "default_inference_args.yaml",
               samples_per_complex: int = 5, batch_size: int = 1) -> Namespace:
    parser = get_inference_parser()
    args = parser.parse_args([])
    args.config = open(config, "r")
    args.protein_path = protein_path
    args.protein_sequence = None
    args.ligand_description = ligand_path
    args.protein_ligand_csv = None
    args.out_dir = out_dir
    args.samples_per_complex = samples_per_complex
    args.batch_size = batch_size
    return args


def run_diffdock(protein_path: str, ligand_path: str, out_dir: str,
                 config: str = "default_inference_args.yaml", samples_per_complex: int = 5,
                 batch_size: int = 1):
    """
    Run DiffDock inference on a single protein/ligand pair.
    Uses the official inference pipeline; downloads weights if missing.
    """
    args = build_args(protein_path, ligand_path, out_dir, config, samples_per_complex, batch_size)
    inference_main(args)


# Optional: minimal custom loop (reuses DiffDock model classes)
def run_diffdock_model_once(protein_path: str, ligand_path: str, out_dir: str,
                            config: str = "default_inference_args.yaml", samples_per_complex: int = 1,
                            batch_size: int = 1):
    """
    Example showing how to call the lower-level APIs without the full CLI loop.
    It mirrors inference.main but trims logging and multi-CSV handling.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    parser = get_inference_parser()
    args = parser.parse_args([])
    args.config = open(config, "r")
    args.protein_ligand_csv = None
    args.protein_path = protein_path
    args.protein_sequence = None
    args.ligand_description = ligand_path
    args.out_dir = out_dir
    args.samples_per_complex = samples_per_complex
    args.batch_size = batch_size

    # Merge yaml config
    import yaml
    if args.config:
        config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__
        for key, value in config_dict.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value

    os.makedirs(args.out_dir, exist_ok=True)
    with open(f'{args.model_dir}/model_parameters.yml') as f:
        score_model_args = Namespace(**yaml.full_load(f))

    t_to_sigma = lambda *ts: t_to_sigma_compl(*ts, args=score_model_args)

    # Build dataset
    complex_name_list = ["complex_0"]
    test_dataset = InferenceDataset(out_dir=args.out_dir, complex_names=complex_name_list,
                                    protein_files=[args.protein_path], ligand_descriptions=[args.ligand_description],
                                    protein_sequences=[None],
                                    lm_embeddings=True,
                                    receptor_radius=score_model_args.receptor_radius,
                                    remove_hs=score_model_args.remove_hs,
                                    c_alpha_max_neighbors=score_model_args.c_alpha_max_neighbors,
                                    all_atoms=score_model_args.all_atoms,
                                    atom_radius=score_model_args.atom_radius,
                                    atom_max_neighbors=score_model_args.atom_max_neighbors,
                                    knn_only_graph=False if not hasattr(score_model_args, 'not_knn_only_graph') else not score_model_args.not_knn_only_graph)

    test_loader = torch.utils.data.DataLoader(dataset=test_dataset, batch_size=1, shuffle=False)

    model = get_model(score_model_args, device, t_to_sigma=t_to_sigma, no_parallel=True)
    state_dict = torch.load(f'{args.model_dir}/{args.ckpt}', map_location=torch.device('cpu'))
    model.load_state_dict(state_dict, strict=True)
    model = model.to(device)
    model.eval()

    tr_schedule = get_t_schedule(inference_steps=args.inference_steps, sigma_schedule='expbeta')

    for idx, orig_complex_graph in tqdm(enumerate(test_loader)):
        if not orig_complex_graph.success[0]:
            continue
        data_list = [copy.deepcopy(orig_complex_graph) for _ in range(args.samples_per_complex)]
        randomize_position(data_list, score_model_args.no_torsion, False, score_model_args.tr_sigma_max,
                           initial_noise_std_proportion=args.initial_noise_std_proportion,
                           choose_residue=args.choose_residue)
        # run reverse diffusion
        data_list, confidence = sampling(data_list=data_list, model=model,
                                         inference_steps=args.actual_steps if args.actual_steps is not None else args.inference_steps,
                                         tr_schedule=tr_schedule, rot_schedule=tr_schedule, tor_schedule=tr_schedule,
                                         device=device, t_to_sigma=t_to_sigma, model_args=score_model_args,
                                         visualization_list=None, confidence_model=None,
                                         confidence_data_list=None, confidence_model_args=None,
                                         batch_size=args.batch_size, no_final_step_noise=args.no_final_step_noise,
                                         temp_sampling=[args.temp_sampling_tr, args.temp_sampling_rot, args.temp_sampling_tor],
                                         temp_psi=[args.temp_psi_tr, args.temp_psi_rot, args.temp_psi_tor],
                                         temp_sigma_data=[args.temp_sigma_data_tr, args.temp_sigma_data_rot, args.temp_sigma_data_tor])
        # Save poses
        lig = orig_complex_graph.mol[0]
        ligand_pos = [cg['ligand'].pos.cpu().numpy() + orig_complex_graph.original_center.cpu().numpy() for cg in data_list]
        out_dir_complex = Path(out_dir) / complex_name_list[idx]
        out_dir_complex.mkdir(parents=True, exist_ok=True)
        from datasets.process_mols import write_mol_with_coords
        for rank, pos in enumerate(ligand_pos):
            mol_pred = copy.deepcopy(lig)
            if score_model_args.remove_hs:
                from rdkit.Chem import RemoveAllHs
                mol_pred = RemoveAllHs(mol_pred)
            write_mol_with_coords(mol_pred, pos, out_dir_complex / f'rank{rank+1}.sdf')
    print('Finished custom loop; results in', out_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--protein', required=False, default='data/1a0q/1a0q_protein_processed.pdb')
    parser.add_argument('--ligand', required=False, default='data/1a0q/1a0q_ligand.sdf')
    parser.add_argument('--out_dir', required=False, default='results/local_diffdock')
    parser.add_argument('--config', default='default_inference_args.yaml')
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--batch', type=int, default=1)
    args = parser.parse_args()
    run_diffdock(args.protein, args.ligand, args.out_dir, args.config, args.samples, args.batch)
