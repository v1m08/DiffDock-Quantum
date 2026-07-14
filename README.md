# Quantum DiffDock

This project is a fork of
[DiffDock](https://github.com/gcorso/DiffDock) exploring physically
informed, quantum-inspired pipeline for structure-based virtual screening. This repo has very minimal research value, and was primarily a fun side-project exploring [a really cool paper on diffusion models in computational biology]([url](https://arxiv.org/abs/2210.01776)).

## Project contribution

This fork adds:

- protonation-state generation and comparison with RDKit;
- AutoDock Vina scoring and re-ranking of DiffDock poses;
- Windows and WSL/Linux runners for the protonation-aware workflow;
- setup, verification, and example scripts for local experiments; and
- a Kaggle notebook and documentation for reproducing the prototype.

The project investigates whether physically informed energy calculations can
complement geometric RMSD-based docking metrics. The evaluation
reported 30.7% of ligand poses below 2 Å RMSD and 72.4% of centroid predictions
below 2 Å on a partial 127-sample time-split evaluation. In individual poses,
the post-processing optimizer improved the reported Vina binding-affinity
metric by as much as 117%, but often moved poses away from their experimental
geometry. 

The proposed next step is to train physical binding-affinity and RMSD objectives
together, then evaluate a simulated-bifurcation/quantum-annealing component in
the learned pipeline. That quantum component is research direction described by
the project; this repository currently contains the classical protonation/Vina
prototype.

## Running the added workflow

For the shortest setup path, see [PROTONATION_README.md](PROTONATION_README.md).
WSL/Linux users can also follow [README_WSL.md](README_WSL.md). Generated model
caches and analysis outputs are intentionally ignored by Git.

## Repository lineage

This project was built from the upstream DiffDock repository and retains its
history, model code, license, and citations. 
