@echo off
REM Simple batch script to run protonation-aware scoring
REM Usage: run_protonation_analysis.bat

echo ========================================
echo PROTONATION-AWARE MOLECULAR DOCKING
echo ========================================
echo.

REM Activate conda environment
call conda activate diffdock-cpu

REM Check if arguments provided
if "%~1"=="" (
    echo Using default example files...
    set PROTEIN=data/1a0q/1a0q_protein_processed.pdbqt
    set LIGAND=data/1a0q/1a0q_ligand.sdf
    set CENTER=0 0 0
    set SIZE=20 20 20
    set OUTPUT=results/protonation_analysis
) else (
    echo Please provide: protein.pdbqt ligand.sdf center_x center_y center_z size_x size_y size_z
    echo Example: run_protonation_analysis.bat protein.pdbqt ligand.sdf 10.5 15.2 8.7 20 20 20
    exit /b 1
)

echo Running protonation analysis...
echo Protein: %PROTEIN%
echo Ligand: %LIGAND%
echo Binding site center: %CENTER%
echo Binding site size: %SIZE%
echo Output directory: %OUTPUT%
echo.

python local_runner/standalone_protonation_scorer.py --protein_pdbqt %PROTEIN% --ligand %LIGAND% --center %CENTER% --size %SIZE% --out_dir %OUTPUT%

echo.
echo ========================================
echo Analysis complete!
echo Results saved to: %OUTPUT%
echo ========================================
pause