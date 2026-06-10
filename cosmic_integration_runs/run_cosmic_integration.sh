#!/bin/bash
#SBATCH --array=0-3
#SBATCH --job-name=CosmicIntegration
#SBATCH --account=sdp140
#SBATCH --output=/home/msantiago/batch_files/outputs/CosmicIntegration_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=0-02:00:0

source activate /home/msantiago/.conda/envs/conda_env

export COMPAS_ROOT_DIR=/home/msantiago/COMPAS

cd /home/msantiago/COMPAS

FILES=(N5e6_MassiveWDWD_NSNS_CE_alpha025 N5e6_MassiveWDWD_NSNS_CE_alpha05 N5e6_MassiveWDWD_NSNS_CE_alpha075 N5e6_MassiveWDWD_NSNS_CE_alpha2)

INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}

python -m compas_python_utils.cosmic_integration.FastCosmicIntegration --path /expanse/lustre/scratch/msantiago/temp_project/COMPAS_DATA/$INPUT/MainRun/COMPAS_Output_wWeights.h5 --dco_type all --weight mixture_weight --maxzdet 8 --zstep 0.01 --sens O3 --m1min 1 --m1max 150 --fbin 0.7 --mu0 0.025 --muz -0.049 --sigma0 1.129 --sigmaz 0.048 --alpha -1.79 --aSF 0.02 --bSF 1.48 --cSF 4.44 --dSF 5.90
