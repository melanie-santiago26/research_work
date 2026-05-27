#!/bin/bash
#SBATCH --job-name=CosmicIntegration
#SBATCH --output=CosmicIntegration_%j.out
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --account=sdp140
#SBATCH --export=ALL

module purge
module load cpu

cd /expanse/lustre/scratch/msantiago/temp_project/COMPAS_DATA/N5e6_MassiveWDWD_NSNS_CE_alpha025/

source /home/msantiago/.conda/env/bin/activate

export COMPAS_ROOT_DIR=/home/msantiago/COMPAS
export PYTHONPATH=$COMPAS_ROOT_DIR:$PYTHONPATH

python -m compas_python_utils.cosmic_integration.FastCosmicIntegration --path /expanse/lustre/scratch/msantiago/temp_project/COMPAS_DATA/N5e6_MassiveWDWD_NSNS_CE_alpha025/MainRun/COMPAS_Output_wWeights.h5 --dco_type NSNS --weight mixture_weight --maxzdet 8 --zstep 0.01 --sens O3 --m1min 5 --m1max 150 --fbin None --mu0 0.025 --muz -0.049 --sigma0 1.129 --sigmaz 0.048 --alpha -1.79 --aSF 0.02 --bSF 1.48 --cSF 4.44 --dSF 5.90
