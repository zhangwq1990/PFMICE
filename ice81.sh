#!/bin/bash
# Example submission script for 'hello world' program
#SBATCH --partition=defq
# This specifies type of node job will use
#SBATCH --nodes=1
# This specifies job uses 1 node
#SBATCH --ntasks-per-node=25
# This specifies job only use 1 core on the node
#SBATCH --mem=4g
# This specifies maximum memory use will be 2 gigabytes
#SBATCH --time=20:00:00
# This specifies job will last no longer than 1 hour

#below use Linux commands, which will run on compute node
module purge
#module load openmpi-uon/gcc6.3.0/2.1.5
module load intel/compiler/2019/5.075
module load intel/impi/2018/
echo "Running on `hostname`"
cd ${SLURM_SUBMIT_DIR}

mpirun ./PFMICE > ./log
sleep 4
echo "Finished job now"

