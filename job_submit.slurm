#!/bin/bash
#SBATCH --job-name=cosmo_test     # create a short name for your job
#SBATCH --nodes=2               # node count
# SBATCH --exclusive              # fully dedicate single node 
#SBATCH -p physics
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=100G               # memory for entire job
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=/scratch/gpfs/sr39/m12ires450000_test/output/%x-%j.out # output file
#SBATCH --mail-type=all          # send email on job start, end, and fail                                                            
#SBATCH --mail-user=sandiproy@princeton.edu                                                                                                            

module purge
module load openmpi/gcc/1.10.2/64 gsl/2.4 fftw/gcc/openmpi-1.10.2/3.3.4 hdf5/gcc/1.8.16

MYPATH=/home/sr39/gizmo
#srun $MYPATH/GIZMO /tigress/sr39/gas_disk_test_Feb2020/gas_spin_analhalo_w1/gas_spin_w1.params 1
#srun $MYPATH/GIZMO_ncrit /tigress/sr39/gizmo_adm_test_v2/starcriteria_test/ncrit_100.params
srun $MYPATH/GIZMO /scratch/gpfs/sr39/m12ires450000_test/gizmo_parameters.txt 2
