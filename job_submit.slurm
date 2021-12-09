#!/bin/bash
#SBATCH --job-name=v1adm    # create a short name for your job
#SBATCH --account=rrg-murray-ac   # Account resource allocation
#SBATCH --nodes=2               # node count
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --output=/scratch/m/murray/sr39/debug/m12ires450000_v1/output/%x-%j.out
#SBATCH --time=23:59:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end, and fail                                                            
#SBATCH --mail-user=sandiproy@princeton.edu                                                                                                            

module purge
#source /home/m/murray/sr39/gizmo_mod.sh
module load NiaEnv/2018a intel/2018.2 gsl/2.4 intelmpi/2018.2 hdf5-mpi/1.10.2

MYPATH=/home/m/murray/sr39/gizmo-debug
#srun $MYPATH/GIZMO /tigress/sr39/gas_disk_test_Feb2020/gas_spin_analhalo_w1/gas_spin_w1.params 1
#srun $MYPATH/GIZMO_ncrit /tigress/sr39/gizmo_adm_test_v2/starcriteria_test/ncrit_100.params
srun $MYPATH/GIZMO_debug1 /scratch/m/murray/sr39/debug/m12ires450000_v1/gizmo_parameters.txt
