#! /bin/bash
#SBATCH -A INF18_teongrav
#SBATCH --partition bdw_usr_prod
##SBATCH --partition knl_usr_prod
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -c 35
#SBATCH -J mkn
#SBATCH -o /marconi_scratch/userexternal/aperego0/macrokilonova_bayes/test/mkn.out
#SBATCH -e /marconi_scratch/userexternal/aperego0/macrokilonova_bayes/test/mkn.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alvin.0983@gmail.com

# here one has to load everything necessary to run the code on the shell
module purge
module load profile/base
module load env-bdw
#module load env-knl  #one of the two...
module load python/3.5.2
module load intel/pe-xe-2017--binary
module load blas
module load mkl/2017--binary                 # not sure it is really needed
module load numpy/1.11.2--python--3.5.2      # not sure it is really needed
source ~/my_venv/bin/activate   # here I have scipy, matplotlib, cython
                                # I have also locally installed cpnest: copy the folder in the home and then compile it with
                                # LDSHARED="icc -shared" CC=icc python setup.py install
echo "Preparing:"
set -x                          # Output commands
set -e                          # Abort on errors
cd /marconi_scratch/userexternal/aperego0/macrokilonova_bayes/mk_source/source/

echo "Checking:"
pwd
hostname
date

export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
export KMP_AFFINITY=disabled
echo "Starting:"
srun -n 1 -c 35 python -u macrokilonova_model.py -o /marconi_scratch/userexternal/aperego0/macrokilonova_bayes/test/ -t 35 -f 1 --nlive 1000 --maxmcmc 100 --poolsize 100

echo "Stopping:"
date
