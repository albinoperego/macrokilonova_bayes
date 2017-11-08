#!/bin/bash
#SBATCH --account=s667
#SBATCH --ntasks=64
#SBATCH --ntasks-per-core=2
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --constraint=mc
#SBATCH --partition=normal
#SBATCH --job-name="lc1"
##SBATCH --mail-type=ALL
##SBATCH --mail-user=albino.perego@unibas.ch
#SBATCH --output=/scratch/snx1600/albino/mk_source/lc.out
#SBATCH --error=/scratch/snx1600/albino/mk_source/lc.err

# source ~/venv-2.7/bin/activate

module load PyExtensions

#set -ex

##cd $SCRATCH
export WRK=/scratch/snx1600/albino/mk_source
##export NLSPATH=/users/albino/pathscale/%N.cat
#export CORE_MAX_PROCESSES=0
#export CORE_ACTION_FIRST=KILL
#export CORE_ACTION_OTHER=KILL
##export F90_DUMP_MAP=YES
#export F90_BOUNDS_CHECK_ABORT=YES
#unset IOBUF_PARAMS    #SS:might be left away for #proc<64
#export MPICH_RANK_REORDER_METHOD=1   #SS:5.11.07
#export MPICH_FAST_MEMCPY=1           #SS:5.11.07
#export MPICH_UNEX_BUFFER_SIZE=260M
#export FILENV=$HOME/.filenv   #SS:speed up writing to files
##assign -F f77.bufsize=4194304 g:su
##export PSC_OMP_AFFINITY=TRUE  #SS: 20.5.09 multi threading pathscale
#
#export OMP_NUM_THREADS=1

cd $WRK  
pwd

python lightcurve_vary_k.py 64

exit
