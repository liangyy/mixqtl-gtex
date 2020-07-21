#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l mem=64gb
#PBS -e 01_extract_qtl_sum_stats.err
#PBS -o 01_extract_qtl_sum_stats.out

source ~/.bash_profile
source ~/.bashrc
conda activate mixqtl
cd ${PBS_O_WORKDIR}

python extract_from_mixqtl.py 
