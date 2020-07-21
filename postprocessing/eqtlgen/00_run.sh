#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l mem=32gb
#PBS -e 00_run.err
#PBS -o 00_run.out


if [[ ! -z ${PBS_O_WORKDIR} ]]
then
  source ~/.bash_profile
  source ~/.bashrc
  conda activate mixqtl
  cd ${PBS_O_WORKDIR}
fi

eqtlpos=/gpfs/data/im-lab/nas40t2/yanyul/eQTLGen/cis-eQTL_significant_20181017_with_GTExV8ID.txt.gz
eqtlneg=/gpfs/data/im-lab/nas40t2/yanyul/eQTLGen/cis-eQTLs_full_20180905_with_GTExV8ID.txt.gz

codepath=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/mixqtl-pipeline
genelist=/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine/Whole_Blood/genes-passed-qc.txt.gz
outdir=samples_from_eqtlgen

mkdir -p $outdir

for seed in `seq 1 10`
do
  
  echo "Seed = $seed: positive samples"
  bash sample_positive_from_eqtlgen.sh \
    $codepath \
    $genelist \
    100000 \
    $seed \
    $eqtlpos \
    $outdir
  
  echo "Seed = $seed: negative samples"
  bash sample_negative_from_eqtlgen.sh \
    $codepath \
    $genelist \
    100000 \
    $seed \
    $eqtlneg \
    $outdir

done
