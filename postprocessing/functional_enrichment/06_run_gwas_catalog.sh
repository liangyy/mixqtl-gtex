#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l mem=16gb
#PBS -e logs/run_gwas_catalog_${tissue}.err
#PBS -o logs/run_gwas_catalog_${tissue}.out


source ~/.bash_profile
source ~/.bashrc


# ARGS1: tissue name
if [[ ! -z $1 ]]
then
  tissue=$1
fi
conda activate mixqtl

if [[ ! -z $PBS_O_WORKDIR ]]
then
  cd $PBS_O_WORKDIR
fi

gwas_catalog='/gpfs/data/im-lab/nas40t2/yanyul/mv_from_scratch/repo_new/rotation-at-imlab/analysis/annotate_gwas_catalog/output/GWAS_SNP_in_GTEx.with_key_and_info.GTExID2rsID__GWAS_catalog_in_rsID.txt'
mixqtl='/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl/{tissue}/mixqtl.{tissue}_GTEx_eGene.cis_qtl_pairs.mixQTL.chr{chr_num}.parquet'
top_qtl='/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/{tissue}.v8.egenes.txt.gz'
strong_gene='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine/{tissue}/genes-passed-qc.txt.gz'
mixfine_pip='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_pip.mixfine.chr{chr_num}.parquet'
nefine_pip='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_pip.nefine.chr{chr_num}.parquet'
mixfine_cs='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_cs.mixfine.chr{chr_num}.parquet'
nefine_cs='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_cs.nefine.chr{chr_num}.parquet'

cache_dir='/scratch/t.cri.yliang/mixQTL-GTExV8/cache/'

mkdir -p $cache_dir
if [[ ! -f enrichment/gwas_catalog_enrichment_$tissue.csv ]]
then
  python functional_enrichment_gwas_catalog.py \
    --tissue $tissue \
    --gwas_catalog_table $gwas_catalog 0 \
    --mixqtl $mixqtl \
    --top_qtl $top_qtl \
    --strong_gene $strong_gene \
    --mixfine_pip $mixfine_pip \
    --nefine_pip $nefine_pip \
    --mixfine_cs $mixfine_cs \
    --nefine_cs $nefine_cs \
    --cache_dir $cache_dir \
    --output enrichment/gwas_catalog_enrichment_$tissue.csv
fi

