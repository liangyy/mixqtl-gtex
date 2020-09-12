# ARGS1: if specify run on the corresponding tissue


mkdir -p logs
mkdir -p /scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl_cleanup/

if [[ ! -z $1 ]]
then

  tissue=$1
  qsub -N cleanup-$tissue -v TISSUE=$tissue,\
INPREFIX=/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl/$tissue/mixqtl.${tissue}_GTEx_eGene.cis_qtl_pairs.mixQTL.chr,\
INSUFFIX=.parquet,\
OUTPREFIX=/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl_cleanup/mixqtl.${tissue}.chr,\
TENSORQTL=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/tensorqtl,\
OUTSUFFIX=.parquet cleanup.qsub

else

  for tissue in `cat ../preprocessing/expression/tissue_metainfo.csv | cut -f 1 -d,`
  do
    
    qsub -N cleanup-$tissue -v TISSUE=$tissue,\
INPREFIX=/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl/$tissue/mixqtl.${tissue}_GTEx_eGene.cis_qtl_pairs.mixQTL.chr,\
INSUFFIX=.parquet,\
OUTPREFIX=/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl_cleanup/mixqtl.${tissue}.chr,\
TENSORQTL=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/tensorqtl,\
OUTSUFFIX=.parquet cleanup.qsub
      
  done
fi

