# ARGS1: code path to mixqtl-pipeline: https://github.com/liangyy/mixqtl-pipeline
# ARGS2: gene list
# ARGS3: number of sample to extract
# ARGS4: random seed
# ARGS5: eQTLGen cis-QTL full summary statistics
# ARGS6: outdir

CODEPATH=$1
GENELIST=$2
nsample=$3
MYSEED=$4
eqtlgen=$5

prefix=$6/eqtlgen_neg

if [ ! -f $prefix.pval_gt_0.5.txt.gz ]; then
  zcat $eqtlgen | awk '{split($1,b,"E"); if($1>0.5 && b[2]=="") print $0}' | gzip > $prefix.pval_gt_0.5.txt.gz
fi

if [ ! -f $prefix.pval_gt_0.5-with-gene-qc.txt.gz ];then
  awk 'FNR==NR{a[$1]=1;next}{ if($8 in a) print $0}' <(zcat $GENELIST) <(zcat $prefix.pval_gt_0.5.txt.gz) | gzip > $prefix.pval_gt_0.5-with-gene-qc.txt.gz
fi

nrow=`zcat $prefix.pval_gt_0.5-with-gene-qc.txt.gz | wc -l`

Rscript $CODEPATH/misc_data/subsample_idx.R --number_total $nrow --number_sub $nsample --output $prefix.idx.subsample$size-with-gene-qc.seed$MYSEED.txt --seed $MYSEED

awk 'FNR==NR{a[$1]=1;next}{ if(FNR in a) print $0}' $prefix.idx.subsample$size-with-gene-qc.seed$MYSEED.txt <(zcat $prefix.pval_gt_0.5-with-gene-qc.txt.gz) | gzip > $prefix.subsample$size-with-gene-qc.seed$MYSEED.txt.gz

