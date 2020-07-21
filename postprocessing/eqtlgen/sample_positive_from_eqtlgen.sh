# ARGS1: code path to mixqtl-pipeline: https://github.com/liangyy/mixqtl-pipeline
# ARGS2: gene list
# ARGS3: number of sample to extract
# ARGS4: random seed
# ARGS5: eQTLGen significant cis-QTL summary statistics
# ARGS6: outdir

CODEPATH=$1
GENELIST=$2
nsample=$3
MYSEED=$4
eqtlgen=$5

prefix=$6/eqtlgen_pos

if [ ! -f $prefix.pos-with-gene-qc.txt.gz ];then
  awk 'FNR==NR{a[$1]=1;next}{ if($8 in a) print $0}' $genelist <(zcat $eqtlgen) | gzip > $prefix.pos-with-gene-qc.txt.gz
fi

nrow=`zcat $prefix.pos-with-gene-qc.txt.gz | wc -l` 

Rscript $CODEPATH/misc_data/subsample_idx.R --number_total $nrow --number_sub $nsample --seed $MYSEED --output $prefix.idx.subsample$size-with-gene-qc.seed$MYSEED.txt

awk 'FNR==NR{a[$1]=1;next}{ if(FNR in a) print $0}' $prefix.idx.subsample$size-with-gene-qc.seed$MYSEED.txt <(zcat $prefix.pos-with-gene-qc.txt.gz) | gzip > $prefix.subsample$size-with-gene-qc.seed$MYSEED.txt.gz
