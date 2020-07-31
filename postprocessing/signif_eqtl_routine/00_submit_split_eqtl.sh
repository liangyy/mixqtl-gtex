# ARGS1: OUTDIR

OUTDIR=$1

mkdir -p logs/
mkdir -p $OUTDIR


for tissue in `cat ../preprocessing/expression/tissue_metainfo.csv|cut -f 1 -d,`
do
  qsub -v TISSUE=$tissue,OUTDIR=$OUTDIR run_split_eqtl.qsub
done 