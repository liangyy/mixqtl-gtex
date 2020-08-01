# ARGS1: OUTDIR

OUTDIR=$1

mkdir -p logs/
mkdir -p $OUTDIR

config='mixqtl eqtl trcqtl'

for tissue in `cat ../../preprocessing/expression/tissue_metainfo.csv|cut -f 1 -d,`
do
  for cc in $config
  do
    if [[ -f logs/signif_${cc}_${tissue}.log ]]
    then
      e=`cat logs/signif_${cc}_${tissue}.out | grep Exit | tail -n 1 | grep 1 | wc -l`
      if [[ $e != 1 ]]
      then
        qsub -v TISSUE=$tissue,OUTDIR=$OUTDIR,CONFIG=$cc run_signif_eqtl.qsub 
      fi
    else
      qsub -v TISSUE=$tissue,OUTDIR=$OUTDIR,CONFIG=$cc run_signif_eqtl.qsub
    fi
  done 
done
