# ARGS1: OUTDIR

OUTDIR=$1

mkdir -p logs/
mkdir -p $OUTDIR


for tissue in `cat ../../preprocessing/expression/tissue_metainfo.csv|cut -f 1 -d,`
do
  if [[ -f logs/split_$tissue.out ]]
  then
    e=`cat logs/split_$tissue.out|grep Exit | tail -n 1|grep 1 | wc -l`
    if [[ $e == 1 ]]
    then
      echo qsub -v TISSUE=$tissue,OUTDIR=$OUTDIR run_split_eqtl.qsub 
    fi
  else
    echo qsub -v TISSUE=$tissue,OUTDIR=$OUTDIR run_split_eqtl.qsub
  fi
done 
