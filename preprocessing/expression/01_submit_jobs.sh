mkdir -p logs

for ii in `seq 1 49`
do
  row=`cat tissue_metainfo.csv | head -n $ii | tail -n 1`
  IFS=',' read -a array <<< "$row"
  tissue="${array[0]}"
  tissuestr="${array[1]}"
  echo qsub -N expr-$tissue -v TISSUE=$tissue,TISSUESTR="$tissuestr" run.qsub
done
