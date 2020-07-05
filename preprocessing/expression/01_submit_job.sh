mkdir -p logs

for row in `cat tissue_metainfo.csv`
do
  IFS=', ' read -a array <<< "'$row'"
  tissue="${array[0]}"
  tissuestr="${array[1]}"
  echo qsub -N expr-$tissue -v TISSUE=$tissue,TISSUESTR="'$tissuestr'" run.qsub
done