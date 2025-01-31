mkdir -p logs

npeerlist=/scratch/t.cri.yliang/mixQTL-GTExV8/misc_data/npeer_list.csv

for ii in `seq 1 49`
do
  row=`cat $npeerlist | head -n $ii | tail -n 1`
  IFS=',' read -a array <<< "$row"
  tissue="${array[0]}"
  nfactor="${array[1]}"
  if [[ ! -f ~/scratch/mixQTL-GTExV8/preprocessing/peer/$tissue/covariate-combined.txt ]]
  then
    qsub -N peer-$tissue -v TISSUE=$tissue,NFACTOR=$nfactor run.qsub
  fi
done
