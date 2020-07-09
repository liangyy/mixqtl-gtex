# ARGS1: number of cores to use per job
# ARGS2: if specify run on the corresponding tissue

ncore=$1

mkdir -p logs

if [[ ! -z $2 ]]
then
  tissue=$2
  qsub -N mixqtl-$tissue -l nodes=1:ppn=$ncore -v TISSUE=$tissue,NCORE=$ncore run.qsub
else

  for tissue in `cat ../preprocessing/expression/tissue_metainfo.csv | cut -f 1 -d,`
  do
    if [[ -f /scratch/t.cri.yliang/mixQTL-GTExV8/preprocessing/prepare_data/$tissue/gtex_v8_library_size.txt.gz ]]
    then
      have=0 
      for kk in `seq 1 22`
      do
        if [[ -f /scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl/$tissue/mixqtl.${tissue}_GTEx_eGene.cis_qtl_pairs.mixQTL.chr$kk.parquet ]]
        then
          have=1
        fi
      done
      if [[ $have == 0 ]]
      then
        qsub -N mixqtl-$tissue -l nodes=1:ppn=$ncore -v TISSUE=$tissue,NCORE=$ncore run.qsub
      fi
    fi
  done
fi

