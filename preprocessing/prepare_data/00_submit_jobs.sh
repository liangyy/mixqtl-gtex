mkdir -p logs

for tt in `cat ../tissue_metainfo.csv | cut -f 1 -d,`
do
  if [[ -f /scratch/t.cri.yliang/mixQTL-GTExV8/preprocessing/peer/$tt/covariate-combined.txt ]]
  then
    if [[ ! -f /scratch/t.cri.yliang/mixQTL-GTExV8/preprocessing/prepare_data/$tt/covariates.txt.gz ]]
    then
      echo qsub -N prepare-$tt -v TISSUE=$tt run.qsub
    fi
  fi
done
