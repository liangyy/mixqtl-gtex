# submit the list of jobs to upload mixQTL sum stats to zenodo

mkdir -p zenodo_meta

for tissue in `cat ../preprocessing/expression/tissue_metainfo.csv | cut -f 1 -d,`
do 
  cat zenodo_meta_pattern.yaml | sed "s#TISSUE-PLACEHOLDER#$tissue#g" > zenodo_meta/$tissue.yaml
done

# split tissues into 49 tissues per batch
cat ../preprocessing/expression/tissue_metainfo.csv | cut -f 1 -d, > zenodo_meta/tissues.txt
split -l 49 -d zenodo_meta/tissues.txt zenodo_meta/batch

for i in `ls zenodo_meta/batch*`
do
  mynum=`echo $i | sed 's#zenodo_meta/batch##g'`
  qsub -v BATCH=$mynum -N zenodo_$mynum run_upload.qsub
done
