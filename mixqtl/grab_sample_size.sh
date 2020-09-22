mydir=`pwd`
cd /scratch/t.cri.yliang/mixQTL-GTExV8/preprocessing/prepare_data 
for i in `ls`; do num=`zcat $i/total_count.bed.gz |head -n 1| awk '{print NF-4}'`; echo $i $num; done > $mydir/sample_size_from_pipeline.txt
