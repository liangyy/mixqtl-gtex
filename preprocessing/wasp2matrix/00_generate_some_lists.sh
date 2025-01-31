OUTDIR=/scratch/t.cri.yliang/mixQTL-GTExV8/misc_data

cd $OUTDIR

zcat /gpfs/data/im-lab/nas40t2/yanyul/ASE/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_matrix.gw_phased.txt.gz |head -n 1|cut -f 5- | tr '\t' '\n' > sample_list_phaser_all.txt
zcat /gpfs/data/im-lab/nas40t2/yanyul/ASE/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_matrix.gw_phased.txt.gz |tail -n +2|cut -f 2 > gene_list_phaser_all.txt
zcat /gpfs/data/im-lab/nas40t2/yanyul/ASE/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_matrix.gw_phased.txt.gz |head -n 1|cut -f 5- |tr '\t' '\n'|awk -F"-" '{print $1"-"$2}' | sort | uniq > indiv_list_phaser_all.txt
