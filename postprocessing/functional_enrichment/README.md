We have two sets of functional annotations. 

1. Regulatory element annotation of all GTEx variants curated by GTEx v8 working group
    - The preparation of the summary table is done by `functional_enrichment.py` 
    - `for i in `cat ../../mixfine/tissue_w_small_sample_size.txt`; do qsub -v tissue=$i 02_run.sh ; done`
2. ENCODE cCRE (merged by tissue) pre-processed by script at [here](https://github.com/liangyy/encode_ccre_extractor) 
    - The preparation of the summary table is done by `functional_enrichment_encode.py` 
    - `for i in `cat ~/labshare/GitHub/encode_ccre_extractor/query.yaml |cut -f 1 -d':'`; do qsub -v tissue=$i 04_run_encode.sh ; done`
    - `for i in `cat ~/labshare/GitHub/encode_ccre_extractor/query.yaml |cut -f 1 -d':'`; do qsub -v tissue=$i 08_run_encode_all.sh ; done` 
3. GWAS catalog variants. The pre-processing script (map GWAS catalog variant to GTEx variant) was done by the script at [here](https://bitbucket.org/yanyul/rotation-at-imlab/src/master/analysis/annotate_gwas_catalog/).
    - The preparation of the summary table is done by `functional_enrichment_gwas_catalog.py`
    - `for i in `cat ../../mixfine/tissue_w_small_sample_size.txt`; do qsub -v tissue=$i 06_run_gwas_catalog.sh ; done` 
