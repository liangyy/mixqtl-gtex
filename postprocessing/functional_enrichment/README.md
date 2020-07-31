We have two sets of functional annotations. 

1. Regulatory element annotation of all GTEx variants curated by GTEx v8 working group
    - The preparation of the summary table is done by `functional_enrichment.py` 
2. ENCODE cCRE (merged by tissue) pre-processed by script at [here](https://github.com/liangyy/encode_ccre_extractor) 
    - The preparation of the summary table is done by `functional_enrichment_encode.py` 
3. GWAS catalog variants. The pre-processing script (map GWAS catalog variant to GTEx variant) was done by the script at [here](https://bitbucket.org/yanyul/rotation-at-imlab/src/master/analysis/annotate_gwas_catalog/).
    - The preparation of the summary table is done by `functional_enrichment_gwas_catalog.py` 
