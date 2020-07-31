# About

Here we implement the simple FDR routine to obtain significant eQTL from nominal p-values.

The script `extract_signif_eqtl.py` will apply `qvalue` 
(R function with Python wrapper provided in `tensorqtl` with default parameters) 
to the sequence of p-values for each gene. 
The pairs exceeding an FDR cutoff will be output.

Run example:

```
bash 00_submit_split_eqtl.sh /scratch/t.cri.yliang/mixQTL-GTExV8/split_eqtl
``` 

