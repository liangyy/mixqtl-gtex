grep_runtime = function(tissue) {
  if(!file.exists(paste0("~/scratch/mixQTL-GTExV8/mixqtl/", tissue, "/mixqtl.", tissue, "_GTEx_eGene.chr", 22, ".log"))) {
    return(NULL)
  }
  df = list()
  for(i in 1 : 22) {
    cmd = paste0("cat ~/scratch/mixQTL-GTExV8/mixqtl/", tissue, "/mixqtl.", tissue, "_GTEx_eGene.chr", i, ".log |grep elapsed|awk '{print $3}'")
    tmp = read.table(pipe(cmd), header = FALSE)
    colnames(tmp) = 'runtime_in_min'
    tmp$tissue = tissue
    tmp$chr = i
    cmd = paste0("cat ~/scratch/mixQTL-GTExV8/mixqtl/", tissue, "/mixqtl.", tissue, "_GTEx_eGene.chr", i, ".log |grep pheno|head -n 2 | tail -n 1| awk '{print $2}'")
    npheno = read.table(pipe(cmd))$V1[1]
    tmp$npheno = npheno
    df[[length(df) + 1]] = tmp
  }
  do.call(rbind, df)
}

grep_walltime = function(tissue) {
  if(!file.exists(paste0("logs/run_", tissue, ".out"))) {
    return(NULL)
  }
  cmd = paste0("cat logs/run_", tissue, ".out | tail -n 8|head -n 1|awk -F',' '{print $3}' | sed 's#walltime=##g' ")
  # print(cmd)
  tmp = read.table(pipe(cmd), header = FALSE, sep = ':')
  colnames(tmp) = c('hour', 'min', 'sec')
  tmp$tissue = tissue
  tmp
}

df_tissue = read.csv('../preprocessing/expression/tissue_metainfo.csv', header = FALSE)
df = list()
for(tissue in df_tissue$V1) {
  df[[length(df) + 1]] = grep_runtime(tissue)
}
df = do.call(rbind, df)

write.table(df, 'mixqtl_runtime.tsv', col = T, row = F, quo = F, sep = '\t')

df = list()
for(tissue in df_tissue$V1) {
  df[[length(df) + 1]] = grep_walltime(tissue)
}
df = do.call(rbind, df)

write.table(df, 'mixqtl_walltime.tsv', col = T, row = F, quo = F, sep = '\t')

