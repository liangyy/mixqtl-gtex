library(optparse)

option_list <- list(
    make_option(c("-g", "--gene_list"), type="character", default=NULL,
                help="gene list to begin with",
                metavar="character"),
    make_option(c("-a", "--asc_cutoff"), type="numeric", default=50,
                help="ASC cutoff (will ignore asc observations below this cutoff)",
                metavar="character"),
    make_option(c("-s", "--asc_nobs_cutoff"), type="numeric", default=15,
                help="Ignore genes with fewer than asc_nobs_cutoff number of ASC observations",
                metavar="character"),
    make_option(c("-c", "--gene_col"), type="character", default=NULL,
                help="the column name of gene in gene_list",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output file",
                metavar="character"),
    make_option(c("-n", "--count_data_rds"), type="character", default=NULL,
                help="Count data RDS (prepared in preprocessing)",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)

trim_dot = function(ss) {
  unlist(lapply(strsplit(ss, '\\.'), function(x) { x[1] }))
}

dl = readRDS(opt$count_data_rds)

# we impose restriction on TRC
trc_cutoff = 100
trc_nobs_cutoff = 50

# restriction on ASC
asc_cutoff = opt$asc_cutoff
asc_nobs_cutoff = opt$asc_nobs_cutoff

trc = rowSums(dl$df_trc >= trc_cutoff)
asc = rowSums(dl$df_ase1 >= asc_cutoff & dl$df_ase2 >= asc_cutoff)
df = data.frame(n_good_trc = trc, n_good_asc = asc, gene = rownames(dl$df_trc))

df = df %>% mutate(pass_trc_qc = n_good_trc >= trc_nobs_cutoff, pass_asc_qc = n_good_asc >= asc_nobs_cutoff)



trc_median_of_indiv_pass_qc = apply(dl$df_trc, 1, function(x) {
  median(x[x >= trc_cutoff])
})
df2 = df %>% filter(pass_asc_qc, pass_trc_qc) %>% select(-pass_asc_qc, -pass_trc_qc)
df2$median_trc = trc_median_of_indiv_pass_qc[match(df2$gene, names(trc_median_of_indiv_pass_qc))]

print(head(df2))
gene_list = read.table(opt$gene_list, header = T, sep = '\t')
df2 = df2[ df2$gene %in% trim_dot(as.character(gene_list[, opt$gene_col])), ]

gz1 <- gzfile(opt$output, "w")
write.table(df2 %>% select(gene, median_trc), gz1, row = F, col = T, quo = F, sep = '\t')
close(gz1)

