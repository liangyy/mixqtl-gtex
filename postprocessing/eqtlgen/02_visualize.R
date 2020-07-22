library(data.table)
library(dplyr)
library(ggplot2)
library(SilverStandardPerformance)
theme_set(theme_bw(base_size = 15))

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

seed_ids = 1 : 10
methods = c('mixqtl', 'eqtl')
sets = c('pos', 'neg')

df = list()

z2p = function(z) {
  2 * exp(pnorm(abs(z), lower.tail = F, log.p = T))
}

p2chisq = function(p) {
  qchisq(p, df=1, lower.tail = F)
}

clean_up_pval = function(df, method) {
  if(method == 'mixqtl') {
    pval = df$pval_meta
    pval[ is.na(pval) ] = df$pval_trc[ is.na(pval) ]
    pval[ is.na(pval) ] = df$pval_asc[ is.na(pval) ]
    pval[ is.na(pval) ] = runif(sum(is.na(pval)))
  }
  if(method == 'eqtl') {
    pval = z2p(abs(df$slope / df$slope_se))
  }
  return(pval)
}

for(mm in methods) {
  df[[mm]] = list()
  for(ss in sets) {
    ff = paste0('samples_from_eqtlgen/eqtlgen_', ss, '.', mm, '.parquet')
    tmp = arrow::read_parquet(ff)
    tmp = tmp[, -which(colnames(tmp) == '__index_level_0__')]
    tmp = tmp %>% distinct()
    tmp$pval = clean_up_pval(tmp, method = mm)
    tmp = tmp[, c('phenotype_id', 'variant_id', 'pval')]
    df[[mm]][[ss]] = tmp
  }
}


df_dat = list()
for(seed in seed_ids) {
  df_dat[[seed]] = list()
  for(ss in sets) {
    df_dat[[seed]][[ss]] = list()
    filename = paste0(paste0('samples_from_eqtlgen/eqtlgen_', ss, '.subsample-with-gene-qc.seed', seed, '.txt.gz'))
    tmp = fread(cmd = paste0('zcat < ', filename), header = F, data.table = F) %>% select(V8, V15) %>% rename(phenotype_id = V8, variant_id = V15)
    for(mm in methods) {
      tmp2 = inner_join(tmp, df[[mm]][[ss]], by = c('phenotype_id', 'variant_id'))
      df_dat[[seed]][[ss]][[mm]] = tmp2
    }
  }
}

# calculate chisq for positive set from pval and visualize
df_chisq = list()
for(seed in seed_ids) {
  for(mm in methods) {
    tmp = df_dat[[seed]][['pos']][[mm]]
    tmp$chisq = p2chisq(tmp$pval)
    df_chisq[[length(df_chisq) + 1]] = tmp %>% mutate(seed = seed, method = mm)
  }
}
df_chisq = do.call(rbind, df_chisq)
p = df_chisq %>% reshape2::dcast(phenotype_id + variant_id + seed ~ method, value.var = 'chisq') %>% 
  ggplot() + geom_bin2d(aes(x = eqtl, y = mixqtl)) + 
  geom_abline(slope = 1, intercept = 0, color = 'lightgray') + 
  xlab(expression('eQTL ' * chi^2)) + 
  ylab(expression('mixQTL ' * chi^2)) +
  facet_wrap(~seed) +
  th2 + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient(limits = c(0, 100), low = 'blue', high = 'red') +
  coord_equal()
ggsave('pairwise_chisq.png', p, width = 10, height = 8)


# roc
breaks = c(0, 1:20/10, 10:15/5, 4:50, Inf)
df_roc = list()
df_pr = list()
for(seed in seed_ids) {
  message('on seed = ', seed)
  tmp_p = inner_join(df_dat[[seed]][['pos']][['mixqtl']], df_dat[[seed]][['pos']][['eqtl']], by = c('phenotype_id', 'variant_id'), suffix = c('.mixqtl', '.eqtl'))
  tmp_n = inner_join(df_dat[[seed]][['neg']][['mixqtl']], df_dat[[seed]][['neg']][['eqtl']], by = c('phenotype_id', 'variant_id'), suffix = c('.mixqtl', '.eqtl'))
  collector = rbind(tmp_p %>% mutate(type = 'positive'), tmp_n %>% mutate(type = 'negative'))
  collector = collector %>% mutate(eqtl = paste(phenotype_id, variant_id))
  true_eqtls = collector %>% filter(type == 'positive')
  
  # ROC
  cur_mix = gen_roc_curve(true_genes = true_eqtls$eqtl, gene = collector$eqtl, score = -log(collector$pval.mixqtl), method = 'gt', cutoff = breaks)
  cur_eqtl = gen_roc_curve(true_genes = true_eqtls$eqtl, gene = collector$eqtl, score = -log(collector$pval.eqtl), method = 'gt', cutoff = breaks)
  e2 = rbind(cur_mix %>% mutate(method = 'mixQTL'), cur_eqtl %>% mutate(method = 'eQTL'))
  df_roc[[length(df_roc) + 1]] = e2 %>% mutate(seed = seed)
  
  # PR
  cur_mix = gen_fdr_power_curve(true_genes = true_eqtls$eqtl, gene = collector$eqtl, score = -log(collector$pval.mixqtl), method = 'gt', cutoff = breaks)
  cur_eqtl = gen_fdr_power_curve(true_genes = true_eqtls$eqtl, gene = collector$eqtl, score = -log(collector$pval.eqtl), method = 'gt', cutoff = breaks)
  e = rbind(cur_mix %>% mutate(method = 'mixQTL') %>% filter(recall != 0) , cur_eqtl %>% mutate(method = 'eQTL') %>% filter(recall != 0))
  df_pr[[length(df_pr) + 1]] = e %>% mutate(seed = seed)
}
df_roc = do.call(rbind, df_roc)
df_pr = do.call(rbind, df_pr)

p = df_roc %>% ggplot() + geom_abline(intercept = 0, slope = 1, color = 'lightgray', group = seed) + 
  geom_path(aes(x = fpr, y = tpr, color = method), size = 1, alpha = .8) + coord_equal() + 
  facet_wrap(~seed) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(legend.position = c(0.8, 0.1), legend.title = element_blank(), legend.text = element_text(size = 12))
p = p + th2
ggsave('roc.png', p, width=10, height=8)

p2 = df_pr %>% ggplot() + 
  geom_path(aes(x = recall, y = precision, color = method), size = 1, alpha = .8) + 
  theme(aspect.ratio = 1) +
  facet_wrap(~seed) 
p2 = p2 + th 
p2 = p2 + xlab('Recall') + ylab('Precision') + 
  theme(legend.position = c(0.8, 0.1), legend.title = element_blank(), legend.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p2
ggsave('pr.png', p2, width=10, height=8)


