library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 15))

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

df = arrow::read_parquet('../eqtlgen/samples_from_eqtlgen/eqtlgen_pos.mixqtl.parquet')
df = as.data.frame(df)
df = df[, -which(colnames(df) == '__index_level_0__')]
df = df %>% distinct()

clean_up_val = function(df, cols) {
  val = df[, cols[1]]
  val[ is.na(val) ] = df[is.na(val), cols[2]]
  val[ is.na(val) ] = df[is.na(val), cols[3]]
  val
}

get_quantile = function(dd, prob) {
  as.numeric(quantile(dd, probs = prob))
}

bin_by_val = function(val, breaks) {
  labels = c()
  o = rep(NA, length(val))
  o[ val < breaks[1] ] = paste0('< ', breaks[1])
  labels = c(labels, paste0('< ', breaks[1]))
  for(i in 1 : (length(breaks) - 1)) {
    o[ val >= breaks[i] & val < breaks[i + 1]] = paste0(breaks[i], ' - ', breaks[i + 1])
    labels = c(labels, paste0(breaks[i], ' - ', breaks[i + 1]))
  }
  o[val >= breaks[length(breaks)]] = paste0('>=', breaks[length(breaks)])
  labels = c(labels, paste0('>=', breaks[length(breaks)]))
  factor(o, levels = labels)
}

df$bhat = clean_up_val(df, c('beta_meta', 'beta_trc', 'beta_asc'))
df$pval = clean_up_val(df, c('pval_meta', 'pval_trc', 'pval_asc'))


tmp = df %>% filter(pval < 0.01) %>% mutate(pval_bin = bin_by_val(pval, breaks = c(1e-14, 1e-7, 1e-4, 2e-3, 1e-2))) 
tmp2 = tmp %>% group_by(pval_bin) %>% 
  summarize(q05 = get_quantile(bhat, 0.05), q95 = get_quantile(bhat, 0.95))

p = tmp %>% ggplot() + geom_histogram(aes(bhat), binwidth = 0.05) + facet_wrap(~pval_bin, ncol = 2, labeller = label_both) +
  geom_vline(data = tmp2, aes(xintercept = q05), linetype = 2) +
  geom_vline(data = tmp2, aes(xintercept = q95), linetype = 2) +
  geom_label(data = tmp2, 
             aes(label = paste0('90% interval \n', round(q05, 3), ' ~ ', round(q95, 3))), 
             x = 1.7, y = 10000) +
  th2
ggsave('effect_size_dist.png', p, width = 7, height = 8)   
message(nrow(tmp))

tmp3 = tmp %>% filter(pval < 1e-4) %>% 
  summarize(q05 = get_quantile(bhat, 0.05), q95 = get_quantile(bhat, 0.95))
p2 = tmp %>% filter(pval < 1e-4) %>% 
  ggplot() + geom_histogram(aes(bhat), binwidth = 0.05) +
  geom_vline(data = tmp3, aes(xintercept = q05), linetype = 2) +
  geom_vline(data = tmp3, aes(xintercept = q95), linetype = 2) +
  geom_label(data = tmp3, 
             aes(label = paste0('90% interval \n', round(q05, 3), ' ~ ', round(q95, 3))), 
             x = 1.7, y = 10000) +
  th
p2
ggsave('effect_size_dist_p_lt_1e-4.png', p2, width = 6, height = 4)   

df = arrow::read_parquet('signif_mixqtl_in_whole_blood.parquet')
tmp4 = df %>% summarize(q025 = get_quantile(bhat, 0.025), q975 = get_quantile(bhat, 0.975))
p3 = df %>% 
  ggplot() + geom_histogram(aes(bhat), binwidth = 0.05) +
  geom_vline(data = tmp4, aes(xintercept = q025), linetype = 2) +
  geom_vline(data = tmp4, aes(xintercept = q975), linetype = 2) +
  geom_label(data = tmp4, 
             aes(label = paste0('95% interval \n', round(q025, 3), ' ~ ', round(q975, 3))), 
             x = 4, y = 3e5) +
  th
p3
ggsave('effect_size_dist_fdr_0.05.png', p3, width = 6, height = 4)   
df %>% summarize(median = median(abs(bhat)), q95 = quantile(abs(bhat), probs = 0.95))

relative_error_exp_approx = function(x) {
  (exp(x) - 1 - x) / x
}
xseq = seq(-1, 1, by = 0.05)
data.frame(xseq, relative_error_exp_approx(xseq))
