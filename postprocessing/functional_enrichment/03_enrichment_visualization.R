library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')



.get_table = function(xtotal, ytotal, xin, yin) {
  matrix(as.numeric(c(xtotal - xin, xin, ytotal - yin, yin)), ncol = 2)
  # matrix(c(yin, ytotal - yin, xin - yin, xtotal - ytotal - xin + yin), ncol = 2)
}

get_fisher = function(dd) {
  tot = dd %>% filter(X == 'total')
  out = list()
  for(i in 1 : (nrow(dd) - 1)) {
    ddi = dd[i, ]
    name = dd$X[i]
    for(cc in c('mixqtl', 'eqtl')) {
      mat = .get_table(tot[1, 'baseline'], tot[1, cc], ddi[1, 'baseline'], ddi[1, cc])
      res = fisher.test(mat)
      out[[length(out) + 1]] = data.frame(pval = res$p.value, or = res$estimate, or_95ci_low = res$conf.int[1], or_95ci_high = res$conf.int[2], category = name, method = cc)
    }
  }
  do.call(rbind, out)
}

include_func = c('TF_binding_site', 'promoter', 'enhancer')  # 'open_chromatin_region', 

tissues = read.table('../../mixfine/tissue_w_small_sample_size.txt')$V1
# tissues = tissues[ -which(tissues == 'Cells_EBV-transformed_lymphocytes')]
# tissues = c('Liver', 'Uterus', 'Ovary', 'Kidney_Cortex')

df_test = list()
df_combine = list()
for(tt in tissues) {
  df = read.csv(paste0('enrichment/enrichment_', tt, '.csv'))
  df_combine[[length(df_combine) + 1]] = df %>% mutate(tissue = tt)
  df_type = df %>% select(gene_list, type) %>% distinct()
  for(i in 1 : nrow(df_type)) {
    gl = df_type$gene_list[i]
    ty = df_type$type[i]
    df_test[[length(df_test) + 1]] = get_fisher(
      df %>% filter(gene_list == gl, type == ty)
    ) %>% mutate(gene_list = gl, type = ty, tissue = tt)
  }
}
df_test = do.call(rbind, df_test)
df_combine = do.call(rbind, df_combine)
df_combine = df_combine %>% group_by(gene_list, type, X) %>% 
  summarize(baseline = sum(baseline), mixqtl = sum(mixqtl), eqtl = sum(eqtl)) %>%
  ungroup()
df_test = df_test %>% mutate(qtl_category = paste(gene_list, type))

df_type = df_combine %>% select(gene_list, type) %>% distinct()
df_test_c = list()
for(i in 1 : nrow(df_type)) {
  gl = df_type$gene_list[i]
  ty = df_type$type[i]
  df_test_c[[length(df_test_c) + 1]] = get_fisher(
    df_combine %>% filter(gene_list == gl, type == ty)
  ) %>% mutate(gene_list = gl, type = ty)
}
df_test_c = do.call(rbind, df_test_c)
df_test_c = df_test_c %>% mutate(qtl_category = paste(gene_list, type))

to_show = c('strong top_qtl', 'strong top_pip', 'strong (in both) top_pip', 'strong (in both) 95%_cs (shared) top snp') # , 'strong 95%_cs (not shared) top snp')
rename_map = list(
  'strong top_qtl' = 'top QTL',
  'strong top_pip' = 'top PIP',
  'strong (in both) top_pip' = 'top PIP \n (common gene)',
  'strong (in both) 95%_cs (shared) top snp' = 'top PIP within 95% CS \n (common CS)'
  # 'strong 95%_cs (not shared) top snp' = 'top PIP within 95% CS \n (distinct CS)'
)

df_test %>% filter(category %in% include_func, qtl_category %in% to_show) %>% ggplot() + 
  geom_bar(aes(x = qtl_category, y = or, fill = method), position = position_dodge(width=0.8), stat = 'identity') + 
  geom_errorbar(aes(x = qtl_category, ymin = or_95ci_low, ymax = or_95ci_high, group = method), width = 0.5, position = position_dodge(width=0.8), stat = 'identity') +
  facet_grid(category~tissue, scales = 'free_y') + 
  theme(axis.text.x = element_text(angle = 90))

rename_col = function(dd, map) {
  for(nn in names(map)) {
    dd[dd == nn] = map[[nn]]
  }
  factor(dd, levels = as.character(unlist(map)))
}


df_test_c = df_test_c %>% mutate(qtl_category_rename = rename_col(qtl_category, rename_map))
p = df_test_c %>% 
  filter(category %in% include_func, qtl_category %in% to_show) %>%
  ggplot() + 
  geom_bar(aes(x = qtl_category_rename, y = or, fill = method), position = position_dodge(width=0.8), stat = 'identity') + 
  geom_errorbar(aes(x = qtl_category_rename, ymin = or_95ci_low, ymax = or_95ci_high, group = method), width = 0.2, position = position_dodge(width=0.8), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~category, scales = 'free_y', ncol = 3) + th2 +
  theme(axis.title.x = element_blank()) +
  ylab('Odds ratio')
df_combine %>% mutate(qtl_category = paste(gene_list, type)) %>% filter(qtl_category %in% to_show, X %in% include_func) %>% print(n=100)
ggsave('functional_enrichment.png', p, width = 7, height = 5)

get_se = function(low_ci, high_ci, signif_level = 0.95) {
  (high_ci - low_ci) / 2 / abs(qnorm((1 - signif_level) / 2))
}
get_diff = function(val, label) {
  val[label == 'mixqtl'] - val[label == 'eqtl']
}
get_diff_se = function(val, label) {
  sqrt(val[label == 'mixqtl'] ^ 2 + val[label == 'eqtl'] ^ 2)
}
inverse_variance_meta = function(bhat, se) {
  w = 1 / se ^ 2
  sum(w * bhat) / sqrt(sum(w))
}

df_test = df_test %>% mutate(or_se = get_se(or_95ci_low, or_95ci_high))
df_test_summary = df_test %>% group_by(category, qtl_category, tissue) %>% summarize(or_diff = get_diff(or, method), or_diff_se = get_diff_se(or_se, method)) %>% ungroup()
df_test_summary = df_test_summary %>% group_by(category, qtl_category) %>% summarize(z = inverse_variance_meta(or_diff, or_diff_se))
df_test_summary %>% 
  filter(category %in% include_func, qtl_category %in% to_show)
