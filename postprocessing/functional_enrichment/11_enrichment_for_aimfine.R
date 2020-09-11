library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))


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
    for(cc in c('combined', 'eqtl')) {
      mat = .get_table(tot[1, 'baseline'], tot[1, cc], ddi[1, 'baseline'], ddi[1, cc])
      res = fisher.test(mat)
      out[[length(out) + 1]] = data.frame(pval = res$p.value, or = res$estimate, or_95ci_low = res$conf.int[1], or_95ci_high = res$conf.int[2], category = name, method = cc)
    }
  }
  do.call(rbind, out)
}


rename_col = function(dd, map) {
  for(nn in names(map)) {
    dd[dd == nn] = map[[nn]]
  }
  factor(dd, levels = as.character(unlist(map)))
}


source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

# df_aim = read.csv('enrichment/gwas_catalog_for_aimfine_enrichment_Kidney_Cortex_aimfine.csv')
# df_mix = read.csv('enrichment/gwas_catalog_for_aimfine_enrichment_Kidney_Cortex_mixfine.csv')

methods = c('mixfine', 'aimfine')
df_test = list()
for(mm in methods) {
  df = read.csv(paste0('enrichment/gwas_catalog_for_aimfine_enrichment_Kidney_Cortex_', mm, '.csv'))
  colnames(df)[colnames(df) == 'mixqtl'] = 'combined'
  df_type = df %>% select(gene_list, type) %>% distinct()
  for(i in 1 : nrow(df_type)) {
    gl = df_type$gene_list[i]
    ty = df_type$type[i]
    df_test[[length(df_test) + 1]] = get_fisher(
      df %>% filter(gene_list == gl, type == ty)
    ) %>% mutate(gene_list = gl, type = ty, method2 = mm)
  }
}
df_test = do.call(rbind, df_test) %>% mutate(qtl_category = paste(gene_list, type))

rename_map = list(
  # 'strong (in both) 95%_cs (shared)' = '95% CS (in both)',
  # 'strong (in both) 95%_cs (shared)' = '95% CS (in both, shared)',
  # 'strong 95%_cs' = '95% CS',
  # 'strong top_qtl' = 'top QTL',
  'strong (in both) 95%_cs (shared)' = '95% CS (in both, shared)',
  'strong top_pip' = 'top PIP',
  'strong (in both) top_pip' = 'top PIP \n (common gene)'
  # 'strong 95%_cs (not shared) top snp' = 'top PIP within 95% CS \n (distinct CS)',
  # 'strong (in both) 95%_cs (shared) top snp' = 'top PIP within 95% CS \n (common CS)',
  # 'strong 95%_cs (not shared)' = '95% CS (not shared)',
  # 
)
to_show = unique(names(rename_map))

p2 = df_test %>% 
  mutate(qtl_category_rename = rename_col(qtl_category, rename_map)) %>%
  filter(qtl_category %in% to_show) %>% 
  ggplot() + 
  geom_bar(aes(x = qtl_category_rename, y = or, fill = method), position = position_dodge(width=0.8), stat = 'identity') + 
  geom_errorbar(aes(x = qtl_category_rename, ymin = or_95ci_low, ymax = or_95ci_high, group = method), width = 0.5, position = position_dodge(width=0.8), stat = 'identity') +
  facet_wrap(~method2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + th2 +
  theme(axis.title.x = element_blank()) +
  ylab('Odds ratio')
p2
ggsave('enrichment_for_aimfine.png', p2, width = 7, height = 5)
