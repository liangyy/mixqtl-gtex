library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

df1 = readRDS('top_qtl.rds')
df2 = readRDS('top_qtl_encode.rds')
df3 = readRDS('top_qtl_gwas_catalog.rds')
df = rbind(
  df1 %>% mutate(annot = 'GTEx annotation') %>% filter(!category %in% c('open_chromatin_region', 'promoter_flanking_region')),
  df2 %>% mutate(annot = 'ENCODE'),
  df3 %>% mutate(annot = 'GWAS catalog')
) %>% mutate(new_annot = paste0(annot, '\n', category))

new = rep('All genes', nrow(df))
new[df$qtl_category == 'strong top_qtl'] = 'Strong genes'
df$new_qtl_annot = new


p = df %>% ggplot() + 
  geom_bar(aes(x = new_qtl_annot, y = or, fill = method), position = position_dodge(width=0.8), stat = 'identity') + 
  geom_errorbar(aes(x = new_qtl_annot, ymin = or_95ci_low, ymax = or_95ci_high, group = method), width = 0.2, position = position_dodge(width=0.8), stat = 'identity') +
  facet_wrap(~new_annot, nrow = 1) +
  ylab('Odds ratio') + th2 + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave('top_qtl.png', p, width = 8, height = 4)
