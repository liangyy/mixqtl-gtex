library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

df_ss = read.csv('https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_sample_counts_by_tissue.csv', stringsAsFactors = F)
df_color = read.csv('https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_tissue_colors.csv', stringsAsFactors = F, skip = '')
df_color$tissue_site_detail_id[df_color$tissue_site_detail_id == 'Cells_Transformed_fibroblasts'] = 'Cells_Cultured_fibroblasts'
color_guide = paste0('#', df_color$tissue_color_hex); names(color_guide) = df_color$tissue_site_detail_id

df = read.csv('signif_summary.csv')

head(df)

df %>% reshape2::dcast(tissue ~ method, value.var = 'unique_pair') %>% ggplot() +
  geom_point(aes(x = eqtl, mixqtl, color = 'mixqtl')) + 
  geom_point(aes(x = eqtl, trcqtl, color = 'trcqtl')) + 
  geom_abline(slope = 1, intercept = 0)


df %>% reshape2::dcast(tissue ~ method, value.var = 'unique_gene') %>% ggplot() +
  geom_point(aes(x = eqtl, mixqtl, color = 'mixqtl')) + 
  geom_point(aes(x = eqtl, trcqtl, color = 'trcqtl')) + 
  geom_abline(slope = 1, intercept = 0)

p1 = df %>% reshape2::dcast(tissue ~ method, value.var = 'unique_pair') %>% ggplot() +
  geom_point(aes(x = eqtl, mixqtl, color = tissue), size = 5, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = color_guide) + 
  theme(legend.position = "none") +
  xlab('Number of significant variant/gene pairs \n (FDR < 0.05) in eQTL') +
  ylab('Number of significant variant/gene pairs \n (FDR < 0.05) in mixQTL') +
  th; p1

p2 = df %>% reshape2::dcast(tissue ~ method, value.var = 'unique_gene') %>% ggplot() +
  geom_point(aes(x = eqtl, mixqtl, color = tissue), size = 5, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = color_guide) + 
  theme(legend.position = "none") +
  xlab('Number of genes with significant variant/gene pairs \n (FDR < 0.05) in eQTL') +
  ylab('Number of genes with significant variant/gene pairs \n (FDR < 0.05) in mixQTL') +
  th; p2  

ggsave('signif_eqtl.png', p1, width = 5, height = 5)
ggsave('signif_gene.png', p2, width = 5, height = 5)

df %>% reshape2::dcast(tissue ~ method, value.var = 'unique_gene') %>% summarize(mean_increase = mean(mixqtl - eqtl)) %>% pander::pander(caption = 'unique_gene')
df %>% reshape2::dcast(tissue ~ method, value.var = 'unique_pair') %>% summarize(mean_increase = mean(mixqtl - eqtl)) %>% pander::pander(caption = 'unique_pair')
