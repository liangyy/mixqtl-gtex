# runtime of mixqtl by tissue and sample size
rm(list=ls())
library(ggplot2)
theme_set(theme_bw(base_size = 12))

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
library(dplyr)
library(ggplot2)
df_time = read.table('mixqtl_runtime.tsv', header = T, stringsAsFactors = F)
df_ss = read.csv('https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_sample_counts_by_tissue.csv', stringsAsFactors = F)

df_sample_size = read.table('sample_size_from_pipeline.txt', header = F)
df_ss = inner_join(df_ss, df_sample_size, by = c('tissue' = 'V1'))

df_wall = read.table('mixqtl_walltime.tsv', header = T, stringsAsFactors = F) %>% mutate(total_walltime_in_sec = hour * 3600 + min * 60 + sec)

df_color = read.csv('https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_tissue_colors.csv', stringsAsFactors = F, skip = '')
df_color$tissue_site_detail_id[df_color$tissue_site_detail_id == 'Cells_Transformed_fibroblasts'] = 'Cells_Cultured_fibroblasts'
color_guide = paste0('#', df_color$tissue_color_hex); names(color_guide) = df_color$tissue_site_detail_id
df_c = df_time %>% group_by(tissue) %>% summarize(nchr = n(), total_time_in_min = sum(runtime_in_min), avg_time_in_sec_per_gene = sum(runtime_in_min) / sum(npheno) * 60, npheno = sum(npheno)) %>% left_join(df_ss %>% select(tissue, V2), by = 'tissue')

df_c = inner_join(df_c, df_wall, by = 'tissue')
df_c = df_c %>% mutate(avg_walltime_in_sec_per_gene = total_walltime_in_sec / npheno)

p = df_c %>% filter(nchr == 22) %>% ggplot() + 
  geom_point(aes(x = V2, y = avg_walltime_in_sec_per_gene, color = tissue), size = 5, alpha = 0.5) + 
  scale_color_manual(values = color_guide) + 
  theme(legend.position = "none") + 
  xlab('Sample size') + 
  ylab('Run time per gene (in second)') +
  th; p
ggsave(filename = 'runtime_vs_sample_size.png', p, width = 5, height = 5)
df_c %>% filter(tissue %in% c('Whole_Blood', 'Kidney_Cortex'))
df_time %>% group_by(tissue) %>% summarize(ngene = sum(npheno)) %>% ungroup() %>% summarize(mean(ngene))
sum(df_c$total_walltime_in_sec) / 3600
