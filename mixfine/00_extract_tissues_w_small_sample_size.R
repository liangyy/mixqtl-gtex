df = read.csv('https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_sample_counts_by_tissue.csv')
df$tissue[ df$tissue == 'Cells_Transformed_fibroblasts'] = 'Cells_Cultured_fibroblasts'
df_to_run = df[ df$v8_all < 260 & df$v8_all > 50, ]
write.table(df_to_run$tissue, 'tissue_w_small_sample_size.txt', quo=F, col=F, row=F)

