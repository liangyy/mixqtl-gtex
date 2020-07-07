# following gtex v8 main paper 
# we use the same number of peer factors as descried 
# in gtex v8 main paper supplementary section 3.5.1
# the rule is
# 15 for n < 150
# 30 for 150 <= n < 250
# 45 for 250 <= n < 350
# 60 for n >= 350

# source of sample size information: https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_sample_counts_by_tissue.csv

wkdir = setwd('/scratch/t.cri.yliang/mixQTL-GTExV8/misc_data/')

df = read.csv('https://bitbucket.org/yanyul/rotation-at-imlab/raw/90f5ada130f89e0a0251eae0034b4687e24b395f/data/gtex_sample_counts_by_tissue.csv')

df = df[ df$v8_all >= 73, ]

assign = function(x) {
  o = rep(15, length(x))
  o[x < 150] = 15
  o[x >= 150 & x < 250] = 30
  o[x >= 250 & x < 350] = 45
  o[x >= 350] = 60
  o
}

out = data.frame(tissue = df$tissue, npeer = assign(df$v8_all))
write.table(out, 'npeer_list.csv', row = F, col = F, sep = ',', quo = F)

