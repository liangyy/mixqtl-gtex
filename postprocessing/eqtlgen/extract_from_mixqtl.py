import pandas as pd
import sys
sys.path.insert(0, '../functional_enrichment/')
import functional_enrichment

def trim_dot(ss):
    return [ i.split('.')[0] for i in ss ]

pos = 'samples_from_eqtlgen/eqtlgen_pos.subsample-with-gene-qc.seed{idx}.txt.gz'
neg = 'samples_from_eqtlgen/eqtlgen_neg.subsample-with-gene-qc.seed{idx}.txt.gz'
mix = '/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl/Whole_Blood/mixqtl.Whole_Blood_GTEx_eGene.cis_qtl_pairs.mixQTL.chr{chr_num}.parquet'
qtl = '/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/Whole_Blood.allpairs.txt.gz'
cache_mix_pos = 'samples_from_eqtlgen/eqtlgen_pos.mixqtl.parquet'
cache_mix_neg = 'samples_from_eqtlgen/eqtlgen_neg.mixqtl.parquet'
cache_qtl_pos = 'samples_from_eqtlgen/eqtlgen_pos.eqtl.parquet'
cache_qtl_neg = 'samples_from_eqtlgen/eqtlgen_neg.eqtl.parquet'

if functional_enrichment.file_exists(cache_mix_pos) and functional_enrichment.file_exists(cache_mix_neg):
    pass
else:
    # aggregate all positive/negative variant/gene pairs from samples
    df_pos = []
    df_neg = []
    print('Load positive/negative pairs', flush=True)
    for i in range(1, 11):
        pp = pd.read_csv(pos.format(idx=i), compression='gzip', sep='\t', header=None)
        pp = pp.iloc[:, [7, 14] ]
        df_pos.append(pp)
        nn = pd.read_csv(neg.format(idx=i), compression='gzip', sep='\t', header=None)
        nn = nn.iloc[:, [7, 14] ]
        df_neg.append(nn)
    df_pos = pd.concat(df_pos, axis=0)
    df_neg = pd.concat(df_neg, axis=0)
    df_pos.columns = ['phenotype_id', 'variant_id']
    df_neg.columns = ['phenotype_id', 'variant_id']

    # load mixQTLs that appear in positive and negative sets
    mix_pos = []
    mix_neg = []
    for i in range(1, 23):
        print(f'mixQTL: working on chr{i}', flush=True)
        dd = pd.read_parquet(mix.format(chr_num=i))
        mix_pos.append(pd.merge(dd, df_pos, on=['phenotype_id', 'variant_id']))
        mix_neg.append(pd.merge(dd, df_neg, on=['phenotype_id', 'variant_id']))

    mix_pos = pd.concat(mix_pos, axis=0)
    mix_neg = pd.concat(mix_neg, axis=0)
    mix_pos.to_parquet(cache_mix_pos)
    mix_neg.to_parquet(cache_mix_neg)

if functional_enrichment.file_exists(cache_qtl_pos) and functional_enrichment.file_exists(cache_qtl_neg):
    pass
else:
    # extract from eQTL
    print('Load eQTL all pairs', flush=True)
    df_qtl = pd.read_csv(qtl, compression='gzip', sep='\t', header=0)
    df_qtl_pos = df_qtl[ df_qtl.variant_id.isin(df_pos.variant_id) ]
    df_qtl_neg = df_qtl[ df_qtl.variant_id.isin(df_neg.variant_id) ]
    df_qtl_pos['phenotype_id'] = trim_dot(df_qtl_pos['gene_id'])
    df_qtl_neg['phenotype_id'] = trim_dot(df_qtl_neg['gene_id'])

    qtl_pos = pd.merge(df_qtl_pos, df_pos, on=['phenotype_id', 'variant_id'])
    qtl_neg = pd.merge(df_qtl_neg, df_neg, on=['phenotype_id', 'variant_id'])
    qtl_pos.to_parquet(cache_qtl_pos)
    qtl_neg.to_parquet(cache_qtl_neg)

