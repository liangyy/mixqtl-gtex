import os
from functional_enrichment import *

def id2pos(ll):
    ch = []
    po = []
    for i in ll:
        c, p, _, _, _ = i.split('_')
        ch.append(c)
        po.append(int(p) - 1)
    return pd.DataFrame({'chr': ch, 'start': po, 'end': np.array(po) + 1, 'var': ll})

def save_bed(ff, filename):
    ff = ff.sort_values(by='start', ascending=True, na_position='first')
    ff.to_csv(filename, sep='\t', index=False, header=False, compression='gzip')
        
def construct_annotation_table_from_variant_list(variant_list, annot_bed, tmp_prefix='test'):
    variant_list = list(set(variant_list))
    df = pd.DataFrame({'variant_id': variant_list})
    df_var = id2pos(variant_list)
    save_bed(df_var, f'{tmp_prefix}.bed.gz')
    sys_call = f'bedtools intersect -a {tmp_prefix}.bed.gz -b {annot_bed} | gzip > {tmp_prefix}.join.bed.gz'
    os.system(sys_call)
    ee = pd.read_csv(f'{tmp_prefix}.join.bed.gz', compression='gzip', sep='\t', header=None)
    annot = df.variant_id.isin(ee.iloc[:, 3]) * 1
    os.remove(f'{tmp_prefix}.join.bed.gz')
    os.remove(f'{tmp_prefix}.bed.gz')
    df = pd.DataFrame({'VARIANT': variant_list, 'cCRE': annot})
    return df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='functional_enrichment_encode.py', description='''
        Form enrichment table using QTL and fine-mapping results and ENCODE cCRE annotation.
    ''')
    parser.add_argument('--tissue', help='''
        Name of tissue
    ''')
    parser.add_argument('--functional_annotation', default='/gpfs/data/im-lab/nas40t2/yanyul/GitHub/encode_ccre_extractor/output/merged_bed_exclude_unclassified.{tissue}.bed.gz', help='''
        Provide the pattern. It is tissue-specific so need include 
        {tissue} as wildcard.
    ''')
    parser.add_argument('--mixqtl', default='/scratch/t.cri.yliang/mixQTL-GTExV8/mixqtl/{tissue}/mixqtl.{tissue}_GTEx_eGene.cis_qtl_pairs.mixQTL.chr{chr_num}.parquet', help='''
        mixQTL results.
    ''')
    parser.add_argument('--top_qtl', default='/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/{tissue}.v8.egenes.txt.gz', help='''
        eQTL results (top SNP only).
    ''')
    parser.add_argument('--strong_gene', default='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine_old/{tissue}/genes-passed-qc.txt.gz', help='''
        Genes passing mixQTL/mixFine QC.
    ''')
    parser.add_argument('--mixfine_pip', default='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine_old/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_pip.mixfine.chr{chr_num}.parquet', help='''
        mixFine pip.
    ''')
    parser.add_argument('--nefine_pip', default='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine_old/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_pip.nefine.chr{chr_num}.parquet', help='''
        neFine pip.
    ''')
    parser.add_argument('--mixfine_cs', default='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine_old/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_cs.mixfine.chr{chr_num}.parquet', help='''
        neFine pip.
    ''')
    parser.add_argument('--nefine_cs', default='/scratch/t.cri.yliang/mixQTL-GTExV8/mixfine_old/{tissue}/mixfine.{tissue}_GTEx_eGene.finemap_cs.nefine.chr{chr_num}.parquet', help='''
        neFine pip.
    ''')
    parser.add_argument('--cache_dir', default='.', help='''
        Directory to save cache.
    ''')
    parser.add_argument('--output', help='''
        output CSV
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    tissue = args.tissue
    
    if os.path.isdir(args.cache_dir):
        cache = args.cache_dir
        cache_top_qtl_mixqtl = f'{cache}/cache_top_qtl_mixqtl.{tissue}.parquet'
        cache_top_qtl_eqtl = f'{cache}/cache_top_qtl_eqtl.{tissue}.parquet'
        cache_top_pip_mixqtl = f'{cache}/cache_top_pip_mixqtl.{tissue}.parquet'
        cache_top_pip_eqtl = f'{cache}/cache_top_pip_eqtl.{tissue}.parquet'
        cache_top_cs_mixqtl = f'{cache}/cache_top_cs_mixqtl.{tissue}.parquet'
        cache_top_cs_eqtl = f'{cache}/cache_top_cs_eqtl.{tissue}.parquet'
        cache_qtl_all = f'{cache}/cache_qtl_all.{tissue}.npy'
        cache_qtl_strong = f'{cache}/cache_qtl_strong.{tissue}.npy'
    else:
        raise ValueError('{} not directory'.format(args.cache_dir))
    
    logging.info('Loading strong gene list')
    df_gene = pd.read_csv(args.strong_gene.format(tissue=tissue), sep='\t', compression='gzip')
    func_annot_bed = args.functional_annotation.format(tissue=tissue)    

    logging.info('Loading QTL/fine-mapping results')
    
    # mixFine CS
    if file_exists(cache_top_cs_mixqtl):
        top_cs_mix = pd.read_parquet(cache_top_cs_mixqtl)
    else:
        top_cs_mix = []
        for i in range(1, 23):
            print(f'Working on chr{i}')
            tmp = pd.read_parquet(args.mixfine_cs.format(chr_num=i, tissue=tissue))
            tmp = tmp[ tmp.variable_prob > 0 ].reset_index(drop=True)
            top_cs_mix.append(tmp[['phenotype_id', 'variant_id', 'cs', 'variable_prob']])
        top_cs_mix = pd.concat(top_cs_mix, axis=0)
        top_cs_mix.to_parquet(cache_top_cs_mixqtl)
    
    # neFine CS
    if file_exists(cache_top_cs_eqtl):
        top_cs_qtl = pd.read_parquet(cache_top_cs_eqtl)
    else:
        top_cs_qtl = []
        for i in range(1, 23):
            print(f'Working on chr{i}')
            tmp = pd.read_parquet(args.nefine_cs.format(chr_num=i, tissue=tissue))
            tmp = tmp[ tmp.variable_prob > 0 ].reset_index(drop=True)
            top_cs_qtl.append(tmp[['phenotype_id', 'variant_id', 'cs', 'variable_prob']])
        top_cs_qtl = pd.concat(top_cs_qtl, axis=0)
        top_cs_qtl.to_parquet(cache_top_cs_eqtl)
    
    # mixFine pip    
    if file_exists(cache_top_pip_mixqtl):
        top_pip_mix = pd.read_parquet(cache_top_pip_mixqtl)
    else:
        top_pip_mix = []
        for i in range(1, 23):
            print(f'Working on chr{i}')
            tmp = pd.read_parquet(args.mixfine_pip.format(chr_num=i, tissue=tissue))
            tmp = tmp[ tmp.variable_prob > 0]
            toppip = tmp.groupby('phenotype_id').apply(get_max_pip)
            top_pip_mix.append(toppip)
        top_pip_mix = pd.concat(top_pip_mix, axis=0)
        top_pip_mix.to_parquet(cache_top_pip_mixqtl)       

    # neFine pip
    if file_exists(cache_top_pip_eqtl):
        top_pip_qtl = pd.read_parquet(cache_top_pip_eqtl)
    else:
        top_pip_qtl = []
        for i in range(1, 23):
            print(f'Working on chr{i}')
            tmp = pd.read_parquet(args.nefine_pip.format(chr_num=i, tissue=tissue))
            tmp = tmp[ tmp.variable_prob > 0]
            toppip = tmp.groupby('phenotype_id').apply(get_max_pip)
            top_pip_qtl.append(toppip)
        top_pip_qtl = pd.concat(top_pip_qtl, axis=0)
        top_pip_qtl.to_parquet(cache_top_pip_eqtl)
    
    # mixQTL
    if file_exists(cache_top_qtl_mixqtl) and file_exists(cache_qtl_all) and file_exists(cache_qtl_strong):
        top_mix = pd.read_parquet(cache_top_qtl_mixqtl)
        df_qtl_all = np.load(cache_qtl_all, allow_pickle=True)
        df_qtl_strong = np.load(cache_qtl_strong, allow_pickle=True)
    else:
        top_mix = []
        df_qtl = []
        for i in range(1, 23):
            print(f'Working on chr{i}')
            tmp = pd.read_parquet(args.mixqtl.format(chr_num=i, tissue=tissue))
            df_qtl.append(tmp[['phenotype_id', 'variant_id']])
            tmp['pval'] = tmp['pval_meta']
            tmp.loc[tmp.pval.isna(), 'pval'] = tmp.pval_trc[tmp['pval'].isna()]
            topqtl = tmp.groupby(['phenotype_id']).apply(get_min)
            top_mix.append(topqtl)
        top_mix = pd.concat(top_mix, axis=0)
        df_qtl = pd.concat(df_qtl, axis=0)
        top_mix.to_parquet(cache_top_qtl_mixqtl)
        df_qtl_strong = df_qtl[ df_qtl.phenotype_id.isin(df_gene.gene) ].variant_id.unique()
        df_qtl_all = df_qtl.variant_id.unique()
        np.save(cache_qtl_all, df_qtl_all)
        np.save(cache_qtl_strong, df_qtl_strong)
    
    # QTL
    if file_exists(cache_top_qtl_eqtl):
        top_eqtl = pd.read_parquet(cache_top_qtl_eqtl)
    else:
        top_eqtl = pd.read_csv(args.top_qtl.format(tissue=tissue), compression='gzip', sep='\t')
        top_eqtl['gene_id'] = trim_dot(top_eqtl['gene_id'].tolist())
        top_eqtl.to_parquet(cache_top_qtl_eqtl)
    gene = in_both(top_eqtl.gene_id, top_mix.index)
    top_mix = top_mix[ top_mix.index.isin(gene) ].reset_index(drop=False)
    top_eqtl = top_eqtl[ top_eqtl.gene_id.isin(gene) ].reset_index(drop=True)
    top_eqtl.rename(columns={'gene_id': 'phenotype_id', 'pval_nominal': 'pmin'}, inplace=True)
    top_eqtl = top_eqtl[['phenotype_id', 'pmin', 'variant_id']]
    
    
    logging.info('Extract baseline')
    # df_baseline = df_annot[ df_annot.VARIANT.isin(df_qtl_all) ].reset_index(drop=True)
    df_baseline = construct_annotation_table_from_variant_list(df_qtl_all, func_annot_bed, tmp_prefix='{}_tmp'.format(tissue))
    df_baseline['total'] = 1
    # df_baseline_strong = df_annot[ df_annot.VARIANT.isin(df_qtl_strong) ].reset_index(drop=True)
    df_baseline_strong = construct_annotation_table_from_variant_list(df_qtl_strong, func_annot_bed, tmp_prefix='{}_strong_tmp'.format(tissue))
    df_baseline_strong['total'] = 1
    
    logging.info('Subset variants')
    # top QTL in strong gene 
    top_mix2 = top_mix[ top_mix.phenotype_id.isin(df_gene.gene) ].reset_index(drop=True)
    top_eqtl2 = top_eqtl[ top_eqtl.phenotype_id.isin(df_gene.gene) ].reset_index(drop=True)
    # shared genes with fine-mapping signals
    gene_pip = in_both(top_pip_mix.index, top_pip_qtl.index)
    top_pip_mix2 = top_pip_mix[ top_pip_mix.index.isin(gene_pip) ].reset_index(drop=False)
    top_pip_qtl2 = top_pip_qtl[ top_pip_qtl.index.isin(gene_pip) ].reset_index(drop=False)
    gene_cs = in_both(top_cs_mix.phenotype_id, top_cs_qtl.phenotype_id)
    top_cs_mix2 = top_cs_mix[ top_cs_mix.phenotype_id.isin(gene_cs) ].reset_index(drop=True)
    top_cs_qtl2 = top_cs_qtl[ top_cs_qtl.phenotype_id.isin(gene_cs) ].reset_index(drop=True)
    # shared CS 
    cs_both = pd.merge(top_cs_mix, top_cs_qtl, on=['phenotype_id', 'variant_id'], suffixes=['_mix', '_qtl'])[['phenotype_id', 'cs_mix', 'cs_qtl']]
    cs_both = cs_both.drop_duplicates()

    cs_mix = glue(cs_both.phenotype_id.values, cs_both.cs_mix.values)
    cs_qtl = glue(cs_both.phenotype_id.values, cs_both.cs_qtl.values)
    top_cs_mix['pheno_cs'] = glue(top_cs_mix.phenotype_id.values, top_cs_mix.cs.values)
    top_cs_qtl['pheno_cs'] = glue(top_cs_qtl.phenotype_id.values, top_cs_qtl.cs.values)

    top_cs_mix3 = top_cs_mix[top_cs_mix.pheno_cs.isin(cs_mix)]
    top_cs_qtl3 = top_cs_qtl[top_cs_qtl.pheno_cs.isin(cs_qtl)]

    top_cs_mix3_topvar = top_cs_mix3.groupby(['phenotype_id', 'cs']).apply(get_max_pip)
    top_cs_qtl3_topvar = top_cs_qtl3.groupby(['phenotype_id', 'cs']).apply(get_max_pip)

    top_cs_mix4 = top_cs_mix[~top_cs_mix.pheno_cs.isin(cs_mix)]
    top_cs_qtl4 = top_cs_qtl[~top_cs_qtl.pheno_cs.isin(cs_qtl)]

    top_cs_mix4_topvar = top_cs_mix4.groupby(['phenotype_id', 'cs']).apply(get_max_pip)
    top_cs_qtl4_topvar = top_cs_qtl4.groupby(['phenotype_id', 'cs']).apply(get_max_pip)
    
    
    logging.info('Construct tables')
    # cols = ['enhancer', 'promoter', 'open_chromatin_region', 'promoter_flanking_region', 'TF_binding_site', 'total']
    cols = ['cCRE', 'total'] 
    results = []
    
    logging.info(' * top QTL in all genes')
    tmp = get_df(df_baseline, top_mix, top_eqtl, cols)
    tmp['type'] = 'top_qtl'
    tmp['gene_list'] = 'all'
    results.append(tmp)
    
    logging.info(' * top QTL in strong genes')
    tmp = get_df(df_baseline_strong, top_mix2, top_eqtl2, cols)
    tmp['type'] = 'top_qtl'
    tmp['gene_list'] = 'strong'
    results.append(tmp)
    
    logging.info(' * top PIP in strong genes')
    tmp = get_df(df_baseline_strong, top_pip_mix, top_pip_qtl, cols)
    tmp['type'] = 'top_pip'
    tmp['gene_list'] = 'strong'
    results.append(tmp)
    
    logging.info(' * top PIP in strong genes (shared genes)')
    tmp = get_df(df_baseline_strong, top_pip_mix2, top_pip_qtl2, cols)
    tmp['type'] = 'top_pip'
    tmp['gene_list'] = 'strong (in both)'
    results.append(tmp)
    
    logging.info(' * 95% CS in strong genes')
    tmp = get_df(df_baseline_strong, top_cs_mix, top_cs_qtl, cols)
    tmp['type'] = '95%_cs'
    tmp['gene_list'] = 'strong'
    results.append(tmp)
    
    logging.info(' * 95% CS in strong genes (shared genes)')
    tmp = get_df(df_baseline_strong, top_cs_mix2, top_cs_qtl2, cols)
    tmp['type'] = '95%_cs'
    tmp['gene_list'] = 'strong (in both)'
    results.append(tmp)
    
    logging.info(' * 95% CS (shared CS) in strong genes (shared genes)')
    tmp = get_df(df_baseline_strong, top_cs_mix3, top_cs_qtl3, cols)
    tmp['type'] = '95%_cs (shared)'
    tmp['gene_list'] = 'strong (in both)'
    results.append(tmp)
    
    logging.info(' * 95% CS (not shared CS) in strong genes (shared genes)')
    tmp = get_df(df_baseline_strong, top_cs_mix4, top_cs_qtl4, cols)
    tmp['type'] = '95%_cs (not shared)'
    tmp['gene_list'] = 'strong'
    results.append(tmp)
    
    logging.info(' * 95% CS top SNP (shared CS) in strong genes (shared genes)')
    tmp = get_df(df_baseline_strong, top_cs_mix3_topvar, top_cs_qtl3_topvar, cols)
    tmp['type'] = '95%_cs (shared) top snp'
    tmp['gene_list'] = 'strong (in both)'
    results.append(tmp)
    
    logging.info(' * 95% CS top SNP (not shared CS) in strong genes (shared genes)')
    tmp = get_df(df_baseline_strong, top_cs_mix4_topvar, top_cs_qtl4_topvar, cols)
    tmp['type'] = '95%_cs (not shared) top snp'
    tmp['gene_list'] = 'strong'
    results.append(tmp)
    
    logging.info('Collect results')
    res = pd.concat(results, axis=0)
    res.to_csv(args.output, index=True)
