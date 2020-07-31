def parse_varid(ll):
    o = []
    for l in ll:
        g, _, _, _, _ = ll.split('_')
        o.append(g)
    return o

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='split_eqtl_by_chr.py', description='''
        Split eQTL txt.gz file by chromosome and filter according to a gene list.
    ''')
    parser.add_argument('--input', nargs='+', help='''
        GTEx v8 all pairs association in txt.gz format.
    ''')
    parser.add_argument('--egene', nargs='+', help='''
        The list of gene to extract.
        Provide the column name of gene id.
    ''')
    parser.add_argument('--output-prefix', help='''
        Prefix of output.
        Will append .chrN.parquet.
    ''')
    args = parser.parse_args()
    
    import logging, time, sys, os
    from extract_signif_eqtl import read_table
    
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    logging.info('Load gene list.')
    df_g = read_table(args.egene[0])
    gene_list = list(df_g[ args.egene[1] ])
    del df_g
    
    logging.info('Load eQTL table.')
    df = read_table(args.input)
    
    logging.info('Add chr.')
    df['chr'] = parse_varid(list(df['variant_id']))
    
    for chr_num in range(1, 23):
        logging.info(f'Working on chr{chr_num}.')
        df_i = df[ df.chr == f'chr{chr_num}' ].reset_index(drop=True)
        df = df[ df.chr != f'chr{chr_num}' ].reset_index(drop=True)
        df_i = df_i[ df_i.gene_id.isin(gene_list) ].reset_index(drop=True)
        df_i.to_parquet(args.output_prefix + f'.chr{chr_num}.parquet')
        del df_i
    
    logging.info('Done')
    