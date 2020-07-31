import pandas as pd
import numpy as np
import os.path
import pdb

def read_table(ff):
    fname, file_extension = os.path.splitext(ff)
    if file_extension == '.parquet':
        return pd.read_parquet(ff)
    elif file_extension == '.gz':
        _, nested_ext = os.path.splitext(fname)
        if nested_ext == '.txt':
            return pd.read_csv(ff, compression='gzip', sep='\t', header=0)
        else:
            raise ValueError('Cannot recognize the file name. Not sure how to load the file.')
    
def load_pvalue(df, mode):
    if mode == 'mixqtl':
        pval = df[ 'pval_meta' ].values
        pval[ np.isnan(pval) ] = df[ 'pval_trc' ].values[ np.isnan(pval) ]
        pval[ np.isnan(pval) ] = df[ 'pval_asc' ].values[ np.isnan(pval) ]
        return pval

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='extract_signif_eqtl.py', description='''
        Apply FDR control to eQTL nominal p-values gene by gene.
    ''')
    parser.add_argument('--input', nargs='+', help='''
        Parquet or txt.gz. 
        Also provide the column names of
        variant and phenotype.
        For instance:
        --input test.parquet variant_id phenotype_id
    ''')
    parser.add_argument('--fdr-cutoff', default=0.05, type=float, help='''
        FDR cutoff.
    ''')
    parser.add_argument('--mode', type=str, help='''
        Modes to load p-value.
        Now it supports: eqtl, mixqtl, trcqtl.
    ''')
    parser.add_argument('--tensorqtl', help='''
        Path to tensorqtl repository code folder.
        Now it needs to load rfunc from tensorqtl.
    ''')
    parser.add_argument('--output', help='''
        Output significant gene/variant pairs along with 
        their q-values in parquet format.
    ''')
    args = parser.parse_args()
    
    import logging, time, sys, os
    try:
        sys.path.insert(0, args.tensorqtl)
        import rfunc
    except:
        print('Cannot load tensorqtl.rfunc')
        os.exit()
    
    from tqdm import tqdm
    
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    # load args.input
    if len(args.input) != 3: 
        raise ValueError('Need 3 args for --input.')
    input_file, variant_col, pheno_col = args.input
    
    # load eqtl summary statistics
    logging.info('Loading eQTL table.')
    df = read_table(input_file)
    
    # pool over gene
    logging.info('Looping over genes.')
    res = []
    genes = df[pheno_col].unique()
    for gene in tqdm(genes):
        df_i = df[ df[pheno_col] == gene ].reset_index(drop=True)
        pval = load_pvalue(df_i, mode=args.mode)
        qval = rfunc.qvalue(pval)
        tmp = pd.DataFrame({
            'phenotype_id': df_i[ pheno_col ],
            'variant_id': df_i[ variant_col ],
            'pval': pval,
            'qval': qval
        })
        res.append(tmp[ tmp.qval < args.fdr ].reset_index(drop=True))
    res = pd.concat(res, axis=0)
    
    # save output
    logging.info('Writing output to disk.')
    res.to_parquet(args.output)
        
    
