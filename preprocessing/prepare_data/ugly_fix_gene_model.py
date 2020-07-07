# adapted from https://github.com/liangyy/mixqtl-pipeline/blob/dev/misc_scripts/prepare_data_for_py_mixqtl_runs/ugly_fix_gene_model.py

import argparse
parser = argparse.ArgumentParser(prog='ugly_fix_gene_model.py', description='''
    Fix the missing suffix in gene ensembl ID.
''')
parser.add_argument('--egene', help='''
    eGene file
''')
parser.add_argument('--gene-model', help='''
    gene annotation
''')
parser.add_argument('--output', help='''
    output
''')
args = parser.parse_args()

import logging, time, sys
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

import pandas as pd
import numpy as np

df1 = pd.read_csv(args.egene, compression='gzip', sep='\t', header=0)
df2 = pd.read_csv(args.gene_model, sep='\t', header=0)
df1['gene_id'] = df1['gene_id'].map(lambda x: x.split('.')[0])
merge = pd.merge(df2, df1, left_on='gene_id', right_on='gene_id', how='left', suffixes=('', '_new'))
merge['start'] = merge[['start', 'gene_start']].apply(lambda x: int(x[0]) if pd.isna(x[1]) else int(x[1]), axis=1)
merge['end'] = merge[['end', 'gene_end']].apply(lambda x: int(x[0]) if pd.isna(x[1]) else int(x[1]), axis=1)
merge['strand'] = merge[['strand', 'strand_new']].apply(lambda x: x[0] if pd.isna(x[1]) else x[1], axis=1)

merge[df2.columns].to_csv(args.output, sep='\t', index=False)


