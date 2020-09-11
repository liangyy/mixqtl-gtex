
'''
This is an ad hoc script to aggregate the p-values from TRC, ASC, META in mixQTL.
The aggregation function has been added to the mixqtl function so it won't need 
in the future.
In other word, this script is one-time business.
'''


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='clean_up_parquet.py', description='''
        Add the cleaned-up p-value, bhat, bhat_se to the old mixQTL output.
    ''')
    parser.add_argument('--input-parquet', help='''
        Input parquet file.
    ''')
    parser.add_argument('--output-parquet', help='''
        Output parquet file.
    ''')
    parser.add_argument('--path-to-tensorqtl', help='''
        Path to tensorqtl to load the function for p-value cleaning.
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
    
    import pandas as pd
    
    sys.path.insert(0, '{}/tensorqtl'.format(args.path_to_tensorqtl))
    import mixqtl
    
    logging.info('Loading input')
    df = pd.read_parquet(args.input_parquet)
    
    logging.info('Converting')
    df = mixqtl.cleanup_pvalue(df)
    
    logging.info('Writing output')
    df = pd.to_parquet(args.output_parquet)
    