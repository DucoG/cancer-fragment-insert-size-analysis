import pandas as pd
import numpy as np
import argparse
import os

# Argument parser for list of parquet files
parser = argparse.ArgumentParser(description='Calculate iwFAF from parquet files.')
parser.add_argument('parquet_files', type=str, nargs='+', help='List of parquet files.')
parser.add_argument('--output', type=str, help='Output file path.', default='iwfaf.csv')
args = parser.parse_args()

# load fragment length normalization table from tsv file
df_size_norm = pd.read_csv('/home/d.gaillard/paired_ovarian/fragment_lengh_distibution/data/iwfaf_normalization_tables/size_normalization_This study.tsv', sep='\t')
df_size_norm.drop(columns=['Aberrant fragments', 'Total fragments'], inplace=True)
df_size_norm.rename(columns={'Fragment length (bp)': 'fragment_length', 'FAF': 'FAF_size'}, inplace=True)

df_gc_norm = pd.read_csv('/home/d.gaillard/paired_ovarian/fragment_lengh_distibution/data/iwfaf_normalization_tables/gc_normalization_This study.tsv', sep='\t')
df_gc_norm.drop(columns=['Aberrant fragments', 'Total fragments'], inplace=True)
df_gc_norm.rename(columns={'Fragment GC-content (%)': 'gc_rounded','FAF': 'FAF_gc'}, inplace=True)

# if df_gc_norm['FAF_gc'] is 0, set to nan
df_gc_norm['FAF_gc'] = df_gc_norm['FAF_gc'].replace(0, np.nan)
df_size_norm['FAF_size'] = df_size_norm['FAF_size'].replace(0, np.nan)

results = {}



for file in args.parquet_files:
    print('processing {}'.format(file))
    id = os.path.basename(file).split('_')[0]

    df_full = pd.read_parquet(file, engine='fastparquet')
    df_RPRs = df_full[(df_full['read1_intersect'] == True) | (df_full['read2_intersect'] == True) | (df_full['RPR_overlap'] == True)]
    df_RPRs.loc[:, 'aberrant'] = (df_RPRs['read1_intersect'] | df_RPRs['read2_intersect'])
    df_RPRs.loc[:, 'gc_rounded'] = df_RPRs['gc_content'] * 100
    df_RPRs.loc[:, 'gc_rounded'] = df_RPRs['gc_rounded'].round(0).astype(int)
    df_RPRs = df_RPRs.loc[:, ['fragment_length','gc_rounded', 'read1_intersect', 'read2_intersect', 'RPR_overlap', 'aberrant']]

    df_RPRs = df_RPRs.merge(df_size_norm, on='fragment_length', how='left')
    df_RPRs = df_RPRs.merge(df_gc_norm, on='gc_rounded', how='left')

    df_RPRs['weight_aberrant'] = np.log2(1/(df_RPRs['FAF_size'] * df_RPRs['FAF_gc'])) * df_RPRs['aberrant']
    df_RPRs['weight_non_aberrant'] = np.log2(1/((1-df_RPRs['FAF_size']) * (1-df_RPRs['FAF_gc']))) * ~df_RPRs['aberrant']

    df_RPRs['weight_aberrant'].fillna(0, inplace=True)
    df_RPRs['weight_non_aberrant'].fillna(0, inplace=True)

    iwFAF = df_RPRs['weight_aberrant'].sum()/(df_RPRs['weight_aberrant'].sum() + df_RPRs['weight_non_aberrant'].sum())
    results[id] = iwFAF

print('saving results to {}'.format(args.output))
pd.DataFrame.from_dict(results, orient='index', columns=['iwFAF']).to_csv(args.output)
