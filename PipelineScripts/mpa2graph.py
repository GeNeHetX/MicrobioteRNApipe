import glob
import argparse

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

#TODO: Add description
parser = argparse.ArgumentParser()

parser.add_argument("-mpa", type=str,
                    help="path to the combine mpa file")
parser.add_argument("-counts", type=str, nargs="+",
                    help="path to the star read counts summary")
parser.add_argument("-o", type=str,
                    help="output directory for the graphs")

args = parser.parse_args()

#TODO: fix the transpose dataframe (simplfy, choose one)


#REPLACED by get read counts
# def import_count_dict(star_sumary_fp, common_suffix):
# # Import star count dict and fix sample name to correspond to the mpa dataframe
#     count_dict = {}
#     with open(star_sumary_fp, 'r') as in_f:
#         for line in in_f:
#             print(line)
#             if 'Number of input reads' in line:
#                 sample = line.split(':')[0].replace(common_suffix,'')
#                 sample = sample.replace(common_suffix,'')
#                 sample = sample.replace('_Log.final.out','')
#                 read_count = line.split('|')[-1].strip()
#                 count_dict[sample] = read_count
                
                
#     return count_dict

def get_read_counts(star_summary_list, common_suffix):
# Extracts numberof reads in the input fastq from stray summary report
    count_dict = {}
    for f in star_summary_list:
        sample = f.split(':')[0].replace(common_suffix,'')
        sample = sample.replace('_Log.final.out','')

        with open(f, 'r') as in_f:
            for line in in_f:
                if 'Number of input reads' in line:
                    count = line.split('\t')[-1].strip()
                    count_dict[sample] = count
                    break
    return count_dict


def clean_sample_name(df):
# Remove common suffix from all the sample name
    df = df.T
    same_suffix = True
    common_suffix = ''
    i = 1
    idx_list = list(df.index)[1:] # Ignore the first colum (#Classification)
    while same_suffix:
        suffix_list = [x[-i:] for x in idx_list]
        if len(set(suffix_list)) == 1:
            common_suffix = suffix_list[0][-i:]
            i+= 1
        else:
            same_suffix = False
    df['sample_id'] = df.apply(lambda x: x.name.replace(common_suffix, ''), axis=1)
    df.set_index('sample_id', inplace=True)
    # _reads.kreport is not present in ohe rfile names
    return df.T, common_suffix.replace('_reads.kreport', '')
    

def count_to_fraction(data_df, count_dict):
# Calculate the fraction of mapped read for each categories (over total original fastq reads)
    for idx, row in data_df.iterrows():
        sample_read_count = float(count_dict[idx])
        for name, elem in row.items():
            elem = elem/sample_read_count
            data_df.at[idx, name]=elem
    data_df.fillna(0, inplace=True)
    return data_df

def select_rank(mpa_df, rank):
# Extract all the lines of the selected rank from the original mpa combine dataframe
    selection_list = []
    for idx, row in mpa_df.iterrows():
        level = row['#Classification'].split('|')[-1]
        if f'{rank}__' in level and 'Eukaryota' not in row['#Classification']:
            selection_list.append(True)
        else:
            selection_list.append(False)
    # TODO: use deep copy
    df =  mpa_df[selection_list]
    df.loc[:,'#Classification']=df['#Classification'].apply(lambda x: x.split('|')[-1].strip(f'{rank}__'))
    df.set_index('#Classification', inplace=True)
    df=df.astype(np.float64)
    return df.T

def filter_abundance(mpa_df, thres, method):
    # Filter the dataframe in function of the rarity of the reads matching a species/phylum...
    # Method in_sample, min_count, min_frac
    if method == 'in_sample':
        return mpa_df

    elif method == 'min_count':
        return mpa_df

    elif method == 'min_frac':
        return mpa_df
    

mpa_df = pd.read_csv(args.mpa, sep = '\t', header = 0).convert_dtypes()

clean_mpa, common_suffix = clean_sample_name(mpa_df)

#count_dict = import_count_dict(args.counts, common_suffix)
count_dict = get_read_counts(args.counts, common_suffix)


rank_list = ['k','p','c','o','f','g','s']
rank_dict = {'k':'Kingdom','p':'Phylum','c':'Class','o':'Order','f':"Family",'g':'Genius','s':'Species'}

rank_data = {}
for rank in rank_list:
    rank_data[rank]=count_to_fraction(select_rank(clean_mpa, rank), count_dict)


rank_data['p'].to_csv('res.tsv', sep = '\t')

for selected_rank in ['p','c','o','f','g','s']:
    data = count_to_fraction(rank_data[selected_rank], count_dict).T
    fig = sns.clustermap(data, norm=LogNorm())
    plt.savefig(f'{args.o}/{rank_dict[selected_rank]}.png')