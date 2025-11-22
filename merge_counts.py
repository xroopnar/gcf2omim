#!/bin/env python
#
import pandas as pd
import numpy as np
import argparse
import sys


def _get_gene_list_from_series(gene_series: pd.Series) -> list:
    if gene_series.empty:
        return []
    all_genes_str = ','.join(gene_series.dropna())
    if not all_genes_str:
        return []
    gene_list = [gene.strip() for gene in all_genes_str.split(',') if gene.strip()]
    return gene_list

def get_total_gene_count(gene_series: pd.Series) -> int:
    return len(_get_gene_list_from_series(gene_series))

def get_unique_gene_string(gene_series: pd.Series) -> str:
    gene_list = _get_gene_list_from_series(gene_series)
    if not gene_list:
        return ""
    
    unique_genes = set(gene_list)
    return ','.join(sorted(list(unique_genes)))

def get_unique_gene_count(gene_series: pd.Series) -> int:
   return len(set(_get_gene_list_from_series(gene_series)))

def process_dataframe(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if 'location' not in df.columns or 'gene' not in df.columns:
        print(f"Warning: DataFrame for prefix '{prefix}' is missing 'location' or 'gene' column.")
        return pd.DataFrame()

    grouped = df.groupby('location')
    agg_df = grouped['gene'].agg(
        total_count=get_total_gene_count,
        unique_genes=get_unique_gene_string,
        unique_count=get_unique_gene_count
    )
    
    # Rename columns to include the specified prefix
    agg_df = agg_df.rename(columns={
        'total_count': f'{prefix}_total_gene_count',
        'unique_genes': f'{prefix}_all_genes',
        'unique_count': f'{prefix}_all_gene_count'
    })
    
    return agg_df

def compare_gene_counts(hg19_df: pd.DataFrame, hg38_df: pd.DataFrame) -> pd.DataFrame:
   
    # Process each dataframe individually
    hg19_results = process_dataframe(hg19_df, 'hg19')
    hg19_df["hg19_location"] = hg19_df.apply(make_loc,axis=1)
    hg19_results = pd.merge(hg19_df[["hg19_location","location"]].drop_duplicates(),hg19_results,on="location",how="right")
    hg38_results = process_dataframe(hg38_df, 'hg38')
    hg38_df["hg38_location"] = hg38_df.apply(make_loc,axis=1)
    hg38_results = pd.merge(hg38_df[["hg38_location","location"]].drop_duplicates(),hg38_results,on="location",how="right")

    output_df = pd.concat([hg19_results, hg38_results], axis=1)
    
    fill_values = {
        'hg19_total_gene_count': 0,
        'hg19_all_genes': '',
        'hg19_all_gene_count': 0,
        'hg38_total_gene_count': 0,
        'hg38_all_genes': '',
        'hg38_all_gene_count': 0
    }
    output_df = output_df.fillna(fill_values)
    
    s19 = output_df['hg19_all_genes'].apply(lambda x: set(x.split(',')) if x else set())
    s38 = output_df['hg38_all_genes'].apply(lambda x: set(x.split(',')) if x else set())

    hg19_only_set = s19 - s38
    output_df['hg19_only_genes'] = hg19_only_set.apply(lambda s: ','.join(sorted(list(s))) if s else '')
    output_df['hg19_only_gene_count'] = hg19_only_set.apply(len)

    hg38_only_set = s38 - s19
    output_df['hg38_only_genes'] = hg38_only_set.apply(lambda s: ','.join(sorted(list(s))) if s else '')
    output_df['hg38_only_gene_count'] = hg38_only_set.apply(len)
    
    output_df['no_genes'] = (output_df['hg19_all_gene_count'] == 0) & \
                            (output_df['hg38_all_gene_count'] == 0)
                            
    output_df['no_diff'] = output_df['hg19_all_genes'] == output_df['hg38_all_genes']
    count_cols = [col for col in output_df.columns if 'count' in col]
    output_df[count_cols] = output_df[count_cols].astype(int)
    output_df = output_df.reset_index().rename(columns={'index': 'location'})
    
    column_order = [
        'hg19_location',
        'hg38_location',
        'location',
        'no_genes',
        'no_diff',
        'hg19_total_gene_count',
        'hg38_total_gene_count',
        'hg19_all_gene_count',
        'hg38_all_gene_count',
        'hg19_all_genes',
        'hg38_all_genes',
        'hg19_only_gene_count',
        'hg19_only_genes',
        'hg38_only_gene_count',
        'hg38_only_genes'
    ]
    final_columns = [col for col in column_order if col in output_df.columns]
    
    return output_df[final_columns].drop("location",axis=1)

def make_loc(x):
    out = x["chrom"]+"-"+str(x["start"])+"-"+str(x["stop"])
    return(out)

def main():

    parser = argparse.ArgumentParser(description='A script that requires that expects combined gene counts files with the columns:\nchrom,start,stop,location,gene_chrom,gene_start,gene_stop,gene,overlap_bp.\n\nThis format is the output of a bedtools intersect with -wao option for overlap_bp to be generated.')
    
    parser.add_argument('--hg19', 
                    type=str, 
                    required=True, 
                    help='Required: Path to the hg19 gene counts.')

    parser.add_argument('--hg38', 
                    type=str, 
                    required=True, 
                    help='Required: Path to the hg38 gene counts.')

    args = parser.parse_args()

    hg19 = args.hg19
    hg38 = args.hg38
    
    hg19_df = pd.read_csv(hg19,sep="\t",header=None).replace(".","")
    hg38_df = pd.read_csv(hg38,sep="\t",header=None).replace(".","")
    cols  = ["chrom","start","stop","location","gene_chrom","gene_start","gene_stop","gene","overlap_bp"]
    hg19_df.columns = cols 
    hg38_df.columns = cols
    output = compare_gene_counts(hg19_df=hg19_df,hg38_df=hg38_df)
    #hg19_df["hg19_location"] = hg19_df.apply(make_loc,axis=1)
    #hg38_df["hg38_location"] = hg38_df.apply(make_loc,axis=1)
    #output = pd.merge(output,hg19_df[["hg19_location","location"]],on="location",how="left")
    #output = pd.merge(output,hg38_df[["hg38_location","location"]],on="location",how="left")
    output.to_csv(sys.stdout,sep="\t",index=False)
    
if __name__=="__main__":
    main()


