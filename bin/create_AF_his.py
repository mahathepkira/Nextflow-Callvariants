#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
#import matplotlib.pyplot as plt




def main(input_file):
    
    data = pd.read_csv(input_file, sep='\t', usecols=['{ALLELE:FREQ}'])

   
    alleles_freqs = data['{ALLELE:FREQ}'].str.extract(r'(?P<Allele>[ACGT<N>])\:(?P<Frequency>[0-1]?[0-9]*(?:\.[0-9]+)?)')
    alleles_freqs['Frequency'] = pd.to_numeric(alleles_freqs['Frequency'], errors='coerce')

   
    bins = np.arange(0, 1.05, 0.05)
    alleles_freqs['AF_binned'] = pd.cut(alleles_freqs['Frequency'], bins)
    counts = alleles_freqs['AF_binned'].value_counts().sort_index()

    
    table = pd.DataFrame({'AF_Range': counts.index.astype(str), 'Count': counts.values})
    table['Log10_Count'] = np.log10(table['Count'].replace(0, np.nan))

    
    output_file = input_file.replace('.frq', '_allele_frequency.csv')
    table.to_csv(output_file, index=False)
    """
    plt.figure(figsize=(10, 6))
    plt.bar(table['AF_Range'], table['Log10_Count'], color='blue', alpha=0.7)
    plt.xlabel('Allele Frequency Range')
    plt.ylabel('Log10_Count')
    plt.title('Histogram of Allele Frequency')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('allele_frequency_histogram.png', dpi=300)
    """
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process allele frequency file.")
    parser.add_argument('input_file', help="Input .frq file containing allele frequencies")
    args = parser.parse_args()
    main(args.input_file)

