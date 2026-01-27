#!/usr/bin/env python


import pandas as pd
import numpy as np
import argparse



def main(input_file):
    
    df = pd.read_csv(input_file, sep='\s+')
   
    bins = np.arange(0, 1.05, 0.05)
    labels = [f"{round(bins[i], 2)}-{round(bins[i+1], 2)}" for i in range(len(bins)-1)]
    counts, _ = np.histogram(df["F_MISS"], bins=bins)
   
    log_counts = np.log10(counts + 1)
    summary_df = pd.DataFrame({"F_MISS_RANGE": labels, "COUNT": counts, "LOG10_COUNT": log_counts})

    output_file = input_file.replace('.lmiss', '_lmiss_count.csv')
    summary_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process lmiss file.")
    parser.add_argument('input_file', help="Input .lmiss file containing lmiss")
    args = parser.parse_args()
    main(args.input_file)

