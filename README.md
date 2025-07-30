# nextflow-Callvariants

## หัวข้อ
1. [บทนำ](#1-บทนำ)
2. [การใช้งาน nextflow-Callvariants](#2-การใช้งาน-nextflow-Callvariants)
3. [การเตรียมเครื่องมือและข้อมูลสำหรับ nextflow-Callvariants](#3-การเตรียมเครื่องมือและข้อมูลสำหรับ-nextflow-Callvariants)
4. [รายละเอียดขั้นตอนใน nextflow-Callvariants](#4-รายละเอียดขั้นตอนใน-nextflow-Callvariants)
5. [การปรับแต่งการ Annotations ใน VEP](#5-การปรับแต่งการ-Annotations-ใน-VEP)
6. [Output](#6-Output)

---

## 1. บทนำ
nextflow-vep เป็น bioinformatics pipline ที่พัฒนาขึ้นสำหรับการทำ Variants Calling โดยจะมีขั้นตอนดังต่อไปนี้ 
1. Quality Control
2. Sequence Alignment
3. Quality Mapped
4. Mark Duplicates
5. Base Recalibrate
6. Variants Calling
7. VCF stats
8. Convert VCF to BED,BIM,FAM and hmp

![ภาพ nextflow](Nextflow-Callvariants.drawio.png)

## 4. รายละเอียดขั้นตอนใน-nextflow-Callvariants
### VCF stats
```bash
process VCFtools_stats {

  tag { "${vcf}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(vcf)

  output:
  file("*.frq")
  file("*.lmiss")
  file("*.TsTv.summary")
  file("*.summary")
  script:
  prefix=vcf

  """
  bash /nbt_main/home/lattapol/nextflow-Callvariants/bin/quality.sh ${vcf}
  """
}
```
```bash

process Histogram {

  tag { "${frq}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(frq)
  file(lmiss)

  output:
  file("*.csv")

  script:

  """
  python /nbt_main/home/lattapol/nextflow-Callvariants/bin/create_AF_his.py ${frq}
  python /nbt_main/home/lattapol/nextflow-Callvariants/bin/create_lmiss_his.py ${lmiss}
  """
}
```
รายละเอียดใน quality.sh 
```bash
#!/bin/bash


echo "===== running ====="

if [ -z "$1" ]; then
    echo "Usage: bash quality.sh <input_vcf.gz>"
    exit 1
fi


vcf="$1"
prefix=$(basename "$vcf" .vcf.gz)

vcftools --gzvcf "$vcf" --freq --out "${prefix}"
vcftools --gzvcf "$vcf" --missing-site --out "${prefix}"
vcftools --gzvcf "$vcf" --TsTv-summary --out "${prefix}"

zcat "$vcf" | vcf-annotate --fill-type | grep -oP "TYPE=\\w+" | sort | uniq -c > "${prefix}.summary"

echo "===== finished ====="
```
รายละเอียดใน create_AF_his.py
```bash
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

```
รายละเอียดใน create_lmiss_his.py 
```bash
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
```
## 6. Output
### ตัวอย่างผลลัพธ์จาก {samples}.frq
### ตัวอย่างผลลัพธ์จาก {samples}.lmiss
### ตัวอย่างผลลัพธ์จาก {samples}.TstV.summary
### ตัวอย่างผลลัพธ์จาก {samples}.summary
### ตัวอย่างผลลัพธ์จาก {samples}_allele_frequency.csv
### ตัวอย่างผลลัพธ์จาก {samples}_lmiss_count.csv


