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
```bash
CHROM   POS     N_ALLELES       N_CHR   {ALLELE:FREQ}
1       2914066 2       552     T:0.713768      G:0.286232
1       2914171 2       556     T:0.751799      A:0.248201
1       3706018 2       554     C:0.519856      T:0.480144
1       4429897 2       554     C:0.819495      T:0.180505
1       4490461 2       554     A:0.830325      T:0.169675
1       8510027 2       560     T:0.585714      A:0.414286
1       9029842 2       554     G:0.790614      T:0.209386
1       9084979 2       554     C:0.819495      G:0.180505
1       9300541 2       560     G:0.514286      A:0.485714
1       10066126        2       552     T:0.753623      C:0.246377
1       10067544        2       552     G:0.695652      A:0.304348
1       10067852        2       556     A:0.938849      T:0.0611511
1       12208400        2       558     G:0.526882      A:0.473118
1       14463522        2       558     A:0.802867      G:0.197133
1       17595070        2       554     T:0.949458      G:0.0505415
1       17601517        2       560     C:0.767857      G:0.232143
1       17680592        2       556     C:0.92446       T:0.0755396
1       21457474        2       554     G:0.509025      T:0.490975
1       22601426        2       554     G:0.68231       A:0.31769
1       23266657        2       552     G:0.778986      A:0.221014
1       26283149        2       556     G:0.863309      C:0.136691
1       28084996        2       560     C:0.885714      T:0.114286
1       28085057        2       558     G:0.81362       T:0.18638
1       34479054        2       558     G:0.921147      T:0.078853
1       36548526        2       556     C:0.723022      A:0.276978
1       36548537        2       552     C:0.905797      T:0.0942029
1       44531771        2       556     G:0.931655      C:0.0683453
1       44531888        2       556     C:0.589928      T:0.410072
1       46255037        2       558     A:0.78853       T:0.21147
1       46810807        2       558     C:0.910394      G:0.0896057
1       49533656        2       552     C:0.865942      A:0.134058
1       51514741        2       552     C:0.913043      A:0.0869565
```
### ตัวอย่างผลลัพธ์จาก {samples}.lmiss
```bash
CHR     POS     N_DATA  N_GENOTYPE_FILTERED     N_MISS  F_MISS
1       2914066 562     0       10      0.0177936
1       2914171 562     0       6       0.0106762
1       3706018 562     0       8       0.0142349
1       4429897 562     0       8       0.0142349
1       4490461 562     0       8       0.0142349
1       8510027 562     0       2       0.00355872
1       9029842 562     0       8       0.0142349
1       9084979 562     0       8       0.0142349
1       9300541 562     0       2       0.00355872
1       10066126        562     0       10      0.0177936
1       10067544        562     0       10      0.0177936
1       10067852        562     0       6       0.0106762
1       12208400        562     0       4       0.00711744
1       14463522        562     0       4       0.00711744
1       17595070        562     0       8       0.0142349
1       17601517        562     0       2       0.00355872
1       17680592        562     0       6       0.0106762
1       21457474        562     0       8       0.0142349
1       22601426        562     0       8       0.0142349
1       23266657        562     0       10      0.0177936
1       26283149        562     0       6       0.0106762
```
### ตัวอย่างผลลัพธ์จาก {samples}.TstV.summary
```bash
MODEL   COUNT
AC      64
AG      188
AT      64
CG      85
CT      226
GT      53
Ts      414
Tv      266
```
### ตัวอย่างผลลัพธ์จาก {samples}.summary
```bash
680 TYPE=snp
```
### ตัวอย่างผลลัพธ์จาก {samples}_allele_frequency.csv
```bash
AF_Range,Count,Log10_Count
"(0.0, 0.05]",0,
"(0.05, 0.1]",408640,5.611340875583322
"(0.1, 0.15]",155040,5.19044375970669
"(0.15, 0.2]",82247,4.915120065804905
"(0.2, 0.25]",50121,4.700019727329084
"(0.25, 0.3]",33086,4.519644265409076
"(0.3, 0.35]",24812,4.394661772492928
"(0.35, 0.4]",22706,4.356140633518128
"(0.4, 0.45]",26107,4.416756969168746
"(0.45, 0.5]",55027,4.740575836289971
"(0.5, 0.55]",0,
"(0.55, 0.6]",0,
"(0.6, 0.65]",0,
"(0.65, 0.7]",0,
"(0.7, 0.75]",0,
"(0.75, 0.8]",0,
"(0.8, 0.85]",0,
"(0.85, 0.9]",0,
"(0.9, 0.95]",0,
"(0.95, 1.0]",0,
```
### ตัวอย่างผลลัพธ์จาก {samples}_lmiss_count.csv
```bash
F_MISS_RANGE,COUNT,LOG10_COUNT
0.0-0.05,25304571,7.403198996076333
0.05-0.1,7310187,6.863928546084085
0.1-0.15,4214411,6.624736989881529
0.15-0.2,2968462,6.472531640091798
0.2-0.25,2240096,6.350266824429099
0.25-0.3,1780134,6.25045293912454
0.3-0.35,1395093,6.144603470911365
0.35-0.4,1233859,6.091265885242811
0.4-0.45,1030921,6.013225807620031
0.45-0.5,873076,5.941052547486134
0.5-0.55,730919,5.863869845549046
0.55-0.6,599437,5.777744270586954
0.6-0.65,480383,5.681588534060951
0.65-0.7,356563,5.552137493100174
0.7-0.75,280559,5.448025752873446
0.75-0.8,199088,5.299047265194292
0.8-0.85,138720,5.1421422108480845
0.85-0.9,86093,4.93497288597527
0.9-0.95,32813,4.516059173758283
0.95-1.0,21080,4.323891208256881
```


