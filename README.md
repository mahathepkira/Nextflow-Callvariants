# nextflow-Callvariants

## หัวข้อ
1. [บทนำ](#1-บทนำ)
2. [การใช้งาน nextflow-Callvariants](#2-การใช้งาน-nextflow-Callvariants)
3. [การเตรียมเครื่องมือและข้อมูลสำหรับ nextflow-Callvariants](#3-การเตรียมเครื่องมือและข้อมูลสำหรับ-nextflow-Callvariants)
4. [รายละเอียดขั้นตอนใน nextflow-Callvariants](#4-รายละเอียดขั้นตอนใน-nextflow-Callvariants)
5. [Output](#6-Output)

---

## 1. บทนำ
Nextflow-Callvariants เป็น bioinformatics pipline ที่พัฒนาขึ้นสำหรับการทำ Variants Calling โดยจะมีขั้นตอนดังต่อไปนี้ 
1. Quality Control
2. Sequence Alignment
3. Quality Mapped
4. Mark Duplicates
5. Base Recalibrate
6. Variants Calling
7. VCF stats
8. Convert VCF to BED,BIM,FAM and hmp

![ภาพ nextflow](Nextflow-Callvariants.drawio.png)


## 2. การใช้งาน nextflow-Callvariants
### การใช้งานแบบไม่ใช้ขั้นตอน Comapare_VCF 
ผู้ใช้งานสามารถใช้คำสั่งต่อไปนี้ในการสั่งใช้งาน nextflow-vep โดยข้อมูลที่อยู่ใน data จะต้องอยู่ในรูป vcf.gz โดย workflow การทำงานจะเป็นไปตามเส้นแดง

```bash
nextflow run main.nf -profile gb --input data --reference <name-species> --outdir results
```
### Options
- `--reference` = reference ที่ต้องการทำ Alignment และ Variants Calling (จำเป็น)
- `--input` = โฟลเดอร์ input (จำเป็น:ค่าเริ่มต้น:data)
- `--outdir` = โฟล์เดอร์ output (จำเป็น:ค่าเริ่มต้น:output)
- `-profile`  = เลือกไฟล์ config ในการรัน Nextflow


## 3. การเตรียมเครื่องมือและข้อมูลสำหรับ nextflow-vep
### เครืองมือ 
Nextflow: version 24
1. Quality Control
   - FastQC
   - Trimmomatric
3. Sequence Alignment
   - BWA
   - samtools
5. Quality Mapped
   - Qualimap
7. Mark Duplicates
   - Picard
9. Variants Calling
   - 
11. VCF stats
12. Convert VCF to BED,BIM,FAM and hmp

### การเตรียม Config
ผู้ใช้งานสามารปรับแต่งเครื่องมือที่ใช้งานในไฟล์ gb.config ให้เหมาะสมกับทรัพยากรในเครื่อง โดย gb.config จะทำงานรวมกับ nextflow.config โดยจะใช้ตัวเลือก `-profile` เพื่อเลือก config ที่จะใช้งาน
```bash
process {
  executor = 'slurm'
  queue = 'memory'
  cache = 'lenient'

  withName: Trimmmomatic {
    module = 'Trimmomatic/0.38-Java-1.8'
    cpus = 8
    memory = '64 GB'
  }

  withName: FastQC {
    module = 'FastQC/0.11.9-Java-11'
    cpus = 8
    memory = '64 GB'
  }

  withName: Alignment_bwa {
    module = 'BWA/0.7.17-GCCcore-11.2.0:SAMtools/1.18-GCC-12.3.0:picard/2.25.1-Java-11:'
    cpus = 8
    memory = '64 GB'
  }

  withName: Mark_duplicates {
  module = 'Java/11.0.16:picard/2.25.1-Java-11'
  cpus = 4
  memory = '32 GB'
  }

  withName: Base_recalibrator {
  module = 'GATK/4.5.0-java-17'
  cpus = 8
  memory = '64 GB'
  }

  withName: Call_GVCF {
  module = 'GATK/4.5.0-java-17:HTSlib/1.19.1-GCC-12.3.0'
  cpus = 10
  memory = '96 GB'
  }

  withName: Combine_GVCF {
  module = 'GATK/3.8-1-Java-1.8.0_144:picard/2.18.27-Java-1.8.0_201:HTSlib/1.9-foss-2018b'
  cpus = 12
  memory = '800 GB'
  }

  withName: VcftoBed {
  module = 'BCFtools/1.17-GCC-12.2.0:PLINK/1.9b_4.1-x86_64:TASSEL/5.2.59'
  cpus = 16
  memory = '32 GB'
  }

  withName: Qualimap {
  container = '/nbt_main/share/singularity/qualimap:2.3--hdfd78af_0'
  cpus = 4
  memory = '32 GB'
  }

  withName: VCFtools_stats {
  module = 'VCFtools/0.1.16-GCC-11.3.0'
  cpus = 8
  memory = '16 GB'
  }

  withName: BCFtools_stats {
  module = 'BCFtools/1.17-GCC-12.2.0'
  cpus = 8
  memory = '16 GB'
  }

  withName: Histogram {
  module = 'Python/3.10.4-GCCcore-11.3.0'
  cpus = 8
  memory = '16 GB'
  }
}

```
 
## 4. รายละเอียดขั้นตอนใน-nextflow-Callvariants
## Trimmmomatic
```bash
process Trimmmomatic {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(read1), file(read2)

  output:
  tuple val(fileId), file("${prefix}_R1_paired.fastq.gz"), file("${prefix}_R2_paired.fastq.gz")
  tuple val(fileId), file("${prefix}_R1_unpaired.fastq.gz"), file("${prefix}_R2_unpaired.fastq.gz")

  script:
  prefix=fileId

  """
  java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 -threads 8 \
  ${read1} ${read2} \
  ${prefix}_R1_paired.fastq.gz ${prefix}_R1_unpaired.fastq.gz \
  ${prefix}_R2_paired.fastq.gz ${prefix}_R2_unpaired.fastq.gz \
  ILLUMINACLIP:/nbt_main/home/lattapol/tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:20 TRAILING:20 \
  SLIDINGWINDOW:4:20 MINLEN:50
  """
}
```
## FastQC
```bash
process FastQC {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(read1), file(read2)

  output:
  tuple val(fileId), file("*.zip"), file("*.html")

  script:
  prefix=fileId

  """
  fastqc ${read1} ${read2} --threads 8

  """
}
```
## Alignment_bwa
```bash
process Alignment_bwa {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(read1), file(read2)

  output:
  tuple val(fileId), file("${fileId}.aln.sorted.bam")

  script:
  def prefix = fileId.tokenize('_')[1]

  """
  bwa mem -t 8 -R "@RG\\tID:${fileId}\\tLB:lib1\\tPL:illumina\\tSM:${prefix}\\tPU:unit1" -M ${params.reference} ${read1} ${read2} | \

  samtools view -bS -@ 8 - | \

  java -XX:ParallelGCThreads=8 -jar \$EBROOTPICARD/picard.jar SortSam I=/dev/stdin O=${fileId}.aln.sorted.bam SORT_ORDER=coordinate

  """
}
```
## Sam_view 
```bash
process Sam_view {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
 //  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(bwa_align)

  output:
  tuple val(fileId), file("${prefix}.aln.bam")

  script:
  prefix=bwa_align.simpleName

  """
  samtools view -bS -@ 12 -o ${prefix}.aln.bam ${bwa_align}
  """
}
```
## Sort_bam
```bash
process Sort_bam {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(samview)

  output:
  tuple val(fileId), file("${prefix}.aln.sorted.bam")

  script:
  prefix=samview.simpleName

  """
  java -XX:ParallelGCThreads=8 -jar \$EBROOTPICARD/picard.jar SortSam I=${samview} O=${prefix}.aln.sorted.bam SORT_ORDER=coordinate
  """
}
```

## Mark_duplicates
```bash
process Mark_duplicates {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
 //  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(sortbam)

  output:
  tuple val(fileId), file("${prefix}.aln.sorted.mrkDup.bam")
  tuple val(fileId), file("${prefix}.dup_metrics.txt")

  script:
  prefix=sortbam.simpleName

  """
  java -XX:ParallelGCThreads=8 -Djava.io.tmpdir=java_temps -jar \$EBROOTPICARD/picard.jar MarkDuplicates I=${sortbam} O=${prefix}.aln.sorted.mrkDup.bam METRICS_FILE=${prefix}.dup_metrics.txt CREATE_INDEX=true
  """
}
```

## Base_recalibrator
```bash
process Base_recalibrator {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
 //  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(markDup)

  output:
  tuple val(fileId), file("${prefix}.recal.bam")

  script:
  prefix=markDup.simpleName

  """
  gatk --java-options "-Xmx16G" BaseRecalibrator -R ${params.reference} -I ${markDup} -known-sites ${params.knownsite} -O ${prefix}.recal.table

  gatk --java-options "-Xmx16G" ApplyBQSR -R ${params.reference} -I ${markDup} --bqsr-recal-file ${prefix}.recal.table -O "${prefix}.recal.bam"
  """
}
```

## Call_GVCF
```bash
process Call_GVCF {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
  // publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(base_recalibrator)

  output:
  tuple val(fileId), file("${fileId}.snps.indels.g.vcf.gz")
  tuple val(fileId), file("${fileId}.snps.indels.g.vcf.gz.tbi")

  script:
  prefix=base_recalibrator.simpleName

  """
  gatk BuildBamIndex -I ${base_recalibrator}
  gatk --java-options "-Xmx32G" HaplotypeCaller -R ${params.reference} -I ${base_recalibrator} -ERC GVCF -O ${prefix}.snps.indels.g.vcf

  bgzip ${prefix}.snps.indels.g.vcf

  tabix -p vcf ${prefix}.snps.indels.g.vcf.gz
  """
}
```

## Combine_GVCF
```bash
process Combine_GVCF {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
 //  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(groupedGVCFs)
  //tuple fileId, file(groupedtbi)

  output:
  file("all.snps.indels.vcf.gz")

  script:

  """
  echo "${fileId.withIndex().collect{ fileId, idx -> "${groupedGVCFs[idx]}" }.join("\n")}" > all_gvcf.list

  getindexvcf.sh

  java -Xmx500G -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ${params.reference} --dbsnp ${params.knownsite} -nt 12 --max_alternate_alleles 6 \$(cat all_gvcf.list | xargs -I {} echo "-V {}") -o all.snps.indels.vcf

  bgzip all.snps.indels.vcf
  """
}
```

## VcftoBed 
```bash
process VcftoBed {

  tag { "${combine_gvcf}" }

  publishDir "${outputPrefixPath(params, task)}"
 //  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(combine_gvcf)

  output:
  file("allsample.bed")
  file("allsample.bim")
  file("allsample.fam")
  file("allsample.hmp.txt")
  file("all.snps.indels_bi.vcf.gz")

  script:

  """
  bcftools view -m2 -M2 -v snps ${combine_gvcf} -o all.snps.indels_bi.vcf.gz
  plink --vcf all.snps.indels_bi.vcf.gz --make-bed --out allsample --double-id --allow-extra-chr
  run_pipeline.pl -Xmx32G -vcf all.snps.indels_bi.vcf.gz -sortPositions -export allsample.hmp.txt -exportType Hapmap
  """
}
```
### VCF stats
สำหรับเครื่องมือชีวสารสนเทศที่ใช้ในขั้นตอนการทำข้อมูลสถิติของ VCF ได้แก่ VCFTools (version 0.1.16) ในการสรุปข้อมูล allele frequency, missing data, Transition/Transversion (Ts/Tv) ratio และ สรุปข้อมูลจำนวน snps
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
ทำการสรุปข้อมูล allele frequency และ missing data สำหรับการทำ histrogram ด้วย Python
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


