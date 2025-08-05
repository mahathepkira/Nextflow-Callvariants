include {outputPrefixPath } from './nbt/utils'

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



process Qualimap {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(bam)

  output:
  path "*"

  script:
  prefix=fileId

  """
  qualimap bamqc -bam ${bam} --java-mem-size=32G

  """
}


process BCFtools_stats {

  tag { "${vcf}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(vcf)

  output:
  path "*"

  script:
  prefix=vcf

  """
  tabix -p vcf ${vcf}
  bcftools stats --threads 8 ${vcf} > ${prefix}.stat
  """
}


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
