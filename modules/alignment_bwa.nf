include {outputPrefixPath } from './nbt/utils'

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


// process Sort_with_picard {

//   tag { "${sampleId}" }

//   input:
//   tuple sampleId, file(alignedbam)

//   output:
//   tuple fileId, file("${prefix}.aln.sorted.bam")

//   script:
//   prefix=alignedbam.simpleName

//   """
//   java  -jar $EBROOTPICARD/picard.jar SortSam I=${alignedbam} O=${prefix}.aln.sorted.bam SORT_ORDER=coordinate
//   """
// }

// process Mark_duplicates {

//   tag { "${sampleId}" }

//   input:
//   tuple sampleId, file(alignedbamsort)

//   output:
//   tuple fileId, file("${prefix}.aln.sorted.mrkDup.bam") optional true
//   tuple fileId, file("${prefix}.dup_metrics.txt") optional true

//   script:
//   prefix=alignedbamsort.simpleName


//   """
//   java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${alignedbamsort} O=${prefix}.aln.sorted.mrkDup.bam METRICS_FILE=${prefix}.dup_metrics.txt CREATE_INDEX=true

//   rm -f ${prefix}.aln.ba* ${prefix}.aln.sorted.ba*

//   """

// }

// process Base_recalibrator {

//   tag { "${sampleId}" }

//   input:
//   tuple sampleId, file(alignedbamsortmrkDup)

//   output:
//   tuple fileId, file("${prefix}.recal.bam") optional true

//   script:
//   prefix=alignedbamsortmrkDup.simpleName

//   """
//   java -jar $EBROOTGATK/GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ${ref_fasta} -I ${prefix}.aln.sorted.mrkDup.bam -knownSites ${known_sites} -o ${prefix}.recal.table

//   java -jar $EBROOTGATK/GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ${ref_fasta} -I ${prefix}.aln.sorted.mrkDup.bam -BQSR ${prefix}.recal.table -o ${prefix}.recal.bam

//   rm -f ${prefix}.aln.sorted.mrkDup.ba* ${prefix}.recal.table

//   """

// }


