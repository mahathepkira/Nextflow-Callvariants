include {outputPrefixPath } from './nbt/utils'

process AlignmentSingle {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(read1)

  output:
  tuple val(fileId), file("${fileId}.aln.sorted.bam")

  script:
  prefix = fileId

  """
  bwa mem -t 8 -R "@RG\\tID:${fileId}\\tLB:lib1\\tPL:illumina\\tSM:${prefix}\\tPU:unit1" -M ${params.reference} ${read1} | \

  samtools view -bS -@ 8 - | samtools view -q ${params.mapQ} -b | \

  java -XX:ParallelGCThreads=8 -jar \$EBROOTPICARD/picard.jar SortSam I=/dev/stdin O=${fileId}.aln.sorted.bam SORT_ORDER=coordinate

  """
}

process AlignmentPaired {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"

  input:
  tuple val(fileId), file(read1), file(read2)

  output:
  tuple val(fileId), file("${fileId}.aln.sorted.bam")
  
  script:
  prefix = fileId
  
  """
  bwa mem -t 8 -R "@RG\\tID:${fileId}\\tLB:lib1\\tPL:illumina\\tSM:${prefix}\\tPU:unit1" -M ${params.reference} ${read1} ${read2} | \
  
  samtools view -bS -@ 8 - | samtools view -q ${params.mapQ} -b | \

  java -XX:ParallelGCThreads=8 -jar \$EBROOTPICARD/picard.jar SortSam I=/dev/stdin O=${fileId}.aln.sorted.bam SORT_ORDER=coordinate
   
  """
}

process Mark_duplicates {

  tag { "${fileId}" }

  publishDir "${outputPrefixPath(params, task)}"

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


process GenomicsDBImport {

  input:
  tuple val(region), path(gvcfs), path(gvcftbi)

  output:
  tuple val(region), path("gendb")

  script:
  """
  gatk --java-options "-Xmx${params.gendb_mem ?: '16g'}" GenomicsDBImport \
    -R ${params.reference} \
    -L ${region.region} \
    --genomicsdb-workspace-path gendb \
    \$(printf -- "--variant %s " ${gvcfs})
  """
}

process GenotypeGVCFs {

  input:
  tuple val(region), path(gendb)

  output:
  path("*.vcf.gz")
  path("*.vcf.gz.tbi")

  script:
  def region_tag = region.region.replace(':','_')
  """
  gatk --java-options "-Xmx16g" GenotypeGVCFs \
    --allow-old-rms-mapping-quality-annotation-data \
    -R ${params.reference} \
    -L ${region.region} \
    -V gendb://gendb \
    -O ${region_tag}.vcf.gz
  """
}


process Combine_finalVCF {

 publishDir "${outputPrefixPath(params, task)}"

 input:
 path(vcf)
 path(vcftbi)
 
 output:
 path("allsample.sorted.vcf.gz")
 path("allsample.sorted.vcf.gz.tbi")

 script:
 """
 gatk MergeVcfs ${vcf.collect { "-I $it" }.join(' ')} -O allsample.vcf.gz
 gatk SortVcf -R ${params.reference} -I allsample.vcf.gz -O allsample.sorted.vcf.gz
 gatk IndexFeatureFile -I allsample.sorted.vcf.gz
 """
}

process VCFtoVCFbi {

  tag { "${vcf}" }

  publishDir "${params.output}/Postprocess/VCF_biallelic/"

  input:
  file(vcf)
  file(vcftbi)

  output:
  file("${prefix}_all.snps.indels_bi.vcf.gz")

  script:
  prefix=vcf.simpleName

  """
  bcftools view -m2 -M2 -v snps ${vcf} -o ${prefix}_all.snps.indels_bi.vcf.gz
  """
}

process VCFtoHMP {

  tag { "$vcf}" }

  publishDir "${params.output}/Postprocess/VCF_conversions/"

  input:
  file(vcf)

  output:
  file("${prefix}.hmp.txt")

  script:
  prefix=vcf.simpleName

  """
  run_pipeline.pl -Xmx128G -vcf ${vcf} -sortPositions -export ${prefix}.hmp.txt -exportType Hapmap
  """
}

process VCFtoBEDBIMFAM {

  tag { "${vcf}" }

  publishDir "${params.output}/Postprocess/VCF_conversions/"

  input:
  file(vcf)

  output:
  file("${prefix}.bed")
  file("${prefix}.bim")
  file("${prefix}.fam")

  script:
  prefix=vcf.simpleName

  """
  plink --vcf ${vcf} --make-bed --out ${prefix} --double-id --allow-extra-chr
  """
}

process VcftoBed {

  tag { "${combine_gvcf}" }

  publishDir "${outputPrefixPath(params, task)}"
 //  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(combine_gvcf)
  file(combine_gvcftbi)

  output:
  file("${prefix}_allsample.bed")
  file("${prefix}_allsample.bim")
  file("${prefix}_allsample.fam")
  file("${prefix}_allsample.hmp.txt")
  file("${prefix}_all.snps.indels_bi.vcf.gz")  
  
  script:
  prefix=combine_gvcf.simpleName

  """
  bcftools view -m2 -M2 -v snps ${combine_gvcf} -o ${prefix}_all.snps.indels_bi.vcf.gz
  plink --vcf ${prefix}_all.snps.indels_bi.vcf.gz --make-bed --out ${prefix}_allsample --double-id --allow-extra-chr
  run_pipeline.pl -Xmx128G -vcf ${prefix}_all.snps.indels_bi.vcf.gz -sortPositions -export ${prefix}_allsample.hmp.txt -exportType Hapmap
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


