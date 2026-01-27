nextflow.enable.dsl=2
/*
================================================================================
=                           Sinonkt Style I N I T                              =
================================================================================
*/

include {readAvroSchema} from './modules/nbt/utils'
include {getDefaultThenResolveParams} from './modules/nbt/utils'
include {groupTupleWithoutCommonKey} from './modules/nbt/utils'

if (params.exportKeySchema) exit 0, printKeySchema()
if (params.exportValueSchema) exit 0, printValueSchema()

params.MAINTAINERS = [
  'Krittin Phornsiricharoenphant (oatkrittin@gmail.com)'
]

def schema = readAvroSchema("${workflow.projectDir}/schemas/value.avsc")
__params = getDefaultThenResolveParams(schema, params)

include {handleErrorMessage} from './modules/nbt/log' params(__params)
include {handleCompleteMessage} from './modules/nbt/log' params(__params)
include {helpMessage} from './modules/nbt/help' params(__params)

if (params.version) exit 0, workflowVersionMessage()
if (params.help) exit 0, helpMessage(schema)

/*
================================================================================
=                   Sinonkt Style Workflows definitions                        =
================================================================================
*/
include { AlignmentSingle }  from './modules/alignment_bwa' params(__params)
include { AlignmentPaired }  from './modules/alignment_bwa' params(__params)
include { Mark_duplicates }  from './modules/alignment_bwa' params(__params)
include { Base_recalibrator }  from './modules/alignment_bwa' params(__params)
include { Call_GVCF }  from './modules/alignment_bwa' params(__params)
include { GenomicsDBImport }  from './modules/alignment_bwa' params(__params)
include { GenotypeGVCFs }  from './modules/alignment_bwa' params(__params)
include { Combine_finalVCF }  from './modules/alignment_bwa' params(__params)
include { VcftoBed }  from './modules/alignment_bwa' params(__params)
include { VCFtoVCFbi }  from './modules/alignment_bwa' params(__params)
include { VCFtoHMP }  from './modules/alignment_bwa' params(__params)
include { VCFtoBEDBIMFAM }  from './modules/alignment_bwa' params(__params)
include { Trimmmomatic_Single } from './modules/quality' params(__params)
include { Trimmmomatic_Paired } from './modules/quality' params(__params)
include { FastqcForPaired as FastqcForPairedBefore } from './modules/quality' params(__params)
include { FastqcForPaired as FastqcForPairedAfter } from './modules/quality' params(__params)
include { FastqcForSingle as FastqcForSingleBefore } from './modules/quality' params(__params)
include { FastqcForSingle as FastqcForSingleAfter } from './modules/quality' params(__params)
include { Qualimap } from './modules/quality' params(__params)
include { VCFtoolsBefore_stats } from './modules/quality' params(__params)
include { VCFtoolsAfter_stats } from './modules/quality' params(__params)
include { BCFtoolsBefore_stats } from './modules/quality' params(__params)
include { BCFtoolsAfter_stats } from './modules/quality' params(__params)
include { HistogramBefore } from './modules/quality' params(__params)
include { HistogramAfter } from './modules/quality' params(__params)
include { FastqcForSingle_visualize } from './modules/quality' params(__params)
include { FastqcForPaired_visualize } from './modules/quality' params(__params)
include { Qualimap_visualize } from './modules/quality' params(__params)
include { MultiQC as multiqc_fastqc_af } from './modules/quality' params(__params)
include { MultiQC as multiqc_fastqc_bf } from './modules/quality' params(__params)
include { MultiQC as multiqc_qualimap } from './modules/quality' params(__params)

def makeChrChannel(fasta) {
    Channel
        .fromPath("${fasta}.fai")
        .splitText()
        .map { line ->
            def f = line.split('\t')
            [
                chr : f[0],
                len : f[1].toLong()
            ]
        }
}


def makeWindowChannel(chr_ch, window_size) {
    chr_ch.flatMap { meta ->
        (1..meta.len).step(window_size).collect { start ->
            def end = Math.min(start + window_size - 1, meta.len)
            [
                chr    : meta.chr,
                start  : start,
                end    : end,
                region : "${meta.chr}:${start}-${end}"
            ]
        }
    }
}



workflow Alignment {
   take:
     reads_trimed

   main:
     if (__params.reads_type == "SE"){
         sortbam = AlignmentSingle(reads_trimed)
     }
     else if (__params.reads_type == "PE"){
         sortbam = AlignmentPaired(reads_trimed)
     }

     (markDup, dup_metrics)=Mark_duplicates(sortbam)

     if (__params.vcfknowsite == "on"){
        final_bam = Base_recalibrator(markDup)
     }
     else if (__params.vcfknowsite == "off"){
        final_bam = markDup
     }
     
     qmap = Qualimap(final_bam).toList()
     Qualimap_visualize(qmap)

   emit:
     final_bam
     qmap
}

workflow Callvariant {
   take:
     bam
     regions_ch

   main:
     (gvcf, gvcftbi) = Call_GVCF(bam)
     gvcf_files = gvcf.map { sid, gvcf_path -> gvcf_path }
     gvcftbi_files = gvcftbi.map { sid, gvcftbi_path -> gvcftbi_path }
     gvcf_list = gvcf_files.toList()
     gvcftbi_list = gvcftbi_files.toList()
     gvcf_by_region = regions_ch.combine(gvcf_list.toList()).combine(gvcftbi_list.toList()).map { region, gvcf, gvcftbi -> [region, gvcf, gvcftbi] }
     gendb = GenomicsDBImport(gvcf_by_region)
     (vcf,vcftbi) = GenotypeGVCFs(gendb)
     vcf_all = vcf.collect()
     vcftbi_all = vcftbi.collect()
     (final_vcf,final_vcftbi) = Combine_finalVCF(vcf_all,vcftbi_all)

   emit:
     final_vcf
     final_vcftbi
}

workflow MultiQC {
   take:
     qc_af
     qc_bf
     qualimap
     

   main:
     qc_af_results = multiqc_fastqc_af(qc_af)
     qc_bf_results = multiqc_fastqc_bf(qc_bf)
     qualimap_results = multiqc_qualimap(qualimap)

   emit:
     qc_af_results
     qc_bf_results
     qualimap_results
}


workflow Postprocess {
   take:
     vcf
     vcftbi

   main:
      
     if (__params.bi_allelic == "filter") {
        (frq1,lmiss1,tstv1,sum1) = VCFtoolsBefore_stats(vcf)
        BCFtoolsBefore_stats(vcf)
        HistogramBefore(frq1,lmiss1)
        vcf_final = VCFtoVCFbi(vcf,vcftbi)
        (frq2,lmiss2,tstv2,sum2) = VCFtoolsAfter_stats(vcf_final)
        BCFtoolsAfter_stats(vcf_final)
        HistogramAfter(frq2,lmiss2)
     }
     else if (__params.bi_allelic == "non-filter") {
        vcf_final = vcf
     }
 
     if (__params.export == "hmp") {
        VCFtoHMP(vcf_final)
     }
     else if (__params.export == "bedbimfam") {
        VCFtoBEDBIMFAM(vcf_final)
     }
     else if (__params.export == "both") {
        VCFtoHMP(vcf_final)
        VCFtoBEDBIMFAM(vcf_final)
     }

   emit:
     vcf_final
}


workflow QC {
   take:
     reads

   main:

     if (__params.reads_type == "SE") {     
        report_qc1 = FastqcForSingleBefore(reads)
        reads_trim = Trimmmomatic_Single(reads)
        report_qc2 = FastqcForSingleAfter(reads_trim)
        qc1 = report_qc1.map {[it[1]]}.flatten().toList()
        qc2 = report_qc2.map {[it[1]]}.flatten().toList()
        combined_qc = qc1+qc2
        sum_qc = FastqcForSingle_visualize(combined_qc)
    }
    else if (__params.reads_type == "PE") {
        report_qc1 = FastqcForPairedBefore(reads)
        reads_trim = Trimmmomatic_Paired(reads)
        report_qc2 = FastqcForPairedAfter(reads_trim)
        qc1 = report_qc1.map {[it[1]]}.flatten().toList()
        qc2 = report_qc2.map {[it[1]]}.flatten().toList()
        combined_qc = qc1+qc2
        sum_qc = FastqcForPaired_visualize(combined_qc)
    }

   emit:
     reads_trim
     qc1
     qc2
}





/*
================================================================================
=                           Sinonkt Style M A I N                              =
================================================================================
*/

workflow {
  println(__params)

  chr_ch     = makeChrChannel(__params.reference)
  regions_ch = makeWindowChannel(chr_ch, __params.window_size)

  if (__params.reads_type == "SE"){
     reads = Channel.fromPath("${__params.input}/*.fastq.gz")
      .map { file -> def sampleID = file.name.replaceFirst(/\.fastq\.gz$/, "")
             [ sampleID, file ]}
     reads.view()
  }

  else if (__params.reads_type == "PE"){
     reads = Channel.fromFilePairs("${__params.input}/*R{1,2}.fastq.gz")
     .map {[it.first(), *it.last()]}
     reads.view()
  }

  (reads_trim,qc1,qc2) = QC(reads)
  (final_bam,qmap) = Alignment(reads_trim)
  MultiQC(qc1,qc2,qmap)
  (vcf, vcftbi) = Callvariant(final_bam,regions_ch)
  Postprocess(vcf, vcftbi)


//  (vcf, vcftbi) = CombineGVCF(gvcf, gvcftbi)
//  (bedd,bimm,famm,hmp) = VCF_Stats(vcf, vcftbi)

//  vcfgz = Channel.fromPath("${__params.input}/*g.vcf.gz")
//    .map { file ->
//        def sample = file.baseName.tokenize('.')[0]
//        return [sample, file]
//     }

}

workflow.onComplete { handleCompleteMessage() }
workflow.onError { handleErrorMessage() }
