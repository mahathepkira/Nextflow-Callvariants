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
include { Alignment_bwa }  from './modules/alignment_bwa' params(__params)
include { Sam_view }  from './modules/alignment_bwa' params(__params)
include { Sort_bam }  from './modules/alignment_bwa' params(__params)
include { Mark_duplicates }  from './modules/alignment_bwa' params(__params)
include { Base_recalibrator }  from './modules/alignment_bwa' params(__params)
include { Call_GVCF }  from './modules/alignment_bwa' params(__params)
include { Combine_GVCF }  from './modules/alignment_bwa' params(__params)
include { VcftoBed }  from './modules/alignment_bwa' params(__params)
include { Trimmmomatic } from './modules/quality' params(__params)
include { FastQC as FastQC_before_trim } from './modules/quality' params(__params)
include { FastQC as FastQC_after_trim } from './modules/quality' params(__params)
include { FastQC as FastQC_unpaired } from './modules/quality' params(__params)
include { Qualimap as Qualimap_before } from './modules/quality' params(__params)
include { Qualimap as Qualimap_after } from './modules/quality' params(__params)
include { VCFtools_stats as VCFtools_before } from './modules/quality' params(__params)
include { VCFtools_stats as VCFtools_after } from './modules/quality' params(__params)
include { Histogram as Histogram_before } from './modules/quality' params(__params)
include { Histogram as Histogram_after } from './modules/quality' params(__params)
include { FastQC_visualize } from './modules/visualize' params(__params)
include { Qualimap_visualize } from './modules/visualize' params(__params)

workflow Callvariantout {
   take:
     //pairs
     vcfgz
   main:
     //report_qc1 = FastQC_before_trim(pairs)
     //(paired, unpaired) = Trimmmomatic(pairs)
     //report_qc2 = FastQC_after_trim(paired)
     //qc1 = report_qc1.map {[it[1]]}.flatten().toList()
     //qc2 = report_qc2.map {[it[1]]}.flatten().toList()
     //combined_qc = qc1+qc2
     //sum_qc = FastQC_visualize(combined_qc)
     //report_qc3 = FastQC_unpaired(unpaired)   
     //sortbam = Alignment_bwa(paired)
     //qmap1 = Qualimap_before(sortbam).toList() 
     //(markDup, dup_metrics)=Mark_duplicates(sortbam)
     //base_recalibrator = Base_recalibrator(markDup)
     //qmap2 = Qualimap_after(base_recalibrator).toList()
     //combine_qmap = qmap1+qmap2
     //Qualimap_visualize(combine_qmap) 
     //(vcfgz, vcfgztbi) = Call_GVCF(base_recalibrator)
     groupedGVCFs = groupTupleWithoutCommonKey(vcfgz, true)
     combine_gvcf=Combine_GVCF(groupedGVCFs)
     (frq1,lmiss1,tstv1,sum1) = VCFtools_before(combine_gvcf)
     Histogram_before(frq1,lmiss1)
     (bedd,bimm,famm,hmp,alle_bi)=VcftoBed(combine_gvcf)
     (frq2,lmiss2,tstv2,sum2) = VCFtools_after(alle_bi)
     Histogram_after(frq2,lmiss2)
  
   emit:
     //vcfgz
     bedd
     bimm
     famm
     hmp
}


/*
================================================================================
=                           Sinonkt Style M A I N                              =
================================================================================
*/

workflow {
  println(__params)

//  pairs = Channel.fromFilePairs("${__params.input}/*_R{1,2}_001.fastq.gz")
//  .map {[it.first(), *it.last()]}
//  pairs.view()
   
//  Callvariantout(pairs)

  vcfgz = Channel.fromPath("${__params.input}/*g.vcf.gz")
    .map { file ->
        def sample = file.baseName.tokenize('.')[0]
        return [sample, file]
     }
  
  vcfgz.view()

  Callvariantout(vcfgz)

}

workflow.onComplete { handleCompleteMessage() }
workflow.onError { handleErrorMessage() }
