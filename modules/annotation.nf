
include {outputPrefixPath } from './nbt/utils'

process ANN_snpEff {

  tag { "${vcfgz}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(vcfgz)

  output:
  path "*"

  script:
  
  prefix=vcfgz.simpleName
  """
  
  snpEff -nodownload -c /nbt_main/home/lattapol/test_tools/snpEff/snpEff.config -v Manihot_esculenta ${vcfgz} | bgzip -c > ${prefix}.ann.vcf.gz  
  
  """
}


process FastQC_plot {


  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"


  input:
  path sum

  output:
  path "*"

  script:

  """ 
  python --version
  ml list
  #python /nbt_main/home/lattapol/mycassava/gwas-cassava2/hh.py
  python /nbt_main/home/lattapol/mycassava/gwas-cassava2/bin/plot_qc_graph.py
  """
}

process Qualimap_visualize {

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  path qmap

  output:
  path "*"

  script:
 
  """
  echo "filename,number_of_reads,number_of_mapped_reads(%),number_of_unmapped_reads(%)" > qualimap_summary.csv

  for dir in *_stats; do
      file="\$dir/genome_results.txt"
      if [[ -f "\$file" ]]; then
         name=\$(echo "\$dir" | sed 's/_stats//')
         total_reads=\$(grep "number of reads =" "\$file" | awk '{print \$NF}' | tr -d ',')
         mapped_percent=\$(grep "number of mapped reads =" "\$file" | awk -F '[()]' '{print \$2}' | tr -d '%')

         unmapped_percent=\$(echo "scale=2; (100 - \$mapped_percent)" | bc)

         echo "\$name,\$total_reads,\$mapped_percent,\$unmapped_percent" >> qualimap_summary.csv
      fi
  done

  #python /nbt_main/home/lattapol/mycassava/gwas-cassava2/bin/plot_qualimap_graph.py
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
  bash /nbt_main/home/lattapol/mycassava/gwas-cassava2/bin/quality2.sh ${vcf}  
  """
}


process Histogram {

  tag { "${frq}" }

  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  file(frq)

  output:
  file("*.csv")
 
  script:
  prefix=frq

  """
  python /nbt_main/home/lattapol/mycassava/gwas-cassava2/bin/create_dataframe.py ${frq}
  """


}
