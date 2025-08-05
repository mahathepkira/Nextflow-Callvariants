include {outputPrefixPath } from './nbt/utils'

process FastQC_visualize {


  publishDir "${outputPrefixPath(params, task)}"
//  publishDir "${s3OutputPrefixPath(params, task)}"

  input:
  path zip

  output:
  path "summary.csv"

  script:

  """
  
  for file in *.zip; do unzip -o "\$file" -d extracted/; done
  echo "file_name,Total_Sequences_BeforeQC,Total_Sequences_AfterQC,Sequence_length_BeforeQC,Sequence_length_AfterQC,%GC_BeforeQC,%GC_AfterQC" > summary.csv
  declare -A before_qc after_qc

  for file in extracted/*/fastqc_data.txt; do
      filename=\$(grep "Filename" "\$file" | cut -f2)
      total_sequences=\$(grep "Total Sequences" "\$file" | cut -f2)
      seq_length=\$(grep "Sequence length" "\$file" | cut -f2)
      gc_content=\$(grep "%GC" "\$file" | cut -f2)

      base_name=\$(echo "\$filename" | sed -E 's/_R[12]_.+//')

    
     if [[ "\${filename}" != "\${filename%_paired.fastq.gz}" ]]; then
          before_qc["\$base_name"]="\$total_sequences,\$seq_length,\$gc_content"
          echo "Filename: \$filename"
     else
          after_qc["\$base_name"]="\$total_sequences,\$seq_length,\$gc_content"
     fi

  done

  for key in "\${!before_qc[@]}"; do
      before="\${before_qc[\$key]:-,,}"
      after="\${after_qc[\$key]:-,,}"

      before_total_sequences=\$(echo \$after | cut -d',' -f1)
      before_seq_length=\$(echo \$after | cut -d',' -f2)
      before_gc_content=\$(echo \$after | cut -d',' -f3)

      after_total_sequences=\$(echo \$before | cut -d',' -f1)
      after_seq_length=\$(echo \$before | cut -d',' -f2)
      after_gc_content=\$(echo \$before | cut -d',' -f3)

      echo "\$key,\$before_total_sequences,\$after_total_sequences,\$before_seq_length,\$after_seq_length,\$before_gc_content,\$after_gc_content" >> summary.csv

  done
  
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

