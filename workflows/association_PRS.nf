/*
 * Run Polygenic Risk Score analysis
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()

Channel.fromPath(["${params.gwas_vcf_file}", "${params.gwas_vcf_file}.tbi"]).into{ch_vcf1; ch_vcf2}

ch_pgs = Channel.of('PGS000027').concat(Channel.fromPath("${params.prs_file}"))
ld_folder = file("${params.ld_folder}")
map_file = file("${params.mapFile}")



if (params.inputTable) {
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true), row[2], row[3]]}
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .set { ch_prs_input}
  } else { exit 1, "--inputTable should be text file with the cohort name and a path to the folder with imputation files" }



process definePRS_R {

  tag "$dataset"
  label 'process_long'
  label 'process_high'

  publishDir "${params.outdir}/associations/PRS/${date}", mode: 'copy'

  input:
  set val(dataset), file(sumstats), val(n_cases), val(n_controls) from ch_prs_input
  path(ld_folder)
  file(map_file)

  output:
  set val(dataset), file("${dataset}.ldpred2_prs.txt.gz") into prs_r
  file("${dataset}.ldpred_auto.Rdata") into prs_rData

  """
  mkdir tmp-data
  definePRS_ldpred2.R $sumstats $n_cases $n_controls ${task.cpus}
  bgzip -c pgs.txt > ${dataset}.ldpred2_prs.txt.gz
  mv ldpred_auto.Rdata ${dataset}.ldpred_auto.Rdata
  """

}

ch_prs = ch_pgs.collect().concat(prs_r)

process processPRS {

  tag "$dataset"

  label 'process_long'

  input:
  set val(dataset),file(prs) from ch_prs
  set file(vcf), file(tbi) from ch_vcf1.collect()

  output:
  set val(dataset), file("prs.plink") into ch_prs_plink


  script:
  """
  bgzip -cd $prs | grep -v '^#' > prs.txt
  cut -f1,2 prs.txt | tail -n +2 > regs.txt
  bcftools query -f '%CHROM  %POS  %ID\n' -R regs.txt $vcf | awk '{print \$1":"\$2"\t"\$3}' | sort -k1,1 > impIds.txt
  awk '{print \$1":"\$2"\t"\$3"\t"\$5}' prs.txt | sort -k1,1 > prs.part
  join -j1 -o 2.2,1.2,1.3 prs.part impIds.txt > prs.plink
  """
}

process computePRS {

  tag "$dataset"

  publishDir "${params.outdir}/associations/PRS/${date}", mode: 'copy'

  input:
  set val(dataset), file(prs) from ch_prs_plink
  set file(vcf), file(tbi) from ch_vcf2.collect()

  output:
  file("${dataset}.profile") into ch_out

  script:
  """
  plink --vcf $vcf --const-fid --score $prs sum --out $dataset
  """
}
