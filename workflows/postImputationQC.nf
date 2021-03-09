/*
 * Prepare data from imputation to association
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()
sample_excel = file("${params.sample_excel}")

if (params.inputTable) {
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true)]}
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .set { ch_Imputation_Folder }
  } else { exit 1, "--inputTable should be text file with the cohort name and a path to the folder with imputation files" }

// Add workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

process processSampleSheet {

  publishDir "${params.outdir}/phenotypes/${date}", mode: 'copy'

  input:
  file(sampleSheet) from sample_excel

  output:
  file("samples.sex.txt") into ch_samps_sex
  file("bad.samples") into ch_exclude_samps
  file("bad.samples.txt")

  script:
  """
  processSampleSheet.R $sampleSheet
  """
}


process collapseChromosomes {

  tag "$cohort"
  label 'process_medium'
  label 'process_long'

  input:
  set val(cohort), path("vcfs") from ch_Imputation_Folder

  output:
  set val(cohort), file("${cohort}.imputation.vcf.gz") into ch_collapsed_chr

  script:
  """
  bcftools concat vcfs/*.dose.vcf.gz --threads ${task.cpus} | bcftools sort -o ${cohort}.imputation.vcf.gz -O z
  """

}

process indexVCF {

  tag "$cohort"

  input:
  set val(cohort), file(vcf) from ch_collapsed_chr

  output:
  set val(cohort), file(vcf), file("${vcf}.tbi") into ch_collapsed_chr_tbi

  script:
  """
  tabix -p vcf $vcf
  """
}

process mergeVCFs {

  publishDir "${params.outdir}/postImputation/merged/${date}", mode: 'copy'

  label 'process_long'

  input:
  file(vcf) from ch_collapsed_chr_tbi.collect()

  output:
  file("merged.imputation.vcf.gz") into ch_merge_vcf

  script:
  """
  bcftools merge --force-samples *.gz -m none -o merged.imputation.vcf.gz -O z
  """
}


process plinkStatistics {

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  file(vcf) from ch_merge_vcf
  file(sex) from ch_samps_sex

  output:
  file("stats*") into ch_plink_out
  set file("stats.genome.gz"), file("stats.sexcheck") into ch_stats

  script:
  """
  plink --vcf $vcf --const-fid --update-sex $sex --genome gz --min 0.1  --check-sex --pca --out stats
  """
}


process makeBadSamplesList {

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  set file(family), file(sex) from ch_stats

  output:
  file("exclude.samps") into ch_bad_samps
  file("exclude.samps.log")
  file("exclude.samps.reason.txt")

  script:
  """
  selectBadSamples.R $family $sex
  """
}


process makeFilteredVCF {

  publishDir "${params.outdir}/postImputation/merged/${date}", mode: 'copy'

  label 'process_long'

  input:
  file(vcf) from ch_merge_vcf
  file(geno) from ch_bad_samps
  file(excel) from ch_exclude_samps

  output:
  file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz")
  file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz.tbi")

  script:
  """
  ## Merge filter
  cat $geno >> $excel

  bcftools view -S ^$excel --force-samples $vcf -o ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz -O z
  tabix -p vcf ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz
  """
}
