/*
 * Prepare data from imputation to association
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()
sample_excel = file("${params.sample_excel}")
obe1_extra = file("${params.obesity1Extra}")

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


  label 'process_long'

  input:
  file(vcf) from ch_collapsed_chr_tbi.collect()

  output:
  file("merged.ini.vcf.gz") into ch_vcf_ini

  script:
  """
  bcftools merge --force-samples *.gz -m none -o merged.ini.vcf.gz -O z
  """
}


process getCommonVariants {


  label 'process_long'

  input:
  file(vcf) from ch_vcf_ini

  output:
  file("merged.imputation.vcf.gz") into ch_vcf_tabix, ch_merge_vcf

  script:
  """
  bcftools filter -i 'N_MISSING < 10' $vcf -o merged.imputation.vcf.gz -O z
  """
}

process indexVCFmerged {

  publishDir "${params.outdir}/postImputation/merged/${date}", mode: 'copy'

  input:
  file(vcf) from ch_vcf_tabix

  output:
  set file(vcf), file("${vcf}.tbi") into ch_vcf_peddy

  script:
  """
  tabix -p vcf $vcf
  """
}



process makePlink {

  label 'process_medium'

  input:
  file(vcf) from ch_merge_vcf
  file(sex) from ch_samps_sex

  output:
  set file("merged.bed"), file("merged.bim"), file("merged.fam") into ch_plink, ch_plink_peddy, ch_plink_ini

  script:
  """
  plink --vcf $vcf --const-fid --update-sex $sex --make-bed --out merged
  """
}


process plinkStatistics {

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  set file("plink.bed"), file("plink.bim"), file("plink.fam") from ch_plink

  output:
  file("stats*") into ch_plink_out
  set file("stats.genome.gz"), file("stats.sexcheck") into ch_stats, ch_sex_plink

  script:
  """
  plink --bfile plink --genome gz --min 0.1  --check-sex --out stats
  """
}


process runPeddy {

  label 'process_medium'

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  set file(bed), file(bim), file(fam) from ch_plink_peddy
  set file(vcf), file(tbi) from ch_vcf_peddy

  output:
  file("merged.het_check.csv") into ch_peddy_out

  script:
  """
  sed -i 's/ /\\t/g' $fam
  python -m peddy -p ${task.cpus} --prefix merged $vcf $fam
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

process combineBadSamples {

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  file(geno) from ch_bad_samps
  file(excel) from ch_exclude_samps
  file(ancestry) from ch_peddy_out
  file(extra) from obe1_extra

  output:
  file("bad.samples.combined.txt") into ch_bad_samps_comb_vcf
  file("bad.samples.plink.txt") into ch_bad_samps_comb_plink

  script:
  """
  ## Merge filter
  cat $geno >> $excel
  cat $extra >> $excel

  ## Filter non-european samples
  awk '\$12 != "EUR" || \$13 <= 0.9 {print \$1}' FS="," $ancestry | grep -v sample_id >> $excel
  sort -u $excel > bad.samples.combined.txt
  sed 's/^/0\t/g' bad.samples.combined.txt >  bad.samples.plink.txt
  """
}

process makeFilteredVCF {

  publishDir "${params.outdir}/postImputation/merged/${date}", mode: 'copy'

  label 'process_long'

  input:
  file(vcf) from ch_merge_vcf
  file(samps) from ch_bad_samps_comb_vcf


  output:
  file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz")
  file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz.tbi")

  script:
  """
  bcftools view -S ^$samps --force-samples $vcf -o ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz -O z
  tabix -p vcf ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz
  """
}


process makePlinkFinal {

  label 'process_medium'
  publishDir "${params.outdir}/postImputation/merged/${date}", mode: 'copy'

  input:
  set file(bed), file(bim), file(fam) from ch_plink_ini
  set file(family), file(sex) from ch_sex_plink
  file(samp) from ch_bad_samps_comb_plink

  output:
  set file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.bed"), file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.bim"), file("ObesityGWAS.imputatedArrays.mergedDataset.filtered.fam") into ch_plink_filtered

  script:
  """
  awk '{print \$1, \$2, \$4}' $sex > plink.sex
  grep -v SB $fam | cut -d ' ' -f1,2 - > cases.txt
  plink --bfile merged --remove $samp --update-sex plink.sex --make-pheno cases.txt '*' --make-bed --out ObesityGWAS.imputatedArrays.mergedDataset.filtered
  """
}

process computePCA {

  label 'process_medium'
  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  set file("data.bed"), file("data.bim"), file("data.fam") from ch_plink_filtered

  output:
  file("mergedDataset.filtered.eigenval")
  file("mergedDataset.filtered.eigenvec")

  script:
  """
  plink --bfile data --pca 10 --out mergedDataset.filtered
  """
}
