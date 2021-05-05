/*
 * Prepare data from imputation to association
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()

ch_imputation_folder = Channel.fromPath("${params.imputation_folder}")
ch_fam_file = Channel.fromPath("${params.fam_file}")
ch_pheno = Channel.fromPath("${params.pheno_file}")

process collapseChromosomes {

  label 'process_medium'
  label 'process_long'

  input:
  path("vcfs") from ch_imputation_folder

  output:
  file("VLF.imputation.vcf.gz") into ch_collapsed_chr, ch_merge_vcf

  script:
  """
  bcftools concat vcfs/*.dose.vcf.gz --threads ${task.cpus} | bcftools sort -o VLF.imputation.vcf.gz -O z
  """

}

process indexVCF {
  publishDir "${params.outdir}/postImputation/VLF/${date}", mode: 'copy'

  input:
  file(vcf) from ch_collapsed_chr

  output:
  set file(vcf), file("${vcf}.tbi") into ch_collapsed_chr_tbi

  script:
  """
  tabix -p vcf $vcf
  """
}

process makePlink {

  label 'process_medium'

  input:
  set file(vcf), file(tbi) from ch_collapsed_chr_tbi
  file(fam) from ch_fam_file

  output:
  set file("VLF.bed"), file("VLF.bim"), file("VLF.fam") into ch_plink, ch_plink_ini

  script:
  """
  plink --vcf $vcf --const-fid --make-bed --out ini
  plink --bfile ini --fam $fam --make-bed --out VLF
  """
}


process plinkStatistics {

  label 'process_medium'

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  set file("plink.bed"), file("plink.bim"), file("plink.fam") from ch_plink

  output:
  file("VLF.stats*") into ch_plink_out
  file("VLF.stats.sexcheck") into ch_sex

  script:
  """
  plink --bfile plink --genome gz --min 0.1  --check-sex --out VLF.stats
  """
}

process getBadSamples {

  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  file(sex) from ch_sex

  output:
  file("VLF.bad.samples.vcf.txt") into ch_bad_samps_comb_vcf
  file("VLF.bad.samples.plink.txt") into ch_bad_samps_comb_plink

  script:
  """
  ## Filter samples with sex problems

  awk '(\$3 == 1 && \$4 != 1) || (\$3 == 2 && \$4 != 0) {print \$2}' $sex > VLF.bad.samples.vcf.txt
  awk '(\$3 == 1 && \$4 != 1) || (\$3 == 2 && \$4 != 0) {print \$1, \$2}' $sex > VLF.bad.samples.plink.txt
  """
}

process makeFilteredVCF {

  publishDir "${params.outdir}/postImputation/VLF/${date}", mode: 'copy'

  label 'process_long'

  input:
  file(vcf) from ch_merge_vcf
  file(samps) from ch_bad_samps_comb_vcf


  output:
  file("VLF.imputatedArrays.filtered.vcf.gz")
  file("VLF.imputatedArrays.filtered.vcf.gz.tbi")

  script:
  """
  bcftools view -S ^$samps --force-samples $vcf -o VLF.imputatedArrays.filtered.vcf.gz -O z
  tabix -p vcf VLF.imputatedArrays.filtered.vcf.gz
  """
}


process makePlinkFinal {

  label 'process_medium'
  publishDir "${params.outdir}/postImputation/VLF/${date}", mode: 'copy'

  input:
  set file(bed), file(bim), file(fam) from ch_plink_ini
  file(samp) from ch_bad_samps_comb_plink
  file(pheno) from ch_pheno

  output:
  set file("VLF.imputatedArrays.filtered.bed"), file("VLF.imputatedArrays.filtered.bim"), file("VLF.imputatedArrays.filtered.fam") into ch_plink_filtered, ch_pca

  script:
  """
  grep -v "#" $pheno > pheno.filt
  join -j2 -o 1.1,1.2,2.4 <(sort -k2 $fam) <(sort -k2 pheno.filt) > tmp.txt
  plink --bfile VLF --remove $samp --pheno tmp.txt --1 --make-bed --out VLF.imputatedArrays.filtered
  """
}


process computePCA {

  label 'process_medium'
  publishDir "${params.outdir}/postImputation/stats/${date}", mode: 'copy'

  input:
  set file("data.bed"), file("data.bim"), file("data.fam") from ch_pca

  output:
  file("VLF.eigenval")
  file("VLF.eigenvec")

  script:
  """
  plink --bfile data --pca 10 --out VLF
  """
}
