/*
 * Run association analysis
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()

plink = [file("${params.plink_prefix}.bed"), file("${params.plink_prefix}.bim"), file("${params.plink_prefix}.fam")]
Channel.fromPath(plink).into{ch_plink_filter; ch_plink_freqs; ch_plink_assoc; ch_plink_depict; ch_plink_assoc_cohort; ch_plink_patho; ch_plink_raw}
Channel.fromPath(["${params.gwas_vcf_file}", "${params.gwas_vcf_file}.tbi"]).into{ch_vcf1; ch_vcf2; ch_vcf3; ch_vcfburden; ch_vcf_inv}

pcs = file("${params.pcs_file}")
prs = file("${params.prs_file}")
phenotype = plink[2]

samples = file("${params.sample_excel}")

Channel.fromPath(["${params.snps_1000G}", "${params.snps_1000G}.tbi"]).set{snps_1000G}
flat_gene = file("${params.flat_gene}")
depictFold = file("${params.depictFold}")
depictFoldData = file("${params.depictFold}/data/")
depictFoldLocus = file("${params.depictFold}/LocusGenerator/")
depictFoldDepict = file("${params.depictFold}/Depict/")

if (params.cohortTable) {
      Channel
            .fromPath(params.cohortTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true)]}
            .ifEmpty { exit 1, "params.cohortTable was empty - no input files supplied" }
            .branch{
              cases: it[0] != "Epicuro"
              controls: it[0] == "Epicuro"
            }
            .set { ch_cohort_fams }
  } else { exit 1, "--cohortTable should be text file with the cohort name and a path to the folder with imputation files" }


// Add workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

process createCovars {

  input:
  file(fam) from phenotype
  file(pcs) from pcs

  output:
  file("covars.txt") into ch_covars, ch_covars_inv, ch_covars_cohorts, ch_covars_patho

  script:
  """
  a="FID IID FATID MATID"
  echo \$a {A..J} > covars.txt
  join -j 2 -o 1.1,1.2,1.3,1.4,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12 $fam $pcs >> covars.txt
  """
}


process runSingleGWASplink {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  label 'process_medium'
  label 'process_long'

  input:
  set file("set.bed"), file("set.bim"), file("set.fam") from ch_plink_assoc.collect()
  file(covars) from ch_covars

  output:
  file("Obesity.GWAS.plink.assoc.logistic") into plink_assoc

  script:
  """
  cut -d ' ' -f 1,2,5-13 $covars > cov.txt
  plink --bfile set --logistic --covar cov.txt --out Obesity.GWAS.plink
  """
}


process filterGWASplink {

  input:
  file(assoc) from plink_assoc

  output:
  file("Obesity.GWAS.plink.assoc.SNPs.logistic") into plink_assoc_snps, plink_assoc_snps_depict, plink_assoc_comp

  script:
  """
  head -n1 $assoc > Obesity.GWAS.plink.assoc.SNPs.logistic
  grep ADD $assoc >> Obesity.GWAS.plink.assoc.SNPs.logistic
  """
}


process compressSingleGWASplink {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  input:
  file(assoc) from plink_assoc_comp

  output:
  file("Obesity.GWAS.locuszoom.input.gz")

  script:
  """
  tail -n +2 $assoc | tr -s ' ' '\t' | cut -f 2-10  > tmp.txt
  cut -f2 tmp.txt | awk -F ":" 'BEGIN {OFS = "\t"} {print \$3, \$4}' > alleles.txt
  paste tmp.txt alleles.txt | bgzip -c  > Obesity.GWAS.locuszoom.input.gz
  """
}

// Based on PMID: 29358691
process defineLoci {

  label 'process_medium'
  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  input:
  set file("set.bed"), file("set.bim"), file("set.fam") from ch_plink_filter.collect()
  file(assoc) from plink_assoc_snps

  output:
  file("Obesity.GWAS.plink.clumped") into ch_loci

  script:
  """
  plink --bfile set --clump $assoc --clump-p1 5e-8 --clump-p2 1e-5 --clump-kb 250 -out Obesity.GWAS.plink
  """
}

pval_threshold = ["assoc", "signif"]

process extractDepictSNPs {

  input:
  file(assoc) from plink_assoc_snps_depict
  set file(vcf), file(tbi) from snps_1000G.collect()
  val(type) from pval_threshold


  output:
  set val(type), file("*.snp.ids") into depict_snps

  script:

  def thres = "${type}" == "1e-5" ? 1e-5 : 5e-8
  """
  LC_ALL=C awk '\$9 < $thres  {print \$2}' $assoc | grep ':' > snp.list
  awk -F ":" 'BEGIN {OFS = "\t"} {print \$1, \$2}' snp.list > snp.pos
  bcftools query -f '%ID\n' -R snp.pos $vcf > ${type}.snp.ids
  """
}


process runDepict {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  input:
  set val(type), file(snps) from depict_snps
  path(depictFoldData)
  path(depictFoldLocus)
  path(depictFoldDepict)

  output:
  file("*loci.txt")

  script:
  """
  mkdir results
  depict.py $type $snps
  mv results/*.* .
  """
}

process splitEpicuro {

  input:
  set val(cohort), file(fam) from ch_cohort_fams.controls


  output:
  set val("Obesity1"), file("subsetaa") into ch_epicuro_sets1
  set val("Obesity2"), file("subsetab") into ch_epicuro_sets2
  set val("Obesity3"), file("subsetac") into ch_epicuro_sets3

  script:
  """
  split -l\$((`wc -l < $fam`/3 + 1)) $fam subset
  """
}
ch_epicuro_sets = ch_epicuro_sets1.concat(ch_epicuro_sets2.concat(ch_epicuro_sets3))
ch_cohorts = ch_epicuro_sets.join(ch_cohort_fams.cases)


process runCohortGWASplink {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  label 'process_medium'
  label 'process_long'

  input:
  set file("set.bed"), file("set.bim"), file("set.fam") from ch_plink_assoc_cohort.collect()
  set val(cohort), file(fam1), file(fam2) from ch_cohorts
  file(covars) from ch_covars_cohorts

  output:
  set val(cohort), file("Obesity.GWAS.plink.${cohort}.assoc.logistic") into plink_assoc_cohort

  script:
  """
  cut -d ' ' -f 1,2,5-13 $covars > cov.txt
  cat $fam1 > selsamps.txt
  cat $fam2 >> selsamps.txt
  cut -d ' ' -f 2 selsamps.txt | sed 's/^/0\t/g'  >  selsamps
  plink --bfile set --keep selsamps --logistic --covar cov.txt --out Obesity.GWAS.plink.${cohort}
  """
}


process filterCohortGWASplink {

  input:
  set val(cohort), file(assoc) from plink_assoc_cohort

  output:
  set val(cohort), file("Obesity.GWAS.plink.${cohort}.assoc.SNPs.logistic") into plink_assoc_cohort_comp

  script:
  """
  head -n1 $assoc > Obesity.GWAS.plink.${cohort}.assoc.SNPs.logistic
  grep ADD $assoc >> Obesity.GWAS.plink.${cohort}.assoc.SNPs.logistic
  """
}


process compressCohortGWASplink {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  input:
  set val(cohort), file(assoc) from plink_assoc_cohort_comp

  output:
  file("Obesity.GWAS.${cohort}.locuszoom.input.gz")

  script:
  """
  tail -n +2 $assoc | tr -s ' ' '\t' | cut -f 2-10  > tmp.txt
  cut -f2 tmp.txt | awk -F ":" 'BEGIN {OFS = "\t"} {print \$3, \$4}' > alleles.txt
  paste tmp.txt alleles.txt | bgzip -c  > Obesity.GWAS.${cohort}.locuszoom.input.gz
  """
}


process computeFreqs {

  input:
  set file("set.bed"), file("set.bim"), file("set.fam") from ch_plink_freqs.collect()

  output:
  file("set.frq") into ch_freq

  script:
  """
  plink --bfile set --freq -out set
  """
}

process callInversions {

  label 'process_medium'

  publishDir "${params.outdir}/associations/inversions/${date}", mode: 'copy'

  input:
  set file(vcf), file(tbi) from ch_vcf_inv.collect()
  file(freq) from ch_freq

  output:
  file("ObesityscoreInvHapClassDF.Rdata") into ch_inv
  file("Obesity.scoreInvHaplist.Rdata")

  script:
  """
  callInversions.R $vcf $freq Obesity
  """
}

process assocInversions {

  publishDir "${params.outdir}/associations/inversions/${date}", mode: 'copy'

  input:
  file(invs) from ch_inv
  file(fam) from phenotype
  file(covars) from ch_covars_inv


  output:
  file("Obesity.inversions.assoc.Rdata")

  """
  assocInversions.R $invs $fam $covars
  """
}



process processPRS {

  label 'process_long'

  input:
  file(prs) from prs
  set file(vcf), file(tbi) from ch_vcf2.collect()

  output:
  file("prs.plink") into ch_prs_plink


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

  publishDir "${params.outdir}/associations/PRS/${date}", mode: 'copy'

  input:
  file(prs) from ch_prs_plink
  set file(vcf), file(tbi) from ch_vcf3.collect()

  output:
  file("obesity.profile") into ch_prs

  script:
  """
  plink --vcf $vcf --const-fid --score $prs sum --out obesity
  """
}

process getPathoSamples {

  input:
  file(samples)

  output:
  file("pathogenicSamps.txt") into ch_path_samps

  script:
  """
  selectSamplesWithVariants.R $samples
  """
}

process runGWASplinkPatho {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  label 'process_medium'
  label 'process_long'

  input:
  set file("set.bed"), file("set.bim"), file("set.fam") from ch_plink_patho.collect()
  file(covars) from ch_covars_patho
  file(bad_samps) from ch_path_samps

  output:
  set val("non_patho"), file("Obesity.GWAS.plink.non_patho.assoc.logistic") into plink_assoc_patho

  script:
  """
  cut -d ' ' -f 1,2,5-13 $covars > cov.txt
  plink --bfile set --logistic --remove $bad_samps --out Obesity.GWAS.plink.non_patho
  """
}

process runGWASplinkRaw {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  label 'process_medium'
  label 'process_long'

  input:
  set file("set.bed"), file("set.bim"), file("set.fam") from ch_plink_raw.collect()

  output:
  set val("raw"), file("Obesity.GWAS.plink.raw.assoc.logistic") into plink_assoc_raw

  script:
  """
  plink --bfile set --logistic --out Obesity.GWAS.plink.raw
  """
}

ch_gwas_comb = plink_assoc_patho.concat(plink_assoc_raw)

process filterGWASpatho {

  input:
  set val(subset), file(assoc) from ch_gwas_comb

  output:
  set val(subset), file("Obesity.GWAS.plink.${subset}.assoc.SNPs.logistic") into plink_assoc_patho_comp

  script:
  """
  head -n1 $assoc > Obesity.GWAS.plink.${subset}.assoc.SNPs.logistic
  grep ADD $assoc >> Obesity.GWAS.plink.${subset}.assoc.SNPs.logistic
  """
}


process compressGWASpatho {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  input:
  set val(subset), file(assoc) from plink_assoc_patho_comp

  output:
  file("Obesity.GWAS.plink.${subset}.locuszoom.input.gz")

  script:
  """
  tail -n +2 $assoc | tr -s ' ' '\t' | cut -f 2-10  > tmp.txt
  cut -f2 tmp.txt | awk -F ":" 'BEGIN {OFS = "\t"} {print \$3, \$4}' > alleles.txt
  paste tmp.txt alleles.txt | bgzip -c  > Obesity.GWAS.plink.${subset}.locuszoom.input.gz
  """
}
