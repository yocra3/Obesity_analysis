/*
 * Run association analysis
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()

plink = [file("${params.plink_prefix}.bed"), file("${params.plink_prefix}.bim"), file("${params.plink_prefix}.fam")]
Channel.fromPath(plink).into{ch_plink_filter; ch_plink_freqs; ch_plink_assoc; ch_plink_depict; ch_plink_assoc_prs}
Channel.fromPath(["${params.gwas_vcf_file}", "${params.gwas_vcf_file}.tbi"]).into{ch_vcf1; ch_vcf2; ch_vcf3; ch_vcfburden; ch_vcf_inv}

pcs = file("${params.pcs_file}")
prs = file("${params.prs_file}")
phenotype = plink[2]


Channel.fromPath(["${params.snps_1000G}", "${params.snps_1000G}.tbi"]).set{snps_1000G}
flat_gene = file("${params.flat_gene}")
depictFold = file("${params.depictFold}")
depictFoldData = file("${params.depictFold}/data/")
depictFoldLocus = file("${params.depictFold}/LocusGenerator/")
depictFoldDepict = file("${params.depictFold}/Depict/")

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
  file("covars.txt") into ch_covars, ch_covars_burden, ch_covars_inv, ch_covars_prs

  script:
  """
  a="FID IID FATID MATID"
  echo \$a {A..J} > covars.txt
  join -j 2 -o 1.1,1.2,1.3,1.4,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12 $fam $pcs >> covars.txt
  """
}

/*
process runSingleGWAS {

  publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'

  label 'process_medium'
  label 'process_long'

  input:
  file(fam) from phenotype
  file(covars) from ch_covars
  set file(vcf), file(tbi) from ch_vcf1.collect()

  output:
  file("Obesity.GWAS.SingleWald.assoc.gz") into single_assoc

  script:
  """
  rvtest --inVcf $vcf --pheno $fam --out Obesity.GWAS --single wald \
  --covar $covars --covar-name A,B,C,D,E,F,G,H,I,J --noweb --sex --numThread ${task.cpus}
  bgzip Obesity.GWAS.SingleWald.assoc
  """
}
*/

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
  cut -d ' ' -f 1,2,5-14 $covars > cov.txt
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






// process runBurdenGWAS {
//
//   publishDir "${params.outdir}/associations/GWAS/${date}", mode: 'copy'
//
//   label 'process_medium'
//   label 'process_long'
//
//   input:
//   file(fam) from phenotype
//   file(covars) from ch_covars_burden
//   set file(vcf), file(tbi) from ch_vcfburden.collect()
//   file(gene) from flat_gene
//
//   output:
//   file("Obesity.GWAS.CMCWald.assoc")
//
//   script:
//   """
//   rvtest --inVcf $vcf --pheno $fam --out Obesity.GWAS --geneFile $flat_gene \
//   --burden cmcWald --covar $covars \
//   --covar-name A,B,C,D,E,F,G,H,I,J --noweb --sex --numThread ${task.cpus}
//   """
// }

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
