/*
 * Prepare data for imputation in Michigan from FullDataTable
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()

ref1000g = file(params.ref1000g)

// Add workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

plink = [file("${params.plink_prefix}.bed"), file("${params.plink_prefix}.bim"), file("${params.plink_prefix}.fam")]
ch_plink = Channel.fromPath(plink)

chain = file("${params.chain}")
liftover = file("${params.liftover}")

// Convert coordinates to hg19
process converthg19 {

  input:
  set file("VLF.bed"), file("VLF.bim"), file("VLF.fam") from ch_plink.collect()
  file(chain)
  file(liftover)

  output:
  set file("VLF_19.bed"), file("VLF_19.bim"), file("VLF_19.fam") into ch_plink_19

  script:
  """
  plink --bfile VLF --output-chr MT --recode ped --out VLF_18
  liftOverPlink.py --map VLF_18.map --out VLF_lifted --chain $chain -e ./$liftover
  cut -f 4 VLF_lifted.bed.unlifted | sed "/^#/d" > VLF_to_exclude.dat
  plink --file VLF_18 --recode ped --out VLF_18_int --exclude VLF_to_exclude.dat
  plink --ped VLF_18_int.ped --map VLF_lifted.map --make-bed --out VLF_19
  """
}



// Add correct allele and remove SNP with low call rate
process excludeBadSNPs {

  input:
  set file("VLF.bed"), file("VLF.bim"), file("VLF.fam") from ch_plink_19.collect()

  output:
  set file("allele.bed"), file("allele.bim"), file("allele.fam") into ch_plink_allele

  script:
  """
  plink --bfile VLF --geno 0.1 --make-bed --out allele
  """
}


// Exclude bad samples
process excludeBadSamples {

  publishDir "${params.outdir}/preImputation/plink/VLF/${date}", mode: 'copy'

  input:
  set file(bed), file(bim), file(fam) from ch_plink_allele

  output:
  set file("VLF_filt.bed"), file("VLF_filt.bim"), file("VLF_filt.fam") into ch_geno_plink, ch_geno_plink2
  file("VLF_filt.log")

  script:
  """
  plink --bfile allele  --mind 0.1 --make-bed --out VLF_filt
  """
}

// Create freq file for HRC check script
process createFreqFile {


  input:
  set file(bed), file(bim), file(fam) from ch_geno_plink

  output:
  file("VLF_filt.frq") into ch_geno_freq

  script:
  """
  plink --freq --bfile VLF_filt  --nonfounders --out VLF_filt
  """
}
ch_plink_merge = ch_geno_plink2.concat(ch_geno_freq)

process correctGenotypes {

  label 'process_medium'

  publishDir "${params.outdir}/preImputation/HRC_check/VLF/${date}", mode: 'copy', pattern: 'LOG*'

  input:
  set file(bed), file(bim), file(fam), file(freq) from ch_plink_merge.collect()
  file(ref1000g)

  output:
  set file("VLF_filt-updated.bed"), file("VLF_filt-updated.bim"), file("VLF_filt-updated.fam") into ch_plink_final
  file("LOG*")

  script:
  """
  HRC-1000G-check-bim.pl -b $bim -f $freq -r $ref1000g -g -p AMR ## Run HRC check script
  sed -i 's,'"\$PWD/"',,g' Run-plink.sh
  bash ./Run-plink.sh
  """
}

process generateVCFs {

  publishDir "${params.outdir}/preImputation/VCF/VLF/${date}", mode: 'copy'


  input:
  set file(bed), file(bim), file(fam) from ch_plink_final

  output:
  file("*.vcf.gz") into ch_vcf

  script:
  """
  for i in \$(seq 1 1 23)
  do
    plink --bfile VLF_filt-updated --real-ref-alleles --recode vcf-iid bgz --chr \$i --out VLF_filt.chr\${i}
  done
  """

}
