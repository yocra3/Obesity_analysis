/*
 * Prepare data for imputation in Michigan from FullDataTable
 */

// Create date variable to avoid overwriting or merging results
date = java.time.LocalDate.now()

HRCref = file(params.HRCref)

// Add workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""


if (params.inputTable) {
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true)]}
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .set { ch_FullDataTable }
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[2], checkIfExists: true)]}
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .set { ch_Annot }

} else { exit 1, "--inputTable should be text file with the cohort name and a path to FullDataTable." }

// Adapt annotation
process adaptAnnotation{

  tag "$cohort"

  input:
  set val(cohort), file(annot) from ch_Annot

  output:
  set val(cohort), file("map.txt") into ch_map

  script:
  """
  awk -F',' 'NF > 10 {print \$2, "A", "B", \$4}' $annot | tail -n +2 > map.txt ## Select annot fields
  sed -i 's/\\// /'g map.txt
  sed -i 's/\\[//'g map.txt
  sed -i 's/\\]//'g map.txt
  """

}

// Convert fullDataTable to tped format
process fulldatatableTotped {

  tag "$cohort"

  input:
  set val(cohort), file(geno) from ch_FullDataTable

  output:
  set val(cohort), file("cohort.tped"), file("cohort.tfam") into ch_geno_tped

  script:

  if( cohort == 'Obesity1' )
    ncol = 4
  else
    ncol = 3
  """
  awk -F'\\t' 'NR == 1 {for (i = 4; i <= NF; i += $ncol) print \$i;}' $geno > samps.txt ## Get sample names
  sed -i 's/.GType//g' samps.txt
  dups=\$(sort samps.txt | uniq -d | sort -nr | grep -w -n -f - samps.txt | cut -d : -f 1) ## Select duplicated samples

  awk -F'\\t' 'BEGIN { OFS = FS } NR > 1 {print \$2, \$1, "0", \$3;}' $geno > annot.txt ## Select annot fields
  awk -F'\\t' 'NR > 1 {for (i = 4; i <= NF; i += $ncol) printf ("%s%c", \$i, i + $ncol <= NF ? "\\t" : "\\n");}' $geno > genomat.txt ## Prepare genotypes

  if [ -z "\$dups" ]
  then
    mv genomat.txt genomat.filt
    mv samps.txt samps.filt
  else

    index=""
    ## Create vector of duplicated samples
    for i in \${dups[@]}
    do
      index=\$index,\$i
    done
    index=\${index:1}

    ## Remove duplicated samples from genomat
    cut --complement -f \$index genomat.txt > genomat.filt
    ## Remove duplicated samples from samps
    sort samps.txt | uniq -d | sort -nr | grep -w -v -f - samps.txt > samps.filt
  fi

  ## Transform genotypes to tped format
  sed -i 's/AA/A\tA/g' genomat.filt
  sed -i 's/AB/A\tB/g' genomat.filt
  sed -i 's/BB/B\tB/g' genomat.filt
  sed -i 's/NC/0\t0/g' genomat.filt
  paste annot.txt genomat.filt > cohort.tped


  awk -F'\\t' 'BEGIN { OFS = FS } {print \$1, \$1, 0, 0, 0, 0}' samps.filt > cohort.tfam



  """
}

// Convert tped to plink
process tpedToPlink {

  tag "$cohort"

  input:
  set val(cohort), file(tped), file(tfam) from ch_geno_tped

  output:
  set val(cohort), file("cohort.bed"), file("cohort.bim"), file("cohort.fam") into ch_plink_ini

  script:
  """
  plink --tfile cohort --make-bed --out cohort
  """
}
ch_plink_variants = ch_plink_ini.join(ch_map)

// Add correct allele and remove SNP with low call rate
process addAlleles {

  tag "$cohort"


  input:
  set val(cohort), file(bed), file(bim), file(fam), file(annot) from ch_plink_variants

  output:
  set val(cohort), file("allele.bed"), file("allele.bim"), file("allele.fam") into ch_plink_allele

  script:
  """
  plink --bfile cohort --update-alleles $annot --geno 0.1 --make-bed --out allele
  """
}


// Exclude bad samples
process excludeBadSamples {

  tag "$cohort"

  publishDir "${params.outdir}/preImputation/plink/${cohort}/${date}", mode: 'copy'

  input:
  set val(cohort), file(bed), file(bim), file(fam) from ch_plink_allele

  output:
  set val(cohort), file("${cohort}.bed"), file("${cohort}.bim"), file("${cohort}.fam") into ch_geno_plink, ch_geno_plink2
  file("${cohort}.log")

  script:
  """
  plink --bfile allele  --mind 0.1 --make-bed --out $cohort
  """
}

// Create freq file for HRC check script
process createFreqFile {

  tag "$cohort"

  input:
  set val(cohort), file(bed), file(bim), file(fam) from ch_geno_plink

  output:
  set val(cohort), file("${cohort}.frq") into ch_geno_freq

  script:
  """
  plink --freq --bfile $cohort --out $cohort
  """
}
ch_plink_merge = ch_geno_plink2.join(ch_geno_freq)

process correctGenotypes {

  tag "$cohort"
  publishDir "${params.outdir}/preImputation/HRC_check/${cohort}/${date}", mode: 'copy', pattern: 'LOG*'

  input:
  set val(cohort), file(bed), file(bim), file(fam), file(freq) from ch_plink_merge
  file(HRCref)

  output:
  set val(cohort), file("${cohort}-updated.bed"), file("${cohort}-updated.bim"), file("${cohort}-updated.fam") into ch_plink_final
  file("LOG*")

  script:
  """
  HRC-1000G-check-bim.pl -b $bim -f $freq -r $HRCref -h ## Run HRC check script
  sed -i 's,'"\$PWD/"',,g' Run-plink.sh
  bash ./Run-plink.sh
  """
}

process generateVCFs {

  tag "$cohort"

  publishDir "${params.outdir}/preImputation/VCF/${cohort}/${date}", mode: 'copy'


  input:
  set val(cohort), file(bed), file(bim), file(fam) from ch_plink_final

  output:
  set val(cohort), file("*.vcf.gz") into ch_vcf

  script:
  """
  for i in \$(seq 1 1 23)
  do
    plink --bfile ${cohort}-updated --real-ref-alleles --recode vcf-iid bgz --chr \$i --out ${cohort}.chr\${i}
  done
  """

}
