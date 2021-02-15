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
} else { exit 1, "--inputTable should be text file with the cohort name and a path to FullDataTable." }


// Convert fullDataTable to tped format
process fulldatatableTotped {

  tag "$cohort"

  input:
  set val(cohort), file(geno) from ch_FullDataTable

  output:
  set val(cohort), file("cohort.tped"), file("cohort.tfam") into ch_geno_tped

  script:
  """
  awk -F'\\t' 'BEGIN { OFS = FS } NR > 1 {print \$2, \$1, "0", \$3;}' $geno > annot.txt ## Select annot fields
  awk -F'\\t' 'NR > 1 {for (i = 4; i <= NF; i += 3) printf ("%s%c", \$i, i + 3 <= NF ? "\\t" : "\\n");}' $geno > genomat.txt ## Prepare genotypes
  ## Transform genotypes to tped format
  sed -i 's/AA/1\t1/g' genomat.txt
  sed -i 's/AB/1\t2/g' genomat.txt
  sed -i 's/BB/2\t2/g' genomat.txt
  sed -i 's/NC/0\t0/g' genomat.txt
  paste annot.txt genomat.txt > cohort.tped

  awk -F'\\t' 'NR == 1 {for (i = 4; i <= NF; i += 3) print \$i;}' $geno > samps.txt ## Get sample names
  sed -i 's/.GType//g' samps.txt
  awk -F'\\t' 'BEGIN { OFS = FS } {print \$1, \$1, 0, 0, 0, 0}' samps.txt > cohort.tfam
  """
}

// Convert tped to plink
process tpedToPlink {

  tag "$cohort"

  input:
  set val(cohort), file(tped), file(tfam) from ch_geno_tped

  output:
  set val(cohort), file("${cohort}.bed"), file("${cohort}.bim"), file("${cohort}.fam") into ch_geno_plink, ch_geno_plink2

  script:
  """
  plink --tfile cohort --make-bed --out $cohort
  """
}

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

process checkGenotypes {

  tag "$cohort"

  input:
  set val(cohort), file(bed), file(bim), file(fam), file(freq) from ch_plink_merge
  file(HRCref)

  script:
  """
  perl HRC-1000G-check-bim.pl -b $bim -f $freq -r $HRCref -h
  """
}
