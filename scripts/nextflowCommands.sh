#'#################################################################################
#'#################################################################################
# Nextflow commands
#'#################################################################################
#'#################################################################################
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAzMzE2fS5kNjdiNWM3OGVkYWUyYWE3ZjM1OWYxOGNlYTU2NDBhYmMzNjdlODY2

## Prepare genotype files for imputation
nextflow run workflows/imputeData.nf  --inputTable data/inputImputation.tsv  -profile docker -resume -with-tower

## Get files postImputation
nextflow run workflows/postImputationQC.nf --inputTable data/inputImputationQC.tsv -profile docker -resume -with-tower

## Run associations
nextflow run workflows/associations.nf --prs_file data/PGS000027.txt.gz \
--gwas_vcf_file results/postImputation/merged/2021-04-06/ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz \
--plink_prefix results/postImputation/merged/2021-04-06/ObesityGWAS.imputatedArrays.mergedDataset.filtered \
--pcs_file results/postImputation/stats/2021-04-06/mergedDataset.filtered.eigenvec -profile docker -resume
