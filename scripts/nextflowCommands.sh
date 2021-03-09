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
