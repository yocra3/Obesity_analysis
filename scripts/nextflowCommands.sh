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
--gwas_vcf_file results/postImputation/merged/2021-06-04/ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz \
--plink_prefix results/postImputation/merged/2021-06-04/ObesityGWAS.imputatedArrays.mergedDataset.filtered \
--cohortTable data/inputAssociation.tsv \
--pcs_file results/postImputation/stats/2021-06-04/mergedDataset.filtered.eigenvec -profile docker -resume

## Prepare VLF files for imputation
nextflow run workflows/imputeData_VLF.nf  --plink_prefix data/VLF  -profile docker -resume -with-tower

## Get files postImputation
nextflow run workflows/postImputationQC_VLF.nf --imputation_folder results/postImputation/MichiganOutput/VLF/ \
--fam_file  results/preImputation/plink/VLF/2021-04-30/VLF_filt.fam \
--pheno_file data/VLF.obesity.txt -profile docker -resume -with-tower

## Run PGS analysis
nextflow run workflows/association_PRS.nf --prs_file data/PGS000027.txt.gz \
--gwas_vcf_file results/postImputation/merged/2021-06-04/ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz \
--inputTable data/inputPGS.tab -profile docker -resume
