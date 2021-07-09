#'#################################################################################
#'#################################################################################
# Set project folders' structure and link files
#'#################################################################################
#'#################################################################################

## Make folders
mkdir data
mkdir results
mkdir bin

## Link raw data
### Genotypes
ln -s /home/SHARED/CNIO//Epicuro/FullDataTable.txt data/Epicuro.FullDataTable.txt
ln -s /home/SHARED/CNIO/1a_cohort/FullDataTable_pennCNV.txt data/Obesity1.FullDataTable.txt
ln -s /home/SHARED/CNIO/2a_cohort/FullDataTable_pennCNV.txt data/Obesity2.FullDataTable.txt
ln -s /home/SHARED/CNIO/JArgente_EOO/pennCNV/FullDataTable_pennCNV.txt data/Obesity3.FullDataTable.txt

### Array annotation
ln -s /home/SHARED/CNIO/1a_cohort/HumanOmni1-Quad_v1-0_H.csv data/Obesity1.annot.csv
ln -s /home/SHARED/CNIO/2a_cohort/HumanOmniExpress-12-v1-0-K.csv data/Obesity2.annot.csv
ln -s /home/SHARED/CNIO/JArgente_EOO/InfiniumOmniExpress-24v1-3_A1.csv data/Obesity3.annot.csv
ln -s /home/SHARED/CNIO//Epicuro/Human1Mv1_C.csv data/Epicuro.annot.csv

## Phenotypes
ln -s /home/SHARED/DISK/CNIO/Taula\ resum\ pacients\ \(def\).xlsx data/Obesity.SamplesSummary.xlsx

### Viva la Familia
ln -s /media/Lacie_1/DATA/dbGaP/phs000616.v2.p2.VIVA_LA_FAMILIA/GenotypeFiles/phg000386.v1.NIDDK_VIVA_LA_FAMILIA-HumanOmni1_Quad_v1-0_B.genotype-calls-matrixfmt.c1.DS-OBS-IRB-RD.update/7_dbGaP_1Mdata_20130522.bed data/VLF.bed
ln -s /media/Lacie_1/DATA/dbGaP/phs000616.v2.p2.VIVA_LA_FAMILIA/GenotypeFiles/phg000386.v1.NIDDK_VIVA_LA_FAMILIA-HumanOmni1_Quad_v1-0_B.genotype-calls-matrixfmt.c1.DS-OBS-IRB-RD.update/7_dbGaP_1Mdata_20130522.bim data/VLF.bim
ln -s /media/Lacie_1/DATA/dbGaP/phs000616.v2.p2.VIVA_LA_FAMILIA/GenotypeFiles/phg000386.v1.NIDDK_VIVA_LA_FAMILIA-HumanOmni1_Quad_v1-0_B.genotype-calls-matrixfmt.c1.DS-OBS-IRB-RD.update/7_dbGaP_1Mdata_20130522.fam data/VLF.fam

ln -s /media/Lacie_1/DATA/dbGaP/phs000616.v2.p2.VIVA_LA_FAMILIA/PhenotypeFiles/phs000616.v2.pht003368.v2.p2.VIVA_LA_FAMILIA_Subject.MULTI.txt data/VLF.obesity.txt
ln -s /media/Lacie_1/DATA/dbGaP/phs000616.v2.p2.VIVA_LA_FAMILIA/PhenotypeFiles/phs000616.v2.pht003371.v2.p2.c1.VIVA_LA_FAMILIA_Subject_Phenotypes.DS-OBS-IRB-RD.txt data/VLF.pheno.txt

## Add liftover binary to bin
cp /home/SHARED/SOFTWARE/liftOverPlink/liftOverPlink.py workflows/bin/liftOverPlink.py

## Download GWAS results
wget http://egg-consortium.org/Childhood_Obesity_2019/CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt.gz -P data/
wget https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz -P data/
#/home/SHARED/DISK/DATA/REFERENCES/1000GP_Imputation -> Path to reference for script
