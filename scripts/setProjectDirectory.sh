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
ln -s /media/Lacie_1/DATA/Epicuro/FullDataTable.txt data/Epicuro.FullDataTable.txt
ln -s /home/SHARED/DISK/CNIO/1a_cohort/FullDataTable_pennCNV.txt data/Obesity1.FullDataTable.txt
ln -s /home/SHARED/DISK/CNIO/2a_cohort/FullDataTable_pennCNV.txt data/Obesity2.FullDataTable.txt
ln -s /home/SHARED/DISK/CNIO/JArgente_EOO/pennCNV/FullDataTable_pennCNV.txt data/Obesity3.FullDataTable.txt

### Array annotation
ln -s /media/Lacie_1/DATA/CNIO_CEGEN/Arrays/HumanOmni1-Quad\ \(1a\ Tanda\)/HumanOmni1-Quad_v1-0_H.csv data/Obesity1.annot.csv
ln -s /media/Lacie_1/DATA/CNIO_CEGEN/Arrays/OmniExpress-12v1\ \(2a\ Tanda\)/HumanOmniExpress-12-v1-0-K.csv data/Obesity2.annot.csv
ln -s /media/Lacie_1/DATA/CNIO_CEGEN/Arrays/OmniExpress-24v1.2\ \(3a\ Tanda\)/InfiniumOmniExpress-24v1-2_A1.csv data/Obesity3.annot.csv
ln -s /home/SHARED/DISK/PROJECTS/GWAS_OBE/DATA/EPICURO/Human1Mv1_C.csv data/Epicuro.annot.csv

## Phenotypes
ln -s /home/SHARED/DISK/CNIO/Taula\ resum\ pacients\ \(def\).xlsx data/Obesity.SamplesSummary.xlsx

#/home/SHARED/DISK/DATA/REFERENCES/1000GP_Imputation -> Path to reference for script
