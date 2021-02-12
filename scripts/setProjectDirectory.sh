#'#################################################################################
#'#################################################################################
# Set project folders' structure and link files
#'#################################################################################
#'#################################################################################

## Make folders
mkdir data
mkdir results

## Link raw data
### Genotypes
ln -s /media/Lacie_1/DATA/Epicuro/FullDataTable.txt data/Epicuro.FullDataTable.txt
ln -s /home/SHARED/DISK/CNIO/1a_cohort/FullDataTable_pennCNV.txt data/Obesity1.FullDataTable.txt
ln -s /home/SHARED/DISK/CNIO/2a_cohort/FullDataTable_pennCNV.txt data/Obesity2.FullDataTable.txt
ln -s /home/SHARED/DISK/CNIO/JArgente_EOO/pennCNV/FullDataTable_pennCNV.txt data/Obesity3.FullDataTable.txt

## Phenotypes
ln -s /home/SHARED/DISK/CNIO/Taula\ resum\ pacients\ \(def\).xlsx data/Obesity.SamplesSummary.xlsx

#/home/SHARED/DISK/DATA/REFERENCES/1000GP_Imputation -> Path to reference for script