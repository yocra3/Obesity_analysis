#'#################################################################################
#'#################################################################################
# Send VCFs for imputation
#'#################################################################################
#'#################################################################################

## Install imputation bot
curl -sL imputationbot.now.sh | bash
mv imputation* bin/

## Configure imputation bot
bin/imputationbot add-instance

# Parameters
## Imputationserver Url [https://imputationserver.sph.umich.edu]:

# Send jobs
## Parameters
##  - Panel: HRC
##  - r2Filter: 0.3
##  - Population: EUR
##  - Build: hg19
vcf=results/preImputation/VCF/
bin/imputationbot impute --files $vcf/Epicuro/2021-02-17/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Epicuro --project Obesity-proj
bin/imputationbot impute --files $vcf/Obesity1/2021-02-17/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity1 --project Obesity-proj
bin/imputationbot impute --files $vcf/Obesity2/2021-02-17/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity2 --project Obesity-proj
bin/imputationbot impute --files $vcf/Obesity3/2021-02-17/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity3 --project Obesity-proj


## Download files
bin/imputationbot download Obesity-proj --output results/postImputation/MichiganOutput/
