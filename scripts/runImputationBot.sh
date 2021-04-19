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
bin/imputationbot impute --files $vcf/Epicuro/2021-03-16/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Epicuro --project Obesity-proj
bin/imputationbot impute --files $vcf/Obesity1/2021-03-16/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity1 --project Obesity-proj
bin/imputationbot impute --files $vcf/Obesity2/2021-03-16/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity2 --project Obesity-proj
bin/imputationbot impute --files $vcf/Obesity3/2021-03-16/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity3 --project Obesity-proj

## Download files
bin/imputationbot download Obesity-proj --output results/postImputation/MichiganOutput/

## Decompress
folds=$(ls results/postImputation/MichiganOutput)
for fold in $folds
do
  new="$(cut -d'-' -f7 <<<$fold)"
  for file in $(ls results/postImputation/MichiganOutput/$fold/local)
  do
    unzip -P `cat results/postImputation/MichiganOutput/$fold/pass` results/postImputation/MichiganOutput/$fold/local/$file -d results/postImputation/MichiganOutput/$new
  done
done
