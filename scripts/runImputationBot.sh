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
bin/imputationbot impute --files $vcf/Epicuro/2021-05-31/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Epicuro --project Obesity
bin/imputationbot impute --files $vcf/Obesity1/2021-05-31/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity1 --project Obesity
bin/imputationbot impute --files $vcf/Obesity2/2021-05-31/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity2 --project Obesity
bin/imputationbot impute --files $vcf/Obesity3/2021-05-31/ --refpanel hrc-r1.1 --r2Filter 0.3 --population eur --name Obesity3 --project Obesity
bin/imputationbot impute --files $vcf/VLF/2021-05-11 --refpanel 1000g-phase-3-v5 --r2Filter 0.3 --population amr --name VLF --project Obesity

## Download files
bin/imputationbot download Obesity --output results/postImputation/MichiganOutput/

## Decompress
folds=$(ls results/postImputation/MichiganOutput)
for fold in $folds
do
  new="$(cut -d'-' -f6 <<<$fold)"
  for file in $(ls results/postImputation/MichiganOutput/$fold/local)
  do
    unzip -P `cat results/postImputation/MichiganOutput/$fold/pass` results/postImputation/MichiganOutput/$fold/local/$file -d results/postImputation/MichiganOutput/$new
  done
  echo $new
done


## Viva la Familia. Use 1000 Genomes
bin/imputationbot impute --files $vcf/VLF/2021-05-11 --refpanel 1000g-phase-3-v5 --r2Filter 0.3 --population amr --name VLF
bin/imputationbot download job-20210430-063418-809 --password 3CebxjK2zm6EAW --output results/postImputation/MichiganOutput/

mkdir results/postImputation/MichiganOutput/VLF
mv results/postImputation/MichiganOutput/job-20210430-063418-809-VLF/local/*.gz results/postImputation/MichiganOutput/VLF/
