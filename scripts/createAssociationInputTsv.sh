#'#################################################################################
#'#################################################################################
# Create input tsv for association step
#'#################################################################################
#'#################################################################################

cohorts=( Obesity1 Obesity2 Obesity3 Epicuro )

for cohort in "${cohorts[@]}"
do
  echo -e "$cohort\tresults/preImputation/plink/${cohort}/2021-03-16/${cohort}.fam" >> data/inputAssociation.tsv
done
