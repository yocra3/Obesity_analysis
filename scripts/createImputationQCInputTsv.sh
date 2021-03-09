#'#################################################################################
#'#################################################################################
# Create input tsv for imputation
#'#################################################################################
#'#################################################################################

cohorts=( Obesity1 Obesity2 Obesity3 Epicuro )

for cohort in "${cohorts[@]}"
do
  echo -e "$cohort\tresults/postImputation/MichiganOutput/$cohort/" >> data/inputImputationQC.tsv
done
