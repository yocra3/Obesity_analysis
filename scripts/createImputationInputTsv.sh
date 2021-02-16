#'#################################################################################
#'#################################################################################
# Create input tsv for imputation
#'#################################################################################
#'#################################################################################

cohorts=( Obesity1 Obesity2 Obesity3 Epicuro )

for cohort in "${cohorts[@]}"
do
  echo -e "$cohort\tdata/${cohort}.FullDataTable.txt\tdata/${cohort}.annot.csv" >> data/inputImputation.tsv
done
