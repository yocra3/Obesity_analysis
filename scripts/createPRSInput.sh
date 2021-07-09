#'#################################################################################
#'#################################################################################
# Create input files for PGS computation
#'#################################################################################
#'#################################################################################

## Locke AE (GWAS which PGS0000027 is based)
bgzip -cd data/25673413-GCST002783-EFO_0004340.h.tsv.gz | \
awk  '{print $3, $4, $5, $6, $17, $18}' | bgzip -c > data/PGS_inputfile/Locke.input38.txt.gz

## Code in R to perform build transformation
library(bigsnpr)
library(bigreadr)
options(scipen=999)
sumstats <- fread2("data/PGS_inputfile/Locke.input38.txt.gz", na.strings = "NULL",
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se"))
sumstats$pos <- as.numeric(sumstats$pos )
sumstats_filt <- subset(sumstats, !is.na(pos) & chr != "NA")
sumstats_19 <- snp_modifyBuild(sumstats_filt, "/home/SHARED/SOFTWARE/liftOver", "hg38", "hg19")
sumstats_19filt <- subset(sumstats_19, !is.na(pos) & chr != "NA")
write.table(sumstats_19filt, file = "data/PGS_inputfile/Locke.input19.txt", quote = FALSE, row.names = FALSE)

bgzip data/PGS_inputfile/Locke.input19.txt

## GIANT
bgzip -cd data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | \
awk  '{print $1, $2, $5, $4, $7, $8}' | bgzip -c > data/PGS_inputfile/GIANT.input.txt.gz

## EGG European
bgzip -cd data/CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt.gz | \
awk  '{print $2, $3, $5, $4, $18, $19}' | bgzip -c > data/PGS_inputfile/EGG.input.txt.gz

## SCOOP obesity
bgzip -cd data/SCOOP_STILTS_ldcorrected.gz | \
awk  '{print $2, $3, $5, $4, $6, $8}' | bgzip -c > data/PGS_inputfile/SCOOP_obesity.input.txt.gz

bgzip -cd data/SCOOP_UKHLS_ldcorrected.gz | \
awk  '{print $2, $3, $5, $4, $6, $8}' | bgzip -c > data/PGS_inputfile/SCOOP_extreme.input.txt.gz
