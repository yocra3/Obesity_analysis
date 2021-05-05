library(VariantAnnotation)
library(snpStats)
library(scoreInvHap)
library(tidyverse)
library(invClust)

range <- inversionGR["inv16_009"]
seqlevelsStyle(range) <- "NCBI"
vcf <- readVcf("results/postImputation/merged/2021-04-06/ObesityGWAS.imputatedArrays.mergedDataset.filtered.vcf.gz",
               genome = "hg19", ScanVcfParam(which = range))

sum <- snpSummary(vcf)
vcf.filt <- vcf[sum$a0Freq > 0.05 & sum$a0Freq < 0.95, ]
snpsVCF <- genotypeToSnpMatrix(vcf.filt)
gen <- as(snpsVCF$genotypes, "numeric")
mds <- cmdscale(dist(gen), k = 3)
colnames(mds) <- paste0("MDS", 1:3)

load("results/associations/inversions/2021-04-20/ObesityscoreInvHapClassDF.Rdata")
load("results/associations/inversions/2021-04-20/Obesity.scoreInvHaplist.Rdata")

mds %>%
  data.frame() %>%
  mutate(id = rownames(.)) %>%
  left_join(mutate(scClassDF, id = rownames(scClassDF)), by = "id") %>%
  ggplot(aes(x = MDS1, y = MDS3, col = inv16_009)) +
  geom_point()


hc <- hclust(dist(mds))
haplos <- cutree(hc, k = 11)

plot(mds, col = factor(scClassDF$inv16_009))
text(mds, labels = haplos)
table(haplos, scClassDF$inv16_009)

genos <- c("NaNa", "II", "NaI", "NaI", "NbNb", "NaNb", "NbI", "II", "NaI", "NbI", "NaNa")
genos2 <- genos[haplos]
table(genos2, scClassDF$inv16_009)

genos2.comb <- ifelse(genos2 == "II", "II", ifelse(grepl("I", genos2), "NI", "NN"))
prop.table(table(scClassDF$inv16_009, grepl("SB", rownames(mds))), margin = 2)
prop.table(table(genos2, grepl("SB", rownames(mds))), margin = 2)
prop.table(table(genos2.comb, grepl("SB", rownames(mds))), margin = 2)

chisq.test(table(scClassDF$inv16_009, grepl("SB", rownames(mds))))
chisq.test(table(genos2, grepl("SB", rownames(mds))))
chisq.test(table(genos2.comb, grepl("SB", rownames(mds))))

genos2.combA <- ifelse(genos2 == "NaNa", "RR", ifelse(grepl("Na", genos2), "RA", "AA"))
genos2.combB <- ifelse(genos2 == "NbNb", "RR", ifelse(grepl("Nb", genos2), "RA", "AA"))

prop.table(table(genos2.combA, grepl("SB", rownames(mds))), margin = 2)
chisq.test(table(genos2.combA, grepl("SB", rownames(mds))))

prop.table(table(genos2.combB, grepl("SB", rownames(mds))), margin = 2)
chisq.test(table(genos2.combB, grepl("SB", rownames(mds))))


genos2.comb2 <- ifelse(genos2 == "II", "II", ifelse(genos2 == "NaI", "NI", "NN"))
genos2.comb3 <- ifelse(genos2 == "II", "II", ifelse(genos2 == "NbI", "NI", "NN"))

prop.table(table(genos2.comb2, grepl("SB", rownames(mds))), margin = 2)
chisq.test(table(genos2.comb2, grepl("SB", rownames(mds))))

prop.table(table(genos2.comb3, grepl("SB", rownames(mds))), margin = 2)
chisq.test(table(genos2.comb3, grepl("SB", rownames(mds))))



par(mfrow = c(1, 2))
plot(mds, col = factor(scClassDF$inv16_009))
plot(mds, col = factor(genos2))
