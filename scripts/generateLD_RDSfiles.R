#!/usr/local/bin/Rscript
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
bgfile <- args[1]

## Rename as gz
bgfile2 <- gsub("bgz", "gz", bgfile)
system(paste("ln -s", bgfile, bgfile2))

lines <- system(paste("bgzip -cd", bgfile2, "| wc -l"), intern = TRUE)
lines <-  as.numeric(strsplit(lines, " ")[[1]][1])
mat <- matrix(scan(gzfile(bgfile2), sep = "\t"), nrow = lines, byrow = TRUE)
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
mat_sp <- Matrix(mat, sparse = TRUE)

rdsfile <- gsub("csv.bgz", "rds", bgfile)
saveRDS(mat_sp, file = rdsfile)