#'#################################################################################
#'#################################################################################
#' Generate LD-matrices from UKbiobank for LDpred
#'#################################################################################
#'#################################################################################

import hail as hl
from hail.linalg import BlockMatrix
ukbb_pan_ancestry import *
import h5py
import pandas as pd
import numpy as np
import math
from os import path
import os

## Load data
ht_idx = hl.read_table('s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht')                                              
ht_score = hl.read_table('s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldscore.ht')                                              

## Import stats data
stats = pd.read_csv('data/UKB_ldref/full_variant_qc_metrics.txt', sep = '\t')
stats = stats[['varid', 'high_quality', 'af_EUR', 'gnomad_genomes_af_EUR'] ]

## Import scores  table
score_pd = ht_score.to_pandas()
score_pd[['a0','a1']] = pd.DataFrame(score_pd.alleles.tolist(), index = score_pd.index)
score_pd = score_pd.rename(columns={'locus.contig':'chr', 'locus.position':'pos'})
score_pd = score_pd[['chr', 'pos', 'a0', 'a1', 'rsid', 'varid', 'AF', 'ld_score'] ]

## Import index table
idx_pd = ht_idx.to_pandas()
idx_pd[['a0','a1']] = pd.DataFrame(idx_pd.alleles.tolist(), index = idx_pd.index)
idx_pd = idx_pd.rename(columns={'locus.contig':'chr', 'locus.position':'pos'})
idx_comb = pd.merge(idx_pd, score_pd, on = ['chr', 'pos', 'a0', 'a1'])
del idx_pd
del score_pd

idx_comb['chr'] = idx_comb['chr'].astype(str)

idx_out = pd.merge(idx_comb, stats, on = ['varid'])
idx_out = idx_out[idx_out.high_quality]
idx_out = idx_out.drop_duplicates(subset=['idx'])

n_rows = 10000
var = []
for chr_id in idx_out.chr.unique():
  rows = sum(idx_out.chr == chr_id)
  g = math.ceil(rows/n_rows)
  
  var = np.append(var, np.concatenate([([i]*n_rows) for i in range(0, g)], axis=0)[:rows])

idx_out['group'] = var
idx_out['group'] = idx_out['group'].astype(str)
idx_out = idx_out[['chr', 'pos', 'a0', 'a1', 'group', 'rsid', 'varid', 'idx', 'AF', 'ld_score'] ]
idx_out.to_csv('data/UKB_ldref/mapLDref.tsv.gz', compression='gzip',  sep = "\t", index = False)
