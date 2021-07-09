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
bm = BlockMatrix.read('s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm')

## Export index table
idx_pd = ht_idx.to_pandas()
idx_pd[['a0','a1']] = pd.DataFrame(idx_pd.alleles.tolist(), index = idx_pd.index)
idx_pd = idx_pd.rename(columns={'locus.contig':'chr', 'locus.position':'pos'})



stats = pd.read_csv('data/UKB_ldref/full_variant_qc_metrics.txt', sep = '\t')
stats = stats[['chrom', 'pos', 'ref', 'alt', 'rsid', 'varid',
       'pass_gnomad_genomes', 'high_quality', 'af_EUR', 'gnomad_genomes_af_EUR'] ]
stats = stats.rename(columns={'chrom':'chr'})
stats['chr'] = stats['chr'].astype(str)

idx_pd['chr'] = idx_pd['chr'].astype(str)

idx_comb = pd.merge(idx_pd, stats, on = ['chr', 'pos'])
idx_comb = idx_comb[idx_comb.high_quality]
idx_comb = idx_comb.drop_duplicates(subset=['idx'])

n_rows = 30000
var = []
for chr_id in idx_comb.chr.unique():
  rows = sum(idx_comb.chr == chr_id)
  g = math.ceil(rows/n_rows)
  
  var = np.append(var, np.concatenate([([i]*n_rows) for i in range(0, g)], axis=0)[:rows])

idx_comb['group'] = var
idx_comb['group'] = idx_comb['group'].astype(str)
idx_out = idx_comb[['chr', 'pos', 'a0', 'a1', 'group', 'rsid', 'varid', 'idx'] ]
idx_out.to_csv('data/UKB_ldref/mapLDref.tsv.gz', compression='gzip',  sep = "\t", index = False)

for chr_id in idx_comb.chr.unique():
  
  idx_comb_chr = idx_comb[idx_comb.chr == chr_id]
  for group_id in idx_comb_chr.group.unique():
    chridx = idx_comb_chr[idx_comb_chr.group == group_id].idx.tolist()
    
    ext = 'data/UKB_ldref/tmp/chr' + chr_id + '.' + group_id
    if path.exists(ext + '.csv.bgz'):
      continue
    os.system("python3 scripts/makeLDfiles.py " + chr_id + ' ' + group_id)

bmchr[:30, :30].write('data/UKB_ldref/tmp/cot.bm', force_row_major=True)

BlockMatrix.export(path_in='data/UKB_ldref/tmp/cot.bm',
  path_out='data/UKB_ldref/tmp/cot.csv.bgz',
  delimiter='\t',
  entries='upper')

