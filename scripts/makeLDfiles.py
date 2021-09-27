import hail as hl
from hail.linalg import BlockMatrix
from os import path
import sys
import pandas as pd

chr_id, group_id =  sys.argv[1], sys.argv[2]
print(chr_id + ' '  + group_id)
idx_comb = pd.read_csv("mapLDref.tsv.gz", sep = "\t") 
idx_comb['chr'] = idx_comb['chr'].astype(str)
idx_comb['group'] = idx_comb['group'].astype(str)

idx_comb_chr = idx_comb[idx_comb.chr == chr_id]
chridx = idx_comb_chr[idx_comb_chr.group == group_id].idx.tolist()
    
ext = 'chr' + chr_id + '.' + group_id
## Load data
bm = BlockMatrix.read('s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm')


bmchr = bm.filter(chridx, chridx)
bmchr.write(ext + '.bm', force_row_major=True)
BlockMatrix.export(ext + '.bm', ext + '.csv.bgz', delimiter='\t')
