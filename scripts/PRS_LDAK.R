
# Download annotation 
## Run in /home/SHARED/DATA/REFERENCES/GRCh37
wget https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip
unzip genetic_map_b37.zip

wget https://genetics.ghpc.au.dk/doug/bld.zip
unzip bld.zip

## Create 1000 GEnomes fror LDAK (run in /home/SHARED/DATA/REFERENCES/1000Genomes/)
for j in {1..22}; do
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$j.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
-P VCF/
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$j.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi \
-P VCF/
done

for j in {1..22}; do
plink --vcf VCF_version_b/ALL.chr$j.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--make-bed --out plink/chr$j --maf 0.0001 --keep eur.keep
done

rm list.txt
for j in {1..22}; do echo plink/chr$j >> list.txt; done
ldak5.1.linux --make-bed LDAK/ref --mbfile list.txt --exclude-dups YES  --exclude-same YES
ldak5.1.linux --make-bed LDAK/ref_clean --bfile LDAK/ref --exclude-dups YES  --exclude-same YES

awk '{print $2, $1":"$4}' LDAK/ref_clean.bim > LDAK/snp.map
plink --bfile LDAK/ref_clean --update-name LDAK/snp.map --make-bed --out LDAK/ref_bp

for j in {1..22} 
do
  plink --bfile LDAK/ref_bp --chr $j --cm-map ../GRCh37/genetic_map_b37/genetic_map_chr@_combined_b37.txt --make-bed --out LDAK/map$j
done

cat LDAK/map{1..22}.bim | awk '{print $2, $3}' > LDAK/map.all
awk '(NR==FNR){arr[$1]=$2;next}{print $1, $2, arr[$2], $4, $5, $6}' LDAK/map.all LDAK/ref_bp.bim > LDAK/1000g.bim
cp LDAK/ref_bp.bed LDAK/1000g.bed
cp LDAK/ref_bp.fam LDAK/1000g.fam

## Create PRS
## Run in home project
awk < /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g.fam '(NR%3==1){print $0 > "results/LDAK/keepa"}(NR%3==2){print $0 > "results/LDAK/keepb"}(NR%3==0){print $0 > "results/LDAK/keepc"}'

#Identify the non-ambiguous eMERGE SNPs
# awk < eMERGE.bim '(($5=="A"&&$6=="C") || ($5=="A"&&$6=="G") || ($5=="C"&&$6=="T") || \
# ($5=="C"&&$6=="A") ||($5=="G"&&$6=="A") || ($5=="T"&&$6=="C")){print $2}' > eMERGE.nonamb

#Format the summary statistics
## SNP name - A1 - A2 - Direction - P - n
## SNP name as chr:bp
gunzip -c data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | sed 's/:.:.//g' | awk '(NR>1){snp=$1":"$2;a1=$4;a2=$5;\
dir=$7;p=$9;n=$10}(NR==1){print "Predictor A1 A2 Direction P n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") \
&& (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, dir, p, n}' - > results/LDAK/bmi.txt


#Check for duplicates (there are none)
awk '{print $1}' results/LDAK/bmi.txt  | sort | uniq -d > results/LDAK/bmi_dups.txt
grep -v NA results/LDAK/bmi.txt | grep -v -f results/LDAK/bmi_dups.txt > results/LDAK/bmi_clean.txt

awk '(NR==FNR){arr[$2]=$2;ars[$2]=$5$6;next}(FNR==1){print $0}($1 in arr && \
($2$3==ars[$1]||$3$2==ars[$1])){$1=arr[$1];print $0}' /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g.bim results/LDAK/bmi_clean.txt > results/LDAK/bmi_clean_good.txt


awk '(NR==FNR && (($2=="A"&&$3=="C") || ($2=="A"&&$3=="G") || ($2=="C"&&$3=="A") || ($2=="C"&&$3=="T") || ($2=="G"&&$3=="A") || ($2=="G"&&$3=="T") || ($2=="T"&&$3=="C") || ($2=="T"&&$3=="G")) ){arr[$1];next}($2 in arr)\
{print $2}' results/LDAK/bmi_clean_good.txt /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g.bim >  results/LDAK/use.snps
ldak5.1.linux --cut-genes results/LDAK/highld --bfile /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g --genefile /home/SHARED/DATA/REFERENCES/GRCh37/LDAK_annot/highld.txt

## Calculate tagging
ldak5.1.linux --calc-tagging results/LDAK/bld.ldak --bfile /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g \
--extract results/LDAK/use.snps --ignore-weights YES --power -.25 --annotation-number 64 \
--annotation-prefix /home/SHARED/DATA/REFERENCES/GRCh37/LDAK_annot/bld --window-cm 1 --save-matrix YES --max-threads 15

## Estimate heritability contributed by each SNP
ldak5.1.linux --sum-hers results/LDAK/bld.ldak --tagfile results/LDAK/bld.ldak.tagging \
--extract results/LDAK/use.snps --summary results/LDAK/bmi_clean_good.txt --matrix results/LDAK/bld.ldak.matrix --check-sums NO

## Pseudo summary statistics
ldak5.1.linux --pseudo-summaries results/LDAK/bmi --bfile /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g \
--extract results/LDAK/use.snps --summary results/LDAK/bmi_clean_good.txt  --training-proportion .9 --keep results/LDAK/keepa

## Predictor predictor correlations
ldak5.1.linux --calc-cors results/LDAK/cors --bfile /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g \
--extract results/LDAK/use.snps --window-cm 3 --keep results/LDAK/keepb

## Estimate effect size
ldak5.1.linux --mega-prs results/LDAK/bayesR --model bayesr --bfile /home/SHARED/DATA/REFERENCES/1000Genomes/LDAK/1000g \
--extract results/LDAK/use.snps --cors results/LDAK/cors \
--ind-hers results/LDAK/bld.ldak.ind.hers --summary results/LDAK/bmi_clean.txt --summary2 results/LDAK/bmi.train.summaries --window-cm 1

## Compute PRS weights
./ldak5.1.linux --calc-scores bayesR --bfile 1000g --keep keepc --scorefile bayesR.effects.train \
--summary bmi.train.summaries --power 0 --final-effects bayesR.effects.final --exclude highld/genes.predictors.used

## Get PRS
sed 's/:[A-Z]:[A-Z]//g' results/postImputation/merged/2021-06-04/ObesityGWAS.imputatedArrays.mergedDataset.filtered.bim > results/LDAK/eoo.bim
cp results/postImputation/merged/2021-06-04/ObesityGWAS.imputatedArrays.mergedDataset.filtered.bed  results/LDAK/eoo.bed
cp results/postImputation/merged/2021-06-04/ObesityGWAS.imputatedArrays.mergedDataset.filtered.fam  results/LDAK/eoo.fam

plink2 --bfile results/LDAK/eoo --rm-dup  --make-bed --out results/LDAK/eoo.clean
plink --bfile results/LDAK/eoo -exclude results/LDAK/eoo.clean.rmdup.mismatch  --make-bed --out results/LDAK/eoo.clean
plink --bfile results/LDAK/eoo.clean --score results/LDAK/bayesR.effects.best 1 2 5 sum --out results/LDAK/eoo.profile