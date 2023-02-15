# MDDSymptoms
Code to reproduce analyses in the MDD symptoms paper: link_to_paper

## GWAS 

All GWAS in the paper are performed using [PLINK2](https://www.cog-genomics.org/plink/2.0/)

```
#GWAS on quantitative phenotypes 
plink2 --bfile $bfile --pheno $pfile --linear hide-covar --variance-standardize --covar $covar --out $outfile 
#GWAS on binary phenotypes 
plink2 --bfile $bfile --pheno $pfile --1 --logistic hide-covar --covar-variance-standardize --ci 0.95 --covar $covar --out $outfile 
```

## Summary statistics 

All summary statistics of GWAS described in the paper (WorstEpisode, dsWorstEpisode, PHQ9, dsPHQ9, PHQ9Skip, LifetimeMDD) are available [here](https://doi.org/10.6084/m9.figshare.22041212). All summary statistics have the following columns:

```
CHR: chromosome
SNP: rsID of SNP
POS: base-pair coordinate
A1: effect allele
A2: other allele
MAF: minor allele frequency
OR: odds ratio of effect allele
SE: standard error of log(OR)
P: P value of association
N: sample size (which may differ between different SNPs due to missing genotypes at SNPs)
```

## Heritability
All heritability estimation in the paper are performed using [ldsc](https://github.com/bulik/ldsc)
```
# munge
ldsc/munge_sumstats.py --sumstats $sumstats --out $outfile --snp SNP --a1 A1 --a2 A2 --p P --frq MAF --signed-sumstats OR,1 --N-cas $ncases --N-con $ncontrols
# h2
ldsc/ldsc.py --h2 $mungefile.sumstats.gz --ref-ld-chr $reffile --w-ld-chr $reffile --samp-prev $prev --pop-prev $prev --out $outfile
```

## Genetic Correlation
All genetic correlation estimation in the paper are performed using [ldsc](https://github.com/bulik/ldsc)
```
ldsc/ldsc.py --rg munge1.sumstats.gz,munge2.sumstats.gz --ref-ld-chr $reffile --w-ld-chr $reffile --samp-prev $prev1,$prev2 --pop-prev $prev1,$prev2 --out $outfile
```

## Mendelian Randomization (MR) 

R scripts for univariate (UVMR) and multivariate MR analyses (MR-BMA) are in the ```MR``` directory:

```read_format_datasets.R```: script for reading in and formatting data for MR analyses 
```uvmr.R```: script for UVMR analysis
```mvmr.R```: script for MR-BMA analysis
```run_mvmr_withpvalue.R```: obtaining P values for MR-BMA analysis 

## PRS predictions

We use the following PRSice-2 command line, switching --binary-target T/F for binary and quantitative target phenotypes. We performed PRS in 10-fold cross validation in UKBiobank for in-sample PRS analysis, so phenotype/genotype files are in folds, to do this without having separate files one can use --remove/--keep options. Confidence intervals for prediction R2 are derived from standard errors between estimates from the 10 folds. 

```
a=$GWAS ## summary statistics 
b=$targetpheno ## 50 non-MDD phenotypes as detailed in Supplementary Table 4 of the paper
for test in {1..10}
do 
PRSice_linux --base $a.ma --cov $covar --num-auto 22 \
--a1 A1 --a2 A2 --pvalue P --snp SNP --stat or --OR --base-maf MAF:0.05 \
--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
--interval 5e-05 --upper 0.5 --lower 5e-08 --num-auto 22 \
--bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
--out $outdir/$b/$a.test$test.prs \
--pheno $phendir/test$test.$b.phen \
--target $genodir/tenfold/allchr.test$test --binary-target F \
--thread 1
done 
```

## PRS Pleiotropy 

R scripts for gathering and summarising PRS prediction R2s across folds are in the ```PRS_pleio``` directory:

## genomicSEM 

R scipts for [genomicSEM](https://github.com/GenomicSEM/GenomicSEM) analyses are in the ```genomicSEM``` directory:

```genomicSEM.MDD.R```: script used to gather R2s from all fold predictions 
```EFACFA.MDD.R```: script used to gather R2s from all fold predictions 
```plotCFA.MDD.R```: plotting scripts   
