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

All summary statistics of GWAS described in the paper (WorstEpisode, dsWorstEpisode, PHQ9, dsPHQ9, PHQ9Skip, LifetimeMDD) are available [here](link_to_figshare)

## Heritability



## Genetic Correlation



## PRS
