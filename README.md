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

All summary statistics of GWAS described in the paper (WorstEpisode, dsWorstEpisode, PHQ9, dsPHQ9, PHQ9Skip, LifetimeMDD) are available [here](https://doi.org/10.6084/m9.figshare.22041212)

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

## PRS Pleiotropy 

R scripts for gathering and summarising PRS prediction R2s across folds are in the ```PRS_pleio``` directory:

### PRS Pleiotropy analysis 

```getdata.R```: script used to gather R2s from all fold predictions 
```utilities.R```: misc functions  

