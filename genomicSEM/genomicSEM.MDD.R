## try genomicSEM
require(GenomicSEM)
require(Matrix)
require(stats)
library(psych)
library(GPArotation)

## LD reference
ld="../ldsc/"
wld="../ldsc/"

## prevalence 
prev=read.table("prevall.txt",header=T)

## FOR iPSYCH2012 and iPSYCH2015i MDD, change "LifetimeMDD" to "iPSYCH2012" or "iPSYCH2015i"

## WorstEpisode symptoms + LifetimeMDD 
traits=c(paste0("WorstEpisode_A",c(1:4,6:9),".sumstats.gz"),"LifetimeMDD.sumstats.gz")
trait.names=c(paste0("A",c(1:4,6:9)),"LifetimeMDD")
sample.prev=c(prev$Prevalence[match(trait.names,prev$Pheno)],0.243)
population.prev=sample.prev
LDSCoutput=ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS=commonfactor(covstruc = LDSCoutput, estimation="DWLS")
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))
save(LDSCoutput, file="WE.LifetimeMDD.genomicSEM.RData")

## PHQ9 symptoms + LifetimeMDD 
traits=c(paste0("PHQ9_A",c(1:9),".sumstats.gz"),"LifetimeMDD.sumstats.gz")
trait.names=c(paste0("A",c(1:9)),"LifetimeMDD")
prev$Pheno=gsub("lencurrent","",prev$Pheno)
sample.prev=c(prev$Prev[match(paste0("A",c(1:9)),prev$Pheno)],0.243)
population.prev=sample.prev
LDSCoutput=ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS=commonfactor(covstruc = LDSCoutput, estimation="DWLS")
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))
save(LDSCoutput, file="PHQ9.LifetimeMDD.genomicSEM.RData")

## PHQ9Skip + PHQA1A2 symptoms + lifetimeMDD 
traits=c(paste0("WorstEpisode_A",c(1:2),".sumstats.gz"),paste0("PHQ9Skip_A",c(3:9),".sumstats.gz"),"LifetimeMDD.sumstats.gz")
trait.names=c(paste0("A",c(1:2)),paste0("SkipA",c(3:9)),"LifetimeMDD")
prev1=prev[which(prev$data=="phq"),]
prev2=prev[which(prev$data=="phqskip"),]
sample.prev=c(prev1$samp_prev[match(paste0("A",c(1:2)),prev1$symp)],prev2$samp_prev[match(paste0("A",c(3:9)),prev2$symp)],0.243)
population.prev=sample.prev
LDSCoutput=ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS=commonfactor(covstruc = LDSCoutput, estimation="DWLS")
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))
save(LDSCoutput, file="PHQ9SkipA1A2.LifetimeMDD.genomicSEM.RData")
