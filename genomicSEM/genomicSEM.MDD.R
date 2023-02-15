## try genomicSEM
require(GenomicSEM)
require(Matrix)
require(stats)
library(psych)
library(GPArotation)

## lifetime symptoms + lifetimeMDD 
setwd("/no-backup/caina/biobank/genomicSEM/mungestats")
traits=c(paste0("Munge_LifetimeSymptoms_A",c(1:4,6:9),".sumstats.gz"),"munge_ukbmdd.sumstats.gz")
trait.names=c(paste0("A",c(1:4,6:9)),"LifetimeMDD")
prev=read.table("Prevalence_Phenos_FullSample.txt",header=T)
sample.prev=c(prev$Prevalence[match(trait.names,prev$Pheno)],0.243)
population.prev=sample.prev
ld="/projects/biobank/release2/ldscores/maf5/"
wld="/projects/biobank/release2/ldscores/maf5/"

LDSCoutput=ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS=commonfactor(covstruc = LDSCoutput, estimation="DWLS")
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))

save(LDSCoutput, file="WE.LifetimeMDD.genomicSEM.RData")

## PHQ9 symptoms + lifetimeMDD 
setwd("/no-backup/caina/biobank/genomicSEM/mungestats")
traits=c(paste0("Munge_LenCurrent_A",c(1:9),".sumstats.gz"),"munge_ukbmdd.sumstats.gz")
trait.names=c(paste0("A",c(1:9)),"LifetimeMDD")
prev=read.table("PrevalenceTable_LenCurrent.txt",header=T)
prev$Pheno=gsub("lencurrent","",prev$Pheno)
sample.prev=c(prev$Prev[match(paste0("A",c(1:9)),prev$Pheno)],0.243)
population.prev=sample.prev
ld="/projects/biobank/release2/ldscores/maf5/"
wld="/projects/biobank/release2/ldscores/maf5/"

LDSCoutput=ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS=commonfactor(covstruc = LDSCoutput, estimation="DWLS")
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))

save(LDSCoutput, file="PHQ9.LifetimeMDD.genomicSEM.RData")

## PHQ9Skip + PHQA1A2 symptoms + lifetimeMDD 
setwd("/no-backup/caina/biobank/genomicSEM/mungestats")
traits=c(paste0("Munge_LenCurrent_A",c(1:2),".sumstats.gz"),paste0("Munge_LenStop_A",c(3:9),".sumstats.gz"),"munge_ukbmdd.sumstats.gz")
trait.names=c(paste0("A",c(1:2)),paste0("SkipA",c(3:9)),"LifetimeMDD")
prev=read.table("prevall.txt",header=T)
prev1=prev[which(prev$data=="phq"),]
prev2=prev[which(prev$data=="phqskip"),]
sample.prev=c(prev1$samp_prev[match(paste0("A",c(1:2)),prev1$symp)],prev2$samp_prev[match(paste0("A",c(3:9)),prev2$symp)],0.243)
population.prev=sample.prev
ld="/projects/biobank/release2/ldscores/maf5/"
wld="/projects/biobank/release2/ldscores/maf5/"

LDSCoutput=ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS=commonfactor(covstruc = LDSCoutput, estimation="DWLS")
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))

save(LDSCoutput, file="PHQ9SkipA1A2.LifetimeMDD.genomicSEM.RData")
