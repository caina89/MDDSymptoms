library("plyr") 
library("tidyverse")
library("stringr") 
library("vroom") 
library("data.table") 
library("MendelianRandomization") 
#library("TwoSampleMR")
library("MRPRESSO") 
library("ieugwasr")
library(ggplot2)
library(reshape2)
library("stringr")

source("/path/softwares/mvmr/demo_AMD/summary_mvMR_SSS.R")
source("/path/softwares/mvmr/demo_AMD/summary_mvMR_BF.R")
source("/path/read_format_datasets.R")
source("/path/run_mvmr_withpvalue.R")




items=c("pgc29","cpds_lencurrentA1","cpds_lencurrentA2","cpds_lencurrentA3","cpds_lencurrentA4","cpds_lencurrentA6","cpds_lencurrentA7","cpds_lencurrentA8","cpds_lencurrentA9")
#items=c("ipsych2012","cpds_lencurrentA1","cpds_lencurrentA2","cpds_lencurrentA3","cpds_lencurrentA4","cpds_lencurrentA6","cpds_lencurrentA7","cpds_lencurrentA8","cpds_lencurrentA9")
#items=c("ipsych2015i","cpds_lencurrentA1","cpds_lencurrentA2","cpds_lencurrentA3","cpds_lencurrentA4","cpds_lencurrentA6","cpds_lencurrentA7","cpds_lencurrentA8","cpds_lencurrentA9")


data_all=NULL 

for (i in items){
    if (i == items[1]) {
    	data_all=mvmr_format(read_data(i),i)
    }
    else {
	data_tmp=mvmr_format(read_data(i),i)
	data_all=dplyr::inner_join(x=data_all, y=data_tmp, by="rsid")    	
    }
}

run_mvmr_withpvalue(data_all,items)
