library("plyr") 
library("tidyverse")
library("stringr") 
library("vroom") 
library("data.table") 
library("MendelianRandomization") 
library("TwoSampleMR")
library("MRPRESSO") 
library("ieugwasr")
library("ggplot2")
library("reshape2")
source("/read_format_datasets.R")



exposures=c("cpds_lencurrentA1","cpds_lencurrentA2","cpds_lencurrentA3","cpds_lencurrentA4","cpds_lencurrentA6","cpds_lencurrentA7","cpds_lencurrentA8","cpds_lencurrentA9")
#exposures=c("cpds_lifetimeA1","cpds_lifetimeA2","cpds_lifetimeA3","cpds_lifetimeA4","cpds_lifetimeA6","cpds_lifetimeA7","cpds_lifetimeA8","cpds_lifetimeA9")
#outcomes="pgc29"
#outcomes="ipsych2012"
outcomes="ipsych2015i"

#exposures="pgc29"
#exposures="ipsych2012"
#exposures="ipsych2015i"
#outcomes=c("cpds_lencurrentA1","cpds_lencurrentA2","cpds_lencurrentA3","cpds_lencurrentA4","cpds_lencurrentA6","cpds_lencurrentA7","cpds_lencurrentA8","cpds_lencurrentA9")
#outcomes=c("cpds_lifetimeA1","cpds_lifetimeA2","cpds_lifetimeA3","cpds_lifetimeA4","cpds_lifetimeA6","cpds_lifetimeA7","cpds_lifetimeA8","cpds_lifetimeA9")





for (exposure in exposures){
	for (outcome in outcomes){
		# read data
		data_exposure=exposure_format(read_data(exposure))
		data_outcome=outcome_format(read_data(outcome))
		# Histogram of p-values showing relationship between SNPs and the exposure
		pdf(paste0(path_output,"/plot/hist_",exposure,"_to_",outcome,".pdf"))
		pl_hist=hist(data_exposure$p_exp,main=paste0(exposure,"_pvalue"))
		print(pl_hist)
		dev.off()
		# Subset sig. SNPs.
		sigsnp = data_exposure[data_exposure$p_exp<5*10^(-6),]
		sigsnp = na.omit(sigsnp)
		sigsnp
		# merge sigsnp_exposure and outcome
		merge_tmp = dplyr::inner_join(x=sigsnp, y=data_outcome, by="rsid")  
		# remove duplicates
		merge_tmp = merge_tmp[!(duplicated(merge_tmp$rsid)), ]
		# swap/delete unmatched alleles
		for (i in 1:nrow(merge_tmp)){
		  if (merge_tmp[i,"a1_exp"]==merge_tmp[i,"a1_out"] & merge_tmp[i,"a2_exp"]==merge_tmp[i,"a2_out"]) {
		    next  
		  }
		  else if (merge_tmp[i,"a1_exp"]==merge_tmp[i,"a2_out"] & merge_tmp[i,"a2_exp"]==merge_tmp[i,"a1_out"]) {
		    merge_tmp[i,"beta_out"]=-1*merge_tmp[i,"beta_out"]
		    merge_tmp[i,"maf_out"]=1-merge_tmp[i,"maf_out"]
		  }
		  else if (merge_tmp[i,"a1_exp"]!=merge_tmp[i,"a1_out"] & merge_tmp[i,"a1_exp"]!=merge_tmp[i,"a2_out"]) {
		    merge_tmp=merge_tmp[-i,]
		  }
		  else if (merge_tmp[i,"a2_exp"]!=merge_tmp[i,"a1_out"] & merge_tmp[i,"a2_exp"]!=merge_tmp[i,"a2_out"]) {
		    merge_tmp=merge_tmp[-i,]
		  }
		}
		# clumpe SNPs
		colnames(merge_tmp)[colnames(merge_tmp)=="p_exp"]="pval"
		tab_clump = ieugwasr::ld_clump(merge_tmp) 
		write.csv(nrow(tab_clump),paste0(path_output,"/snp/nsnps_clumping_",exposure,"_to_",outcome,".csv"))  
		# F_Stats
		tab_clump$F = (tab_clump$beta_exp/tab_clump$se_exp)^2
		write.csv(ave(tab_clump$F)[1],paste0(path_output,"/f_stats/fstats_",exposure,"_to_",outcome,".csv"))    
		tab_clump
		# MR
		print("MR ...")
		mr_input = mr_input(
		  bx = tab_clump$beta_exp, 
		  bxse = tab_clump$se_exp, 
		  by = tab_clump$beta_out, 
		  byse = tab_clump$se_out, 
		  exposure = exposure, 
		  outcome = outcome, 
		  snps = tab_clump$rsid
		  )		

		# try MR to avoid stop looping
		controller_mr = try(mr_allmethods(mr_input),silent = TRUE) 
		  
		# run MR
		if ( !inherits(controller_mr, "try-error") ) {
		  write.csv(mr_allmethods(mr_input)@Values,paste0(path_output,"/uvmr/uvmr_allmethods_",exposure,"_to_",outcome,".csv"))
		  
		  # IVW plot
		  pdf(paste0(path_output,"/plot/ivw_",exposure,"_to_",outcome,".pdf"))
		  p2=mr_plot(mr_input, interactive = F)
		  print(p2)
		  dev.off() 
		  mr_plot(mr_input, interactive = F)

		  # main plot
		  pdf(paste0(path_output,"/plot/main_",exposure,"_to_",outcome,".pdf"))
		  p3=mr_plot(mr_allmethods(mr_input, method = "main"))
		  print(p3)
		  dev.off()
		  mr_plot(mr_allmethods(mr_input, method = "main"))

		  # MR egger plot
		  pdf(paste0(path_output,"/plot/egger_",exposure,"_to_",outcome,".pdf"))
		  p4=mr_plot(mr_allmethods(mr_input, method = "egger"))
		  print(p4)
		  dev.off()
		  mr_plot(mr_allmethods(mr_input, method = "egger"))
		  
		  tmp=mr_allmethods(mr_input)@Values
		}
	}
}
