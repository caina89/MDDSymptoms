run_mvmr_withpvalue=function(data_all,items){
    for (i in items){
    print(i)  

    # set one item as outcome
    data_exp = data_all[!colnames(data_all) %in% c(paste("beta_",i,sep=""), paste("se_",i,sep=""), paste("p_",i,sep=""), paste("a1_",i,sep=""), paste("a2_",i,sep=""))]
    data_out = data_all[colnames(data_all) %in% c("rsid",paste("beta_",i,sep=""), paste("se_",i,sep=""), paste("p_",i,sep=""), paste("a1_",i,sep=""), paste("a2_",i,sep=""))]
    colnames(data_out) = str_replace_all(colnames(data_out), i, "out") # change colnames for output 

    # subset iv-related rows, assign the smallest p-value in each row among all exposures
    p_tmp = select(data_all,contains("p_")) 
    p_tmp = p_tmp[!colnames(p_tmp) %in% paste("p_",i,sep="")]
    p_min = apply(p_tmp, FUN=min, MARGIN=1) # find the smallest p-value in each row among all exposures
    data_exp$pval = p_min # name as pval for furture steps
    iv = data_exp[data_exp$pval<5*10^(-6),] # select potential ivs by minimum p-value
    data_merge = dplyr::inner_join(x=data_out, y=iv, by="rsid") # merge iv-exposure and outcome
    data_merge = data_merge[!(duplicated(data_merge$rsid)), ] # remove duplicate variants
    #data_merge$pval=as.numeric(data_merge$pval)

    # clump and get clumped snp list
    hist(data_merge$pval)
    print("summary pval:")
    print(summary(data_merge$pval))
    data_clump = ieugwasr::ld_clump(data_merge)
    
    # define variables
    tmp = select(data_clump,contains("beta_"))
    beta_x = tmp[!colnames(tmp)%in%"beta_out"]
    beta_x = as.matrix(beta_x)
    hist(c(beta_x))
    list_exposure = colnames(beta_x)
    list_snp = data_clump$rsid
    beta_y = data_clump$beta_out
    se_y = data_clump$se_out
    beta_x_ivw = beta_x/se_y # ivw variable
    beta_y_ivw = beta_y/se_y # ivw variable
    cor_beta = cor(beta_x_ivw) # the correlation structure between exposures  

    # Visualising the correlation structure
    melted_cormat = melt(cor_beta)
    plt = ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", 
      midpoint = 0, limit = c(-1,1), space = "Lab", 
      name="Correlation \n between betas"
       ) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
    coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
     )
    print(plt)  

    ## Univariable MR
    ivw = matrix(0,ncol=4, nrow= ncol(beta_x_ivw))
    rownames(ivw) = list_exposure
    colnames(ivw) = c("Estimate","StdError","tvalue","pval")
    for(c in 1:ncol(beta_x_ivw)){
      ivw[c,] = summary(lm(beta_y_ivw ~ beta_x_ivw[,c] -1))$coefficients
    }
    print("ivw:")
    print(ivw[order(ivw[,4], decreasing = F),]) 

    ## MR-BMA prior = 0.5
    # MR-BMA
    bma_input=new("mvMRInput", betaX = as.matrix(beta_x_ivw), betaY = as.matrix(beta_y_ivw), snps=list_snp, exposure=list_exposure, outcome = "mdd")
    bma_output=summarymvMR_SSS(bma_input,kmin=length(items)-1,kmax=length(items)-1, prior_prob=0.5)   

    # What are the next best individual models?
    bma_bestmodel = sss.report.best.model(bma_output, top = 20,  write.out = TRUE, csv.file.name="best.model_sigma05.csv")
    print(i)
    print(bma_bestmodel)    

    # MR-BMA output: Marginal inclusion probabilities and model-averaged effect (MACE) for each risk factor 
    bma_marginal = sss.report.mr.bma(bma_output, top = length(items)-1, write.out = FALSE)
    print(bma_marginal)
    bma_marginal_tmp=data.frame(bma_marginal)
    colnames(bma_marginal_tmp)[colnames(bma_marginal_tmp)=="rf"]="exposure"
    bma_marginal_tmp$outcome=i
    
    # empirical p-values
    bma_permute = create.permutations(bma_output, nrepeat = 100000, save.matrix=TRUE, file.name = "permutation_mrBMA.csv")
    empirical_p = calculate.p(bma_output, bma_permute)
    empirical_p_tmp=data.frame(empirical_p)
    empirical_p_tmp$exposure=row.names(empirical_p_tmp)
    rownames(empirical_p_tmp)=c()
    
    comb=merge(bma_marginal_tmp,empirical_p_tmp,by="exposure")
    comb=comb[order(-comb$mip),]
    final=comb[,c("outcome","exposure","marginal.inclusion","average.effect","pval","fdr")]
    final$exposure=str_replace_all(final$exposure,"beta_","")
    
    write.table(final,paste0(path_output,"/mvmr/out_",i,"_exp_",str_replace_all(paste(list_exposure,collapse="_"), "beta_", ""),".csv"),row.names=FALSE,quote=FALSE,sep="\t")     
  }
}



