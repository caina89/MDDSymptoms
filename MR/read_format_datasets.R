# ------
# read from different path
# ------
read_lifetime=function(symp_tmp){
  symp_mediate=substring(symp_tmp,9,11)
  data_tmp=vroom(paste0("/path/GwasSumstats_FullSample_Pheno_",symp_mediate))
  data_tmp$symp=paste0("lifetime",symp_mediate)
  return(data_tmp)
}

read_lencurrent=function(symp_tmp){
  symp_mediate=substring(symp_tmp,11,12)
  data_tmp=vroom(paste0("/path/GwasSumstats_LenCurrent_",symp_mediate))
  data_tmp$symp=paste0("lencurrent",symp_mediate)
  return(data_tmp)
}


read_anxiety=function(symp_tmp){
  data_tmp=vroom(paste0("/path/allchr.autocomplete.anxiety.recent.",symp_tmp,".glm.maf05.logistic"))
  data_tmp$symp=paste0("anxiety_",symp_tmp)
  return(data_tmp)
}


read_neuroticism=function(symp_tmp){
  if (symp_tmp=="subjectivewellbeing.lifemeaningful"){
    data_tmp=vroom(paste0("/path/allchr.subjectivewellbeing.lifemeaningful.glm.maf05.linear"))
    data_tmp$symp=symp_tmp
  }
  else if(symp_tmp=="selfharm.notworthliving"){
    data_tmp=vroom(paste0("/path/allchr.selfharm.notworthliving.glm.maf05.logistic"))
    data_tmp$symp=symp_tmp
  }
  else {
    data_tmp=vroom(paste0("/path/allchr.",symp_tmp,".glm.maf05.logistic"))
    data_tmp$symp=symp_tmp    
  }
  
  return(data_tmp)
}


read_stress=function(symp_tmp){
  data_tmp=vroom(paste0("/path/allchr.",symp_tmp,".glm.maf05.logistic"))
  data_tmp$symp=paste0("stress_",symp_tmp)
  return(data_tmp)
}


read_externalsymp=function(symp_tmp){
  if (symp_tmp %in% c("bmi","educationage.baseline","income.baseline","insomnia.baseline","neuroticismscore.baseline","subjectivewellbeing.general")){
    data_tmp=vroom(paste0("/path/allchr.original.",symp_tmp,".glm.maf05.linear"))
    data_tmp$symp=paste0("original_",symp_tmp)
  }
  else {
    data_tmp=vroom(paste0("/path/allchr.original.",symp_tmp,".glm.maf05.logistic"))
    data_tmp$symp=paste0("original_",symp_tmp)    
  }

  return(data_tmp)
}


read_lifetimemdd=function(symp_tmp){
  data_tmp=vroom(paste0("/path/gwassumstats_matched_ukbmdd.txt"))
  data_tmp$symp=paste0("lifetimeMDD")
  return(data_tmp)
}

read_pgc29=function(symp_tmp){
  data_tmp=vroom(paste0("/path/gwassumstats_matched_pgc29.txt"))
  data_tmp$symp=paste0("pgc29")
  return(data_tmp)
}

read_ipsych2012=function(symp_tmp){
  data_tmp=vroom(paste0("/path/gwassumstats_matched_ipsych2012_maf.txt"))
  data_tmp$symp=paste0("ipsych2012")
  return(data_tmp)
}

read_ipsych2015i=function(symp_tmp){
  data_tmp=vroom(paste0("path/gwassumstats_matched_ipsych2015i_maf.txt"))
  data_tmp$symp=paste0("ipsych2015i")
  return(data_tmp)
}


# read in downsample, a bit tricky
# read_dslifetime_orig including A3 A4 A6 A7 A8, where lencurrent was downsampled, 
read_dslifetime_orig=function(symp_tmp){
  symp_mediate=substring(symp_tmp,12,13)
  data_tmp=vroom(paste0("/path/GwasSumstats_FullSample_Pheno_",symp_mediate))
  data_tmp$symp=paste0("ds_lifetime",symp_mediate)
  return(data_tmp)
}

read_dslifetime_ds=function(symp_tmp){
  symp_mediate=substring(symp_tmp,12,13)
  data_tmp=vroom(paste0("/path/WorstEpisode",symp_mediate,".logistic"))
  data_tmp$symp=paste0("ds_lifetime",symp_mediate)
  return(data_tmp)
}


read_dslencurrent_orig=function(symp_tmp){
  symp_mediate=substring(symp_tmp,14,15)
  data_tmp=vroom(paste0("/path/GwasSumstats_LenCurrent_",symp_mediate))
  data_tmp$symp=paste0("ds_lencurrent",symp_mediate)
  return(data_tmp)
}


read_dslencurrent_ds=function(symp_tmp){
  symp_mediate=substring(symp_tmp,14,15)
  data_tmp=vroom(paste0("/path/PHQ9",symp_mediate,".logistic"))
  data_tmp$symp=paste0("ds_lencurrent",symp_mediate)
  return(data_tmp)
}




# read in downsample_constentprevalence, a bit tricky
# worst: A1 A2 A9 downsample_constentprev, others original
# phq: A3 A4 A6 A7 A8 downsample_constentprev, others original
read_cpdslifetime_orig=function(symp_tmp){
  symp_mediate=substring(symp_tmp,14,15)
  data_tmp=fread(paste0("/path/GwasSumstats_FullSample_Pheno_",symp_mediate))
  data_tmp=data.frame(data_tmp)
  data_tmp$symp=paste0("cpds_lifetime",symp_mediate)
  return(data_tmp)
}

read_cpdslifetime_cpds=function(symp_tmp){
  symp_mediate=substring(symp_tmp,14,15)
  data_tmp=fread(paste0("/path/constprev.WorstEpisode",symp_mediate,".logistic"))
  data_tmp=data.frame(data_tmp)
  data_tmp$symp=paste0("cpds_lifetime",symp_mediate)
  return(data_tmp)
}


read_cpdslencurrent_orig=function(symp_tmp){
  symp_mediate=substring(symp_tmp,16,17)
  data_tmp=fread(paste0("/path/GwasSumstats_LenCurrent_",symp_mediate))
  data_tmp=data.frame(data_tmp)
  data_tmp$symp=paste0("cpds_lencurrent",symp_mediate)
  return(data_tmp)
}


read_cpdslencurrent_cpds=function(symp_tmp){
  symp_mediate=substring(symp_tmp,16,17)
  data_tmp=fread(paste0("/path/constprev.PHQ9",symp_mediate,".logistic"))
  data_tmp=data.frame(data_tmp)
  data_tmp$symp=paste0("cpds_lencurrent",symp_mediate)
  return(data_tmp)
}









# ------
# decide which read() to be used
# ------
read_data=function(symp_tmp){
  if (symp_tmp %in% c("lifetimeA1","lifetimeA2","lifetimeA3","lifetimeA3a","lifetimeA3b","lifetimeA3c","lifetimeA4","lifetimeA4a","lifetimeA4b","lifetimeA4c","lifetimeA6","lifetimeA7","lifetimeA8","lifetimeA9")) {
    tmp = read_lifetime(symp_tmp)
  }
  else if (symp_tmp %in% c("lencurrentA1","lencurrentA2","lencurrentA3","lencurrentA4","lencurrentA5","lencurrentA6","lencurrentA7","lencurrentA8","lencurrentA9")) {
    tmp = read_lencurrent(symp_tmp)
  }
  else if (symp_tmp %in% c("anxiety","cantstopworrying","foreboding","irritability","restlessness","troublerelaxing","worrytoomanythings")) {
    tmp = read_anxiety(symp_tmp)
  }
  else if (symp_tmp %in% c("neuroticism.anxious","neuroticism.fedup","neuroticism.guilt","neuroticism.irritability","neuroticism.loneliness","neuroticism.miserableness","neuroticism.moodswings","neuroticism.nerves","neuroticism.nervous","neuroticism.risktaking","neuroticism.sensitivity","neuroticism.tense","neuroticism.worrytoolong","selfharm.notworthliving","subjectivewellbeing.lifemeaningful")) {
    tmp = read_neuroticism(symp_tmp)
  }
  else if (symp_tmp %in% c("stressbinary","trauma")) {
    tmp = read_stress(symp_tmp)
  }
  else if (symp_tmp %in% c("bmi","educationage.baseline","income.baseline","insomnia.baseline","neuroticismscore.baseline","subjectivewellbeing.general","anxiety.diagnosis","eversamesexpartner.baseline","eversmoked.baseline","mentalhealth.distress","selfharm.everselfharm","snoring.baseline")) {
    tmp = read_externalsymp(symp_tmp)
  }
  else if (symp_tmp %in% c("lifetimemdd","lifetimeMDD")) {
    tmp = read_lifetimemdd(symp_tmp)
  }
  else if (symp_tmp %in% c("ds_lifetimeA1","ds_lifetimeA2","ds_lifetimeA9")) {
    tmp = read_dslifetime_ds(symp_tmp)
  }
  else if (symp_tmp %in% c("ds_lifetimeA3","ds_lifetimeA4","ds_lifetimeA6","ds_lifetimeA7","ds_lifetimeA8")) {
    tmp = read_dslifetime_orig(symp_tmp)
  }
  else if (symp_tmp %in% c("ds_lencurrentA3","ds_lencurrentA4","ds_lencurrentA6","ds_lencurrentA7","ds_lencurrentA8")) {
    tmp = read_dslencurrent_ds(symp_tmp)
  }
  else if (symp_tmp %in% c("ds_lencurrentA1","ds_lencurrentA2","ds_lencurrentA9")) {
    tmp = read_dslencurrent_orig(symp_tmp)
  }
  else if (symp_tmp %in% c("cpds_lifetimeA3","cpds_lifetimeA4","cpds_lifetimeA6","cpds_lifetimeA7","cpds_lifetimeA8")){
    tmp = read_cpdslifetime_orig(symp_tmp)
  }
  else if (symp_tmp %in% c("cpds_lifetimeA1","cpds_lifetimeA2","cpds_lifetimeA9")){
    tmp = read_cpdslifetime_cpds(symp_tmp)
  } 
  else if (symp_tmp %in% c("cpds_lencurrentA1","cpds_lencurrentA2","cpds_lencurrentA9")) {
    tmp = read_cpdslencurrent_orig(symp_tmp)
  } 
  else if (symp_tmp %in% c("cpds_lencurrentA3","cpds_lencurrentA4","cpds_lencurrentA6","cpds_lencurrentA7","cpds_lencurrentA8")){
    tmp = read_cpdslencurrent_cpds(symp_tmp)
  }
  else if (symp_tmp %in% c("pgc29")) {
    tmp = read_pgc29(symp_tmp)
  }
  else if (symp_tmp %in% c("ipsych2012")) {
    tmp = read_ipsych2012(symp_tmp)
  }
  else if (symp_tmp %in% c("ipsych2015i")) {
    tmp = read_ipsych2015i(symp_tmp)
  }
  return(tmp)
}





# ------
# format data based on whether it is exposure/ outcome/ or no specific format for UVMR
# ------
exposure_format=function(data_tmp) {
  # rsid  
  if ("ID" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ID"]="rsid"}
  if ("SNP" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SNP"]="rsid"}
  # or,beta
  if ("OR" %in% colnames(data_tmp)){data_tmp$beta_exp=log(data_tmp$OR)}
  if ("BETA" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETA"]="beta_exp"}
  # se
  if ("SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SE"]="se_exp"}
  if ("LOG(OR)_SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOG(OR)_SE"]="se_exp"}
  if ("BETASE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETASE"]="se_exp"}
  if ("LOGORSE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOGORSE"]="se_exp"}
  # allele1
  if ("A1" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1"]="a1_exp"} 
  if ("ALT" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ALT"]="a1_exp"}
  # allele2
  if ("A2" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A2"]="a2_exp"} 
  if ("REF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="REF"]="a2_exp"}
  if ("A0" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A0"]="a2_exp"}
  # maf
  if ("MAF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="MAF"]="maf_exp"}
  if ("A1FREQ" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1FREQ"]="maf_exp"}
  if ("FREQA1" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="FREQA1"]="maf_exp"}
  if ("maf" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="maf"]="maf_exp"}
  # p
  if ("P" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="P"]="p_exp"}
  if ("PVALUE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="PVALUE"]="p_exp"}
  # symp
  if ("symp" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="symp"]="symp_exp"}

  data_tmp=data_tmp %>% select("rsid","beta_exp","se_exp","p_exp","a1_exp","a2_exp","maf_exp","symp_exp")

  return(data_tmp)
}

outcome_format=function(data_tmp) {
  # rsid  
  if ("ID" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ID"]="rsid"}
  if ("SNP" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SNP"]="rsid"}
  # or,beta
  if ("OR" %in% colnames(data_tmp)){data_tmp$beta_out=log(data_tmp$OR)}
  if ("BETA" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETA"]="beta_out"}
  # se
  if ("SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SE"]="se_out"}
  if ("LOG(OR)_SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOG(OR)_SE"]="se_out"}
  if ("BETASE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETASE"]="se_out"}
  if ("LOGORSE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOGORSE"]="se_out"}
  # allele1
  if ("A1" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1"]="a1_out"} 
  if ("ALT" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ALT"]="a1_out"}
  # allele2
  if ("A2" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A2"]="a2_out"} 
  if ("REF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="REF"]="a2_out"}
  if ("A0" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A0"]="a2_out"}
  # maf
  if ("MAF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="MAF"]="maf_out"}
  if ("A1FREQ" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1FREQ"]="maf_out"}
  if ("FREQA1" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="FREQA1"]="maf_out"}
  if ("maf" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="maf"]="maf_out"}
  # p
  if ("P" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="P"]="p_out"}
  if ("PVALUE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="PVALUE"]="p_out"}
  # symp
  if ("symp" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="symp"]="symp_out"}

  data_tmp=data_tmp %>% select("rsid","beta_out","se_out","p_out","a1_out","a2_out","maf_out","symp_out")

  return(data_tmp)
}


mvmr_format=function(data_tmp,suffix) {
  # rsid  
  if ("ID" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ID"]="rsid"}
  if ("SNP" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SNP"]="rsid"}
  # or,beta
  if ("OR" %in% colnames(data_tmp)){data_tmp[,paste0("beta_",suffix)]=log(data_tmp$OR)}
  if ("BETA" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETA"]=paste0("beta_",suffix)}
  # se
  if ("SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SE"]=paste0("se_",suffix)}
  if ("LOG(OR)_SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOG(OR)_SE"]=paste0("se_",suffix)}
  if ("BETASE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETASE"]=paste0("se_",suffix)}
  if ("LOGORSE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOGORSE"]=paste0("se_",suffix)}
  # allele1
  if ("A1" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1"]=paste0("a1_",suffix)} 
  if ("ALT" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ALT"]=paste0("a1_",suffix)}
  # allele2
  if ("A2" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A2"]=paste0("a2_",suffix)} 
  if ("REF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="REF"]=paste0("a2_",suffix)}
  if ("A0" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A0"]=paste0("a2_",suffix)}
  # p
  if ("P" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="P"]=paste0("p_",suffix)}
  if ("PVALUE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="PVALUE"]=paste0("p_",suffix)}

  data_tmp=data_tmp %>% select("rsid",paste0("beta_",suffix),paste0("se_",suffix),paste0("p_",suffix),paste0("a1_",suffix),paste0("a2_",suffix))

  return(data_tmp)
}



general_format=function(data_tmp) {
  # rsid  
  if ("ID" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ID"]="rsid"}
  if ("SNP" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SNP"]="rsid"}
  # or,beta
  if ("OR" %in% colnames(data_tmp)){data_tmp$beta=log(data_tmp$OR)}
  if ("BETA" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETA"]="beta"}
  # se
  if ("SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="SE"]="se"}
  if ("LOG(OR)_SE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOG(OR)_SE"]="se"}
  if ("BETASE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="BETASE"]="se"}
  if ("LOGORSE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="LOGORSE"]="se"}
  # allele1
  if ("A1" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1"]="a1"} 
  if ("ALT" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="ALT"]="a1"}
  # allele2
  if ("A2" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A2"]="a2"} 
  if ("REF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="REF"]="a2"}
  if ("A0" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A0"]="a2"}
  # maf
  if ("MAF" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="MAF"]="maf"}
  if ("A1FREQ" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="A1FREQ"]="maf"} 
  if ("maf" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="maf"]="maf"} 
  # p
  if ("P" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="P"]="p"}
  if ("PVALUE" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="PVALUE"]="p"}
  # symp
  if ("symp" %in% colnames(data_tmp)){colnames(data_tmp)[colnames(data_tmp)=="symp"]="symp"}

  data_tmp=data_tmp %>% select("rsid","beta","se","p","a1","a2","maf","symp")

  return(data_tmp)
}










# ------
# ------
# ------
create_subdir=function(doit){
    if(file.exists(file.path(paste0(path_output,"/f_stats")))==FALSE){
      dir.create(file.path(paste0(path_output,"/f_stats")), showWarnings = FALSE)}
    if(file.exists(file.path(paste0(path_output,"/mr_presso")))==FALSE){
      dir.create(file.path(paste0(path_output,"/mr_presso")), showWarnings = FALSE)}
    if(file.exists(file.path(paste0(path_output,"/plot")))==FALSE){
      dir.create(file.path(paste0(path_output,"/plot")), showWarnings = FALSE)}  
    if(file.exists(file.path(paste0(path_output,"/rimage")))==FALSE){
      dir.create(file.path(paste0(path_output,"/rimage")), showWarnings = FALSE)}
    if(file.exists(file.path(paste0(path_output,"/snp")))==FALSE){
      dir.create(file.path(paste0(path_output,"/snp")), showWarnings = FALSE)}
    if(file.exists(file.path(paste0(path_output,"/uvmr")))==FALSE){
      dir.create(file.path(paste0(path_output,"/uvmr")), showWarnings = FALSE)}
    if(file.exists(file.path(paste0(path_output,"/mvmr")))==FALSE){
      dir.create(file.path(paste0(path_output,"/mvmr")), showWarnings = FALSE)}
}

