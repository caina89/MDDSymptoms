library(data.table)
library(tidyr)

### ------
# paired t-test for rG 4 MDD symptom vs. 4 external phenotype
rm(list=ls())

# read data
data=fread("/no-backup/lianyun/20220222_plot_mddsymptomspaper/data/fig3_rg_lifetime_lencurrent_4in1.csv")
data=data.frame(data)

df=cbind(
    c("item","subjectivewellbeing.general","insomnia.baseline","selfworth.subjectivewellbeing.lifemeaningful","selfharm.notworthliving"),
    c("rg_worst",data[1,2],data[3,2],data[5,2],data[7,2]),
    c("rg_phq",data[2,2],data[4,2],data[6,2],data[8,2])
)

colnames(df)=df[1,]
df=df[-1,]
df=data.table(df)
df$rg_worst=as.numeric(df$rg_worst)
df$rg_phq=as.numeric(df$rg_phq)

# shift A1-subjectivewellbeing and A7-meaninginlife into positive, 
# it is possible, because, A1 sad mood correspond to bad subjective wellbeing. Therefore, good mood would mean good subjective wellbeing. 
# so i can shift the direction of numbers to positive.
# the same for A7 and lifemeaningful: feeling of worthlessness means life is not meaningful -> therefore, not feeling worthless means life is meaningful
# so i can shift this number as well
# note that this does not apply to all situations
df$rg_worst=abs(df$rg_worst)
df$rg_phq=abs(df$rg_phq)

t.test(x=df$rg_worst,y=df$rg_phq,alternative="less",mu=0,paired=TRUE,var.equal=FALSE,conf.level=0.95)  
### ------


### ------
# paired t-test for rG neuro, anxiety, trauma, stress

rm(list=ls())

item="stress" # neuro, anxiety, trauma, stress

if (item=="neuro") {
    data=fread(paste0("/path/","fig3c_heatmap_neuroticism",".csv"))
} else if (item == "anxiety") {
    data=fread(paste0("/path/","fig3b_heatmap_anxiety",".csv"))
} else if (item == "trauma") {
    data=fread(paste0("/path/","fig3d_heatmap_trauma",".csv"))
} else if (item == "stress") {
    data=fread(paste0("/path/","fig3d_heatmap_stress",".csv"))
}

data=data[-which(data$symp %like% "A5"),]
data=data[,c("symp","rg","externalsymp")]
data$tool=data$symp
data$tool=gsub("A.*","",data$tool)
data$symp=gsub("WorstEpisode","",data$symp)
data$symp=gsub("PHQ9","",data$symp)
data=data %>% spread(tool,rg) 

t.test(x=data$PHQ9,y=data$WorstEpisode,alternative="greater",mu=0,paired=TRUE,var.equal=FALSE,conf.level=0.95)  
### ------



### ------
# t-test UVMR 4 external phenotype to MDD symptom
rm(list=ls())

data=fread("/path/summary_uvmr_sympvs4.csv")
data=data.frame(data)
data=data[-which(data$exp %like% "4a"),]
data=data[-which(data$exp %like% "4b"),]
data=data[-which(data$exp %like% "4c"),]
data=data[-which(data$out %like% "4a"),]
data=data[-which(data$out %like% "4b"),]
data=data[-which(data$out %like% "4c"),]

df=NULL
items=c("subjectivewellbeing.general", "insomnia.baseline", "subjectivewellbeing.lifemeaningful", "selfharm.notworthliving")
#items=c("subjectivewellbeing.general", "subjectivewellbeing.lifemeaningful")
#items=c("insomnia.baseline", "selfharm.notworthliving")

use="beta" # c("beta","or")
for (exp in items) {
    tmp=data[which(data$exp==exp & data$method=="IVW"),]
    df=rbind(df,tmp)
}

df$data=gsub("A.*","",df$out)
df$symp=gsub(".*A","",df$out)
df$symp=paste0("A",df$symp)

if (use=="beta") {
    df=df[,c("exp","symp","data","estimate_beta")]
    df=df %>% spread(data,estimate_beta) 
} else if (use=="or") {
    df$or=exp(df$estimate_beta)
    df=df[,c("exp","symp","data","or")]
    df=df %>% spread(data,or)
}

dfpositive=df
dfpositive$lencurrent=abs(dfpositive$lencurrent)
dfpositive$lifetime=abs(dfpositive$lifetime)
t.test(x=dfpositive$lencurrent,y=dfpositive$lifetime,alternative="greater",mu=0,paired=TRUE,var.equal=FALSE,conf.level=0.95)  

### ------




### ------
# uvmr symp to/ from neuro, anxiety, trauma, stress


path_input="/no-backup/lianyun/20220222_plot_mddsymptomspaper/tmp"
title_vec=c("neuro","anxiety","trauma","stress")
#title_vec=c("trauma") 

# for trauma, it looks like many significant arrows for phq in the plot, but only a few for worst.
# but it's not significant in the t-test
# because i do t-test for all estimated beta A1-A9, while in the plot only those that survived bf qvalue are shown

# filename patterns:
# paste0("/","causalplot_phq_","neuro","_uvmr.csv")  # phq9 symptoms
# paste0("/","causalplot_worst_","neuro","_uvmr.csv") # worstepi symptoms

# function: from the corresponding summary table, select useful rows and columns
# in_title: neuro, anxiety or trauma
# in_mode: mdd2symp or symp2mdd
# in_worst: worst symptoms
# in_phq: phq symptoms

get_mr_stats = function(in_title, in_mode, in_worst, in_phq){
    
    if(in_mode %in% c("mdd2symp", "symp2mdd")){
        if(in_mode == "mdd2symp"){
            worst_tmp=in_worst[in_worst$exp %like% in_title,] # search title in exp, so it's mdd2symp 
            phq_tmp=in_phq[in_phq$exp %like% in_title,]

            worst_tmp$definition="worst"
            worst_tmp=worst_tmp[,c("estimate_beta","definition")] # we want to compare beta
            
            phq_tmp$definition="phq"
            phq_tmp=phq_tmp[-c(5),] # worst has no symp A5 so we remove it
            phq_tmp=phq_tmp[,c("estimate_beta","definition")]

            return(list(worst = worst_tmp, phq = phq_tmp))
            
        }else if(in_mode == "symp2mdd"){
            worst_tmp=in_worst[in_worst$out %like% in_title,] # search title in out, so it's symp2mdd
            phq_tmp=in_phq[in_phq$out %like% in_title,]
            
            worst_tmp$definition="worst"
            worst_tmp=worst_tmp[,c("estimate_beta","definition")]
            
            phq_tmp$definition="phq"
            phq_tmp=phq_tmp[-c(5),] 
            phq_tmp=phq_tmp[,c("estimate_beta","definition")]            

            return(list(worst = worst_tmp, phq = phq_tmp))
        }
    }else{
        return(NA)
    }    
}


# functions 
get_list_model = function(in_title, in_alternative, in_mode){
    # read in data
    worst=fread(paste0(path_input,"/","causalplot_worst_",in_title,"_uvmr.csv"))
    worst=data.frame(worst)
    #worst=worst[which(worst$bf_qvalue<0.05),]
    phq=fread(paste0(path_input,"/","causalplot_phq_",in_title,"_uvmr.csv"))
    phq=data.frame(phq)
    #phq=phq[which(phq$bf_qvalue<0.05),]

    # select useful rows and columns using get_mr_stats function
    mr_stats = get_mr_stats(in_title = in_title, in_mode = in_mode, in_worst = worst, in_phq = phq)
    worst_tmp = mr_stats$worst
    phq_tmp = mr_stats$phq

    # put worst and phq data in the same table
    combine=rbind(worst_tmp,phq_tmp)
    #combine=combine[which(combine$bf_qvalue<0.05),] # filter out non-significant results
    combine

    # run paired t-test
    #t_model = t.test(formula=combine$estimate_beta~combine$definition,
    #          alternative=in_alternative, mu=0,
    #          paired=TRUE,var.equal=FALSE,conf.level=0.95)
    
    # using this, you dont really know x is greater or lesser than y, so i changed to the function below
    
    worst_tmp=worst_tmp[,c("estimate_beta")]
    phq_tmp=phq_tmp[,c("estimate_beta")]
    
    t_model = t.test(x=worst_tmp,y=phq_tmp,
              alternative=in_alternative, mu=0,
              paired=TRUE,var.equal=FALSE,conf.level=0.95)    
    # to test whether worst is greater or lesser than phq
    
    return(t_model)
}



# compare three parameters: less, greater and two.sided

##########
# mdd to symp
list_tmodel_less_mdd2symp = sapply(title_vec, FUN = function(title){
   return(get_list_model(in_title = title, in_alternative = "less", in_mode = "mdd2symp"))
},simplify = F, USE.NAMES = T)

list_tmodel_greater_mdd2symp = sapply(title_vec, FUN = function(title){
   return(get_list_model(in_title = title, in_alternative = "greater", in_mode = "mdd2symp"))
},simplify = F, USE.NAMES = T)

list_tmodel_2sided_mdd2symp = sapply(title_vec, FUN = function(title){
   return(get_list_model(in_title = title, in_alternative = "two.sided", in_mode = "mdd2symp"))
},simplify = F, USE.NAMES = T)
##########


##########
# symp to mdd
list_tmodel_less_symp2mdd = sapply(title_vec, FUN = function(title){
   return(get_list_model(in_title = title, in_alternative = "less", in_mode = "symp2mdd"))
},simplify = F, USE.NAMES = T)

list_tmodel_greater_symp2mdd = sapply(title_vec, FUN = function(title){
   return(get_list_model(in_title = title, in_alternative = "greater", in_mode = "symp2mdd"))
},simplify = F, USE.NAMES = T)

list_tmodel_2sided_symp2mdd = sapply(title_vec, FUN = function(title){
   return(get_list_model(in_title = title, in_alternative = "two.sided", in_mode = "symp2mdd"))
},simplify = F, USE.NAMES = T)
##########


# generate results table, showing overall t-test results
t_test_df_mdd2symp = do.call(rbind, lapply(title_vec, FUN = function(x){
    
    col_cohort = rep(x, 4)
    
    col_parameter = c('mean_diff', 'df', 'tstats', 'pvalue')
    
    col_greater = c(list_tmodel_greater_mdd2symp[[x]]$estimate, list_tmodel_greater_mdd2symp[[x]]$parameter, 
                                list_tmodel_greater_mdd2symp[[x]]$statistic, list_tmodel_greater_mdd2symp[[x]]$p.value)
    col_2sided = c(list_tmodel_2sided_mdd2symp[[x]]$estimate, list_tmodel_2sided_mdd2symp[[x]]$parameter, 
                            list_tmodel_2sided_mdd2symp[[x]]$statistic, list_tmodel_2sided_mdd2symp[[x]]$p.value)
    col_less = c(list_tmodel_less_mdd2symp[[x]]$estimate, list_tmodel_less_mdd2symp[[x]]$parameter, 
                            list_tmodel_less_mdd2symp[[x]]$statistic, list_tmodel_less_mdd2symp[[x]]$p.value)
    
    
    curr_df = cbind(cohort = col_cohort, 
                    parameter = col_parameter,
                    greater = col_greater, 
                    two_sided = col_2sided,
                    less = col_less)
    
    rownames(curr_df) = col_parameter
    
    return(curr_df)
}))



t_test_df_symp2mdd = do.call(rbind, lapply(title_vec, FUN = function(x){
    
    col_cohort = rep(x, 4)
    
    col_parameter = c('mean_diff', 'df', 'tstats', 'pvalue')
    
    col_greater = c(list_tmodel_greater_symp2mdd[[x]]$estimate, list_tmodel_greater_symp2mdd[[x]]$parameter, 
                                list_tmodel_greater_symp2mdd[[x]]$statistic, list_tmodel_greater_symp2mdd[[x]]$p.value)
    col_2sided = c(list_tmodel_2sided_symp2mdd[[x]]$estimate, list_tmodel_2sided_symp2mdd[[x]]$parameter, 
                            list_tmodel_2sided_symp2mdd[[x]]$statistic, list_tmodel_2sided_symp2mdd[[x]]$p.value)
    col_less = c(list_tmodel_less_symp2mdd[[x]]$estimate, list_tmodel_less_symp2mdd[[x]]$parameter, 
                            list_tmodel_less_symp2mdd[[x]]$statistic, list_tmodel_less_symp2mdd[[x]]$p.value)
    
    
    curr_df = cbind(cohort = col_cohort, 
                    parameter = col_parameter,
                    greater = col_greater, 
                    two_sided = col_2sided,
                    less = col_less)
    
    rownames(curr_df) = col_parameter
    
    return(curr_df)
}))


print("mdd2symp") # we look into mdd2symp, which is to test if items lead to symptoms or not
t_test_df_mdd2symp

print("symp2mdd")
t_test_df_symp2mdd

### ------


### ------
# uvmr symptoms from/to pgc, ipsych2012, ipsych2015

datasets=c("pgcmdd","ipsych2012","ipsych2015i")

worst=NULL
for (title in datasets){
    tmp=fread(paste0(path_input,"/causalplot_lifetimesymp_",title,"_bothmr_uvmr.csv"))
    tmp=data.frame(tmp)
    tmp=tmp[tmp$exp %like% title,]
    tmp$definition="worst"
    tmp$out=gsub(".*A","",tmp$out)
    tmp$item=paste0(tmp$exp,"_A",tmp$out)
    tmp=tmp[,c("item","definition","estimate_beta")]
    worst=rbind(worst,tmp)
}


phq=NULL
for (title in datasets){
    tmp=fread(paste0(path_input,"/causalplot_lencurrentsymp_",title,"_bothmr_uvmr.csv"))
    tmp=data.frame(tmp)
    tmp=tmp[tmp$exp %like% title,]
    tmp$definition="phq"
    tmp$out=gsub(".*A","",tmp$out)
    tmp$item=paste0(tmp$exp,"_A",tmp$out)
    tmp=tmp[,c("item","definition","estimate_beta")]
    phq=rbind(phq,tmp)
}
phq=phq[-which(phq$item %like% "A5"),]

data=rbind(worst,phq)
data= data%>% spread(definition,estimate_beta)

t.test(x=data$worst,y=data$phq,alternative="less",mu=0,paired=TRUE,var.equal=FALSE,conf.level=0.95)  
### ------





