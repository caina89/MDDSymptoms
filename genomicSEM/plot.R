library(semPlot)
library(lavaan)

semPlotModel_GSEM=function(gsem.object=GWISoutput , est.label="STD_All"){ 
       object=gsem.object$results
       object$free=0
       numb=1:length(which(object$op!="~~"))
       object$free[which(object$op!="~~")]=numb
       varNames <- lavaanNames(object, type = "ov")
       factNames <- lavaanNames(object, type = "lv")
       factNames <- factNames[!factNames %in% varNames]
       n <- length(varNames)
       k <- length(factNames)
       if (is.null(object$label)) 
         object$label <- rep("", nrow(object))
       semModel <- new("semPlotModel")
       object$est <- object[,est.label]
       if (is.null(object$group)) 
         object$group <- ""
       semModel@Pars <- data.frame(label = object$label, lhs = ifelse(object$op == "~" | object$op == "~1", object$rhs, object$lhs), edge = "--", 
                                   rhs = ifelse(object$op == "~" | object$op == "~1", object$lhs, object$rhs), est = object$est, std = NA, group = object$group, 
                                   fixed = object$free==0, par = object$free, stringsAsFactors = FALSE)
       semModel@Pars$edge[object$op == "~~"] <- "<->"
       semModel@Pars$edge[object$op == "~*~"] <- "<->"
       semModel@Pars$edge[object$op == "~"] <- "~>"
       semModel@Pars$edge[object$op == "=~"] <- "->"
       semModel@Pars$edge[object$op == "~1"] <- "int"
       semModel@Pars$edge[grepl("\\|", object$op)] <- "|"
       semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge), 
                                            -(3:4)]
       semModel@Pars <- semModel@Pars[!object$op %in% c(":=", "<", 
                                                        ">", "==", "|", "<", ">"), ]
       semModel@Vars <- data.frame(name = c(varNames, factNames), 
                                   manifest = c(varNames, factNames) %in% varNames, exogenous = NA, 
                                   stringsAsFactors = FALSE)
       semModel@ObsCovs <- list()
       semModel@ImpCovs <- list()
       semModel@Computed <- FALSE
       semModel@Original <- list(object)
       return(semModel)
       
}
## plot CFA for EFA1
symptoms=c("WE","PHQ9","PHQ9SkipA1A2")
cohorts=c("LifetimeMDD","PGC29","iPSYCH2012","iPSYCH2015i")
for (symptom in symptoms){
       print(symptom)
       for (cohort in cohorts){
              print(cohort)
              file=paste0(symptom,".",cohort,".genomicSEM.gCFAofEFA1")
              load(paste0(file,".RData"))
              pdf(paste0(file,".pdf"),height=3,width=3,useDingbats=FALSE)
              semPaths(semPlotModel_GSEM(gCFAofEFA1),layout="tree2",style="OpenMx",whatLabels = "est")
              dev.off()        
       }
}
## plot CFA for EFA2
symptoms=c("WE","PHQ9","PHQ9SkipA1A2")
cohorts=c("LifetimeMDD","PGC29","iPSYCH2012","iPSYCH2015i")
for (symptom in symptoms){
       print(symptom)
       for (cohort in cohorts){
              print(cohort)
              file=paste0(symptom,".",cohort,".genomicSEM.gCFAofEFA2")
              load(paste0(file,".RData"))
              pdf(paste0(file,".pdf"),height=3,width=3,useDingbats=FALSE)
              semPaths(semPlotModel_GSEM(gCFAofEFA2),layout="tree2",style="OpenMx",whatLabels = "est")
              dev.off()        
       }
}
## plot CFA for EFA3
symptoms=c("WE","PHQ9","PHQ9SkipA1A2")
cohorts=c("LifetimeMDD","PGC29","iPSYCH2012","iPSYCH2015i")
for (symptom in symptoms){
       print(symptom)
       for (cohort in cohorts){
              print(cohort)
              file=paste0(symptom,".",cohort,".genomicSEM.gCFAofEFA3")
              if (file.exists(paste0(file,".RData"))){
                     load(paste0(file,".RData"))
                     pdf(paste0(file,".pdf"),height=3,width=3,useDingbats=FALSE)
                     semPaths(semPlotModel_GSEM(gCFAofEFA3),layout="tree2",style="OpenMx",whatLabels = "est")
                     dev.off()
              }
       }
}


