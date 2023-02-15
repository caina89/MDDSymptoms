## get EFA 
library(psych)
library(GPArotation)
setwd("/no-backup/caina/biobank/genomicSEM/mungestats")

## WorstEpisode symptoms 
LDSCoutput=load("WE.LifetimeMDD.genomicSEM.RData")

EFA1<-fa(r=as.matrix(LDSCoutput$S), nfactors=1, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA2<-fa(r=as.matrix(LDSCoutput$S), nfactors=2, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA3<-fa(r=as.matrix(LDSCoutput$S), nfactors=3, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA4<-fa(r=as.matrix(LDSCoutput$S), nfactors=4, rotate="promax", SMC=FALSE, fm="minres",covar=T)

CFAofEFA1 <- 'F1=~NA*LifetimeMDD+A1+A2+A3+A4+A6+A7+A8+A9'
gCFAofEFA1=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA1, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA1$modelfit
#      chisq df      p_chisq     AIC       CFI      SRMR
# df 177.911 27 3.660276e-24 213.911 0.9696789 0.1271246
save(gCFAofEFA1, file="WE.LifetimeMDD.genomicSEM.gCFAofEFA1.RData")

CFAofEFA2 <- 	'F1=~NA*LifetimeMDD+A1+A2
				F2=~NA*A3+A4+A6+A7+A8+A9
				F1~~1*F1
				F2~~1*F2
				F1~~F2
				A1~~c*A1
				A8~~d*A8
				A9~~e*A9
				c>0.001
				d>0.001
				e>0.001'
gCFAofEFA2=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA2, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA2$modelfit
#       chisq df      p_chisq      AIC      CFI       SRMR
# df 91.82711 26 2.820498e-09 129.8271 0.986774 0.08462389
save(gCFAofEFA2, file="WE.LifetimeMDD.genomicSEM.gCFAofEFA2.RData")

CFAofEFA3 <- 'F1=~NA*LifetimeMDD+A1+A2
				F2=~NA*A4+A6+A7+A8+A9
				F3=~NA*LifetimeMDD+A3+A4
				## variances
				F1~~1*F1
				F2~~1*F2
				F3~~1*F3
				## covariances
				F1~~F2
				F2~~F3
				F1~~F3
				## residuals 
				LifetimeMDD~~a*LifetimeMDD
				A1~~b*A1
				A2~~c*A2
				A3~~d*A3
				A4~~e*A4
				A6~~f*A6
				A7~~g*A7
				A8~~h*A8
				A9~~i*A9
				a>0.001
				b>0.001
				c>0.001
				d>0.001
				e>0.001
				f>0.001
				g>0.001
				h>0.001
				i>0.001'

gCFAofEFA3=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA3, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA3$modelfit
#      chisq df     p_chisq     AIC       CFI       SRMR
# df 44.1436 22 0.003402846 90.1436 0.9955509 0.07665413
save(gCFAofEFA3, file="WE.LifetimeMDD.genomicSEM.gCFAofEFA3.RData")

## PHQ9 symptoms 
LDSCoutput=load("PHQ9.LifetimeMDD.genomicSEM.RData")

EFA1<-fa(r=as.matrix(LDSCoutput$S), nfactors=1, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA2<-fa(r=as.matrix(LDSCoutput$S), nfactors=2, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA3<-fa(r=as.matrix(LDSCoutput$S), nfactors=3, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA4<-fa(r=as.matrix(LDSCoutput$S), nfactors=4, rotate="promax", SMC=FALSE, fm="minres",covar=T)

CFAofEFA1 <- 'F1=~NA*LifetimeMDD+A1+A2+A3+A4+A5+A6+A7+A8+A9'
gCFAofEFA1=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA1, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA1$modelfit
#       chisq df      p_chisq      AIC       CFI       SRMR
# df 78.29752 35 3.721039e-05 118.2975 0.9852366 0.05882035
save(gCFAofEFA1, file="PHQ9.LifetimeMDD.genomicSEM.gCFAofEFA1.RData")

CFAofEFA2 <- 	'F1=~NA*LifetimeMDD+A1+A2+A3+A5+A7+A8+A9
				F2=~NA*A1+A2+A3+A4+A5+A6+A8
				F1~~1*F1
				F2~~1*F2
				F1~~F2
				A4~~c*A4
				A7~~d*A7
				A9~~e*A9
				c>0.001
				d>0.001
				e>0.001'
gCFAofEFA2=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA2, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA2$modelfit
#       chisq df   p_chisq      AIC       CFI       SRMR
# df 34.28931 29 0.2288274 86.28931 0.9981965 0.04056499
save(gCFAofEFA2, file="PHQ9.LifetimeMDD.genomicSEM.gCFAofEFA2.RData")

CFAofEFA3 <- 'F1=~NA*LifetimeMDD+A1+A2+A3+A7+A8+A9
				F2=~NA*A1+A2+A3+A4+A5+A6+A8
				F3=~NA*LifetimeMDD+A5
				## variances
				F1~~1*F1
				F2~~1*F2
				F3~~1*F3
				## covariances
				F1~~F2
				F2~~F3
				F1~~F3
				## residuals 
				LifetimeMDD~~a*LifetimeMDD
				A7~~g*A7
				A8~~h*A8
				A9~~i*A9
				a>0.001
				g>0.001
				h>0.001
				i>0.001'

gCFAofEFA3=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA3, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
#The primary model did not converge!


## PHQ9Skip symptoms
LDSCoutput=load("PHQ9SkipA1A2.LifetimeMDD.genomicSEM.RData")

EFA1<-fa(r=as.matrix(LDSCoutput$S), nfactors=1, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA2<-fa(r=as.matrix(LDSCoutput$S), nfactors=2, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA3<-fa(r=as.matrix(LDSCoutput$S), nfactors=3, rotate="promax", SMC=FALSE, fm="minres",covar=T)
EFA4<-fa(r=as.matrix(LDSCoutput$S), nfactors=4, rotate="promax", SMC=FALSE, fm="minres",covar=T)

CFAofEFA1 <- 'F1=~NA*LifetimeMDD+A1+A2+SkipA3+SkipA4+SkipA5+SkipA6+SkipA7+SkipA8+SkipA9'
gCFAofEFA1=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA1, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA1$modelfit
#      chisq df      p_chisq      AIC       CFI      SRMR
# df 66.6886 35 0.0009816996 106.6886 0.9850723 0.1619797
save(gCFAofEFA1, file="PHQ9SkipA1A2.LifetimeMDD.genomicSEM.gCFAofEFA1.RData")

CFAofEFA2 <- 	'F1=~NA*LifetimeMDD+A1+A2+SkipA3+SkipA4+SkipA5+SkipA6+SkipA8
				F2=~NA*LifetimeMDD+A1+A2+SkipA7+SkipA8+SkipA9
				F1~~1*F1
				F2~~1*F2
				F1~~F2
				LifetimeMDD~~a*LifetimeMDD
				SkipA4~~b*SkipA4
				SkipA5~~c*SkipA5
				SkipA6~~d*SkipA6
				SkipA7~~e*SkipA7
				SkipA9~~f*SkipA9
				c>0.001
				b>0.001
				c>0.001
				d>0.001
				e>0.001
				f>0.001'
gCFAofEFA2=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA2, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# > gCFAofEFA2$modelfit
#       chisq df     p_chisq      AIC       CFI      SRMR
# df 56.35188 30 0.002487156 106.3519 0.9875863 0.1170016
save(gCFAofEFA2, file="PHQ9SkipA1A2.LifetimeMDD.genomicSEM.gCFAofEFA2.RData")

CFAofEFA3 <- 'F1=~NA*LifetimeMDD+A1+A2+SkipA3+SkipA4+SkipA5+SkipA6
				F2=~NA*LifetimeMDD+A1+A2+SkipA7+SkipA8+SkipA9
				F3=~NA*SkipA3+SkipA4+SkipA8+SkipA9
				## variances
				F1~~1*F1
				F2~~1*F2
				F3~~1*F3
				## covariances
				F1~~F2
				F2~~F3
				F1~~F3
				## residuals 
				A1~~a*A1
				SkipA3~~b*SkipA3
				SkipA4~~c*SkipA4
				SkipA5~~d*SkipA5
				SkipA6~~e*SkipA6
				SkipA7~~f*SkipA7
				SkipA8~~g*SkipA8
				SkipA9~~h*SkipA9
				a>0.001
				b>0.001
				c>0.001
				d>0.001
				e>0.001
				f>0.001
				g>0.001
				h>0.001'
gCFAofEFA3=usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA3, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# >  gCFAofEFA3$modelfit
#       chisq df    p_chisq      AIC       CFI       SRMR
# df 41.61361 25 0.01976817 101.6136 0.9921737 0.08636107
save(gCFAofEFA3, file="PHQ9SkipA1A2.LifetimeMDD.genomicSEM.gCFAofEFA3.RData")


