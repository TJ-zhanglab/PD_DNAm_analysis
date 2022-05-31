library(survival)
library(survminer)
library(coxme)
library(ggplot2)

#########################fig1A######baseline_all####################################
AOA1_AA1_blood <- read.csv('../DNAm_file/AOA1-AO1.csv')
AOA1_AA1.cox <- coxph(Surv(AOA1, Status1) ~ AA +sex+relatedness +CD8T + CD4T +Bcell +Gran, data = AOA1_AA1_blood)
AOA1_AA1_with <- with(AOA1_AA1_blood, Surv(AOA1, Status1 == 1))
AOA1_AA1_sf <- survfit(AOA1_AA1_with ~ AA, data = AOA1_AA1_blood, conf.type = "log-log")
AOA1_AA1_pic<-ggsurvplot( AOA1_AA1_sf, data=AOA1_AA1_blood, conf.int=FALSE, legend.labs=c("Slow (AA<-3)","Normal (-3<AA<3)","Fast (AA>3)"),xlab="Age of symptom onset", ylab="Cumulative incidence", legend.title=c("DNAm-age acceleration"), palette = c("#009E73", "#56B4E9","#D55E00"), risk.table=TRUE,legend=c(0.2,0.9), fun = function(x) {1 - x}, break.x.by=10, xlim = c(30,80))


########################fig1B###########3-year-follow-up####################################
AOA1_AA1_blood <- read.csv('AOA2-AO2.csv')
AOA1_AA1.cox <- coxph(Surv(AOA2, Status2) ~ AA +sex+relatedness +CD8T + CD4T +Bcell +Gran, data = AOA1_AA1_blood)
OA1_AA1_with <- with(AOA1_AA1_blood, Surv(AOA2, Status2 == 1))
AOA1_AA1_sf <- survfit(AOA1_AA1_with ~ AA, data = AOA1_AA1_blood, conf.type = "log-log")
AOA1_AA1_pic<-ggsurvplot( AOA1_AA1_sf, data=AOA1_AA1_blood, conf.int=FALSE, legend.labs=c("Slow (AA<-3)","Normal (-3<AA<3)","Fast (AA>3)"),xlab="Age of symptom onset", ylab="Cumulative incidence", legend.title=c("DNAm-age acceleration"), palette = c("#009E73", "#56B4E9","#D55E00"), risk.table=TRUE,legend=c(0.2,0.9), fun = function(x) {1 - x}, break.x.by=10, xlim = c(30,80),ylim = c(0,0.6))


#####################fig2A#############baseline-91-manifesting#########################################################
AOA1_AA1_blood <- read.csv('LRRK2_91_baseline.csv')
AOA1_AA1.cox <- coxph(Surv(AO,Status1) ~ AA+interval +sex +CD8T + CD4T +Bcell +Gran, data = AOA1_AA1_blood)
AOA1_AA1_with <- with(AOA1_AA1_blood, Surv(AO,Status1==1))
AOA1_AA1_sf <- survfit(AOA1_AA1_with ~ AA, data = AOA1_AA1_blood, conf.type = "log-log")
AOA1_AA1_pic<-ggsurvplot( AOA1_AA1_sf, data=AOA1_AA1_blood, conf.int=FALSE, legend.labs=c("Slow (AA<-3)","Normal (-3<AA<3)","Fast (AA>3)"),xlab="Age of symptom onset", ylab="Cumulative incidence", legend.title=c("DNAm-age acceleration"), palette = c("#009E73", "#56B4E9","#D55E00"), risk.table=TRUE,legend=c(0.2,0.9), fun = function(x) {1 - x}, break.x.by=10, xlim = c(30,80))


#####################fig2B############96-idioPD###############################################################################
AOA1_AA1_blood <- read.csv('AO-AA-96ipd.csv')
AOA1_AA1.cox <- coxph(Surv(AO,status) ~ AA+interval+sex+CD8T + CD4T +Bcell +Gran, data = AOA1_AA1_blood)
AOA1_AA1_with <- with(AOA1_AA1_blood, Surv(AO,status==1))
AOA1_AA1_sf <- survfit(AOA1_AA1_with ~ AA, data = AOA1_AA1_blood, conf.type = "log-log")
AOA1_AA1_pic<-ggsurvplot( AOA1_AA1_sf, data=AOA1_AA1_blood, conf.int=FALSE, legend.labs=c("Slow (AA<-3)","Normal (-3<AA<3)","Fast (AA>3)"),xlab="Age of symptom onset", ylab="Cumulative incidence", legend.title=c("DNAm-age acceleration"), palette = c("#009E73", "#56B4E9","#D55E00"), risk.table=TRUE,legend=c(0.2,0.9), fun = function(x) {1 - x}, break.x.by=10, xlim = c(20,80))

#######################fig3A#############AO-AA-baseline-linear##########################################
AO_AA1_blood <- read.csv('AO-AA1.csv')
AO_AA1 <- lm(Age_symptom_onset ~ AA1+GENDER+interval,data = AO_AA1_blood)
AO_AA1_pic1 <- ggplot(data=AO_AA1_blood, aes(x=AA1, y=Age_symptom_onset)) +
  geom_point()+labs(x='DNAm Age acceleration at baseline of sample collection',y='Age of onset')
AO_AA1_pic2<-geom_abline(intercept = coef(AO_AA1)[1],slope = coef(AO_AA1)[2], color = "blue") 
AO_AA1_pic<-AO_AA1_pic1+AO_AA1_pic2+theme_classic()+scale_x_continuous(breaks=seq(-20,30,5))+scale_y_continuous(breaks=seq(40,100,10))

#########################fig3B#############AO-AA-3-year-follow-up#########################################
AO_AA2_blood <- read.csv('AO-AA2.csv')
AO_AA2 <- lm(Age_symptom_onset ~ AA2+GENDER+interval,data = AO_AA2_blood)
AO_AA2_pic1 <- ggplot(data=AO_AA2_blood, aes(x=AA2, y=Age_symptom_onset)) +
  geom_point()+labs(x='DNAm Age acceleration at the 3-year time point',y='Age of onset')
AO_AA2_pic2<-geom_abline(intercept = coef(AO_AA2)[1],slope = coef(AO_AA2)[2], color = "blue") 
AO_AA2_pic<-AO_AA2_pic1+AO_AA2_pic2+theme_classic()+scale_x_continuous(breaks=seq(-20,20,5))+scale_y_continuous(breaks=seq(10,90,10))

##########################fig3C########AO-AA-96-idiopathic PD####################################################
AO_AA3_blood <- read.csv('AO-AA-96ipd.csv')
AO_AA3 <- lm(AO ~ AA1+sex+interval,data = AO_AA3_blood)
AO_AA3_pic1 <- ggplot(data=AO_AA3_blood, aes(x=AA1, y=AO)) +
  geom_point()+labs(x='DNAm Age acceleration in idiopathic PD cases',y='Age of onset')
AO_AA3_pic2<-geom_abline(intercept = coef(AO_AA3)[1],slope = coef(AO_AA3)[2], color = "blue")
AO_AA3_pic<-AO_AA3_pic1+AO_AA3_pic2+theme_classic()+scale_x_continuous(limits=c(-15,20),breaks=seq(-20,30,5))+scale_y_continuous(limits=c(10,90),breaks=seq(10,90,10))

