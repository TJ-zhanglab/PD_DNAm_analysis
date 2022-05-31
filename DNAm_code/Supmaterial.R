library(ggplot2)
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
############################figS1A############################
four_na <- read.csv('fourtimepoint_nonmanifesting.csv')
sample <- four[1:115,3:6]
time_point<-c(1,2,3,4)
plot(time_point,sample[1,],type='b',pch=16, col='#372B2F', xlim=c(1,4),ylim=c(-15,20),xaxt='n',xlab='time points of non-manifesting carriers ',ylab='DNAm-age acceleration')
axis(1,c(1:4))
for (i in 2:115)
  {
    lines(time_point,sample[i,],type='b',pch=16,col='#372B2F')
   }

############################figS1B############################
four <- read.csv('fourtimepoint_manifesting.csv')
sample <- four[1:82,3:6]
time_point<-c(1,2,3,4)
plot(time_point,sample[1,],type='b',pch=16, col='#372B2F', xlim=c(1,4),ylim=c(-15,20),xaxt='n',xlab='time points of manifesting carriers ',ylab='DNAm-age acceleration')
axis(1,c(1:4))
for (i in 2:82)
  {
    lines(time_point,sample[i,],type='b',pch=16,col='#372B2F')
  }

############################figS2############################
cpg_1<-read.csv('./beta_baseline.csv') #please replace beta-value file u saved from minfi.R 
cpg_353<-read.csv('./353CpGlist.csv')
colnames(cpg_1)[1]<-'CpGSNP'
cpg_1_1<-semi_join(cpg_1,cpg_353,by="CpGSNP")
AO_AA1<-read.csv('AO-AA1.csv')
AO_AA1<-AO_AA1[,c(1,2,3)]
cpg<-cpg_1_1$CpGSNP
myCombat<-cpg_1_1[,-1]
row.names(myCombat)<-cpg
data_list<-apply(myCombat,1, list)
data_list<-lapply(data_list,cbind,AO_AA1)

for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[1] <- 'CpG'
}

fit <- lapply(data_list, lm, formula = CpG~Age_symptom_onset+GENDER)
fun2 <- function(x){
  y = data.frame(summary(x)[4])[2,4]
  return(y)}
p <- lapply(fit, fun2)

p_value <- data.frame(matrix(unlist(p), ncol = 1))
rownames(p_value) <- names(fit)
colnames(p_value)[1] <- 'pval'
p_value$CpG <- rownames(p_value)

p_sort <- sort(p_value$pval, index.return =T)
p_value_sort <- data.frame(p_sort = p_sort$x)
p_value_sort$cpg <- rownames(p_value)[p_sort$ix]
p_value_sort$bonferroni <- p.adjust(p_value_sort$p_sort, method = 'bonferroni')
p_value_sort$fdr <- p.adjust(p_value_sort$p_sort, method = 'fdr')
par(mfrow=c(2,3))
plot(data_list[["cg23124451"]]$Age_symptom_onset,data_list[["cg23124451"]]$CpG,xlab='AO',ylab='Beta',main='cg23124451',pch=20)+abline(fit[["cg23124451"]],col='blue')
plot(data_list[["cg06493994"]]$Age_symptom_onset,data_list[["cg06493994"]]$CpG,xlab='AO',ylab='Beta',main='cg06493994',pch=20)+abline(fit[["cg06493994"]],col='blue')
plot(data_list[["cg01820374"]]$Age_symptom_onset,data_list[["cg01820374"]]$CpG,xlab='AO',ylab='Beta',main='cg01820374',pch=20)+abline(fit[["cg01820374"]],col='blue')
plot(data_list[["cg09809672"]]$Age_symptom_onset,data_list[["cg09809672"]]$CpG,xlab='AO',ylab='Beta',main='cg09809672',pch=20)+abline(fit[["cg09809672"]],col='blue')
plot(data_list[["cg22736354"]]$Age_symptom_onset,data_list[["cg22736354"]]$CpG,xlab='AO',ylab='Beta',main='cg22736354',pch=20)+abline(fit[["cg22736354"]],col='blue')

############################figS3############################
cpg_1<-read.csv('./beta_3_year.csv')
cpg_353<-read.csv('./353CpGlist.csv')
colnames(cpg_1)[1]<-'CpGSNP'
cpg_1_1<-semi_join(cpg_1,cpg_353,by="CpGSNP")
AO_AA1<-read.csv('./AO-AA2.csv')
AO_AA1<-AO_AA1[,c(1,2,3)]
cpg<-cpg_1_1$CpGSNP
myCombat<-cpg_1_1[,-1]
row.names(myCombat)<-cpg
data_list<-apply(myCombat,1, list)
data_list<-lapply(data_list,cbind,AO_AA1)

for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[1] <- 'CpG'
}

fit <- lapply(data_list, lm, formula = CpG~Age_symptom_onset+GENDER)
fun2 <- function(x){
  y = data.frame(summary(x)[4])[2,4]
  return(y)}
p <- lapply(fit, fun2)

p_value <- data.frame(matrix(unlist(p), ncol = 1))
rownames(p_value) <- names(fit)
colnames(p_value)[1] <- 'pval'
p_value$CpG <- rownames(p_value)

p_sort <- sort(p_value$pval, index.return =T)
p_value_sort <- data.frame(p_sort = p_sort$x)
p_value_sort$cpg <- rownames(p_value)[p_sort$ix]
p_value_sort$bonferroni <- p.adjust(p_value_sort$p_sort, method = 'bonferroni')
p_value_sort$fdr <- p.adjust(p_value_sort$p_sort, method = 'fdr')
par(mfrow=c(2,3))
plot(data_list[["cg23124451"]]$Age_symptom_onset,data_list[["cg23124451"]]$CpG,xlab='AO',ylab='Beta',main='cg23124451',pch=20)+abline(fit[["cg23124451"]],col='blue')
plot(data_list[["cg06493994"]]$Age_symptom_onset,data_list[["cg06493994"]]$CpG,xlab='AO',ylab='Beta',main='cg06493994',pch=20)+abline(fit[["cg06493994"]],col='blue')
plot(data_list[["cg01820374"]]$Age_symptom_onset,data_list[["cg01820374"]]$CpG,xlab='AO',ylab='Beta',main='cg01820374',pch=20)+abline(fit[["cg01820374"]],col='blue')
plot(data_list[["cg09809672"]]$Age_symptom_onset,data_list[["cg09809672"]]$CpG,xlab='AO',ylab='Beta',main='cg09809672',pch=20)+abline(fit[["cg09809672"]],col='blue')
plot(data_list[["cg22736354"]]$Age_symptom_onset,data_list[["cg22736354"]]$CpG,xlab='AO',ylab='Beta',main='cg22736354',pch=20)+abline(fit[["cg22736354"]],col='blue')


############################figS4############################
cpg_1<-read.csv('./96ipd.csv')
cpg_353<-read.csv('./353CpGlist.csv')
colnames(cpg_1)[1]<-'CpGSNP'
cpg_1_1<-semi_join(cpg_1,cpg_353,by="CpGSNP")
AO_AA1<-read.csv('../DNAm_file/AO-AA-96ipd.csv')
AO_AA1<-AO_AA1[,c(1,6,8)]
cpg<-cpg_1_1$CpGSNP
myCombat<-cpg_1_1[,-1]
row.names(myCombat)<-cpg
data_list<-apply(myCombat,1, list)
data_list<-lapply(data_list,cbind,AO_AA1)

for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[1] <- 'CpG'
}

fit <- lapply(data_list, lm, formula = CpG~AO+sex)
fun2 <- function(x){
  y = data.frame(summary(x)[4])[2,4]
  return(y)}
p <- lapply(fit, fun2)

p_value <- data.frame(matrix(unlist(p), ncol = 1))
rownames(p_value) <- names(fit)
colnames(p_value)[1] <- 'pval'
p_value$CpG <- rownames(p_value)

p_sort <- sort(p_value$pval, index.return =T)
p_value_sort <- data.frame(p_sort = p_sort$x)
p_value_sort$cpg <- rownames(p_value)[p_sort$ix]
p_value_sort$bonferroni <- p.adjust(p_value_sort$p_sort, method = 'bonferroni')
p_value_sort$fdr <- p.adjust(p_value_sort$p_sort, method = 'fdr')
par(mfrow=c(2,3))
plot(data_list[["cg23124451"]]$AO,data_list[["cg23124451"]]$CpG,xlab='AO',ylab='Beta',main='cg23124451',pch=20)+abline(fit[["cg23124451"]],col='blue')
plot(data_list[["cg06493994"]]$AO,data_list[["cg06493994"]]$CpG,xlab='AO',ylab='Beta',main='cg06493994',pch=20)+abline(fit[["cg06493994"]],col='blue')
plot(data_list[["cg01820374"]]$AO,data_list[["cg01820374"]]$CpG,xlab='AO',ylab='Beta',main='cg01820374',pch=20)+abline(fit[["cg01820374"]],col='blue')
plot(data_list[["cg09809672"]]$AO,data_list[["cg09809672"]]$CpG,xlab='AO',ylab='Beta',main='cg09809672',pch=20)+abline(fit[["cg09809672"]],col='blue')
plot(data_list[["cg22736354"]]$AO,data_list[["cg22736354"]]$CpG,xlab='AO',ylab='Beta',main='cg22736354',pch=20)+abline(fit[["cg22736354"]],col='blue')


############################figS5############################
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain, features = c("Cbx7", "Scgn","Lag3","Edaradd","Nhlrc1"))
brain <- LoadData("stxBrain", type = "posterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain, features = c("Cbx7", "Scgn","Lag3","Edaradd","Nhlrc1"))
############################figS6############################
