load_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
  if(length(new_pkgs)) install.packages(new_pkgs)
  for(pkg in pkgs){
    suppressWarnings(suppressMessages(library(pkg, character.only = T)))
  }
}
pkgs <- c('dplyr', 'minfi', 'ggplot2', 'survival','survminer',
          'coxme','Seurat','SeuratData','patchwork')
load_pkgs(pkgs)

###############################################################################
#codes for preprocessing raw data and getting the file of DNAm beta value 
###############################################################################
#read idat files according to the sample sheet 
targets<-read.metharray.sheet("./",pattern = "baseline_annotable.csv")
rgSet<-read.metharray.exp(targets = targets,force=TRUE)

#calculate probes' p-value and the probes with p-value > 0.01 are filtered out
probep<-detectionP(rgSet)
keep<-apply(probep,1,function(t){all(t<0.01)})
rgSet<-rgSet[keep,]

#probes with mean p-value > 0.05 are filtered out
keep<-apply(probep, 2, mean) < 0.05   
rgSet<-rgSet[,keep]

#normalize p-value
GRSet.Quantile <- preprocessQuantile(rgSet) 

#save the file with beta values
beta.Quantile <- getBeta(GRSet.Quantile)
write.csv(beta.Quantile,file = './beta_baseline.csv')

# produce a file with beta values that contains all probes 
#in datMiniAnnotation3.csv which the DNAm clock caculator website provides
annote<-read.csv('datMiniAnnotation3.csv')
input<-as.data.frame(beta.Quantile)
input$ProbeID<-rownames(input)
output<-left_join(annote,input,by='ProbeID')
outputs <- output[,-(2:7)]
write.csv(outputs,file = 'baseline_output.csv',row.names = FALSE) 
#the file will be uploaded to the DNAm clock caculator website (http://dnamage.genetics.ucla.edu)

###############################################################################
#codes used to generate figures in the main text
###############################################################################
#AO=age of symptom onset,AOA=age of assessment,AA= DNAm age acceleration, Status=1:manifesting

#fig1-2:cox proportional hazard regression analysis and Kaplan-Meier curves#
AOA_AA_blood <- read.csv('AOA-AO.csv') #the file with AO,AA,AOA and other information for analysis
AOA_AA.cox <- coxph(Surv(AOA, Status1) ~ AA +gender+relatedness +CD8T + CD4T +Bcell +Gran, data = AOA_AA_blood)
AOA_AA_with <- with(AOA_AA_blood, Surv(AOA, Status1 == 1))
AOA_AA_sf <- survfit(AOA_AA_with ~ AA, data = AOA_AA_blood, conf.type = "log-log")
AOA_AA_pic<-ggsurvplot( AOA_AA_sf, data=AOA_AA_blood, conf.int=FALSE, legend.labs=c("Slow (AA<-3)","Normal (-3<AA<3)","Fast (AA>3)"),xlab="Age of symptom onset", ylab="Cumulative incidence", legend.title=c("DNAm-age acceleration"), palette = c("#009E73", "#56B4E9","#D55E00"), risk.table=TRUE,legend=c(0.2,0.9), fun = function(x) {1 - x}, break.x.by=10, xlim = c(30,80))


#fig3:scatter plot and linear regression analysis of AO and AA#
AO_AA_blood <- read.csv('AO-AA.csv') #file with AO, AA and other information for analysis
AO_AA <- lm(AO ~ AA+gender+interval,data = AO_AA_blood)
AO_AA_pic1 <- ggplot(data=AO_AA_blood, aes(x=AA, y=AO)) +
  geom_point()+labs(x='DNAm Age acceleration at baseline sample collection',y='Age of onset')
AO_AA_pic2<-geom_abline(intercept = coef(AO_AA)[1],slope = coef(AO_AA)[2], color = "blue") 
AO_AA_pic<-AO_AA_pic1+AO_AA_pic2+theme_classic()+scale_x_continuous(breaks=seq(-20,30,5))+scale_y_continuous(breaks=seq(40,100,10))

###############################################################################
#codes used to generate figures in the supplementary file
###############################################################################
#figS1:AA changes at four time-points#
four <- read.csv('fourtimepoint.csv') #the file with AA values at four time-points
sample <- four[1:115,3:6] #115 samples' AA values at four time-points 
time_point<-c(1,2,3,4)
plot(time_point,sample[1,],type='b',pch=16, col='#372B2F', xlim=c(1,4),ylim=c(-15,20),xaxt='n',xlab='time-points of non-manifesting carriers ',ylab='DNAm-age acceleration')
axis(1,c(1:4))
for (i in 2:115)
{
  lines(time_point,sample[i,],type='b',pch=16,col='#372B2F')
}

#figS2-S5: scatter plot and linear regression analysis of 353 CpG sites and AO#
cpg_1<-read.csv('./beta_value.csv') #the beta value file which is saved from the preprocessing raw data step 
cpg_353<-read.csv('./353CpGlist.csv')  #353 age-related CpGs file
colnames(cpg_1)[1]<-'CpGSNP'
cpg_1_1<-semi_join(cpg_1,cpg_353,by="CpGSNP")
AO_AA<-read.csv('AO-AA.csv') #the file with patientID, age of symptom onset and gender
AO_AA<-AO_AA[,c(1,2,3)]
cpg<-cpg_1_1$CpGSNP
myCombat<-cpg_1_1[,-1]
row.names(myCombat)<-cpg
data_list<-apply(myCombat,1, list)
data_list<-lapply(data_list,cbind,AO_AA) #a list with beta value, probeID, patientID, age of symptom onset and gender

for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[1] <- 'CpG'
}

fit <- lapply(data_list, lm, formula = CpG~AO+gender) #linear regression analysis of 353 CpGs and AO
fun2 <- function(x){
  y = data.frame(summary(x)[4])[2,4]
  return(y)}
p <- lapply(fit, fun2) #get p-value

p_value <- data.frame(matrix(unlist(p), ncol = 1))
rownames(p_value) <- names(fit)
colnames(p_value)[1] <- 'pval'
p_value$CpG <- rownames(p_value) 

p_sort <- sort(p_value$pval, index.return =T) 
p_value_sort <- data.frame(p_sort = p_sort$x)
p_value_sort$cpg <- rownames(p_value)[p_sort$ix]
p_value_sort$bonferroni <- p.adjust(p_value_sort$p_sort, method = 'bonferroni') #use bonferroni and fdr to adjust p-value
p_value_sort$fdr <- p.adjust(p_value_sort$p_sort, method = 'fdr')
write.csv(p_value_sort,'p_value.csv',row.names = FALSE) #save the file with CpGs and p-values

par(mfrow=c(2,3))
#plot 5 shared age-related CpGs in LRRK2 carriers, idiopathic PD patients and sporadic PD patients
plot(data_list[["cg23124451"]]$AO,data_list[["cg23124451"]]$CpG,xlab='AO',ylab='Beta',main='cg23124451',pch=20)+abline(fit[["cg23124451"]],col='blue')
plot(data_list[["cg06493994"]]$AO,data_list[["cg06493994"]]$CpG,xlab='AO',ylab='Beta',main='cg06493994',pch=20)+abline(fit[["cg06493994"]],col='blue')
plot(data_list[["cg01820374"]]$AO,data_list[["cg01820374"]]$CpG,xlab='AO',ylab='Beta',main='cg01820374',pch=20)+abline(fit[["cg01820374"]],col='blue')
plot(data_list[["cg09809672"]]$AO,data_list[["cg09809672"]]$CpG,xlab='AO',ylab='Beta',main='cg09809672',pch=20)+abline(fit[["cg09809672"]],col='blue')
plot(data_list[["cg22736354"]]$AO,data_list[["cg22736354"]]$CpG,xlab='AO',ylab='Beta',main='cg22736354',pch=20)+abline(fit[["cg22736354"]],col='blue')


#figS6:spatial gene expression in the mouse brain#
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain, features = c("Cbx7", "Scgn","Lag3","Edaradd","Nhlrc1"))
brain <- LoadData("stxBrain", type = "posterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain, features = c("Cbx7", "Scgn","Lag3","Edaradd","Nhlrc1"))

#figS7:scatter plot and linear regression analysis of AA and DaTScan striatal binding ratio#
table <- read.csv("./LRRK2_Methylation_DAT_Full_clinical.csv") #the file with AA, gender,DaTScan striatal binding ratio and other information for analysis
colnames(table)[3] <- "Diagnosis"

a <- filter(table,Diagnosis == "1") #only keep manifesting LRRK2 carriers
a <- filter(a,AA != "#VALUE!")
a <- filter(a,Left_caudate_at_baseline != "NA")
a$Left_caudate_at_baseline <- as.numeric(a$Left_caudate_at_baseline)#replace it with Right_caudate_at_baseline, Left_putamen_at_baseline and Right_putamen_at_baseline to plot other 3 figures
a$AA<-as.numeric(a$AA)
a$Disease_duration<-as.numeric(a$Disease_duration)
cau_l <- lm(Left_caudate_at_baseline ~ AA + gender + family_info + Disease_duration, data = a)
plot<- ggplot(data=a,aes(x=AA,y= Left_caudate_at_baseline))+geom_point() + 
  geom_abline(intercept = coef(cau_l)[1],slope = coef(cau_l)[2],color="blue",size=1) + 
  scale_y_continuous(name = "DatScan striatal binding ratio of Left caudate",limits=c(0,4)) + 
  theme_bw() + theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line=element_line(size = 1,colour = 'black'))