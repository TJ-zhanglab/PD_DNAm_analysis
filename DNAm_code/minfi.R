library(minfi)
library(dplyr)
targets<-read.metharray.sheet("./",pattern = "baseline_annotable.csv")
rgSet<-read.metharray.exp(targets = targets,force=TRUE)
pd <- pData(rgSet)

probep<-detectionP(rgSet)
keep<-apply(probep,1,function(t){all(t<0.01)})
rgSet<-rgSet[keep,]
keep<-apply(probep, 2, mean) < 0.05
rgSet<-rgSet[,keep]

GRSet.Quantile <- preprocessQuantile(rgSet) 

annotation <- getAnnotation(GRSet.Quantile)
sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX","chrY")]
keep <- !(featureNames(GRSet.Quantile) %in% sex_probe)
GRSet.Quantile <- GRSet.Quantile[keep,]

beta.Quantile <- getBeta(GRSet.Quantile)
write.csv(beta.Quantile,file = './beta_baseline.csv')