#This code used ttScreening method in a two step screening process to screening high dimensional epigenetic biomarkers
#IOW data used for demonstration
#after performing step 1 and step 2 of ttscreening, return the ttscreening passing rate for identified DNAm
#outcome variables are ASCON 1 18months 2 4 10 18 (children' asthma observed at 6 consecutive measures)

library("MASS")
library(dplyr)
library(lme4)
library(lmerTest)
library(geepack)
library(reshape2)
library(data.table)
library("haven")
setwd("C:/Users/lxie3/OneDrive - The University of Memphis/Data_Mediation/OneDrive_1_10-13-2020")
setwd("C:\\Users\\Xie\\OneDrive - The University of Memphis\\Data_Mediation\\OneDrive_1_10-13-2020")
load("DNAm Guthrie 12112019.RData")
#load("/home/lxie3/DNAm Guthrie 12112019.RData")



##Read in phenotype data
PhenoType=as.matrix(read_sas("C:\\Users\\lxie3\\OneDrive\\mediation\\Data\\OneDrive_1_3-17-2020\\iow_v36a.sas7bdat"))
#PhenoType=as.matrix(read_sas("/home/lxie3/iow_v36a.sas7bdat"))
dim(PhenoType)
#[1] 1536 2243
Pool=paste("X_",as.numeric(PhenoType[,"STUDYid"]),sep="")
rownames(PhenoType)=Pool


#prepare the exposure and outcome variables
PhenoType1=PhenoType[,c("MSMK_03","ASCON_1", "AST_2","ASCON_2", "ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10", "M_INVESTIGATORDIAGNOSEDASTHMA_18",   "SEX_18","ASTHMY1_0","GESTAGE")] #dim 1536 6
PhenoType1_complete <- PhenoType1[complete.cases(PhenoType1), ]
dim(PhenoType1_complete)
#[1] 924 10


##Read in the DNAm data
dim(DNAm)
#[1] 551710    797
#Note that for DNAm of size 551710, please submit job to HPC for computation
#For demonstration only, here I choose the first 500 DNAm as DNAm_sample
DNAm_ID=colnames(DNAm)
sample_ID=intersect(DNAm_ID,rownames(PhenoType1_complete))
length(sample_ID)
#591
rownames(DNAm)=DNAm[,1]
DNAm_sample=DNAm[1:500,sample_ID]

#dim 500 591
PhenoType.data=PhenoType1[sample_ID,c("ASCON_1", "AST_2","ASCON_2", "ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10","M_INVESTIGATORDIAGNOSEDASTHMA_18")]
#dim 591 6 as repeated measures
table(rownames(PhenoType.data)==colnames(DNAm_sample))
#TRUE
# 591



############################################################################################################################################
##Start two-step ttScreening process########################################################################################################
##Users need to prepare two datasets, PhenoType.data, DNAm.data#############################################################################
############################################################################################################################################


DNAm.data<-DNAm_sample

### P value using geeglm for binary repeated outcomes
### change family, corstr based on real data structure

geeglm_p=function(datanew){
  datanew2=as.data.frame(datanew)
  long <- melt(setDT(datanew2), id.vars = c("ID","M"),
               measure.vars=c("Y1","Y2","Y3","Y4", "Y5","Y6"),
               variable.name = "Y",variable.factor=TRUE)
 pvalue =coef(summary(geeglm(value~ M + Y,id=ID,
                              family = "binomial",corstr = "AR-1",
                              data =  long)))[2, c("Pr(>|W|)")]
  return(pvalue)
}


p=dim(DNAm.data)[1]#number of DNAm for screening
n=dim(DNAm.data)[2]#sample size
class(PhenoType.data)="numeric"
pvalue.all=rep(NA, p)


### step 1 :screening each cg
for (i in 1:p){
  ID=1:n
  tempdata <- cbind(ID, Y=PhenoType.data[,c("ASCON_1","AST_2","ASCON_2" ,"ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10","M_INVESTIGATORDIAGNOSEDASTHMA_18")], M=t(DNAm.data[i,]))
  colnames(tempdata)=c("ID","Y1","Y2","Y3","Y4","Y5","Y6","M")
  pvalue.all[i]=geeglm_p(tempdata)
}


#step 2: use 0.1 as cutoff, use ttscreening for those DNAm's which has p value less than 0.1 in the first step
Final=matrix(NA,nrow=p,ncol=2)
colnames(Final)<-c("pass","tt_rate")
rownames(Final)<-rownames(DNAm.data)

s=100
for (i in 1:p){
  if (pvalue.all[i] >=0.1) {
    Final[i]=0
  }else{
    ID=1:n
    temp<- cbind(ID,Y=PhenoType.data[,c("ASCON_1","AST_2","ASCON_2" ,"ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10","M_INVESTIGATORDIAGNOSEDASTHMA_18")], M=t(DNAm.data[i,]))
    colnames(temp)=c("ID","Y1","Y2","Y3","Y4","Y5","Y6","M")
  #start ttScreening for 100 times of iteration
    sig_path=rep(NA,s)
    for (k in 1:s) {
      set.seed(k)

      sample <- sample.int(n = n, size = floor(.67*n), replace = F)## Every loop generate IDs for train and test datasets

      ### we need to try each
      train <- temp[sample, ]
      test  <- temp[-sample, ]
      Path.out.train=geeglm_p(train)
      Path.out.test=geeglm_p(test)

      sig_path[k] <- ifelse (Path.out.train<0.05 &Path.out.test<0.1,1,0)
    }


    Final[i,1]=ifelse (sum(sig_path)>10,1,0) #Passing ttScreening for 10 times, change this number to 50 or higher, for high dimensional data
    Final[i,2]=sum(sig_path)/s
  }

}


#Final[,1]==1 will give selected DNAm from the two-step screening process

#List the selected DNAm and their passing rate
selected<-subset(Final,Final[,1]==1)
No.selected<-dim(selected)[1] #number of identified DNAm after ttscreening

selected_ID=rownames(selected)

DNAm_selected<-DNAm.data[selected_ID,]

selected


#CpG site names and passing rate will be listed here
#in the end we want the function to be in this way, repeated outcome should be in terms of the wide format
ttscreening_repeated(RepeatedOutcome,m=6,DNAm,family = "binomial",corstr = "AR-1", Cutoff.Joint=0.1,Iterations=100,Train.SigLevel=0.05,Test.SigLevel=0.1,Percentage=0.1)











