library(survminer) 
library(survival)
data<-readRDS("/home/bio/Documents/glioma_HR/CGGA.Rds")

phe=data$pData
phe=phe[,c("survival","status")]
phe=na.omit(phe)
exp<-data$expr
exp<-exp[,9:ncol(exp)]
exp<-as.data.frame(t(exp))
exp<-exp[,rownames(phe)]
exp[1:3,1:3]
# colsd <- apply(exp, 1, sd)
# exp<-exp[colsd!=0,]

mySurv <- with(phe, Surv(survival, status))

gene_keep<-data.frame(gene=NULL,
                      id=NULL)

for (i in rownames(exp)) {
  gene= as.numeric(exp[i,])
  group=ifelse(gene>median(gene),'high','low')
  
  if( length(table(group))==2 & table(group)[1]>20){
    tmp<-data.frame(gene=i,id=2)
  }else{
    tmp<-data.frame(gene=i,id=1)
  }
  
  gene_keep<-rbind(gene_keep,tmp)
}

table(gene_keep$id)

exp<-exp[gene_keep$gene[which(gene_keep$id==2)],]

cox_results <-apply(exp, 1 , function(gene){
  
  group=ifelse(gene>median(gene),'high','low')
  
  survival_dat <- data.frame(group=group,# stage=phe$stage,
                             stringsAsFactors = F)
  m=coxph(mySurv ~ group,
          data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
})

dim(cox_results)
cox_results<-as.data.frame(t(cox_results))
cox_results$gene<-rownames(cox_results)

cox_results2<-cox_results[which(cox_results$HR<1 & cox_results$p<0.05),]
cox_results2$gene<-rownames(cox_results2)
# Rembrandt_hr<-cox_results
# save(Rembrandt_hr,file = "/home/bio/Documents/glioma_HR/Rembrandt_HR.Rda")
#cgga_hr<-cox_results2
#tcga_hr<-cox_results2
#Rembrandt_hr<-cox_results2
# survival_dat = phe
# gene = as.numeric(exp["IGF2BP3",])
# survival_dat$gene = ifelse(gene > median(gene),'high','low')
# table(survival_dat$gene)
# 
# fit <- survfit(Surv(survival, status) ~ gene,
#                data = survival_dat)
# 
# ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                    #legend.title = i, #定义图例的名称
#                    # legend = "top",#图例位置
#                    # legend.labs = c('High', 'Low'),
#                    pval = T, #在图上添加log rank检验的p值
#                    # pval.method = TRUE,#添加p值的检验方法
#                    risk.table = TRUE,
#                    risk.table.y.text = F,
#                    xlab = "Time in years", #x轴标题
#                    # xlim = c(0, 10), #展示x轴的范围
#                    # break.time.by = 1, #x轴间隔
#                    size = 1.5, #线条大小
#                    #ggtheme = theme_ggstatsplot(),
#                    palette="nejm", #配色
# )

gene<-intersect(rownames(cgga_hr),rownames(tcga_hr))
gene<-intersect(gene,rownames(Rembrandt_hr))

# 1013+667+286


library(tidyverse)
Rembrandt<-data.table::fread("/home/bio/Documents/glioma_HR/Rembrandt_mRNA_array_475.txt",header = T,data.table = F)
Rembrandt[1:3,1:3]
Rembrandt<-Rembrandt %>% column_to_rownames(var = "gene")

clin=data.table::fread("/home/bio/Documents/glioma_HR/Rembrandt_mRNA_array_475_clinical.txt",header = T,data.table = F)
clin=as.data.frame(clin)
rownames(clin)<-clin$gene
library(survival)
library(survminer)
df<-data.frame(row.names = rownames(clin),
               Histology=clin$Histology)
df$survival<-as.numeric(clin[rownames(df),"OS"])/30
df$status<-clin[rownames(df),"Censor"]
df=df[which(df$Histology!="NON_TUMOR"),]
df=na.omit(df)
Rembrandt=Rembrandt[,rownames(df)]

phe=df
exp=Rembrandt

deg<-data.table::fread("/home/bio/Documents/glioma_HR/Glioma_GSE131928_Smartseq2_AllDiffGenes_table.tsv",header = T,data.table = F)

deg2<-deg[deg$log2FC>1,]

table(deg2$`Celltype (malignancy)`)

gene2=deg2$Gene[which(deg2$`Celltype (malignancy)`=="Malignant cells")]

gliomaHR=intersect(gene,gene2)

save(gliomaHR,file = "/home/bio/Documents/glioma_HR/gliomaHR.Rda")


marker<-data.table::fread("Glioma_GSE131928_Smartseq2_AllDiffGenes_table (1).tsv",header = T)
marker_cancer<-marker[which(marker$`Celltype (malignancy)`=="Malignant cells"),]
marker_cancer<-marker_cancer[marker_cancer$log2FC>=1 & marker_cancer$`Adjusted p-value`<0.05,]

gene_glioma=intersect(marker_cancer$Gene,gene)
geneAll<-gene

df<-data.frame(geneAll)
rownames(df)<-geneAll
df$group<-1
df$group<-ifelse(df$geneAll %in% gene_glioma,1,2)

TMEgene<-rownames(df)[df$group==2]

save(gene_glioma,geneAll,TMEgene,file="glioma_Sig.Rda")



