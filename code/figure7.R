library(openxlsx)
library(survival)
library(survminer)
library(meta)
load("glioma_Sig.Rda")

tcga<-readRDS("/Users/mac/Downloads/shiny_GlioVis-master/data/datasets/Rembrandt.Rds")
clin=tcga$pData
clin=as.data.frame(clin)
#unique(clin$Histology)
#clin=clin[which(clin$Histology!="GBM"),]
table(clin$Histology)
clin=clin[which(clin$Histology %in% c("Oligodendroglioma",
                                      "Astrocytoma")),]
table(clin$Recurrence)

clin=clin[which(clin$Recurrence!="Primary"),]

rownames(clin)=clin$Sample

colnames(clin)
exp<-tcga$expr
exp<-exp[,9:ncol(exp)]
exp<-as.data.frame(t(exp))
exp[1:3,1:20]
library(survival)
library(survminer)
library(Seurat)
exp<-exp[,rownames(clin)]

seno<-read.xlsx("41467_2022_32552_MOESM4_ESM (1).xlsx",sheet = 1)

gsva_ex<-GSVA::gsva(expr=as.matrix(exp), 
                    gset.idx.list = list(seno$`Gene(human)`), 
                    method="ssgsea",verbose= FALSE)

clin$sig<-as.numeric(gsva_ex)
clin$sig<-as.numeric(exp["ITGAV",])
res.cut <- surv_cutpoint(clin, #数据集
                         time = "survival", #生存状态
                         event = "status", #生存时间
                         variables = c("sig") #需要计算的数据列名
)

plot(res.cut, "sig", palette = "npg")

#clin$cluster <-ifelse(clin$sig>median(clin$sig),"High","Low")

clin$cluster <-ifelse(clin$sig>res.cut$cutpoint[1,1],"High","Low")

fit <- survfit(Surv(survival,status) ~ cluster,  data = clin) 

ggsurvplot(fit, # 创建的拟合对象
           data = clin,  # 指定变量数据来源
           conf.int = F, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = TRUE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           palette = "nejm")  # 自定义调色板

table(clin$cluster,clin$status)


df<-openxlsx::read.xlsx("meta.xlsx")

m1 <- metabin(e1, n1, e2, n2,
              data = df,
              studlab = study,
              sm="OR",
              comb.random=FALSE)
m1

pdf("forest.pdf",height = 5,width = 12)
forest(m1)
dev.off()