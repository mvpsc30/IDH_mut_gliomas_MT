load("scemk.Rda")
scemk2<-scemk[which(scemk$cluster=="positive"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

diff<-read.csv("IDHi_treated.csv",header = T)
library(clusterProfiler)
library(org.Hs.eg.db)

diff<-diff[order(diff$log2.ratio,decreasing = T),]
gene_list<-diff$log2.ratio
names(gene_list)<-diff$Gene
gene_list
marker<-readRDS("Neftel_et_al_2019_four_state_neoplastic_markers.rds")
marker<-marker$NPC
load("glioma_Sig.Rda")
#marker<-sc.marker$gene[which(sc.marker$cluster=="MESlike" & sc.marker$avg_log2FC>=1)]

geneset<-data.frame(term="SenTAM",gene=sc.marker$gene)

y <- GSEA(gene_list,TERM2GENE = geneset,pvalueCutoff = 0.99)

yd<-data.frame(y)

library(GseaVis)
library(ggplot2)
gseaNb(object = y,
       geneSetID = 'SenTAM',
       subPlot = 2,
       addPval= T,
       pvalX= 0.95,
       pvalY= 0.8,
       # addGene = T,
       # markTopgene = T,
       # topGeneN=15,
       curveCol= c("#C3DD9B","#E63D70"),
       htCol= c("#C3DD9B","#E63D70"),
       rankCol= c("#C3DD9B","white","#E63D70"))

#ggsave("IDH_treated_senTAM.pdf",height = 4,width = 6)
#ggsave("IDHi_NPC.pdf",height=4,width = 6)
files <- list.files("/Users/mac/Documents/scRNA/glioma/glass/cell_sensence/")

senSet<-list()

for (i in 1:15) {
  gst<-read.gmt(paste0("cell_sensence/",files[i]))
  senSet[[i]]<-gst$gene
}

names(senSet)<-gsub(files, pattern = ".v2023.2.Hs.gmt", replacement = "")

SenMayo<-openxlsx::read.xlsx("/Users/mac/Documents/scRNA/glioma/glass/41467_2022_32552_MOESM4_ESM (1).xlsx")

senSet$SenMayo<-SenMayo$`Gene(human)`

signature_immport_df<-stack(senSet)
colnames(signature_immport_df)<-c("gene","term")
signature_immport_df<-signature_immport_df[,c("term","gene")]

alldiff <- diff[order(diff$log2.ratio,decreasing = T),]

id <- alldiff$log2.ratio
names(id) <- alldiff$Gene
id

y <- GSEA(id,
          TERM2GENE = signature_immport_df,
          pvalueCutoff = 0.99,
)
yd<-data.frame(y)

setid<-rownames(yd)

gseaNb(object = y,
       geneSetID = setid[1],
       subPlot = 2,
       addPval= T,
       pvalX= 0.95,
       pvalY= 0.8,
       addGene = T,
       markTopgene = T,
       topGeneN=15,
       curveCol= c("#C3DD9B","#E63D70"),
       htCol= c("#C3DD9B","#E63D70"),
       rankCol= c("#C3DD9B","white","#E63D70"))


library(GseaVis)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(patchwork)
library(ggpubr)
library(fgsea)
gs_list<-readRDS("/Users/mac/Downloads/difg_glass_2019/gene_set_list.Rds")
gs_list<-gs_list$LM22
library(clusterProfiler)
signature_immport_df<-stack(gs_list)
colnames(signature_immport_df)<-c("gene","term")
signature_immport_df<-signature_immport_df[,c("term","gene")]
################################################################################
load("marker_scRNA.Rda")
sce.marker$pt<-sce.marker$pct.1-sce.marker$pct.2

marker<-sce.marker[which(sce.marker$avg_log2FC>=1 & sce.marker$pt>0.2),]
table(marker$cluster)

signature_immport_df<-data.frame(term=marker$cluster,gene=marker$gene)


alldiff <- diff
alldiff$log2.ratio<- -(alldiff$log2.ratio)
alldiff <- alldiff[order(alldiff$log2.ratio,decreasing = T),]

id <- alldiff$log2.ratio
names(id) <- alldiff$Gene
id

y <- GSEA(id,
          TERM2GENE = signature_immport_df,
          pvalueCutoff = 1,
)
yd<-data.frame(y)

setid<-c("Diff-like","Stem-like","Prolife.stem-like",
         "Myeloid","Oligodendrocyte","Endothelial",
         "Fibroblast","T cell")
#setid<-setid[c(2,3,5)]
gseaNb(object = y,
       geneSetID = setid[8],
       subPlot = 2,
       addPval= T,
       pvalX= 0.95,
       pvalY= 0.8,
       topGeneN=15,
       curveCol= c("#C3DD9B","#E63D70"),
       htCol= c("#C3DD9B","#E63D70"),
       rankCol= c("#C3DD9B","white","#E63D70"))


# my_theme<-theme(panel.background = element_rect(fill = "black", colour = "white",linewidth=1),
#                 plot.background = element_rect(fill = "black", colour = "white"),
#                 panel.grid.major = element_line(color = "black"),
#                 panel.grid.minor = element_line(color = "black"),
#                 text = element_text(color = "white"),
#                 axis.text.y = element_text(color = "white"),
#                 axis.title = element_text(color = "white"))

pdf("IDHtreated_gsea.pdf",height = 15,width = 7)

gseaNb(object = y,
       geneSetID = setid,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX= 0.5,
       pvalY = 0.2,
       pCol="white",
       termWidth=15)+black_theme
dev.off()