library(Seurat)
library(clusterProfiler)
# install the CARD package
library(CARD)
library(Seurat)
library(nichenetr)
library(readxl)
library(dplyr)
library(readr)
library(nichenetr)
library(dplyr)
library(hdf5r)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

load("scemk.Rda")
scemk2<-scemk[which(scemk$cluster=="positive"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

scemk2<-scemk[which(scemk$cluster=="Bg"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker2<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

senescence_marker= openxlsx::read.xlsx("41467_2022_32552_MOESM4_ESM (1).xlsx",sheet = 1)
senescence_marker=senescence_marker$`Gene(human)`

list<-list(gs=sc.marker$gene,
           bg=sc.marker2$gene,
           sen=senescence_marker)

pt<-Load10X_Spatial(data.dir = "/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM4_spaceranger_out/",
                    filename = "filtered_feature_bc_matrix.h5",
                    slice ="gbm")

pt <- SCTransform(pt, assay = "Spatial", verbose = FALSE)
## Assign Ivy GAP defined niches
DefaultAssay(pt) <- "SCT"

pt <-AddModuleScore(pt,features = list,
                    ctrl = min(c(vapply(X = list, 
                                        FUN = length, 
                                        FUN.VALUE = numeric(length = 1))), 
                               200),
                    nbin=10)

colnames(pt@meta.data)[6:8]<-c("pos","bg","sen")

median(pt$pos)
summary(pt$pos)

pt$region<-ifelse(pt$pos>pt$bg,"pos","bg")

pt$region2<-ifelse(pt$sen>0.2,"sen","bg")
table(pt$region)
table(pt$region2)

p1=SpatialDimPlot(pt,group.by = "region",  stroke = 0,   image.alpha = 0, alpha = 1)+
  scale_fill_manual(values = c("lightblue","darkolivegreen3"))

p2=SpatialDimPlot(pt,group.by = "region2",  stroke = 0,   image.alpha = 0, alpha = 1) + 
  scale_fill_manual(values = c( "lightblue","violetred1"))



p3=SpatialDimPlot(pt,group.by = "region2",  stroke = 0,   image.alpha = 1, alpha =0) + 
  scale_fill_manual(values = c( "darkolivegreen3","violetred1"))

p1/p2/p3

ggsave("spatial.pdf",height = 10,width = 4)

theme(panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"))



pt<-Load10X_Spatial(data.dir = "/Users/mac/Downloads/GSE237183_RAW/bwh35/",
                    filename = "filtered_feature_bc_matrix.h5",
                    slice ="gbm")

pt <- SCTransform(pt, assay = "Spatial", verbose = FALSE)

pt <-AddModuleScore(pt,features = list,
                    ctrl = min(c(vapply(X = list, 
                                        FUN = length, 
                                        FUN.VALUE = numeric(length = 1))), 
                               100))

colnames(pt@meta.data)[6:8]<-c("pos","bg","sen")

pt$region<-ifelse(pt$pos>pt$bg,"pos","bg")
pt$region<-ifelse(pt$pos==0,NULL,pt$region)

pt$region2<-ifelse(pt$sen>0.2,"sen","bg")
table(pt$region)
table(pt$region2)

SpatialDimPlot(pt,group.by = "region",  stroke = 0,   image.alpha = 0, alpha = 1)+
  SpatialDimPlot(pt,group.by = "region2",  stroke = 0,   image.alpha = 0, alpha = 1) + 
  scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) +
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"))


source("/Users/mac/Documents/scRNA/GBM/cell_2024_colocal/Kloosterman_et_al_Cell_2024-main/unique functions/assignLocation.R")
IvyGAP <- read_excel("/Users/mac/Documents/scRNA/GBM/cell_2024_colocal/Kloosterman_et_al_Cell_2024-main/signatures/IVY_gap_signatures.xlsx")
colnames(IvyGAP)[1]<-"Gene"
ct.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "CT")
pan.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "CTpan")
le.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "LE")
mvp.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "CTmvp")

pt <- assignLocation(pt, 
                     ct.features = ct.all, 
                     pan.features = pan.all, 
                     le.features = le.all, 
                     mvp.features = mvp.all)
SpatialDimPlot(pt, 
               group.by = "Location" ,  stroke = 0, image.alpha = 0, alpha = 1) +
  scale_fill_manual(values = c("#C1DB82", "#7776AD", "#D58B3F","#3070B0", "grey")) 


load("mye.Rda")
DimPlot(mye)
mye$celltype<-as.character(mye$cell) %>% {
  .[.=='P2RY12+(Microglia)'] <- "Microglia"
  .[.=='CD87+'] <- "TAM3"
  .[.=='GPM6B+'] <- 'TAM2'
  .[.=='FTL+'] <- "TAM1"
  ;.}

Idents(mye)<-"celltype"
DimPlot(mye)
marker<-FindAllMarkers(mye,only.pos = T,logfc.threshold = 1)

mye$Rec<-as.character(mye$id) %>% {
  .[.%in% c("RV18002","RV19005","RV18003")] <- "Recurrence"
  .[.%in% c("RV18001","RV19007","RV19001")] <- "Initial"
}

meta<-as.data.frame(mye@meta.data)

meta$Rec<-as.character(meta$id) %>% {
  .[.%in% c("RV18002","RV19005","RV18003")] <- "Recurrence"
  .[.%in% c("RV18001","RV19007","RV19001")] <- "Initial"
}
table(mye$sampleid)
table(mye$Rec)

mye$Rec<-ifelse(mye$id=="RV19005","Recurrence",mye$Rec)


p1=scRNAtoolVis::cellRatioPlot(mye,sample.name = "Rec",celltype.name = "celltype",fill.col = c("#E1A39C","#9C623D","#4C74A7"))
library(Seurat)


p2=VlnPlot(mye,group.by = "Rec",features = "Fridman senescence",cols = c("#C3DD9B","#C58177"),pt.size = 0)+
  xlab("")+NoLegend()+ylab("Senescence scores")+
  stat_compare_means()
p1+p2+patchwork::plot_layout(widths = c(1,1.5))
ggsave("Tam_propor.pdf",height = 6,width = 8)

DimPlot(mye,split.by = "Rec",group.by = "celltype")

tam1 <-  subset(marker$gene, subset = marker$cluster == "TAM2 & TAM1")
tam2 <-  subset(marker$gene, subset = marker$cluster == "TAM2 & TAM1")
tam3 <-  subset(marker$gene, subset = marker$cluster == "TAM3")

pt <- assignLocation2(pt, 
                      tam1 = tam1, 
                      tam2 = tam2, 
                      tam3 = tam3)

SpatialDimPlot(pt, 
               group.by = "Location2" ,  stroke = 0, image.alpha = 0, alpha = 1) +
  scale_fill_manual(values = c("#C1DB82", "#7776AD", "#D58B3F","#3070B0", "grey")) 


#### figure 3
library(openxlsx)
data<-read.xlsx("GSE153508_EtOH_vs_4OHT_DESEQ.xlsx")
data<-na.omit(data)

library(clusterProfiler)
library(org.Mm.eg.db)
colnames(data)
load("scemk.Rda")
scemk2<-scemk[which(scemk$cluster=="positive"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

load("/Users/mac/Documents/scRNA/glioma/glass/mye_back.Rda")
library(Seurat)
marker_mye<-FindAllMarkers(mye,only.pos = T,logfc.threshold = 0.25)
marker_mye_gene<-data.frame(gene=marker_mye$gene,term=marker_mye$cluster)

library(nichenetr)
sc.marker<-marker_mye_gene[which(marker_mye_gene$term=="TAM3"),]
gs<-convert_human_to_mouse_symbols(sc.marker$gene)
gs<-as.character(gs)
gs<-na.omit(gs)


signature_immport_df<-data.frame(term="sentam",gene=gs)


alldiff <- data
alldiff <- alldiff[order(alldiff$`log2FoldChange.(4-OHT/Control)`,decreasing = T),]

id <- alldiff$`log2FoldChange.(4-OHT/Control)`
names(id) <- alldiff$Gene
id

y <- GSEA(id,
          TERM2GENE = signature_immport_df,
          pvalueCutoff = 2
)
yd<-data.frame(y)

library(GseaVis)
gseaNb(object = y,
       geneSetID ="sentam",
       addPval = T,
       rmHt = T,
       pvalX= 0.95,
       pvalY = 0.7,subPlot=2,
       termWidth=15)
dev.off()
