library(Seurat)
library(harmony)
library(clusterProfiler)
library(SCP)
library(ggplot2)
library(ggrastr)
library(patchwork)

gene_sets<-list(AllSig=geneAll,
                MTSig=gene_glioma,
                TMESig=TMEgene)
sce<-AddModuleScore(sce,gene_sets)
colnames(sce@meta.data)[29:31]<-c("AllSig","MTSig","TMESig")

# FeaturePlot(sce,features = c("AllSig","MTSig","TMESig"),
#             reduction = "tsne",
#             cols = c("blue","yellow","red"))+
#   plot_layout(nrow = 1,ncol = 3)

p1=FeatureDimPlot(sce,features = c("AllSig"),
                  reduction = "UMAP",bg_cutoff = -1,
                  show_stat = F,theme_use = "theme_classic")
p1

p2=FeatureDimPlot(sce,features = c("MTSig"),
                  lower_cutoff = NULL,
                  upper_cutoff = NULL,
                  reduction = "UMAP",bg_cutoff = -1,
                  show_stat = F,theme_use = "theme_classic")
p2

p3=FeatureDimPlot(sce,features = c("TMESig"),
                  reduction = "UMAP",bg_cutoff = -1,
                  show_stat = F,theme_use = "theme_classic")
p3

pdf("signature1.pdf",height = 3,width = 3.4)
rasterise(p1,dpi = 400)
dev.off()

pdf("signature2.pdf",height = 3,width = 3.4)
rasterise(p2,dpi = 400)
dev.off()

pdf("signature3.pdf",height = 3,width = 3.4)
rasterise(p3,dpi = 400)
dev.off()


mye<-subset(sce,name %in% rownames(cellmeta))
mye$cell<-cellmeta$cell

mye <- mye %>% 
  NormalizeData(., 
                normalization.method = "LogNormalize", 
                scale.factor = 10000) %>%
  FindVariableFeatures(., selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(., features = rownames(.)) %>%
  RunPCA(.,  npcs = 25) %>% 
  RunHarmony(.,reduction = "pca",
             dims.use = 1:25,
             group.by.vars = "id",
             plot_convergence = TRUE) %>% 
  RunUMAP(.,  reduction = "harmony", dims = 1:15) %>% 
  RunTSNE(.,  reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(., reduction = "harmony", dims = 1:15) %>% 
  FindClusters(., resolution = 0.5)

Idents(mye)<-"cell"
DimPlot(mye)
# save(mye,file = "mye_back.Rda")
DotPlot(mye,features = c("PLAUR"))

DimPlot(mye)


CellDimPlot(mye,group.by = "celltype",label = T,theme_use = "theme_blank",
            add_mark = F,mark_expand = unit(0.1, "mm"),show_stat = F,
            mark_type = c("rect"),palcolor=c("#E1A39C","#9C623D","#4C74A7","#481F72"))


mye$celltype<-as.character(mye$cell) %>% {
  .[.=='P2RY12+(Microglia)'] <- "Microglia"
  .[.=='CD87+'] <- "TAM3"
  .[.=='GPM6B+'] <- 'TAM2'
  .[.=='FTL+'] <- "TAM1"
  ;.}

Idents(mye)<-"celltype"
DimPlot(mye)

mye<-AddModuleScore(mye,gene_sets)

colnames(mye@meta.data)[31:33]<-c("AllSig","MTSig","TMESig")

library(ggpubr)

VlnPlot(mye,features = "TMESig",cols = c("#E1A39C","#9C623D","#4C74A7","#481F72"),pt.size = 0)+
  xlab("")+NoLegend()+scale_y_continuous(limits = c(-0.1,0.25))+ylab("TME signature scores")+
  stat_compare_means(comparisons = list(c("TAM3","TAM1"),
                                        c("TAM3","TAM2"),
                                        c("TAM3","Microglia")),
                     label.y =0.2)


mye.marker<-FindAllMarkers(mye,logfc.threshold = 0)

library(scRNAtoolVis)
# plot
markerVocalno(markers = mye2.marker,
              topn = 5,
              labelCol = c("#E1A39C","#9C623D","#4C74A7","#481F72"))

# mymarkerVocalno<-function (markers = NULL, ownGene = NULL, topn = 5, log2FC = 0.25, 
#               labelCol = NULL, hlineSize = 1, hlineColor = "grey50", pforce = 5, 
#               nforce = 2.5, nudge_x = 0.8, pnudge_y = 0.25, nnudge_y = 0, 
#               base_size = 14, facetColor = NA, facetFill = "white", ylab = "Log2-Fold Change", 
#               nrow = 1) 
# {
#   if (is.null(ownGene)) {
#     markers$dif <-markers$pct.1 - markers$pct.2
#     toppos <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = topn, 
#                                                                     wt = dif)
#     topnegtive <- markers %>% dplyr::group_by(cluster) %>% 
#       dplyr::top_n(n = -topn, wt = dif)
#     topgene <- rbind(toppos, topnegtive)
#   }
#   else {
#     topgene <- markers %>% dplyr::filter(gene %in% ownGene)
#     toppos <- topgene %>% dplyr::filter(avg_log2FC > 0)
#     topnegtive <- topgene %>% dplyr::filter(avg_log2FC < 
#                                               0)
#   }
#   ggplot2::ggplot(markers, ggplot2::aes(x = pct.1 - pct.2, 
#                                         y = avg_log2FC)) + ggplot2::geom_point(color = "grey80") + 
#     ggplot2::geom_hline(yintercept = c(-log2FC, log2FC), 
#                         lty = "dashed", size = hlineSize, color = hlineColor) + 
#     ggrepel::geom_text_repel(data = toppos, ggplot2::aes(x = pct.1 - 
#                                                            pct.2, y = avg_log2FC, label = gene, color = cluster), 
#                              show.legend = F, direction = "y", hjust = 1, nudge_y = pnudge_y, 
#                              force = pforce, nudge_x = -nudge_x - (toppos$pct.1 - 
#                                                                      toppos$pct.2)) + ggrepel::geom_text_repel(data = topnegtive, 
#                                                                                                                ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC, label = gene, 
#                                                                                                                             color = cluster), show.legend = F, direction = "y", 
#                                                                                                                hjust = 0, nudge_y = nnudge_y, force = nforce, nudge_x = nudge_x - 
#                                                                                                                  (topnegtive$pct.1 - topnegtive$pct.2)) + ggplot2::geom_point(data = topgene, 
#                                                                                                                                                                               show.legend = F, ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC, 
#                                                                                                                                                                                                             color = cluster)) + ggplot2::scale_color_manual(name = "", 
#                                                                                                                                                                                                                                                             values = labelCol) + ggplot2::theme_bw(base_size = base_size) + 
#     ggplot2::theme(panel.grid = ggplot2::element_blank(), 
#                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
#                    strip.background = ggplot2::element_rect(color = facetColor, 
#                                                             fill = facetFill)) + ggplot2::xlab(expression(Delta ~ 
#                                                                                                             "Percentage Difference")) + ggplot2::ylab(ylab) + ggplot2::facet_wrap(~cluster, 
#                                                                                                                                                                                   nrow = nrow, scales = "fixed")
# }
ggsave("marker_mye.pdf",height = 4,width = 10)

load("glioma_Sig.Rda")
list<-list(sig=sc.marker$gene,
           top=gene_glioma)

library(GSVA)
library(Seurat)
gsva_matrix<- gsva(as.matrix(exp),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix<-as.data.frame(t(gsva_matrix))

clin$Sig<-as.numeric(gsva_matrix$sig)
clin$MTS<-as.numeric(gsva_matrix$top)

table(clin$Recurrence)
clin$Rec<-ifelse(clin$Recurrence %in% c("Recurrent","Secondary"),"Recurrent",as.character(clin$Recurrence))
table(clin$Rec)
clin<-clin[!is.na(clin$Rec),]

p1=ggboxplot(clin[!is.na(clin$Rec),],
             x="Rec",
             y="Sig",add = "jitter",
             fill = "Rec",
             font.label = list(size = 20, color = "black"),
             palette = c("#5BA4AE","#951F41"),
             add.params = list(size=0.5))+
  stat_compare_means(comparisons = list(c("Primary","Recurrent")))+
  #scale_fill_manual(values = c("#BCDE93","#D17D74","#947A97"))+
  xlab("")+ylab("TAM3 Signature Scores")+NoLegend()+
  theme(axis.text.x = element_text(size = 12,colour = "black",angle = 45,hjust = 1))+
  scale_y_continuous(limits = c(0,1.1))
p1

clin$MTS<-as.numeric(gsva_matrix$top)

table(clin$Grade)
clin$Grade2021<-ifelse(clin$Grade=="II",2,ifelse(clin$Grade=="III",3,ifelse(clin$Grade=="IV",4,as.character(clin$Grade))))

p2=ggboxplot(clin[!is.na(clin$Grade2021),],
             x="Grade2021",y="Sig",
             add = "jitter",
             fill ="Grade2021",
             #facet.by = "Rec",
             ylab = "TAM3 Signature scores",
             font.label = list(size = 20, color = "black"),
             add.params = list(size=0.5),
             xlab = "",
             palette = c("#BED99A","#C28177","#8F7B92"))+NoLegend()+
  stat_compare_means(comparisons = list(c("2","3"),
                                        c("3","4"),
                                        c("2","4")),
                     method = "wilcox.test")+
  theme(axis.text.x = element_text(size = 12,colour = "black",angle = 45,hjust = 1))+
  scale_y_continuous(limits = c(0,1.1))
p2

p1+p2

p3=ggscatter(clin,x="MTS",y="Sig",
             add = "reg.line", conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             ylab = "TAM3 Signature scores",
             xlab="MTS scores",
             shape = "Rec",
             color = "Rec",
             palette = c("#6596A5","#7B263A")
)+ stat_cor(method = "pearson",color='black',
            label.y.npc = "top")+
  theme(legend.position = "bottom")+
  scale_y_continuous(limits = c(0,1.1))
p3
library(patchwork)
p1+p2+p3+plot_layout(widths = c(1,1.5,1.5))



library(Seurat)
rm(sc)
sc<-readRDS("GSM6261341_astrocytoma_paired_data_merged_seurat_object.rds")

DimPlot(sc2)

library(SCP)
load("/Users/mac/Documents/scRNA/glioma/glass/sce_after.Rda")
sce$cell<-as.character(sce$main_cluster)
table(sce$cell)

max(sc$nCount_RNA);min(sc$nCount_RNA)
max(sc$nFeature_RNA);min(sc$nFeature_RNA)

# minGene=500
# maxGene=4000
# pctMT=15
# 
# table(sc$nFeature_RNA<4000)
# 
# sc<-subset(sc,nFeature_RNA<4000)

library(harmony)
sc@meta.data[["gem_id"]]
DimPlot(sc)
meta<-sc@meta.data
sc<-CreateSeuratObject(sc@assays$RNA@counts,
                       meta.data = meta)

sc <- sc %>% 
  NormalizeData(., 
                normalization.method = "LogNormalize", 
                scale.factor = 10000) %>%
  FindVariableFeatures(., selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(., features = VariableFeatures(.)) %>%
  RunPCA(.,  npcs = 25) %>% 
  RunHarmony(.,reduction = "pca",
             dims.use = 1:25,
             group.by.vars = "gem_id",
             plot_convergence = TRUE) %>% 
  RunUMAP(.,  reduction = "harmony", dims = 1:15) %>% 
  RunTSNE(.,  reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(., reduction = "harmony", dims = 1:15) %>% 
  FindClusters(., resolution = 0.5)

sc@meta.data[["Subclusters"]]
DimPlot(sc,group.by = "Subclusters")
sc <- RunKNNPredict(
  srt_query = sc, srt_ref = sce,
  ref_group = "cell", filter_lowfreq = 20
)
colnames(sc@meta.data)


table(sc$gem_id)

sc$patient <- as.character(sc$gem_id) %>% {
  .[.%in% c('NCH557','NCH758A')] <- 'p1';
  .[.%in% c('NCH511B','NCH678')] <- 'p2';
  .[.%in% c('NCH302','NCH645')]  <- 'p3';  # MT
  .[.%in% c('NCH988','NCH2375')] <- 'p4';
  .[.%in% c('NCH740W','NCH2367')] <- 'p5'; # low quality
  .[.%in% c('NCH673D','NCH2260')] <- 'p6'; # MT
  .}

sc$MT <- as.character(sc$patient) %>% {
  .[.%in% c('p1','p2','p4','p5')] <- 'Not-MT';
  .[.%in% c('p3','p6')] <- 'MT';
  .}

sc$id<-paste0(sc$patient,"_",sc$time)

black_theme <- theme(panel.background = element_rect(fill = "black", colour = "white",linewidth=1),
                     plot.background = element_rect(fill = "black", colour = "white"),
                     panel.grid.major = element_line(color = "black"),
                     panel.grid.minor = element_line(color = "black"),
                     text = element_text(color = "white"),
                     axis.text = element_text(color = "white"),
                     axis.ticks = element_line(color = "white"),
                     axis.title = element_text(color = "white"),
                     legend.background = element_rect(fill = "black", colour = "black"),
                     legend.text = element_text(color = "white"),
                     legend.title = element_text(color = "white"))

cols<-colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", 
                  "#fdbf6f", "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", 
                  "#cab2d6", "#ffff99")


library(ggrastr)
pdf("umap2.pdf",height = 5,width = 6.5)
rasterize(p, layers = 'Point', dpi = 600) #先调一个低分辨率，便于直观对比
dev.off()

p<-DimPlot(sc,group.by = "id",cols = cols)+black_theme
p<-DimPlot(sc,group.by = "KNNPredict_classification")+black_theme


load("scemk.Rda")
scemk2<-scemk[which(scemk$cluster=="positive"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

FeaturePlot(sc,features = "PLAUR",split.by = "time")

sc<-AddModuleScore(sc,list(sc.marker$gene))

FeatureDimPlot(sc,features = "Cluster1",split.by = "MT")
colnames(sc@meta.data)

mye<-subset(sc,Subclusters=="Microglia")
mye <- mye %>% 
  NormalizeData(., 
                normalization.method = "LogNormalize", 
                scale.factor = 10000) %>%
  FindVariableFeatures(., selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(., features = VariableFeatures(.)) %>%
  RunPCA(.,  npcs = 25) %>% 
  RunHarmony(.,reduction = "pca",
             dims.use = 1:25,
             group.by.vars = "gem_id",
             plot_convergence = TRUE) %>% 
  RunUMAP(.,  reduction = "harmony", dims = 1:15) %>% 
  RunTSNE(.,  reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(., reduction = "harmony", dims = 1:15) %>% 
  FindClusters(., resolution = 0.5)

DimPlot(mye,label = T)
mye<-subset(mye,seurat_clusters %in% c(0:8))
mye<-AddModuleScore(mye,list(sc.marker$gene))
colnames(mye@meta.data)[24]<-"SenTAM.sig"

FeatureDimPlot(mye,features = "SenTAM.sig",
               split.by = "time",reduction = "TSNE")


mye$group<-paste0(mye$MT,"_",mye$time)

FeatureDimPlot(mye,features = "SenTAM.sig",split.by = "group",reduction = "TSNE")
df <- VlnPlot(mye,features = "SenTAM.sig",
              group.by = "group",pt.size = 0)

dput(unique(mye$group))
c("MT_relapse", "Not-MT_relapse", "MT_primary", "Not-MT_primary")

df+xlab("Myeloid subsets")+stat_compare_means(label.y = 0.25,
                                              color = "white",
                                              comparisons = list(c("MT_relapse","MT_primary"),
                                                                 c("Not-MT_relapse","Not-MT_primary"),
                                                                 c("MT_relapse","Not-MT_relapse")))+
  black_theme+
  NoLegend()

