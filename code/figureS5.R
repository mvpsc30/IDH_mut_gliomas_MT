```R
library(tidyverse)
library(ggpubr)
library(data.table)
library(openxlsx)
library(GSVA)
library(ggpubr)
library(clusterProfiler)
load("scemk.Rda")
scemk2<-scemk[which(scemk$cluster=="positive"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

scemk2<-scemk[which(scemk$cluster=="Bg"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker2<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

phe<-read.xlsx("ClinicalInformation.xlsx")

data1<-read.table("panglioma544.initial.RPKM.txt",header = T,row.names = 1)
data2<-read.table("panglioma544.recurrence.RPKM.txt",header = T,row.names = 1)

all<-cbind(data1,data2)

all_phe<-data.frame(id=colnames(all),
                    id2=substring(colnames(all),1,5),
                    group=substring(colnames(all),7,7))

id<-read.csv("idh_regrading.csv",header = T)
id$name<-paste0(id$id,"_",id$group)
table(id$MT)

exp <-all[,id$name]

max(exp);min(exp)

exp<-log2(exp+1)

gs_list<-readRDS("/Users/mac/Downloads/difg_glass_2019/gene_set_list.Rds")

gs_list<-gs_list$charoentong

sig<-openxlsx::read.xlsx("/Users/mac/Downloads/difg_glass_2019/sig_immune.xlsx")

list<-list(MDSC=gs_list$MDSC,
           SenTAM=sc.marker$gene,
           otherTAM=sc.marker2$gene,
           Exhausted=c(na.omit(sig$Exhausted)),
           Treg=gs_list$`Regulatory T cell`
)

gsva_ex<-GSVA::gsva(expr=as.matrix(exp), 
                    gset.idx.list = list, 
                    method="ssgsea",verbose= FALSE)

gsva_ex<-as.data.frame(t(gsva_ex))

id$SenTAM<-gsva_ex$SenTAM
id$MDSC<-gsva_ex$MDSC
id$OtherTAM<-gsva_ex$otherTAM
id$Exhausted<-gsva_ex$Exhausted
id$Treg<-gsva_ex$Treg

id2<-id[65:128,]

p1=ggscatter(id2,x="OtherTAM",y="Exhausted",
             add = "reg.line", 
             conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             xlab = "OtherTAM.sig",
             ylab = "Exhausted T.sig"
)+ stat_cor(method = "pearson",
            label.y.npc = "top")

p2=ggscatter(id2,x="OtherTAM",y="Treg",
             add = "reg.line", 
             conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             xlab = "OtherTAM.sig",
             ylab = "Regulatory T.sig"
)+ stat_cor(method = "pearson",
            label.y.npc = "top")

p3=ggscatter(id2,x="OtherTAM",y="MDSC",
             add = "reg.line", 
             conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             xlab = "OtherTAM.sig",
             ylab = "MDSC.sig"
)+ stat_cor(method = "pearson",
            label.y.npc = "top")

f1=p1/p2/p3

library(patchwork)

p4=ggscatter(id2,x="SenTAM",y="Exhausted",
             add = "reg.line", 
             conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             xlab = "SenTAM.sig",
             ylab = "Exhausted T.sig"
)+ stat_cor(method = "pearson",
            label.y.npc = "top")

p5=ggscatter(id2,x="SenTAM",y="Treg",
             add = "reg.line", 
             conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             xlab = "SenTAM.sig",
             ylab = "Regulatory T.sig"
)+ stat_cor(method = "pearson",
            label.y.npc = "top")

p6=ggscatter(id2,x="SenTAM",y="MDSC",
             add = "reg.line", 
             conf.int = TRUE,
             add.params = list(fill = "lightgray"),
             xlab = "SenTAM.sig",
             ylab = "MDSC.sig"
)+ stat_cor(method = "pearson",
            label.y.npc = "top")

f2=p4/p5/p6
```
```R
TME_data <- read.csv("CIBERSORTx_Job1_Results (1).csv",header = T,row.names = 1)
head(TME_data)

TME_data = select(TME_data, stemcell_tumor:b_cell)

# 
TME_New = TME_data[id$name,]

id$cell<-TME_New$prolif_stemcell_tumor*100

phe_sub<-TME_New[,c(1,3,4)]*100
colnames(phe_sub)<-c("Stem-like","Diff-like","Prolife.Stem-like")

phe_sub$name<-rownames(phe_sub)
phe_sub$group<-id$group

phe_sub<-melt(phe_sub)
phe_sub$id<-substring(phe_sub$name,1,5)
colnames(phe_sub)
id$Grade.of.Recurrent.Tumor
phe_sub2<-left_join(phe_sub,id[,c("name","MT",
                                  "Grade.of.Initial.Tumor",
                                  "Grade.of.Recurrent.Tumor")],by="name")

df<-TME_New
df<-df[order(df$`Prolife.Stem-like`,decreasing = T),]

df<-df[,c("Prolife.Stem-like","Stem-like", "Diff-like",
          "Granulocyte", "T cells",  "Pericyte", 
          "Oligodendrocyte","Myeloid", "Endothelial", 
          "DCs", "Fibroblast", "B cells")]

dd1 <- df %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols=2:13,
               names_to= "celltype",
               values_to = "Proportion")

dd1$sample<-factor(dd1$sample,levels = rownames(df))

library(ggplot2)
library(RColorBrewer)

plot_order<-c("Prolife.Stem-like","Stem-like", "Diff-like", 
              "Granulocyte", "T cells",  "Pericyte", 
              "Oligodendrocyte","Myeloid", "Endothelial", 
              "DCs", "Fibroblast", "B cells")

dd1$celltype<-factor(dd1$celltype,levels = rev(plot_order))

dput(colorRampPalette(brewer.pal(15,"Set3"))(15))

Cols <-c("#EEF0FC", "#BBD7E4", "#66AED5", "#2A80BC", "#0C4D9A", "#24A15A", 
         "#FCFED4", "#FFE08D", "#FABB9F", "#F9CDE4", "#F9684A", "#E41A1C", 
         "#C2ADC0", "#D6EBB2", "#FFED6F")

ggplot(dd1,aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity")+
  scale_fill_manual(values=Cols) +
  theme_classic()+
  guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("")

x= aggregate(dd1$Proportion,by=list(dd1$sample),sum)
x$sample<-x$Group.1

dd1<-left_join(dd1,x,by="sample")
dd1$Proportion<-dd1$Proportion/dd1$x


p1=ggplot(dd1[1:576,],aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity")+
  scale_fill_manual(values=Cols) +
  theme_classic()+
  guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("")+ylab("Primary")


p2=ggplot(dd1[577:1152,],aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity")+
  scale_fill_manual(values=Cols) +
  theme_classic()+
  guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("")+ylab("Recurrent")
library(patchwork)
p1+p2+plot_layout(guides = "collect",ncol = 1)

```
```R
p1=ggpaired(phe_sub2, 
            x="group", 
            y="value", 
            fill="group",
            id = "id",
            add="jitter",
            line.color = "gray", line.size = 0.5,
            palette=c("#BCDE93","#D17D74"),
            facet.by = c("MT","variable"),
            #scales = "free_y",
            title = "Tumor cell subsets",
            xlab=" ",
            ylab = "Proportion(%)",
            legend.title=" ",show.legend = F,
            ncol=6) + 
  stat_compare_means(paired = TRUE) +#配对t检验
  theme(legend.position = 'none')

p2=ggpaired(id[which(id$MT!="No"),], 
            x="group", 
            y="cell", 
            fill="group",
            id = "id",
            add="jitter",
            line.color = "gray", line.size = 0.5,
            palette=c("#BCDE93","#D17D74"),
            facet.by = "grade",
            #scales = "free_y",
            title = "Grade change(Prolife.Stem-like)",
            xlab=" ",
            ylab = "Proportion(%)",
            legend.title=" ",
            show.legend = F,
            ncol=6) + 
  stat_compare_means(paired = TRUE) +#配对t检验
  theme(legend.position = 'none')

p1+(p3/p2)

ggsave("box_combine_fig1.pdf",height =5,width = 12)

ggboxplot(data = id[which(id$group!="Primary"),],
          x="Grade.of.Recurrent.Tumor",
          y="sig",fill="Grade.of.Recurrent.Tumor",
          title = "Grade of Recurrent Glioma",
          xlab=" ", ylab = "Proportion(%)",
          legend.title=" ",show.legend = F,
          add = "jitter",
          palette = c("#BCDE93","#D17D74","#846E89"))+
  stat_compare_means(comparisons = list(c("2","3"),
                                        c("2","4"),
                                        c("3","4"))) +
  theme(legend.position = 'none')
```
```R
library(semla)
library(tibble)
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)
library(patchwork)
library(ggplot2)
library(RcppML)
load("scemk.Rda")
scemk2<-scemk[which(scemk$cluster=="positive"),]
scemk2$difference <- scemk2$pct.1 - scemk2$pct.2
scemk2_sig <- scemk2[which(scemk2$p_val_adj<0.05 & abs(scemk2$avg_log2FC) >=1),]
sc.marker<-scemk2_sig[which(scemk2_sig$difference>=0.2 & scemk2_sig$avg_log2FC>=1),]

# Assemble spaceranger output files
samples <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM4_spaceranger_out/filtered_feature_bc_matrix.h5")
imgs <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM4_spaceranger_out/spatial/tissue_hires_image.png")
spotfiles <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM4_spaceranger_out/spatial/tissue_positions_list.csv")
json <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM4_spaceranger_out/spatial/scalefactors_json.json")


samples <- Sys.glob("/Users/mac/Downloads/ST_matrix/Pt7/filtered_feature_bc_matrix.h5")
imgs <- Sys.glob("/Users/mac/Downloads/ST_matrix/Pt7/spatial/tissue_hires_image.png")
spotfiles <- Sys.glob("/Users/mac/Downloads/ST_matrix/Pt7/spatial/tissue_positions_list.csv")
json <- Sys.glob("/Users/mac/Downloads/ST_matrix/Pt7/spatial/scalefactors_json.json")


samples <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM1_spaceranger_out/filtered_feature_bc_matrix.h5")
imgs <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM1_spaceranger_out/spatial/tissue_hires_image.png")
spotfiles <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM1_spaceranger_out/spatial/tissue_positions_list.csv")
json <- Sys.glob("/Users/mac/Documents/scRNA/data_glioma/GSE194329_RAW/GBM1_spaceranger_out/spatial/scalefactors_json.json")


# Create a tibble/data.frame with file paths
infoTable <- tibble(samples, imgs, spotfiles, json)

# Create Seurat object
se <- ReadVisiumData(infoTable = infoTable)

se <- LoadImages(se)

se <- se |>
  NormalizeData() |>
  ScaleData() |>
  FindVariableFeatures() |>
  RunPCA() |>
  FindNeighbors(reduction = "pca", dims = 1:30) |>
  FindClusters()

MapLabels(se, column_name = "seurat_clusters", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

se$cluster_9 <- ifelse(se$seurat_clusters %in% "3", "3", NA)

MapLabels(se, column_name = "cluster_9", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

se <- RadialDistance(se, column_name = "seurat_clusters", selected_groups = "3")

MapFeatures(se, features = "r_dist_3", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

se$r_dist_3_sqrt <- sign(se$r_dist_3)*sqrt(abs(se$r_dist_3))
MapFeatures(se, features = "r_dist_3_sqrt", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

library(tidyr)

sel_genes <- c("TOP2A", "IL1B")

se[[]] |> 
  bind_cols(FetchData(se, vars = sel_genes)) |> 
  filter(r_dist_3 < 1e3) |> 
  pivot_longer(all_of(sel_genes), names_to = "variable", values_to = "value") |> 
  ggplot(aes(r_dist_3, value, color = variable)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed")


MapFeatures(SubsetSTData(se, expression = r_dist_3 < 1e3), 
            features = sel_genes, scale_alpha = TRUE, ncol = 3,
            image_use = "raw", pt_size = 1.5, colors = viridis::viridis(n = 9),
            max_cutoff = 0.99)


se <- RegionNeighbors(se, 
                      column_name = "seurat_clusters", 
                      column_labels = "3")

MapLabels(se, 
          #crop_area = c(0.45, 0.6, 0.8, 0.93),
          column_name = "nb_to_3", drop_na = TRUE,
          image_use = "raw", pt_size = 3)

se <- RegionNeighbors(se, column_name = "seurat_clusters", 
                      column_labels = "3", mode = "inner_outer")

MapLabels(se, 
          #crop_area = c(0.45, 0.6, 0.8, 0.93),
          column_name = "nb_to_3",
          image_use = "raw", pt_size = 3, drop_na = TRUE)

VlnPlot(se,features ="SPP1",group.by = "nb_to_3")
VlnPlot(se,features ="Cluster1",group.by = "nb_to_3")


se<-AddModuleScore(se,list(sc.marker$gene))

cols <- viridis::rocket(11, direction = -1)

MapFeatures(se, 
            features = c("TOP2A","PLAUR","CD8A"),
            colors = cols)


p1=MapLabels(se, 
             #crop_area = c(0.45, 0.6, 0.8, 0.93),
             column_name = "nb_to_3",
             image_use ="raw", pt_size = 3, drop_na = TRUE)
p+p1

FeaturePlot(sce,"PDGFB",reduction = "tsne")
FeaturePlot(sce,"PDGFRA",reduction = "tsne")

FeaturePlot(mye,"LGALS9")
# Compute spatial autocorrelation for variable features
spatgenes <- CorSpatialFeatures(se)

loc_meta<-se@meta.data
loc_df$Sig<-loc_meta$Cluster1

cor_matrix <- loc_df |> 
  mutate_all(~ if_else(.x<0.01, 0, .x)) |>  # Filter lowest values (-> set as 0)
  cor()

diag(cor_matrix) <- NA
max_val <- max(cor_matrix, na.rm = T)
cols <- RColorBrewer::brewer.pal(7, "RdYlBu") |> rev(); cols[4] <- "white"
pheatmap::pheatmap(cor_matrix, 
                   breaks = seq(-max_val, max_val, length.out = 100),
                   color=colorRampPalette(cols)(100),
                   cellwidth = 14, cellheight = 14, 
                   treeheight_col = 10, treeheight_row = 10, 
                   main = "Cell type correlation\nwithin spots")



DefaultAssay(se) <- "Spatial"

ti <- Sys.time()
unique(sce$main_cluster)
unique(sce$celltype)

# sce2<-subset(sce,main_cluster %in% c("Stem-like", "Diff-like", "Prolife.stem-like", "Oligodendrocyte", 
#                                      "Endothelial", "T cell", "CD87+", 
#                                      "Fibroblast","P2RY12+(Microglia)"))
sce <-sce |>
  NormalizeData() |>
  ScaleData() |>
  FindVariableFeatures() |>
  RunPCA()
sce <-sce |> RunUMAP(dims = 1:15)

DimPlot(sce,reduction = "pca")

se@assays$celltypeprops<-NULL

se <- RunNNLS(object = se, 
              singlecell_object = sce, 
              groups = "celltype")

# Plot selected cell types
DefaultAssay(se) <- "celltypeprops"
dput(unique(sce$celltype))

df<-as.data.frame(t(se@assays$celltypeprops@data))
as.numeric(df$TAM3)

selected_celltypes <-c("Stem-like", "Diff-like", "Prolife.stem-like", "TAM1", 
                       "TAM3","TAM2")

se<- LoadImages(se, image_height = 1e3)

plots <- lapply(seq_along(selected_celltypes), function(i) {
  MapFeatures(se, pt_size = 1.3,
              features = selected_celltypes[i], image_use = "raw",
              arrange_features = "row", scale = "shared", 
              override_plot_dims = TRUE,
              colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
              scale_alpha = TRUE)  +
    plot_layout(guides = "collect") & 
    theme(legend.position = "right", legend.margin = margin(b = 50),
          legend.text = element_text(angle = 0),
          plot.title = element_blank())
}) |> setNames(nm = selected_celltypes)

plots$TAM3

MapFeatures(se, pt_size = 1.3,
            features = c("Prolife.stem-like",
                         "TAM3"),
            scale_alpha = F)  +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())

# Plot multiple features
MapMultipleFeatures(se, 
                    image_use = "raw", 
                    pt_size = 2, max_cutoff = 0.99,
                    override_plot_dims = T, 
                    colors = c("#ADCCE0", "#3B76AF", "#BBDE93", "#D7A69E", "#946342", "red"),
                    features = selected_celltypes) +
  plot_layout(guides = "collect")&
  theme(panel.grid.major = element_line(linetype = "dashed"), axis.text = element_text())


ggsave("sGBM1.pdf",height = 6,width = 10)

data<-readRDS("data/TCGA_GBMLGG.Rds")
data<-data$pData
table(data$IDH.status)

MapMultipleFeatures(se, 
                    image_use = "raw", 
                    pt_size = 4, max_cutoff = 0.99,
                    crop_area = c(0.3, 0.7, 0.5, 0.85),
                    colors = c("#ADCCE0", "#3B76AF", "#BBDE93", "#D7A69E", "#946342", "red"),
                    features = selected_celltypes) +
  plot_layout(guides = "collect")

ggsave("sGBM1_sub.pdf",height = 4,width = 5)

cols_he <- viridis::viridis(11)

p <- MapFeatures(se, 
                 features = selected_celltypes[6], 
                 image_use = "raw", 
                 color = cols_he, 
                 pt_alpha = 0.5) &
  theme(panel.grid.major = element_line(linetype = "dashed"), axis.text = element_text())
p


MapMultipleFeatures(se, 
                    image_use = "raw", 
                    pt_size = 2.5, max_cutoff = 0.99,
                    crop_area = c(0.5, 0.5, 0.7, 0.7),
                    colors = c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
                    features = selected_celltypes) +
  plot_layout(guides = "collect") &
  theme(plot.title = element_blank()) &
  theme(panel.grid.major = element_line(linetype = "dashed"), axis.text = element_text())

MapFeatures(se, 
            features = c("PLAUR","ITGAV"), 
            image_use = "raw",
            color = cols_he)


se$group<-ifelse(is.na(se$nb_to_3),"Bg",se$nb_to_3)


cor_matrix <- FetchData(se, selected_celltypes) |> 
  mutate_all(~ if_else(.x<0.0000001,0, .x)) |>  # Filter lowest values (-> set as 0)
  cor()
df<-na.omit(df)
cor_matrix <-cor(df)
diag(cor_matrix) <- NA
max_val <- max(cor_matrix, na.rm = T)
cols <- RColorBrewer::brewer.pal(7, "RdYlBu") |> rev(); cols[4] <- "white"

pdf("cell_correlation.pdf",height = 5,width = 5)
pheatmap::pheatmap(cor_matrix, 
                   breaks = seq(-max_val, max_val, length.out = 100),
                   color=colorRampPalette(cols)(100),
                   cellwidth = 14, cellheight = 14, 
                   treeheight_col = 10, treeheight_row = 10, 
                   main = "Cell type correlation\nwithin spots")
dev.off()
```
  
  