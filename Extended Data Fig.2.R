# heatmap
load("/mnt/wutan/data/9.CMSPlus/Figure1/gsva.all.42.RData")
load("/mnt/wutan/data/9.CMSPlus/code/basement/CRCSC.expression.clindata.RData")
colannoA = HeatmapAnnotation(
    CMS= clindata$cms_label, 
    col = list(CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"))
)
col_fun = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#368ABF", "#68C3A7", "#AED8A3", "#E3EA9B", "white", "#FDE18B", "#F9AF62", "#F36D46", "#D54150"))(35))
col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))
dendA = cluster_within_group(gsva.all.42, clindata$cms_label)
PPCRCSC2 <- Heatmap(gsva.all.42[1:15, ], col = col_fun, cluster_rows = FALSE, cluster_columns = dendA, top_annotation = colannoA, column_labels=rep("", dim(gsva.all.42)[2]))
PPCRCSC3 <- Heatmap(gsva.all.42[16:42, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = dendA, column_labels=rep("", dim(gsva.all.42)[2]))
ht_list = PPCRCSC2%v%PPCRCSC3
draw(ht_list)
# cluster
mat <- gsva.all.42
factor=clindata$cms_label
if (!is.factor(factor)) {
     factor = factor(factor, levels = unique(factor))
}
le <- "CMS1"
m = mat[, factor == le, drop = FALSE]
dim(m)
hc1 = hclust(dist(1:ncol(m)))
plot(hc1)
