### Figure1a
library(ComplexHeatmap)
library(circlize)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/model/intoeder.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")
Datasetcolor=RColorBrewer::brewer.pal(10, "Set3")[-9]
names(Datasetcolor) <- names(table(CMSclin2$dataset))
col_fun = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#368ABF", "#68C3A7", "#AED8A3", "#E3EA9B", "white", "#FDE18B", "#F9AF62", "#F36D46", "#D54150"))(35))
col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))
colannoA4 = HeatmapAnnotation(
		CMSPlus=CMSclin2$NewCluster,
		CMS= CMSclin2$cms_label, 
		Dataset=CMSclin2$dataset,
		MSI=CMSclin2$msi, 
		CIMP=CMSclin2$cimp, 
		KRAS=CMSclin2$kras_mut, 
		BRAF=CMSclin2$braf_mut, 
		Age=CMSclin2$age,
		Gender=CMSclin2$gender,
		Grade=CMSclin2$grade,
		Pathologic_Tstage=CMSclin2$pt, 
		Pathologic_Nstage=CMSclin2$pn, 
		Pathologic_Mstage=CMSclin2$pm, 
		pathologic_TNMstage=CMSclin2$tnm, 
    col = list(CMSPlus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', `CMS4-TME+`='#3C5488', `CMS4-TME-`='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20'), Dataset=Datasetcolor, Gender=c(female="#E69F00", male="#56B4E9"), Grade=c('1'=RColorBrewer::brewer.pal(10, "Blues")[2], '2'=RColorBrewer::brewer.pal(10, "Blues")[4], '3'=RColorBrewer::brewer.pal(10, "Blues")[6], '4'=RColorBrewer::brewer.pal(10, "Blues")[8]), Pathologic_Tstage=c('1'=RColorBrewer::brewer.pal(10, "YlGn")[2], '2'=RColorBrewer::brewer.pal(10, "YlGn")[4], '3'=RColorBrewer::brewer.pal(10, "YlGn")[6], '4'=RColorBrewer::brewer.pal(10, "YlGn")[8]), Pathologic_Nstage=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(10, "Paired")[2]), Pathologic_Mstage=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(12, "Paired")[12]), pathologic_TNMstage=c("I"=RColorBrewer::brewer.pal(10, "BuPu")[2], "II"=RColorBrewer::brewer.pal(10, "BuPu")[4], "III"=RColorBrewer::brewer.pal(10, "BuPu")[6], "IV"=RColorBrewer::brewer.pal(10, "BuPu")[8]), Age=colorRamp2(c(20, 100), c("white", "red")))
)
row.subsections1 <- c(3, 4, 5, 3)
row_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4"), row.subsections1))
col.subsections1 <- as.numeric(table(CMSclin2$NewCluster))
col_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4-TME+", "CMS4-TME-"), col.subsections1))
PPCRCSC2 <- Heatmap(gsvaall2[1:15, ], col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = colannoA4, column_labels=rep("", dim(CMSclin2)[1]), row_split = row_split1, column_split = col_split1)
row.subsections2 <- c(8, 19)
row_split2 = data.frame(factor(rep(c("Fibroblastes", "Immune"), row.subsections2), levels=c("Immune", "Fibroblastes")))
PPCRCSC3 <- Heatmap(gsvaall2[16:42, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, column_labels=rep("", dim(CMSclin2)[1]), row_split = row_split2, column_split = col_split1)
ht_list = PPCRCSC2%v%PPCRCSC3
draw(ht_list)

# Figure1b
library(EPIC)
setwd("/mnt/wutan/data/9.CMSPlus/data")
load("CMSclin2.gsvaall2.RData")
load("CRCSC.expression.clindata.RData")
expcrctcga <- expall.corr[, CMSclin2$sample]
out <- EPIC(bulk = expcrctcga)
CRCEPIC.pst <- PST(as.matrix(out$cellFractions))
CRCEPIC.pst <- data.frame(Sample=CRCEPIC.pst[, 1], ImmuneInfiltration=CRCEPIC.pst[, 2], Score=as.numeric(CRCEPIC.pst[, 3]))
CRCEPIC.pst$CMSPlus <- CMSclin2$NewCluster[match(CRCEPIC.pst$Sample, CMSclin2$sample)]
CRCEPIC.pst$CMSPlus <- gsub("CMS4_LowIF", "CMS4-TME-", CRCEPIC.pst$CMSPlus)
CRCEPIC.pst$CMSPlus <- gsub("CMS4_HighIF", "CMS4-TME+", CRCEPIC.pst$CMSPlus)
anno_CRCEPIC <- compare_means(Score ~ CMSPlus, method="wilcox.test", group.by = "ImmuneInfiltration", data = CRCEPIC.pst, p.adjust.method = "BH") %>%
 mutate(y_pos = 0.25, p.adj = format.pval(p.adj, digits = 2))
anno_CRCEPIC <- anno_CRCEPIC[anno_CRCEPIC$group1 %in% c("CMS4-TME+", "CMS4-TME-") & anno_CRCEPIC$group2 %in% c("CMS4-TME+", "CMS4-TME-"), ]
CMSPlusSubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
CRCEPIC.pst$CMSPlus <- factor(CRCEPIC.pst$CMSPlus, levels=c("CMS1", "CMS2", "CMS3", "CMS4-TME+", "CMS4-TME-"))
ggboxplot(CRCEPIC.pst, x = "CMSPlus", y = "Score", fill = "CMSPlus",
          ggtheme = theme_bw(), palette = CMSPlusSubtype,
          outlier.shape = NA) +
  facet_wrap(~ImmuneInfiltration, ncol = 5, scales = "free_y") +
  stat_compare_means(
    comparisons = list(c("CMS4-TME+", "CMS4-TME-")),
    method = "wilcox.test",
    label = "p.format"
  ) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

# Figure1c
p <- as.data.frame(t(out$cellFractions))
shannon <- apply(p, 2, function(x) {
  x <- x[x > 0]  # 移除 0，避免 log(0)
  if (length(x) == 0) return(0)
  -sum(x * log(x))
})
# Pielou's evenness：标准化到 [0,1]
n_features <- nrow(gsvaall2)
evenness <- shannon / log(n_features)
tme_entropy <- data.frame(
  sample = colnames(gsvaall2),
  shannon = shannon,
  evenness = evenness
)
identical(tme_entropy$sample, CMSclin2$sample)
tme_entropy$NewCluster <- CMSclin2$NewCluster
my_comparisons = list( c("CMS4_HighIF", "CMS4_LowIF"), c("CMS4_HighIF", "CMS1"), c("CMS4_HighIF", "CMS2"), c("CMS4_HighIF", "CMS3"), c("CMS1", "CMS4_LowIF"), c("CMS2", "CMS4_LowIF"), c("CMS3", "CMS4_LowIF") )
bpc3 <- ggplot(tme_entropy, aes(x=NewCluster, y=shannon, fill=NewCluster)) + geom_boxplot(width=0.8, outlier.shape = NA) + 
      labs(title="Shannon entropy",x="Subtype", y = "Shannon entropy") + 
      scale_fill_manual(values=CMSPlusSubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "wilcox.test", step_increase = .1, tip_length = 0, vjust=0)

# Figure1d
library(maftools)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
TCGAclindata <- CRC.clin.exp.he_v2
TCGAclindata$Tumor_Sample_Barcode <- CRC.clin.exp.he_v2$sample
CRCMaf <- read.maf(maf = "/data/home2/wutan/Task/9.CMS.SDI/Procedures/StartFrom20210310/CRC_maftools.maf", clinicalData = TCGAclindata, verbose = FALSE, isTCGA = TRUE)
CRCMaf.subtype = clinicalEnrichment(maf = CRCMaf, clinicalFeature = 'CMSPlus')
CMSPlusSubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
MutationCountSample <- table(CRCMaf@data$Tumor_Sample_Barcode)
MutationCountDF <- data.frame(Sample=names(MutationCountSample), MutationCount=as.numeric(MutationCountSample), Subtype=TCGAclindata$CMSPlus[match(names(MutationCountSample), TCGAclindata$Tumor_Sample_Barcode)])
t.test(MutationCountDF$MutationCount[which(MutationCountDF$Subtype %in% "CMS4-TME+")], MutationCountDF$MutationCount[which(MutationCountDF$Subtype %in% "CMS4-TME-")])
# p-value = 0.01493
my_comparisons = list( c("CMS4-TME+", "CMS4-TME-") )
MutationCountDF <- MutationCountDF[!is.na(MutationCountDF$Subtype), ]
bp1 <- ggplot(MutationCountDF, aes(x=Subtype, y=MutationCount, fill=Subtype)) + geom_boxplot(width=0.8) + 
      labs(x="Subtype", y = "Mutation Count") + 
      scale_fill_manual(values=CMSPlusSubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(0, 1500)) +
      geom_signif(comparisons=my_comparisons, annotations=c("1.49e-02"), y_position = c(700), tip_length = 0, vjust=0)

# Figure1e
load("/data/home2/wutan/Task/MultiomicsSubtyping/data/TCGAbiolinks/data/CRCGistic.rda")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
TCGAclindata <- CRC.clin.exp.he_v2
TCGAclindata$Tumor_Sample_Barcode <- CRC.clin.exp.he_v2$sample
thresholedbygene <- gistic.thresholedbygene[, -(1:3)]
rownames(thresholedbygene) <- gistic.thresholedbygene[, 1]
thresholedbygene <- thresholedbygene[, grep("01", substr(colnames(thresholedbygene), 14, 15))]
colnames(thresholedbygene) <- substr(colnames(thresholedbygene), 1, 12)
colnames(thresholedbygene) <- gsub("\\.", "-", colnames(thresholedbygene))
for(i in 1:length(colnames(thresholedbygene))){
	thresholedbygene[, i] <- as.numeric(thresholedbygene[, i])
}
copynumberload <- apply(thresholedbygene, 2, function(x) sum(abs(x)))
SCNACountDF <- data.frame(Sample=names(copynumberload), SCNACount=as.numeric(copynumberload), Subtype=TCGAclindata$CMSPlus[match(names(copynumberload), TCGAclindata$Tumor_Sample_Barcode)])
##
t.test(SCNACountDF$SCNACount[which(SCNACountDF$Subtype %in% "CMS4-TME+")], SCNACountDF$SCNACount[which(SCNACountDF$Subtype %in% "CMS4-TME-")])
# p-value = 0.04313
my_comparisons = list( c("CMS4-TME+", "CMS4-TME-") )
SCNACountDF <- SCNACountDF[!is.na(SCNACountDF$Subtype), ]
bp2 <- ggplot(SCNACountDF, aes(x=Subtype, y=SCNACount, fill=Subtype)) + geom_boxplot(width=0.8) + 
      labs(x="Subtype", y = "SCNA Count") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(0, 25000)) +
      geom_signif(comparisons=my_comparisons, annotations=c("4.31e-02"), y_position = c(23000), tip_length = 0, vjust=0)

# Figure1f
setwd("/mnt/wutan/data/9.CMSPlus/Figure3/model")
load("CMSclin2.gsvaall2.RData")
load("CRC.clin.exp.he_v2.RData")
library(survival)
library(gfplot)
SurvivalFGPlotSimplify <- function(DataFrame, Time, Status, Lable, Color, type="days"){
	if(type=="days"){
		DataFrame[, Time] <- DataFrame[, Time]/30
	}
	clin <- Surv(DataFrame[, Time], DataFrame[, Status])
	labs <- factor(DataFrame[, Lable])
	ylable <- unlist(sapply(strsplit(Time, "\\_"), function(x) x[1]))
	parameters <- list(clin=clin, labs=labs, ylable=ylable, Color=Color)
	return(parameters)
}
CRC.clin.exp.he_v2$CMSPlus <- CMSclin2$NewCluster[match(CRC.clin.exp.he_v2$sample, CMSclin2$sample)]
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "old.rfs.delay", "old.rfs.event", "CMSPlus", Color=c('#E69F24','#0273B3','#CC79A7', '#3C5488', '#8491B4'), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
