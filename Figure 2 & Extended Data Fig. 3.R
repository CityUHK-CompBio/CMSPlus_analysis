####### Figure2A
options(stringsAsFactors=FALSE)
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/3.differential/exp")
library(HTSanalyzeR2)
load("CMS.DEIF.GSEA.result.RData")
DEIS.GSEA.rslt <- GSEA.rslt
DEIS.heat.dat.p <- CCSEnrichmentHeatmapData(DEIS.GSEA.rslt, P="Pvalue", terms=TRUE)
DEIS.heat.dat.p$Genesetname <- sapply(strsplit(rownames(DEIS.heat.dat.p), "\\."), function(x) x[1])
labelcolors <- list(Subtype = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4', 'CMS4-TME- vs CMS4-TME+'='#1E8BA9'))
aaa <- DEIS.heat.dat.p[, 1:6]
rownames(aaa) <- gsub("CMS.gs.", "", rownames(aaa))
CCSEnrichmentHeatmapRename <- function(heat.dat, labelcolors, cancerintrinsic=FALSE, colorheatmap=colorRampPalette(c("#4A7F4B", "#FFFFFF", "#B41600"))(100), breaks = c(seq(-4, 4, length.out = 100))){
  Signature.heat <- heat.dat[match(c("EPITH_LOBODA",
                                     "WNT_FLIER",
                                     "MYC_TARGETS_ZELLER",
                                     "MESENCH_LOBODA",
                                     "EMT_CORE_GENES",
                                     "TGFB_KEGG",
                                     "MATRIX_REMODEL_REACTOME",
                                     "WOUND_RESPONSE_GO_BP",
                                     "CSC_BATLLE"), rownames(heat.dat)), ]
  
  rownames(Signature.heat) <- c("Epithelial", "WNT targets", "MYC targets", "Mesenchymal",
                                "EMT activation", "TGF-beta activation",
                                "Matrix remodeling", "Wound response", "Cancer stem cell")
  ## Pathways.heat
  Pathways.heat <- heat.dat[match(c("MAPK_KEGG",
                                    "PI3K_ACT_REACTOME",
                                    "SRC_ACT_BILD",
                                    "JAK_STAT_KEGG",
                                    "CASPASE_BIOCARTA",
                                    "PROTEASOME_KEGG",
                                    "KEGG_CELL_CYCLE",
                                    "TRANSLATION_RIBOS_REACTOME",
                                    "INTEGRIN_BETA3_CP",
                                    "VEGF_VEGFR_REACTOME"), rownames(heat.dat)), ]
  rownames(Pathways.heat) <- c("MAPK", "PI3K", "SRC",
                               "JAK-STAT", "Caspases", "Proteosome", "Cell cycle",
                               "Translation ribosome", "Integrin-beta3", "VEGF VEGFR")
  
  ## Estimate.heat
  Estimate.heat <- heat.dat[match(c("IMMUNE_ESTIMATE",
                                    "STROMAL_ESTIMATE"), rownames(heat.dat)), ]
  rownames(Estimate.heat) <- c("Immune infiltration", "Stromal infiltration")
  
  ## Immune.heat
  Immune.heat <- heat.dat[match(c("IMMUNE_RESP_GO_BP",
                                  "PD1_REACTOME",
                                  "IMMUNE_NKC_BREAST",
                                  "IMMUNE_TH1_GALON",
                                  "IMMUNE_THF_BREAST",
                                  "IMMUNE_TH17_GOUNARI",
                                  "IMMUNE_TREG_GALON",
                                  "COMPLEMENT_COAG_KEGG"), rownames(heat.dat)), ]
  
  rownames(Immune.heat) <- c("Immune response","PD1 activation", "NK cell infiltration",
                             "TH1 infiltration","TFH infiltration", "TH17 activation",
                             "Treg activation", "Complement activation")
  
  ## Metabolism.heat
  Metabolism.heat <- heat.dat[match(c("AMINO_SUGAR_NUCLEO_METAB_KEGG",
                                      "PENTOSE_GLUC_METAB_KEGG",
                                      "FRUTOSE_MANNOSE_METAB_KEGG",
                                      "GALACTOSE_METAB_KEGG",
                                      "GLUTAMINE_GO_BP",
                                      "GLUTATHIONE_KEGG",
                                      "NITROGEN_METAB_KEGG",
                                      "GLYCEROPHOSPHOLIPID_GO_BP",
                                      "LYSOPHOSPHOLIPID_PID",
                                      "FATTY_ACID_METAB_KEGG"), rownames(heat.dat)), ]
  rownames(Metabolism.heat) <- c("Sugar aa nucleotide", "Glucose pentose",
                                 "Fructose mannose", "Galactose", "Glutamine", "Glutathione",
                                 "Nitrogen", "Glycerophospholipid", "Lysophospholipid", "Fatty acid")
  if(cancerintrinsic == TRUE){
    Intrinsic.heat <- heat.dat[match(c("CRIS", "Eschrich", "Jorissen", "Kennedy", "PDAC", "Popovici", "Sadanandam", "The_30_gene"), rownames(heat.dat)), ]
  }
  ## pheatmap
  library(pheatmap)
  col.annotation <- data.frame("group" = colnames(heat.dat),
                               "Subtype" = colnames(heat.dat))
  annotation_col <- data.frame("Subtype" = as.factor(col.annotation$Subtype))
  rownames(annotation_col) <- colnames(Metabolism.heat)
  ann_colors = labelcolors
  ## Signature.heatmap
  #breaks <- c(seq(-4, 4, length.out = 100))
  #breaks <- c(seq(-4, -1.5, length.out = 10), seq(-1.49, 0, length.out = 40), seq(0.01, 1.5, length.out = 40), seq(1.51, 4, length.out = 10))
  Signature.p <- pheatmap(Signature.heat,
                          filename = "Signature.pdf",
                          fontsize =10, breaks = breaks,
                          show_rownames=T,annotation_col = annotation_col, annotation_colors = ann_colors,
                          color = colorheatmap,
                          # border_color = F,
                          annotation_names_col=F,show_colnames=F,
                          cluster_cols =F,cluster_rows=F,fontsize_col=12,
                          cellwidth = 15,cellheight=15,
                          main = "Signatures")
  ## Pathways.heatmap
  Pathways.p <- pheatmap(Pathways.heat,
                         filename = "Pathways.pdf",
                         fontsize =10, breaks = breaks,
                         show_rownames=T,annotation_col = annotation_col, annotation_colors = ann_colors,
                         color = colorheatmap,
                         # border_color = F,
                         annotation_names_col=F,show_colnames=F,
                         cluster_cols =F,cluster_rows=F,fontsize_col=12,
                         cellwidth = 15,cellheight=15,
                         legend = F, annotation_legend = F,
                         main = "Pathways")
  ## Estimate.heatmap
  Estimate.p <- pheatmap(Estimate.heat,
                         filename = "Estimate.pdf",
                         fontsize =10, breaks = breaks,
                         show_rownames=T,annotation_col = annotation_col, annotation_colors = ann_colors,
                         color = colorheatmap,
                         # border_color = F,
                         annotation_names_col=F,show_colnames=F,
                         cluster_cols =F,cluster_rows=F,fontsize_col=12,
                         cellwidth = 15,cellheight=15,
                         legend = F, annotation_legend = F,
                         main = "Estimate")
  ## Immune.heatmap
  Immune.p <- pheatmap(Immune.heat,
                       filename = "Immune.pdf",
                       fontsize =10, breaks = breaks,
                       show_rownames=T,annotation_col = annotation_col, annotation_colors = ann_colors,
                       color = colorheatmap,
                       # border_color = F,
                       annotation_names_col=F,show_colnames=F,
                       cluster_cols =F,cluster_rows=F,fontsize_col=12,
                       cellwidth = 15,cellheight=15,
                       legend = F, annotation_legend = F,
                       main = "Immune")
  
  ## Metabolism.heatmap
  # breaks <- c(seq(-5, 5, length.out = 100))
  Metabolism.p <- pheatmap(Metabolism.heat,
                           filename = "Metabolism.pdf",
                           fontsize =10, breaks = breaks,
                           show_rownames=T,annotation_col = annotation_col, annotation_colors = ann_colors,
                           color = colorheatmap,
                           # border_color = F,
                           annotation_names_col=F,show_colnames=F,
                           cluster_cols =F,cluster_rows=F,fontsize_col=12,
                           cellwidth = 15,cellheight=15,
                           legend = F, annotation_legend = F,
                           main = "Metabolism")
  save(Signature.p, Pathways.p, Immune.p, Metabolism.p, Estimate.p,
       file = "Figure.RData")
  ## cancer intrinsic
  if(cancerintrinsic == TRUE){
    Intrinsic.p <- pheatmap(Intrinsic.heat,
                            filename = "Intrinsic.pdf",
                            fontsize =10, breaks = breaks,
                            show_rownames=T,annotation_col = annotation_col, annotation_colors = ann_colors,
                            color = colorheatmap,
                            # border_color = F,
                            annotation_names_col=F,show_colnames=F,
                            cluster_cols =F,cluster_rows=F,fontsize_col=12,
                            cellwidth = 15,cellheight=15,
                            legend = F, annotation_legend = F,
                            main = "Cancer signatures")
    save(Signature.p, Pathways.p, Immune.p, Metabolism.p, Estimate.p, Intrinsic.p,
         file = "Figure.RData")
  }
}
CCSEnrichmentHeatmapRename(aaa, labelcolors=labelcolors, colorheatmap=colorRampPalette(c("#532788", "#FFFFFF", "#B35807"))(100), cancerintrinsic=TRUE, breaks = c(seq(-2, 2, length.out = 100)))
load("Figure.RData")
pdf("pathway.pdf", height = 9.1, width = 10)
	r1 <- plot_grid(Signature.p$gtable, Pathways.p$gtable, labels = NULL,
	                nrow = 1, align = "h")
	r2 <- plot_grid(Immune.p$gtable, Estimate.p$gtable, labels = NULL,
	                nrow = 2, align = "v")
	r3 <- plot_grid(Metabolism.p$gtable, r2, labels = NULL,
	                nrow = 1, align = "h")
	plot_grid(r1, r3, labels = NULL,
	          nrow = 2, align = "none")
	dev.off()


####### Figure2B
library(maftools)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
TCGAclindata <- CRC.clin.exp.he_v2
TCGAclindata$Tumor_Sample_Barcode <- CRC.clin.exp.he_v2$sample
CRCMaf <- read.maf(maf = "/data/home2/wutan/Task/9.CMS.SDI/Procedures/StartFrom20210310/CRC_maftools.maf", clinicalData = TCGAclindata, verbose = FALSE, isTCGA = TRUE)
CRCMaf.subtype = clinicalEnrichment(maf = CRCMaf, clinicalFeature = 'CMSPlus')
CMSPlusSubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
MutationCountSample <- table(CRCMaf@data$Tumor_Sample_Barcode)
MutationCountDF <- data.frame(Sample=names(MutationCountSample), MutationCount=as.numeric(MutationCountSample), Subtype=TCGAclindata$CMSPlus[match(names(MutationCountSample), TCGAclindata$Tumor_Sample_Barcode)])
##
t.test(MutationCountDF$MutationCount[which(MutationCountDF$Subtype %in% "CMS4-TME+")], MutationCountDF$MutationCount[which(MutationCountDF$Subtype %in% "CMS4-TME-")])
# p-value = 0.01493
my_comparisons = list( c("CMS4-TME+", "CMS4-TME-") )
MutationCountDF <- MutationCountDF[!is.na(MutationCountDF$Subtype), ]
bp1 <- ggplot(MutationCountDF, aes(x=Subtype, y=MutationCount, fill=Subtype)) + geom_boxplot(width=0.8) + 
      labs(x="Subtype", y = "Mutation Count") + 
      scale_fill_manual(values=CMSPlusSubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(0, 1500)) +
      geom_signif(comparisons=my_comparisons, annotations=c("1.49e-02"), y_position = c(700), tip_length = 0, vjust=0)


####### Figure2C
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


####### Figure2D
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")
# age
CMSclin2$Age <- "50-"
CMSclin2$Age[CMSclin2$age > 50 & CMSclin2$age<= 65] <- "51-65"
CMSclin2$Age[CMSclin2$age > 65 & CMSclin2$age <= 75] <- "66-75"
CMSclin2$Age[CMSclin2$age > 75] <- "76+"

fisher.test(table(CMSclin2$CMSTME, CMSclin2$Age), workspace = 2e7, simulate.p.value=TRUE)
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$Age), workspace = 2e7, simulate.p.value=TRUE)
SDI1.1_Age <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="Age", Sample="sample", adjust=FALSE)
age1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$Age))
pp1 <- ggplot(age1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("#1E8BA9", "#FECB24", "#DA1F27", "#704122")) + 
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Age") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) 

# Gender 

fisher.test(table(CMSclin2$CMSTME, CMSclin2$gender), workspace = 2e7, simulate.p.value=TRUE)
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$gender), workspace = 2e7, simulate.p.value=TRUE)
SDI1.1_Gender <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="gender", Sample="sample", adjust=FALSE)
Gender1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$gender))
pp2 <- ggplot(Gender1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Gender") + ggtitle("p-value = 0.005") + theme(axis.text.x=element_text(angle=60,hjust=1))

# Stage

fisher.test(table(CMSclin2$CMSTME, CMSclin2$stage), workspace = 2e7, simulate.p.value=TRUE)
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$stage), workspace = 2e7, simulate.p.value=TRUE)
SDI1.1_TumorLocation1 <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="stage", Sample="sample", adjust=FALSE)
TumorLocation11.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$stage))
pp4 <- ggplot(TumorLocation11.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("1"=RColorBrewer::brewer.pal(9, "YlGn")[2], "2"=RColorBrewer::brewer.pal(9, "YlGn")[4], "3"=RColorBrewer::brewer.pal(9, "YlGn")[6], "4"=RColorBrewer::brewer.pal(9, "YlGn")[8])) + 
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Stage") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) 

# tnm
CMSclin2$TNM <- CMSclin2$tnm
CMSclin2$TNM <- gsub("A", "", CMSclin2$TNM)
CMSclin2$TNM <- gsub("B", "", CMSclin2$TNM)
CMSclin2$TNM <- gsub("C", "", CMSclin2$TNM)
CMSclin2$TNM <- gsub("/II", "", CMSclin2$TNM)
table(CMSclin2$TNM, useNA="always")

# pn
fisher.test(table(CMSclin2$CMSTME, CMSclin2$pn), workspace = 2e7, simulate.p.value=TRUE)
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$pn), workspace = 2e7, simulate.p.value=TRUE)
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="pn", Sample="sample", adjust=FALSE)
TNM1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$pn))
ppn1 <- ggplot(TNM1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(10, "Paired")[2])) + 
		xlab("Subtype") + ylab("Percentage") + labs(fill = "pn") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) 

# pm
fisher.test(table(CMSclin2$CMSTME, CMSclin2$pm), workspace = 2e7, simulate.p.value=TRUE)
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$pm), workspace = 2e7, simulate.p.value=TRUE)
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="pm", Sample="sample", adjust=FALSE)
TNM1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$pm))
ppm1 <- ggplot(TNM1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(12, "Paired")[12])) + 
		xlab("Subtype") + ylab("Percentage") + labs(fill = "pm") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1))

# Location

load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
ccc <- CRC.clin.exp.he_v2
colnames(ccc)[colnames(ccc) %in% "Immune.Subtype"] <- "Subtype"
colnames(ccc)[colnames(ccc) %in% "tumor.location"] <- "TumorLocation"

fisher.test(table(ccc$Subtype, ccc$TumorLocation), workspace = 2e7, simulate.p.value=TRUE)
ccc_CMS4 <- ccc[grep("CMS4", ccc$Subtype), ]
fisher.test(table(ccc_CMS4$Subtype, ccc_CMS4$TumorLocation), workspace = 2e7, simulate.p.value=TRUE)
SDI1.1_TumorLocation1 <- HyperGeoForTwoGroups(T90Clin=ccc, GroupName1="Subtype", GroupName2="TumorLocation", Sample="sample", adjust=FALSE)
TumorLocation11.1box <- as.data.frame(table(ccc$Subtype, ccc$TumorLocation))
pp5 <- ggplot(TumorLocation11.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Dark2")[3:4]) + 
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Tumor Location") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) 

figure <- ggarrange(pp1, pp2, pp4, pp5, ppn1, ppm1,
                    labels = c("A", "B", "C", "D", "E", "F", "G"),
                    ncol = 4, nrow = 2)
figure

####### Figure2E & FigureS3

load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/5.otherclassifier")
CRIStcga <- read.table("CRIS.label.txt", header=T, sep="\t", quote="", fill=T)  ## download from paper
CRC.clin.exp.he_v2$CRIS.TSP <- CRIStcga$CRIS.TSP.Class.assignment[match(CRC.clin.exp.he_v2$sample, CRIStcga$Sample.ID)]
CRC.clin.exp.he_v2$CRIS.NTP80 <- CRIStcga$CRIS.NTP80.Class.Assignment[match(CRC.clin.exp.he_v2$sample, CRIStcga$Sample.ID)]
clin3 <- CRC.clin.exp.he_v2[!is.na(CRC.clin.exp.he_v2$CRIS.TSP) & CRC.clin.exp.he_v2$CRIS.TSP != "", ]

### FigureS3
library(survival)
library(gfplots)
source("~/R/functions/functions.R")
cmssdirfs <- SurvivalFGPlotSimplify(clin3, "rfs.delay", "rfs.event", "CRIS.TSP", Color=c("#EF4C2A", "#D42127", "#1C275C", "#018647", "#00AD9B"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

### Figure2E
GetHGPlot <- function(atable, colname1, colname2, samplename){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 0.05, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow=length((table(atable[, colname1]))))
	bk = c(seq(0, 0.05, length.out=20), seq(0.051, 1, length.out=20))
	my_palette <- c(rev(colorpanel(n=20, low="#BFD8EA",high="#F0F3FA")), rev(colorpanel(n=20, low="#02519D",high="#BFD8EA")))
	pheatmap(SIAGCadj, scale="none", cluster_row = FALSE, cluster_col = FALSE, color = rev(my_palette), breaks = bk, display_numbers = SIAGCshow, fontsize =13, cellwidth = 30, cellheight=30, gaps_col=c(1, 2, 3, 4), gaps_row=c(1, 2, 3, 4))
}
GetHGPlot(atable=clin3, colname1="CMSPlus", colname2="CRIS.TSP", samplename="sample")


####### Figure2F & FigureS3
load("Normalized_Count_TCGA622.RData")
options(stringsAsFactors  = F)
suppressMessages(library(tidyverse))
suppressMessages(library(getopt))
suppressMessages(library(caret))
suppressMessages(library(xgboost))
suppressMessages(library(tibble))
args <- list(input.data= Normalized_Count_TCGA622,   #expression profile file
          output.data = "CCCRC_output.txt", #output file
          log2 = "F",        # Whether the data needs to be Log2 processed:T/F,Default is T
          scale = "T",      # character","Whether the data needs to be scaled:T/F,Default is T
          xgboost.model = "CCCRC_classifier-model.Rdata", # model file  ：CCCRC_classifier-model.Rdata
          genelist = "CCCRC_classifier-genelist.Rdata")    # gene list：CCCRC_classifier-genelist.Rdata
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
cccrctcga <- read.table("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/5.otherclassifier/CCCRC_output.txt", header=T, sep="\t", quote="")
identical(CRC.clin.exp.he_v2$sample, cccrctcga$id)
CRC.clin.exp.he_v2$CCCRC <- cccrctcga$CCCRC_predict

### FigureS3
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "rfs.delay", "rfs.event", "CCCRC", Color=c("#72598D", "#FB9B03", "#FB0206", "#0765F7"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "rfs.delay", "rfs.event", "RF.nearestCMS.mod", Color=c('#E69F24','#0273B3', '#CC79A7', '#009F73'), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

### Figure2F
GetHGPlot <- function(atable, colname1, colname2, samplename){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 0.05, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow=length((table(atable[, colname1]))))
	bk = c(seq(0, 0.05, length.out=20), seq(0.051, 1, length.out=20))
	my_palette <- c(rev(colorpanel(n=20, low="#BFD8EA",high="#F0F3FA")), rev(colorpanel(n=20, low="#02519D",high="#BFD8EA")))
	pheatmap(SIAGCadj, scale="none", cluster_row = FALSE, cluster_col = FALSE, color = rev(my_palette), breaks = bk, display_numbers = SIAGCshow, fontsize =13, cellwidth = 30, cellheight=30, gaps_col=c(1, 2, 3, 4), gaps_row=c(1, 2, 3, 4))
}
GetHGPlot(atable=CRC.clin.exp.he_v2, colname1="CMSPlus", colname2="CCCRC", samplename="sample")

####### Figure2GH
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/5.otherclassifier")
CCCRCids <- read.table("CCCRC_offical_label.txt", header=T, sep="\t", quote="")
CRISids <- read.table("CRIS.label.txt", header=T, sep="\t", quote="", fill=T)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")  ## CMSclin2, gsvaall2
CCCRCids$id <- gsub("\\.", "-", CCCRCids$id)
length(intersect(CMSclin2$sample, CCCRCids$id))  ## 1244
length(intersect(CMSclin2$sample, CRISids$Sample.ID)) ## 1351
SampleInt <- Reduce(intersect, list(CMSclin2$sample, CCCRCids$id, CRISids$Sample.ID))  ## 1239
CMSclin3 <- CMSclin2[CMSclin2$sample %in% SampleInt, ]
CMSclin3$CCCRC <- CCCRCids$CCCRC[match(CMSclin3$sample, CCCRCids$id)]
CMSclin3$CRIS.TSP <- CRISids$CRIS.TSP.Class.assignment[match(CMSclin3$sample, CRISids$Sample.ID)]
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.all.45.v2.RData")
gsva.all.45.clin3 <- gsva.all.45[, CMSclin3$sample]
col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))

### Figure2H
CMSclin32 <- CMSclin3[order(CMSclin3$CCCRC), ]
CCCRCha = HeatmapAnnotation(
    CCCRC = CMSclin32$CCCRC, 
    MSI = CMSclin32$msi,
    CIMP = CMSclin32$cimp,
    KRAS = CMSclin32$kras_mut,
    BRAF = CMSclin32$braf_mut,
    col = list(CMSplus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', `CMS4-TME+`='#3C5488', `CMS4-TME-`='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), CCCRC =c(C1="#72598D", C2="#FB9B03", C3="#FB0206", C4="#0765F7"), CRIS =c(Nolabel="#DEDEDE", CRISA="#EF4C2A", CRISB="#D42127", CRISC="#1C275C", CRISD="#018647", CRISE="#00AD9B"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20')),
    na_col = "white"
)
CCCRChaP <- Heatmap(gsva.all.45.clin3[43:45, CMSclin32$sample], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = CCCRCha, column_labels=rep("", 1239))
#draw(CCCRCha)

### Figure2G
CMSclin33 <- CMSclin3[order(CMSclin3$CRIS.TSP), ]
CMSclin33 <- CMSclin33[CMSclin33$CRIS.TSP != "" & !is.na(CMSclin33$CRIS.TSP), ]
CRISha = HeatmapAnnotation(
    CRIS = CMSclin33$CRIS.TSP, 
    MSI = CMSclin33$msi,
    CIMP = CMSclin33$cimp,
    KRAS = CMSclin33$kras_mut,
    BRAF = CMSclin33$braf_mut,
    col = list(CMSplus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', `CMS4-TME+`='#3C5488', `CMS4-TME-`='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), CCCRC =c(C1="#72598D", C2="#FB9B03", C3="#FB0206", C4="#0765F7"), CRIS =c(Nolabel="#DEDEDE", CRISA="#EF4C2A", CRISB="#D42127", CRISC="#1C275C", CRISD="#018647", CRISE="#00AD9B"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20')),
    na_col = "white"
)
CRIShaP <- Heatmap(gsva.all.45.clin3[43:45, CMSclin33$sample], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = CRISha, column_labels=rep("", 1115))
#draw(CRISha)

