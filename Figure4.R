### SYSU
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
CMSPlusColor=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', `CMS4-TME+`='#3C5488', `CMS4-TME-`='#8491B4')
cmsclin.m <- cmsclin[!is.na(cmsclin$DFS_days), ]
cmsclin.m[cmsclin.m$Sample_RNA == "14TS0Z0425TR1", "Relapse_Death"] <- 0
cmsclin.m[cmsclin.m$DFS_days > 2400, "DFS_days"] <- 2400
cmsclin.m[cmsclin.m$DFS_days > 2400, "Relapse_Death"] <- 0

## Figure4a-b
MidStageSamples <- cmsclin.m
MidStageSamplesCMS1 <- MidStageSamples[MidStageSamples$CMSPlus2 == "CMS1" & !is.na(MidStageSamples$chemotherapy), ]  ## 81 32
MidStageSamplesCMS2 <- MidStageSamples[MidStageSamples$CMSPlus2 == "CMS2" & !is.na(MidStageSamples$chemotherapy), ]  ## 152
MidStageSamplesCMS3 <- MidStageSamples[MidStageSamples$CMSPlus2 == "CMS3" & !is.na(MidStageSamples$chemotherapy), ]  ## 54
MidStageSamplesCMS4_HighIF <- MidStageSamples[MidStageSamples$CMSPlus2 == "CMS4-TME+" & !is.na(MidStageSamples$chemotherapy), ]  ## 74
MidStageSamplesCMS4_LowIF <- MidStageSamples[MidStageSamples$CMSPlus2 == "CMS4-TME-" & !is.na(MidStageSamples$chemotherapy), ]  ## 61
MidStageSamplesCMS4 <- MidStageSamples[intersect(grep("CMS4", MidStageSamples$CMSPlus2), which(!is.na(MidStageSamples$chemotherapy))), ]  ## 122 114

### CMS1
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS1, "DFS_days", "Relapse_Death", "chemotherapy", Color=c("#E69F24", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS2
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS2, "DFS_days", "Relapse_Death", "chemotherapy", Color=c("#0273B3", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS3
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS3, "DFS_days", "Relapse_Death", "chemotherapy", Color=c("#CC79A7", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS4
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS4_HighIF, "DFS_days", "Relapse_Death", "chemotherapy", Color=c("#3C5488", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS4_LowIF, "DFS_days", "Relapse_Death", "chemotherapy", Color=c("#8491B4", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS4, "DFS_days", "Relapse_Death", "chemotherapy", Color=c("#009F73", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

## Figure4c
load("/data1/wutan/14.amc/3.TimFFPE/3.Expression/PREDICATECount.RData")
library(GSVA)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/CMSsignatures42.v2.RData")
gsva.predicate <- gsva(as.matrix(PREDICATENCount), CMSsignatures42.v2, method=c("gsva"), kcdf = "Gaussian", min.sz = 1, parallel.sz=40)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/epitme.rf.model.Rdata")
predicate.vals <- as.data.frame(t(gsva.predicate))
colnames(predicate.vals) <- gsub("-", "_", colnames(predicate.vals))
colnames(predicate.vals) <- gsub(" ", "_", colnames(predicate.vals))
probabilitiespredicate = predict(epitme.rf.model, data = predicate.vals)$predictions
table(probabilitiespredicate )
#probabilitiespredicate
#    CMS1     CMS2     CMS3 CMS4-IF+ CMS4-IF- 
#      12       58       16       25       18
library(haven)
clindata3 <- read_sav('/data0/tan/Task/14.amc/3.TimFFPE/1.rawdata/infor/Data NKR levering6 PREDICATE only SPSS-003.sav')
clindata3 <- as.data.frame(clindata3)
clindata3$org_id <- gsub("R", "", clindata3$org_id)
clindata3$CMSTME <- probabilitiespredicate[match(clindata3$org_id, colnames(gsva.predicate))]
clindata3$CMS <- substr(clindata3$CMSTME, 1, 4)
clindata3 <- clindata3[!is.na(clindata3$CMSTME), ]
CMSTMESubtype=c('#E69F24','#0273B3', '#CC79A7', '#3C5488', '#8491B4')
clindata3.cms4l <- clindata3[clindata3$CMSTME=="CMS4-TME-" , ]
cmssdirfs <- SurvivalFGPlotSimplify(clindata3.cms4l, "vit_stat_int", "vit_stat", "chemo_binair", c("#8491B4", "#6D3E33"), type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
							
## Figure4d
CMSPlusColor=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
PD1GeneSig <- c("CD2", "CTSW", "GZMB", "HLA-E", "HSPA8", "ID2", "ITM2B", "LSP1", "MS4A4A", "MS4A4E", "S100A10", "SERINC3", "SERPINE2", "SHISA")
cpexp2 <- t(scale(t(SYSU.exp[rownames(SYSU.exp) %in% PD1GeneSig, ])))
cpexp2.score <- apply(cpexp2, 2, mean)
cpexp2.df <- data.frame(score=cpexp2.score, subtype=SYSU.clin$CMSPlusLabel)
my_comparisons = list( c("CMS4-TME+", "CMS1"), c("CMS4-TME+", "CMS2"), c("CMS4-TME+", "CMS3"), c("CMS4-TME+", "CMS4-TME-"))
bpis <- ggplot(cpexp2.df, aes(x=subtype, y=score, fill=subtype)) + geom_boxplot(width=0.8) + 
      labs(x="Subtype", y = "PD1 response signature score") + 
      scale_fill_manual(values=CMSPlusColor) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)
#
CMSPlusColor=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
CD8GeneSig <- signature_tme$CD_8_T_effector
cpexp3 <- t(scale(t(SYSU.exp[rownames(SYSU.exp) %in% CD8GeneSig, ])))
cpexp3.score <- apply(cpexp3, 2, mean)
cpexp3.df <- data.frame(score=cpexp3.score, subtype=SYSU.clin$CMSPlusLabel)
my_comparisons = list( c("CMS4-TME+", "CMS1"), c("CMS4-TME+", "CMS2"), c("CMS4-TME+", "CMS3"), c("CMS4-TME+", "CMS4-TME-"))
bpis3 <- ggplot(cpexp3.df, aes(x=subtype, y=score, fill=subtype)) + geom_boxplot(width=0.8) + 
    labs(x="Subtype", y = "CD8 T cell response signature score") + 
    scale_fill_manual(values=CMSPlusColor) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)

## Figure4e
load("/mnt/wutan/data/9.CMSPlus/Figure3/SYSU.clin.CMSPlus.RData")
SYSU.clin$CMSPlusLabel <- factor(SYSU.clin$CMSPlusLabel, levels=c("CMS1", "CMS2", "CMS3", "CMS4-TME+", "CMS4-TME-"))
library(IOBR)
PD1Gene <- c("PDCD1", "CD274", "PDCD1LG2")
PD1GeneSig <- c("CD2", "CTSW", "GZMB", "HLA-E", "ITM2B", "LSP1", "MS4A4A", "SERPINE2")
cpexp3 <- SYSU.exp[rownames(SYSU.exp) %in% PD1GeneSig, ] ##
cpexp4 <- SYSU.exp[rownames(SYSU.exp) %in% signature_tme$CD_8_T_effector, ]
cpexp34 <- t(scale(t(rbind(cpexp3, cpexp4))))
cpexp35 <- t(scale(t(rbind(cpexp4, rbind(cpexp3, SYSU.exp[rownames(SYSU.exp) %in% PD1Gene, ])))))
# heatmap
col_fun3 = colorRamp2(unique(c(seq(-1.5, 0, length.out=18), seq(0, 1.5, length.out=18))), colorRampPalette(c("#3660A5", "#000000", "#CB1A1E"))(35))
dend34 = cluster_within_group(cpexp34, SYSU.clin$CMSPlusLabel)
colannoA34 = HeatmapAnnotation(
		CMSPlus=SYSU.clin$CMSPlusLabel,
    col = list(CMSPlus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4'))
)
row_split34 = data.frame(rep(c("CD8 T effector", "PD1 Gene", "Check point"), c(dim(cpexp4)[1], dim(cpexp3)[1], 4)))
PGARGOIT <- Heatmap(cpexp35, col = col_fun3, cluster_rows = FALSE, cluster_columns = dend34, top_annotation = colannoA34, column_labels=rep("", dim(SYSU.clin)[1]), row_split = row_split34)
set.seed(1234)
columnOrder34 <- c(intersect(column_order(PGARGOIT), grep("CMS1", SYSU.clin$CMSPlusLabel)), intersect(sample(column_order(PGARGOIT)), grep("CMS2", SYSU.clin$CMSPlusLabel)), intersect(sample(column_order(PGARGOIT)), grep("CMS3", SYSU.clin$CMSPlusLabel)), intersect(sample(column_order(PGARGOIT)), grep("CMS4-TME\\+", SYSU.clin$CMSPlusLabel)), intersect(sample(column_order(PGARGOIT)), grep("CMS4-TME-", SYSU.clin$CMSPlusLabel)))
colannoA342 = HeatmapAnnotation(
		CMSPlus=SYSU.clin$CMSPlusLabel[columnOrder34],
    col = list(CMSPlus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4'))
)
col.subsections34 <- as.numeric(table(SYSU.clin$CMSPlusLabel))
col_split34 = data.frame(rep(names(table(SYSU.clin$CMSPlusLabel)), col.subsections34))
PGARGOIT2 <- Heatmap(cpexp35[, columnOrder34], col = col_fun3, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = colannoA342, column_labels=rep("", dim(SYSU.clin)[1]), row_split = row_split34, column_split = col_split34)
PGARGOIT2

