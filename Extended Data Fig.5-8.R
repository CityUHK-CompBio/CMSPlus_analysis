load("/mnt/wutan/data/9.CMSPlus/FigureS5/C_ARGO_clinical_SDI.rdata")
# survival
CMScolor <- c('#E69F24','#0273B3', '#CC79A7', '#009F73')
CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4')
CMSSDIcolor =c('#E69F24','#0273B3', '#CC79A7', '#C2CB1E', '#537A34')
ARGO_sdiused_label[ARGO_sdiused_label$sample == "14TS0Z0425TR1", "dfs.event"] <- 0
ARGO_sdiused_label[ARGO_sdiused_label$sample == "14TS0Z0425TR1", "rfs.event"] <- 0
ARGO_sdiused_label[ARGO_sdiused_label$rfs.delay > 80, "rfs.delay"] <- 80
ARGO_sdiused_label[ARGO_sdiused_label$rfs.delay > 80, "rfs.event"] <- 0
ARGO_sdiused_label[ARGO_sdiused_label$dfs.delay > 80, "dfs.delay"] <- 80
ARGO_sdiused_label[ARGO_sdiused_label$dfs.delay > 80, "dfs.event"] <- 0
#####################################################################
################### ARGO
ARGO_sdiused_label$rfs.event[grep("复发", ARGO_sdiused_label$生存现况)] <- 1
ARGOcutoff <- median(ARGO_sdiused_label$Gau_innerAverSDI)
#[1] 0.737966
#ARGOcutoff <- 0.7476921
ARGO_sdiused_label$SDI <- "HighSDI"
ARGO_sdiused_label$SDI[ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "LowSDI"
## CMS1
ARGO_sdiused_label$CMS1.SDI <- as.character(ARGO_sdiused_label$cms)
ARGO_sdiused_label$CMS1.SDI[ARGO_sdiused_label$cms == "CMS1" & ARGO_sdiused_label$Gau_innerAverSDI > ARGOcutoff] <- "CMS1_HighSDI"
ARGO_sdiused_label$CMS1.SDI[ARGO_sdiused_label$cms == "CMS1" & ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "CMS1_LowSDI"
## CMS2
ARGO_sdiused_label$CMS2.SDI <- as.character(ARGO_sdiused_label$cms)
ARGO_sdiused_label$CMS2.SDI[ARGO_sdiused_label$cms == "CMS2" & ARGO_sdiused_label$Gau_innerAverSDI > ARGOcutoff] <- "CMS2_HighSDI"
ARGO_sdiused_label$CMS2.SDI[ARGO_sdiused_label$cms == "CMS2" & ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "CMS2_LowSDI"
## CMS3
ARGO_sdiused_label$CMS3.SDI <- as.character(ARGO_sdiused_label$cms)
ARGO_sdiused_label$CMS3.SDI[ARGO_sdiused_label$cms == "CMS3" & ARGO_sdiused_label$Gau_innerAverSDI > ARGOcutoff] <- "CMS3_HighSDI"
ARGO_sdiused_label$CMS3.SDI[ARGO_sdiused_label$cms == "CMS3" & ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "CMS3_LowSDI"
## CMS4

####### two groups
### SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "rfs.delay", "rfs.event", "SDI", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "dfs.delay", "dfs.event", "SDI", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS1.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS1", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "CMS1.SDI", Color=c("#FCD7A2", "#DE9B2A"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS1", ARGO_sdiused_label$cms), ], "dfs.delay", "dfs.event", "CMS1.SDI", Color=c("#FCD7A2", "#DE9B2A"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS2.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS2", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "CMS2.SDI", Color=c("#A6D2E8", "#106DA8"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS2", ARGO_sdiused_label$cms), ], "dfs.delay", "dfs.event", "CMS2.SDI", Color=c("#A6D2E8", "#106DA8"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS3.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS3", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "CMS3.SDI", Color=c("#E8B5D3", "#C2759F"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS3", ARGO_sdiused_label$cms), ], "dfs.delay", "dfs.event", "CMS3.SDI", Color=c("#E8B5D3", "#C2759F"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS4.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS4", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "self_all_inAverSDI_5", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS4", ARGO_sdiused_label$cms), ], "dfs.delay", "dfs.event", "self_all_inAverSDI_5", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
######## more than two groups
### CMS-SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "rfs.delay", "rfs.event", "self_all_inAverSDI_5", Color=CMSSDIcolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "dfs.delay", "dfs.event", "self_all_inAverSDI_5", Color=CMSSDIcolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "rfs.delay", "rfs.event", "cms", Color=CMScolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "dfs.delay", "dfs.event", "cms", Color=CMScolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

####### SDI
load("/mnt/wutan/data/9.CMSPlus/FigureS5/ARGO_sdiused_label2.v5.RData")
ARGO_sdiused_label2[ARGO_sdiused_label2$RNAseqID == "14TS0Z0425TR1", "rfs.event"] <- 0
ARGO_sdiused_label2[ARGO_sdiused_label2$RNAseqID == "14TS0Z0425TR1", "dfs.event"] <- 0
ARGO_sdiused_label2[ARGO_sdiused_label2$rfs.delay > 80, "rfs.delay"] <- 80
ARGO_sdiused_label2[ARGO_sdiused_label2$rfs.delay > 80, "rfs.event"] <- 0
ARGO_sdiused_label2[ARGO_sdiused_label2$dfs.delay > 80, "dfs.delay"] <- 80
ARGO_sdiused_label2[ARGO_sdiused_label2$dfs.delay > 80, "dfs.event"] <- 0
ARGOcutoff <- median(ARGO_sdiused_label2$Gau_innerAverSDI)
ARGO_sdiused_label2$SDI <- "SDI+"
ARGO_sdiused_label2$SDI[ARGO_sdiused_label2$Gau_innerAverSDI <= ARGOcutoff] <- "SDI-"
ARGOSDIHigh <- ARGO_sdiused_label2[ARGO_sdiused_label2$SDI == "SDI+" & !is.na(ARGO_sdiused_label2$chemotherapy), ]  ## 194
ARGOSDILow <- ARGO_sdiused_label2[ARGO_sdiused_label2$SDI == "SDI-" & !is.na(ARGO_sdiused_label2$chemotherapy), ]  ## 186
ARGOSDIHigh$chemotherapy <- factor(ARGOSDIHigh$chemotherapy, levels=c("Surgery alone", "Adjuvant chemotherapy"))
ARGOSDILow$chemotherapy <- factor(ARGOSDILow$chemotherapy, levels=c("Surgery alone", "Adjuvant chemotherapy"))
cmssdirfs <- SurvivalFGPlotSimplify(ARGOSDIHigh, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#C2CB1E", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGOSDIHigh, "dfs.delay", "dfs.event", "chemotherapy", Color=c("#C2CB1E", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGOSDILow, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#537A34", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGOSDILow, "dfs.delay", "dfs.event", "chemotherapy", Color=c("#537A34", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
