######### Figure3
library(CMSPlus)
library(GSVA)
library(CMSclassifier)
library(ComplexHeatmap)
library(circlize)
library(randomForest)
library(ranger)
library(parallel)
load("/mnt/wutan/data/9.CMSPlus/Figure3/SYSU.data.RData")
# Figure3D
res <- CMSPlus(exp2symbol=SYSU.exp, CMSlabel=SYSU.clin$Final_CMS)
SYSU.clin$CMSPlusLabel <- res$CMSPlusLabel$nearest[match(rownames(SYSU.clin), rownames(res$CMSPlusLabel))]
save(res, SYSU.exp, SYSU.clin, file="/mnt/wutan/data/9.CMSPlus/Figure3/SYSU.clin.CMSPlus.RData")
write.table(data.frame(Model="testing", Dataset="SYSU", Sample=rownames(SYSU.clin), `Clustered label`=SYSU.clin$CMSPlus, `Predicted label`=SYSU.clin$CMSPlusLabel), file="/mnt/wutan/data/9.CMSPlus/Figure3/SYSU.label.txt", col.names=T, row.names=T, sep="\t", quote=F)
# Figure3B
SYSU.clin$CMSPredict <- substr(SYSU.clin$CMSPlusLabel, 1, 4)
tab <- table(SYSU.clin$Final_CMS, SYSU.clin$CMSPredict)
tab
#       CMS1 CMS2 CMS3 CMS4
#  CMS1   73    2    3   10
#  CMS2    4  292   22   10
#  CMS3   13    8  113    5
#  CMS4    7    5    0  271
acc.SYSU <- sum(diag(tab)) / sum(tab)
# 0.8937947
# Figure3B
# python
from mlxtend.plotting import plot_confusion_matrix
import matplotlib.pyplot as plt
import numpy as np
multiclass = np.array([[73, 2, 3, 10],
                       [4, 292, 22, 10],
                       [13, 8, 113, 5],
                       [7, 5, 0, 271]])
class_names = ['CMS1', 'CMS2', 'CMS3', 'CMS4']
fig, ax = plot_confusion_matrix(conf_mat=multiclass,
                                colorbar=True,
                                show_absolute=True,
                                show_normed=True,
                                class_names=class_names)
plt.savefig("/mnt/wutan/data/9.CMSPlus/Figure3/Figure3B.pdf", format="pdf", bbox_inches="tight")
plt.show()

# Figure3F
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
SYSU.clin.m <- SYSU.clin[!is.na(SYSU.clin$DFS_days), ]
SYSU.clin.m[SYSU.clin.m$Sample_RNA == "14TS0Z0425TR1", "Relapse_Death"] <- 0
SYSU.clin.m[SYSU.clin.m$DFS_days > 2400, "DFS_days"] <- 2400
SYSU.clin.m[SYSU.clin.m$DFS_days > 2400, "Relapse_Death"] <- 0
SYSU.clin.m[SYSU.clin.m$rfs.delay > 80, "rfs.delay"] <- 80
SYSU.clin.m[SYSU.clin.m$rfs.delay > 80, "rfs.event"] <- 0

SYSU.clin.m$CMSPlusLabel <- factor(SYSU.clin.m$CMSPlusLabel, levels=names(CMSPlusColor))
cmssdirfs <- SurvivalFGPlotSimplify(SYSU.clin.m, "DFS_days", "Relapse_Death", "CMSPlusLabel", Color=CMSPlusColor, type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
# FigureS4
cmssdirfs <- SurvivalFGPlotSimplify(SYSU.clin.m, "rfs.delay", "rfs.event", "CMSPlusLabel", Color=CMSPlusColor, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
# CMS
CMSColor = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73")
cmssdirfs <- SurvivalFGPlotSimplify(SYSU.clin.m, "DFS_days", "Relapse_Death", "Final_CMS", Color=CMSColor, type="days")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(SYSU.clin.m, "rfs.delay", "rfs.event", "Final_CMS", Color=CMSColor, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

# CMSPlus clustering
SYSU.clin$CMSPredict <- substr(SYSU.clin$CMSPlusLabel, 1, 4)
SYSU.clin$CMSPlusLabel <- factor(SYSU.clin$CMSPlusLabel, levels=c(names(CMSPlusColor)))
SYSU.clin$CMSPlus <- factor(SYSU.clin$CMSPlus, levels=names(table(SYSU.clin$CMSPlus))[c(1:3, 5, 4)])
tab <- table(SYSU.clin$CMSPlus, SYSU.clin$CMSPlusLabel)
tab
#              CMS1 CMS2 CMS3 CMS4-TME+ CMS4-TME-
#  CMS1          73    2    3         3         7
#  CMS2           4  292   22         4         6
#  CMS3          13    8  113         2         3
#  CMS4-TME(+)    6    0    0       113         4
#  CMS4-TME(-)    1    5    0        47       107
acc.SYSU <- sum(diag(tab)) / sum(tab)
# 0.8329356
# Figure3B
# python
from mlxtend.plotting import plot_confusion_matrix
import matplotlib.pyplot as plt
import numpy as np
multiclass = np.array([[73, 2, 3, 3, 7],
                       [4, 292, 22, 4, 6],
                       [13, 8, 113, 2, 3],
                       [6, 0, 0, 113, 4], 
                       [1, 5, 0, 47, 107]])
class_names = ['CMS1', 'CMS2', 'CMS3', 'CMS4-TME(+)', 'CMS4-TME(-)']
fig, ax = plot_confusion_matrix(conf_mat=multiclass,
                                colorbar=True,
                                show_absolute=True,
                                show_normed=True,
                                class_names=class_names)
plt.savefig("/mnt/wutan/data/9.CMSPlus/Figure3/Figure3B2.pdf", format="pdf", bbox_inches="tight")
plt.show()


### MVRM
load("/mnt/wutan/data/9.CMSPlus/Figure3/publicdata/download/public.datasets.allgenes.MVRM.RData")  ## MVRM
library(CMSPlus);
library(GSVA);
library(CMSclassifier);
library(ComplexHeatmap);
library(circlize);
library(randomForest);
library(ranger);
library(parallel)

# substr data
cmslabepub <- read.table("/mnt/wutan/data/9.CMSPlus/Figure3/publicdata/cms_labels_public_all.txt", header=T, sep="\t", quote="")
MVRM$clindata$CMS <- factor(cmslabepub$CMS_network[match(rownames(MVRM$clindata), cmslabepub$sample)], levels=c("CMS1", "CMS2", "CMS3", "CMS4"))
#Figure3E
res <- CMSPlus(exp2symbol=MVRM$exp, CMSlabel=MVRM$clindata$CMS)
CMSPlusColor=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', `CMS4-TME+`='#3C5488', `CMS4-TME-`='#8491B4')
MVRM$clindata$CMSPlus <- factor(res$CMSPlusLabel$nearest, levels=names(CMSPlusColor))
save(res, MVRM, file="/mnt/wutan/data/9.CMSPlus/Figure3/MVRM.CMSPlus.RData")
write.table(data.frame(Model="testing", Dataset=MVRM$clindata$datasets, Sample=rownames(MVRM$clindata), `Clustered label`=NA, `Predicted label`=MVRM$clindata$CMSPlus), file="/mnt/wutan/data/9.CMSPlus/Figure3/MVRM.label.txt", col.names=T, row.names=T, sep="\t", quote=F)
# Figure3C
MVRM$clindata$CMSPredict <- factor(substr(as.character(MVRM$clindata$CMSPlus), 1, 4), levels=c("CMS1", "CMS2", "CMS3", "CMS4"))
tab <- table(MVRM$clindata$CMS, MVRM$clindata$CMSPredict)
tab 
#       CMS1 CMS2 CMS3 CMS4
#  CMS1   41    1    3    8
#  CMS2    0   95    4   15
#  CMS3    3    2   34    0
#  CMS4    0    4    0   66
acc.MVRM <- sum(diag(tab)) / sum(tab)
# [1] 0.8550725
# python
from mlxtend.plotting import plot_confusion_matrix
import matplotlib.pyplot as plt
import numpy as np
multiclass = np.array([[41, 1, 3, 8],
                       [0, 95, 4, 15],
                       [3, 2, 34, 0],
                       [0, 4, 0, 66]])
class_names = ['CMS1', 'CMS2', 'CMS3', 'CMS4']
fig, ax = plot_confusion_matrix(conf_mat=multiclass,
                                colorbar=True,
                                show_absolute=True,
                                show_normed=True,
                                class_names=class_names)
plt.savefig("/mnt/wutan/data/9.CMSPlus/Figure3/Figure3C.pdf", format="pdf", bbox_inches="tight")
plt.show()
#Figure3G
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
cmssdirfs <- SurvivalFGPlotSimplify(MVRM$clindata, "dss.delay", "dss.event", "CMSPlus", Color=CMSPlusColor, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
#FigureS4
cmssdirfs <- SurvivalFGPlotSimplify(MVRM$clindata, "rfs.delay", "rfs.event", "CMSPlus", Color=CMSPlusColor, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
CMSColor = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73")
cmssdirfs <- SurvivalFGPlotSimplify(MVRM$clindata, "dss.delay", "dss.event", "CMS", Color=CMSColor, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(MVRM$clindata, "rfs.delay", "rfs.event", "CMS", Color=CMSColor, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
