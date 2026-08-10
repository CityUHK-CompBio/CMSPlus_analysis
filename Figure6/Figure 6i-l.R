############################################################
#  Figure 6i to Figure 6l  Kaplan–Meier DFS curves
#  Chemotherapy benefit analysis by CMS4 subtype
#
#  Instructions:
#  1. Choose dataset:
#       file_path <- "table/SYSU_CT_info.xlsx"
#       or
#       file_path <- "table/Liaoning_CT_info.xlsx"
#
#  2. Choose subtype:
#       subtype <- "CMS4-TME+"   # for Figure 6i (SYSU-CT) or Figure 6j (Liaoning-CT)
#       subtype <- "CMS4-TME-"   # for Figure 6k (SYSU-CT) or Figure 6l (Liaoning-CT)
#
#  Figure titles are generated automatically,
#  and outputs (curve + risk table) are saved
#  as PNG / PDF in folder "img_plot/".
############################################################

library(survival)
library(survminer)
library(coxphf)
library(openxlsx)

# ==========================================================
#                Set dataset and subtype
# ==========================================================

file_path <- "table/SYSU_CT_info.xlsx" # or "table/Liaoning_CT_info.xlsx"
subtype   <- "CMS4-TME+"      # or "CMS4-TME-"

# ==========================================================
#                Load and preprocess data
# ==========================================================

data <- read.xlsx(file_path, sheet = 1)
cat("Loaded:", basename(file_path), "\n")

# subset one subtype
data <- data[data$CMSPlus_predlabel == subtype, ]
cat("Subtype:", subtype, "n =", nrow(data), "\n")

# truncate follow-up at 80 months
data$dfs.event <- ifelse(data$dfs.delay > 80, 0, data$dfs.event)
data$dfs.delay <- ifelse(data$dfs.delay > 80, 80, data$dfs.delay)

# define survival object
surv_object <- Surv(time = data$dfs.delay, event = data$dfs.event)

# ==========================================================
#                Fit models
# ==========================================================

fit     <- survfit(surv_object ~ Chemotherapy.adjuvant, data = data)
cox_fit <- coxph(surv_object ~ Chemotherapy.adjuvant, data = data)
summary_cox <- summary(cox_fit)

# hazard ratio
HR        <- summary_cox$coefficients[1, 2]
lower_CI  <- summary_cox$conf.int[1, 3]
upper_CI  <- summary_cox$conf.int[1, 4]
p_val     <- summary_cox$coefficients[1, 5]

cat(sprintf("HR = %.2f (95%% CI %.2f–%.2f), p = %.4f\n", HR, lower_CI, upper_CI, p_val))

# p-value for KM curve
surv_pval  <- surv_pvalue(fit, data = data)$pval
pval_text  <- sprintf("p = %.2e", surv_pval)

# ==========================================================
#                Plot
# ==========================================================

title_text <- switch(
  paste0(basename(file_path), "_", subtype),
  "SYSU_CT_info.xlsx_CMS4-TME+" = "SYSU-CT CMS4-TME+",
  "SYSU_CT_info.xlsx_CMS4-TME-" = "SYSU-CT CMS4-TME-",
  "Liaoning_CT_info.xlsx_CMS4-TME+" = "Liaoning-CT CMS4-TME+",
  "Liaoning_CT_info.xlsx_CMS4-TME-" = "Liaoning-CT CMS4-TME-",
  paste("Chemotherapy effect |", basename(file_path), subtype)
)

p <- ggsurvplot(
  fit,
  data = data,
  conf.int = FALSE,
  pval = pval_text,
  pval.coord = c(0, 0.05),
  pval.size = 5,
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  risk.table.height = 0.25,
  risk.table.fontsize = 5,
  legend.labs = c("Surgery alone", "Adjuvant chemotherapy"),
  legend.title = "",
  legend = c(0.4, 0.4),
  font.legend = c(14, "plain"),
  surv.median.line = "none",
  palette = c("#8491B5", "#674035"),
  ylim = c(0, 1),
  break.x.by = 20,
  font.x = c(16, "plain", "black"),
  font.y = c(16, "plain", "black"),
  font.xtickslab = c(14, "plain", "black"),
  font.ytickslab = c(14, "plain", "black"),
  xlab = "Follow-up (Months)",
  ylab = "DFS",
  ggtheme = theme_survminer() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10))
    ),
  tables.theme = theme_survminer() +
    theme(
      plot.title = element_blank(),
      plot.margin = margin(0, 2.1, 0, 5.8, unit = "mm"),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14, margin = margin(t = 8))
    )
)

p$plot <- p$plot + theme(axis.title.x = element_blank()) + ggtitle(title_text)
print(p)

# ==========================================================
#                Save outputs
# ==========================================================

if (!dir.exists("img_plot")) dir.create("img_plot")

combined_plot <- arrange_ggsurvplots(
  list(p), 
  print = FALSE,
  ncol = 1, nrow = 1
)

prefix <- paste0(gsub(".xlsx", "", basename(file_path)), "_", gsub("\\+", "plus", gsub("-", "minus", subtype)))
png_file <- paste0("img_plot/", prefix, "_KM.png")
pdf_file <- paste0("img_plot/", prefix, "_KM.pdf")

# save everything (plot + risk table)
ggsave(png_file, combined_plot, width = 6, height = 5.5, dpi = 600)
ggsave(pdf_file, combined_plot, width = 6, height = 5.5)

cat("Saved full figure with risk table:", title_text, "→ img_plot/\n")
