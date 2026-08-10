############################################################
#  Kaplan–Meier DFS curves by CMSPlus subtype
#
#  Instructions:
#  1. Choose dataset:
#       file_path <- "table/SYSU_CT_info.xlsx"
#       or
#       file_path <- "table/Liaoning_CT_info.xlsx"
#
#  2. Choose subtype:
#       subtype <- "CMS1"
#       subtype <- "CMS2"
#       subtype <- "CMS3"
#       subtype <- "CMS4-TME+"
#       subtype <- "CMS4-TME-"
#
#  Output includes survival curve + risk table,
#  saved as PNG / PDF in folder "img_plot/".
############################################################

library(survival)
library(survminer)
library(coxphf)
library(openxlsx)
library(ggplot2)

############################################################
#                 Parameter settings
############################################################
file_path <- "table/SYSU_CT_info.xlsx"
subtype   <- "CMS1"   # or CMS1 / CMS2 / CMS3 / CMS4-TME-

############################################################
#                 Load and preprocess data
############################################################
data <- read.xlsx(file_path, sheet = 1)
cat("Loaded:", basename(file_path), "\n")

data <- data[data$CMSPlus_pred == subtype, ]
cat("Subtype:", subtype, " n =", nrow(data), "\n")

# Truncate follow-up beyond 80 months
data$dfs.event <- ifelse(data$dfs.delay > 80, 0, data$dfs.event)
data$dfs.delay <- ifelse(data$dfs.delay > 80, 80, data$dfs.delay)

surv_object <- Surv(time = data$dfs.delay, event = data$dfs.event)

############################################################
#                 Fit survival and Cox models
############################################################
fit <- survfit(surv_object ~ Chemotherapy.adjuvant, data = data)
cox_fit <- coxph(surv_object ~ Chemotherapy.adjuvant, data = data)
summary_cox <- summary(cox_fit)

HR <- summary_cox$coefficients[1, 2]
lower_CI <- summary_cox$conf.int[1, 3]
upper_CI <- summary_cox$conf.int[1, 4]
p_val <- summary_cox$coefficients[1, 5]
cat(sprintf("HR = %.2f (95%% CI %.2f - %.2f), p = %.4f\n", HR, lower_CI, upper_CI, p_val))

surv_pval <- surv_pvalue(fit, data = data)$pval
pval_text <- sprintf("p = %.2e", surv_pval)

############################################################
#                 Dynamic color palette
############################################################
palette_map <- list(
  "CMS1" = c("#E69F00", "#674035"),
  "CMS2" = c("#377EB8", "#674035"),
  "CMS3" = c("#CC7AA7", "#674035"),
  "CMS4-TME+" = c("#3D5588", "#674035"),
  "CMS4-TME-" = c("#8491B5", "#674035")
)

# Default palette if not in the map
if (!(subtype %in% names(palette_map))) {
  palette_colors <- c("#469D76", "#674035")
} else {
  palette_colors <- palette_map[[subtype]]
}

############################################################
#                 Figure title
############################################################
cohort_name <- ifelse(grepl("SYSU", file_path), "SYSU-CT", "Liaoning-CT")
figure_title <- paste0(cohort_name, " ", subtype)

############################################################
#                 Kaplan–Meier plot + risk table
############################################################
p <- ggsurvplot(
  fit,
  data = data,
  conf.int = FALSE,
  pval = pval_text,
  pval.coord = c(0, 0.05),
  pval.size = 6,
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  risk.table.height = 0.25,
  risk.table.fontsize = 5.5,
  legend.labs = c("Surgery alone", "Adjuvant chemotherapy"),
  legend.title = "",
  legend = c(0.5, 0.5),
  font.legend = c(16, "plain"),
  surv.median.line = "none",
  palette = palette_colors,
  ylim = c(0, 1),
  break.x.by = 20,
  font.x = c(18, "plain", "black"),
  font.y = c(18, "plain", "black"),
  font.xtickslab = c(16, "plain", "black"),
  font.ytickslab = c(16, "plain", "black"),
  xlab = "Follow-up (Months)",
  ylab = "DFS",
  ggtheme = theme_survminer() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10))
    ),
  tables.theme = theme_survminer() +
    theme(
      plot.title = element_blank(),
      plot.margin = margin(0, 2.1, 0, 7.2, unit = "mm"),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 16, margin = margin(t = 8))
    )
)

p$plot <- p$plot + theme(axis.title.x = element_blank()) + ggtitle(figure_title)
print(p)

############################################################
#                 Combine and save (includes risk table)
############################################################
if (!dir.exists("img_plot")) dir.create("img_plot")

combined_plot <- arrange_ggsurvplots(list(p), print = FALSE, ncol = 1, nrow = 1)

output_prefix <- paste0(
  gsub(".xlsx", "", basename(file_path)),
  "_",
  gsub("\\+", "plus", gsub("-", "minus", subtype))
)

png_file <- paste0("img_plot/", output_prefix, "_KM.png")
pdf_file <- paste0("img_plot/", output_prefix, "_KM.pdf")

ggsave(png_file, combined_plot, width = 6, height = 5.5, dpi = 600)
ggsave(pdf_file, combined_plot, width = 6, height = 5.5)
ggsave(svg_file, combined_plot, width = 6, height = 5.5)

cat("Saved Figure S7 with full curve and risk table to:", basename(png_file), "\n")
