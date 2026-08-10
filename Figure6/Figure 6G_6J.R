############################################################
#  Figure 6G and Figure 6H  Kaplan–Meier DFS curves
#
#  Instructions:
#  1. Choose dataset:
#       file_path <- "table/SYSU_CT_info.xlsx"
#       or
#       file_path <- "table/Liaoning_CT_info.xlsx"
#
#  2. Set figure title:
#       figure_title <- "SYSU-CT Cohort"      # for Figure 6G
#       figure_title <- "Liaoning-CT Cohort"  # for Figure 6H
#
#  Output includes survival curve + risk table,
#  saved as PNG / PDF in folder "img_plot/".
############################################################

# 1. Install packages if needed
# install.packages(c("survival", "survminer", "extrafont", "openxlsx", "svglite"))

# 2. Load required libraries
library(survival)
library(survminer)
library(extrafont)
library(openxlsx)

# 3. Set environment
loadfonts(device = "win", quiet = TRUE)
par(family = "Arial")

if (!dir.exists("img_plot")) dir.create("img_plot")

# ==========================================================
#                 Select file and title here
# ==========================================================

file_path <- "table/SYSU_CT_info.xlsx"      # or "table/Liaoning_CT_info.xlsx"
figure_title <- "SYSU-CT"  # or "Liaoning‑CT"

# ==========================================================
#                 Data loading and processing
# ==========================================================

data <- read.xlsx(file_path, sheet = "Sheet1")
cat("Data loaded successfully:", dim(data)[1], "rows ×", dim(data)[2], "columns\n")

# Truncate follow‑up beyond 80 months
data$dfs.event <- ifelse(data$dfs.delay > 80, 0, data$dfs.event)
data$dfs.delay <- ifelse(data$dfs.delay > 80, 80, data$dfs.delay)

# Define survival object
surv_object <- Surv(time = data$dfs.delay, event = data$dfs.event)
fit <- survfit(surv_object ~ CMSPlus_predlabel, data = data)

# Compute p-value
surv_pval <- surv_pvalue(fit, data = data)$pval
pval_text <- sprintf("p = %.2e", surv_pval)

# ==========================================================
#                 Kaplan–Meier Plot
# ==========================================================

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
  risk.table.fontsize = 5,
  legend.labs = c("CMS1", "CMS2", "CMS3", "CMS4-TME+", "CMS4-TME-"),
  legend.title = "",
  legend = c(0.7, 0.4),
  font.legend = c(14, "plain"),
  surv.median.line = "none",
  palette = c("#E69F00", "#377EB8", "#CC7AA7", "#3D5588", "#8491B5"),
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
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10))
    ),
  tables.theme = theme_survminer() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_blank(),
      plot.margin = margin(0, 2.1, 0, 5.8, unit = "mm"),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14, margin = margin(t = 8))
    )
)

p$plot <- p$plot +
  theme(axis.title.x = element_blank()) +
  ggtitle(figure_title)

# Display plot
print(p)

# ==========================================================
#                 Save Figures (PNG, PDF, SVG)
# ==========================================================
if (!dir.exists("img_plot")) dir.create("img_plot")

combined_plot <- arrange_ggsurvplots(
  list(p),
  print = FALSE,
  ncol = 1,
  nrow = 1
)

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

cat("Full KM figure with risk table saved to 'img_plot/' →",
    basename(png_file), "\n")
