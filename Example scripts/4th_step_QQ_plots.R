# QQ plots
## load packages
library(ggplot2)
library(ggrepel)
library(CMplot)

data <- $path_to_TWAS

setwd("/Users/leofanfever/Documents/Biostatistics/TIGAR/Thesis/Manuscript/Figures/Supplementary/QQplots/")
columns <- c("p_DPR", "p_EN", "p_FUSION", "p_ACAT")

for (i in 1:4) {
  p_values <- as.numeric(na.omit(data[[columns[i]]]))
  plot_data <- data.frame(Gene = 1:length(p_values) ,CHR = 1:length(p_values),BP = 1:length(p_values),P = p_values)
  plot_data <- plot_data[which(plot_data[,4]!=0),]
  
  CMplot(plot_data,
         plot.type = "q",
         conf.int = TRUE,
         box = TRUE,
         file = "pdf",
         main = "",
         file.name = paste0("$Celltype_", columns[i]),
         dpi = 300,
         file.output = TRUE,
         width = 5,
         height = 5)
}