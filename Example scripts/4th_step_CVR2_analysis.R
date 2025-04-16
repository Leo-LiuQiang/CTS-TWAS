DPR <- $path_to_DPR_CVR2
EN <- $path_to_EN_CVR2
FUSION <- $path_to_FUSION_CVR2

FUSION$CVR2 <- apply(FUSION[, c("top1", "lasso", "enet", "blup")], 1, max, na.rm = TRUE)

dpr_df <- DPR[, c("TargetID", "CVR2")]
colnames(dpr_df)[2] <- "DPR_CVR2"

en_df <- EN[, c("TargetID", "CVR2")]
colnames(en_df)[2] <- "EN_CVR2"

fusion_df <- FUSION[, c("TargetID", "CVR2")]
colnames(fusion_df)[2] <- "FUSION_CVR2"

merged_df <- merge(dpr_df, en_df, by = "TargetID", all = TRUE)
merged_df <- merge(merged_df, fusion_df, by = "TargetID", all = TRUE)
merged_df[is.na(merged_df)] <- 0

library(ggplot2)
library(ggExtra)
library(patchwork)

create_scatter_plot <- function(df, x_var, y_var, x_label, y_label, color_x, color_y, label_x, label_y) {
  df$color <- with(df, ifelse(get(x_var) > 0.005 & get(y_var) <= 0.005, "x_only",
                              ifelse(get(y_var) > 0.005 & get(x_var) <= 0.005, "y_only",
                                     ifelse(get(x_var) > 0.005 & get(y_var) > 0.005, "both", "neither"))))
  
  p <- ggplot(df, aes_string(x = x_var, y = y_var)) +
    geom_point(aes(color = color), alpha = 0.7) +
    scale_color_manual(
      values = c("x_only" = color_x, "y_only" = color_y, "both" = '#00BFC4', "neither" = 'grey'),
      labels = c("x_only" = paste(label_x, "Only"), 
                 "y_only" = paste(label_y, "Only"), 
                 "both" = "Both", 
                 "neither" = "Neither")
    ) +
    theme_minimal() +
    labs(x = x_label, y = y_label) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    theme(panel.grid = element_blank(), 
          axis.line = element_line(color = "black", size = 0.5), 
          axis.ticks = element_line(color = "black", size = 0.5), 
          axis.text = element_text(color = "black"),  
          legend.position = "none",
          plot.title = element_blank()) +
    guides(color = guide_legend(title = "R2 > 0.005"))
  
  ggMarginal(p, type = "histogram", fill = "grey", bins = 250, size = 5)
}

plot1 <- create_scatter_plot(merged_df, "DPR_CVR2", "EN_CVR2", 
                             "CVR2 of TIGAR Model", "CVR2 of Elastic_Net Model", 
                             '#F8766D', '#7CAE00', "TIGAR Model", "Elastic_Net Model")

plot2 <- create_scatter_plot(merged_df, "DPR_CVR2", "FUSION_CVR2", 
                             "CVR2 of TIGAR Model", "CVR2 of FUSION Model", 
                             '#F8766D', '#619CFF', "TIGAR Model", "FUSION Model")

plot3 <- create_scatter_plot(merged_df, "EN_CVR2", "FUSION_CVR2", 
                             "CVR2 of Elastic_Net Model", "CVR2 of FUSION Model", 
                             '#7CAE00', '#619CFF', "Elastic_Net Model", "FUSION Model")

wrapped_plot1 <- wrap_elements(plot1)
wrapped_plot2 <- wrap_elements(plot2)
wrapped_plot3 <- wrap_elements(plot3)

legend_df <- data.frame(Model = c("x_only", "y_only", "z_only", "both", "neither"))
legend_plot <- ggplot(legend_df, aes(x = Model, y = Model, color = Model)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("x_only" = '#F8766D', "y_only" =  '#7CAE00', "z_only" = '#619CFF', "both" = "#00BFC4", "neither" = "grey"),
    labels = c("x_only" = "TIGAR Model Only", 
               "y_only" = "Elastic_Net Model Only",
               "z_only" = "FUSION Model Only",
               "both" = "Both", 
               "neither" = "Neither")) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(
    title = "R2 > 0.005",
    override.aes = list(size = 4),
    nrow = 1
  ))

extract_legend <- function(plot) {
  gtable <- ggplot_gtable(ggplot_build(plot))
  legend <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend_only <- extract_legend(legend_plot)

combined_plot <- (wrapped_plot1 | wrapped_plot2 | wrapped_plot3) / wrap_elements(grid::grobTree(legend_only)) + plot_layout(heights = c(10, 1))
