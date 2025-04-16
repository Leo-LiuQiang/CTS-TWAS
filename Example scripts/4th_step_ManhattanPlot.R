# Manhattan plot
options(stringsAsFactors=FALSE)

library(ggplot2)
library(ggrepel)
library(CMplot)

manhattan_plot <- function(data, # dataframe (columns: 'CHR','POS','Pvalue','label_text') 
                           ## genes: 
                           point_colors=c('#44B5AD','#3A948E','#36807A','#2f615d'), # color for non-sig genes 
                           sig_color1='#FD923F', # color for sig. genes without labels 
                           sig_color2='#D92B26', # color for sig. genes with labels 
                           point_size = 2.5, # from second function
                           point_alpha=0.9, # transparency for genes 
                           ## significance level: 
                           sig_level=0.05, # significance level value 
                           sig_level_line_col='black', # line color 
                           sig_linetype='dashed', # linetype 
                           sig_line_size=1, # line size 
                           ## plot theme, other plot variables: 
                           chr_vec=1:22, # chromosomes to plot 
                           chr_gap=80, # gap between chromosomes in x-axis from second function 
                           theme=theme_grey(), # ggplot2 theme (can use custom theme) from second function 
                           plot_bg_col=NULL, # background color if different from theme 
                           panel_border=element_blank(), # from second function
                           text_size=13, # text size from second function 
                           ## point labelling: 
                           geom_label_size=2.5, # label text size from second function
                           label_fill='white', # label background color 
                           label_col='black', # label border color 
                           label_seg_col='black', # color of line from label to point 
                           min_segment_length=0.01, # minimum length of line from label to point 
                           segment_size=0.2, # line from label to point 
                           label_force=2, # force of repulsion between overlapping text labels 
                           point_padding=1e-06, # padding around genes 
                           box_padding=0.3, # box padding from second function
                           seed=NA, # optional seed for generating label positions 
                           max_iter=20000, # number of iterations to use to generate label positions from second function
                           max.overlaps=10 # max overlaps from second function
){ 
  # setup dataframe for plotting; get plot positions from chromosome positions 
  plot_data <- NULL # dataframe 
  endPos <- 0  # place on x-axis where last chromosome ended 
  x_axis_chr_breaks <- NULL # chromosome label positions 
  x_axis_chr_labels <- NULL # list of chromosomes to label 
  for (chr in chr_vec) { 
    # get data for chr 
    temp <- data[data$CHR==chr, ] 
    if (nrow(temp) > 0) { 
      # append chromosome to list of chromosomes to label 
      x_axis_chr_labels <- c(x_axis_chr_labels, chr) 
      # get unique positions for this chr 
      uniq_pos <- sort(unique(temp$POS)) 
      uniq_pos <- setNames(order(uniq_pos), uniq_pos) 
      # set POS to order value of the unique positions 
      temp$POS <- uniq_pos[as.character(temp$POS)] 
      # get plot positions for genes on this chr 
      temp$plotPos <- (temp$POS - min(temp$POS, na.rm=TRUE) ) + endPos + 1 
      # get new end position based on max position 
      endPos <- max(temp$plotPos, na.rm=TRUE) + chr_gap 
      # append label position 
      x_axis_chr_breaks <- c(x_axis_chr_breaks, mean(temp$plotPos, na.rm=TRUE) ) 
      # add rows to plot_data 
      plot_data <- rbind(plot_data, temp) 
    } 
  } 
  # set min, max values for axes 
  min_x <- min(plot_data$plotPos) 
  max_x <- max(plot_data$plotPos) 
  max_y <- max(-log10(plot_data$Pvalue)) 
  # plot 
  p <- ggplot(data=plot_data,  
              aes(x=plotPos, y=-log10(Pvalue), label=label_text)) +  
    # non-sig. genes: 
    geom_point(data=subset(plot_data, Pvalue >= sig_level), 
               aes(x=plotPos, y=-log10(Pvalue), color=factor(CHR)), 
               size=point_size, alpha=point_alpha) +  
    scale_color_manual(values=rep(point_colors, 22)) + 
    # sig. genes 
    geom_point(data=subset(plot_data, Pvalue < sig_level), 
               aes(x=plotPos, y=-log10(Pvalue), fill=factor(CHR)), 
               size=ifelse(subset(plot_data, Pvalue < sig_level)$label_text=='', point_size + 0.25, point_size + 0.5), 
               color=ifelse(subset(plot_data, Pvalue < sig_level)$label_text=='', sig_color1, sig_color2), 
               alpha=point_alpha) + 
    # add labels 
    geom_label_repel(data=subset(plot_data, Pvalue < sig_level),  
                     min.segment.length=min_segment_length, 
                     segment.size=segment_size, 
                     segment.color=label_seg_col, 
                     box.padding=box_padding, 
                     size=geom_label_size,  
                     alpha=1, 
                     ylim=c(-log10(sig_level), max_y), 
                     xlim=c(min_x, max_x), 
                     force=label_force, 
                     point.padding=point_padding, 
                     max.iter=max_iter, 
                     colour=label_col, 
                     fill=label_fill, 
                     seed=seed) + 
    # significance level line 
    geom_hline(yintercept=-log10(sig_level),  
               linetype=sig_linetype,  
               linewidth=sig_line_size,  
               color=sig_level_line_col) + 
    # remove legend 
    guides(color='none', fill='none') +  
    # set axes titles 
    labs(x='Chromosome', y=bquote(-"log"[10]("q-value"))) +  
    # x-axis labels, breaks 
    scale_x_continuous(breaks=x_axis_chr_breaks,  
                       labels=x_axis_chr_labels,  
                       expand=c(0.01, 0)) +  
    coord_cartesian(clip='off') + 
    # pad y-axis 
    scale_y_continuous(expand=c(0.05, 0)) + 
    theme + 
    theme( 
      text=element_text(size=text_size, face='bold'), 
      axis.title.y=element_text(size=text_size), 
      axis.text.x=element_text(size=text_size-1, angle=-90, vjust=0.5, hjust=0), 
      panel.grid.major.x=element_blank(), 
      panel.grid.minor.x=element_blank(), 
      panel.border=panel_border, 
      plot.background=element_rect(fill=plot_bg_col)) 
  
  return(p) 
}
