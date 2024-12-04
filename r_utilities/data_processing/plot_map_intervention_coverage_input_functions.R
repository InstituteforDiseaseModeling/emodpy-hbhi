# plot_map_intervention_coverage_input_functions.R


library(raster)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(tidyverse)
library(sf)
library(reshape2)
library(pals)
library(prettyGraphs)



##########################################################
# functions for ggplot version
##########################################################

# function to combine multiple plots for same intervention that share a legend
grid_arrange_shared_legend_plotlist =function(..., 
                                              plotlist=NULL,
                                              ncol = length(list(...)),
                                              nrow = NULL,
                                              position = c("bottom", "right")) {
  
  plots <- c(list(...), plotlist)
  
  if (is.null(nrow)) nrow = ceiling(length(plots)/ncol)
  
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position) + guides(fill = guide_legend(reverse=T)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}

change_legend_size <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

# function to create and save maps of intervention coverages
create_coverage_input_maps = function(inter_input, inter_years, output_filename, colorscale, min_value, max_value, num_colors){
  plot_list = list()
  for(yy in 1:length(inter_years)){
    inter_input_cur = inter_input[inter_input$year == inter_years[yy],]
    inter_input_cur$output_value = inter_input_cur[[coverage_colname]]
    admin_cur = admin_shapefile %>%
      left_join(inter_input_cur, by=c('NOMDEP' = 'admin_name')) %>%
      mutate(binned_values=cut(output_value,
                               breaks=round(seq((min_value), (max_value), length.out = (num_colors+1)),2)))
    plot_list[[yy]] = ggplot(admin_cur) +
      geom_sf(aes(fill=binned_values), size=0.1) +
      scale_fill_manual(values=setNames(colorscale, levels(admin_cur$binned_values)), drop=FALSE, name='coverage', na.value='grey96') + 
      ggtitle(inter_years[yy]) +
      guides(fill = guide_legend(reverse=T)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5)) 
    plot_list[[yy]] = change_legend_size(plot_list[[yy]], pointSize=10, textSize=10, spaceLegend=1)
  }
  
  gg = grid_arrange_shared_legend_plotlist(plotlist=plot_list, ncol=length(plot_list), position='right')
  ggsave(output_filename, gg, width = (length(inter_years)+1)*1.8, height=2.3, units='in', dpi=800)
}



# function to create and save maps of intervention coverages
create_maps_state_groups = function(admin_shapefile_filepath, shapefile_admin_colname, admin_group_df, column_name, output_filename, colorscale){
  admin_shapefile = shapefile(admin_shapefile_filepath)
  # standardize shapefile names
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=archetype_info$LGA, origin_names=admin_shapefile$NOMDEP)
  
  admin_group_df$output_value = admin_group_df[[column_name]]
  admin_cur = admin_shapefile %>%
    left_join(inter_input_cur, by=c('NOMDEP' = 'admin_name'))
  gg = ggplot(admin_cur) +
    geom_sf(aes(fill=binned_values), size=0.1) +
    scale_fill_manual(values=setNames(colorscale, levels(admin_cur$output_value)), drop=FALSE, name='group', na.value='grey96') + 
    ggtitle(inter_years[yy]) +
    guides(fill = guide_legend(reverse=T)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) 

  ggsave(output_filename, gg, width = (2)*1.8, height=2.3, units='in', dpi=800)
}
