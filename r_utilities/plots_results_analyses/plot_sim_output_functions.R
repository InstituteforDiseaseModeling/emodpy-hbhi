# plot_sim_output_functions.R

# library(rgdal)
library(raster)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(sf)
library(reshape2)
library(data.table)
library(dplyr)
library(geofacet)
library(ggpattern)


separate_plot_text_size=12
text_size = 15
save_plots = TRUE


####################################################################################
# barplots for burden relative to BAU: percent reduction
####################################################################################

plot_relative_burden_barplots = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                         barplot_start_year, barplot_end_year, 
                                         pyr, chw_cov,
                                         scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, separate_plots_flag=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                         include_to_present=TRUE, burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality', 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all', 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }

  # first comparison name is to-present (skip it), second is BAU (use as reference), comparison scenarios start at the third index
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2  
  }
  # iterate through comparison scenarios, calculating the burden reduction of all metrics relative to BAU (seedwise comparisons, so one output for each run). Combine all scenario reductions into a dataframe (each scenario set in separate rows)
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, 
                                             start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    # only save relevant columns for plotting
    relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    }else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_names[comparison_start_index:length(scenario_names)])
  
  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_x = 0
  standard_max_x = 0.1
  cur_min = min(relative_burden_all_df[,2:(1+length(burden_colnames))])
  cur_max = max(relative_burden_all_df[,2:(1+length(burden_colnames))])
  if(cur_min < standard_min_x) standard_min_x = cur_min
  if(cur_max > standard_max_x) standard_max_x = cur_max
  
  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name = burden_metric_names[bb]
    select_col_names = c(current_burden_name, 'scenario')
    # get mean, min, and max among all runs for this burden metric
    rel_burden_agg = as.data.frame(relative_burden_all_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                       max_rel = max(get(current_burden_name)),
                       min_rel = min(get(current_burden_name)))
    
    gg_list[[bb]] = ggplot(rel_burden_agg) + 
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity') +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_x, standard_max_x)) +   # turn into percent reduction
      ylab('Percent reduction') + 
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values = scenario_palette) + 
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size), 
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg_list[[bb]] = gg_list[[bb]] +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour="black", alpha=0.9, size=1) 
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] + 
        ylab('Percent reduction in burden \n ((Current - Plan) / Current) * 100') + 
        theme(legend.position='none', plot.title = element_blank(), text=element_text(size =separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/','barplot_percent_reduction_', burden_metric_name,'_',district_subset,'.png'), separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }
  # for each burden type, 
  # get mean, min, and max among all runs for each burden metric, each saved as a separate column
  # create barplot for each burden type (using columns of dataframe, separate bar for each scenario)
  
  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  # remove legend from main plots
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position = "none")  + theme(text = element_text(size = text_size))   
  }
  
  if(save_plots){
    gg_saved = grid.arrange(grobs = gg_list[-1], layout_matrix = matrix(c(1:(length(burden_colnames))), nrow=2, byrow=FALSE))
    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_percent_reduction_burden_', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg_saved, dpi=600, width=14, height=7, units='in')
  }
  
  # ----- combine all burden plots ----- #
  # gg = grid.arrange(grobs = gg_list, layout_matrix = matrix(c(1,1,2:(length(burden_colnames)+1)), ncol=2, byrow=TRUE))
  gg = grid.arrange(grobs = gg_list, layout_matrix = rbind(matrix(rep(1, length(burden_colnames)/2), nrow=1), matrix(2:(length(burden_colnames)+1), nrow=2, byrow=FALSE)))
  
  return(gg)
}






####################################################################################
# barplots for burden relative to BAU at the state level, displayed in a state grid
####################################################################################

plot_relative_burden_barplots_by_state = function(sim_future_output_dir, pop_filepath,grid_layout_state_locations,
                                                 barplot_start_year, barplot_end_year, 
                                                 pyr, chw_cov,
                                                 scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                                 burden_metric_subset=c(), include_to_present=TRUE, file_suffix=''){
  
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality', 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all', 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }
  
  # first comparison name is to-present (skip it), second is BAU (use as reference), comparison scenarios start at the third index
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2  
  }
  # iterate through comparison scenarios, calculating the burden reduction of all metrics relative to BAU (seedwise comparisons, so one output for each run). Combine all scenario reductions into a dataframe (each scenario set in separate rows)
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden_by_state(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, 
                                                   start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    # only save relevant columns for plotting
    relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', 'State', burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    }else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_names[comparison_start_index:length(scenario_names)])
  
  # # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  # standard_min_y = 0
  # standard_max_y = 0.1
  # cur_min = min(relative_burden_all_df[,2:(1+length(burden_colnames))])
  # cur_max = max(relative_burden_all_df[,2:(1+length(burden_colnames))])
  # if(cur_min < standard_min_y) standard_min_y = cur_min
  # if(cur_max > standard_max_y) standard_max_y = cur_max
  

  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name = burden_metric_names[bb]
    select_col_names = c(current_burden_name, 'scenario', 'State')
    # get mean, min, and max among all runs for this burden metric
    rel_burden_agg = as.data.frame(relative_burden_all_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
      dplyr::group_by(scenario, State) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                       max_rel = max(get(current_burden_name)),
                       min_rel = min(get(current_burden_name)))
    
    rel_burden_agg$code = rel_burden_agg$State
    gg = ggplot(rel_burden_agg) + 
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity') +
      scale_y_continuous(labels=percent_format(), n.breaks=4) + #,limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
      ylab('Percent reduction in burden \n ((Baseline - Plan) / Baseline) * 100') + 
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values = scenario_palette) + 
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size), 
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm')) + 
      facet_geo(~code, grid = grid_layout_state_locations, label="name") #, scales='free') 
    
    if(show_error_bar){
      gg = gg +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour="black", alpha=0.9, size=1) 
    }
    ggsave(paste0(sim_future_output_dir, '/_plots/','barplot_percent_reduction_', burden_metric_name,'_stateGrid',file_suffix,'.png'), gg, dpi=600, width=18*.6, height=12*.6, units='in')  # , width=18, height=12, units='in'
  }
}




####################################################################################
# barplots for burden reduction relative to BAU: difference
####################################################################################

plot_difference_burden_barplots = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                         barplot_start_year, barplot_end_year, 
                                         pyr, chw_cov,
                                         scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, separate_plots_flag=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                         include_to_present=TRUE, burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality', 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all', 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }
  
  # first comparison name is to-present (skip it), second is BAU (use as reference), comparison scenarios start at the third index
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2  
  }
  # iterate through comparison scenarios, calculating the burden reduction of all metrics relative to BAU (seedwise comparisons, so one output for each run). Combine all scenario reductions into a dataframe (each scenario set in separate rows)
  difference_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    difference_burden_df = get_difference_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, 
                                             start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    # only save relevant columns for plotting
    difference_burden_df = difference_burden_df[,which(colnames(difference_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(difference_burden_all_df) == 0){
      difference_burden_all_df = difference_burden_df
    }else{
      difference_burden_all_df = rbind(difference_burden_all_df, difference_burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  difference_burden_all_df$scenario = factor(difference_burden_all_df$scenario, levels=scenario_names[comparison_start_index:length(scenario_names)])
  
  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_x = 0
  standard_max_x = 0.1
  cur_min = min(difference_burden_all_df[,2:(1+length(burden_colnames))])
  cur_max = max(difference_burden_all_df[,2:(1+length(burden_colnames))])
  if(cur_min < standard_min_x) standard_min_x = cur_min
  if(cur_max > standard_max_x) standard_max_x = cur_max
  
  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name = burden_metric_names[bb]
    select_col_names = c(current_burden_name, 'scenario')
    # get mean, min, and max among all runs for this burden metric
    rel_burden_agg = as.data.frame(difference_burden_all_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                       max_rel = max(get(current_burden_name)),
                       min_rel = min(get(current_burden_name)))
    
    gg_list[[bb]] = ggplot(rel_burden_agg) + 
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity') +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_x, standard_max_x)) +   # turn into percent reduction
      ylab('Burden averted') + 
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values = scenario_palette) + 
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size), 
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg_list[[bb]] = gg_list[[bb]] +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour="black", alpha=0.9, size=1) 
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] + 
        ylab('Burden averted \n ((Baseline - Plan)') + 
        theme(legend.position='none', plot.title = element_blank(), text=element_text(size =separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/','barplot_burden_averted_', burden_metric_name,'_',district_subset,'.png'), separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }
  # for each burden type, 
  # get mean, min, and max among all runs for each burden metric, each saved as a separate column
  # create barplot for each burden type (using columns of dataframe, separate bar for each scenario)
  
  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  # remove legend from main plots
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position = "none")  + theme(text = element_text(size = text_size))   
  }
  
  if(save_plots){
    gg_saved = grid.arrange(grobs = gg_list[-1], layout_matrix = matrix(c(1:(length(burden_colnames))), nrow=2, byrow=FALSE))
    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_burden_averted_', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg_saved, dpi=600, width=14, height=7, units='in')
  }
  
  # ----- combine all burden plots ----- #
  # gg = grid.arrange(grobs = gg_list, layout_matrix = matrix(c(1,1,2:(length(burden_colnames)+1)), ncol=2, byrow=TRUE))
  gg = grid.arrange(grobs = gg_list, layout_matrix = rbind(matrix(rep(1, length(burden_colnames)/2), nrow=1), matrix(2:(length(burden_colnames)+1), nrow=2, byrow=FALSE)))
  
  return(gg)
}






####################################################################################################################################
# barplot of the impact a specific intervention has in relevant admins 
#  (percent reduction when intervention is included versus matched simulation without the intervention)
####################################################################################################################################


plot_barplot_impact_specific_intervention = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                              barplot_start_year, barplot_end_year, 
                                              pyr, chw_cov,
                                              experiment_names_without, experiment_names_with, scenario_palette, scenario_barfill=NA, intervention_name='PMC', age_group = 'U1', LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                              burden_metric_subset=c(), default_ylim_max=0.03){
  admin_pop = read.csv(pop_filepath)
  comparison_scenario_name = intervention_name
  
  # iterate through the matched pairs of experiments without / with the intervention
  rel_burden_agg = data.frame()
  if(length(experiment_names_without) == length(experiment_names_with)){
    for(ii in 1:length(experiment_names_without)){
      # first experiment is without interventions, second experiment is with intervention
      reference_experiment_name = experiment_names_without[ii]
      # calculating the burden reduction of all metrics relative to no intervention (seedwise comparisons, so one output for each run). 
      comparison_experiment_name = experiment_names_with[ii]  
      
      # set which burden metrics are relevant and get relative burden between simulations
      burden_metrics = c('PfPR', 'incidence', 'directMortality', 'allMortality')
      if(age_group=='U1'){
        burden_colnames = c('average_PfPR_U1', 'incidence_U1', 'direct_death_rate_mean_U1', 'all_death_rate_mean_U1')
        burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'direct mortality (U1)', 'mortality (U1)')
        relative_burden_df = get_relative_U1_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      } else if (age_group=='U5'){
        burden_colnames = c('average_PfPR_U5', 'incidence_U5', 'direct_death_rate_mean_U5', 'all_death_rate_mean_U5')
        burden_metric_names = c('PfPR (U5)', 'incidence (U5)', 'direct mortality (U5)', 'mortality (U5)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }else{
        burden_metrics = c(burden_metrics, 'mLBW_deaths', 'MiP_stillbirths')
        burden_colnames = c('average_PfPR_all', 'incidence_all', 'direct_death_rate_mean_all', 'all_death_rate_mean_all', 'annual_num_mLBW', 'annual_num_mStill')
        burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'direct mortality (all ages)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }
      
  
      # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
      if((length(burden_metric_subset)>=1)){
        burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
        burden_colnames = burden_colnames[burden_metrics_subset_indices]
        burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
      }
      
      # only save relevant columns for plotting
      relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
      
      for(bb in 1:length(burden_colnames)){
        current_burden_name = burden_colnames[bb]
        burden_metric_name = burden_metric_names[bb]
        select_col_names = c(current_burden_name, 'scenario')
        # get mean, min, and max among all runs for this burden metric
        rel_burden_agg_bb = as.data.frame(relative_burden_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
          dplyr::group_by(scenario) %>%
          dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                           max_rel = max(get(current_burden_name)),
                           min_rel = min(get(current_burden_name)))
        rel_burden_agg_bb$burden_metric = burden_metric_name
        rel_burden_agg_bb$scenario_name = experiment_names_with[ii]
        if(nrow(rel_burden_agg)<1){
          rel_burden_agg = rel_burden_agg_bb
        } else{
          rel_burden_agg = merge(rel_burden_agg, rel_burden_agg_bb, all=TRUE)
        }
      }
    }
  }

  rel_burden_agg$burden_metric = gsub('\\(births\\)', '', rel_burden_agg$burden_metric)
  rel_burden_agg$burden_metric = factor(rel_burden_agg$burden_metric, levels=gsub('\\(births\\)', '', burden_metric_names))
  rel_burden_agg$scenario_name = factor(rel_burden_agg$scenario_name, levels=experiment_names_with)

  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_x = 0
  standard_max_x = default_ylim_max
  cur_min = min(rel_burden_agg[,2:4])
  cur_max = max(rel_burden_agg[,2:4])
  if(cur_min < standard_min_x) standard_min_x = cur_min
  if(cur_max > standard_max_x) standard_max_x = cur_max
  if(any(is.na(scenario_barfill))){
    scenario_barfill = rep('none', length(unique(rel_burden_agg$scenario_name)))
    names(scenario_barfill) = unique(rel_burden_agg$scenario_name)
  }
  # original without shading:
  # gg = ggplot(rel_burden_agg) +
  #   geom_bar(aes(x=burden_metric, y=mean_rel, fill=scenario_name), stat='identity', position="dodge") +
  #   scale_y_continuous(labels=percent_format(), limits=c(standard_min_x, standard_max_x)) +   # turn into percent reduction
  #   ylab(paste0('Percent reduction in burden \n ((without ', intervention_name, ' - with ', intervention_name, ') / without ', intervention_name, ') * 100')) +
  #   geom_hline(yintercept=0, color='black') +
  #   ggtitle(paste0('Comparison of burden in proposed ', intervention_name, ' districts')) +
  #   scale_fill_manual(values = scenario_palette) +
  #   theme_classic()+
  #   theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size),
  #         axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
  #         plot.margin=unit(c(0,1,1,0), 'cm'))
  
  gg = ggplot(rel_burden_agg) + 
    # geom_bar(aes(x=burden_metric, y=mean_rel, fill=scenario_name, pattern=scenario_name), stat='identity', position="dodge") +
    scale_y_continuous(labels=percent_format(), limits=c(standard_min_x, standard_max_x)) +   # turn into percent reduction
    ylab(paste0('Percent reduction in burden \n ((without ', intervention_name, ' - with ', intervention_name, ') / without ', intervention_name, ') * 100')) + 
    geom_hline(yintercept=0, color='black') +
    ggtitle(paste0('Comparison of burden in proposed ', intervention_name, ' districts')) + 
    scale_fill_manual(values = scenario_palette) + 
    geom_bar_pattern(aes(x=burden_metric, y=mean_rel, fill=scenario_name, pattern=scenario_name), stat='identity', position="dodge", #position = position_dodge(preserve = "single"),
                     # color = "white",
                     pattern_fill = "white",
                     pattern_linetype=0,
                     pattern_angle = 45,
                     pattern_density = 0.35,
                     pattern_spacing = 0.06,# 0.025,
                     pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = scenario_barfill) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme_classic()+ 
    theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size), 
          axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
          plot.margin=unit(c(0,1,1,0), 'cm'))
  
  if(show_error_bar){
    gg = gg +
      geom_errorbar(aes(x=burden_metric, ymin=min_rel, ymax=max_rel, group=scenario_name), position='dodge',  colour="black", alpha=0.9, size=1) # width=0.4,
  }

  if(save_plots){
    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_', intervention_name, '_percent_reduction_burden_', age_group, '', barplot_start_year, '_', barplot_end_year, '.png'), gg, dpi=600, width=4.8, height=4.8, units='in')
  }
  
  return(gg)
}






####################################################################################################################################
# barplot of the impact a specific intervention has in relevant admins 
#  (burden reduction when intervention is included versus matched simulation without the intervention)
####################################################################################################################################


table_difference_impact_specific_intervention = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                                     barplot_start_year, barplot_end_year, 
                                                     pyr, chw_cov,
                                                     experiment_names_without, experiment_names_with, intervention_name='PMC', age_group = 'U1', LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE,
                                                     burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)
  comparison_scenario_name = intervention_name
  
  # iterate through the matched pairs of experiments without / with the intervention
  rel_burden_agg = data.frame()
  if(length(experiment_names_without) == length(experiment_names_with)){
    for(ii in 1:length(experiment_names_without)){
      # first experiment is without interventions, second experiment is with intervention
      reference_experiment_name = experiment_names_without[ii]
      # calculating the burden reduction of all metrics relative to no intervention (seedwise comparisons, so one output for each run). 
      comparison_experiment_name = experiment_names_with[ii]  
      
      # set which burden metrics are relevant and get relative burden between simulations
      burden_metrics = c('PfPR', 'incidence', 'directMortality', 'allMortality')
      if(age_group=='U1'){
        warning('Have not yet added support for burden reduction for U1... need to add relevant functions')
        # burden_colnames = c('average_PfPR_U1', 'incidence_U1', 'direct_death_rate_mean_U1', 'all_death_rate_mean_U1')
        # burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'direct mortality (U1)', 'mortality (U1)')
        # relative_burden_df = get_relative_U1_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      } else if (age_group=='U5'){
        burden_colnames = c('average_PfPR_U5', 'incidence_U5', 'direct_death_rate_mean_U5', 'all_death_rate_mean_U5')
        burden_metric_names = c('PfPR (U5)', 'incidence (U5)', 'direct mortality (U5)', 'mortality (U5)')
        relative_burden_df = get_difference_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }else{
        burden_metrics = c(burden_metrics, 'mLBW_deaths', 'MiP_stillbirths')
        burden_colnames = c('average_PfPR_all', 'incidence_all', 'direct_death_rate_mean_all', 'all_death_rate_mean_all', 'annual_num_mLBW', 'annual_num_mStill')
        burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'direct mortality (all ages)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
        relative_burden_df = get_difference_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }
      
      
      # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
      if((length(burden_metric_subset)>=1)){
        burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
        burden_colnames = burden_colnames[burden_metrics_subset_indices]
        burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
      }
      
      # only save relevant columns for plotting
      relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
      
      for(bb in 1:length(burden_colnames)){
        current_burden_name = burden_colnames[bb]
        burden_metric_name = burden_metric_names[bb]
        select_col_names = c(current_burden_name, 'scenario')
        # get mean, min, and max among all runs for this burden metric
        rel_burden_agg_bb = as.data.frame(relative_burden_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
          dplyr::group_by(scenario) %>%
          dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                           max_rel = max(get(current_burden_name)),
                           min_rel = min(get(current_burden_name)))
        rel_burden_agg_bb$burden_metric = burden_metric_name
        rel_burden_agg_bb$scenario_name = experiment_names_with[ii]
        if(nrow(rel_burden_agg)<1){
          rel_burden_agg = rel_burden_agg_bb
        } else{
          rel_burden_agg = merge(rel_burden_agg, rel_burden_agg_bb, all=TRUE)
        }
      }
    }
  }
  
  rel_burden_agg$burden_metric = gsub('\\(births\\)', '', rel_burden_agg$burden_metric)
  rel_burden_agg$burden_metric = factor(rel_burden_agg$burden_metric, levels=gsub('\\(births\\)', '', burden_metric_names))
  rel_burden_agg$scenario_name = factor(rel_burden_agg$scenario_name, levels=experiment_names_with)
  
  return(rel_burden_agg)
}







plot_barplot_impact_two_specific_interventions = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                                                     barplot_start_year, barplot_end_year,
                                                     pyr, chw_cov,
                                                     experiment_names_without, experiment_names_with, scenario_palette, intervention_name='PMC', age_group = 'U1', LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                                     intervention_strings = c('Vacc')){
  admin_pop = read.csv(pop_filepath)
  comparison_scenario_name = intervention_name

  # iterate through the matched pairs of experiments without / with the intervention
  rel_burden_agg = data.frame()
  if(length(experiment_names_without) == length(experiment_names_with)){
    for(ii in 1:length(experiment_names_without)){
      # first experiment is without interventions, second experiment is with intervention
      reference_experiment_name = experiment_names_without[ii]
      # calculating the burden reduction of all metrics relative to no intervention (seedwise comparisons, so one output for each run).
      comparison_experiment_name = experiment_names_with[ii]

      # set which burden metrics are relevant and get relative burden between simulations
      if(age_group=='U1'){
        burden_colnames = c('average_PfPR_U1', 'incidence_U1', 'direct_death_rate_mean_U1', 'all_death_rate_mean_U1')
        burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'direct mortality (U1)', 'mortality (U1)')
        relative_burden_df = get_relative_U1_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      } else if (age_group=='U5'){
        burden_colnames = c('average_PfPR_U5', 'incidence_U5', 'direct_death_rate_mean_U5', 'all_death_rate_mean_U5')
        burden_metric_names = c('PfPR (U5)', 'incidence (U5)', 'direct mortality (U5)', 'mortality (U5)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }else{
        burden_colnames = c('average_PfPR_all', 'incidence_all', 'direct_death_rate_mean_all', 'all_death_rate_mean_all')
        burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'direct mortality (all ages)', 'mortality (all ages)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }

      # only save relevant columns for plotting
      relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]

      for(bb in 1:length(burden_colnames)){
        current_burden_name = burden_colnames[bb]
        burden_metric_name = burden_metric_names[bb]
        select_col_names = c(current_burden_name, 'scenario')
        # get mean, min, and max among all runs for this burden metric
        rel_burden_agg_bb = as.data.frame(relative_burden_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
          dplyr::group_by(scenario) %>%
          dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                           max_rel = max(get(current_burden_name)),
                           min_rel = min(get(current_burden_name)))
        rel_burden_agg_bb$burden_metric = burden_metric_name
        rel_burden_agg_bb$scenario_name = experiment_names_with[ii]
        if(length(intervention_strings)>1){
          rel_burden_agg_bb$intervention_info = NA
          for(jj in 1:length(intervention_strings)){
            if(grepl(intervention_strings[jj], experiment_names_without[ii])){
              rel_burden_agg_bb$intervention_info = intervention_strings[jj]
            }
          }
        }
        if(nrow(rel_burden_agg)<1){
          rel_burden_agg = rel_burden_agg_bb
        } else{
          rel_burden_agg = merge(rel_burden_agg, rel_burden_agg_bb, all=TRUE)
        }
      }
    }
  }

  rel_burden_agg$burden_metric = factor(rel_burden_agg$burden_metric, levels=burden_metric_names)

  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_x = 0
  standard_max_x = 0.2
  cur_min = min(rel_burden_agg[,2:4])
  cur_max = max(rel_burden_agg[,2:4])
  if(cur_min < standard_min_x) standard_min_x = cur_min
  if(cur_max > standard_max_x) standard_max_x = cur_max

  
  # create list where each element is a barplot corresponding to one of the interventions in intervention_strings
  gg_list = list()
  for(jj in 1:length(intervention_strings)){
    gg = ggplot(rel_burden_agg[rel_burden_agg$intervention_info == intervention_strings[jj],]) +
      geom_bar(aes(x=burden_metric, y=mean_rel, fill=scenario_name), stat='identity', position="dodge") +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_x, standard_max_x)) +   # turn into percent reduction
      ylab(paste0('Percent reduction in burden \n ((without ', intervention_strings[jj], ' - with ', intervention_name, ') / without ', intervention_strings[jj], ') * 100')) +
      geom_hline(yintercept=0, color='black') +
      ggtitle(paste0('Comparison of burden in proposed ', intervention_name, ' districts')) +
      scale_fill_manual(values = scenario_palette) +
      theme_classic()+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size),
            axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    
    if(show_error_bar){
      gg = gg +
        geom_errorbar(aes(x=burden_metric, ymin=min_rel, ymax=max_rel, group=scenario_name), position='dodge',  colour="black", alpha=0.9, size=1) # width=0.4,
    }
    gg_list[[jj]] = gg
  }
  return(gg_list)
}







######################################################################
# create plot panel with all burden metrics, no intervention info (either showing burden or burden relative to burden in specified year)
######################################################################

plot_simulation_output_burden_all = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                             plot_by_month, min_year, max_year, sim_end_years, relative_year=NA,
                                             pyr='', chw_cov='',
                                             scenario_filepaths, scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, 
                                             separate_plots_flag=FALSE, extend_past_timeseries_year=NA, scenario_linetypes=NA, plot_CI=TRUE, include_U1=FALSE,
                                             burden_metric_subset=c()){
  
  if (!is.na(relative_year)){ if(relative_year<min_year){
    warning('specified minimum year must be <= relative year. Setting min_year to relative_year.')
    min_year = relative_year
  }}
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  pop_sizes = read.csv(pop_filepath)
  pop_sizes = pop_sizes[,c('admin_name','pop_size')]
  # if we include all admins, get list of names from population size dataframe
  if(cur_admins[1] == 'all'){
    cur_admins = unique(pop_sizes$admin_name)
  }
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  
  # ----- malaria burden ----- #
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality', 'mLBW_deaths', 'MiP_stillbirths')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
  burden_colnames = c('PfPR_U5', 'PfPR_MiP_adjusted', 'New_clinical_cases_U5', 'New_Clinical_Cases', 'direct_mortality_nonMiP_U5_mean', 'direct_mortality_nonMiP_mean', 'total_mortality_U5_mean', 'total_mortality_mean', 'mLBW_deaths', 'MiP_stillbirths')    
  if(include_U1){
    burden_metrics = c('PfPR', 'PfPR', 'PfPR', 'incidence','incidence', 'incidence', 'directMortality', 'directMortality', 'directMortality', 'allMortality', 'allMortality', 'allMortality', 'mLBW_deaths', 'MiP_stillbirths')
    burden_metric_names = c('PfPR (U1)', 'PfPR (U5)', 'PfPR (all ages)', 'incidence (U1)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U1)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U1)', 'mortality (U5)', 'mortality (all ages)', 'mLBW mortality (births)', 'stillbirths (births)')
    burden_colnames = c('PfPR_U1', 'PfPR_U5', 'PfPR_MiP_adjusted', 'New_clinical_cases_U1', 'New_clinical_cases_U5', 'New_Clinical_Cases', 'direct_mortality_nonMiP_U1_mean', 'direct_mortality_nonMiP_U5_mean', 'direct_mortality_nonMiP_mean', 'total_mortality_U1_mean', 'total_mortality_U5_mean', 'total_mortality_mean', 'mLBW_deaths', 'MiP_stillbirths')
  }
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
    burden_metrics = burden_metrics[burden_metrics_subset_indices]
  }
    
  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    burden_metric_name = burden_metric_names[bb]
    burden_colname = burden_colnames[bb]
    burden_metric = burden_metrics[bb]
    
    if(grepl('U1', burden_metric_name)){
      age_plotted = 'U1'
    } else if(grepl('U5', burden_metric_name)){
      age_plotted = 'U5'
    } else if(grepl('births', burden_metric_name)){
      age_plotted = 'births'
    } else age_plotted = 'all'
    
    
    # check whether burden output already exists for this comparison
    if(LLIN2y_flag){
      llin2y_string = '_2yLLIN'
    } else{
      llin2y_string = ''
    }
    
    # iterate through scenarios, storing relevant output
    burden_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_sim_output_agg = get_burden_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                     exp_name = scenario_names[ee], district_subset=district_subset,
                                                     cur_admins=cur_admins, pop_sizes = pop_sizes, min_year=min_year, max_year=max_year, burden_colname=burden_colname, age_plotted=age_plotted, plot_by_month=plot_by_month)
      if(nrow(burden_df)==0){
        burden_df = cur_sim_output_agg
      } else{
        burden_df = rbind(burden_df, cur_sim_output_agg)
      }
    }
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # Connect the 'to-present' and 'future-projection' simulations in the plot
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # Two alternatives for how this is done, controlled by extend_past_timeseries:
    #   - (FALSE) extend the 'future-projection' lines all back to the end of the 'to-present' simulations, which is desirable if the future projection scenarios separate right away
    #   - (TRUE) extend the end of the 'to-present' line up to the specified point in the 'future-projections' timeseries. This is only desirable if all 'future-projections' are 
    #            identical up to that point (e.g., 'to-present' simulations only run to 2020 and we are currently in 2023, so 2021-2022 are identical in all 'future projection' scenarios)
    if('to-present' %in% burden_df$scenario){
      connect_future_with_past = TRUE
      similarity_threshold = 0.15
      if(!is.na(extend_past_timeseries_year) & (extend_past_timeseries_year %in% burden_df$year[burden_df$scenario != 'to-present']) & time_string=='annual'){
        # check whether the future projections are all nearly identical (minus stochasticity) for the initial year (otherwise, use version that extends future-projection lines back to past)
        future_df = burden_df[burden_df$scenario != 'to-present',]
        earliest_future_year = min(future_df$year)
        compare_burdens = future_df$mean_burden[future_df$year == earliest_future_year]
        if(all(compare_burdens<(compare_burdens[1]*(1+similarity_threshold))) & all(compare_burdens>(compare_burdens[1]*(1-similarity_threshold)))){
          connect_future_with_past = FALSE
          merge_years = earliest_future_year
          if(extend_past_timeseries_year > earliest_future_year){
            # check which years (up to a maximum of extend_past_timeseries_year) should be included in the to-present line
            yy = earliest_future_year + 1
            while(yy <= extend_past_timeseries_year){
              compare_burdens = future_df$mean_burden[future_df$year == (yy)]
              if(all(compare_burdens<(compare_burdens[1]*1.05)) & all(compare_burdens>(compare_burdens[1]*0.95))){
                merge_years = c(merge_years, yy)
                yy = yy+1
              } else{  # as soon as they don't match for a year, stop trying to match any future years
                yy=99999999
              }
            }
          }
          # get the mean value from the 'future-projection' rows so that it can be added to the 'to-present' scenario
          past_from_future_df = future_df[future_df$year %in% merge_years,]
          past_from_future_df_means = past_from_future_df %>% dplyr::select(-scenario) %>% group_by(year) %>%
            summarise_all(mean) %>% ungroup()
          past_from_future_df_means$scenario = 'to-present'
          # delete the old 'future-projection' rows for all but the final of these years
          delete_future_years = merge_years[merge_years != max(merge_years)]
          if(length(delete_future_years)>0) burden_df = burden_df[-which(burden_df$year %in% merge_years),]
          # add the rows to the 'to-present' scenario in the data frame
          burden_df = merge(burden_df, past_from_future_df_means, all=TRUE)
        }else{
          connect_future_with_past = TRUE
        }
      } 
      if(connect_future_with_past){
        # add the final 'to-present' row to all future simulations for a continuous plot
        to_present_df = burden_df[burden_df$scenario == 'to-present',]
        if(plot_by_month){
          final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
          for(ss in 2:length(scenario_names)){
            final_to_present_row$scenario = scenario_names[ss]
            burden_df = rbind(burden_df, final_to_present_row)
          }
        } else{
          final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
          for(ss in 2:length(scenario_names)){
            final_to_present_row$scenario = scenario_names[ss]
            burden_df = rbind(burden_df, final_to_present_row)
          }
        }
      }
    }
      
    
    
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # if plotting burden relative to specified year, calculate relative values
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    if (!is.na(relative_year) & time_string=='annual'){
      # if the reference year is in the to-present simulation, use the same reference for all scenarios
      # if the reference year is not in the to-present simulation, use average value across scenarios and check that all scenarios have similar values for the reference year (if they do not, send a warning)
      if(('to-present' %in% burden_df$scenario) & (relative_year %in% unique(burden_df$year[burden_df$scenario=='to-present']))){
        # get the burden in the reference year
        reference_burden_cur = burden_df$mean_burden[burden_df$scenario=='to-present' & burden_df$year == relative_year]
      } else{
        similarity_threshold = 0.1
        all_ref_year_burdens = burden_df$mean_burden[burden_df$year == relative_year]
        reference_burden_cur = mean(all_ref_year_burdens)
        if(any(all_ref_year_burdens>(reference_burden_cur*(1+similarity_threshold))) | any(all_ref_year_burdens<(reference_burden_cur*(1-similarity_threshold)))){
          warning(paste0('in the reference year, some scenarios have different burdens for ',burden_metric_name))
        }
      }
      # calculate all relative burden values as mean_burden / reference_burden_cur: this will be referred to as 'burden relative to burden in relative_year'
      burden_df$mean_burden = burden_df$mean_burden / reference_burden_cur
      burden_df$max_burden = NA
      burden_df$min_burden = NA
        
      ylab_add_component = paste0('\n relative to ', relative_year)
      relative_string = paste0('_relativeTo', relative_year)
    } else{
      ylab_add_component = ''
      relative_string = ''
    }
      
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # create scenario-comparison plots
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    
    # subset to relevant scenarios currently being compared
    burden_df = burden_df[burden_df$scenario %in% scenario_names,]
    # get factors in the correct order (rather than alphabetical)
    burden_df$scenario = factor(burden_df$scenario, levels=rev(scenario_names))
    
    if(is.na(scenario_linetypes[1])){
      scenario_linetypes = rep(1, length(unique(burden_df$scenario)))
      names(scenario_linetypes) = unique(burden_df$scenario)
    }
    # ----- malaria burden ----- #
    
    if(plot_by_month){
      gg_list[[bb]] = ggplot(burden_df, aes(x=as.Date(date), y=mean_burden, color=scenario)) +
        geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
        scale_fill_manual(values = rev(scenario_palette)) + 
        geom_line(size=1) + 
        scale_color_manual(values = rev(scenario_palette)) + 
        xlab('date') + 
        ylab(paste0(gsub('\\(births\\)', '', burden_metric_name),ylab_add_component)) + 
        xlim(as.Date(paste0(min_year, '-01-01')), as.Date(paste0(max_year, '-01-01'))) +
        coord_cartesian(ylim=c(0, NA)) +
        theme_classic()+ 
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size))
    } else{
      gg_list[[bb]] = ggplot(burden_df, aes(x=year, y=mean_burden, color=scenario, linetype=scenario))  +
        geom_line(size=1) + 
        scale_linetype_manual(values=rev(scenario_linetypes)) +
        scale_color_manual(values = rev(scenario_palette)) + 
        xlab('year') + 
        ylab(paste0(gsub('\\(births\\)', '', burden_metric_name), ylab_add_component)) + 
        xlim(min_year, max_year) + 
        scale_x_continuous(breaks= pretty_breaks()) +
        coord_cartesian(ylim=c(0, NA)) +
        theme_classic()+ 
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size))
    }
    if(plot_CI){
      gg_list[[bb]] =  gg_list[[bb]] +
        geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
        scale_fill_manual(values = scenario_palette)
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] + theme(legend.position='none', text=element_text(size =separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_', relative_string, burden_metric_name,'_',district_subset,'.png'), separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }
  # gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]])), (ggplot() + theme_void())), gg_list)
  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  # remove legend from main plots
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position = "none")  + theme(text = element_text(size = text_size))   
  }
  # ----- combine all burden plots ----- #
  # gg = grid.arrange(grobs = gg_list, layout_matrix = matrix(c(1,1,2:(length(burden_colnames)+1)), ncol=2, byrow=TRUE))  # other orientation
  nrow_plot = 2
  if(include_U1) nrow_plot = 3
  num_in_matrix = ceiling(length(burden_colnames)/nrow_plot)*nrow_plot
  gg = grid.arrange(grobs = gg_list, layout_matrix = rbind(matrix(rep(1, ceiling(length(burden_colnames)/nrow_plot)), nrow=1), matrix(2:(num_in_matrix+1), nrow=nrow_plot, byrow=FALSE)))
  
  if(save_plots){
    ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_burden', relative_string,'_pyr', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg, dpi=600, width=9, height=3*nrow_plot, units='in')
  }
  
  return(gg)
}





######################################################################
# create state grid plot with timeseries of burden, no intervention info
######################################################################

plot_simulation_output_burden_by_state = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations,
                                             min_year, max_year, sim_end_years, relative_year=NA,
                                             scenario_filepaths, scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, 
                                             extend_past_timeseries_year=NA, scenario_linetypes=NA,filename_suffix=''){
  if (!is.na(relative_year)){ if(relative_year<min_year){
    warning('specified minimum year must be <= relative year. Setting min_year to relative_year.')
    min_year = relative_year
  }}
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))


  
  # iterate through scenarios, storing relevant output
  burden_df_all = data.frame()
  for(ee in 1:length(scenario_filepaths)){
    exp_filepath = scenario_filepaths[ee]
    exp_name = scenario_names[ee]
    cur_sim_output_agg = get_burden_timeseries_by_state(exp_filepath=exp_filepath, exp_name=exp_name, pop_filepath=pop_filepath, overwrite_files=overwrite_files)
    if(nrow(burden_df_all)==0){
      burden_df_all = cur_sim_output_agg
    } else{
      burden_df_all = rbind(burden_df_all, cur_sim_output_agg)
    }
  }
  
  # subset to relevant scenarios currently being compared
  burden_df_all = burden_df_all[burden_df_all$scenario %in% scenario_names,]

  # ----- malaria burden ----- #
  burden_metrics = c( 'PfPR_U5', 'PfPR_all', 'incidence_pp_U5', 'incidence_pp_all',  'direct_mortality_pp_U5', 'direct_mortality_pp_all',  'total_mortality_pp_U5', 'total_mortality_pp_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence per person (U5)', 'incidence per person (all ages)',  'direct mortality per person (U5)','direct mortality per person (all ages)', 'total mortality per person (U5)','total mortality per person (all ages)')
  
  for(bb in 1:length(burden_metrics)){
    burden_metric = burden_metrics[bb]
    burden_metric_name = burden_metric_names[bb]
    burden_df = burden_df_all[,c('State','year', 'scenario', burden_metric)]
    burden_df$mean_burden = burden_df[,burden_metric]

    # connect the 'to-present' and 'future-projection' simulations in the plot. Two alternatives for how this is done, controlled by extend_past_timeseries:
    #   - (FALSE) extend the 'future-projection' lines all back to the end of the 'to-present' simulations, which is desirable if the future projection scenarios separate right away
    #   - (TRUE) extend the end of the 'to-present' line up to the specified point in the 'future-projections' timeseries. This is only desirable if all 'future-projections' are 
    #            identical up to that point (e.g., 'to-present' simulations only run to 2020 and we are currently in 2023, so 2021-2022 are identical in all 'future projection' scenarios)
    if('to-present' %in% burden_df$scenario){
      connect_future_with_past = TRUE
      similarity_threshold = 0.15
      if(!is.na(extend_past_timeseries_year) & (extend_past_timeseries_year %in% burden_df$year[burden_df$scenario != 'to-present'])){
        # check whether the future projections are all nearly identical (minus stochasticity) for the initial year (otherwise, use version that extends future-projection lines back to past)
        future_df = burden_df[burden_df$scenario != 'to-present',]
        earliest_future_year = min(future_df$year)
        # compare_burdens = future_df$mean_burden[future_df$year == earliest_future_year]
        # if(all(compare_burdens<(compare_burdens[1]*(1+similarity_threshold))) & all(compare_burdens>(compare_burdens[1]*(1-similarity_threshold)))){
        #   connect_future_with_past = FALSE
        #   merge_years = earliest_future_year
        #   if(extend_past_timeseries_year > earliest_future_year){
        #     # check which years (up to a maximum of extend_past_timeseries_year) should be included in the to-present line
        #     yy = earliest_future_year + 1
        #     while(yy <= extend_past_timeseries_year){
        #       compare_burdens = future_df$mean_burden[future_df$year == (yy)]
        #       if(all(compare_burdens<(compare_burdens[1]*1.05)) & all(compare_burdens>(compare_burdens[1]*0.95))){
        #         merge_years = c(merge_years, yy)
        #         yy = yy+1
        #       } else{  # as soon as they don't match for a year, stop trying to match any future years
        #         yy=99999999
        #       }
        #     }
        #   }
        #   # get the mean value from the 'future-projection' rows so that it can be added to the 'to-present' scenario
        #   past_from_future_df = future_df[future_df$year %in% merge_years,]
        #   past_from_future_df_means = past_from_future_df %>% dplyr::select(-scenario) %>% group_by(year) %>%
        #     summarise_all(mean) %>% ungroup()
        #   past_from_future_df_means$scenario = 'to-present'
        #   # delete the old 'future-projection' rows for all but the final of these years
        #   delete_future_years = merge_years[merge_years != max(merge_years)]
        #   if(length(delete_future_years)>0) burden_df = burden_df[-which(burden_df$year %in% merge_years),]
        #   # add the rows to the 'to-present' scenario in the data frame
        #   burden_df = merge(burden_df, past_from_future_df_means, all=TRUE)
        # }else{
          connect_future_with_past = TRUE
        # }
      } 
      if(connect_future_with_past){
        # add the final 'to-present' row to all future simulations for a continuous plot
        to_present_df = burden_df[burden_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          burden_df = rbind(burden_df, final_to_present_row)
        }
      }
    }
    
    
    
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # if plotting burden relative to specified year, calculate relative values
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    if (!is.na(relative_year)){
      # if the reference year is in the to-present simulation, use the same reference for all scenarios
      # if the reference year is not in the to-present simulation, use average value across scenarios and check that all scenarios have similar values for the reference year (if they do not, send a warning)
      if(('to-present' %in% burden_df$scenario) & (relative_year %in% unique(burden_df$year[burden_df$scenario=='to-present']))){
        # get the burden in the reference year
        reference_burden_cur = burden_df[burden_df$scenario=='to-present' & burden_df$year == relative_year, c('State','mean_burden')]
      } else{
        similarity_threshold = 0.1
        all_ref_year_burdens = burden_df$mean_burden[burden_df$year == relative_year, c('State','mean_burden')]
        reference_burden_cur = all_ref_year_burdens %>% group_by(State) %>%
          summarise(mean_burden = mean(mean_burden)) %>% ungroup()
        # if(any(all_ref_year_burdens>(reference_burden_cur*(1+similarity_threshold))) | any(all_ref_year_burdens<(reference_burden_cur*(1-similarity_threshold)))){
        #   warning(paste0('in the reference year, some scenarios have different burdens for ',burden_metric_name))
        # }
      }
      # calculate all relative burden values as mean_burden / reference_burden_cur: this will be referred to as 'burden relative to burden in relative_year'
      colnames(reference_burden_cur)[colnames(reference_burden_cur)=='mean_burden'] = 'ref_mean_burden'
      burden_df = merge(burden_df, reference_burden_cur, all=TRUE)
      burden_df$mean_burden = burden_df$mean_burden / burden_df$ref_mean_burden
      
      burden_ylab = paste0(gsub(' per person', '', burden_metric_name), ' relative to ', relative_year)
      relative_string = paste0('_relativeTo', relative_year)
    } else{
      burden_ylab = burden_metric_name
      relative_string = ''
    }
    
      
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # create scenario-comparison plots
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    
    # get factors in the correct order (rather than alphabetical)
    burden_df$scenario = factor(burden_df$scenario, levels=rev(scenario_names))
    
    if(is.na(scenario_linetypes[1])){
      scenario_linetypes = rep(1, length(unique(burden_df$scenario)))
      names(scenario_linetypes) = unique(burden_df$scenario)
    }

    burden_df$code = burden_df$State
    gg = ggplot(burden_df, aes(x=year, y=mean_burden, color=scenario, linetype=scenario))+
      geom_line(size=1) + 
      scale_linetype_manual(values=scenario_linetypes) +
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(burden_ylab) + 
      coord_cartesian(ylim=c(0, ifelse(!is.na(relative_year),2,NA)), xlim=c(min_year, max_year)) + 
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name")#, scales='free') 

      ggsave(paste0(sim_future_output_dir, '/_plots/Timeseries_burden',relative_string,'_state_grid_',burden_metric,filename_suffix,'.png'), gg, dpi=600, width=12*0.7, height=10*0.7, units='in')
  }
  return(gg)
}









######################################################################
# create plot panel with selected burden metric and intervention info
######################################################################
# note: plot of ITN use rates through time is for the entire population (always shows all-age, even when burden plot shows U5)

plot_simulation_intervention_output = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                               plot_by_month, min_year, max_year, sim_end_years, 
                                               burden_metric, age_plotted, 
                                               pyr, chw_cov,
                                               scenario_filepaths, scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                               indoor_protection_fraction=0.75, LLIN2y_flag=FALSE, overwrite_files=FALSE){
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  pop_sizes = read.csv(pop_filepath)
  pop_sizes = pop_sizes[,c('admin_name','pop_size')]
  # if we include all admins, get list of names from population size dataframe
  if(cur_admins[1] == 'all'){
    cur_admins = unique(pop_sizes$admin_name)
  }
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  
  # ----- malaria burden ----- #
  
  # Get output column name for specified burden metric
  # Note: need to divide by pop size and multiply by 1000 if not PfPR
  burden_colname = NA
  if(burden_metric == 'PfPR'){
    if(age_plotted == 'U5'){
      burden_colname = 'PfPR_U5'
    } else if(age_plotted == 'all'){
      burden_colname = 'PfPR_MiP_adjusted'
    }
  } else if(burden_metric == 'incidence'){
    if(age_plotted == 'U5'){
      burden_colname = 'New_clinical_cases_U5'
    } else if(age_plotted == 'all'){
      burden_colname = 'New_Clinical_Cases'
    }
  } else if(burden_metric == 'directMortality'){
    if(age_plotted == 'U5'){
      burden_colname = 'direct_mortality_nonMiP_U5_mean'
    } else if(age_plotted == 'all'){
      burden_colname = 'direct_mortality_nonMiP_mean'
    }
  } else if(burden_metric == 'allMortality'){
    if(age_plotted == 'U5'){
      burden_colname = 'total_mortality_U5_mean'
    } else if(age_plotted == 'all'){
      burden_colname = 'total_mortality_mean'
    }
  } 
  if(is.na(burden_colname)){
    warning('PROBLEM DETECTED: name of burden metric or age group not currently supported')
  }
  
  # check whether burden output already exists for this comparison
  if(LLIN2y_flag){
    llin2y_string = '_2yLLIN'
  } else{
    llin2y_string = ''
  }
  burden_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_burden_',time_string,'Timeseries_', burden_metric, '_', age_plotted, '_pyr', pyr, '_', chw_cov, 'CHW_',district_subset, llin2y_string,'.csv')
  if(file.exists(burden_df_filepath) & !overwrite_files){
    burden_df = read.csv(burden_df_filepath)
  } else{
    # iterate through scenarios, storing relevant output
    burden_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_sim_output_agg = get_burden_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                     exp_name = scenario_names[ee],  district_subset=district_subset,
                                                     cur_admins=cur_admins, pop_sizes = pop_sizes, min_year=min_year, max_year=max_year, burden_colname=burden_colname, age_plotted=age_plotted, plot_by_month=plot_by_month)
      if(nrow(burden_df)==0){
        burden_df = cur_sim_output_agg
      } else{
        burden_df = merge(burden_df, cur_sim_output_agg, all=TRUE)
      }
    }
    
    # add the final 'to-present' row to all future simulations for a continuous plot
    to_present_df = burden_df[burden_df$scenario == 'to-present',]
    if(nrow(to_present_df)>0){
      if(plot_by_month){
        final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          burden_df = rbind(burden_df, final_to_present_row)
        }
      } else{
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          burden_df = rbind(burden_df, final_to_present_row)
        }
      }
    }
    write.csv(burden_df, burden_df_filepath, row.names=FALSE)
  }
  
  
  
  
  # ----- LLIN, vaccine, and IRS intervention coverage ----- #
  
  # check whether LLIN/IRS output already exists for this comparison
  llin_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_llin_irs_',time_string,'Timeseries', '_pyr', pyr, '_', chw_cov, 'CHW_',district_subset, llin2y_string,'.csv')
  if(file.exists(llin_df_filepath) & !overwrite_files){
    net_use_df = read.csv(llin_df_filepath)
  } else{
    # iterate through scenarios, storing relevant output
    net_use_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_net_agg = get_intervention_use_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                        exp_name = scenario_names[ee], 
                                                        cur_admins=cur_admins, pop_sizes=pop_sizes, min_year=min_year, max_year=max_year, indoor_protection_fraction=indoor_protection_fraction, plot_by_month=plot_by_month)
      if(nrow(net_use_df)==0){
        net_use_df = cur_net_agg
      } else{
        net_use_df = merge(net_use_df, cur_net_agg, all=TRUE)
      }
    }
    
    # first, remove the final 'to-present' month or year - it should was overwritten in the pick-up from burn-in
    # then, add the final 'to-present' row to all future simulations for a continuous plot
    if(plot_by_month){
      # remove excess month from to-present simulation
      max_to_present_date = max(net_use_df$date[net_use_df$scenario == 'to-present'])
      row_to_remove = intersect(which(net_use_df$scenario == 'to-present'), which(net_use_df$date == max_to_present_date))
      net_use_df = net_use_df[-row_to_remove,]
      
      # join past and future simulation trajectories
      to_present_df = net_use_df[net_use_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        net_use_df = rbind(net_use_df, final_to_present_row)
      }
    } else{
      # remove excess year from to-present simulation
      if(any(net_use_df$scenario == 'to-present')){
        max_to_present_date = max(net_use_df$year[net_use_df$scenario == 'to-present'])
        row_to_remove = intersect(which(net_use_df$scenario == 'to-present'), which(net_use_df$year == max_to_present_date))
        net_use_df = net_use_df[-row_to_remove,]
        
        # join past and future simulation trajectories
        to_present_df = net_use_df[net_use_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          net_use_df = rbind(net_use_df, final_to_present_row)
        }
      }
    }
    write.csv(net_use_df, llin_df_filepath, row.names=FALSE)
  }
  
  
  
  
  # ----- Case management ----- #
  
  # check whether CM output already exists for this comparison
  cm_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_cm_',time_string,'Timeseries', '_pyr', pyr, '_', chw_cov, 'CHW_',district_subset, llin2y_string,'.csv')
  if(file.exists(cm_df_filepath)){
    cm_df = read.csv(cm_df_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    cm_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      cm_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      
      cur_cm_agg = get_cm_timeseries_exp(cm_filepath=cm_filepath, pop_sizes=pop_sizes, end_year=end_year, exp_name = scenario_names[ee], 
                                         cur_admins=cur_admins, min_year=min_year, plot_by_month=plot_by_month)
      
      if(nrow(cm_df)==0){
        cm_df = cur_cm_agg
      } else{
        cm_df = rbind(cm_df, cur_cm_agg)
      }
    }
    
    # add the final 'to-present' row to all future simulations for a continuous plot
    if(plot_by_month){
      # join past and future simulation trajectories
      to_present_df = cm_df[cm_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        cm_df = rbind(cm_df, final_to_present_row)
      }
    } else{
      # join past and future simulation trajectories
      to_present_df = cm_df[cm_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        cm_df = rbind(cm_df, final_to_present_row)
      }
    }
    write.csv(cm_df, cm_df_filepath, row.names=FALSE)
  }
  
  
  
  
  
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  burden_df$scenario = factor(burden_df$scenario, levels=scenario_names)
  
  # ----- malaria burden ----- #
  
  if(plot_by_month){
    g_burden = ggplot(burden_df, aes(x=as.Date(date), y=mean_burden, color=scenario)) +
      geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0(burden_metric, ' - ', age_plotted)) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size))
  } else{
    g_burden = ggplot(burden_df, aes(x=year, y=mean_burden, color=scenario)) +
      geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0(burden_metric, ' - ', age_plotted)) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size))
  }
  
  inter_plot_list = list()
  # ----- LLIN use and distribution ----- #
  # plot net use through time
  if(plot_by_month){
    g_net_use = ggplot(net_use_df, aes(x=as.Date(date), y=coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('LLIN use (all ages)')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  } else{
    g_net_use = ggplot(net_use_df, aes(x=year, y=coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      # geom_hline(yintercept=0.22, alpha=0.1)+
      # geom_hline(yintercept=0.39, alpha=0.1)+
      xlab('year') + 
      ylab(paste0('LLIN use (all ages)')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
    
    g_all_inter = ggplot() +
      geom_line(data=net_use_df, aes(x=year, y=coverage), color=rgb(1,0.6,1), size=1) + 
      xlab('year') + 
      ylab(paste0('coverage metric')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  }
  inter_plot_list = append(inter_plot_list, list(g_net_use))
  # # plot net distribution numbers through time (how many nets distributed in each month or year per person?)
  # if(plot_by_month){
  #   g_net_dist = ggplot(net_use_df, aes(x=as.Date(date), y=new_net_per_cap, color=scenario)) +
  #     geom_point(size=1) + 
  #     scale_color_manual(values = scenario_palette) + 
  #     xlab('date') + 
  #     ylab(paste0('LLINs distributed per person')) + 
  #     theme_classic()+ 
  #     theme(legend.position = "none", text = element_text(size = text_size))
  # } else{
  #   g_net_dist = ggplot(net_use_df, aes(x=year, y=new_net_per_cap, color=scenario)) +
  #     geom_point(size=2) + 
  #     geom_line(alpha=0.2, size=2) +
  #     scale_color_manual(values = scenario_palette) + 
  #     xlab('year') + 
  #     ylab(paste0('LLINs distributed per person')) + 
  #     theme_classic()+ 
  #     theme(legend.position = "none", text = element_text(size = text_size))
  # }
  
  
  # ----- vaccine ----- #
  if('vacc_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_vacc = ggplot(net_use_df, aes(x=as.Date(date), y=vacc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('Vaccines (primary series + booster) per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_vacc = ggplot(net_use_df, aes(x=year, y=vacc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, size=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Vaccines (primary series + booster) per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=vacc_per_cap), color=rgb(0,0.3,0), size=1)
        
    }
    inter_plot_list = append(inter_plot_list, list(g_vacc))
  }
  
  
  # ----- PMC ----- #
  if('pmc_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_pmc = ggplot(net_use_df, aes(x=as.Date(date), y=pmc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('PMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_pmc = ggplot(net_use_df, aes(x=year, y=pmc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, size=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('PMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=pmc_per_cap), color=rgb(0.0,0.4,1), size=1)
    }
    inter_plot_list = append(inter_plot_list, list(g_pmc))
  }
  
  
  # ----- SMC ----- #
  if('smc_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_smc = ggplot(net_use_df, aes(x=as.Date(date), y=smc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('SMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_smc = ggplot(net_use_df, aes(x=year, y=smc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, size=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=smc_per_cap), color=rgb(0.0,0.4,1), size=1)
    }
    inter_plot_list = append(inter_plot_list, list(g_smc))
  }
  
  
  # ----- IRS ----- #
  if('irs_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_irs = ggplot(net_use_df, aes(x=as.Date(date), y=irs_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('IRS rounds per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_irs = ggplot(net_use_df, aes(x=year, y=irs_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, size=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('IRS per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=irs_per_cap), color=rgb(1,0,1), size=1)
    }
    inter_plot_list = append(inter_plot_list, list(g_irs))
  }
  
  
  # ----- Case management ----- #
  if(plot_by_month){
    g_cm = ggplot(cm_df, aes(x=as.Date(date), y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  } else{
    g_cm = ggplot(cm_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
    
    g_all_inter = g_all_inter +
      geom_line(data=cm_df, aes(x=year, y=mean_coverage), color=rgb(0.1,0.9,0.4), size=1)
  }
  inter_plot_list = append(inter_plot_list, list(g_cm))
  
  
  
  # ----- combine burden and intervention plots ----- #
  gg_leg = ggpubr::as_ggplot(ggpubr::get_legend(g_burden))
  g_burden = g_burden + theme(legend.position = "none")
  gg = plot_grid(plotlist=append(list(gg_leg, g_burden), inter_plot_list), ncol=1, nrow=(2+length(inter_plot_list)), align='vh', axis='lrtb')  # (gg_leg, g_burden, plot_list)
  
  if(save_plots){
    ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_', burden_metric, '_', age_plotted, '_versusInterventions_pyr', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg, dpi=600, width=7, height=4*(2+length(inter_plot_list)), units='in')
    
    if(!plot_by_month){
      ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_interventions_pyr', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), g_all_inter, dpi=600, width=7, height=5, units='in')
    }
  }
  return(gg)
}





######################################################################
# create plot panel with CM timeseries for included scenarios in each state
######################################################################

plot_state_grid_cm = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                               plot_by_month, min_year, max_year, sim_end_years, 
                               scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                               overwrite_files=FALSE){

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]

  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  # check whether CM output already exists for this comparison
  cm_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_cm_state_',time_string,'Timeseries.csv')
  if(file.exists(cm_df_filepath)){
    cm_df = read.csv(cm_df_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    cm_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      cm_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      
      cur_cm_agg = get_cm_timeseries_by_state(cm_filepath=cm_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                              min_year=min_year, plot_by_month=plot_by_month)
      
      if(nrow(cm_df)==0){
        cm_df = cur_cm_agg
      } else{
        cm_df = rbind(cm_df, cur_cm_agg)
      }
    }
    
    if(any(grepl('to-present', cm_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      if(plot_by_month){
        # join past and future simulation trajectories
        to_present_df = cm_df[cm_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          cm_df = rbind(cm_df, final_to_present_row)
        }
      } else{
        # join past and future simulation trajectories
        to_present_df = cm_df[cm_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          cm_df = rbind(cm_df, final_to_present_row)
        }
      }
    }
    write.csv(cm_df, cm_df_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  cm_df$scenario = factor(cm_df$scenario, levels=rev(scenario_names))
  cm_df$code = cm_df$State
  
  if(plot_by_month){
    g_cm = ggplot(cm_df, aes(x=as.Date(date), y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year))+
      theme_bw()+ 
      theme(legend.position = "none", text = element_text(size = text_size))+
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  } else{
    g_cm = ggplot(cm_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_CM_by_state.png'), g_cm, dpi=600, width=12, height=10, units='in')
}














#####################################################################
# plot map of admin subsets
#####################################################################
plot_included_admin_map = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, admin_shapefile_filepath, shapefile_admin_colname){
  admin_pop = read.csv(pop_filepath)
  admin_shapefile = st_read(admin_shapefile_filepath)
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=admin_pop$admin_name, origin_names=admin_shapefile[[shapefile_admin_colname]])
  
  admin_in_map = data.frame(admin_name = admin_pop$admin_name, admin_included='no')
  admin_in_map$admin_included[admin_in_map$admin_name %in% cur_admins] = 'yes'
  included_colors = c('#006692', 'grey96')
  names(included_colors) = c('yes', 'no')

  admin_cur = admin_shapefile %>%
    dplyr::left_join(admin_in_map, by=c('NOMDEP' = 'admin_name'))
  
  gg_map = ggplot(admin_cur) +
    geom_sf(aes(fill=admin_included), size=0.5, color='black') +
    scale_fill_manual(values=included_colors, drop=FALSE, na.value='grey96') + 
    theme_void() +
    theme(legend.position = 'none') 
  ggsave(paste0(sim_future_output_dir, '/_plots/map_admins_included_', district_subset, '.png'), gg_map, dpi=600, width=4.8, height=4.8, units='in')
}



#####################################################################
# plot maps of burden with and without the intervention
#####################################################################
plot_burden_maps = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                            barplot_start_year, barplot_end_year,
                            pyr, chw_cov,
                            scenario_names, experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', LLIN2y_flag=FALSE,
                            filename_sufffix='', overwrite_files=FALSE){


  admin_pop = read.csv(pop_filepath)
  if(!(cur_admins[1] == 'all')){
    admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),]
  }
  admin_shapefile = shapefile(admin_shapefile_filepath)
  # standardize shapefile names
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=archetype_info$LGA, origin_names=admin_shapefile$NOMDEP)
  
  years_included = barplot_end_year - barplot_start_year + 1

  # burden metrics
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'death_rate_mean_U5', 'death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  # burden_colnames_for_map = c('pfpr_all', 'incidence_all', 'mortality_rate_all')
  # burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'mortality (all)')
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  burden_colnames_for_map = c('pfpr_u5', 'pfpr_all', 'incidence_u5', 'incidence_all', 'mortality_rate_u5',  'mortality_rate_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #   iterate through scenarios, creating dataframe including all burden metrics
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  num_scenarios = length(experiment_names)
  burden_df_all = data.frame()
  for(ee in 1:num_scenarios){
    experiment_name = experiment_names[ee]
    cur_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop, comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
    cur_burden_df$scenario_name = scenario_names[ee]
    if(nrow(burden_df_all) == 0){
      burden_df_all = cur_burden_df
    } else{
      burden_df_all = rbind(burden_df_all, cur_burden_df)
    }
  }



  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #      create maps showing each burden metric for all scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  if(LLIN2y_flag){
    llin2y_string = '_2yLLIN'
  } else{
    llin2y_string = ''
  }
  num_colors = 40
  # colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
  colorscale = colorRampPalette(brewer.pal(9, 'YlOrRd'))(num_colors)
  


  # iterate through burden metrics, creating plots for each
  for(cc in 1:length(burden_colnames_for_map)){

    if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_', burden_colnames_for_map[cc], '_', pyr, '_', chw_cov, 'CHW_', district_subset, llin2y_string, filename_sufffix, '.png'), res=600, width=(num_scenarios*3+2)*3/4, height=3, units='in')
    par(mar=c(0,1,2,0))
    # set layout for panel of maps
    layout_matrix = matrix(rep(c(rep(1:num_scenarios, each=3),rep((num_scenarios+1),2)),2), nrow=2, byrow=TRUE)
    layout(mat = layout_matrix)

    cur_colname = burden_colnames_for_map[cc]
    min_value = min(min(burden_df_all[[cur_colname]], na.rm=TRUE), 0)
    max_value = max(max(burden_df_all[[cur_colname]], na.rm=TRUE), 0.65)

    # iterate through scenarios
    for(ee in 1:num_scenarios){
      cur_burden_df = burden_df_all[burden_df_all$scenario_name == scenario_names[ee],]
      vals_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
      for (i_ds in 1:length(vals_ordered$ds_ordered)){
        cur_ds = vals_ordered$ds_ordered[i_ds]
        if(toupper(cur_ds) %in% toupper(cur_burden_df$admin_name)){
          vals_ordered$value[i_ds] = cur_burden_df[which(toupper(cur_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
        }
      }

      col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
      col_cur[is.na(col_cur)] = 'grey'
      plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=scenario_names[ee])
    }
    # legend
    legend_label_vals = seq(min_value, max_value, length.out=5)
    legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
    text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
    rasterImage(legend_image, 0, 0, 1,1)
    # fourth blank plot
    # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
    par(mfrow=c(1,1), mar=c(5,4,4,2))
    if(save_plots) dev.off()
  }

}





# plot map with the reduction in burden relative to the first scenario
plot_burden_relative_reduction_maps = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                            barplot_start_year, barplot_end_year,
                            pyr, chw_cov,
                            scenario_names, experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', LLIN2y_flag=FALSE,
                            overwrite_files=FALSE){
  
  
  admin_pop = read.csv(pop_filepath)
  if(!(cur_admins[1] == 'all')){
    admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),]
  }
  admin_shapefile = shapefile(admin_shapefile_filepath)
  # standardize shapefile names
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=archetype_info$LGA, origin_names=admin_shapefile$NOMDEP)
  
  years_included = barplot_end_year - barplot_start_year + 1
  
  # burden metrics
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'death_rate_mean_U5', 'death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  # burden_colnames_for_map = c('pfpr_all', 'incidence_all', 'mortality_rate_all')
  # burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'mortality (all)')
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  burden_colnames_for_map = c('pfpr_u5', 'pfpr_all', 'incidence_u5', 'incidence_all', 'mortality_rate_u5',  'mortality_rate_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #   iterate through scenarios, creating dataframe including all burden metrics
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  num_scenarios = length(experiment_names)
  burden_df_all = data.frame()
  for(ee in 1:num_scenarios){
    experiment_name = experiment_names[ee]
    cur_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop, comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
    cur_burden_df$scenario_name = scenario_names[ee]
    if(nrow(burden_df_all) == 0){
      burden_df_all = cur_burden_df
    } else{
      burden_df_all = rbind(burden_df_all, cur_burden_df)
    }
  }
  
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #      create maps showing each burden metric for all scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  if(LLIN2y_flag){
    llin2y_string = '_2yLLIN'
  } else{
    llin2y_string = ''
  }
  num_colors = 40
  colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
  
  
  
  # iterate through burden metrics, creating plots for each
  for(cc in 1:length(burden_colnames_for_map)){
    
    if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_rel_reduction_', burden_colnames_for_map[cc], '_', pyr, '_', chw_cov, 'CHW_', district_subset, llin2y_string, '.png'), res=600, width=(num_scenarios*3+2)*3/4, height=3, units='in')
    par(mar=c(0,1,2,0))
    # set layout for panel of maps
    layout_matrix = matrix(rep(c(rep(1:num_scenarios, each=3),rep((num_scenarios+1),2)),2), nrow=2, byrow=TRUE)
    layout(mat = layout_matrix)
    
    cur_colname = burden_colnames_for_map[cc]
    
    # get reference column for this burden metric
    burden_df_ref = burden_df_all[burden_df_all$scenario_name==scenario_names[1],c('admin_name',cur_colname)]
    colnames(burden_df_ref)[which(colnames(burden_df_ref)==cur_colname)]='reference_value'
    burden_df_relative = merge(burden_df_all, burden_df_ref, all=TRUE)
    burden_df_relative$rel_reduction = (burden_df_relative$reference_value - burden_df_relative[[cur_colname]]) / burden_df_relative$reference_value
    
    # set minimum and maximum plotted in legend
    abs_max = max(abs( burden_df_relative$rel_reduction), na.rm=TRUE)
    min_value = -1 * abs_max
    max_value = abs_max
    
    # iterate through scenarios
    for(ee in 2:num_scenarios){
      cur_burden_df = burden_df_relative[burden_df_relative$scenario_name == scenario_names[ee],]
      vals_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
      for (i_ds in 1:length(vals_ordered$ds_ordered)){
        cur_ds = vals_ordered$ds_ordered[i_ds]
        if(toupper(cur_ds) %in% toupper(cur_burden_df$admin_name)){
          vals_ordered$value[i_ds] = cur_burden_df[which(toupper(cur_burden_df$admin_name) == toupper(cur_ds)), 'rel_reduction']
        }
      }
      
      col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
      col_cur[is.na(col_cur)] = 'grey'
      plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=scenario_names[ee])
    }
    # legend
    legend_label_vals = seq(min_value, max_value, length.out=5)
    legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
    text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
    rasterImage(legend_image, 0, 0, 1,1)
    # fourth blank plot
    # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
    par(mfrow=c(1,1), mar=c(5,4,4,2))
    if(save_plots) dev.off()
  }
  
}
# 
# 
# 
# 
# 
# #####################################################################
# # plot maps of burden with and without IPTi
# #####################################################################
# 
# line2user <- function(line, side) {
#   lh <- par('cin')[2] * par('cex') * par('lheight')
#   x_off <- diff(grconvertX(0:1, 'inches', 'user'))
#   y_off <- diff(grconvertY(0:1, 'inches', 'user'))
#   switch(side,
#          `1` = par('usr')[3] - line * y_off * lh,
#          `2` = par('usr')[1] - line * x_off * lh,
#          `3` = par('usr')[4] + line * y_off * lh,
#          `4` = par('usr')[2] + line * x_off * lh,
#          stop("side must be 1, 2, 3, or 4", call.=FALSE))
# }
# 
# 
# 
# plot_IPTi_burden_maps = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
#                                  barplot_start_year, barplot_end_year, 
#                                  experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', overwrite_files=FALSE){
#   
#   admin_pop = read.csv(pop_filepath)
#   admin_shapefile = shapefile(admin_shapefile_filepath)
#   
#   years_included = barplot_end_year - barplot_start_year + 1 
#   
#   # burden metrics
#   burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'mortality (U1)')
#   burden_colnames_for_map = c('pfpr_u1', 'incidence_u1', 'mortality_rate_u1')
#   
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      read in and format malaria burden simulation output in IPTi admins
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   experiment_names_descriptions = c('noIPTi', 'IPTi')
#   pfpr_u1_df = data.frame(admin_name=cur_admins)
#   deaths_u1_df = data.frame(admin_name=cur_admins)
#   clinical_cases_u1_df = data.frame(admin_name=cur_admins)
#   pop_u1_df = data.frame(admin_name=admin_pop$admin_name)
#   
#   # no-IPTi burden df
#   experiment_name = experiment_names[1]
#   option_name = experiment_names_descriptions[1]
#   noIPTi_burden_df = get_total_U1_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   # IPTi burden df
#   experiment_name = experiment_names[2]
#   option_name = experiment_names_descriptions[2]
#   IPTi_burden_df = get_total_U1_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   
#   
#   
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      create panel of maps showing all burden metrics
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_IPTi_burden_', pyr, '_', chw_cov, 'CHW.png'), res=600, width=6, height=3*length(burden_metric_names), units='in')
#   par(mar=c(0,1,2,0))
#   
#   num_colors = 40
#   colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
#   
#   # set layout for panel of maps
#   base_matrix = matrix(c(1,1,1,2,2,2,3,3, 1,1,1,2,2,2,3,3), nrow=2, byrow=TRUE)
#   # add rows for each burden metric
#   layout_matrix = base_matrix
#   for(cc in 2:length(burden_colnames_for_map)){
#     layout_matrix = rbind(layout_matrix, base_matrix + 3*(cc-1))
#   }
#   # add row for title
#   layout_matrix = layout_matrix + 1
#   layout_matrix = rbind(rep(1, ncol(layout_matrix)), layout_matrix)
#   layout(mat = layout_matrix)
#   
#   # title
#   plot.new()
#   text(0.5,0.5,"Malaria burden in each health district",cex=2.5,font=1)
#   # text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
#   #      line2user(line=4, side=3), "Malaria burden in each health district", xpd=NA, cex=2, font=2)
#   
#   # iterate through burden metrics, creating plots of each
#   for(cc in 1:length(burden_colnames_for_map)){
#     cur_colname = burden_colnames_for_map[cc]
#     vals_ordered_noipti = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     vals_ordered_ipti = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     for (i_ds in 1:length(vals_ordered_noipti$ds_ordered)){
#       cur_ds = vals_ordered_noipti$ds_ordered[i_ds]
#       if(toupper(cur_ds) %in% toupper(noIPTi_burden_df$admin_name)){
#         vals_ordered_noipti$value[i_ds] = noIPTi_burden_df[which(toupper(noIPTi_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#         vals_ordered_ipti$value[i_ds] = IPTi_burden_df[which(toupper(IPTi_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#       }
#     }
#     min_value = min(c(vals_ordered_noipti$value, vals_ordered_ipti$value), na.rm=TRUE)
#     max_value = max(c(vals_ordered_noipti$value, vals_ordered_ipti$value), na.rm=TRUE)
#     # without IPTi
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_noipti$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - without IPTi'))
#     # with IPTi
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_ipti$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - with IPTi'))
#     # legend
#     legend_label_vals = seq(min_value, max_value, length.out=5)
#     legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
#     plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
#     text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
#     rasterImage(legend_image, 0, 0, 1,1)
#     # fourth blank plot
#     # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
#     
#   }
#   par(mfrow=c(1,1), mar=c(5,4,4,2))
#   if(save_plots) dev.off()
#   
#   
# }
# 
# 
# 
# #####################################################################
# # plot maps of burden with and without the intervention
# #####################################################################
# 
# plot_burden_maps_with_without_inter = function(sim_future_output_dir, pop_filepath, cur_admins, district_subset,
#                                                barplot_start_year, barplot_end_year,
#                                                experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', inter_name='PBO',
#                                                overwrite_files=FALSE){
#   
#   admin_pop = read.csv(pop_filepath)
#   admin_shapefile = shapefile(admin_shapefile_filepath)
#   
#   years_included = barplot_end_year - barplot_start_year + 1
#   
#   # burden metrics
#   # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'death_rate_mean_U5', 'death_rate_mean_all')
#   # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
#   burden_colnames_for_map = c('pfpr_all', 'incidence_all', 'mortality_rate_all')
#   burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'mortality (all)')
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      read in and format malaria burden simulation output in specified set of admins
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   experiment_names_descriptions = paste0(c('no', ''), inter_name)
#   
#   # no-intervention burden df
#   experiment_name = experiment_names[1]
#   option_name = experiment_names_descriptions[1]
#   noPBO_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   # with intervention burden df
#   experiment_name = experiment_names[2]
#   option_name = experiment_names_descriptions[2]
#   PBO_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   
#   # # increase in burden in each admin: (without intervention - with intervention)  /  with intervention
#   # rel_burd_increase = (noPBO_burden_df$incidence_all - PBO_burden_df$incidence_all ) / PBO_burden_df$incidence_all
#   # min(rel_burd_increase)
#   # max(rel_burd_increase)
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      create panel of maps showing all burden metrics
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_withWithout', inter_name,'_burden_', pyr, '_', chw_cov, 'CHW.png'), res=600, width=6, height=3*length(burden_metric_names), units='in')
#   par(mar=c(0,1,2,0))
#   
#   num_colors = 40
#   colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
#   
#   # set layout for panel of maps
#   base_matrix = matrix(c(1,1,1,2,2,2,3,3, 1,1,1,2,2,2,3,3), nrow=2, byrow=TRUE)
#   # add rows for each burden metric
#   layout_matrix = base_matrix
#   for(cc in 2:length(burden_colnames_for_map)){
#     layout_matrix = rbind(layout_matrix, base_matrix + 3*(cc-1))
#   }
#   # add row for title
#   layout_matrix = layout_matrix + 1
#   layout_matrix = rbind(rep(1, ncol(layout_matrix)), layout_matrix)
#   layout(mat = layout_matrix)
#   
#   # title
#   plot.new()
#   text(0.5,0.5,"Malaria burden in each health district",cex=2.5,font=1)
#   # text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
#   #      line2user(line=4, side=3), "Malaria burden in each health district", xpd=NA, cex=2, font=2)
#   
#   # iterate through burden metrics, creating plots of each
#   for(cc in 1:length(burden_colnames_for_map)){
#     cur_colname = burden_colnames_for_map[cc]
#     vals_ordered_noPBO = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     vals_ordered_PBO = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     for (i_ds in 1:length(vals_ordered_noPBO$ds_ordered)){
#       cur_ds = vals_ordered_noPBO$ds_ordered[i_ds]
#       if(toupper(cur_ds) %in% toupper(noPBO_burden_df$admin_name)){
#         vals_ordered_noPBO$value[i_ds] = noPBO_burden_df[which(toupper(noPBO_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#         vals_ordered_PBO$value[i_ds] = PBO_burden_df[which(toupper(PBO_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#       }
#     }
#     min_value = min(c(vals_ordered_noPBO$value, vals_ordered_PBO$value), na.rm=TRUE)
#     max_value = max(c(vals_ordered_noPBO$value, vals_ordered_PBO$value), na.rm=TRUE)
#     # without intervention
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_noPBO$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - without ', inter_name))
#     # with intervention
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_PBO$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - with ', inter_name))
#     # legend
#     legend_label_vals = seq(min_value, max_value, length.out=5)
#     legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
#     plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
#     text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
#     rasterImage(legend_image, 0, 0, 1,1)
#     # fourth blank plot
#     # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
#     
#   }
#   par(mfrow=c(1,1), mar=c(5,4,4,2))
#   if(save_plots) dev.off()
#   
#   
# }








