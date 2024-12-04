# plot_sim_input_functions.R

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




#####################################################################################################
# ================================================================================================= #
# create state-grid-arranged plot panel with intervention timeseries for included scenarios
# ================================================================================================= #
#####################################################################################################

##############
# CM
##############
plot_state_grid_cm = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                              min_year, max_year, sim_end_years, 
                              scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                              separate_admin_lines_flag = FALSE,  act_adherence_effective_multiplier=1, overwrite_files=FALSE){

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]

  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether CM output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_cm_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(file.exists(timeseries_filepath)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_cm_timeseries_by_state(cm_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_cm_timeseries_by_state(cm_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                        min_year=min_year, get_state_level=TRUE)
      }

      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  # adjust for ACT adherence to get coverage of individuals taking ACTs (not reduced effective coverage due to not completing the regimen)
  if(act_adherence_effective_multiplier<1){
    timeseries_df$mean_coverage = timeseries_df$mean_coverage / act_adherence_effective_multiplier
    timeseries_df$min_coverage = timeseries_df$min_coverage / act_adherence_effective_multiplier
    timeseries_df$max_coverage = timeseries_df$max_coverage / act_adherence_effective_multiplier
  }
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(aes(group_by=admin_name), size=0.8) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Effective treatment rate (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 

  } else{
      gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
        scale_fill_manual(values = scenario_palette) + 
        geom_line(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Effective treatment rate (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_CM_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}



##############
# ITN ANC
##############
plot_state_grid_itn_anc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                    min_year, max_year, sim_end_years, 
                                    scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                    separate_admin_lines_flag = FALSE,  act_adherence_effective_multiplier=1, overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether CM output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_itn_anc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(file.exists(timeseries_filepath)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input ITN ANC coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$ANC_ITN_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_itn_anc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                             min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_itn_anc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                             min_year=min_year, get_state_level=TRUE)
      }

      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }

  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group_by=admin_name), size=0.5) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ANC ITN coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ANC ITN coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_ITN_ANC_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}






##############
# SMC
##############

plot_state_grid_smc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                               min_year, max_year, sim_end_years, 
                               scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                               separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_smc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(file.exists(timeseries_filepath)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$SMC_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_smc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_smc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(aes(group_by=admin_name), size=0.8) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC coverage (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
        # scale_fill_manual(values = scenario_palette) + 
        geom_line(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC coverage (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_SMC_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}





##############
# ITN mass campaign
##############

plot_state_grid_itn = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                               min_year, max_year, sim_end_years, 
                               scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                               separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_itn_mass_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(file.exists(timeseries_filepath)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$ITN_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_itn_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_itn_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group_by=admin_name), size=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ITN initial use rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      # scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ITN initial use rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_ITN_mass_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}








##############
# Vaccines
##############

plot_state_grid_vacc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                min_year, max_year, sim_end_years, 
                                scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_vacc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(file.exists(timeseries_filepath)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$vacc_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_vacc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_vacc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group_by=admin_name), size=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('EPI malaria vaccine coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      # scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('EPI malaria vaccine coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_vacc_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}







##############
# PMC
##############

plot_state_grid_pmc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                min_year, max_year, sim_end_years, 
                                scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_pmc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(file.exists(timeseries_filepath)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$PMC_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_pmc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                          min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_pmc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                          min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group_by=admin_name), size=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('PMC coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      # scale_fill_manual(values = scenario_palette) + 
      geom_line(size=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('PMC coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_PMC_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
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


