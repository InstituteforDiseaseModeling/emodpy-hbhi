# process_sim_input_functions.R


library(data.table)
library(dplyr)
library(tidyverse)
library(lubridate)
# library(rgdal)
library(lubridate)




#####################################################################################################
# ================================================================================================= #
# calculate state- or LGA-level intervention  timeseries for included scenarios, based on inputs
# ================================================================================================= #
#####################################################################################################

##############
# CM
##############

# ----- admin- or state-level CM intervention coverage ----- #
get_cm_timeseries_by_state = function(cm_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(cm_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  input_df = merge(input_df, admin_info, by='admin_name')
  
  # CM is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_cm_years = unique(input_df$year)
    for(yy in cur_cm_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  if(any(input_df$duration==-1) & (max(input_df$year)<end_year)){
    cm_repeated = input_df[input_df$duration == -1,]
    for(rr in 1:(end_year - cm_repeated$year[1])){
      temp_year = cm_repeated
      temp_year$year = cm_repeated$year + rr
      temp_year$simday = cm_repeated$simday + rr*365
      input_df = rbind(input_df, temp_year)
    }
  }
  
  if(get_state_level){
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise_all(mean) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    # get population-weighted CM coverage across admins
    input_df$multiplied_U5_cm = input_df$U5_coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied CM coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cm, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cm / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    cm_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise_all(mean) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    input_df_agg_admin = input_df
    
    # take average across simulation seeds
    cm_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  cm_agg$scenario = exp_name
  return(cm_agg)
}


# ----- national-level CM intervention coverage ----- #
get_cm_timeseries_exp = function(cm_filepath, pop_sizes, end_year, exp_name, cur_admins, min_year, plot_by_month=TRUE){
  
  input_df = read.csv(cm_filepath)
  # subset to appropriate admins
  input_df = input_df[input_df$admin_name %in% cur_admins,] 
  if(!('seed' %in% input_df)) input_df$seed = 1
  
  # CM is sometimes repeated for several years but only listed once; change to repeate the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_cm_years = unique(input_df$year)
    for(yy in cur_cm_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  if(any(input_df$duration==-1) & (max(input_df$year)<end_year)){
    cm_repeated = input_df[input_df$duration == -1,]
    for(rr in 1:(end_year - cm_repeated$year[1])){
      temp_year = cm_repeated
      temp_year$year = cm_repeated$year + rr
      temp_year$simday = cm_repeated$simday + rr*365
      input_df = rbind(input_df, temp_year)
    }
  }
  # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
  input_df <- input_df %>% group_by(year, admin_name, seed) %>%
    summarise_all(mean) %>% ungroup()
  
  input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
  
  
  # get population-weighted CM coverage across admins
  input_df = merge(input_df, pop_sizes, by='admin_name')
  input_df$multiplied_U5_cm = input_df$U5_coverage * input_df$pop_size
  
  # get sum of population sizes and multiplied CM coverage across included admins
  input_df_agg_admin <- input_df %>% dplyr::select(year, seed, multiplied_U5_cm, pop_size) %>% group_by(year, seed) %>%
    summarise_all(sum) %>% ungroup()
  # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
  input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cm / input_df_agg_admin$pop_size
  
  
  # take average, max, and min burdens across simulation seeds
  if(plot_by_month){
    # subdivide year values and add dates
    # date dataframe
    included_years = unique(input_df_agg_admin$year)
    all_months = as.Date(paste0(rep(included_years, each=12),'-',c('01','02','03','04','05','06','07','08','09','10','11','12'), '-01' ))
    date_df = data.frame(year=rep(included_years, each=12), date=all_months)
    input_df_agg_admin_monthly = merge(input_df_agg_admin, date_df, by='year', all.x=TRUE, all.y=TRUE)
    cm_agg = as.data.frame(input_df_agg_admin_monthly) %>% dplyr::select(year, date, U5_coverage) %>%
      dplyr::group_by(year, date) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    cm_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, U5_coverage) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  
  cm_agg$scenario = exp_name
  return(cm_agg)
  
}




##############
# ANC ITN
##############

# ----- admin- or state-level CM intervention coverage ----- #
get_itn_anc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  input_df = merge(input_df, admin_info, by='admin_name')
  
  # ANC ITN is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_years = unique(input_df$year)
    for(yy in cur_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  if(any(input_df$duration==-1) & (max(input_df$year)<end_year)){
    df_repeated = input_df[input_df$duration == -1,]
    for(rr in 1:(end_year - df_repeated$year[1])){
      temp_year = df_repeated
      temp_year$year = df_repeated$year + rr
      temp_year$simday = df_repeated$simday + rr*365
      input_df = rbind(input_df, temp_year)
    }
  }
  
  if(get_state_level){
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise(coverage=mean(coverage), pop_size=mean(pop_size)) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    # get population-weighted CM coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied CM coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise(coverage=mean(coverage)) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    input_df_agg_admin = input_df
    
    # take average across simulation seeds
    input_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_agg$scenario = exp_name
  return(input_agg)
}



##############
# SMC
##############
get_smc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  input_df$u5_coverage_total = input_df$coverage_high_access_U5 * input_df$high_access_U5 + input_df$coverage_low_access_U5 * (1-input_df$high_access_U5)
  
  # get mean across rounds within a year
  input_df <- input_df %>% group_by(year, admin_name, seed) %>%
    summarise_all(mean) %>% ungroup()
  
  # create data frame with all admin-years as zeros if SMC isn't specified as given in input file
  admin_years = merge(admin_info, data.frame('year'=seq(min_year, end_year)), all=TRUE)
  input_df = merge(input_df, admin_years, by=c('admin_name','year'), all=TRUE)
  input_df$u5_coverage_total[is.na(input_df$u5_coverage_total)] = 0
  input_df$U5_coverage = input_df$u5_coverage_total
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_U5_cov = input_df$U5_coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  input_ave_df$scenario = exp_name

  return(input_ave_df)
}





##############
# ITN mass campaign
##############
get_itn_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1

  # get mean across rounds within a year
  input_df <- input_df %>% group_by(year, admin_name, seed) %>%
    summarise(itn_u5 = mean(itn_u5)) %>% ungroup()
  
  # create data frame with all admin-years as zeros if SMC isn't specified as given in input file
  admin_years = merge(admin_info, data.frame('year'=seq(min_year, end_year)), all=TRUE)
  input_df = merge(input_df, admin_years, by=c('admin_name','year'), all=TRUE)
  input_df$itn_u5[is.na(input_df$itn_u5)] = 0
  input_df$U5_coverage = input_df$itn_u5
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_U5_cov = input_df$U5_coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  input_ave_df$scenario = exp_name
  
  return(input_ave_df)
}






##############
# EPI vaccine
##############
get_vacc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  input_df = merge(input_df, admin_info, by='admin_name')
  
  # look just at coverage for primary series
  input_df = input_df[input_df$vaccine == 'primary',]
  # true start day of administration is RTSS_day+rtss_touchpoint (due to how things are setup in campaign for birth-triggered delayed vaccine)
  input_df$simday = input_df$RTSS_day + input_df$rtss_touchpoint
  input_df$year = round(input_df$simday/365 + min_year)
  
  # EPI vaccine coverage is specified only once (and then continues for remainder of simulation after the start day)
  # NOTE: assumes same delivery dates across LGAs
  df_repeated = input_df
  for(rr in 1:(end_year - df_repeated$year[1])){
    temp_year = df_repeated
    temp_year$year = df_repeated$year + rr
    temp_year$simday = df_repeated$simday + rr*365
    input_df = rbind(input_df, temp_year)
  }
  
  # # create data frame with all admin-years as zeros before start of vaccine
  # admin_years = merge(admin_info, data.frame('year'=seq(min_year, vacc_start_year-1), coverage=0), all=TRUE)
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_ave_df$scenario = exp_name
  
  return(input_ave_df)
}





##############
# PMC
##############
get_pmc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  
  # true start day of administration is simday+pmc_touchpoints (due to how things are setup in campaign for birth-triggered delayed vaccine)
  input_df$simday = input_df$simday + input_df$pmc_touchpoints
  
  # get average coverage across touchpoints; for simday, take date of earliest touchpoint simday
  input_df <- input_df %>% group_by(admin_name, seed) %>%
    summarise(coverage = mean(coverage), simday=min(simday)) %>% 
    ungroup()  
  input_df$year = round(input_df$simday/365 + min_year)
  input_df = merge(input_df, admin_info, by='admin_name')
  
  # PMC coverage is specified only once (and then continues for remainder of simulation after the start day)
  # NOTE: assumes same delivery dates across LGAs
  df_repeated = input_df
  for(rr in 1:(end_year - df_repeated$year[1])){
    temp_year = df_repeated
    temp_year$year = df_repeated$year + rr
    temp_year$simday = df_repeated$simday + rr*365
    input_df = rbind(input_df, temp_year)
  }
  
  # # create data frame with all admin-years as zeros before start of vaccine
  # admin_years = merge(admin_info, data.frame('year'=seq(min_year, vacc_start_year-1), coverage=0), all=TRUE)
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_ave_df$scenario = exp_name
  
  return(input_ave_df)
}
