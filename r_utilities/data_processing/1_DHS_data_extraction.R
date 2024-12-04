# DHS/MIS output for Sierra Leone
# first, use DHS_code_examination2.R to determine which datasets and code ids to use for each variable

library(foreign)
# if("plyr" %in% (.packages())){
#   if("reshape2" %in% (.packages())){
#     detach("package:reshape2", unload=TRUE) 
#   }
#   detach("package:plyr", unload=TRUE) 
# }
library(dplyr)
# library(rgdal)  # not available for newer version of R
library(raster)
library(sp)
library(pals)
library(prettyGraphs) 
library(stringr)
library(haven)


# function to read in relevant dta file and extract the number of positive and negative results, along with the number tested in each cluster and cluster locations
get_cluster_level_outputs = function(dta_dir, cur_dta, DHS_file_recode_df, var_index, MIS_outputs, include_itn_weight=FALSE, alternate_positive_patterns = c()){
  cur_dta$pos = NA
  cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$pos_pattern[var_index])] = 1
  cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$neg_pattern[var_index])] = 0
  if(length(alternate_positive_patterns)>0){
    for(alt_pos_pattern in alternate_positive_patterns){
      cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == alt_pos_pattern)] = 1
    }
  }
  if(include_itn_weight){
    colnames(cur_dta)[which(colnames(cur_dta) =='hv005')] = 'clust_weight'
    colnames(cur_dta)[which(colnames(cur_dta) =='v005')] = 'clust_weight'
    cur_dta$clust_weight = cur_dta$clust_weight/1000000
    dta_cluster_0 = cur_dta  %>%
      filter(!is.na(cur_dta$pos)) %>%
      group_by_at(DHS_file_recode_df$cluster_id_code[var_index]) %>%
      summarize(rate = mean(pos, na.rm = TRUE),
                num_pos = sum(pos),
                num_tested = n(),
                itn_weights = mean(clust_weight, na.rm=TRUE)) 
  } else{
    dta_cluster_0 = cur_dta  %>%
      filter(!is.na(cur_dta$pos)) %>%
      group_by_at(DHS_file_recode_df$cluster_id_code[var_index]) %>%
      summarize(rate = mean(pos, na.rm = TRUE),
                num_pos = sum(pos),
                num_tested = n()) 
  }

  # match 'hv001' with 'clusterid'
  MIS_outputs = merge(MIS_outputs, dta_cluster_0, by.y=DHS_file_recode_df$cluster_id_code[var_index], by.x='clusterid', all=TRUE)
  
  # add date, if relevant
  if(!('mean_date' %in% colnames(MIS_outputs)) & !is.na(DHS_file_recode_df$month_code[var_index]) & !is.na(DHS_file_recode_df$year_code[var_index])){
    # if month is given as a string instead of numeric, convert
    if(cur_dta[which(colnames(cur_dta) ==DHS_file_recode_df$month_code[var_index])][1,] %in% month.name){
      # cur_dta[which(colnames(cur_dta) == DHS_file_recode_df$month_code[var_index])] = as.character(cur_dta[which(colnames(cur_dta) == DHS_file_recode_df$month_code[var_index])])
      cur_dta$temp_month = sapply(cur_dta[which(colnames(cur_dta) == DHS_file_recode_df$month_code[var_index])], as.character)
      cur_dta[which(colnames(cur_dta) == DHS_file_recode_df$month_code[var_index])] = match(cur_dta$temp_month, month.name)
    }
    cur_dta$date = as.Date(paste0(as.vector(cur_dta[which(colnames(cur_dta) ==DHS_file_recode_df$year_code[var_index])][,1]),'-',
                                  as.vector(cur_dta[which(colnames(cur_dta) ==DHS_file_recode_df$month_code[var_index])][,1]),'-01'),
                           tryFormats = c("%Y-%m-%d","%y-%m-%d"))
    dta_cluster_date = cur_dta  %>%
      filter(!is.na(cur_dta$date)) %>%
      group_by_at(DHS_file_recode_df$cluster_id_code[var_index]) %>%
      summarize(mean_date = mean(date, na.rm = TRUE)) 
    MIS_outputs = merge(MIS_outputs, dta_cluster_date, by.y=DHS_file_recode_df$cluster_id_code[var_index], by.x='clusterid', all=TRUE)
  }
  
  return(MIS_outputs)
}





# function to read in relevant dta file and extract the number of positive and negative results, along with the number tested in each cluster and cluster locations
get_cluster_level_vacc_outputs = function(dta_dir, cur_dta, DHS_file_recode_df, var_index, MIS_outputs, include_itn_weight=FALSE, alternate_positive_patterns = c('vaccination date on card','vaccination marked on card', 'reported by mother'),
                                          survey_month_code = 'v006', survey_year_code='v007', age_month_code='b1', age_year_code='b2', min_age_months_included=12){
  cur_dta$pos = NA
  cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$pos_pattern[var_index])] = 1
  cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$neg_pattern[var_index])] = 0
  if(length(alternate_positive_patterns)>0){
    for(alt_pos_pattern in alternate_positive_patterns){
      cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == alt_pos_pattern)] = 1
    }
  }
  
  # calculate age at time of survey
  cur_dta$survey_date = as.Date(paste0(cur_dta[,which(colnames(cur_dta)==survey_year_code)],'-',cur_dta[,which(colnames(cur_dta)==survey_month_code)],'-28'))
  cur_dta$birth_date = as.Date(paste0(cur_dta[,which(colnames(cur_dta)==age_year_code)],'-',cur_dta[,which(colnames(cur_dta)==age_month_code)],'-01'))
  cur_dta$age = as.numeric(cur_dta$survey_date - cur_dta$birth_date)
  # only include entries above the minimum age
  cur_dta$over_min_age = cur_dta$age > (min_age_months_included*30.4)
  
  # # compare rates of positivity between all ages and ages over min (for debugging/checking)
  # sum(cur_dta$pos, na.rm=T)/sum(!is.na(cur_dta$pos))  # fraction with the vaccine among all ages
  # sum(cur_dta$pos[(cur_dta$over_min_age)], na.rm=T)/sum(!is.na(cur_dta$pos[(cur_dta$over_min_age)]))  # fraction with the vaccine among those over the minimum age

  # remove entries for individuals under the age cutoff (i.e., who are to young to have received the vaccine yet)
  cur_dta$pos[!(cur_dta$over_min_age)] = NA
  
  dta_cluster_0 = cur_dta  %>%
    filter(!is.na(cur_dta$pos)) %>%
    group_by_at(DHS_file_recode_df$cluster_id_code[var_index]) %>%
    summarize(rate = mean(pos, na.rm = TRUE),
              num_pos = sum(pos),
              num_tested = n()) 
  
  # match 'hv001' with 'clusterid'
  MIS_outputs = merge(MIS_outputs, dta_cluster_0, by.y=DHS_file_recode_df$cluster_id_code[var_index], by.x='clusterid', all=TRUE)
  return(MIS_outputs)
}






# function to read in relevant dta file and extract the number of positive and negative results, along with the number tested in each cluster and cluster locations
get_cluster_level_vacc_mics_outputs = function(dta_dir, cur_dta, DHS_file_recode_df, var_index, MIS_outputs, include_itn_weight=FALSE, alternate_positive_patterns = c('vaccination date on card','vaccination marked on card', 'reported by mother'),
                                               min_age_months_included=12){
  cur_dta$pos = NA
  cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$pos_pattern[var_index])] = 1
  cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$neg_pattern[var_index])] = 0
  if(length(alternate_positive_patterns)>0){
    for(alt_pos_pattern in alternate_positive_patterns){
      cur_dta$pos[which(cur_dta[,which(colnames(cur_dta) == DHS_file_recode_df$code[var_index])] == alt_pos_pattern)] = 1
    }
  }

  # only include entries above the minimum age
  cur_dta$over_min_age = cur_dta[[DHS_file_recode_df$age_code[var_index]]] > (min_age_months_included*30.4)
  
  # # compare rates of positivity between all ages and ages over min (for debugging/checking)
  # sum(cur_dta$pos, na.rm=T)/sum(!is.na(cur_dta$pos))  # fraction with the vaccine among all ages
  # sum(cur_dta$pos[(cur_dta$over_min_age)], na.rm=T)/sum(!is.na(cur_dta$pos[(cur_dta$over_min_age)]))  # fraction with the vaccine among those over the minimum age
  
  # remove entries for individuals under the age cutoff (i.e., who are to young to have received the vaccine yet)
  cur_dta$pos[!(cur_dta$over_min_age)] = NA
  
  dta_cluster_0 = cur_dta  %>%
    filter(!is.na(cur_dta$pos)) %>%
    group_by_at(DHS_file_recode_df$cluster_id_code[var_index]) %>%
    summarize(rate = mean(pos, na.rm = TRUE),
              num_pos = sum(pos),
              num_tested = n()) 
  
  # match 'hv001' with 'clusterid'
  MIS_outputs = merge(MIS_outputs, dta_cluster_0, by.y=DHS_file_recode_df$cluster_id_code[var_index], by.x='clusterid', all=TRUE)
  return(MIS_outputs)
}




match_lga_names = function(lga_name){
  lga_name = toupper(lga_name)
  lga_name = str_replace_all(lga_name, pattern=' ', replacement='-')
  lga_name = str_replace_all(lga_name, pattern='/', replacement='-')
  return(lga_name)
}




# do any values in a single df row match one of the positive codes?
any_matches = function(df_row, pos_codes){
  return(any(df_row %in% pos_codes))
}


# create new column describing whether an individual received effective treatment given they received any antimalarial. 
#   Column value will be 1 if an individual received rectal artesunate, IV artesunate, or an ACT and 0 if they received a different antimalarial
#    note: currently does not include individuals who reported country-specific antimalarial ("ml13g", "ml13f"), since it's not clear whether or not those are ACT
received_art_antimalarial = function(dta_dir, DHS_file_recode_df, var_index, art_codes = c("ml13e", "ml13ab"), non_art_codes = c("ml13a", "ml13b", "ml13c", "ml13d", "ml13aa", "ml13da", "ml13h")){
  cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
  # change columns to strings (the any_matches function does not work as expected if there are factors)
  cur_dta = data.frame(lapply(cur_dta, as.character), stringsAsFactors=FALSE)
  # find rows where at least one of the antimalarials matches the positive code
  antimalarial_codes = c(art_codes[which(art_codes %in% colnames(cur_dta))], non_art_codes[which(non_art_codes %in% colnames(cur_dta))])

  if(length(antimalarial_codes)>1){
    cur_dta$any_antimalarial = as.numeric(apply(cur_dta[,antimalarial_codes], 1, any_matches, pos_codes=c('yes', 'Yes', 'YES', 'Y', 1, '1', 'T', 'TRUE', 'True') ))
    cur_dta$any_antimalarial[cur_dta$any_antimalarial != 1] = NA
    # find rows where at least one of the artesunate/ACTs matches the positive code
    art_codes = art_codes[which(art_codes %in% colnames(cur_dta))]
    if(length(art_codes)>1){
      cur_dta$art_antimalarial = as.numeric(apply(cur_dta[,art_codes], 1, any_matches, pos_codes=c('yes', 'Yes', 'YES', 'Y', 1, '1', 'T', 'TRUE', 'True') )) 
      cur_dta$art_given_antimalarial = cur_dta$art_antimalarial / cur_dta$any_antimalarial
    } else if(length(art_codes)==1){
      cur_dta$art_antimalarial = as.numeric(cur_dta[,art_codes] %in% c('yes', 'Yes', 'YES', 'Y', 1, '1', 'T', 'TRUE', 'True'))
      cur_dta$art_given_antimalarial = cur_dta$art_antimalarial / cur_dta$any_antimalarial
    } else{
      cur_dta$art_given_antimalarial = NA
    }
  } else{
    cur_dta$art_given_antimalarial = NA
  }
  return(cur_dta)
}





# create new column describing whether an individual sought any type of treatment (public or private or pharmacy)
sought_any_treatment = function(dta_dir, year, DHS_file_recode_df, var_index, treat_codes=NA){
  if(is.na(treat_codes)){
    if(year !=2013){
      treat_codes = c("h32a", "h32b", "h32c", "h32d", "h32e", "h32f", "h32g", "h32h", "h32i", "h32j", "h32k", "h32l", "h32m", "h32n", "h32o", "h32p", "h32q", "h32r", "h32na", "h32nb", "h32nc", "h32nd", "h32ne", "h32s",  "h32u", "h32v", "h32w", "h32x")  # currently not including traditional practitioner "h32t",
    } else{
      treat_codes = c("h32a", "h32b", "h32c", "h32d", "h32e", "h32f", "h32g", "h32h", "h32i", "h32j", "h32k", "h32l", "h32m", "h32n", "h32o", "h32p", "h32q", "h32r", "h32na", "h32nb", "h32nc", "h32nd", "h32ne", "h32s",  "h32u", "h32t", "h32w", "h32x")  # currently not including traditional practitioner "h32v",
    }
  }
  cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
  # change columns to strings (the any_matches function does not work as expected if there are factors)
  cur_dta = data.frame(lapply(cur_dta, as.character), stringsAsFactors=FALSE)
  # find rows where at least one of the treatment seeking behaviors matches the positive code
  treat_codes = c(treat_codes[which(treat_codes %in% colnames(cur_dta))])
  
  if(length(treat_codes)>1){
    cur_dta$sought_treatment = as.numeric(apply(cur_dta[,treat_codes], 1, any_matches, pos_codes=c('yes', 'Yes', 'YES', 'Y', 1, '1', 'T', 'TRUE', 'True') ))
    cur_dta$responded_sought_treatment = as.numeric(apply(cur_dta[,treat_codes], 1, any_matches, pos_codes=c('yes', 'Yes', 'YES', 'Y', 1, '1', 'T', 'TRUE', 'True', 'no', 'No', 'NO', 'N', 0, '0', 'F', 'FALSE', 'False') ))
    cur_dta$sought_treatment[cur_dta$responded_sought_treatment != 1] = NA
  } else{
    cur_dta$sought_treatment = NA
  }
  return(cur_dta)
}


#################################################################################################################
# main function to extract DHS data and plot maps of the results
#################################################################################################################
extract_DHS_data = function(hbhi_dir, dta_dir, years, admin_shape, ds_pop_df_filename, min_num_total=30, variables=c('mic', 'rdt','itn_all','itn_u5','itn_5_10','itn_10_15','itn_15_20','itn_o20','iptp','cm','blood_test', 'art_given_antimal')){
  
  ####=========================================================================================================####
  # iterate through years, creating csvs with cluster-level and admin-level counts and rates for all variables
  ####=========================================================================================================####
  
  for(yy in 1:length(years)){
    year = years[yy]
    DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',year,'_files_recodes_for_sims.csv'))
    location_index = which(DHS_file_recode_df$variable == 'locations')
    locations_shp = shapefile(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[location_index], '/', DHS_file_recode_df$filename[location_index]))
    locations = data.frame(clusterid = locations_shp$DHSCLUST, latitude=locations_shp$LATNUM, longitude=locations_shp$LONGNUM)
    MIS_outputs = locations
    

    ### - - - - - - - - - - - - - - - - - - ###
    # PfPR (microscopy)
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'mic')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'mic_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'mic_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'mic_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # PfPR (RDT)
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'rdt')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'rdt_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'rdt_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'rdt_num_total'
    }
    
    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - all ages
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs, include_itn_weight=TRUE)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_all_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_all_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_all_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - U5
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      cur_dta = cur_dta[cur_dta[[DHS_file_recode_df$age_code[var_index]]]<=5,]
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_u5_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_u5_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_u5_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - O5
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      cur_dta = cur_dta[cur_dta[[DHS_file_recode_df$age_code[var_index]]]>5,]
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_o5_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_o5_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_o5_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - 5-10
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      cur_dta = cur_dta[intersect(which(cur_dta[[DHS_file_recode_df$age_code[var_index]]]>5), which(cur_dta[[DHS_file_recode_df$age_code[var_index]]]<=10)),]
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_5_10_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_5_10_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_5_10_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - 10-15
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      cur_dta = cur_dta[intersect(which(cur_dta[[DHS_file_recode_df$age_code[var_index]]]>10), which(cur_dta[[DHS_file_recode_df$age_code[var_index]]]<=15)),]
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_10_15_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_10_15_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_10_15_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - 15-20
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      cur_dta = cur_dta[intersect(which(cur_dta[[DHS_file_recode_df$age_code[var_index]]]>15), which(cur_dta[[DHS_file_recode_df$age_code[var_index]]]<=20)),]
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_15_20_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_15_20_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_15_20_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ITNs - >20
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'itn')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      cur_dta = cur_dta[cur_dta[[DHS_file_recode_df$age_code[var_index]]]>20,]
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'itn_o20_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'itn_o20_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'itn_o20_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # Case management - sought treatment at a facility
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'cm')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'cm_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'cm_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'cm_num_total'
    }
    
    ### - - - - - - - - - - - - - - - - - - ###
    # Case management - receive treatment (1-no treatment or advice sought)
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'received_treatment')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'received_treatment_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'received_treatment_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'received_treatment_num_total'
    }

    
    
    ### - - - - - - - - - - - - - - - - - - ###
    # Case management - sought advice or treatment (sum across all types of treatment except traditional healers)
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'sought_treatment')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = sought_any_treatment(dta_dir, year, DHS_file_recode_df, var_index, treat_codes=NA)
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'sought_treatment_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'sought_treatment_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'sought_treatment_num_total'
    }
    
    
    ### - - - - - - - - - - - - - - - - - - ###
    # Case management - heel prick or blood test
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'blood_test')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'blood_test_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'blood_test_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'blood_test_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # IPTp
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'iptp')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'iptp_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'iptp_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'iptp_num_total'
    }

    ### - - - - - - - - - - - - - - - - - - ###
    # ACT or artesunate given antimalarial
    ### - - - - - - - - - - - - - - - - - - ###
    var_index = which(DHS_file_recode_df$variable == 'art_given_antimal')
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = received_art_antimalarial(dta_dir=dta_dir, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index)
      MIS_outputs=get_cluster_level_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs)
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = 'art_rate'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = 'art_num_true'
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = 'art_num_total'
    }

    

    

    

    ####=========================================================================================================####
    # determine which clusters are in which admins
    ####=========================================================================================================####
    # turn MIS output data frame into spatial points data frame
    points_crs = crs(admin_shape)
    MIS_shape = SpatialPointsDataFrame(MIS_outputs[,c('longitude', 'latitude')],
                                       MIS_outputs,
                                       proj4string = points_crs)
    # find which admins each cluster belongs to
    MIS_locations = over(MIS_shape, admin_shape)
    if(nrow(MIS_locations) == nrow(MIS_shape)){
      MIS_shape$NOMDEP = MIS_locations$NOMDEP
      MIS_shape$NOMREGION = MIS_locations$NOMREGION
      # MIS_shape$NAME_1 = MIS_locations$NAME_1
    }

    write.csv(as.data.frame(MIS_shape), paste0(hbhi_dir, '/estimates_from_DHS/DHS_cluster_outputs_', years[yy], '.csv'))

    
    ####=========================================================================================================####
    #   get values for each admin
    ####=========================================================================================================####
    # aggregate across clusters to get total number tested and positive within each admin and within each region for all variables
    MIS_shape = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_cluster_outputs_', years[yy], '.csv'))[,-1]
    # remove rows without known admin/region
    MIS_shape = MIS_shape[!is.na(MIS_shape$NOMREGION),]

    # admin level values
    include_cols = c(which(names(MIS_shape) %in% c('NOMREGION','NOMDEP')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
    admin_sums = MIS_shape[,include_cols] %>% 
      group_by( NOMREGION, NOMDEP) %>%
      summarise_all(sum, na.rm = TRUE)
    
    for(var in variables){
      if(paste0(var,'_num_true') %in% colnames(admin_sums)){
        admin_sums[[paste0(var, '_rate')]] = admin_sums[[paste0(var,'_num_true')]] / admin_sums[[paste0(var,'_num_total')]]
      }
    }
    
    # for case management, there are several ways we could attempt to estimate the treatment-seeking and effective treatment rates. 
    # Some relevant components are 
    #    a) the treatment-seeking rate, but this does not include whether the individual received effective/appropriate treatment
    #    b) the blood-test rate, but this may ignore individuals who were given ACTs without receiving a test
    #    c) the fraction of individual who received various antimalarials for fever or difficulty breathing in past two weeks, but the denominator will include individuals without malaria who SHOULDN'T have received antimalarials
    # A previous approach was to assume that the rate of receiving effective treatment is halfway in between a) and b) - i.e., half of the people who seek treatment but don't receive a blood test are nonetheless given ACTs
    # The current approach is to use the national fraction of individuals who are given ACTs/artesunate among those who receive any antimalarial, multiplied by the local treatment-seeking rate
    use_art_probs = TRUE
    use_cm_probs = TRUE  # set which of the variables for treatment-seeking is used: any facility (TRUE) or any treatment/advice aside from traditional (FALSE)
    if(use_art_probs){  # new approach using probability someone who gets an antimalarial gets an ACT
      # use the country-wide, cluster-weighted probability of ACT/artesunate use given some type of antimalarial use if cluster weights are available
      if('art_rate' %in% colnames(MIS_shape) & 'itn_weights' %in% colnames(MIS_shape)){
        national_act_rate = stats::weighted.mean(MIS_shape$art_rate, MIS_shape$itn_weights, na.rm=TRUE)
      } else if('art_num_true' %in% colnames(MIS_shape) & 'art_num_total' %in% colnames(MIS_shape) ){
        national_act_rate = sum(MIS_shape$art_num_true, na.rm=TRUE) / sum(MIS_shape$art_num_total, na.rm=TRUE)
      } else{
        warning('ACT rates not available, using mean between treatment-seeking and blood-test rates instead.')
        national_act_rate = NA
        use_art_probs = FALSE
      }
      admin_sums$national_act_rate = national_act_rate
      if(use_cm_probs){
        admin_sums$effective_treatment_rate = admin_sums$national_act_rate * admin_sums$cm_rate
        admin_sums$effective_treatment_num_total = admin_sums$cm_num_total
      } else{
        admin_sums$effective_treatment_rate = admin_sums$national_act_rate * admin_sums$sought_treatment_rate
        admin_sums$effective_treatment_num_total = admin_sums$sought_treatment_num_total
      }
      admin_sums$effective_treatment_num_true = admin_sums$effective_treatment_num_total * admin_sums$effective_treatment_rate
    }
    if(!use_art_probs){  # old approach using an average of treatment-seeking and blood-test rates
      if('cm_num_total' %in% colnames(admin_sums) & 'blood_test_num_total' %in% colnames(admin_sums)){
        admin_sums$effective_treatment_num_total = (admin_sums$cm_num_total + admin_sums$blood_test_num_total)/2
        admin_sums$effective_treatment_num_true = (admin_sums$cm_num_true + admin_sums$blood_test_num_true)/2
        admin_sums$effective_treatment_rate = (admin_sums$cm_rate + admin_sums$blood_test_rate)/2
      }
    }
   
    
    # for date, take average dates across all clusters in the admin
    if('mean_date' %in% colnames(MIS_shape)){
      MIS_shape$mean_date = as.Date(MIS_shape$mean_date)
      admin_dates = MIS_shape[,c('NOMREGION','NOMDEP', 'mean_date')] %>% 
        group_by( NOMREGION, NOMDEP) %>%
        summarise_all(mean, na.rm = TRUE)
      admin_sums = merge(admin_sums, admin_dates, all.x=TRUE)
    }

    # add any admins that did not have any DHS clusters (with all NAs and 0s) and record which state and archetype each cluster belongs to
    archetype_info = read.csv(ds_pop_df_filename)
    colnames(archetype_info)[colnames(archetype_info)=='LGA'] = 'NOMDEP'
    colnames(archetype_info)[colnames(archetype_info)=='DS'] = 'NOMDEP'
    colnames(archetype_info)[colnames(archetype_info)=='admin_name'] = 'NOMDEP'
    colnames(archetype_info)[colnames(archetype_info)=='State'] = 'NOMREGION'
    # colnames(archetype_info)[colnames(archetype_info)=='NOMDEP'] = 'NOMDEP_target'
    archetype_info$name_match = sapply(archetype_info$NOMDEP, match_lga_names)
    archetype_info = archetype_info[,c('name_match', 'NOMDEP', 'NOMREGION', 'Archetype')]
    admin_sums$name_match = sapply(admin_sums$NOMDEP, match_lga_names)
    colnames(admin_sums)[colnames(admin_sums) == 'NOMDEP'] = 'NOMDEP_dhs_orig'
    admin_sums_expanded = merge(admin_sums, archetype_info, all=TRUE)
    # admin_sums_expanded$NOMDEP[is.na(admin_sums_expanded$NOMDEP)] = admin_sums_expanded$NOMDEP_backup[is.na(admin_sums_expanded$NOMDEP)]
    # check that all names have been matched successfully
    if(length(admin_sums$name_match[which(!(admin_sums$name_match %in% archetype_info$name_match))])>0) warning('Not all LGA names from the shapefile were matched with archetype file')
    if(length(admin_sums$NOMREGION[which(!(admin_sums$NOMREGION %in% archetype_info$NOMREGION))])>0) warning('Not all state names from the shapefile were matched with archetype file')
    # remove extra columns
    admin_sums_expanded = admin_sums_expanded[,-c(which(colnames(admin_sums_expanded) %in% c('NOMDEP_dhs_orig', 'name_match')))]
    # change to zero sample size for admins that were not included in DHS
    admin_sums_expanded[,grep('num_total', names(admin_sums_expanded))][is.na(admin_sums_expanded[,grep('num_total', names(admin_sums_expanded))])] = 0
    write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_', years[yy], '.csv'))
    admin_sums = admin_sums_expanded
    
    #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # when the number surveyed in a admin is lower than the threshold, use the region value instead
    #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # region level values
    include_cols = c(which(names(MIS_shape) %in% c('NOMREGION')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
    region_sums = MIS_shape[,include_cols] %>% 
      group_by(NOMREGION) %>%
      summarise_all(sum, na.rm = TRUE)
    for(var in variables){
      if(paste0(var,'_num_true') %in% colnames(region_sums)){
        region_sums[[paste0(var, '_rate')]] = region_sums[[paste0(var,'_num_true')]] / region_sums[[paste0(var,'_num_total')]]
      }
    }
    
    # for date, take average dates across all clusters in the admin
    if('mean_date' %in% colnames(MIS_shape)){
      MIS_shape$mean_date = as.Date(MIS_shape$mean_date)
      region_dates = MIS_shape[,c('NOMREGION', 'mean_date')] %>% 
        group_by( NOMREGION) %>%
        summarise_all(mean, na.rm = TRUE)
      region_sums = merge(region_sums, region_dates, all.x=TRUE)
    }
    
    for(var in variables){
      if(paste0(var,'_num_true') %in% colnames(admin_sums)){
        for(i_admin in 1:nrow(admin_sums)){
          if(admin_sums[[paste0(var,'_num_total')]][i_admin]<min_num_total){
            admin_sums[[paste0(var,'_num_total')]][i_admin] = region_sums[[paste0(var,'_num_total')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
            admin_sums[[paste0(var,'_num_true')]][i_admin] = region_sums[[paste0(var,'_num_true')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
            admin_sums[[paste0(var,'_rate')]][i_admin] = region_sums[[paste0(var,'_rate')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
            if('mean_date' %in% colnames(admin_sums) & grepl('itn', var)){
              admin_sums$mean_date[i_admin] = region_sums$mean_date[region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
            }
          }
        }
        # par(mfrow=c(1,2))
        # hist(admin_sums[[paste0(var,'_rate')]], main=paste0(var, ' - raw values'))
        # hist(admin_sums[[paste0(var,'_rate2')]], main=paste0(var,' - using region when small sample size'))
      }
    }

    # calculate estimate for effective treatment-seeking rate
    if(use_art_probs){  # new approach using probability someone who gets an antimalarial gets an ACT
      admin_sums$national_act_rate = national_act_rate
      if(use_cm_probs){
        admin_sums$effective_treatment_rate = admin_sums$national_act_rate * admin_sums$cm_rate
        admin_sums$effective_treatment_num_total = admin_sums$cm_num_total
      } else{
        admin_sums$effective_treatment_rate = admin_sums$national_act_rate * admin_sums$sought_treatment_rate
        admin_sums$effective_treatment_num_total = admin_sums$sought_treatment_num_total
      }
      admin_sums$effective_treatment_num_true = admin_sums$effective_treatment_num_total * admin_sums$effective_treatment_rate
    } else{  # old approach using an average of treatment-seeking and blood-test rates
      if('cm_num_total' %in% colnames(admin_sums) & 'blood_test_num_total' %in% colnames(admin_sums)){
        admin_sums$effective_treatment_num_total = (admin_sums$cm_num_total + admin_sums$blood_test_num_total)/2
        admin_sums$effective_treatment_num_true = (admin_sums$cm_num_true + admin_sums$blood_test_num_true)/2
        admin_sums$effective_treatment_rate = (admin_sums$cm_rate + admin_sums$blood_test_rate)/2
      }
    }
    

    #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # when the number surveyed in a admin1 is lower than the threshold, use the archetype value instead
    #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # archetype level values
    MIS_shape = merge(MIS_shape, distinct(archetype_info[,c('NOMREGION', 'Archetype')]), all=TRUE)
    include_cols = c(which(names(MIS_shape) %in% c('Archetype')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
    arch_sums = MIS_shape[,include_cols] %>% 
      group_by(Archetype) %>%
      summarise_all(sum, na.rm = TRUE)
    for(var in variables){
      if(paste0(var,'_num_true') %in% colnames(arch_sums)){
        arch_sums[[paste0(var, '_rate')]] = arch_sums[[paste0(var,'_num_true')]] / arch_sums[[paste0(var,'_num_total')]]
      }
    }
    
    # for date, take average dates across all clusters in the admin
    if('mean_date' %in% colnames(MIS_shape)){
      MIS_shape$mean_date = as.Date(MIS_shape$mean_date)
      arch_dates = MIS_shape[,c('Archetype', 'mean_date')] %>% 
        group_by( Archetype) %>%
        summarise_all(mean, na.rm = TRUE)
      arch_sums = merge(arch_sums, arch_dates, all.x=TRUE)
    }
    
    for(var in variables){
      if(paste0(var,'_num_true') %in% colnames(admin_sums)){
        for(i_admin in 1:nrow(admin_sums)){
          if(admin_sums[[paste0(var,'_num_total')]][i_admin]<min_num_total){
            admin_sums[[paste0(var,'_num_total')]][i_admin] = arch_sums[[paste0(var,'_num_total')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
            admin_sums[[paste0(var,'_num_true')]][i_admin] = arch_sums[[paste0(var,'_num_true')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
            admin_sums[[paste0(var,'_rate')]][i_admin] = arch_sums[[paste0(var,'_rate')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
            if('mean_date' %in% colnames(admin_sums) & grepl('itn', var)){
              admin_sums$mean_date[i_admin] = arch_sums$mean_date[arch_sums$Archetype == admin_sums$Archetype[i_admin]]
            }
          }
        }
        # par(mfrow=c(1,2))
        # hist(admin_sums[[paste0(var,'_rate')]], main=paste0(var, ' - raw values'))
        # hist(admin_sums[[paste0(var,'_rate2')]], main=paste0(var,' - using region when small sample size'))
      }
    }
    
    
    # calculate estimate for effective treatment-seeking rate
    if(use_art_probs){  # new approach using probability someone who gets an antimalarial gets an ACT
      admin_sums$national_act_rate = national_act_rate
      if(use_cm_probs){
        admin_sums$effective_treatment_rate = admin_sums$national_act_rate * admin_sums$cm_rate
        admin_sums$effective_treatment_num_total = admin_sums$cm_num_total
      } else{
        admin_sums$effective_treatment_rate = admin_sums$national_act_rate * admin_sums$sought_treatment_rate
        admin_sums$effective_treatment_num_total = admin_sums$sought_treatment_num_total
      }
      admin_sums$effective_treatment_num_true = admin_sums$effective_treatment_num_total * admin_sums$effective_treatment_rate
    } else{  # old approach using an average of treatment-seeking and blood-test rates
      if('cm_num_total' %in% colnames(admin_sums) & 'blood_test_num_total' %in% colnames(admin_sums)){
        admin_sums$effective_treatment_num_total = (admin_sums$cm_num_total + admin_sums$blood_test_num_total)/2
        admin_sums$effective_treatment_num_true = (admin_sums$cm_num_true + admin_sums$blood_test_num_true)/2
        admin_sums$effective_treatment_rate = (admin_sums$cm_rate + admin_sums$blood_test_rate)/2
      }
    }


    # save different version with specified minimum sample size in admin before using region value
    write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_minN',min_num_total,'_', years[yy], '.csv'))
    # write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_minN',min_num_total,'_includeArch_', years[yy], '.csv'))
    
  }
}




# extract cluster-level data for EPI vaccination coverages for single year
extract_vaccine_DHS_data = function(hbhi_dir, dta_dir, year, admin_shape, ds_pop_df_filename, min_num_total=30, vaccine_variables=c('vacc_dpt1', 'vacc_dpt2', 'vacc_dpt3'), vaccine_alternate_positive_patterns=c('reported by mother', 'vaccination marked on card'),
                                    survey_month_code = 'v006', survey_year_code='v007', age_month_code='b1', age_year_code='b2', min_age_months_included=12){
  
  ####=========================================================================================================####
  # create csvs with cluster-level and admin-level counts and rates for all vaccination variables
  ####=========================================================================================================####
  
  DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',year,'_files_recodes_for_sims.csv'))
  location_index = which(DHS_file_recode_df$variable == 'locations')
  locations_shp = shapefile(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[location_index], '/', DHS_file_recode_df$filename[location_index]))
  locations = data.frame(clusterid = locations_shp$DHSCLUST, latitude=locations_shp$LATNUM, longitude=locations_shp$LONGNUM)
  MIS_outputs = locations
  
  for(vv in 1:length(vaccine_variables)){
    var_index = which(DHS_file_recode_df$variable == vaccine_variables[vv])
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_vacc_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs, alternate_positive_patterns=vaccine_alternate_positive_patterns,
                                                 survey_month_code=survey_month_code, survey_year_code=survey_year_code, age_month_code=age_month_code, age_year_code=age_year_code, min_age_months_included=DHS_file_recode_df$min_age_months_to_include[var_index])
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = paste0(vaccine_variables[vv], '_rate')
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = paste0(vaccine_variables[vv], '_num_true')
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = paste0(vaccine_variables[vv], '_num_total')
    }
    
  }
  
  ####=========================================================================================================####
  # determine which clusters are in which admins
  ####=========================================================================================================####
  # turn MIS output data frame into spatial points data frame
  points_crs = crs(admin_shape)
  MIS_shape = SpatialPointsDataFrame(MIS_outputs[,c('longitude', 'latitude')],
                                     MIS_outputs,
                                     proj4string = points_crs)
  # find which admins each cluster belongs to
  MIS_locations = over(MIS_shape, admin_shape)
  if(nrow(MIS_locations) == nrow(MIS_shape)){
    MIS_shape$NOMDEP = MIS_locations$NOMDEP
    MIS_shape$NOMREGION = MIS_locations$NOMREGION
    # MIS_shape$NAME_1 = MIS_locations$NAME_1
  }
  write.csv(as.data.frame(MIS_shape), paste0(hbhi_dir, '/estimates_from_DHS/DHS_vaccine_cluster_outputs_', year, '.csv'))
  
  
  ####=========================================================================================================####
  #   get values for each admin
  ####=========================================================================================================####
  # aggregate across clusters to get total number tested and positive within each admin and within each region for all variables
  MIS_shape = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_vaccine_cluster_outputs_', year, '.csv'))[,-1]
  # remove rows without known admin/region
  MIS_shape = MIS_shape[!is.na(MIS_shape$NOMREGION),]
  
  # admin level values
  include_cols = c(which(names(MIS_shape) %in% c('NOMREGION','NOMDEP')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
  admin_sums = MIS_shape[,include_cols] %>% 
    group_by(NOMREGION, NOMDEP) %>%
    summarise_all(sum, na.rm = TRUE)
  
  for(var in vaccine_variables){
    if(paste0(var,'_num_true') %in% colnames(admin_sums)){
      admin_sums[[paste0(var, '_rate')]] = admin_sums[[paste0(var,'_num_true')]] / admin_sums[[paste0(var,'_num_total')]]
    }
  }
  
  # add any admins that did not have any DHS clusters (with all NAs and 0s) and record which state and archetype each cluster belongs to
  archetype_info = read.csv(ds_pop_df_filename)
  colnames(archetype_info)[colnames(archetype_info)=='LGA'] = 'NOMDEP'
  colnames(archetype_info)[colnames(archetype_info)=='DS'] = 'NOMDEP'
  colnames(archetype_info)[colnames(archetype_info)=='admin_name'] = 'NOMDEP'
  colnames(archetype_info)[colnames(archetype_info)=='State'] = 'NOMREGION'
  # colnames(archetype_info)[colnames(archetype_info)=='NOMDEP'] = 'NOMDEP_target'
  archetype_info$name_match = sapply(archetype_info$NOMDEP, match_lga_names)
  archetype_info = archetype_info[,c('name_match', 'NOMDEP', 'NOMREGION', 'Archetype')]
  admin_sums$name_match = sapply(admin_sums$NOMDEP, match_lga_names)
  colnames(admin_sums)[colnames(admin_sums) == 'NOMDEP'] = 'NOMDEP_dhs_orig'
  admin_sums_expanded = merge(admin_sums, archetype_info, all=TRUE)
  # admin_sums_expanded$NOMDEP[is.na(admin_sums_expanded$NOMDEP)] = admin_sums_expanded$NOMDEP_backup[is.na(admin_sums_expanded$NOMDEP)]
  # check that all names have been matched successfully
  if(length(admin_sums$name_match[which(!(admin_sums$name_match %in% archetype_info$name_match))])>0) warning('Not all LGA names from the shapefile were matched with archetype file')
  if(length(admin_sums$NOMREGION[which(!(admin_sums$NOMREGION %in% archetype_info$NOMREGION))])>0) warning('Not all state names from the shapefile were matched with archetype file')
  # remove extra columns
  admin_sums_expanded = admin_sums_expanded[,-c(which(colnames(admin_sums_expanded) %in% c('NOMDEP_dhs_orig', 'name_match')))]
  # change to zero sample size for admins that were not included in DHS
  admin_sums_expanded[,grep('num_total', names(admin_sums_expanded))][is.na(admin_sums_expanded[,grep('num_total', names(admin_sums_expanded))])] = 0
  write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/DHS_vaccine_admin_', year, '.csv'))
  admin_sums = admin_sums_expanded
  
  
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # when the number surveyed in a admin is lower than the threshold, use the region value instead
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # region level values
  include_cols = c(which(names(MIS_shape) %in% c('NOMREGION')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
  region_sums = MIS_shape[,include_cols] %>% 
    group_by(NOMREGION) %>%
    summarise_all(sum, na.rm = TRUE)
  for(var in vaccine_variables){
    if(paste0(var,'_num_true') %in% colnames(region_sums)){
      region_sums[[paste0(var, '_rate')]] = region_sums[[paste0(var,'_num_true')]] / region_sums[[paste0(var,'_num_total')]]
    }
  }
  
  for(var in vaccine_variables){
    if(paste0(var,'_num_true') %in% colnames(admin_sums)){
      for(i_admin in 1:nrow(admin_sums)){
        if(admin_sums[[paste0(var,'_num_total')]][i_admin]<min_num_total){
          admin_sums[[paste0(var,'_num_total')]][i_admin] = region_sums[[paste0(var,'_num_total')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
          admin_sums[[paste0(var,'_num_true')]][i_admin] = region_sums[[paste0(var,'_num_true')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
          admin_sums[[paste0(var,'_rate')]][i_admin] = region_sums[[paste0(var,'_rate')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
        }
      }
    }
  }
  

  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # when the number surveyed in a admin1 is lower than the threshold, use the archetype value instead
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # archetype level values
  MIS_shape = merge(MIS_shape, distinct(archetype_info[,c('NOMREGION', 'Archetype')]), all=TRUE)
  include_cols = c(which(names(MIS_shape) %in% c('Archetype')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
  arch_sums = MIS_shape[,include_cols] %>% 
    group_by(Archetype) %>%
    summarise_all(sum, na.rm = TRUE)
  for(var in vaccine_variables){
    if(paste0(var,'_num_true') %in% colnames(arch_sums)){
      arch_sums[[paste0(var, '_rate')]] = arch_sums[[paste0(var,'_num_true')]] / arch_sums[[paste0(var,'_num_total')]]
    }
  }
  
  for(var in vaccine_variables){
    if(paste0(var,'_num_true') %in% colnames(admin_sums)){
      for(i_admin in 1:nrow(admin_sums)){
        if(admin_sums[[paste0(var,'_num_total')]][i_admin]<min_num_total){
          admin_sums[[paste0(var,'_num_total')]][i_admin] = arch_sums[[paste0(var,'_num_total')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
          admin_sums[[paste0(var,'_num_true')]][i_admin] = arch_sums[[paste0(var,'_num_true')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
          admin_sums[[paste0(var,'_rate')]][i_admin] = arch_sums[[paste0(var,'_rate')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
        }
      }
    }
  }
  
  # save different version with specified minimum sample size in admin before using region value
  write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/DHS_vaccine_admin_minN',min_num_total,'_', year, '.csv'))
}


 













# extract cluster-level data for EPI vaccination coverages for single year
extract_vaccine_MICS_data = function(hbhi_dir, dta_dir, year, admin_shape, ds_pop_df_filename, min_num_total=30, vaccine_variables=c('vacc_dpt1', 'vacc_dpt2', 'vacc_dpt3'), vaccine_alternate_positive_patterns=c('reported by mother', 'vaccination marked on card'),
                                    min_age_months_included=12, use_admin_level=FALSE){
  
  ####=========================================================================================================####
  # create csvs with cluster-level and admin-level counts and rates for all vaccination variables
  ####=========================================================================================================####
  
  DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/MICS_',year,'_files_recodes_for_sims.csv'))
  location_index = which(DHS_file_recode_df$variable == 'locations')
  locations_shp = shapefile(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[location_index], '/', DHS_file_recode_df$filename[location_index]))
  locations = data.frame(clusterid = locations_shp$HH1, latitude=locations_shp$LATITUDE, longitude=locations_shp$LONGITUDE)
  MIS_outputs = locations
  
  for(vv in 1:length(vaccine_variables)){
    var_index = which(DHS_file_recode_df$variable == vaccine_variables[vv])
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read_sav(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      MIS_outputs=get_cluster_level_vacc_mics_outputs(dta_dir=dta_dir, cur_dta=cur_dta, DHS_file_recode_df=DHS_file_recode_df, var_index=var_index, MIS_outputs=MIS_outputs, alternate_positive_patterns=vaccine_alternate_positive_patterns,
                                                      min_age_months_included=DHS_file_recode_df$min_age_months_to_include[var_index])
      colnames(MIS_outputs)[colnames(MIS_outputs)=='rate'] = paste0(vaccine_variables[vv], '_rate')
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_pos'] = paste0(vaccine_variables[vv], '_num_true')
      colnames(MIS_outputs)[colnames(MIS_outputs)=='num_tested'] = paste0(vaccine_variables[vv], '_num_total')
    }
  }
  
  ####=========================================================================================================####
  # determine which clusters are in which admins
  ####=========================================================================================================####
  # turn MIS output data frame into spatial points data frame
  points_crs = crs(admin_shape)
  MIS_shape = SpatialPointsDataFrame(MIS_outputs[,c('longitude', 'latitude')],
                                     MIS_outputs,
                                     proj4string = points_crs)
  # find which admins each cluster belongs to
  MIS_locations = over(MIS_shape, admin_shape)
  MIS_locations$NOMDEP = MIS_locations$GEONAMET
  MIS_locations$NOMREGION = MIS_locations$GEONAMES
  if(nrow(MIS_locations) == nrow(MIS_shape)){
    MIS_shape$NOMDEP = MIS_locations$NOMDEP
    MIS_shape$NOMREGION = MIS_locations$NOMREGION
    # MIS_shape$NAME_1 = MIS_locations$NAME_1
  }
   write.csv(as.data.frame(MIS_shape), paste0(hbhi_dir, '/estimates_from_DHS/MICS_vaccine_cluster_outputs_', year, '.csv'))
  
  
  ####=========================================================================================================####
  #   get values for each admin
  ####=========================================================================================================####
  # aggregate across clusters to get total number tested and positive within each admin and within each region for all variables
  MIS_shape = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/MICS_vaccine_cluster_outputs_', year, '.csv'))[,-1]
  # remove rows without known admin/region
  MIS_shape = MIS_shape[!is.na(MIS_shape$NOMREGION),]
  
  # standardize admin and state names
  archetype_info = read.csv(ds_pop_df_filename)
  if(use_admin_level){
    MIS_shape$admin_name = standardize_admin_names_in_vector(target_names=archetype_info$admin_name, origin_names=MIS_shape$NOMDEP)
    MIS_shape = standardize_admin_names_in_df(target_names_df=archetype_info, origin_names_df=MIS_shape, target_names_col='admin_name', origin_names_col='NOMDEP', additional_id_col='State', possible_suffixes=c(1,2, 3))
  }
  MIS_shape$State = MIS_shape$NOMREGION
  MIS_shape = standardize_state_names_in_df(target_names_df=archetype_info, origin_names_df=MIS_shape, target_names_col='State', origin_names_col='State')
  
  if(use_admin_level){
    # admin level values
    include_cols = c(which(names(MIS_shape) %in% c('NOMREGION','NOMDEP')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
    admin_sums = MIS_shape[,include_cols] %>% 
      group_by(NOMREGION, NOMDEP) %>%
      summarise_all(sum, na.rm = TRUE)
    
    for(var in vaccine_variables){
      if(paste0(var,'_num_true') %in% colnames(admin_sums)){
        admin_sums[[paste0(var, '_rate')]] = admin_sums[[paste0(var,'_num_true')]] / admin_sums[[paste0(var,'_num_total')]]
      }
    }
    
    # add any admins that did not have any DHS clusters (with all NAs and 0s) and record which state and archetype each cluster belongs to
    colnames(archetype_info)[colnames(archetype_info)=='LGA'] = 'NOMDEP'
    colnames(archetype_info)[colnames(archetype_info)=='DS'] = 'NOMDEP'
    colnames(archetype_info)[colnames(archetype_info)=='admin_name'] = 'NOMDEP'
    colnames(archetype_info)[colnames(archetype_info)=='State'] = 'NOMREGION'
    # colnames(archetype_info)[colnames(archetype_info)=='NOMDEP'] = 'NOMDEP_target'
    archetype_info$name_match = sapply(archetype_info$NOMDEP, match_lga_names)
    archetype_info = archetype_info[,c('name_match', 'NOMDEP', 'NOMREGION', 'Archetype')]
    admin_sums$name_match = sapply(admin_sums$NOMDEP, match_lga_names)
    colnames(admin_sums)[colnames(admin_sums) == 'NOMDEP'] = 'NOMDEP_dhs_orig'
    admin_sums_expanded = merge(admin_sums, archetype_info, all=TRUE)
    # admin_sums_expanded$NOMDEP[is.na(admin_sums_expanded$NOMDEP)] = admin_sums_expanded$NOMDEP_backup[is.na(admin_sums_expanded$NOMDEP)]
    # check that all names have been matched successfully
    if(length(admin_sums$name_match[which(!(admin_sums$name_match %in% archetype_info$name_match))])>0) warning('Not all LGA names from the shapefile were matched with archetype file')
    if(length(admin_sums$NOMREGION[which(!(admin_sums$NOMREGION %in% archetype_info$NOMREGION))])>0) warning('Not all state names from the shapefile were matched with archetype file')
    # remove extra columns
    admin_sums_expanded = admin_sums_expanded[,-c(which(colnames(admin_sums_expanded) %in% c('NOMDEP_dhs_orig', 'name_match')))]
    # change to zero sample size for admins that were not included in DHS
    admin_sums_expanded[,grep('num_total', names(admin_sums_expanded))][is.na(admin_sums_expanded[,grep('num_total', names(admin_sums_expanded))])] = 0
    write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/MICS_vaccine_admin_', year, '.csv'))
    admin_sums = admin_sums_expanded
  }
  
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # when the number surveyed in a admin is lower than the threshold, use the region value instead
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # region level values
  include_cols = c(which(names(MIS_shape) %in% c('NOMREGION', 'State')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
  region_sums = MIS_shape[,include_cols] %>% 
    group_by(NOMREGION, State) %>%
    summarise_all(sum, na.rm = TRUE)
  for(var in vaccine_variables){
    if(paste0(var,'_num_true') %in% colnames(region_sums)){
      region_sums[[paste0(var, '_rate')]] = region_sums[[paste0(var,'_num_true')]] / region_sums[[paste0(var,'_num_total')]]
    }
  }
  write.csv(as.data.frame(region_sums), paste0(hbhi_dir, '/estimates_from_DHS/MICS_vaccine_state_', year, '.csv'))
  
  
  if(use_admin_level){
    for(var in vaccine_variables){
      if(paste0(var,'_num_true') %in% colnames(admin_sums)){
        for(i_admin in 1:nrow(admin_sums)){
          if(admin_sums[[paste0(var,'_num_total')]][i_admin]<min_num_total){
            admin_sums[[paste0(var,'_num_total')]][i_admin] = region_sums[[paste0(var,'_num_total')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
            admin_sums[[paste0(var,'_num_true')]][i_admin] = region_sums[[paste0(var,'_num_true')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
            admin_sums[[paste0(var,'_rate')]][i_admin] = region_sums[[paste0(var,'_rate')]][region_sums$NOMREGION == admin_sums$NOMREGION[i_admin]]
          }
        }
      }
    }

  
    #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # when the number surveyed in a admin1 is lower than the threshold, use the archetype value instead
    #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # archetype level values
    MIS_shape = merge(MIS_shape, distinct(archetype_info[,c('NOMREGION', 'Archetype')]), all=TRUE)
    include_cols = c(which(names(MIS_shape) %in% c('Archetype')), grep('num_total', names(MIS_shape)), grep('num_true', names(MIS_shape)))
    arch_sums = MIS_shape[,include_cols] %>% 
      group_by(Archetype) %>%
      summarise_all(sum, na.rm = TRUE)
    for(var in vaccine_variables){
      if(paste0(var,'_num_true') %in% colnames(arch_sums)){
        arch_sums[[paste0(var, '_rate')]] = arch_sums[[paste0(var,'_num_true')]] / arch_sums[[paste0(var,'_num_total')]]
      }
    }
    
    for(var in vaccine_variables){
      if(paste0(var,'_num_true') %in% colnames(admin_sums)){
        for(i_admin in 1:nrow(admin_sums)){
          if(admin_sums[[paste0(var,'_num_total')]][i_admin]<min_num_total){
            admin_sums[[paste0(var,'_num_total')]][i_admin] = arch_sums[[paste0(var,'_num_total')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
            admin_sums[[paste0(var,'_num_true')]][i_admin] = arch_sums[[paste0(var,'_num_true')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
            admin_sums[[paste0(var,'_rate')]][i_admin] = arch_sums[[paste0(var,'_rate')]][arch_sums$Archetype == admin_sums$Archetype[i_admin]]
          }
        }
      }
    }
    
    # save different version with specified minimum sample size in admin before using region value
    write.csv(as.data.frame(admin_sums), paste0(hbhi_dir, '/estimates_from_DHS/MICS_vaccine_admin_minN',min_num_total,'_', year, '.csv'))
  }
}






















# get country-level estimates using 1) unweighted, aggregated values from entire country and 2) DHS-household weightings, which I believe is meant to be representative within a province
extract_country_level_DHS_data = function(hbhi_dir, dta_dir, years){
  
  for(yy in 1:length(years)){
    year = years[yy]
    DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',year,'_files_recodes_for_sims.csv'))
    country_weighted_rates_df = data.frame('admin' = 'country')
    country_unweighted_rates_df = data.frame('admin' = 'country')
    for (var_index in 1:length(DHS_file_recode_df$variable)){
      if((!is.na(DHS_file_recode_df$filename[var_index])) & !(DHS_file_recode_df$variable[var_index] %in% c('locations', 'iptp_doses'))){
        cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
        if(grepl('itn',DHS_file_recode_df$variable[var_index])){
          # repeat for several age groups
          group_names = c('u5','5_10','10_15','15_20','o20','all')
          lower_age = c(-0.1,5,10,15,20,-0.1)
          upper_age = c(5,10,15,20,200,200)
          for(i_age in 1:length(group_names)){
            cur_dta_age = cur_dta
            var_rate_name = paste0('itn_',group_names[i_age],'_rate')
            cur_dta_age = cur_dta_age[intersect(which(cur_dta_age[[DHS_file_recode_df$age_code[var_index]]]>lower_age[i_age]), which(cur_dta_age[[DHS_file_recode_df$age_code[var_index]]]<=upper_age[i_age])),]
            cur_dta_age$pos = NA
            cur_dta_age$pos[which(cur_dta_age[,grep(colnames(cur_dta_age), pattern = paste0('\\b',DHS_file_recode_df$code[var_index],'\\b'))] == DHS_file_recode_df$pos_pattern[var_index])] = 1
            cur_dta_age$pos[which(cur_dta_age[,grep(colnames(cur_dta_age), pattern = paste0('\\b',DHS_file_recode_df$code[var_index],'\\b'))] == DHS_file_recode_df$neg_pattern[var_index])] = 0
            colnames(cur_dta_age)[which(colnames(cur_dta_age) =='hv005')] = 'clust_weight'
            colnames(cur_dta_age)[which(colnames(cur_dta_age) =='v005')] = 'clust_weight'
            cur_dta_age$clust_weight = cur_dta_age$clust_weight/1000000    
            # tabulate weighted mean for entire country
            country_weighted_rates_df[[var_rate_name]] = weighted.mean(cur_dta_age$pos, cur_dta_age$clust_weight, na.rm=TRUE)
            # # tabulate unweighted mean for entire country
            country_unweighted_rates_df[[var_rate_name]] = mean(cur_dta_age$pos, na.rm=TRUE)
          }
        } else{
          var_rate_name = paste0(DHS_file_recode_df$variable[var_index],'_rate')
          cur_dta$pos = NA
          cur_dta$pos[which(cur_dta[,grep(colnames(cur_dta), pattern = paste0('\\b',DHS_file_recode_df$code[var_index],'\\b'))] == DHS_file_recode_df$pos_pattern[var_index])] = 1
          cur_dta$pos[which(cur_dta[,grep(colnames(cur_dta), pattern = paste0('\\b',DHS_file_recode_df$code[var_index],'\\b'))] == DHS_file_recode_df$neg_pattern[var_index])] = 0
          colnames(cur_dta)[which(colnames(cur_dta) =='hv005')] = 'clust_weight'
          colnames(cur_dta)[which(colnames(cur_dta) =='v005')] = 'clust_weight'
          cur_dta$clust_weight = cur_dta$clust_weight/1000000  
          # tabulate weighted mean for entire country
          country_weighted_rates_df[[var_rate_name]] = weighted.mean(cur_dta$pos, cur_dta$clust_weight, na.rm=TRUE)
          country_weighted_rates_df$year = years[yy]
          # # tabulate unweighted mean for entire country
          country_unweighted_rates_df[[var_rate_name]] = mean(cur_dta$pos, na.rm=TRUE)
          country_unweighted_rates_df$year = years[yy]
        }
      }
    }
    if(yy==1){
      country_weighted_rates_df_all = country_weighted_rates_df
      country_unweighted_rates_df_all = country_unweighted_rates_df
    } else{  # combine current year entries with existing data frame, adding new columns as needed
      country_weighted_rates_df_all = merge(country_weighted_rates_df_all, country_weighted_rates_df, by=intersect(colnames(country_weighted_rates_df), colnames(country_weighted_rates_df_all)), all.x=TRUE, all.y=TRUE)
      country_unweighted_rates_df_all = merge(country_unweighted_rates_df_all, country_unweighted_rates_df, by=intersect(colnames(country_unweighted_rates_df), colnames(country_unweighted_rates_df_all)), all.x=TRUE, all.y=TRUE)
    }
    # write.csv(as.data.frame(country_weighted_rates_df), paste0(hbhi_dir, '/estimates_from_DHS/DHS_country_weighted_rates_', years[yy], '.csv'))
    # write.csv(as.data.frame(country_unweighted_rates_df), paste0(hbhi_dir, '/estimates_from_DHS/DHS_country_unweighted_rates_', years[yy], '.csv'))
  }
  write.csv(as.data.frame(country_weighted_rates_df_all), paste0(hbhi_dir, '/estimates_from_DHS/DHS_country_weighted_rates.csv'))
  write.csv(as.data.frame(country_unweighted_rates_df_all), paste0(hbhi_dir, '/estimates_from_DHS/DHS_country_unweighted_rates.csv'))
}




# get archetype-level estimates using aggregated raw counts within boundaries
extract_archetype_level_DHS_data = function(hbhi_dir, dta_dir, ds_pop_df_filename, years){
  # read in admin-level DHS outputs, then aggregate counts up to archetype level (using archetype-admin relations from pop file) and find archetype-level rates (not weighted)
  # archetype - admin associations
  archetype_info = read.csv(ds_pop_df_filename)
  archetype_names = unique(archetype_info$seasonality_archetype)

  for(yy in 1:length(years)){
    year = years[yy]
    # read in admin-level counts from DHS
    DHS_counts = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_',year,'.csv'))
    # get rid of the rate columns
    DHS_counts = DHS_counts[,-grep('rate',colnames(DHS_counts))]
    # get the DHS variables for this year
    variables0 = colnames(DHS_counts)[grep('_num_total', colnames(DHS_counts))]
    variables = gsub('_num_total','', variables0)
    for(i_arch in 1:length(archetype_names)){
      admins_in_arch = archetype_info$admin_name[which(archetype_info$seasonality_archetype == archetype_names[i_arch])]
      DHS_counts_arch = DHS_counts[which(DHS_counts$NOMDEP %in% admins_in_arch),]
      DHS_counts_arch = DHS_counts_arch[,grep('_num', colnames(DHS_counts_arch))] 
      # take column sums to get totals in archetype, then calculate archetype-level rate and save in data frame
      DHS_colsums = colSums(DHS_counts_arch, na.rm=TRUE)
      DHS_aggregate_arch_cur = data.frame('archetype' = archetype_names[i_arch], 'year' = year)
      for(i_var in 1:length(variables)){
        DHS_aggregate_arch_cur[[paste0(variables[i_var], '_rate')]] = DHS_colsums[[paste0(variables[i_var], '_num_true')]] / DHS_colsums[[paste0(variables[i_var], '_num_total')]]
      }

      if((i_arch == 1) & (yy == 1)){
        DHS_arch_aggregates = DHS_aggregate_arch_cur
      } else{
        DHS_arch_aggregates = merge(DHS_arch_aggregates, DHS_aggregate_arch_cur, by=intersect(colnames(DHS_arch_aggregates), colnames(DHS_aggregate_arch_cur)), all.x=TRUE, all.y=TRUE) # use merge instead of rbind because sometimes later years have more columns
      }
    }
  }
  write.csv(as.data.frame(DHS_arch_aggregates), paste0(hbhi_dir, '/estimates_from_DHS/DHS_archetype_rates.csv'), row.names=FALSE)
}








# plot coverage/prevalence values in each cluster and admin, as extracted from DHS
plot_extracted_DHS_data = function(hbhi_dir, years, admin_shape, min_num_total=30, variables=c('mic', 'rdt','itn_all','itn_u5','itn_5_10','itn_10_15','itn_15_20','itn_o20','iptp','cm','received_treatment','sought_treatment','blood_test'), colors_range_0_to_1=NA, all_years_int_plot_panel=TRUE, separate_plots_for_each_var=TRUE, plot_separate_pdfs=FALSE, plot_vaccine=FALSE, plot_suffix=''){
  
  if(any(is.na(colors_range_0_to_1))){
    colors_range_0_to_1 = add.alpha(pals::parula(101), alpha=0.5)
  }
  
  if(plot_vaccine){
    vacc_string='vaccine_'
  } else vacc_string=''

  #=========================================================================================================##
  # plots of cluster-level DHS results
  #=========================================================================================================##
  if(!dir.exists(paste0(hbhi_dir,'/estimates_from_DHS/plots'))) dir.create(paste0(hbhi_dir,'/estimates_from_DHS/plots'))
  if(plot_separate_pdfs){
    for(yy in 1:length(years)){
      pdf(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_',vacc_string, 'cluster_observations_', years[yy], '.pdf'), width=7, height=5, useDingbats = FALSE)
      par(mfrow=c(1,1))
      for(i_var in 1:length(variables)){
        var = variables[i_var]
        cluster_obs = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'cluster_outputs_', years[yy], '.csv'))[,-1]
        cluster_obs$latitude[which(cluster_obs$latitude == 0)] = NA
        cluster_obs$longitude[which(cluster_obs$longitude == 0)] = NA
        if(paste0(var,'_num_total') %in% colnames(cluster_obs)){
          max_survey_size = max(cluster_obs[[paste0(var,'_num_total')]], na.rm=TRUE)
          layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
          plot(admin_shape, main=var)
          points(cluster_obs$longitude, cluster_obs$latitude, col=colors_range_0_to_1[1+round(cluster_obs[[paste0(var,'_rate')]]*100)], pch=20, cex=cluster_obs[[paste0(var,'_num_total')]]/round(max_survey_size/5))#, xlim=c(min(cluster_obs$longitude), max(cluster_obs$longitude)), ylim=c(min(cluster_obs$latitude), max(cluster_obs$latitude)))
          # legend - colorbar
          legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
          plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
          text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
          rasterImage(legend_image, 0, 0, 1,1)
          # legend - survey size
          plot(rep(0,5), seq(1, max_survey_size, length.out=5), cex=seq(1,max_survey_size, length.out=5)/round(max_survey_size/5), pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)
        }
      }
      dev.off()
    }
  }



  ##=========================================================================================================##
  # plots of cluster-level DHS results - separated by interventions, each plot panel showing across years
  ##=========================================================================================================##
  if(all_years_int_plot_panel){
    if(!dir.exists(paste0(hbhi_dir,'/estimates_from_DHS/plots'))) dir.create(paste0(hbhi_dir,'/estimates_from_DHS/plots'))
    # use subset of variables
    if(!plot_vaccine) {
      variables2 = variables[variables %in% c('mic', 'rdt', 'itn_all', 'itn_u5', 'iptp','cm','received_treatment','sought_treatment','blood_test', 'art')]
    } else variables2 = variables
    nyears = length(years)

    if(separate_plots_for_each_var){
      for(i_var in 1:length(variables2)){
        var = variables2[i_var]
        png(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_',var,'_',vacc_string, 'cluster_observations_all_years',plot_suffix,'.png'), width=0.5*7*nyears, height=0.5*6*1, units='in', res=900)
        base_layout_matrix = matrix(c(rep(1:nyears, each=3), nyears+1, rep(1:nyears, each=3),nyears+2),nrow=2, byrow=TRUE)
        layout(base_layout_matrix)
        par(mar=c(0,0,1,0))
        for(yy in 1:nyears){
          cluster_obs = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'cluster_outputs_', years[yy], '.csv'))[,-1]
          cluster_obs$latitude[which(cluster_obs$latitude == 0)] = NA
          cluster_obs$longitude[which(cluster_obs$longitude == 0)] = NA
          if(paste0(var,'_num_total') %in% colnames(cluster_obs)){
            max_survey_size = max(cluster_obs[[paste0(var,'_num_total')]], na.rm=TRUE)
            plot(admin_shape, main=paste0(var, ' - ', years[yy]), border=rgb(0.5,0.5,0.5,0.5))
            points(cluster_obs$longitude, cluster_obs$latitude, col=colors_range_0_to_1[1+round(cluster_obs[[paste0(var,'_rate')]]*100)], pch=20, cex=cluster_obs[[paste0(var,'_num_total')]]/round(max_survey_size/5))#, xlim=c(min(cluster_obs$longitude), max(cluster_obs$longitude)), ylim=c(min(cluster_obs$latitude), max(cluster_obs$latitude)))
          }else{
            plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
          }
        }
        # legend - colorbar
        legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
        text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        rasterImage(legend_image, 0, 0, 1,1)
        # legend - survey size
        plot(rep(0,5), seq(1, max_survey_size, length.out=5), cex=seq(1,max_survey_size, length.out=5)/round(max_survey_size/5), pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)
        dev.off()
      }
    } else{
      # pdf(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_cluster_observations_all_years2.pdf'), width=28, height=6*length(variables), useDingbats = FALSE)
      png(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_',vacc_string, 'cluster_observations_all_years',plot_suffix,'.png'), width=0.5*7*nyears, height=0.5*6*length(variables2), units='in', res=900)
      base_layout_matrix = matrix(c(rep(1:nyears, each=3), nyears+1, rep(1:nyears, each=3),nyears+2),nrow=2, byrow=TRUE)
      layout_matrix = base_layout_matrix
      if (length(variables2)>1){
        for(vv in 2:length(variables2)){
          layout_matrix = rbind(layout_matrix, base_layout_matrix + (vv-1)*max(base_layout_matrix))
        }
      }
      layout(layout_matrix)
      par(mar=c(0,0,1,0))
      for(i_var in 1:length(variables2)){
        var = variables2[i_var]
        for(yy in 1:nyears){
          cluster_obs = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'cluster_outputs_', years[yy], '.csv'))[,-1]
          cluster_obs$latitude[which(cluster_obs$latitude == 0)] = NA
          cluster_obs$longitude[which(cluster_obs$longitude == 0)] = NA
          if(paste0(var,'_num_total') %in% colnames(cluster_obs)){
            max_survey_size = max(cluster_obs[[paste0(var,'_num_total')]], na.rm=TRUE)
            plot(admin_shape, main=paste0(var, ' - ', years[yy]), border=rgb(0.5,0.5,0.5,0.5))
            points(cluster_obs$longitude, cluster_obs$latitude, col=colors_range_0_to_1[1+round(cluster_obs[[paste0(var,'_rate')]]*100)], pch=20, cex=cluster_obs[[paste0(var,'_num_total')]]/round(max_survey_size/5))#, xlim=c(min(cluster_obs$longitude), max(cluster_obs$longitude)), ylim=c(min(cluster_obs$latitude), max(cluster_obs$latitude)))
          }else{
            plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
          }
        }
        # legend - colorbar
        legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
        text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        rasterImage(legend_image, 0, 0, 1,1)
        # legend - survey size
        plot(rep(0,5), seq(1, max_survey_size, length.out=5), cex=seq(1,max_survey_size, length.out=5)/round(max_survey_size/5), pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)
      }
      dev.off()
    }
  }
  #
  #
  #
  # ####=========================================================================================================####
  # # map of LGA-level DHS results, allowing for aggregation to admin1 level when sample sizes too small
  # ####=========================================================================================================####
  if(plot_separate_pdfs){
    for(yy in 1:length(years)){
      pdf(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_',vacc_string, 'admin_minN', min_num_total,'_', years[yy], '.pdf'), width=7, height=5, useDingbats = FALSE)
      for(i_var in 1:length(variables)){
        var = variables[i_var]
        admin_sums0 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'admin_minN', min_num_total,'_', years[yy], '.csv'))[,-1]
        reorder_admins = match(sapply(admin_shape$NOMDEP, match_lga_names), sapply(admin_sums0$NOMDEP, match_lga_names))
        admin_sums = admin_sums0[reorder_admins,]
        if(all(sapply(admin_shape$NOMDEP, match_lga_names) == sapply(admin_sums$NOMDEP, match_lga_names))){
          if(paste0(var,'_num_total') %in% colnames(admin_sums)){
            layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
            admin_colors = colors_range_0_to_1[1+round(admin_sums[[paste0(var,'_rate')]]*100)]
            plot(admin_shape, main=var, col=admin_colors)

            # legend - colorbar
            legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
            plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
            text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
            rasterImage(legend_image, 0, 0, 1,1)
          }
        } else warning('during plot generation, order of districts in shapefile and data frame did not match, skipping plotting.')
      }
      dev.off()
    }
  }






  
  ##=========================================================================================================##
  # map of LGA-level DHS results, allowing for aggregation to archetype level when sample sizes too small
  #     separated by interventions, each plot panel showing across years
  ##=========================================================================================================##  
  if(all_years_int_plot_panel){
    if(!dir.exists(paste0(hbhi_dir,'/estimates_from_DHS/plots'))) dir.create(paste0(hbhi_dir,'/estimates_from_DHS/plots'))
    # use subset of variables
    if(!plot_vaccine) {
      variables2 =  c(variables[variables %in% c('mic', 'rdt', 'itn_all', 'itn_u5', 'iptp','cm','received_treatment','sought_treatment','blood_test','art')], 'effective_treatment')
    } else variables2 = variables
    nyears = length(years)
    
    if(separate_plots_for_each_var){
      for(i_var in 1:length(variables2)){
        var = variables2[i_var]
        png(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_',var,'_',vacc_string, 'admin_minN',min_num_total,'_observations_acrossYears',plot_suffix,'.png'), width=0.5*7*nyears, height=0.5*6*length(variables2), units='in', res=900)
        base_layout_matrix = matrix(c(rep(1:nyears, each=3), nyears+1, rep(1:nyears, each=3),nyears+2),nrow=2, byrow=TRUE)
        layout(base_layout_matrix)
        par(mar=c(0,0,1,0))
        for(yy in 1:nyears){
          
          # admin_sums0 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_minN',min_num_total,'_includeArch_', years[yy], '.csv'))[,-1]
          admin_sums0 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'admin_minN', min_num_total,'_', years[yy], '.csv'))[,-1]
          reorder_admins = match(sapply(admin_shape$NOMDEP, match_lga_names), sapply(admin_sums0$NOMDEP, match_lga_names))
          admin_sums = admin_sums0[reorder_admins,]
          if(all(sapply(admin_shape$NOMDEP, match_lga_names) == sapply(admin_sums$NOMDEP, match_lga_names))){
            if(paste0(var,'_num_total') %in% colnames(admin_sums)){
              admin_colors = colors_range_0_to_1[1+round(admin_sums[[paste0(var,'_rate')]]*100)]
              plot(admin_shape, main=paste0(var, ' - ', years[yy]), border=rgb(0,0,0,0.5), col=admin_colors)
              
            }else {
              plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
            }
          } else {
            plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
            warning('during plot generation, order of districts in shapefile and data frame did not match, skipping plotting.')
          }
        }
        
        # legend - colorbar
        legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
        text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        rasterImage(legend_image, 0, 0, 1,1)
        plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
        dev.off()
      }
    } else{
      # pdf(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_cluster_observations_all_years2.pdf'), width=28, height=6*length(variables), useDingbats = FALSE)
      png(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_',vacc_string, 'admin_minN',min_num_total,'_observations_acrossYears',plot_suffix,'.png'), width=0.5*7*nyears, height=0.5*6*length(variables2), units='in', res=900)
      base_layout_matrix = matrix(c(rep(1:nyears, each=3), nyears+1, rep(1:nyears, each=3),nyears+2),nrow=2, byrow=TRUE)
      layout_matrix = base_layout_matrix
      if (length(variables2)>1){
        for(vv in 2:length(variables2)){
          layout_matrix = rbind(layout_matrix, base_layout_matrix + (vv-1)*max(base_layout_matrix))
        }
      }
      layout(layout_matrix)
      par(mar=c(0,0,1,0))
      for(i_var in 1:length(variables2)){
        var = variables2[i_var]
        for(yy in 1:nyears){
          
          # admin_sums0 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_minN',min_num_total,'_includeArch_', years[yy], '.csv'))[,-1]
          admin_sums0 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'admin_minN', min_num_total,'_', years[yy], '.csv'))[,-1]
          reorder_admins = match(sapply(admin_shape$NOMDEP, match_lga_names), sapply(admin_sums0$NOMDEP, match_lga_names))
          admin_sums = admin_sums0[reorder_admins,]
          if(all(sapply(admin_shape$NOMDEP, match_lga_names) == sapply(admin_sums$NOMDEP, match_lga_names))){
            if(paste0(var,'_num_total') %in% colnames(admin_sums)){
              admin_colors = colors_range_0_to_1[1+round(admin_sums[[paste0(var,'_rate')]]*100)]
              plot(admin_shape, main=paste0(var, ' - ', years[yy]), border=rgb(0,0,0,0.5), col=admin_colors)
  
            }else {
              plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
            }
          } else {
            plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
            warning('during plot generation, order of districts in shapefile and data frame did not match, skipping plotting.')
          }
        }
        
        # legend - colorbar
        legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
        text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        rasterImage(legend_image, 0, 0, 1,1)
        plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
      }
      dev.off()
    }
  }
}





# plot coverage/prevalence values in each cluster and admin, as extracted from DHS
plot_extracted_DHS_clusters_within_state = function(hbhi_dir, years, admin_shape, state, variables=c('mic', 'rdt','itn_all','itn_u5','itn_5_10','itn_10_15','itn_15_20','itn_o20','iptp','cm','received_treatment','sought_treatment','blood_test'), colors_range_0_to_1=NA, plot_vaccine=FALSE, plot_suffix=''){
  admin_shape = admin_shape[admin_shape$State==state,]
  if(any(is.na(colors_range_0_to_1))){
    colors_range_0_to_1 = add.alpha(pals::parula(101), alpha=0.5)
  }
  if(plot_vaccine){
    vacc_string='vaccine_'
  } else vacc_string=''
  if(!dir.exists(paste0(hbhi_dir,'/estimates_from_DHS/plots'))) dir.create(paste0(hbhi_dir,'/estimates_from_DHS/plots'))
  
  ##=========================================================================================================##
  # plots of cluster-level DHS results - separated by interventions, each plot panel showing across years
  ##=========================================================================================================##
  if(!dir.exists(paste0(hbhi_dir,'/estimates_from_DHS/plots'))) dir.create(paste0(hbhi_dir,'/estimates_from_DHS/plots'))
  # use subset of variables
  if(!plot_vaccine) {
    variables2 = variables[variables %in% c('mic', 'rdt', 'itn_all', 'itn_u5', 'iptp','cm','received_treatment','sought_treatment','blood_test', 'art')]
  } else variables2 = variables
  nyears = length(years)
  
  # pdf(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_cluster_observations_all_years2.pdf'), width=28, height=6*length(variables), useDingbats = FALSE)
  png(paste0(hbhi_dir, '/estimates_from_DHS/plots/',state,'_DHS_',vacc_string, 'cluster_observations_all_years',plot_suffix,'.png'), width=0.5*7*nyears, height=0.5*6*length(variables2), units='in', res=900)
  base_layout_matrix = matrix(c(rep(1:nyears, each=3), nyears+1, rep(1:nyears, each=3),nyears+2),nrow=2, byrow=TRUE)
  layout_matrix = base_layout_matrix
  if (length(variables2)>1){
    for(vv in 2:length(variables2)){
      layout_matrix = rbind(layout_matrix, base_layout_matrix + (vv-1)*max(base_layout_matrix))
    }
  }
  layout(layout_matrix)
  par(mar=c(0,0,1,0))
  for(i_var in 1:length(variables2)){
    var = variables2[i_var]
    for(yy in 1:nyears){
      cluster_obs = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',vacc_string, 'cluster_outputs_', years[yy], '.csv'))[,-1]
      cluster_obs = cluster_obs[cluster_obs$NOMREGION == state,]
      cluster_obs$latitude[which(cluster_obs$latitude == 0)] = NA
      cluster_obs$longitude[which(cluster_obs$longitude == 0)] = NA
      if(paste0(var,'_num_total') %in% colnames(cluster_obs)){
        max_survey_size = max(cluster_obs[[paste0(var,'_num_total')]], na.rm=TRUE)
        plot(st_geometry(admin_shape), main=paste0(var, ' - ', years[yy]), border=rgb(0.5,0.5,0.5,0.5))
        points(cluster_obs$longitude, cluster_obs$latitude, col=colors_range_0_to_1[1+round(cluster_obs[[paste0(var,'_rate')]]*100)], pch=20, cex=cluster_obs[[paste0(var,'_num_total')]]/round(max_survey_size/5))#, xlim=c(min(cluster_obs$longitude), max(cluster_obs$longitude)), ylim=c(min(cluster_obs$latitude), max(cluster_obs$latitude)))
      }else{
        plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
      }
    }
    # legend - colorbar
    legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = var)
    text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
    rasterImage(legend_image, 0, 0, 1,1)
    # legend - survey size
    plot(rep(0,5), seq(1, max_survey_size, length.out=5), cex=seq(1,max_survey_size, length.out=5)/round(max_survey_size/5), pch=20, axes=FALSE, xlab='', ylab='sample size'); axis(2)
  }
  dev.off()
    
  
}









################################################################################################
# get rates of ITN use in different age groups at country level relative to all-age use rate
################################################################################################
get_relative_itn_use_by_age = function(hbhi_dir, years){
  if(!dir.exists(paste0(hbhi_dir,'/estimates_from_DHS/plots'))) dir.create(paste0(hbhi_dir,'/estimates_from_DHS/plots'))
  pdf(paste0(hbhi_dir, '/estimates_from_DHS/plots/DHS_relative_ITN_use_by_age.pdf'), width=10, height=8, useDingbats = FALSE)
  par(mfrow=c(2,2))
  
  save_weighted_ratios = matrix(NA, nrow=0,ncol=5)
  save_equal_weights = matrix(NA, nrow=0,ncol=5)
  for(yy in 1:length(years)){
    clust_sums = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_cluster_outputs_', years[yy], '.csv'))[,-1]
    clust_sums$rel_u5 = clust_sums$itn_u5_rate / clust_sums$itn_all_rate
    clust_sums$rel_5_10 = clust_sums$itn_5_10_rate / clust_sums$itn_all_rate
    clust_sums$rel_10_15 = clust_sums$itn_10_15_rate / clust_sums$itn_all_rate
    clust_sums$rel_15_20 = clust_sums$itn_15_20_rate / clust_sums$itn_all_rate
    clust_sums$rel_o20 = clust_sums$itn_o20_rate / clust_sums$itn_all_rate
    
    
    # remove any clusters with Inf
    clust_sums_2 = clust_sums[which(clust_sums$rel_5_10<Inf),]
    
    # unweighted mean of all cluster values
    unweighted_cluster_means = c(
      mean(clust_sums_2$rel_u5, na.rm=TRUE),
      mean(clust_sums_2$rel_5_10, na.rm=TRUE),
      mean(clust_sums_2$rel_10_15, na.rm=TRUE),
      mean(clust_sums_2$rel_15_20, na.rm=TRUE),
      mean(clust_sums_2$rel_o20, na.rm=TRUE)
    )
    # weighted mean of cluster values
    weighted_cluster_means = c(
      weighted.mean(clust_sums_2$rel_u5, clust_sums_2$itn_weights, na.rm=TRUE),
      weighted.mean(clust_sums_2$rel_5_10, clust_sums_2$itn_weights, na.rm=TRUE),
      weighted.mean(clust_sums_2$rel_10_15, clust_sums_2$itn_weights, na.rm=TRUE),
      weighted.mean(clust_sums_2$rel_15_20, clust_sums_2$itn_weights, na.rm=TRUE),
      weighted.mean(clust_sums_2$rel_o20, clust_sums_2$itn_weights, na.rm=TRUE)
    )
    # without aggregating or weighting by cluster, equal weight to all individuals included in survey
    unweighted_survey_means = c(
      (sum(clust_sums$itn_u5_num_true) / sum(clust_sums$itn_u5_num_total) ) / (sum(clust_sums$itn_all_num_true) / sum(clust_sums$itn_all_num_total)),
      (sum(clust_sums$itn_5_10_num_true) / sum(clust_sums$itn_5_10_num_total) ) / (sum(clust_sums$itn_all_num_true) / sum(clust_sums$itn_all_num_total)),
      (sum(clust_sums$itn_10_15_num_true) / sum(clust_sums$itn_10_15_num_total) ) / (sum(clust_sums$itn_all_num_true) / sum(clust_sums$itn_all_num_total)),
      (sum(clust_sums$itn_15_20_num_true) / sum(clust_sums$itn_15_20_num_total) ) / (sum(clust_sums$itn_all_num_true) / sum(clust_sums$itn_all_num_total)),
      (sum(clust_sums$itn_o20_num_true) / sum(clust_sums$itn_o20_num_total) ) / (sum(clust_sums$itn_all_num_true) / sum(clust_sums$itn_all_num_total))
    )
    save_weighted_ratios = rbind(save_weighted_ratios, weighted_cluster_means)
    save_equal_weights = rbind(save_equal_weights, unweighted_survey_means)
    ratios = matrix(c(unweighted_cluster_means, weighted_cluster_means, unweighted_survey_means), byrow=TRUE, nrow=3)
    barplot(ratios, beside=TRUE, col=c(rgb(0.3,0.3,1), rgb(0.6,0.6,1), rgb(0.9,0.9,1)), names.arg=c('U5','5-10','10-15','15-20','>20'), xlab='age group',ylab='coverage relative to all-age', main= paste0('net usage relative to all-age usage - ', years[yy]))
  }
  plot(NA,ylab='', xlab='', xlim=c(0,1), ylim=c(0,1),axes=FALSE)
  legend('left', c('unweighted mean of clusters', 'DHS-weighted mean of clusters', 'equal weights for all individuals'), col=c(rgb(0.3,0.3,1), rgb(0.6,0.6,1), rgb(0.9,0.9,1)), pch=15, bty='n')
  par(mfrow=c(1,1))
  dev.off()
  
  rownames(save_weighted_ratios) = years
  colnames(save_weighted_ratios) = paste0('relative_use_', c('u5','5_10','10_15','15_20','o20'))
  save_weighted_ratios_df = as.data.frame(save_weighted_ratios)
  
  rownames(save_equal_weights) = years
  colnames(save_equal_weights) = paste0('relative_use_', c('u5','5_10','10_15','15_20','o20'))
  save_equal_weights_df = as.data.frame(save_equal_weights)
  
  write.csv(save_weighted_ratios_df, paste0(hbhi_dir, '/estimates_from_DHS/DHS_weighted_relative_ITN_use_by_age.csv'))
  write.csv(save_weighted_ratios_df, paste0(hbhi_dir, '/estimates_from_DHS/DHS_unweighted_relative_ITN_use_by_age.csv'))
  

  # =====================
  # save the 2016 use rates relative to age U5 from 2016, with cluster-weighted values
  clust_sums = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_cluster_outputs_', years[length(years)], '.csv'))[,-1]
  clust_sums$rel_u5 = clust_sums$itn_u5_rate / clust_sums$itn_u5_rate
  clust_sums$rel_5_10 = clust_sums$itn_5_10_rate / clust_sums$itn_u5_rate
  clust_sums$rel_10_15 = clust_sums$itn_10_15_rate / clust_sums$itn_u5_rate
  clust_sums$rel_15_20 = clust_sums$itn_15_20_rate / clust_sums$itn_u5_rate
  clust_sums$rel_o20 = clust_sums$itn_o20_rate / clust_sums$itn_u5_rate

  # remove any clusters with Inf
  clust_sums_2 = clust_sums[which(clust_sums$rel_5_10<Inf),]
  clust_sums_2 = clust_sums_2[which(clust_sums_2$rel_10_15<Inf),]
  clust_sums_2 = clust_sums_2[which(clust_sums_2$rel_15_20<Inf),]
  clust_sums_2 = clust_sums_2[which(clust_sums_2$rel_o20<Inf),]

  # weighted mean of cluster values
  weighted_cluster_means_rel_u5 = c(
    weighted.mean(clust_sums_2$rel_u5, clust_sums_2$itn_weights, na.rm=TRUE),
    weighted.mean(clust_sums_2$rel_5_10, clust_sums_2$itn_weights, na.rm=TRUE),
    weighted.mean(clust_sums_2$rel_10_15, clust_sums_2$itn_weights, na.rm=TRUE),
    weighted.mean(clust_sums_2$rel_15_20, clust_sums_2$itn_weights, na.rm=TRUE),
    weighted.mean(clust_sums_2$rel_o20, clust_sums_2$itn_weights, na.rm=TRUE)
  )
  weighted_cluster_means_rel_u5[weighted_cluster_means_rel_u5>1]=1
  weighted_cluster_means_rel_u5_mat = matrix(weighted_cluster_means_rel_u5, nrow=1)
  colnames(weighted_cluster_means_rel_u5_mat) = paste0(c('u5','5_10','10_15','15_20','o20'),'_use_relative_to_u5')
  write.csv(weighted_cluster_means_rel_u5_mat, paste0(hbhi_dir, '/estimates_from_DHS/DHS_weighted_ITN_use_by_age_relative_to_u5.csv'))
  
  
  # # =====================  
  # # save the 2016 use rates relative to age >20 from 2016, with cluster-weighted values
  # clust_sums = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_cluster_outputs_', years[length(years)], '.csv'))[,-1]
  # clust_sums$rel_u5 = clust_sums$itn_u5_rate / clust_sums$itn_o20_rate
  # clust_sums$rel_5_10 = clust_sums$itn_5_10_rate / clust_sums$itn_o20_rate
  # clust_sums$rel_10_15 = clust_sums$itn_10_15_rate / clust_sums$itn_o20_rate
  # clust_sums$rel_15_20 = clust_sums$itn_15_20_rate / clust_sums$itn_o20_rate
  # clust_sums$rel_o20 = clust_sums$itn_o20_rate / clust_sums$itn_o20_rate
  # 
  # # remove any clusters with Inf
  # clust_sums_2 = clust_sums[which(clust_sums$rel_5_10<Inf),]
  # clust_sums_2 = clust_sums_2[which(clust_sums_2$rel_10_15<Inf),]
  # clust_sums_2 = clust_sums_2[which(clust_sums_2$rel_15_20<Inf),]
  # clust_sums_2 = clust_sums_2[which(clust_sums_2$rel_u5<Inf),]
  # 
  # # weighted mean of cluster values
  # weighted_cluster_means_rel_o20 = c(
  #   weighted.mean(clust_sums_2$rel_u5, clust_sums_2$itn_weights, na.rm=TRUE),
  #   weighted.mean(clust_sums_2$rel_5_10, clust_sums_2$itn_weights, na.rm=TRUE),
  #   weighted.mean(clust_sums_2$rel_10_15, clust_sums_2$itn_weights, na.rm=TRUE),
  #   weighted.mean(clust_sums_2$rel_15_20, clust_sums_2$itn_weights, na.rm=TRUE),
  #   weighted.mean(clust_sums_2$rel_o20, clust_sums_2$itn_weights, na.rm=TRUE)
  # )
  # weighted_cluster_means_rel_o20[weighted_cluster_means_rel_o20>1]=1
  # weighted_cluster_means_rel_o20_mat = matrix(weighted_cluster_means_rel_o20, nrow=1)
  # colnames(weighted_cluster_means_rel_o20_mat) = paste0(c('u5','5_10','10_15','15_20','o20'),'_use_relative_to_o20')
  # write.csv(save_weighted_ratios_df, paste0(hbhi_dir, '/estimates_from_DHS/DHS_weighted_ITN_use_by_age_relative_to_o20.csv'))
  # 
  
  
  

}












# ####################    old     #################################
# 
# 
# 
# 
# par(mfrow=c(3,4))
# # for each of the interventions/measures, count number of individuals surveyed in each admin
# num_surveyed = as.data.frame(MIS_2016_shape[c('NOMDEP','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
#   group_by(NOMDEP) %>%
#   summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
#             num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
#             num_w_fever=sum(num_w_fever, na.rm = TRUE),
#             num_preg=sum(num_preg, na.rm = TRUE))
# num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
# hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,800, length.out=80))
# hist(num_surveyed$num_use_itn_all,  main='ITN', breaks=seq(0,1600, length.out=80))
# hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,200, length.out=80))
# hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,500, length.out=80))
# 
# 
# # look at survey numbers by region
# # for each of the interventions/measures, count number of individuals surveyed in each admin
# num_surveyed = as.data.frame(MIS_2016_shape[c('NOMREGION','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
#   group_by(NOMREGION) %>%
#   summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
#             num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
#             num_w_fever=sum(num_w_fever, na.rm = TRUE),
#             num_preg=sum(num_preg, na.rm = TRUE))
# num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
# hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,800, length.out=80))
# hist(num_surveyed$num_use_itn_all,  main='ITN', breaks=seq(0,1600, length.out=80))
# hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,200, length.out=80))
# hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,500, length.out=80))
# 
# # look at survey numbers by region
# # for each of the interventions/measures, count number of individuals surveyed in each admin
# num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_1','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
#   group_by(NAME_1) %>%
#   summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
#             num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
#             num_w_fever=sum(num_w_fever, na.rm = TRUE),
#             num_preg=sum(num_preg, na.rm = TRUE))
# num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
# hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,3000, length.out=80))
# hist(num_surveyed$num_use_itn_all,  main='ITN',breaks=seq(0,6000, length.out=80))
# hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,700, length.out=80))
# hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,2000, length.out=80))
# par(mfrow=c(1,1))
# 
# 
# 
# 
