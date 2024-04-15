##########################################################################
#############                Survival analysis               ############# 
###########################################################################
# written by Tak Karasaki (takahiro.karasaki@crick.ac.uk)

# Description:
# Script to create Figure 4 & Extended figure 10 of the manuscript "The evolution of lung cancer and impact of subclonal selection in TRACERx"

############################################################################
### Important notes for the survival analysis (DFS) in TRACERx421 cohort ###
############################################################################
# As described in the Methods of the above manuscript, 
#
# 1. One patient (CRUK0682) with synchronous primary lung cancers (LUAD & LUSC) whose tumour with the highest stage (LUAD) 
#    was not sequenced was excluded from the survival analysis. [line 131]
# 2. As for the patients who harboured synchronous multiple primary lung cancers, when associating genomic data from the tumours 
#   with patient level survival information, we used only data from the tumour of the highest pathological TMN stage.
#   "tumour_id_per_patient" column in "per patient" object gives you the tumour_id of the tumour of the highest stage. 
# 3. The tumour of the highest stage could not be selected for the genomically confirmed collision tumours (CRUK0704 and CRUK0881)*.  
#    For these cases, the genomic variable of interest need to be linked case-by-case (e.g. lines 175 - 181). 
# 4. During the follow-up, 4 patients (CRUK0512, CRUK0373, CRUK0428, CRUK0511) developed new primary cancer and subsequent to 
#    this a recurrence from either the 1st primary lung cancer or the new primary cancer diagnosed during the follow up. 
#    These 4 cases were censored at the time of the diagnosis of new primary cancer for DFS analysis, due to the uncertainty 
#    of the origin of the third tumour. [lines 135 - 138]
# *CRUK0039 was also a genomically identified collision tumour but only one "tumour" (genomically identified cluster) could run 
# through the whole WES pipeline and therefore that information was used for the downstream analysis.
 

####################################### 
### libraries and working directory ###
#######################################

library(tidyverse)
library(data.table)
library(readxl)
library(fst)
library(RColorBrewer)
library(ggfittext)
library(cowplot)
library(ggpubr)
library(survival)
library(survminer)
library(car)
library(muhaz)
library(survRM2)
library(rstatix)
library(ggforestplot)

setwd("./")

maindir <- ("./Fig4.and.Supp/")
dir.create(maindir)

########################
## path to input data ##
########################

# clinical data
all_patient_df_path <- '../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_patient_df.rds'
all_tumour_df_path <- '../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_tumour_df.rds'

# genomic data outputs
evo_metrics_path <- "../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_evolutionary_metrics.tsv"
gd_df_tumour_path <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_gd_table_per_tumour.tsv"

# clone diversity metrics (default tree version)
diversity_df1_path <- "../20221109_Tumour_evo_histories_DATA/20221102_TRACERx421_tumourDivCorrectedtree.csv"

# clone diversity metrics (tree with minimum score version)
diversity_df2_path <- "../20221109_Tumour_evo_histories_DATA/20221102_TRACERx421_tumourDivMintree_subclonalExpScore.csv"

#############################
#### Run all analyses !! ####
#############################

## Main figure: use default tree (i.e. tree likely to be less branching than alternative trees)
 # Supp figure: use tree that yields lowest subclonal expansion score (probably more conservative way to define the subclonal expansion score)

analysis_version_list <- c("main", "supp")

for (analysis_version in analysis_version_list){ 
  
  ###############################
  ## load clone diversity data ##
  ###############################
  
  if(analysis_version == "main"){
    
    workdir <- paste0(maindir, analysis_version,"/")
    dir.create(workdir)
    
    diversity_df <- fread(diversity_df1_path) %>%
      mutate(min_simpson_tum = min_simpson_tum_meanccf) %>%
      mutate(max_CCF_terminal_tum = maxCCF_terminalonly_tum_meanccf) %>%
      dplyr::select(patient_id, tumour_id,min_simpson_tum,max_CCF_terminal_tum  )
    
  } else if (analysis_version == "supp"){
    
    workdir <- paste0(maindir, analysis_version,"/")
    dir.create(workdir)
    
    diversity_df1 <- fread(diversity_df1_path)
    diversity_df2 <- fread(diversity_df2_path)
    
    # glimpse(diversity_df1)
    # glimpse(diversity_df2)
    
    diversity_df <- diversity_df1 %>% 
      left_join(diversity_df2[,c("tumour_id", "mintree_maxCCF_terminalonly_tum_meanccf")], by = "tumour_id") %>%
      dplyr::select(-min_simpson_tum) %>%
      mutate(min_simpson_tum = min_simpson_tum_meanccf) %>%
      mutate(max_CCF_terminal_tum = mintree_maxCCF_terminalonly_tum_meanccf)  %>%
      dplyr::select(patient_id, tumour_id,min_simpson_tum,max_CCF_terminal_tum  )
  }
    
  ##########################################  
  ### create input for survival analysis ###
  ##########################################  
  
  ## load & organise clinical data 
  
  all_patient_df <- readRDS(all_patient_df_path )
  all_tumour_df <- readRDS(all_tumour_df_path )
  
  clinical_data_df <- all_patient_df  %>%
    dplyr::select(cruk_id,tumour_id_per_patient, tx100, age, pathologyTNM ,packyears, smoking_status_merged, 
                  histology_lesion1,adjuvant_treatment_YN,
                  cens_os, os_time, cens_dfs, dfs_time, cens_dfs_any_event, dfs_time_any_event, cens_lung_event, lung_event_time, 
                  Recurrence_time_use, Relapse_cat_new
    ) %>%
    filter(!cruk_id %in% c("CRUK0682")) %>%  # no lesion 1 sampled : CRUK0682
    mutate(tumour_id_per_patient = case_when(cruk_id == "CRUK0704" ~  "CRUK0704_Tumour1.2.3",  # these are collision tumours
                                             cruk_id == "CRUK0881" ~  "CRUK0881_Tumour1.2",
                                             TRUE ~ tumour_id_per_patient)) %>%
    mutate(cens_dfs = case_when(cruk_id %in% c("CRUK0512", "CRUK0373","CRUK0428","CRUK0511") ~ 0, 
                                TRUE ~ as.numeric(cens_dfs)),
           dfs_time = case_when(cruk_id %in% c("CRUK0512", "CRUK0373","CRUK0428","CRUK0511") ~ dfs_time_any_event, 
                                TRUE ~ dfs_time )) %>%
    mutate(Relapse_cat_new = case_when(cruk_id %in% c("CRUK0173","CRUK0249", "CRUK0781","CRUK0129") ~ "Unknown Site", # "CRUK0173" : relapsed before death but no CT reports available. "CRUK0249", "CRUK0781","CRUK0129": lung cancer death as first event 
                                       TRUE ~ as.character(Relapse_cat_new) )) %>%
    mutate(Relapse_cat_new = factor(Relapse_cat_new, levels = c("No rec", "Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site"))) %>%
    mutate(is.rec = case_when(Relapse_cat_new %in% c("Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site") ~ TRUE, 
                              TRUE ~ FALSE)) %>%
    mutate(time.to.rec = case_when(is.rec == TRUE ~ Recurrence_time_use)) %>%
    mutate(Relapse.Site = case_when(Relapse_cat_new %in% c("Intrathoracic") ~ "Intrathoracic",
                                    Relapse_cat_new %in% c("Intra & Extra", "Extrathoracic") ~ "Extrathoracic")) %>%
    mutate(Relapse.Site = factor(Relapse.Site , levels = c("Intrathoracic", "Extrathoracic"))) %>%
    
    mutate(TNM_stage = case_when(pathologyTNM %in% c("IA","IB") ~ "I",
                                 pathologyTNM %in% c("IIA","IIB") ~ "II",
                                 pathologyTNM %in% c("IIIA","IIIB") ~ "III")) %>%
    mutate(pTNM_stage = case_when(pathologyTNM %in% c("IIIA","IIIB") ~ "III",
                                  TRUE ~ pathologyTNM)) %>%
    mutate(adjuvant_tx = factor(adjuvant_treatment_YN, levels = c("No adjuvant", "Adjuvant"))) %>%
    mutate(pack_years = packyears,
           smoking_history = factor(smoking_status_merged, levels = c("Never Smoked","Ex-Smoker",  "Smoker" ))) %>%
    mutate(histology = case_when(histology_lesion1 == "Invasive adenocarcinoma" ~ "LUAD",
                                 histology_lesion1 == "Squamous cell carcinoma" ~ "LUSC",
                                 TRUE ~ "Other")) %>%
    mutate(histology_LUAD = case_when(histology_lesion1 == "Invasive adenocarcinoma" ~ "LUAD",
                                      TRUE ~ "Non-LUAD")) %>%
    mutate(histology = factor(histology, levels = c("LUAD", "LUSC", "Other")),
           histology_LUAD = factor(histology_LUAD, levels = c("LUAD","Non-LUAD")),
           adjuvant_tx = factor(adjuvant_treatment_YN, levels = c("No adjuvant", "Adjuvant")),
           smoking_history = factor(smoking_status_merged, levels = c("Never Smoked","Ex-Smoker",  "Smoker" )))  %>%
    dplyr::select(-pathologyTNM, -histology_lesion1, -adjuvant_treatment_YN, -smoking_status_merged )
  
  ### load various ITH metrics
  
  ## SCNA-ITH  
  evo_metrics <- fread(evo_metrics_path) %>%
    dplyr::select(tumour_id, frac_abberant_genom_subcl, frac_gain_genom_subcl, frac_loss_genom_subcl, perc_subclonal ,clonal_mixing_index,
                  clonal_tmb, subclonal_tmb) 
  
  # collision tumours:  CRUK0704, CRUK0881 (CRUK0039 also had collision tumour, though only 1 of them was analysed in primary cohort )
  # use max SCNA-ITH (frac_abberant_genom_subcl) to represent the patient's SCNA-ITH, based on NEJM2017 that high SCNA-ITH is associated with poor prognosis
  evo_metrics_collision <- evo_metrics %>%  
    mutate(patient_id = str_split_fixed(tumour_id, "_", n =2)[,1]) %>%
    filter(patient_id %in%  c("CRUK0704", "CRUK0881")) %>% 
    group_by(patient_id) %>%
    summarise(frac_abberant_genom_subcl = max(frac_abberant_genom_subcl))
  
  ## load subclonal GD data 
  gd_df_tumour_ori <- fread(gd_df_tumour_path)
  
  gd_df_tumour <- gd_df_tumour_ori %>%
    dplyr::select(tumour_id, num_clonal_gds, num_subclonal_gds) 
  
  gds <- data.table(gd_df_tumour)[,c('tumour_id', 'num_clonal_gds', 'num_subclonal_gds')]
  gds[ num_clonal_gds == 0 & num_subclonal_gds == 0, gd_status := 'No GD']
  gds[ num_clonal_gds == 1  & num_subclonal_gds == 0,  gd_status := '1 Truncal GD']
  gds[ num_clonal_gds == 2 & num_subclonal_gds == 0,  gd_status := '2 Truncal GD']
  gds[ num_clonal_gds == 0 & num_subclonal_gds == 1, gd_status := '1 Subclonal GD']
  gds[ num_clonal_gds == 0 & num_subclonal_gds > 1,  gd_status := 'Multi-Subclonal GD']
  gds[ num_clonal_gds == 1 & num_subclonal_gds == 1,  gd_status := '1 Truncal GD & Subclonal GD']
  gds[ num_clonal_gds == 1 & num_subclonal_gds > 1,   gd_status := '1 Truncal GD & Multi-Subclonal GD']
  
  gd_df_tumour <- gd_df_tumour %>%
    left_join(gds[,c("tumour_id", "gd_status")])
  
  
  ## organise ITH metrics data , deal with collision tumours
  
  tumour_ITH_df <- diversity_df %>%
    filter(tumour_id != "CRUK0721_Tumour1", tumour_id != "CRUK0555_Tumour1") %>%  # these are single region tumours in multiple-tumour patients
    distinct(tumour_id , .keep_all = T) %>%  # no need this but no need to remove this
    dplyr::select(tumour_id , 
                  min_simpson_tum, 
                  max_CCF_terminal_tum, 
    ) %>% 
    left_join(evo_metrics, by = "tumour_id") %>%
    left_join(gd_df_tumour, by = "tumour_id") 
  
  
  tumour_ITH_df_collision <- tumour_ITH_df  %>%
    mutate(tumour_id = case_when(grepl("CRUK0704", tumour_id) ~ "CRUK0704_Tumour1.2.3", 
                                 grepl("CRUK0881", tumour_id) ~ "CRUK0881_Tumour1.2", 
                                 TRUE ~ tumour_id)) %>%
    group_by(tumour_id, gd_status) %>%  # gd_status   CRUK0704: 1 Truncal GD x3   CRUK0881: Multi-Subclonal GD x2 
    summarise(min_simpson_tum = min(min_simpson_tum, na.rm = T),
              max_CCF_terminal_tum = max(max_CCF_terminal_tum, na.rm = T),
              frac_abberant_genom_subcl = max(frac_abberant_genom_subcl, na.rm = T),
              perc_subclonal = max(perc_subclonal, na.rm = T),
              num_subclonal_gds = max(num_subclonal_gds, na.rm = T),
              num_clonal_gds = max(num_clonal_gds, na.rm = T))  %>%
    mutate(subclonal_gd = case_when(is.na(num_subclonal_gds) ~ NA_character_,
                                    num_subclonal_gds > 0 ~ "Yes",
                                    num_subclonal_gds == 0 ~ "No" ),
           any_gd = case_when(is.na(num_subclonal_gds) ~ NA_character_,
                              num_subclonal_gds > 0 | num_clonal_gds > 0 ~ "Yes",
                              num_subclonal_gds == 0 & num_clonal_gds == 0 ~ "No" ))
  
  
  # this will be the analytical cohort  n=420 (CRUK0682 excluded from survival analysis)
  in_df_surv_pre <- clinical_data_df %>%
    left_join(tumour_ITH_df_collision, by = c("tumour_id_per_patient" = "tumour_id") )
  
  ### set median and quartile cutoff  per analytical cohort
  
  # median cutoff
  min_simpson_tum_median <- median(in_df_surv_pre$min_simpson_tum, na.rm = T)
  max_CCF_terminal_tum_median <- median(in_df_surv_pre$max_CCF_terminal_tum, na.rm = T)
  frac_abberant_genom_subcl_median <- median(in_df_surv_pre$frac_abberant_genom_subcl, na.rm = T)
  perc_subclonal_median <- median(in_df_surv_pre$perc_subclonal, na.rm = T)
  
  # tertile cutoff 
  max_CCF_terminal_tum_Q3 <- quantile(in_df_surv_pre$max_CCF_terminal_tum, probs = c(1/3,2/3), na.rm = T)
  
  # quartile cutoff 
  max_CCF_terminal_tum_Q4 <- quantile(in_df_surv_pre$max_CCF_terminal_tum, probs = c(1/4,2/4, 3/4), na.rm = T)
  frac_abberant_genom_subcl_Q4 <- quantile(in_df_surv_pre$frac_abberant_genom_subcl, probs = c(1/4,2/4, 3/4), na.rm = T)
  perc_subclonal_Q4 <- quantile(in_df_surv_pre$perc_subclonal, probs = c(1/4,2/4, 3/4), na.rm = T)
  
  
  ## categorise each variables
  in_df_surv <- in_df_surv_pre %>%
    mutate(min_simpson_tum_cat = case_when(min_simpson_tum >= min_simpson_tum_median  ~ "high",
                                           min_simpson_tum < min_simpson_tum_median  ~ "low"),
           max_CCF_terminal_tum_cat = case_when(max_CCF_terminal_tum >= max_CCF_terminal_tum_median  ~ "high",
                                                max_CCF_terminal_tum < max_CCF_terminal_tum_median  ~ "low"),
           frac_abberant_genom_subcl_cat = case_when(frac_abberant_genom_subcl >= frac_abberant_genom_subcl_median  ~ "high",
                                                     frac_abberant_genom_subcl < frac_abberant_genom_subcl_median  ~ "low"),
           perc_subclonal_cat = case_when(perc_subclonal >= perc_subclonal_median  ~ "high",
                                          perc_subclonal < perc_subclonal_median  ~ "low"))%>%
    mutate(max_CCF_terminal_tum_3=  case_when(max_CCF_terminal_tum >= max_CCF_terminal_tum_Q3[2]  ~ "high",
                                              max_CCF_terminal_tum >= max_CCF_terminal_tum_Q3[1] &  max_CCF_terminal_tum < max_CCF_terminal_tum_Q3[2]  ~ "mid",
                                              max_CCF_terminal_tum < max_CCF_terminal_tum_Q3[1]  ~ "low"),
    )%>%
    mutate(max_CCF_terminal_tum_4=  case_when(max_CCF_terminal_tum >= max_CCF_terminal_tum_Q4[3]  ~ "highest",
                                              max_CCF_terminal_tum >= max_CCF_terminal_tum_Q4[2] &  max_CCF_terminal_tum < max_CCF_terminal_tum_Q4[3]  ~ "high",
                                              max_CCF_terminal_tum >= max_CCF_terminal_tum_Q4[1] &  max_CCF_terminal_tum < max_CCF_terminal_tum_Q4[2]  ~ "low",
                                              max_CCF_terminal_tum < max_CCF_terminal_tum_Q4[1]  ~ "lowest"),
           frac_abberant_genom_subcl_4=  case_when(frac_abberant_genom_subcl >= frac_abberant_genom_subcl_Q4[3]  ~ "highest",
                                                   frac_abberant_genom_subcl >= frac_abberant_genom_subcl_Q4[2] &  frac_abberant_genom_subcl < frac_abberant_genom_subcl_Q4[3]  ~ "high",
                                                   frac_abberant_genom_subcl >= frac_abberant_genom_subcl_Q4[1] &  frac_abberant_genom_subcl < frac_abberant_genom_subcl_Q4[2]  ~ "low",
                                                   frac_abberant_genom_subcl < frac_abberant_genom_subcl_Q4[1]  ~ "lowest"),
           perc_subclonal_4=  case_when(perc_subclonal >= perc_subclonal_Q4[3]  ~ "highest",
                                        perc_subclonal >= perc_subclonal_Q4[2] &  perc_subclonal < perc_subclonal_Q4[3]  ~ "high",
                                        perc_subclonal >= perc_subclonal_Q4[1] &  perc_subclonal < perc_subclonal_Q4[2]  ~ "low",
                                        perc_subclonal < perc_subclonal_Q4[1]  ~ "lowest"))%>%
    mutate(min_simpson_tum_cat  = factor(min_simpson_tum_cat ,levels = c("low","high")),
           max_CCF_terminal_tum_cat  = factor(max_CCF_terminal_tum_cat ,levels = c("low","high")),
           frac_abberant_genom_subcl_cat  = factor(frac_abberant_genom_subcl_cat ,levels = c("low","high")),
           max_CCF_terminal_tum_3 = factor(max_CCF_terminal_tum_3, levels =c("low","mid","high")),
           max_CCF_terminal_tum_4 = factor(max_CCF_terminal_tum_4, levels =c("lowest","low","high","highest")),
           frac_abberant_genom_subcl_4 = factor(frac_abberant_genom_subcl_4, levels =c("lowest","low","high","highest")),
           perc_subclonal_4 = factor(perc_subclonal_4, levels =c("lowest","low","high","highest")))
  
  
  # For DFS, censor 4 patients who are currently marked as "recurrence" but uncertain whether the recurrence is from 1st primary or 2nd primary (if from 2nd primary, these cases should NOT be marked as recurrence in TRACERx protocol)
  #  "CRUK0512", "CRUK0373","CRUK0428","CRUK0511" : currently marked as recurrence, after 2nd primary cancer was confirmed
  
  in_data_surv <- in_df_surv %>%
    mutate(subclonal_gd.3 = case_when(subclonal_gd == "Yes" ~ "Subclonal GD",
                                      any_gd == "Yes" ~ "Clonal GD only", 
                                      any_gd == "No" ~ "No GD") ) %>%
    mutate(subclonal_gd.3 = factor(subclonal_gd.3, levels = c("No GD", "Subclonal GD","Clonal GD only"))) %>%
    mutate(SCNA.ITH = frac_abberant_genom_subcl_cat,
           mut.ITH = perc_subclonal_cat,
           RecentSubclExpansionScore_2 =  max_CCF_terminal_tum_cat , 
           RecentSubclExpansionScore_3 =  max_CCF_terminal_tum_3, 
           RecentSubclExpansionScore_4 =  max_CCF_terminal_tum_4 , 
           RecentSubclExpansionScore0.3 =  max_CCF_terminal_tum*(10/3), # to get OR per 0.3 increase (1 sd = 0.27)
           RecentSubclExpansionScore =  max_CCF_terminal_tum*5, # to get OR per 0.2 increase - but no reason to use "per 0.2"
           minSimpson_2 = min_simpson_tum_cat, 
           minSimpson = min_simpson_tum*5, # to get OR per 0.2 increase (1 sd = 0.22)
           Age = age/10,  # to get OR per 10 increase
           PackYears = packyears/10)  # to get OR per 10 increase
  
  
  
  
  ########################################################
  ## check median survival time & median follow up time ##
  ########################################################
  
  in_df <- all_patient_df %>%
    mutate(cens_os_rev = case_when(cens_os == 1 ~ 0,
                                   cens_os == 0 ~ 1))
  
  # median survival time
  survfit(Surv(os_time, cens_os) ~ 1, data = in_df)
  #        n events median 0.95LCL 0.95UCL
  # [1,] 421    177   2197    1767      NA
  2197/365.25
  # [1] 6.015058 # years
  2197/365.25 * 12 
  # [1] 72.1807 # months
  
  # min 
  min(in_df$os_time) / 365.25 * 12 
  # [1] 0.03285421
  # max
  max(in_df$os_time) / 365.25 * 12 
  # [1] 79.86858
  
  # median follow-up time  by reverse KM method
  survfit(Surv(os_time, cens_os_rev) ~ 1, data = in_df)
  #        n events median 0.95LCL 0.95UCL
  # [1,] 421    244   1702    1649    1784
  1702/365.25
  # [1] 4.659822  # years
  1702/365.25 * 12 
  # [1] 55.91786 # months
  
  # min disease free
  min(in_df$os_time[ in_df$cens_dfs_any_event == 0 ]) / 365.25 * 12 
  # [1] 2.168378
  # max disease free
  max(in_df$os_time[ in_df$cens_dfs_any_event == 0 ]) / 365.25 * 12 
  # [1] 73.72485

  
  ################
  ################
  ###   main   ###  
  ################
  ################
  
  time_vec <- c( "dfs_time")
  status_vec <- c("cens_dfs")
  title_vec <- c( "Disease-free survival")
  
  #########################
  ### survival analysis ###  uni-Cox  & KM plot  for 2 groups
  #########################
  
  
  group_vec <- c("RecentSubclExpansionScore_2","SCNA.ITH","any_gd","subclonal_gd", "mut.ITH")
  out_pdf <- paste0(workdir,"/kmplot_all.patients_by2.pdf")
  pdf(out_pdf, width = 5, height = 7)
  
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      # out_f <- paste0(group_vec[n],"_",time_vec[m],".kmplot.png")
      # png(out_f)
      
      survival_data.frame <- in_data_surv 
      colnames_group <-  group_vec[n]
      colnames_time <-  time_vec[m]
      colnames_status <-  status_vec[m]
      title_phrase <- title_vec[m]
      xlim_value <- 2000
      xbreak_value <- 500
      
      surv_df <- survival_data.frame
      surv_df$group <- surv_df[[colnames_group]]
      surv_df$time <- surv_df[[colnames_time]]
      surv_df$status <- surv_df[[colnames_status]]  
      
      if(colnames_group %in% c("any_gd", "subclonal_gd")){
        surv_df <- surv_df %>%
          dplyr::filter(!is.na(group), !is.na(status)) %>%
          mutate(group = factor(group, levels = c("No", "Yes")))
      } else {
        surv_df <- surv_df %>%
          dplyr::filter(!is.na(group), !is.na(status)) %>%
          mutate(group = factor(group, levels = c("low", "high")))
      }
      
      fit <- survfit(Surv(time, status) ~ group, data= surv_df)
      cox_fit <- coxph( Surv(time, status) ~ group, 
                        data = surv_df )
      CI_tab <- summary(cox_fit)$conf.int
      round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
      
      if(colnames_group == "any_gd"){
        label <- paste0('Hazard ratio:\n',
                        'Yes vs No - ', round_cust(CI_tab[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab[, 3][1]), '-', round_cust(CI_tab[, 4][1]),
                        '\n (P=', signif(summary(cox_fit)$coefficients[,5][1],2),')')
        labs_vec <-  c('No', 'Yes')
      } else {
        label <- paste0('Hazard ratio:\n',
                        'High vs Low - ', round_cust(CI_tab[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab[, 3][1]), '-', round_cust(CI_tab[, 4][1]),
                        # '\nPolycl. & Polyphyl. vs Monocl. - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
                        '\n (P=', signif(summary(cox_fit)$coefficients[,5][1],2),')')
        labs_vec <-  c('Low', 'High')
      }
      
      colours <- brewer.pal(3, 'Set1')[1:2]
      
      g <- ggsurvplot(
        fit = fit,
        pval = FALSE, #TRUE,
        pval.method = FALSE, #TRUE,
        #test.for.trend = TRUE, 
        risk.table = TRUE,
        #risk.table.height = 0.16,
        legend.title = '',
        font.legend = 16,
        palette = rev(colours),
        xlab = "Days", 
        ylab = title_phrase,
        xlim = c(0,xlim_value), 
        font.tickslab = c(18),
        font.x = c(20),
        font.y = c(20),
        break.x.by = xbreak_value,
        legend.labs = labs_vec,
        title = paste(title_vec[m],"by",group_vec[n]) ) 
      
      g$plot <- g$plot +
        annotate("text", hjust = 0,
                 x = 10, y = 0.15, # x and y coordinates of the text
                 label = label, size = 4)
      g$table <- g$table + 
        labs(title = '', y = '', x = '')  #theme(plot.title = element_blank())
      
      print(g)
      #dev.off()
    }
  }
  dev.off()
  
  
  #########################
  ### survival analysis ###   uni-Cox  &  KM plot  for 3 clonal dominance 
  #########################
  
  group_vec <- c("RecentSubclExpansionScore_3" )
  out_pdf <- paste0(workdir,"/kmplot_all.patients_by3.pdf")
  
  pdf(out_pdf, width = 5, height = 7)
  
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      
      survival_data.frame <-in_data_surv
      colnames_group <-  group_vec[n]
      colnames_time <-  time_vec[m]
      colnames_status <-  status_vec[m]
      title_phrase <- title_vec[m]
      xlim_value <- 2000
      xbreak_value <- 500
      
      surv_df <- survival_data.frame
      surv_df$group <- surv_df[[colnames_group]]
      surv_df$time <- surv_df[[colnames_time]]
      surv_df$status <- surv_df[[colnames_status]]  
      
      surv_df <- surv_df %>%
        dplyr::filter(!is.na(group), !is.na(status)) %>%
        mutate(group = factor(group, levels = c("low","mid","high")))%>%
        mutate(group_2 = case_when(group == "high" ~ "high",
                                   TRUE ~ "all other")) %>%
        mutate(group_2 = factor(group_2, levels = c("all other", "high")))
      
      round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
      
      fit <- survfit(Surv(time, status) ~ group, data= surv_df)
      cox_fit <- coxph( Surv(time, status) ~ group, 
                        data = surv_df )
      CI_tab <- summary(cox_fit)$conf.int
      
      cox_fit2 <- coxph( Surv(time, status) ~ group_2, 
                         data = surv_df )
      CI_tab2 <- summary(cox_fit2)$conf.int
      
      label <- paste0('Hazard ratio:\n',
                      'High vs Low - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
                      ' (P=', signif(summary(cox_fit)$coefficients[,5][2],2),')',
                      '\nHigh vs all other - ', round_cust(CI_tab2[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab2[, 3][1]), '-', round_cust(CI_tab2[, 4][1]),
                      ' (P=', signif(summary(cox_fit2)$coefficients[,5][1],2),')'
      )
      labs_vec <-  c('Low','Mid', 'High')
      
      colours <- brewer.pal(3, 'Set1')[c(1,3,2)]
      
      
      g <- ggsurvplot(
        fit = fit,
        pval = FALSE, #TRUE,
        pval.method = FALSE, #TRUE,
        #test.for.trend = TRUE, 
        risk.table = TRUE,
        #risk.table.height = 0.16,
        legend.title = '',
        font.legend = 16,
        palette = rev(colours),
        xlab = "Days", 
        ylab = title_phrase,
        xlim = c(0,xlim_value), 
        font.tickslab = c(18),
        font.x = c(20),
        font.y = c(20),
        break.x.by = xbreak_value,
        legend.labs = labs_vec,
        title = paste(title_vec[m],"by",group_vec[n]) ) 
      
      g$plot <- g$plot +
        annotate("text", hjust = 0,
                 x = 10, y = 0.15, # x and y coordinates of the text
                 label = label, size = 3)
      g$table <- g$table + 
        labs(title = '', y = '', x = '')  #theme(plot.title = element_blank())
      
      print(g)
      #dev.off()
    }
  }
  dev.off()
  
  #########################
  ### survival analysis ###  subclonal GD (by 3 groups)
  #########################
  
  group_vec <- c("subclonal_gd.3","subclonal_gd.3.plot")
  out_pdf <- paste0(workdir,"/kmplot_all.patients_subclonalGD_by.3.pdf")
  
  pdf(out_pdf, width = 5, height = 7)
  
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      # out_f <- paste0(group_vec[n],"_",time_vec[m],".kmplot.png")
      # png(out_f)
      
      survival_data.frame <- in_data_surv  %>% 
        mutate(subclonal_gd.3.plot = case_when(subclonal_gd.3 == "No GD" ~ "low",
                                               subclonal_gd.3 == "Clonal GD only" ~ "mid",
                                               subclonal_gd.3 == "Subclonal GD" ~ "high")) %>%
        mutate(subclonal_gd.3.plot = factor(subclonal_gd.3.plot , levels = c("low","mid","high")))
      
      colnames_group <-  group_vec[n]
      colnames_time <-  time_vec[m]
      colnames_status <-  status_vec[m]
      title_phrase <- title_vec[m]
      xlim_value <- 2000
      xbreak_value <- 500
      
      surv_df <- survival_data.frame
      surv_df$group <- surv_df[[colnames_group]]
      surv_df$time <- surv_df[[colnames_time]]
      surv_df$status <- surv_df[[colnames_status]]  
      
      if(colnames_group == "subclonal_gd.3"){
        surv_df <- surv_df %>%
          dplyr::filter(!is.na(group), !is.na(status)) %>%
          mutate(group = factor(group, levels = c("No GD","Clonal GD only", "Subclonal GD")))
        
        surv_df2 <- surv_df %>%
          mutate(group = factor(group, levels = c("Clonal GD only", "Subclonal GD","No GD")))
        
      } else {
        surv_df <- surv_df %>%
          dplyr::filter(!is.na(group), !is.na(status)) %>%
          mutate(group = factor(group, levels = c("low", "mid","high")))
        surv_df2 <- surv_df %>%
          mutate(group = factor(group, levels = c("mid", "high","low")))
      }
      
      fit <- survfit(Surv(time, status) ~ group, data= surv_df)
      cox_fit <- coxph( Surv(time, status) ~ group, 
                        data = surv_df )
      CI_tab <- summary(cox_fit)$conf.int
      
      cox_fit2 <- coxph( Surv(time, status) ~ group, 
                         data = surv_df2 )
      CI_tab2 <- summary(cox_fit2)$conf.int
      
      round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
      
      
      if(colnames_group == "subclonal_gd.3"){
        label <- paste0('Hazard ratio:\n',
                        'Clonal vs No GD - ', round_cust(CI_tab[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab[, 3][1]), '-', round_cust(CI_tab[, 4][1]),
                        ' (P=', signif(summary(cox_fit)$coefficients[,5][1],2),')',
                        '\nSublonal vs No GD - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
                        ' (P=', signif(summary(cox_fit)$coefficients[,5][2],2),')',
                        '\nSublonal vs Clonal - ', round_cust(CI_tab2[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab2[, 3][1]), '-', round_cust(CI_tab2[, 4][1]),
                        ' (P=', signif(summary(cox_fit2)$coefficients[,5][1],2),')')
        
        labs_vec <-  c("No GD","Clonal GD only", "Subclonal GD")
        colours <- c("No GD"='#bababa', "Clonal GD only"='#0000CC', "Subclonal GD"='#990000')
      } else {
        label <- paste0('Hazard ratio:\n',
                        'Clonal vs No GD - ', round_cust(CI_tab[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab[, 3][1]), '-', round_cust(CI_tab[, 4][1]),
                        ' (P=', signif(summary(cox_fit)$coefficients[,5][1],2),')',
                        '\nSublonal vs No GD - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
                        ' (P=', signif(summary(cox_fit)$coefficients[,5][2],2),')',
                        '\nSublonal vs Clonal - ', round_cust(CI_tab2[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab2[, 3][1]), '-', round_cust(CI_tab2[, 4][1]),
                        ' (P=', signif(summary(cox_fit2)$coefficients[,5][1],2),')')
        
        labs_vec <-  c('low','mid', 'high')
        colours <- c('low'='#bababa', 'mid'='#0000CC', 'high'='#990000')
      }
      
      g <- ggsurvplot(
        fit = fit,
        pval = FALSE, #TRUE,
        pval.method = FALSE, #TRUE,
        #test.for.trend = TRUE, 
        risk.table = TRUE,
        #risk.table.height = 0.16,
        legend.title = '',
        font.legend = 16,
        palette = rev(colours),
        xlab = "Days", 
        ylab = title_phrase,
        xlim = c(0,xlim_value), 
        font.tickslab = c(18),
        font.x = c(20),
        font.y = c(20),
        break.x.by = xbreak_value,
        legend.labs = labs_vec,
        title = paste(title_vec[m],"by",group_vec[n]) ) 
      
      g$plot <- g$plot +
        annotate("text", hjust = 0,
                 x = 10, y = 0.15, # x and y coordinates of the text
                 label = label, size = 3)
      g$table <- g$table + 
        labs(title = '', y = '', x = '')  #theme(plot.title = element_blank())
      
      print(g)
      #dev.off()
    }
  }
  dev.off()
  
  #########################
  ### survival analysis ###  uni-Cox  &  KM plot  for 4 clonal dominance 
  #########################
  
  group_vec <- c("RecentSubclExpansionScore_4" )
  out_pdf <- paste0(workdir,"/kmplot_all.patients_by4.pdf")
  
  pdf(out_pdf, width = 5, height = 7)
  
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      
      survival_data.frame <- in_data_surv 
      colnames_group <-  group_vec[n]
      colnames_time <-  time_vec[m]
      colnames_status <-  status_vec[m]
      title_phrase <- title_vec[m]
      xlim_value <- 2000
      xbreak_value <- 500
      
      surv_df <- survival_data.frame
      surv_df$group <- surv_df[[colnames_group]]
      surv_df$time <- surv_df[[colnames_time]]
      surv_df$status <- surv_df[[colnames_status]]  
      
      surv_df <- surv_df %>%
        dplyr::filter(!is.na(group), !is.na(status)) %>%
        mutate(group = factor(group, levels = c("lowest", "low","high", "highest"))) %>%
        mutate(group_2 = case_when(group == "highest" ~ "highest",
                                   TRUE ~ "all other")) %>%
        mutate(group_2 = factor(group_2, levels = c("all other", "highest")))
      
      round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
      
      fit <- survfit(Surv(time, status) ~ group, data= surv_df)
      
      cox_fit <- coxph( Surv(time, status) ~ group, 
                        data = surv_df )
      CI_tab <- summary(cox_fit)$conf.int
      
      cox_fit2 <- coxph( Surv(time, status) ~ group_2, 
                         data = surv_df )
      CI_tab2 <- summary(cox_fit2)$conf.int
      
      label <- paste0('Hazard ratio:\n',
                      'Highest vs Lowest - ', round_cust(CI_tab[, 1][3]), ' ; 95% CI: ',round_cust(CI_tab[, 3][3]), '-', round_cust(CI_tab[, 4][3]),
                      # '\nPolycl. & Polyphyl. vs Monocl. - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
                      '\n (P=', signif(summary(cox_fit)$coefficients[,5][3],2),')',
                      '\nHighest vs all other - ', round_cust(CI_tab2[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab2[, 3][1]), '-', round_cust(CI_tab2[, 4][1]),
                      ' (P=', signif(summary(cox_fit2)$coefficients[,5][1],2),')'
      )
      labs_vec <-  c('Lowest','Low', 'High', "'Highest")
      
      colours <-c("#E31A1C", "#FB9A99", "#A6CEE3", "#1F78B4")
      
      
      g <- ggsurvplot(
        fit = fit,
        pval = FALSE, #TRUE,
        pval.method = FALSE, #TRUE,
        #test.for.trend = TRUE, 
        risk.table = TRUE,
        #risk.table.height = 0.16,
        legend.title = '',
        font.legend = 16,
        palette = rev(colours),
        xlab = "Days", 
        ylab = title_phrase,
        xlim = c(0,xlim_value), 
        font.tickslab = c(18),
        font.x = c(20),
        font.y = c(20),
        break.x.by = xbreak_value,
        legend.labs = labs_vec,
        title = paste(title_vec[m],"by",group_vec[n]) ) 
      
      g$plot <- g$plot +
        annotate("text", hjust = 0,
                 x = 10, y = 0.15, # x and y coordinates of the text
                 label = label, size = 4)
      g$table <- g$table + 
        labs(title = '', y = '', x = '')  #theme(plot.title = element_blank())
      
      print(g)
      #dev.off()
    }
  }
  dev.off()
  
  #########################
  ### survival analysis ###   KM plot  all gd_status
  #########################
  
  group_vec <- c("gd_status" )
  out_pdf <- paste0(workdir,"/kmplot_all.patients_gd_status.pdf")
  pdf(out_pdf, width = 7, height = 7)
  # out_pdf <- paste0(workdir,"/kmplot_all.patients_gd_status.large.pdf")
  # pdf(out_pdf, width = 12, height = 12)
  
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      # out_f <- paste0(group_vec[n],"_",time_vec[m],".kmplot.png")
      # png(out_f)
      
      colours <- c('No GD'='#bababa', '1 Truncal GD'='#9999FF', '2 Truncal GD'='#0000CC', 
                   '1 Truncal GD & Subclonal GD'='#CC66FF', '1 Truncal GD & Multi-Subclonal GD'='#330066', 
                   '1 Subclonal GD'='#FF6666', 'Multi-Subclonal GD'='#990000')
      
      survival_data.frame <- in_data_surv  %>%
        mutate(gd_status = factor(gd_status, levels = names(colours)))
      
      colnames_group <-  group_vec[n]
      colnames_time <-  time_vec[m]
      colnames_status <-  status_vec[m]
      title_phrase <- title_vec[m]
      xlim_value <- 2000
      xbreak_value <- 500
      
      surv_df <- survival_data.frame
      surv_df$group <- surv_df[[colnames_group]]
      surv_df$time <- surv_df[[colnames_time]]
      surv_df$status <- surv_df[[colnames_status]]  
      
      fit <- survfit(Surv(time, status) ~ group, data= surv_df)
      cox_fit <- coxph( Surv(time, status) ~ group, 
                        data = surv_df )
      CI_tab <- summary(cox_fit)$conf.int
      round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
      
      labs_vec <- names(colours)
      
      
      
      g <- ggsurvplot(
        fit = fit,
        pval = FALSE, #TRUE,
        pval.method = FALSE, #TRUE,
        #test.for.trend = TRUE, 
        risk.table = TRUE,
        #risk.table.height = 0.16,
        legend.title = '',
        font.legend = 16,
        palette = colours,
        xlab = "Days", 
        ylab = title_phrase,
        xlim = c(0,xlim_value), 
        font.tickslab = c(18),
        font.x = c(20),
        font.y = c(20),
        break.x.by = xbreak_value,
        legend.labs = labs_vec,
        title = paste(title_vec[m],"by",group_vec[n]) ) 
      g$table <- g$table + 
        labs(title = '', y = '', x = '')  #theme(plot.title = element_blank())
      
      print(g)
    }
  }
  dev.off()
  
  ##############################
  ## Uni-Cox -> .csv ##  
  ##############################
  
  group_vec <- c("RecentSubclExpansionScore0.3","RecentSubclExpansionScore","RecentSubclExpansionScore_2", "SCNA.ITH", "subclonal_gd", "subclonal_gd.3")  # subclonal_gd.3: ref = No GD
  
  # Uni-Cox
  all_cox_res <- vector()                          
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      survival_data.frame <-in_data_surv 
      colnames_group <-  group_vec[n]
      colnames_time <-  time_vec[m]
      colnames_status <-  status_vec[m]
      
      surv_df <- survival_data.frame
      surv_df$group <- surv_df[[colnames_group]]
      surv_df$time <- surv_df[[colnames_time]]
      surv_df$status <- surv_df[[colnames_status]]  
      
      cox_fit <- coxph( Surv(time, status) ~   group, 
                        data = surv_df )
      CI_tab <- summary(cox_fit)$conf.int
      cox_res <- c(colnames_group,colnames_time,
                   signif(CI_tab[, 1][1],3),signif(CI_tab[, 3][1],3), signif(CI_tab[, 4][1],3),
                   signif(summary(cox_fit)$coefficients[,5][1],3) )
      cox_res_mat <- matrix(cox_res, nrow = 1)
      all_cox_res <- rbind(all_cox_res, cox_res)
    }
  }
  all_cox_df <- data.frame(all_cox_res, stringsAsFactors = F)
  names(all_cox_df) <- c("variable", "time","HR","lower95CI","upper95CI", "pval")
  all_cox_df <- all_cox_df %>% 
    mutate(HR = as.numeric(HR),
           lower95CI= as.numeric(lower95CI),
           upper95CI= as.numeric(upper95CI),
           pval= as.numeric(pval))
  
  write.csv(all_cox_df, paste0(workdir, "/all_unicox_df.csv"), row.names = F)
  
  
  ##########################################
  ## Multi-Cox full model -> forest plot  ## 
  ##########################################
  
  in_df <- in_data_surv
  
  # as numerical
  variables_cox <- c("Age","pTNM_stage", "PackYears","histology_LUAD","adjuvant_tx", "SCNA.ITH", "subclonal_gd","RecentSubclExpansionScore0.3" ) 
  out_pdf <- paste0(workdir,"/forestplot_all.patients_ITHx3.num.pdf")
  
  pdf(out_pdf, width = 10, height = 12)
  
  for(var in c(variables_cox) ) {
    in_df <- in_df[!is.na(in_df[[var]]),]
  }
  f <- as.formula(paste("Surv(time, status)",
                        paste(variables_cox, collapse = " + "), 
                        sep = " ~ "))
  for (m in 1:length(time_vec)) {
    
    survival_data.frame <-in_df
    colnames_time <-  time_vec[m]
    colnames_status <-  status_vec[m]
    title_phrase <- title_vec[m]
    
    surv_df <- survival_data.frame
    surv_df$time <- surv_df[[colnames_time]]
    surv_df$status <- surv_df[[colnames_status]]  
    
    cox_out <- coxph(f,data = surv_df)
    
    fp <- ggforest(cox_out, data = NULL, main = paste0("Hazard ratio: ",title_phrase), fontsize = 1.5, refLabel = "reference", noDigits = 2) +
      ggtitle(title_phrase)
    print(fp)
    #dev.off()
  }
  dev.off()
  
  
  # as categorical
  variables_cox <- c("Age","pTNM_stage", "PackYears","histology_LUAD","adjuvant_tx", "SCNA.ITH", "subclonal_gd","RecentSubclExpansionScore_2" ) 
  out_pdf <- paste0(workdir,"/forestplot_all.patients_ITHx3.cat.pdf")
  
  pdf(out_pdf, width = 10, height = 12)
  
  for(var in c(variables_cox) ) {
    in_df <- in_df[!is.na(in_df[[var]]),]
  }
  f <- as.formula(paste("Surv(time, status)",
                        paste(variables_cox, collapse = " + "), 
                        sep = " ~ "))
  for (m in 1:length(time_vec)) {
    
    survival_data.frame <-in_df
    colnames_time <-  time_vec[m]
    colnames_status <-  status_vec[m]
    title_phrase <- title_vec[m]
    
    surv_df <- survival_data.frame
    surv_df$time <- surv_df[[colnames_time]]
    surv_df$status <- surv_df[[colnames_status]]  
    
    cox_out <- coxph(f,data = surv_df)
    
    fp <- ggforest(cox_out, data = NULL, main = paste0("Hazard ratio: ",title_phrase), fontsize = 1.5, refLabel = "reference", noDigits = 2) +
      ggtitle(title_phrase)
    print(fp)
    #dev.off()
  }
  dev.off()
  
  
  ##########################################################
  ## Multi-Cox : SCNA-ITH vs subclonal GD -> forest plot  ## 
  ##########################################################
  
  in_df <- in_data_surv
  
  variables_cox <- c("SCNA.ITH", "subclonal_gd") 
  out_pdf <- paste0(workdir,"/forestplot_all.patients_subclGD.SCNAITH.pdf")
  
  pdf(out_pdf, width = 8, height = 5)
  
  for(var in c(variables_cox) ) {
    in_df <- in_df[!is.na(in_df[[var]]),]
  }
  f <- as.formula(paste("Surv(time, status)",
                        paste(variables_cox, collapse = " + "), 
                        sep = " ~ "))
  for (m in 1:length(time_vec)) {
    
    survival_data.frame <-in_df 
    colnames_time <-  time_vec[m]
    colnames_status <-  status_vec[m]
    title_phrase <- title_vec[m]
    
    surv_df <- survival_data.frame
    surv_df$time <- surv_df[[colnames_time]]
    surv_df$status <- surv_df[[colnames_status]]  
    
    cox_out <- coxph(f,data = surv_df)
    
    fp <- ggforest(cox_out, data = NULL, main = paste0("Hazard ratio: ",title_phrase), fontsize = 1.5, refLabel = "reference", noDigits = 2) +
      ggtitle(title_phrase)
    print(fp)
    
  }
  dev.off()
  
  
  ########################
  # hazard function plot #
  ########################
  
  library(muhaz) 
  in_df_h <- in_data_surv %>%
    filter(SCNA.ITH == "high")
  
  in_df_l <- in_data_surv %>%
    filter(SCNA.ITH == "low")
  
  # muhaz (hazard estimates using kernel-based method)
  mod.h <- muhaz(in_df_h$dfs_time, in_df_h$cens_dfs ) # status must be 1 for failure and 0 for censored
  mod.l <- muhaz(in_df_l$dfs_time, in_df_l$cens_dfs ) # status must be 1 for failure and 0 for censored
  
  out_pdf <- paste0(workdir,"/hazard.plot_DFS_SCNA.ITH.muhaz.pdf")
  pdf(out_pdf, width = 5, height = 5)
  
  plot.muhaz(mod.h, ylim = c(0, 0.0013), xlim = c(0,2000) , col = "#E41A1C" )
  lines(mod.l, ylim = c(0, 0.0013), xlim = c(0,2000), col = "#377EB8" )
  legend(x=1500, y=0.0012, legend = c("high", "low"),  col=c("#E41A1C","#377EB8"),lty=c(1,1), cex=0.8)
  
  dev.off()
  
  #################################
  ### analyse RMST sequentially ### 
  #################################
  library(survRM2)
  
  ### split the cohort into Tx100 and post-tx100 to assess SCNA-ITH
  
  ## let's create figure x = truncation time, y = RMTL ratio (with 95CI) 
  time_vec <- c("dfs_time") 
  status_vec <- c("cens_dfs") 
  title_vec <- c( "Disease-free survival") 
  file_vec <- c("DFS")
  
  in_data_surv_tx100 <- in_data_surv %>%
    filter(tx100 == TRUE)
  in_data_surv_post.tx100 <- in_data_surv %>%
    filter(tx100 != TRUE)
  
  table(in_data_surv_tx100$SCNA.ITH) # tx100 -> 92 patients
  # low high 
  # 43   49 
  
  table(in_data_surv_post.tx100$SCNA.ITH) # post tx100 -> 300 patients
  # low high 
  # 153  147 
  
  table(in_data_surv$SCNA.ITH) # tx421 -> 392 patients
  # low high 
  # 196  196 
  
  group_vec <- c("SCNA.ITH")
  cov_vec <- c("Age","pTNM_stage","histology_LUAD", "PackYears", "adjuvant_tx")
  
  tau_vec_months <- c(6,12,18,24,30,36,42,48,54,60)
  
  # Tx421 cohort
  workdir_sub <- workdir
  
  in_data_surv_sub <- in_data_surv
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      
      all_unadj.res <- vector()
      all_adj.res <- vector()
      
      for(j in 1:length(tau_vec_months)){
        
        colnames_time <-  time_vec[m]
        colnames_status <-  status_vec[m]
        colnames_group <- group_vec[n]
        
        surv_df  <- in_data_surv_sub[!is.na(in_data_surv_sub[[colnames_group]]),] %>%
          mutate(subclonal_gd = case_when(subclonal_gd == "Yes" ~ "high",
                                          subclonal_gd == "No" ~ "low")) %>%
          mutate(subclonal_gd = factor(subclonal_gd, levels = c("low", "high") )) 
        surv_df$time <- surv_df[[colnames_time]]/(365.25/12) # transform days -> months
        surv_df$cens <- surv_df[[colnames_status]]  
        surv_df$group <- surv_df[[colnames_group]]  
        
        surv_df <- surv_df %>%
          filter(group %in% c("low","high"))
        
        surv_df$group_bin <- ifelse(as.numeric(factor(surv_df[[colnames_group]])) == 1, 0, 1)
        
        time <- surv_df$time
        status <- surv_df$cens
        arm <- surv_df$group_bin 
        
        # unadjusted RMST
        obj = rmst2(time, status, arm, tau=tau_vec_months[j])
        #print(obj)
        obj_table <- obj$unadjusted.result
        unadj.res_vec <- c(tau=tau_vec_months[j], obj_table[3,1], obj_table[3,2],obj_table[3,3], obj_table[3,4])
        
        all_unadj.res <- rbind(all_unadj.res, unadj.res_vec)
        
        # adjusted RMST
        
        for(var in c(cov_vec, colnames_group) ) {
          surv_df <- surv_df[!is.na(surv_df[[var]]),]
        }
        
        cov_df <- surv_df[,cov_vec ]  
        cov_model.matrix <- model.matrix(~ ., cov_df) [,-1]  # remove intercept (1st column)
        
        obj.adj = rmst2(time, status, arm, tau=tau_vec_months[j], covariates= cov_model.matrix)
        obj.adj_table <- obj.adj$adjusted.result
        adj.res_vec <- c(tau=tau_vec_months[j], obj.adj_table[3,1], obj.adj_table[3,2],obj.adj_table[3,3], obj.adj_table[3,4])
        
        all_adj.res <- rbind(all_adj.res, adj.res_vec)
      }
      
      
      #
      
      all_unadj.res_df <- as.data.frame(all_unadj.res)
      all_adj.res_df <- as.data.frame(all_adj.res)
      
      colnames(all_unadj.res_df) <- c("TruncationTime", "RMTLratio", "lower95CI","upper95CI","pval")
      colnames(all_adj.res_df) <- c("TruncationTime", "RMTLratio", "lower95CI","upper95CI","pval")
      
      write.csv(all_unadj.res_df, paste0(workdir_sub,"/all_unadj.res_df_",file_vec[m],"_",group_vec[n],".csv"), row.names = F )
      write.csv(all_adj.res_df, paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".csv"), row.names = F )
      
      # plot! 
      
      in_df <- all_adj.res_df %>%
        filter(!is.nan(RMTLratio),!is.nan(lower95CI)) #%>%
      #pivot_longer(cols = RMTLratio:upper95CI, names_to = "RMTLcat", values_to = "RMTL") %>%
      #filter(!is.nan(RMTL))
      
      g <- ggplot(data = in_df) +
        geom_errorbar(aes(x=TruncationTime, y = RMTLratio, ymin=lower95CI, ymax=upper95CI), width=2, size=0.5, color="black") + 
        geom_point(aes(x=TruncationTime, y=RMTLratio), size=4, shape=21, fill="white") +
        geom_hline(yintercept=1, lty=2) +
        # theme_classic()+
        scale_y_log10() +
        scale_x_continuous(breaks = c(12,24,36,48,60,72)) + 
        #geom_hline(yintercept = 1,  color = "red") + 
        theme_cowplot() +
        labs(x = "Truncation time (months)", y = "RMTL ratio (95%CI)")+
        theme(axis.title.x = element_text(size=18),  
              axis.title.y = element_text(size=18),
              axis.text.x = element_text(size=14),  
              axis.text.y = element_text(size=14),
              legend.position = "none",
              legend.title = element_blank(),
              legend.text = element_text(size=14),
              title = element_text(size=20)) 
      
      out_pdf <- paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".pdf")
      
      pdf(out_pdf, width = 4.5, height = 4.5)
      print(g)
      dev.off()
      
      out_pdff <- paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".eps")
      fp
      ggsave(out_pdff,  width = 4.5, height = 4.5)
      
    }
  }
  
  # Tx100 cohort
  workdir_sub <- paste0(workdir,"/tx100/")
  dir.create(workdir_sub)
  
  in_data_surv_sub <- in_data_surv_tx100
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      
      all_unadj.res <- vector()
      all_adj.res <- vector()
      
      for(j in 1:length(tau_vec_months)){
        
        colnames_time <-  time_vec[m]
        colnames_status <-  status_vec[m]
        colnames_group <- group_vec[n]
        
        surv_df  <- in_data_surv_sub[!is.na(in_data_surv_sub[[colnames_group]]),] %>%
          mutate(subclonal_gd = case_when(subclonal_gd == "Yes" ~ "high",
                                          subclonal_gd == "No" ~ "low")) %>%
          mutate(subclonal_gd = factor(subclonal_gd, levels = c("low", "high") )) 
        surv_df$time <- surv_df[[colnames_time]]/(365.25/12) # transform days -> months
        surv_df$cens <- surv_df[[colnames_status]]  
        surv_df$group <- surv_df[[colnames_group]]  
        
        surv_df <- surv_df %>%
          filter(group %in% c("low","high"))
        
        surv_df$group_bin <- ifelse(as.numeric(factor(surv_df[[colnames_group]])) == 1, 0, 1)
        
        time <- surv_df$time
        status <- surv_df$cens
        arm <- surv_df$group_bin 
        
        # unadjusted RMST
        obj = rmst2(time, status, arm, tau=tau_vec_months[j])
        #print(obj)
        obj_table <- obj$unadjusted.result
        unadj.res_vec <- c(tau=tau_vec_months[j], obj_table[3,1], obj_table[3,2],obj_table[3,3], obj_table[3,4])
        
        all_unadj.res <- rbind(all_unadj.res, unadj.res_vec)
        
        # adjusted RMST
        
        for(var in c(cov_vec, colnames_group) ) {
          surv_df <- surv_df[!is.na(surv_df[[var]]),]
        }
        
        cov_df <- surv_df[,cov_vec ]  
        cov_model.matrix <- model.matrix(~ ., cov_df) [,-1]  # remove intercept (1st column)
        
        obj.adj = rmst2(time, status, arm, tau=tau_vec_months[j], covariates= cov_model.matrix)
        obj.adj_table <- obj.adj$adjusted.result
        adj.res_vec <- c(tau=tau_vec_months[j], obj.adj_table[3,1], obj.adj_table[3,2],obj.adj_table[3,3], obj.adj_table[3,4])
        
        all_adj.res <- rbind(all_adj.res, adj.res_vec)
      }
      
      
      #
      
      all_unadj.res_df <- as.data.frame(all_unadj.res)
      all_adj.res_df <- as.data.frame(all_adj.res)
      
      colnames(all_unadj.res_df) <- c("TruncationTime", "RMTLratio", "lower95CI","upper95CI","pval")
      colnames(all_adj.res_df) <- c("TruncationTime", "RMTLratio", "lower95CI","upper95CI","pval")
      
      write.csv(all_unadj.res_df, paste0(workdir_sub,"/all_unadj.res_df_",file_vec[m],"_",group_vec[n],".csv"), row.names = F )
      write.csv(all_adj.res_df, paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".csv"), row.names = F )
      
      # plot! 
      
      in_df <- all_adj.res_df %>%
        filter(!is.nan(RMTLratio),!is.nan(lower95CI)) #%>%
      #pivot_longer(cols = RMTLratio:upper95CI, names_to = "RMTLcat", values_to = "RMTL") %>%
      #filter(!is.nan(RMTL))
      
      g <- ggplot(data = in_df) +
        geom_errorbar(aes(x=TruncationTime, y = RMTLratio, ymin=lower95CI, ymax=upper95CI), width=2, size=0.5, color="black") + 
        geom_point(aes(x=TruncationTime, y=RMTLratio), size=4, shape=21, fill="white") +
        geom_hline(yintercept=1, lty=2) +
        # theme_classic()+
        scale_y_log10() +
        scale_x_continuous(breaks = c(12,24,36,48,60,72)) + 
        #geom_hline(yintercept = 1,  color = "red") + 
        theme_cowplot() +
        labs(x = "Truncation time (months)", y = "RMTL ratio (95%CI)")+
        theme(axis.title.x = element_text(size=18),   
              axis.title.y = element_text(size=18),
              axis.text.x = element_text(size=14),  
              axis.text.y = element_text(size=14),
              legend.position = "none",
              legend.title = element_blank(),
              legend.text = element_text(size=14),
              title = element_text(size=20)) 
      
      out_pdf <- paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".pdf")
      
      pdf(out_pdf, width = 4.5, height = 4.5)
      print(g)
      dev.off()
      
      out_pdff <- paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".eps")
      fp
      ggsave(out_pdff,  width =  4.5, height = 4.5)
      
    }
  }
  
  # Post-Tx100 cohort
  workdir_sub <- paste0(workdir,"/post.tx100/")
  dir.create(workdir_sub)
  
  in_data_surv_sub <- in_data_surv_post.tx100
  for (n in 1:length(group_vec)) {
    for (m in 1:length(time_vec)) {
      
      all_unadj.res <- vector()
      all_adj.res <- vector()
      
      for(j in 1:length(tau_vec_months)){
        
        colnames_time <-  time_vec[m]
        colnames_status <-  status_vec[m]
        colnames_group <- group_vec[n]
        
        surv_df  <- in_data_surv_sub[!is.na(in_data_surv_sub[[colnames_group]]),] %>%
          mutate(subclonal_gd = case_when(subclonal_gd == "Yes" ~ "high",
                                          subclonal_gd == "No" ~ "low")) %>%
          mutate(subclonal_gd = factor(subclonal_gd, levels = c("low", "high") )) 
        surv_df$time <- surv_df[[colnames_time]]/(365.25/12) # transform days -> months
        surv_df$cens <- surv_df[[colnames_status]]  
        surv_df$group <- surv_df[[colnames_group]]  
        
        surv_df <- surv_df %>%
          filter(group %in% c("low","high"))
        
        surv_df$group_bin <- ifelse(as.numeric(factor(surv_df[[colnames_group]])) == 1, 0, 1)
        
        time <- surv_df$time
        status <- surv_df$cens
        arm <- surv_df$group_bin 
        
        # unadjusted RMST
        obj = rmst2(time, status, arm, tau=tau_vec_months[j])
        #print(obj)
        obj_table <- obj$unadjusted.result
        unadj.res_vec <- c(tau=tau_vec_months[j], obj_table[3,1], obj_table[3,2],obj_table[3,3], obj_table[3,4])
        
        all_unadj.res <- rbind(all_unadj.res, unadj.res_vec)
        
        # adjusted RMST
        
        for(var in c(cov_vec, colnames_group) ) {
          surv_df <- surv_df[!is.na(surv_df[[var]]),]
        }
        
        cov_df <- surv_df[,cov_vec ]  
        cov_model.matrix <- model.matrix(~ ., cov_df) [,-1]  # remove intercept (1st column)
        
        obj.adj = rmst2(time, status, arm, tau=tau_vec_months[j], covariates= cov_model.matrix)
        obj.adj_table <- obj.adj$adjusted.result
        adj.res_vec <- c(tau=tau_vec_months[j], obj.adj_table[3,1], obj.adj_table[3,2],obj.adj_table[3,3], obj.adj_table[3,4])
        
        all_adj.res <- rbind(all_adj.res, adj.res_vec)
      }
      
      
      #
      
      all_unadj.res_df <- as.data.frame(all_unadj.res)
      all_adj.res_df <- as.data.frame(all_adj.res)
      
      colnames(all_unadj.res_df) <- c("TruncationTime", "RMTLratio", "lower95CI","upper95CI","pval")
      colnames(all_adj.res_df) <- c("TruncationTime", "RMTLratio", "lower95CI","upper95CI","pval")
      
      write.csv(all_unadj.res_df, paste0(workdir_sub,"/all_unadj.res_df_",file_vec[m],"_",group_vec[n],".csv"), row.names = F )
      write.csv(all_adj.res_df, paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".csv"), row.names = F )
      
      # plot! 
      
      in_df <- all_adj.res_df %>%
        filter(!is.nan(RMTLratio),!is.nan(lower95CI)) #%>%
      #pivot_longer(cols = RMTLratio:upper95CI, names_to = "RMTLcat", values_to = "RMTL") %>%
      #filter(!is.nan(RMTL))
      
      g <- ggplot(data = in_df) +
        geom_errorbar(aes(x=TruncationTime, y = RMTLratio, ymin=lower95CI, ymax=upper95CI), width=2, size=0.5, color="black") + 
        geom_point(aes(x=TruncationTime, y=RMTLratio), size=4, shape=21, fill="white") +
        geom_hline(yintercept=1, lty=2) +
        # theme_classic()+
        scale_y_log10() +
        scale_x_continuous(breaks = c(12,24,36,48,60,72)) + 
        #geom_hline(yintercept = 1,  color = "red") + 
        theme_cowplot() +
        labs(x = "Truncation time (months)", y = "RMTL ratio (95%CI)")+
        theme(axis.title.x = element_text(size=18),  
              axis.title.y = element_text(size=18),
              axis.text.x = element_text(size=14), 
              axis.text.y = element_text(size=14),
              legend.position = "none",
              legend.title = element_blank(),
              legend.text = element_text(size=14),
              title = element_text(size=20)) 
      
      out_pdf <- paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".pdf")
      
      pdf(out_pdf, width =  4.5, height = 4.5)
      print(g)
      dev.off()
      
      out_pdff <- paste0(workdir_sub,"/all_adj.res_df_",file_vec[m],"_",group_vec[n],".eps")
      fp
      ggsave(out_pdff,  width = 4.5, height = 4.5)
    }
  }
  
  
  #################################################
  ## barplot : SCNA-ITH vs relapse site & timing ##
  #################################################
  
  in_data <-  in_data_surv %>%
    mutate(rec_1yr = case_when(lung_event_time < 365 & cens_lung_event == 1 ~ "rec < 1yr",
                               lung_event_time >= 365 & cens_lung_event == 1 ~ "rec >= 1yr")) %>%
    mutate(rec_1yr = factor(rec_1yr , levels = c("rec < 1yr", "rec >= 1yr"))) %>%
    mutate(rec_1.2.3yr = case_when(lung_event_time < 365 & cens_lung_event == 1 ~ "rec < 1y",
                                   lung_event_time < 365*2 & lung_event_time >= 365 & cens_lung_event == 1 ~ "rec 1-2y",
                                   lung_event_time < 365*3 & lung_event_time >= 365*2 & cens_lung_event == 1 ~ "rec 2-3y",
                                   lung_event_time >= 365*3 & cens_lung_event == 1 ~ "rec >= 3y")) %>%
    mutate(rec_1.2.3yr = factor(rec_1.2.3yr , levels = c("rec < 1y", "rec 1-2y","rec 2-3y", "rec >= 3y"))) %>%
    mutate(Relapse_Site = case_when(Relapse_cat_new =="No rec" ~ "No rec before death or new primary", 
                                    Relapse_cat_new %in% c("Intra & Extra", "Extrathoracic") ~ "Extrathoracic (any)",
                                    Relapse_cat_new == "Intrathoracic" ~ "Intrathoracic (only)")) %>%
    mutate(Relapse_Site = factor(Relapse_Site, levels = c("No rec before death or new primary", "Intrathoracic (only)", "Extrathoracic (any)"))) %>%
    mutate(RelapseSite_Timing = case_when(Relapse_Site == "No rec before death or new primary" ~ "No rec before death or new primary",
                                          Relapse_Site %in% c("Intrathoracic (only)", "Extrathoracic (any)") ~ paste(Relapse_Site, rec_1yr))) %>%
    mutate(RelapseSite_Timing_rec = factor(RelapseSite_Timing, levels = c("Intrathoracic (only) rec >= 1yr",
                                                                          "Intrathoracic (only) rec < 1yr",
                                                                          "Extrathoracic (any) rec >= 1yr",
                                                                          "Extrathoracic (any) rec < 1yr"))) 
  
  
  ## cochrane-armitage test for trend
  library(rstatix)
  
  table(in_data$SCNA.ITH, in_data$rec_1.2.3yr)
  # rec < 1y rec 1-2y rec 2-3y rec >= 3y
  # low        20       24       16        12
  # high       50       13        7         4
  
  prop_trend_test(table(in_data$SCNA.ITH, in_data$rec_1.2.3yr), score = NULL)  #default: group score = NULL
  #         n statistic          p p.signif    df method               
  #   1   146      19.3 0.0000111 ****         1 Chi-square trend test
  
  
  ###  barplot - stratify only by time 
  x_var_list <- c("SCNA.ITH") 
  y_var_list <- c("rec_1.2.3yr") 
  out_f <- paste0(workdir,"/barplot_SCNAITH_time.3yr.pdf")
  pdf(out_f, width = 6, height = 6)
  
  for(x_var in x_var_list){
    for(y_var in y_var_list){
      
      df <- in_data %>%
        mutate(rec_1.2.3yr = factor(rec_1.2.3yr , levels = rev(c("rec < 1y", "rec 1-2y","rec 2-3y", "rec >= 3y"))))
      df$x_var <- df[[x_var]]
      df$y_var <- as.factor(df[[y_var]])
      df <- df %>%
        filter(!is.na(x_var),!is.na(y_var))
      
      df_count <- df %>%
        group_by(y_var, x_var) %>%
        summarise(count = n()) %>%
        arrange(desc(y_var))
      
      g <-ggplot(data=df_count, aes(x=x_var, y=count, fill=y_var, label = count)) +
        geom_bar(stat="identity", position = "fill", color = "black") +
        geom_bar_text( position = "fill", reflow = TRUE) + 
        scale_fill_manual(values =rev(c("rec < 1y"= "#B15928" ,"rec 1-2y"= "#FF7F00",  "rec 2-3y"="#FDBF6F", "rec >= 3y"= "#FFFF99")) )+
        theme_cowplot() + 
        labs(x = x_var, y = "Proportion of patients (%)", fill = y_var)   + 
        theme(axis.title.x = element_text(size=24), 
              axis.title.y = element_text(size=24),
              axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
              axis.text.y = element_text(size=16), 
              legend.position = "right",
              legend.title = element_text(size=22), 
              legend.text = element_text(size=18),
              strip.text = element_text(size=16)) 
      print(g)
    }
  }
  dev.off()
  
  ### barplot - site & timing (1yr)
  x_var_list <- c("SCNA.ITH") 
  y_var_list <- c("RelapseSite_Timing_rec")
  out_f <- paste0(workdir,"/barplot_SCNAITH_site.time.pdf")
  pdf(out_f, width = 8, height = 6)
  
  for(x_var in x_var_list){
    for(y_var in y_var_list){
      
      df <- in_data
      df$x_var <- df[[x_var]]
      df$y_var <- as.factor(df[[y_var]])
      df <- df %>%
        filter(!is.na(x_var),!is.na(y_var))
      
      df_count <- df %>%
        group_by(y_var, x_var) %>%
        summarise(count = n())
      
      g <-ggplot(data=df_count, aes(x=x_var, y=count, fill=y_var, label = count)) +
        geom_bar(stat="identity", position = "fill", color = "black") +
        geom_bar_text( position = "fill", reflow = TRUE) +
        scale_fill_brewer(palette = "Paired") +
        theme_cowplot() + 
        labs(x = x_var, y = "Proportion of patients (%)", fill = y_var)   + 
        theme(axis.title.x = element_text(size=24), 
              axis.title.y = element_text(size=24),
              axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
              axis.text.y = element_text(size=16), 
              legend.position = "right",
              legend.title = element_text(size=22), #legend size
              legend.text = element_text(size=18),
              strip.text = element_text(size=16)) 
      print(g)
    }
  }
  dev.off()
  
  
  ###########################################
  ## linear regression for time-to-relapse ##
  ###########################################
  in_df <- in_data_surv 
  
  variables <- c("Age","pTNM_stage", "PackYears","histology_LUAD","adjuvant_tx", "SCNA.ITH", "subclonal_gd","RecentSubclExpansionScore_2" ) 
  
  in_df_nona1 <- in_df %>%
    dplyr::select(time.to.rec, Age, pTNM_stage , PackYears, adjuvant_tx, histology_LUAD,SCNA.ITH,RecentSubclExpansionScore_2, subclonal_gd  ) %>%
    na.omit()
  
  in_df_nona2 <- in_df %>%
    dplyr::select(Relapse.Site, Age, pTNM_stage , PackYears, adjuvant_tx, histology_LUAD,SCNA.ITH, RecentSubclExpansionScore_2, subclonal_gd  ) %>%
    na.omit()
  
  f1 <- as.formula(paste("time.to.rec",
                         paste(variables, collapse = " + "), 
                         sep = " ~ "))
  
  f2 <-  as.formula(paste("Relapse.Site",
                          paste(variables, collapse = " + "), 
                          sep = " ~ "))
  
  res_lm <- lm(f1 , data = in_df_nona1)
  f1
  summary(res_lm)
  
  res_glm <- glm(f2 , data = in_df_nona2, family = "binomial")
  f2
  summary(res_glm)
  
  
  ## forestplot 
  library(ggforestplot)
  
  res_list <- list(res_lm, res_glm)
  
  title_vec <- c("Time to relapse (linear regression)","Relapse site (logistic regression)")
  x_lab_list <- c("Coefficient (days)", "Odds ratio")
  logodds_list <- c(FALSE, TRUE)
  
  out_pdf <- paste0(workdir,"/forestplot_linear.and.logistic.pdf")
  pdf(out_pdf, width = 6, height = 5)
  
  for(i in 1:2){
    title_phrase <- title_vec[i]
    res <- res_list[[i]]
    x_lab <- x_lab_list[i]
    logodds <- logodds_list[[i]]
    
    df <- summary(res)$coef[-1,] %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      dplyr::rename(variables = rowname,
                    se = "Std. Error") 
    
    fp <- ggforestplot::forestplot(
      df = df,
      name = variables,
      estimate = Estimate,
      se = se,
      logodds = logodds,
    ) + 
      ggtitle(title_phrase) +
      labs(x = x_lab)
    print(fp)
    
    # make dataframe to save the coef, 95%CI and pvalues 
    
    if(i == 1){
      df_out <- summary(res)$coef[-1,] %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(variables = rowname,
                      se = "Std. Error")
      CI_table <- Confint(res)[,c(2,3)][-1,] %>%
        as.data.frame()
      df_out <- cbind(df_out,  CI_table, title_phrase)
      df_out1 <- df_out[,c(1,2,6,7,5,8)]
      colnames(df_out1) <- c("variables", "Estimate", "lower95CI", "upper95CI", "Pval", "analysis" )
    } else if (i == 2) {
      df_out <- summary(res)$coef[-1,]%>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(variables = rowname,
                      se = "Std. Error") %>%
        mutate(Estimate = exp(Estimate))
      CI_table <- exp(Confint(res)[,c(2,3)][-1,]) %>%
        as.data.frame()
      df_out <- cbind(df_out,  CI_table, title_phrase)
      df_out2 <- df_out[,c(1,2,6,7,5,8)]
      colnames(df_out2) <- c("variables", "Estimate", "lower95CI", "upper95CI", "Pval", "analysis" )
    }
  }
  dev.off()
  
  df_out_all <- rbind(df_out1, df_out2)
  
  out_f <- paste0(workdir, "/linear.and.logistic.csv")
  write.csv(df_out_all, out_f, row.names = F)
  
  #################################################
}





