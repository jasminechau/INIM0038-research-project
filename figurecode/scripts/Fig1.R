################################################################################################################################################
#############               Swimmers plot for Tx421 patients split into new lesion and no new lesion per cancer type               ############# 
################################################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk)

# Description:
# Script to create Figure 1 of the manuscript "The natural history of NSCLC in TRACERx"

#options and libraries
options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(ggnewscale)
library(fst)
library(gtable)

# setwd('/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/TRACERx_421/TRACERx421_data_code/20220808_Tumour_evoutionary_histories')

#parameters
output_dir  <- './'
crukid_size <- 9

#load patient data
patient_df <- readRDS('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_patient_df.rds')
chemo_df   <- readRDS('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_chemo_df.rds')

#load evo metrics 
metrics <- read.table('../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_evolutionary_metrics.tsv', header = T, sep = '\t')

#load mutation data
muttable_df <- read_fst('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst')



#############################
#####     Functions     #####
#############################

#Extract the legend from a ggplot object with this function
g_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) == 0){
    return(NULL)
  }
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Then use the following function on a list of legends to return a single grob that will contain all the legends aligned
align.legends.fun  <- function(legend.list){
  library(gtable)
  aligned.leg <- legend.list[[1]]
  for(i in 2:length(legend.list)){
    leg1        <- legend.list[[i]]$grobs[[1]]
    leg_a       <- gtable_add_rows(aligned.leg, pos = nrow(aligned.leg) - 1, heights = sum(leg1$heights))
    leg_final   <- gtable_add_grob(leg_a, leg1, t = nrow(leg_a) - 1, l = 3)
    aligned.leg <- leg_final
  }
  return(aligned.leg)
}


#function to plot swimmers with annotations as lower trinagle
lowerTriangle_swimmers_plot <- function(plot_data, plot_chemo_data, legend = T, max_time = 1000, label_size = 7){
  
  #prepare chemo data
  plot_chemo_data <- plot_chemo_data %>%
    left_join(plot_data[,c('patient_id', 'cruk_id')]) %>%
    mutate(agents = paste0('Platinum:', CHMOthDgName_cleaned)) %>%
    filter(!is.na(CHMPlatDgName_cleaned) | !is.na(CHMOthDgName_cleaned))
  
  plot_chemo_range <- plot_chemo_data %>%
    group_by(cruk_id) %>%
    summarise(start = min(CHMCycle_time), 
              stop = max(CHMCycle_time))
  
  #plot names based on relapse category
  plot_data$label_colour <- 'black'
  plot_data$label_colour[!is.na(plot_data$Relapse_cat)] <- '#737373'
  plot_data$label_colour[grep('primary', plot_data$Relapse_cat)] <- '#bdbdbd'
  
  #swimmers plot
  p_swimmer <- ggplot(plot_data, aes(x = os_time, y = cruk_id)) + 
    geom_segment(aes(x = 0, xend = os_time, y = cruk_id, yend = cruk_id), size = 0.1, linetype = 'dotted') +
    geom_segment(aes(x = AdjRadStartTime_manual, xend = AdjRadEndTime_manual, y = cruk_id, yend = cruk_id, colour = 'Radiotherapy'), size = 1.3) +
    scale_colour_manual(name = '', values = '#f781bf') +
    new_scale_colour() +
    geom_segment(data = plot_chemo_range, aes(x = start, xend = stop, y = cruk_id, yend = cruk_id), size = 0.5, colour = '#08306b') +
    geom_point(data = plot_chemo_data, aes(x = CHMCycle_time, y = cruk_id, colour = agents), size = 1.5) +
    scale_colour_manual(name = 'Adjuvant Therapy', values = c('Platinum:Etoposide' = '#a6cee3', 'Platinum:Pemetrexed' = '#8c6bb1', 'Platinum:Vinorelbine' = '#1f78b4')) +
    new_scale_colour() +
    geom_point(aes(x = Recurrence_time_use, y = cruk_id, colour = Relapse_cat_relapse, shape = Relapse_cat), shape = 18, size = 3) +
    geom_point(aes(x = newPrim_time_use, y = cruk_id, colour = Relapse_cat_newPrim, shape = Relapse_cat), shape = 8, size = 3) +
    scale_colour_manual(name = 'New Lesion', values = c('Intrathoracic' = '#fed976', 'Extrathoracic' = '#fd8d3c', 'Intra & Extra' = '#bd0026', 'Missing data' = '#969696', 'Second primary lung' = '#74c476', 'Other second primary' = '#006d2c'),
                        labels = c('Intrathoracic' = 'Intrathoracic Recurrence', 'Extrathoracic' = 'Extrathoracic Recurrence', 'Intra & Extra' = 'Intra- & Extrathoracic Recurrence', 'Missing data' = 'Unknown Recurrence', 'Second primary lung' = 'New Primary Lung', 'New Primary Other' = 'Second primary not lung')) +
    geom_point(data = plot_data[plot_data$cens_os %in% 0,], aes(x = os_time, y = cruk_id), size = 2, shape = 'x') + 
    geom_point(data = plot_data[plot_data$cens_os %in% 1,], aes(x = os_time, y = cruk_id), size = 2, shape = '>') + 
    ylab('') + xlab('Days post surgery') +
    ggtitle('CRUK ID') +
    scale_x_continuous(expand = c(0,0), limits = c(0,max_time), position = 'bottom') +
    theme_classic() + theme(axis.line.y = element_blank(), strip.background = element_rect(fill = '#d9d9d9'), 
                            axis.text.y = element_text(colour = plot_data$label_colour, size = label_size), panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), axis.title = element_blank(), 
                            plot.title.position = 'plot', plot.title = element_text(size = label_size + 1, vjust = -1)) 
  
  #stage
  p_stage <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = pathologyTNM)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = 'Stage', values = c('IA' = '#f5f5f5', 'IB' = '#c7eae5', 'IIA' = '#80cdc1', 'IIB' = '#35978f', 'IIIA' = '#01665e', 'IIIB' = '#003c30'), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Stage') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  
  
  #smoking status
  p_smoking <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = smoking_status_merged)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = 'Smoking Status', values = c("Never Smoked" = "#187FC3", "Ex-Smoker"  = "grey40","Smoker" =  "black"), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Smoking') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  
  #number of regions
  p_nRegions <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = nRegions)) + 
    geom_tile(size = 0) + 
    scale_fill_gradient(name = ' # Regions', low = '#e6f5d0', high = '#276419', limits=c(2,11), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = '#Regions') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())
  
  #Mutation ITH
  p_mutITH <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = subclonal_frac)) + 
    geom_tile(size = 0) + 
    scale_fill_gradient(name = 'SNV ITH', low = '#fde0ef', high = '#c51b7d', limits=c(0,1), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'SNV ITH') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  #Copy Number ITH
  p_cnITH <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = frac_abberant_genom_subcl)) + 
    geom_tile(size = 0) + 
    scale_fill_gradient(name = 'SCNA ITH', low = '#e7d4e8', high = '#3f007d', limits=c(0,1), na.value = "#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'SCNA ITH') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  
  #clonal WGS
  p_clonalWGD <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = num_clonal_gds)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = ' # Clonal GD', values = c('0' = '#d1e5f04D', '1' = '#92c5de', '2' = '#2166ac', 'not applicable' = '#bababa'), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Clonal GD') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())
  
  #subclonal WGS
  p_subclonalWGD <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = num_subclonal_gds)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = ' # Sublonal GD', values = c('0' = '#fddbc74D', '1' = '#f4a582', '2' = '#b2182b', 'not applicable' = '#bababa'), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Subclonal GD') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())
  
  
  #plot to extract legend for death and censored
  p_events <- ggplot(plot_data, aes(x = os_time, y = cruk_id)) + 
    geom_point(aes(shape = factor(cens_os)), size = 2) + 
    scale_shape_manual(name = 'Event', values = c('0' = 'x', '1' = '>'),
                       labels = c('0' = 'Death', '1' = 'Censored')) +
    theme_bw()
  
  
  #align legends
  legends         <- list(g_legend(p_swimmer), g_legend(p_events), g_legend(p_stage), g_legend(p_smoking), g_legend(p_nRegions), g_legend(p_mutITH), g_legend(p_cnITH), g_legend(p_clonalWGD), g_legend(p_subclonalWGD))
  aligend.legends <- align.legends.fun(legends)
  
  
  #exclude legends from original plots
  p_swimmer      <- p_swimmer + theme(legend.position = 'none')
  p_stage        <- p_stage + theme(legend.position = 'none')
  p_smoking      <- p_smoking + theme(legend.position = 'none')
  p_nRegions     <- p_nRegions + theme(legend.position = 'none')
  p_mutITH       <- p_mutITH + theme(legend.position = 'none')
  p_cnITH        <- p_cnITH + theme(legend.position = 'none')
  p_clonalWGD    <- p_clonalWGD + theme(legend.position = 'none')
  p_subclonalWGD <- p_subclonalWGD + theme(legend.position = 'none')
  
  #set margins
  p_swimmer      <- ggplotGrob(p_swimmer + theme(plot.margin = unit(c(0.5, 0, 0.5, 0.1), "cm")))
  p_stage        <- ggplotGrob(p_stage + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_smoking      <- ggplotGrob(p_smoking + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_nRegions     <- ggplotGrob(p_nRegions + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_mutITH       <- ggplotGrob(p_mutITH + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_cnITH        <- ggplotGrob(p_cnITH + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_clonalWGD    <- ggplotGrob(p_clonalWGD + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_subclonalWGD <- ggplotGrob(p_subclonalWGD + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  
  #align heights
  plot_list <- list(p_swimmer, p_stage, p_smoking, p_nRegions, p_mutITH, p_cnITH, p_clonalWGD, p_subclonalWGD)
  all_heights <- lapply(plot_list, function(x) {x$heights})
  plot_list_alignedHeights <- lapply(plot_list, function(x){
    x$heights <- do.call(unit.pmax, all_heights)
    return(x)
  })
  
  if(legend){
    return(list(plot = plot_list_alignedHeights, legend = aligend.legends))
  } else {
    return(plot_list_alignedHeights)
  }
  
  
}



#function to plot swimmers with annotations as uper trinagle
upperTriangle_swimmers_plot <- function(plot_data, plot_chemo_data, max_time = 1000, label_size = 5){
  
  #prepare chemo data
  plot_chemo_data <- plot_chemo_data %>%
    left_join(plot_data[,c('patient_id', 'cruk_id')]) %>%
    mutate(agents = paste0('Platinum:', CHMOthDgName_cleaned)) %>%
    filter(!is.na(CHMPlatDgName_cleaned) | !is.na(CHMOthDgName_cleaned))
  
  plot_chemo_range <- plot_chemo_data %>%
    group_by(cruk_id) %>%
    summarise(start = min(CHMCycle_time), 
              stop = max(CHMCycle_time))
  
  #plot names based on relapse category
  plot_data$label_colour <- 'black'
  plot_data$label_colour[!is.na(plot_data$Relapse_cat)] <- '#737373'
  plot_data$label_colour[grep('primary', plot_data$Relapse_cat)] <- '#bdbdbd'
  
  #swimmers plot
  p_swimmer <- ggplot(plot_data, aes(x = os_time, y = cruk_id)) + 
    geom_segment(aes(x = 0, xend = os_time, y = cruk_id, yend = cruk_id), size = 0.1, linetype = 'dotted') +
    geom_segment(aes(x = AdjRadStartTime_manual, xend = AdjRadEndTime_manual, y = cruk_id, yend = cruk_id, colour = 'Radiotherapy'), size = 1.3) +
    scale_colour_manual(name = '', values = '#f781bf') +
    new_scale_colour() +
    geom_segment(data = plot_chemo_range, aes(x = start, xend = stop, y = cruk_id, yend = cruk_id), size = 0.5, colour = '#08306b') +
    geom_point(data = plot_chemo_data, aes(x = CHMCycle_time, y = cruk_id, colour = agents), size = 1.5) +
    scale_colour_manual(name = 'Adjuvant Therapy', values = c('Platinum:Etoposide' = '#a6cee3', 'Platinum:Pemetrexed' = '#8c6bb1', 'Platinum:Vinorelbine' = '#1f78b4')) +
    new_scale_colour() +
    geom_point(aes(x = Recurrence_time_use, y = cruk_id, colour = Relapse_cat_relapse, shape = Relapse_cat), shape = 18, size = 3) +
    geom_point(aes(x = newPrim_time_use, y = cruk_id, colour = Relapse_cat_newPrim, shape = Relapse_cat), shape = 8, size = 3) +
    scale_colour_manual(name = 'New Lesion', values = c('Intrathoracic' = '#fed976', 'Extrathoracic' = '#fd8d3c', 'Intra & Extra' = '#bd0026', 'Missing data' = '#969696', 'Second primary lung' = '#74c476', 'Other second primary' = '#006d2c'),
                        labels = c('Intrathoracic' = 'Intrathoracic Recurrence', 'Extrathoracic' = 'Extrathoracic Recurrence', 'Intra & Extra' = 'Intra- & Extrathoracic Recurrence', 'Missing data' = 'Unknown Recurrence', 'Second primary lung' = 'New Primary Lung', 'New Primary Other' = 'Second primary not lung')) +
    geom_point(data = plot_data[plot_data$cens_os %in% 0,], aes(x = os_time, y = cruk_id), size = 2, shape = 'x') + 
    geom_point(data = plot_data[plot_data$cens_os %in% 1,], aes(x = os_time, y = cruk_id), size = 2, shape = '<') + 
    ylab('') + xlab('Days post surgery') +
    scale_x_reverse(expand = c(0,0), position = 'top', limits = c(max_time,0)) +
    scale_y_discrete(position = 'right') +
    theme_classic() + theme(axis.line.y = element_blank(), strip.background = element_rect(fill = '#d9d9d9'), 
                            axis.text.y = element_text(colour = plot_data$label_colour, size = label_size), panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), legend.position = 'none', axis.title.y = element_blank()) 
  
  #stage
  p_stage <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = pathologyTNM)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = 'Stage', values = c('IA' = '#f5f5f5', 'IB' = '#c7eae5', 'IIA' = '#80cdc1', 'IIB' = '#35978f', 'IIIA' = '#01665e', 'IIIB' = '#003c30'), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Stage') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), legend.position = 'none')
  
  
  #smoking status
  p_smoking <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = smoking_status_merged)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = 'Smoking Status', values = c("Never Smoked" = "#187FC3", "Ex-Smoker"  = "grey40","Smoker" =  "black"), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Smoking') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none', panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  
  #number of regions
  p_nRegions <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = nRegions)) + 
    geom_tile(size = 0) + 
    scale_fill_gradient(name = ' # Regions', low = '#e6f5d0', high = '#276419', limits=c(2,11), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = '#Regions') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank(), legend.position = 'none')
  
  #Mutation ITH
  p_mutITH <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = subclonal_frac)) + 
    geom_tile(size = 0) + 
    scale_fill_gradient(name = 'SNV ITH', low = '#fde0ef', high = '#c51b7d', limits=c(0,1), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'SNV ITH') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none', panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  
  #Copy Number ITH
  p_cnITH <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = frac_abberant_genom_subcl)) + 
    geom_tile(size = 0) + 
    scale_fill_gradient(name = 'SCNA ITH', low = '#e7d4e8', high = '#3f007d', limits=c(0,1), na.value = "#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'SCNA ITH') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none', panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA))
  
  #clonal WGS
  p_clonalWGD <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = num_clonal_gds)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = ' # Clonal GD', values = c('0' = '#d1e5f04D', '1' = '#92c5de', '2' = '#2166ac', 'not applicable' = '#bababa'), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Clonal GD') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank(), legend.position = 'none')
  
  #subclonal WGS
  p_subclonalWGD <- ggplot(plot_data, aes(x = 1, y = cruk_id, fill = num_subclonal_gds)) + 
    geom_tile(size = 0) + 
    scale_fill_manual(name = ' # Sublonal GD', values = c('0' = '#fddbc74D', '1' = '#f4a582', '2' = '#b2182b', 'not applicable' = '#bababa'), na.value="#bababa") +
    scale_x_continuous(expand = c(0,0), breaks = 1, labels = 'Subclonal GD') +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_rect(fill = "transparent",colour = NA),
                       plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank(), legend.position = 'none')
  
  #set margins
  p_swimmer      <- ggplotGrob(p_swimmer + theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0), "cm")))
  p_stage        <- ggplotGrob(p_stage + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_smoking      <- ggplotGrob(p_smoking + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_nRegions     <- ggplotGrob(p_nRegions + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_mutITH       <- ggplotGrob(p_mutITH + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_cnITH        <- ggplotGrob(p_cnITH + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_clonalWGD    <- ggplotGrob(p_clonalWGD + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  p_subclonalWGD <- ggplotGrob(p_subclonalWGD + theme(plot.margin = unit(c(0.5, 0, 0.5, -0.1), "cm")))
  
  #align heights
  plot_list <- list(p_swimmer, p_stage, p_smoking, p_nRegions, p_mutITH, p_cnITH, p_clonalWGD, p_subclonalWGD)
  all_heights <- lapply(plot_list, function(x) {x$heights})
  plot_list_alignedHeights <- lapply(plot_list, function(x){
    x$heights <- do.call(unit.pmax, all_heights)
    return(x)
  })
  
  return(plot_list_alignedHeights)
  
}




########################
#####     Main     #####
########################

patient_df$patient_id <- patient_df$cruk_id
plot_data_all <- plot_data <- patient_df

#set Recurrence_time_use for CRUK0249, CRUK0781 and CRUK0129 to time of death with unknown recurrence
plot_data_all$Recurrence_time_use[plot_data_all$cruk_id %in% c('CRUK0249', 'CRUK0781', 'CRUK0129')] <- plot_data_all$os_time[plot_data_all$cruk_id %in% c('CRUK0267', 'CRUK0781', 'CRUK0129')]

#filter chemo treatment timings
chemo_df <- chemo_df %>%
  filter(patient_id %in% plot_data_all$patient_id)
chemo_df <- chemo_df[!is.na(chemo_df$CHMCycle_time),]


#calculate fraction of subclonal mutation using most advanced tumour
#--> use all tumour for LTX799 and LTX829 because we don't know which one is most advanced
clonality_mutBurden <- lapply(plot_data_all$patient_id, function(x){
  print(x)
  tumour <- unlist(strsplit(plot_data_all$tumour_id_per_patient[plot_data_all$patient_id == x], ';'))
  counts <- table(muttable_df$combTiming_SC[muttable_df$tumour_id %in% tumour])
  
  #check NA's
  na_clonal <- na_subclonal_shared <- na_subclonal_private <- 0
  if ('NA' %in% names(counts)){
    na_counts <- table(muttable_df$ITHState[muttable_df$tumour_id %in% tumour & muttable_df$combTiming_SC == 'NA'])
    na_clonal <- ifelse('1' %in% names(na_counts), na_counts['1'], 0)
    na_subclonal_shared <- ifelse('2' %in% names(na_counts), na_counts['2'], 0)
    na_subclonal_private <- ifelse('3' %in% names(na_counts), na_counts['3'], 0)
  }
  
  
  #divide subclonal into private and shared
  subclonal_shared <- subclonal_private <- 0
  if('subclonal' %in% names(counts)){
    subclonal_muts    <- table(muttable_df$ITHState[muttable_df$tumour_id %in% tumour & muttable_df$combTiming_SC == 'subclonal'])
    subclonal_shared  <- sum(subclonal_muts[names(subclonal_muts) %in% 1:2]) 
    subclonal_private <- ifelse('3' %in% names(subclonal_muts), subclonal_muts['3'], 0)
  }
  
  #combine counts
  data.frame(patient_id = x,
             Clonal_early = ifelse('early' %in% names(counts), counts['early'], 0),
             Clonal_late = ifelse('late' %in% names(counts), counts['late'], 0),
             Clonal_untimed = ifelse('Unknown' %in% names(counts), counts['Unknown'], 0) + na_clonal,
             Subclonal_shared = subclonal_shared + na_subclonal_shared,
             Subclonal_private = subclonal_private + na_subclonal_private,
             total = sum(muttable_df$tumour_id %in% tumour))
})
clonality_mutBurden                <- Reduce(rbind, clonality_mutBurden)
clonality_mutBurden$clonal_frac    <- (clonality_mutBurden$Clonal_early + clonality_mutBurden$Clonal_late + clonality_mutBurden$Clonal_untimed) / clonality_mutBurden$total
clonality_mutBurden$subclonal_frac <- (clonality_mutBurden$Subclonal_shared + clonality_mutBurden$Subclonal_private) / clonality_mutBurden$total

plot_data_all <- plot_data_all %>% left_join(clonality_mutBurden[,c('patient_id', 'subclonal_frac')], by = 'patient_id')


#use copy number metrics for most advanced tumour 
#--> use all tumour for LTX799 and LTX829 because we don't know which one is most advanced
SCNA_metrics <- lapply(plot_data_all$patient_id, function(x){
  print(x)
  tumour <- unlist(strsplit(plot_data_all$tumour_id_per_patient[plot_data_all$patient_id == x], ';'))
  
  sub_metrics <- metrics %>%
    dplyr::select(tumour_id, perc_subclonal, frac_abberant_genom_subcl, num_clonal_gds, num_subclonal_gds) %>%
    filter(tumour_id %in% tumour)
  
  if(nrow(sub_metrics) > 1){
    out <- data.frame(patient_id = x, matrix(colMeans(sub_metrics[,-1], na.rm = T), nrow = 1))
  } else if(nrow(sub_metrics) == 0){
    out <- data.frame(patient_id = x, matrix(NA, nrow = 1, ncol = 4))
  } else {
    out <- data.frame(patient_id = x, sub_metrics[,-1])
  }
  colnames(out) <- c('patient_id', 'perc_subclonal', 'frac_abberant_genom_subcl', 'num_clonal_gds', 'num_subclonal_gds')
  
  return(out)
})
SCNA_metrics <- Reduce(rbind, SCNA_metrics) 

plot_data_all <- plot_data_all %>% 
  left_join(SCNA_metrics[,c('patient_id', 'frac_abberant_genom_subcl', 'num_clonal_gds', 'num_subclonal_gds')], by = 'patient_id') %>%
  mutate(num_subclonal_gds = ifelse(num_subclonal_gds > 2, 2, num_subclonal_gds))
plot_data_all <- plot_data_all %>%
  mutate(num_clonal_gds = ifelse(is.na(num_clonal_gds), 'not applicable', round(num_clonal_gds)),
         num_subclonal_gds = ifelse(is.na(num_subclonal_gds), 'not applicable', round(num_subclonal_gds)))


#count number of regions
plot_data_all$nRegions <- sapply(plot_data_all$patient_id, function(x){
  tumour   <- unlist(strsplit(plot_data_all$tumour_id_per_patient[plot_data_all$patient_id == x], ';'))
  nRegions <- sapply(tumour, function(i){ 
    sub      <- muttable_df[muttable_df$tumour_id %in% i, c('Is.present'), drop = F]
    length(unlist(strsplit(as.character(sub$Is.present[1]), ';')))
  })
  return(max(nRegions))
})

#for patient CRUK0682 the most advance tumour (LUAD) was not samples which is why tumour specific variables should be set to NA
plot_data_all$nRegions[plot_data_all$cruk_id %in% "CRUK0682"]                  <- NA
plot_data_all$subclonal_frac[plot_data_all$cruk_id %in% "CRUK0682"]            <- NA
plot_data_all$frac_abberant_genom_subcl[plot_data_all$cruk_id %in% "CRUK0682"] <- NA
plot_data_all$num_clonal_gds[plot_data_all$cruk_id %in% "CRUK0682"]            <- NA
plot_data_all$num_subclonal_gds[plot_data_all$cruk_id %in% "CRUK0682"]         <- NA

#create two seperate columns for relapse and new primaries to account for patients which have both
recurrence_and_newPrim <- plot_data_all %>%
  filter(!is.na(Recurrence_time_use) & !is.na(newPrim_time_use)) %>%
  filter(Recurrence_time_use != newPrim_time_use) %>%
  dplyr::select(cruk_id)
recurrence_and_newPrim <- as.character(recurrence_and_newPrim[,1])

plot_data_all$Relapse_cat_relapse <- plot_data_all$Relapse_cat_newPrim <- plot_data_all$Relapse_cat
plot_data_all$Relapse_cat_relapse[grep('primary', plot_data_all$Relapse_cat_relapse)] <- NA
plot_data_all$Relapse_cat_newPrim[grep('primary', plot_data_all$Relapse_cat_newPrim, invert = T)] <- NA
plot_data_all$Relapse_cat_relapse[plot_data_all$cruk_id %in% recurrence_and_newPrim] <- "Missing data"

#exckude CRUK from CRUK-id
plot_data_all$cruk_id <- sub('CRUK', '', plot_data_all$patient_id )


############
### LUAD 
############

#no new lesions
plot_data <- plot_data_all %>%
  filter(histology_lesion1_merged == 'Invasive adenocarcinoma') %>%
  filter(is.na(Relapse_cat)) %>%
  arrange(desc(os_time)) %>%
  mutate(cruk_id = factor(cruk_id, levels = cruk_id))

plot_chemo_data <- chemo_df[chemo_df$patient_id %in% plot_data$patient_id,]

plot_LUAD_noRelapse_list <- lowerTriangle_swimmers_plot(plot_data, plot_chemo_data, legend = T, max_time = max(plot_data_all$os_time) + 5, label_size = crukid_size)

#new lesions
plot_data <- plot_data_all %>%
  filter(histology_lesion1_merged == 'Invasive adenocarcinoma') %>%
  filter(!is.na(Relapse_cat)) %>%
  arrange(os_time) %>%
  mutate(cruk_id = factor(cruk_id, levels = cruk_id))

plot_chemo_data <- chemo_df[chemo_df$patient_id %in% plot_data$patient_id,]

plot_LUAD_Relapse_list <- upperTriangle_swimmers_plot(plot_data, plot_chemo_data, max_time = max(plot_data_all$os_time) + 5, label_size = crukid_size)



############
### LUSC 
############

#no new lesions
plot_data <- plot_data_all %>%
  filter(histology_lesion1_merged == 'Squamous cell carcinoma') %>%
  filter(is.na(Relapse_cat)) %>%
  arrange(desc(os_time)) %>%
  mutate(cruk_id = factor(cruk_id, levels = cruk_id))

plot_chemo_data <- chemo_df[chemo_df$patient_id %in% plot_data$patient_id,]

plot_LUSC_noRelapse_list <- lowerTriangle_swimmers_plot(plot_data, plot_chemo_data, legend = T, max_time = max(plot_data_all$os_time) + 5, label_size = crukid_size)

#new lesions
plot_data <- plot_data_all %>%
  filter(histology_lesion1_merged == 'Squamous cell carcinoma') %>%
  filter(!is.na(Relapse_cat)) %>%
  arrange(os_time) %>%
  mutate(cruk_id = factor(cruk_id, levels = cruk_id))

plot_chemo_data <- chemo_df[chemo_df$patient_id %in% plot_data$patient_id,]

plot_LUSC_Relapse_list <- upperTriangle_swimmers_plot(plot_data, plot_chemo_data, max_time = max(plot_data_all$os_time) + 5, label_size = crukid_size)


############
### Other 
############

#no new lesions
plot_data <- plot_data_all %>%
  filter(!histology_lesion1_merged %in% c('Squamous cell carcinoma', 'Invasive adenocarcinoma')) %>%
  filter(is.na(Relapse_cat)) %>%
  arrange(desc(os_time)) %>%
  mutate(cruk_id = factor(cruk_id, levels = cruk_id))

plot_chemo_data <- chemo_df[chemo_df$patient_id %in% plot_data$patient_id,]

plot_Other_noRelapse_list <- lowerTriangle_swimmers_plot(plot_data, plot_chemo_data, legend = T, max_time = max(plot_data_all$os_time) + 5, label_size = crukid_size)

#new lesions
plot_data <- plot_data_all %>%
  filter(!histology_lesion1_merged %in% c('Squamous cell carcinoma', 'Invasive adenocarcinoma')) %>%
  filter(!is.na(Relapse_cat)) %>%
  arrange(os_time) %>%
  mutate(cruk_id = factor(cruk_id, levels = cruk_id))

plot_chemo_data <- chemo_df[chemo_df$patient_id %in% plot_data$patient_id,]

plot_Other_Relapse_list <- upperTriangle_swimmers_plot(plot_data, plot_chemo_data, max_time = max(plot_data_all$os_time) + 5, label_size = crukid_size)


############
### Plots 
############

full_list <- c(rev(plot_LUAD_noRelapse_list$plot))
widths    <- c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.6) / 2
xpos      <- c(0, cumsum(widths)) + 0.01
full_plot <- ggdraw()
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0, width = widths[x], height = 0.95)
}

full_list <- c(plot_LUAD_Relapse_list)
widths    <- rev(c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.6)) / 2
xpos      <- c(0, cumsum(widths)) + 0.5-sum(widths) - 0.01
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0.38, width = widths[x], height = 0.62)
}


full_list <- c(rev(plot_LUSC_noRelapse_list$plot))
widths    <- c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.6) / 2
xpos      <- c(0, cumsum(widths)) + 0.51
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0.47, width = widths[x], height = 0.45)
}

full_list <- c(plot_LUSC_Relapse_list)
widths    <- rev(c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.6)) / 2
xpos      <- c(0, cumsum(widths)) + 1-sum(widths) - 0.01
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0.53, width = widths[x], height = 0.47)
}

full_list <- c(rev(plot_Other_noRelapse_list$plot))
widths    <- c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.6) / 2
xpos      <- c(0, cumsum(widths)) + 0.51
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0.23, width = widths[x], height = 0.2)
}

full_list <- c(plot_Other_Relapse_list)
widths    <- rev(c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.6)) / 2
xpos      <- c(0, cumsum(widths)) + 1-sum(widths) - 0.01
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = xpos[x], y = 0.25, width = widths[x], height = 0.22)
}

full_plot <- full_plot + draw_plot_label(c('LUAD', 'LUSC', 'Other'), x = c(0, 0.5, 0.5), y = c(1,1,0.46), size = 20)

pdf(paste0(output_dir, '/swimmersPlot_splitNewLesion.pdf'), width = 12, height = 18, useDingbats = F)
plot(full_plot)
dev.off()


#Legend
pdf(paste0(output_dir, '/Legend_swimmersPlot_splitNewLesion.pdf'), width = 5, height = 15, useDingbats = F)
plot(plot_LUAD_noRelapse_list$legend)
dev.off()

