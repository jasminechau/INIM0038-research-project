#####################################################################################################################
#############                         Genomic overview for Tx421 patients                               ############# 
#####################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk)

# Description:
# Script to create Extended figure 5 of the manuscript "The natural history of NSCLC in TRACERx"

#libraries and options
options(stringsAsFactors = F)
options(useDingbats = F)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(RColorBrewer)
library(grid)
library(fst)
library(gtable)

#setwd('/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/TRACERx_421/TRACERx421_data_code/20220808_Tumour_evoutionary_histories/')

#parameters
output_dir  <- './'
gab_size    <- 0
text_size   <- 12
legend_text_size  <- 10
legend_title_size <- 12

#load clinial data
clinicalData <- readRDS('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_patient_df.rds')

#load evo metrics
evo_metrics <- read.table('../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_evolutionary_metrics.tsv', header = T, sep = '\t')

#load mutation data
muttable_df <- read_fst('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst')

#load mutational signatures
signature_weights <- readRDS('../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_mutationSignature_weights.rds')
mutsigs_clonal    <- signature_weights$signature_weights_clonalMuts_perTumour
mutsigs_subclonal <- signature_weights$signature_weights_subclonalMuts_perTumour

#germline drivers
germline_df <- read.table('../20221109_Tumour_evo_histories_DATA/20220808_TRACERx421_GermlineDrivers.tsv', header = T, sep = '\t')

#load fusions
fusions <- read.table('../20221109_Tumour_evo_histories_DATA/20221116_TRACERx421_fusions.txt', header = T, sep = '\t')

#load subclonalExpansionScore
diversity_df <- read.csv('../20221109_Tumour_evo_histories_DATA/20221102_TRACERx421_tumourDivCorrectedtree.csv')

#################################### 
#####        Functions         ##### 
#################################### 

#Extract the legend from a ggplot object with this function:
g_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) == 0){
    return(NULL)
  }
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Then use the following function on a list of legends to return a single grob that will contain all the legends aligned:
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


###############################
#####        Main         ##### 
###############################

#Histologies --> use most advanced tumour for multiple tumour cases 
histology_data       <- clinicalData
histology_data$group <- 'Other'
histology_data$group[histology_data$histology_lesion1_merged %in% 'Invasive adenocarcinoma'] <- 'LUAD'
histology_data$group[histology_data$histology_lesion1_merged %in% 'Squamous cell carcinoma'] <- 'LUSC'
histology_data$group <- factor(histology_data$group, levels = c('LUAD', 'LUSC', 'Other'))

histology_data$plot_hist <- as.character(histology_data$histology_lesion1_merged)
histology_data$plot_hist[histology_data$plot_hist == 'Invasive adenocarcinoma']   <- paste0('LUAD ', ifelse(is.na(histology_data$LUAD_pred_subtype[histology_data$plot_hist == 'Invasive adenocarcinoma']), 'other', histology_data$LUAD_pred_subtype[histology_data$plot_hist == 'Invasive adenocarcinoma']))
histology_data$plot_hist[histology_data$plot_hist == 'Squamous cell carcinoma']   <- 'LUSC'
histology_data$plot_hist[histology_data$plot_hist == '(Invasive adenocarcinoma)'] <- 'LUAD'
histology_data$plot_hist <- factor(histology_data$plot_hist, levels = c('LUAD lepidic', 'LUAD papillary', 'LUAD acinar', 'LUAD cribriform', 'LUAD micropapillary', 'LUAD solid', 'LUAD invasive_mucinous', 'LUAD other',
                                                                        'LUSC','LCNEC', "Adenosquamous carcinoma", "Pleomorphic carcinoma", "Large cell carcinoma", 'Other'))

histology_data$facet_group <- paste(histology_data$group, histology_data$pathologyTNM, sep = ':')
histology_data$facet_group[histology_data$facet_group %in% c('LUAD:IIIA', 'LUAD:IIIB')] <- 'LUAD:IIIA+B'
histology_data$facet_group <- factor(histology_data$facet_group, levels = sort(unique(histology_data$facet_group)))

histology_data$patient_id <- histology_data$cruk_id

#order based on mutation heterogeneity
clonality_mutBurden <- lapply(histology_data$patient_id, function(x){
  print(x)
  tumour <- unlist(strsplit(histology_data$tumour_id_per_patient[histology_data$patient_id == x], ';'))
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
  
  #if multiple tumours (LTX799 & LTX829) --> use median total mutation count
  if(length(tumour) > 1){
    total <- median(table(muttable_df$tumour_id[muttable_df$tumour_id %in% tumour]))
  } else {
    total <- sum(muttable_df$tumour_id %in% tumour)
  }
  
  #combine counts
  data.frame(SampleID = x,
             Clonal_early = ifelse('early' %in% names(counts), counts['early'], 0),
             Clonal_late = ifelse('late' %in% names(counts), counts['late'], 0),
             Clonal_untimed = ifelse('Unknown' %in% names(counts), counts['Unknown'], 0) + na_clonal,
             Subclonal_shared = subclonal_shared + na_subclonal_shared,
             Subclonal_private = subclonal_private + na_subclonal_private,
             nMuts = total,
             total = sum(muttable_df$tumour_id %in% tumour))
})
clonality_mutBurden <- Reduce(rbind, clonality_mutBurden)
clonality_mutBurden_freq <- clonality_mutBurden
clonality_mutBurden_freq[,2:6] <- clonality_mutBurden_freq[,2:6] / clonality_mutBurden_freq$total
clonality_mutBurden_freq  <- clonality_mutBurden_freq[order((clonality_mutBurden_freq$Subclonal_private + clonality_mutBurden_freq$Subclonal_shared),
                                                            clonality_mutBurden_freq$Subclonal_private),]
order_samples <- clonality_mutBurden_freq$SampleID
clonality_mutBurden_freq$SampleID <- factor(clonality_mutBurden_freq$SampleID, levels = order_samples)


#plot nRegions
histology_data$nRegions <- sapply(histology_data$patient_id, function(x){
  tumour   <- unlist(strsplit(histology_data$tumour_id_per_patient[histology_data$patient_id == x], ';'))
  nRegions <- sapply(tumour, function(i){ 
    sub <- muttable_df[muttable_df$tumour_id %in% i, c('Is.present'), drop = F]
    length(unlist(strsplit(as.character(sub$Is.present[1]), ';')))
  })
  return(max(nRegions))
})
histology_data$patient_id <- factor(histology_data$patient_id, levels = order_samples)

p_nRegions <- ggplot(histology_data, aes(x = patient_id, y = nRegions, fill = plot_hist)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(name = 'Histology', values = c(brewer.pal(n = 9, 'OrRd')[-1], '#053061', rev(brewer.pal(n = 7, 'Greys')[-1]))) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,4,8), labels = c('0','4','8')) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('# Regions') +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

#Smoking
smoking_data <- histology_data %>%
  dplyr::select(patient_id, smoking_status_merged, facet_group) %>%
  mutate(patient_id = factor(patient_id, levels = order_samples))

p_smkstatus <- ggplot(smoking_data, aes(x = patient_id, y = 1, fill = smoking_status_merged)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_manual(name = 'Smoking Status', values = c('Never Smoked' = '#187fc3', 'Ex-Smoker'  = 'grey40','Smoker' =  "black")) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('Smoking Status') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), 
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())


#Treatment
treatment_data <- histology_data %>%
  dplyr::select(patient_id, adjuvant_treatment_given, facet_group) %>%
  mutate(patient_id = factor(patient_id, levels = order_samples),
         adjuvant_treatment_given = ifelse(is.na(adjuvant_treatment_given), 'none', adjuvant_treatment_given),
         adjuvant_treatment_given= factor(adjuvant_treatment_given, levels = c('none', 'Platinum chemo', 'Radiotherapy', 'Platinum chemo/radiotherapy'))) 

p_treatment <- ggplot(treatment_data, aes(x = patient_id, y = 1, fill = adjuvant_treatment_given)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_manual(name = 'Adjuvant Therapy', values = c('#f7fcf0', '#c7e9b4', '#41b6c4', '#225ea8'),
                    labels = c('No Treatment', 'Platinum Chemotherapy', 'Radiotherapy', 'Both')) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('Adjuvant Therapy') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), 
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())


#Recurrence
recurrence_data <- histology_data %>%
  dplyr::select(patient_id, Relapse_cat, facet_group) %>%
  mutate(patient_id = factor(patient_id, levels = order_samples),
         Relapse_cat = ifelse(is.na(Relapse_cat), 'no relapse', Relapse_cat))

p_recurrence <- ggplot(recurrence_data, aes(x = patient_id, y = 1, fill = Relapse_cat)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_manual(name = 'New Lesions', values = c('Extrathoracic' = '#fed976', 'Intra & Extra' = '#fd8d3c', 'Intrathoracic' = '#bd0026', 'Second primary lung' = '#74c476', 'Other second primary' = '#006d2c', 'Missing data' = '#bababa', 'no relapse' = '#f0f0f0'),
                    labels = c('Extrathoracic' = 'Extrathoracic Recurrence', 'Intra & Extra' = 'Intra- & Extrathoracic Recurrence', 'Intrathoracic' = 'Intrathoracic Recurrence', 'Second primary lung' = 'New Primary Lung', 'New Primary Other' = 'Second primary not lung', 'Missing data' = 'Unknown', 'no relapse' = 'No new lesion'),
                    na.value = 'white') +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('New Lesions') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), 
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())



### Genomic variables ###
#subclonal expansion score
diversity_data <- histology_data %>%
  dplyr::select(patient_id, facet_group) %>%
  mutate(patient_id = factor(patient_id, levels = order_samples)) %>%
  left_join(diversity_df, by = c('patient_id'))

p_diversity <- ggplot(diversity_data, aes(x = patient_id, y = 1, fill = maxCCF_terminalonly_tum_meanccf)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_gradientn(name = 'Subclonal Expansion Score', colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                       na.value = 'white') +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('Subclonal Expansion') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), 
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())

#mutation clonality
plot_freq_mutBurden <- clonality_mutBurden_freq %>% 
  left_join(histology_data[,c('patient_id', 'facet_group')], by = c('SampleID' = 'patient_id')) 
plot_freq_mutBurden <- reshape2::melt(plot_freq_mutBurden, id.vars = c('SampleID', 'facet_group', 'total', 'nMuts'))
plot_freq_mutBurden$SampleID <- factor(plot_freq_mutBurden$SampleID, levels = order_samples)
plot_freq_mutBurden$variable <- factor(plot_freq_mutBurden$variable, levels = c('Clonal_early', 'Clonal_untimed', 'Clonal_late', 'Subclonal_shared', 'Subclonal_private'))

p_mutBurden_freq <- ggplot(plot_freq_mutBurden, aes(x = SampleID, y = value * 100, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = F)) + 
  scale_fill_manual(name = 'Mutation Timing', values = rev(c('#b2182b', '#f4a582', '#92c5de', '#4393c3', '#2166ac')), drop = F) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100), labels = c('0','50','100')) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('% Mutations') + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

#plot frequency of cn clonality
SCNA_metrics <- lapply(as.character(histology_data$patient_id), function(x){
  print(x)
  tumour <- unlist(strsplit(histology_data$tumour_id_per_patient[histology_data$patient_id == x], ';'))
  
  sub_metrics <- evo_metrics %>%
    dplyr::select(tumour_id, perc_subclonal, frac_abberant_genom_subcl, num_clonal_gds, num_subclonal_gds) %>%
    filter(tumour_id %in% tumour)
  
  if(nrow(sub_metrics) > 1){
    out <- data.frame(SampleID = x, matrix(colMeans(sub_metrics[,-1], na.rm = T), nrow = 1))
  } else if(nrow(sub_metrics) == 0){
    out <- data.frame(SampleID = x, matrix(NA, nrow = 1, ncol = 4))
  } else {
    out <- data.frame(SampleID = x, sub_metrics[,-1])
  }
  colnames(out) <- c('SampleID', 'perc_subclonal', 'frac_abberant_genom_subcl', 'num_clonal_gds', 'num_subclonal_gds')
  
  return(out)
})
SCNA_metrics <- Reduce(rbind, SCNA_metrics) 
SCNA_metrics <- SCNA_metrics %>% left_join(histology_data[,c('patient_id', 'facet_group')], by = c('SampleID' = 'patient_id'))

plot_data <- SCNA_metrics[,c('SampleID', 'facet_group', 'frac_abberant_genom_subcl')]
plot_data$clonal <- 1 - plot_data$frac_abberant_genom_subcl
colnames(plot_data) <- c('sample_id', 'facet_group', 'subclonal', 'clonal')
plot_data <- reshape2::melt(plot_data, id.vars = c('sample_id', 'facet_group'))
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable[is.na(plot_data$value)] <- 'not applicable'
plot_data$value[is.na(plot_data$value)]    <- 0.5
plot_data$variable  <- factor(plot_data$variable, levels = c('subclonal', 'clonal', 'not applicable'))
plot_data$sample_id <- factor(plot_data$sample_id, levels = order_samples)

p_cnClonality <- ggplot(plot_data, aes(x = sample_id, y = value * 100, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + 
  scale_fill_manual(name = 'CN Clonality', values = c('clonal' = '#4393c3', 'subclonal' =  '#b2182b', 'not applicable' = '#bababa'), drop = F) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100), labels = c('0','50','100')) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('% CN') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))


#number of WGD events
plot_data <- SCNA_metrics[,c('SampleID', 'facet_group', 'num_clonal_gds', 'num_subclonal_gds')]
plot_data$num_subclonal_gds <- ifelse(plot_data$num_subclonal_gds > 2, 2, plot_data$num_subclonal_gds)
colnames(plot_data) <- c('sample_id', 'facet_group', 'clonal', 'subclonal')
plot_data <- reshape2::melt(plot_data, id.vars = c('sample_id', 'facet_group'))
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable[is.na(plot_data$value)] <- 'not applicable'
plot_data$value[is.na(plot_data$value)]    <- 1.5
plot_data$variable <- factor(plot_data$variable, levels = c('subclonal', 'clonal', 'not applicable'))
plot_data$sample_id <- factor(plot_data$sample_id, levels = order_samples)

p_WGD <- ggplot(plot_data, aes(x = sample_id, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = F)) + 
  scale_fill_manual(name = 'WGD', values = c('clonal' = '#4393c3', 'subclonal' =  '#b2182b', 'not applicable' = '#bababa'), drop = F) +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('# WGD') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))


# mutation burden #
plot_data <- clonality_mutBurden
plot_data <- plot_data %>% left_join(histology_data[, c('patient_id', 'facet_group')], by = c('SampleID' = 'patient_id')) 
plot_data$SampleID <- factor(plot_data$SampleID, levels = order_samples)
plot_data$nMuts[plot_data$nMuts > 5000] <- 5000

p_mutBurden <- ggplot(plot_data, aes(x = SampleID, y = nMuts)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(expand = c(0,0), breaks = seq(0,5000, 1000), labels = c(seq(0,4000, 1000), '>=5000')) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  ylab('# Mutations') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))


# mutational signatures #
#clonal
plot_mutsigs_clonal_data <- lapply(as.character(histology_data$patient_id), function(x){
  print(x)
  tumour <- unlist(strsplit(histology_data$tumour_id_per_patient[histology_data$patient_id == x], ';'))
  
  sub <- mutsigs_clonal[mutsigs_clonal$tumour_id %in% tumour,,drop = F]
  df  <- data.frame(sample_id = x, matrix(apply(sub[,grep('SBS', colnames(sub))], 2, function(i) mean(i)), nrow = 1))
  colnames(df) <- c('sample_id', grep('SBS', colnames(sub), value = T))
  return(df)
})
plot_mutsigs_clonal_data <- Reduce(rbind, plot_mutsigs_clonal_data)
plot_mutsigs_clonal_data$Unknown   <- 1 - rowSums(plot_mutsigs_clonal_data[,grep('SBS', colnames(plot_mutsigs_clonal_data))])
plot_mutsigs_clonal_data           <- plot_mutsigs_clonal_data %>%
  left_join(histology_data[, c('patient_id', 'facet_group')], by = c('sample_id' = 'patient_id')) %>%
  data.frame()
plot_mutsigs_clonal_data           <- reshape2::melt(plot_mutsigs_clonal_data, id.vars = c('sample_id', 'facet_group'))
plot_mutsigs_clonal_data$variable  <- factor(plot_mutsigs_clonal_data$variable, levels = c('SBS1', 'SBS5', 'SBS2', 'SBS13', 'SBS4', 'SBS92', 'SBS17b', 'SBS44', 'Unknown'))
plot_mutsigs_clonal_data$sample_id <- factor(plot_mutsigs_clonal_data$sample_id, levels = order_samples)

p_mutsigs_clonal <- ggplot(plot_mutsigs_clonal_data, aes(x = sample_id, y = value * 100, fill = variable)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  scale_fill_manual(name = 'Signatures', values = c(rev(brewer.pal(n = length(levels(plot_mutsigs_clonal_data$variable)) - 1, 'Paired')), '#bababa'),
                    labels = c('SBS1 (Aging)', 'SBS5', 'SBS2 (APOBEC)', 'SBS13 (APOBEC)', 'SBS4 (Smoking)', 'SBS92 (Smoking)', 'SBS17b', 'SBS44 (MSI)', 'unknown'),
                    drop = F) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100), labels = c('0', '50', '100')) +  
  ylab('% Clonal\nSignatures') +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), legend.position = 'none')

#subclonal
plot_mutsigs_subclonal_data <- lapply(as.character(histology_data$patient_id), function(x){
  print(x)
  tumour <- unlist(strsplit(histology_data$tumour_id_per_patient[histology_data$patient_id == x], ';'))
  
  sub <- mutsigs_subclonal[mutsigs_subclonal$tumour_id %in% tumour,,drop = F]
  df  <- data.frame(sample_id = x, matrix(apply(sub[,grep('SBS', colnames(sub))], 2, function(i) mean(i)), nrow = 1))
  colnames(df) <- c('sample_id', grep('SBS', colnames(sub), value = T))
  return(df)
})
plot_mutsigs_subclonal_data <- Reduce(rbind, plot_mutsigs_subclonal_data)
plot_mutsigs_subclonal_data$Unknown   <- 1 - rowSums(plot_mutsigs_subclonal_data[,grep('SBS', colnames(plot_mutsigs_subclonal_data))])
plot_mutsigs_subclonal_data           <- plot_mutsigs_subclonal_data %>%
  left_join(histology_data[, c('patient_id', 'facet_group')], by = c('sample_id' = 'patient_id')) %>%
  data.frame()
plot_mutsigs_subclonal_data           <- reshape2::melt(plot_mutsigs_subclonal_data, id.vars = c('sample_id', 'facet_group'))
plot_mutsigs_subclonal_data$variable  <- factor(plot_mutsigs_subclonal_data$variable, levels = c('SBS1', 'SBS5', 'SBS2', 'SBS13', 'SBS4', 'SBS92', 'SBS17b', 'SBS44', 'Unknown'))
plot_mutsigs_subclonal_data$sample_id <- factor(plot_mutsigs_subclonal_data$sample_id, levels = order_samples)

p_mutsigs_subclonal <- ggplot(plot_mutsigs_subclonal_data, aes(x = sample_id, y = value * 100, fill = variable)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) +
  facet_grid(.~facet_group, scales = 'free_x', space = 'free_x') +
  scale_fill_manual(name = 'Signatures', values = c(rev(brewer.pal(n = length(levels(plot_mutsigs_clonal_data$variable)) - 1, 'Paired')), '#bababa'),
                    labels = c('SBS1 (Aging)', 'SBS5', 'SBS2 (APOBEC)', 'SBS13 (APOBEC)', 'SBS4 (Smoking)', 'SBS92 (Smoking)', 'SBS17b', 'SBS44 (MSI)', 'unknown'),
                    drop = F) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100), labels = c('0', '50', '100')) +  
  ylab('% Subclonal\nSignatures') +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))


# count clonal and subclonal drivers #
driver_muttable <- muttable_df[muttable_df$DriverMut,]
count_drivers <- lapply(as.character(histology_data$patient_id), function(x){
  tumour <- unlist(strsplit(histology_data$tumour_id_per_patient[histology_data$patient_id == x], ';'))
  
  if(length(tumour) > 1){
    n_clonal <- sapply(tumour, function(i){
      sub    <- driver_muttable[driver_muttable$tumour_id %in% i,]
      clonal <- sum(sub$combTiming_SC %in% c('early', 'late', 'Unknown'))
      clonal <- clonal + sum(sub$combTiming_SC == 'NA' & sub$ITHState == 1)
      return(clonal)
    })
    tumour <- names(n_clonal)[n_clonal == max(n_clonal)]
  }
  
  sub    <- driver_muttable[driver_muttable$tumour_id %in% tumour,]
  clonal <- sum(sub$combTiming_SC %in% c('early', 'late', 'Unknown'))
  clonal <- clonal + sum(sub$combTiming_SC == 'NA' & sub$ITHState == 1)
  data.frame(SampleID = x, tumour = tumour, total = nrow(sub), clonal = clonal, subclonal = nrow(sub) - clonal)
})
count_drivers <- Reduce(rbind, count_drivers)

plot_driver <- count_drivers
plot_driver <- plot_driver %>% left_join(histology_data[, c('patient_id', 'facet_group')], by = c('SampleID' = 'patient_id')) 
plot_driver$SampleID <- factor(plot_driver$SampleID, levels = order_samples)

# plot_driver$clonal   <- count_drivers$clonal[match(plot_driver$tumour, count_drivers$tumour)]
plot_driver$clonal[plot_driver$clonal > 10] <- 10
# plot_driver$subclonal <- count_drivers$subclonal[match(plot_driver$tumour, count_drivers$tumour)]
plot_driver$subclonal[plot_driver$subclonal > 10] <- 10

p_driver_clonal <-  ggplot(plot_driver, aes(x = SampleID, y = clonal)) + 
  geom_bar(stat = 'identity', fill = '#4393c3') + 
  facet_grid(. ~ facet_group, scales = 'free_x', space = 'free_x', drop = FALSE) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,10,2), labels = c(seq(0,8,2), '>=10')) +  
  ylab('# Clonal\nDrivers') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

p_driver_subclonal <-  ggplot(plot_driver, aes(x = SampleID, y = subclonal)) + 
  geom_bar(stat = 'identity', fill = '#b2182b') + 
  facet_grid(. ~ facet_group, scales = 'free_x', space = 'free_x', drop = FALSE) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,10,2), labels = c(seq(0,8,2), '>=10')) +
  ylab('# Subclonal\nDrivers') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))


# driver heatmap of top 15 most ferquently mutated driver genes in cohort #
# top_genes <- c(names(sort(table(driver_muttable$Hugo_Symbol), decreasing = T)[1:15]), 'EGFR')
top_genes <-c('KRAS', 'TP53',   'KEAP1',  'NFE2L2', 'STK11',  'B2M', 'PTEN',   'SMAD4', 'EGFR')
driver_muttable <- driver_muttable[driver_muttable$tumour_id %in% count_drivers$tumour,]
plot_driver_heatmap <- lapply(top_genes, function(x){
  sub <- driver_muttable[driver_muttable$Hugo_Symbol == x,]
  sub$newTiming <- sub('Unknown' , 'untimed', sub$combTiming_SC)
  
  #check NA's
  if (sum(sub$newTiming == 'NA') > 0){
    sub$newTiming[sub$newTiming == 'NA' & sub$ITHState == 1] <- 'untimed'
    sub$newTiming[sub$newTiming == 'NA' & sub$ITHState == 2] <- 'Subclonal_shared'
    sub$newTiming[sub$newTiming == 'NA' & sub$ITHState == 3] <- 'Subclonal_private'
  }
  
  #divide subclonal into private and shared
  if(sum(sub$newTiming == 'subclonal') > 0){
    sub$newTiming[sub$newTiming == 'subclonal' & sub$ITHState %in% c(1,2)] <- 'Subclonal_shared'
    sub$newTiming[sub$newTiming == 'subclonal' & sub$ITHState == 3] <- 'Subclonal_private'
  }
  
  sub$newTiming[grep('Subclonal', sub$newTiming, invert = T)] <- paste('Clonal', sub$newTiming[grep('Subclonal', sub$newTiming, invert = T)], sep = '_')
  
  #mark duplicated hits
  duplicated_hits <- sub$tumour_id[duplicated(sub$tumour_id)]
  sub$newTiming[sub$tumour_id %in% unique(duplicated_hits)] <- 'Multiple_hits'
  sub <- sub[!duplicated(sub$tumour_id),]
  
  data.frame(gene = x, SampleID = sub$patient_id, tumour = sub$tumour_id, timing = sub$newTiming)
  
})
plot_driver_heatmap <- Reduce(rbind, plot_driver_heatmap)

add_patients <- order_samples[!order_samples %in% unique(plot_driver_heatmap$SampleID)]
plot_driver_heatmap <- rbind(plot_driver_heatmap, 
                             data.frame(gene = 'TP53', SampleID = add_patients, tumour = add_patients, timing = NA))
plot_driver_heatmap <- plot_driver_heatmap %>%  left_join(histology_data[, c('patient_id', 'facet_group')], by = c('SampleID' = 'patient_id')) 
plot_driver_heatmap$SampleID  <- factor(plot_driver_heatmap$SampleID, levels = order_samples)
plot_driver_heatmap$timing    <- factor(plot_driver_heatmap$timing, levels = c('Clonal_early', 'Clonal_untimed', 'Clonal_late', 'Subclonal_shared', 'Subclonal_private', 'Multiple_hits'))
plot_driver_heatmap$gene      <- factor(plot_driver_heatmap$gene, levels = rev(top_genes))

p_driver_heatmap <- ggplot(plot_driver_heatmap, aes(x = SampleID, y = gene, fill = timing)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_manual(name = 'Driver Mutation', values = rev(c('#999999', '#b2182b', '#f4a582', '#92c5de', '#4393c3', '#2166ac')), drop = F, na.translate = F) +
  geom_hline(yintercept = (1:length(top_genes)) - 0.5, size = 0.3) +
  scale_y_discrete(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ facet_group, scales = 'free_x', space = 'free_x', drop = FALSE) +
  theme_bw() + 
  ylab('Driver\nMutations') +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5), axis.text.y = element_text(face = "italic"),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))


# germline drivers #
plot_germline                 <- histology_data
plot_germline$GLdriver_binary <- 'NO'
plot_germline$GLdriver_binary[plot_germline$patient_id %in% germline_df$Patient] <- 'YES'
plot_germline$patient_id  <- factor(plot_germline$patient_id, levels = order_samples)

p_GLdriver <- ggplot(plot_germline, aes(x = patient_id, y = 1, fill = GLdriver_binary)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_manual(name = 'Germline Driver', values = c('NO' = '#fff7f3', 'YES' = '#ae017e'), drop = F) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ facet_group, scales = 'free_x', space = 'free_x', drop = FALSE) +
  ylab('Germline Drivers') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), 
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())


# fusion events #
plot_fusions          <- histology_data
fusions               <- fusions[!duplicated(paste(fusions$patient_id, fusions$Genes, sep = ':')),]
plot_fusions          <- plot_fusions %>% left_join(fusions[, c('patient_id', 'Genes')], by = 'patient_id')
plot_fusions$patient_id  <- factor(plot_fusions$patient_id, levels = order_samples)

p_fusions <- ggplot(plot_fusions, aes(x = patient_id, y = 1, fill = Genes)) + 
  geom_tile(colour = 'white', size = gab_size) +
  scale_fill_manual(name = 'Fusions', values = brewer.pal(n = 6, 'Accent'), drop = F, na.translate=FALSE) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ facet_group, scales = 'free_x', space = 'free_x', drop = FALSE) +
  ylab('Fusions') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = text_size), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), 
        strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing = unit(0.2, "lines"), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), panel.grid = element_blank())



##########################################
###          combine plots             ###
##########################################

#align legends
legends         <- list(g_legend(p_nRegions), g_legend(p_smkstatus), g_legend(p_treatment), g_legend(p_recurrence), g_legend(p_diversity),
                        g_legend(p_mutBurden_freq), g_legend(p_cnClonality), g_legend(p_WGD), g_legend(p_mutsigs_subclonal), 
                        g_legend(p_driver_heatmap), g_legend(p_GLdriver), g_legend(p_fusions))
aligend.legends <- align.legends.fun(legends)

p_nRegions          <- p_nRegions + theme(legend.position = 'none')
p_smkstatus         <- p_smkstatus + theme(legend.position = 'none')
p_treatment         <- p_treatment + theme(legend.position = 'none')
p_recurrence        <- p_recurrence + theme(legend.position = 'none')
p_diversity         <- p_diversity + theme(legend.position = 'none')
p_mutBurden_freq    <- p_mutBurden_freq+ theme(legend.position = 'none')
p_cnClonality       <- p_cnClonality + theme(legend.position = 'none')
p_WGD               <- p_WGD + theme(legend.position = 'none')
p_mutsigs_subclonal <- p_mutsigs_subclonal + theme(legend.position = 'none')
p_driver_heatmap    <- p_driver_heatmap + theme(legend.position = 'none')
p_GLdriver          <- p_GLdriver + theme(legend.position = 'none')
p_fusions           <- p_fusions + theme(legend.position = 'none')


#set margins
p_nRegions       <- p_nRegions + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
p_smkstatus      <- ggplotGrob(p_smkstatus + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_treatment      <- ggplotGrob(p_treatment + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_recurrence     <- ggplotGrob(p_recurrence + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_diversity      <- ggplotGrob(p_diversity + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_mutBurden_freq <- ggplotGrob(p_mutBurden_freq + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_cnClonality    <- ggplotGrob(p_cnClonality + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_WGD            <- ggplotGrob(p_WGD + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_mutBurden      <- ggplotGrob(p_mutBurden + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_mutsigs_clonal    <- ggplotGrob(p_mutsigs_clonal + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_mutsigs_subclonal <- ggplotGrob(p_mutsigs_subclonal + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_driver_clonal     <- ggplotGrob(p_driver_clonal + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_driver_subclonal  <- ggplotGrob(p_driver_subclonal + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_driver_heatmap    <- ggplotGrob(p_driver_heatmap + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_GLdriver          <- ggplotGrob(p_GLdriver + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))
p_fusions           <- ggplotGrob(p_fusions + theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")))


#change facet colours and labels
colours            <- matrix(unlist(strsplit(levels(histology_data$facet_group), ':')), ncol = 2, byrow = T)[,1]
colours[colours == 'LUAD']  <- '#67001f'
colours[colours == 'LUSC']  <- '#053061'
colours[colours == 'Other'] <- '#4d4d4d'
strip_name_colours <- rep('white', length(levels(histology_data$facet_group)))
strip_name         <- matrix(unlist(strsplit(levels(histology_data$facet_group), ':')), ncol = 2, byrow = T)[,2]

g       <- ggplot_gtable(ggplot_build(p_nRegions))
striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
  
  t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
  g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
  
  k <- k+1
}

#align widths
plot_list <- list(g, p_smkstatus, p_treatment, p_recurrence, p_diversity, p_mutBurden_freq, p_cnClonality, p_WGD, p_mutBurden, p_mutsigs_clonal, p_mutsigs_subclonal,
                  p_driver_clonal, p_driver_subclonal, p_driver_heatmap, p_GLdriver, p_fusions)
all_widths <- lapply(plot_list, function(x) {x$widths})
plot_list_alignedWidths <- lapply(plot_list, function(x){
  x$widths <- do.call(unit.pmax, all_widths)
  return(x)
})


#plot
full_list <- rev(plot_list_alignedWidths)
heights    <- c(0.02, 0.02, 0.17, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.02, 0.02, 0.02, 0.02, 0.1)
ypos      <- c(0, cumsum(heights)) + 0.03
full_plot <- ggdraw()
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = 0, y = ypos[x], width = 1, height = heights[x])
}

pdf(paste0(output_dir, 'overviewFigure_clinical_genomics.pdf'), width = 11, height = 10, useDingbats = F)
plot(full_plot)
dev.off()

#plot legend
pdf(paste0(output_dir, 'Legend_overviewFigure_clinical_genomics.pdf'), width = 5, height = 23, useDingbats = F)
plot(aligend.legends)
dev.off()
