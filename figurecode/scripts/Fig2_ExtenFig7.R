######################################################################################################################
#############               Analysis of Tx421 Smokers without detectable smoking mutagenesis             ############# 
######################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk)

# Description:
# Script to create Figure 2 of the manuscript "The natural history of NSCLC in TRACERx"

#options and libraries
options(stringsAsFactors = F)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(rjson)
library(tidyverse)
library(fst)
library(Biostrings)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(MASS)
library(car)
library(data.table)

# setwd('/Volumes/proj-tracerx-lung/tctProjects/frankella/repos/shared_repos/TRACERx421_data_code/20220808_Tumour_evoutionary_histories/')

#parameters 
output_dir  <- '.'
# output_dir  <- '../../../../Tx_exome/Plots/Fig2_suppfig8/20221116/'

problematic_neverSmokers <- c('CRUK0054', 'CRUK0093', 'CRUK0485', 'CRUK0788')

#thresholds 
minMut_active   <- 20
maxMut_inactive <- 50
active_cutoff   <- 0.3
inactive_cutoff <- 0.1
active   <- 'and'
inactive <- 'and'
# --> active:   mut > 20 and weight > 0.3
#     inactive: mut < 50 and weight < 0.1

#load clinical per patient
clinicalData <- readRDS('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_patient_df.rds')

#load clinical data per tumour
all_combined_df <- readRDS('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_tumour_df.rds')
all_combined_df$tumour_id  <- all_combined_df$tumour_id_muttable_cruk
all_combined_df$patient_id <- all_combined_df$cruk_id

#load mutational signatures
signature_weights <- readRDS('../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_mutationSignature_weights.rds')
total_weights     <- signature_weights$signature_weights_perTumour
clonal_weights    <- signature_weights$signature_weights_clonalMuts_perTumour
subclonal_weights <- signature_weights$signature_weights_subclonalMuts_perTumour

signature_counts <- readRDS('../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_mutationSignature_counts.rds')
total_counts     <- signature_counts$signature_counts_perTumour
clonal_counts    <- signature_counts$signature_counts_clonalMuts_perTumour
subclonal_counts <- signature_counts$signature_counts_subclonalMuts_perTumour

#load mutation data
muttable_df <- read_fst('../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst')

#load fusions
fusions <- read.table('../20221109_Tumour_evo_histories_DATA/20221116_TRACERx421_fusions.txt', header = T, sep = '\t')

#load evolutionary metrics
evo_metrics <- read.table('../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_evolutionary_metrics.tsv', header = T, sep = '\t')

#load eCOMSIC signature profiles (v3.2)
comsic_sigs <- read.table('../20221109_Tumour_evo_histories_DATA/COSMIC_v3.2_SBS_GRCh37.txt', header = T, sep = '\t')



###################################
######       Functions       ######   
###################################

#function to plot SBS profile
plot_SBS <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$group <- substr(as.character(signatures_df$channel), start = 3, stop = 5)
  signatures_df <- signatures_df %>%
    mutate(channel = factor(channel, levels = signatures_df$channel),
           group = factor(group, levels = unique(signatures_df$group)),
           value = value * 100)
  
  #set colours
  colours   <- setNames(c('#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'), unique(signatures_df$group))
  strip_name_colours <- c('black','white','white','black','black','black')
  xlabels   <- paste0(substr(as.character(signatures_df$channel), start = 1, stop = 1),
                      substr(as.character(signatures_df$channel), start = 3, stop = 3),
                      substr(as.character(signatures_df$channel), start = 7, stop = 7))
  
  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~ group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours, guide = 'none') +
    xlab('') +
    ylab('% SBS') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(breaks = as.character(signatures_df$channel), labels = xlabels) +
    ggtitle(title) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                            strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    
    k <- k+1
  }
  
  return(g)
}



#function to calculate fractions of SBS4 detection groups for increasing smoke variables
frac_SBS4detection_valueGroups <- function(group_df, gap = 5){
  
  #exclude NA's
  group_df <- group_df[!is.na(group_df[,2]),]
  
  #create sequence-vector
  min_value <- ifelse(min(group_df[,2])%%gap == 0, min(group_df[,2]), min(group_df[,2]) + (gap - min(group_df[,2])%%gap))
  max_value <- ifelse(max(group_df[,2])%%gap == 0, max(group_df[,2]), max(group_df[,2]) + (gap - max(group_df[,2])%%gap))
  seq_vector <- seq(min_value, max_value, gap)
  
  #calculate fraction for value groups
  frac <- lapply(seq_vector, function(x){
    sub <- group_df[group_df[,2] < x,]
    counts <- sub %>%
      group_by(sub[,1]) %>%
      summarise(n = n()) %>%
      mutate(n_group = nrow(sub))
    return(data.frame(max_value = x, counts))
  })
  frac <- Reduce(rbind, frac)
  colnames(frac) <- c('max_groupValue', 'SBS4detected', 'n', 'n_group')
  frac$frac_group <- frac$n / frac$n_group
  
  #calculate fraction for SBS4detected groups per value group
  temp <- group_df %>%
    group_by(group_df[,1]) %>%
    summarise(n_SBS4group = n()) 
  colnames(temp)[1] <- 'group'
  frac <- frac %>% left_join(temp, by = c('SBS4detected' = 'group'))
  frac$frac_SBS4group <- frac$n / frac$n_SBS4group
  
  return(frac)
}



##############################
######       Main       ######   
##############################

#----------------------------------#
# Investigating of SBS4 and SBS92  #
#----------------------------------#

#SBS4 and SBS92 signature profile
SBS4profile  <- plot_SBS(comsic_sigs[,c('Type', 'SBS4')], title = 'SBS4')
SBS92profile <- plot_SBS(comsic_sigs[,c('Type', 'SBS92')], title = 'SBS92')

pdf(paste0(output_dir, "SBS4_SBS92_profiles.pdf"), width = 11, height = 12)
ggdraw() + 
  draw_plot(SBS4profile, x = 0, y = 0.75, width = 1, height = 0.25) + 
  draw_plot(SBS92profile, x = 0, y = 0.5, width = 1, height = 0.25) 
dev.off()

#correlation of clonal SBS4 and SBS92 weights with pack years in LUAD and LUSC
all_combined_df$histology <- 'Other'
all_combined_df$histology[all_combined_df$Histology_per_tumour_id_muttable %in% 'Invasive adenocarcinoma'] <- 'LUAD'
all_combined_df$histology[all_combined_df$Histology_per_tumour_id_muttable %in% 'Squamous cell carcinoma'] <- 'LUSC'

all_combined_df$tumour_id <- all_combined_df$tumour_id_muttable

all_combined_df <- full_join(all_combined_df, clinicalData)
smoke_vars <- all_combined_df[, c('tumour_id', "years_smoking","cigs_perday","age", 'histology', 'smoking_status_merged', 'packyears')]

plot_data_weights <- smoke_vars %>% 
  dplyr::select(tumour_id, histology, packyears) %>%
  left_join(clonal_weights[,c('tumour_id', 'SBS4', 'SBS92')], by = 'tumour_id') %>%
  dplyr::filter(histology %in% c('LUAD', 'LUSC'))
plot_data_weights <- reshape2::melt(plot_data_weights, id.vars = c('tumour_id', 'histology', 'packyears'))

pdf(paste0(output_dir, 'cor_packyears_clonalSmokeSigs.pdf'), width = 6, height = 5)
ggplot(plot_data_weights, aes(x = packyears, y = value)) + 
  geom_point(aes(colour = histology), alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  stat_cor() +
  scale_y_continuous( breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20 )) +
  scale_colour_manual(name = 'Cancer Type', values = c('LUAD' = '#67001f', 'LUSC' = '#053061')) +
  facet_grid(histology ~ variable, scales = 'free_x') + 
  xlab('Pack Years') + ylab('Clonal Exposure Weights') +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))
dev.off()



#-----------------------------------#
# Detections of Smoking Muatgenesis #
#-----------------------------------#

#subset for SBS4 muttaions and combined signature weights and counts in ine data frame
total_df      <- total_weights[,c('tumour_id', 'SBS4')] %>% left_join(total_counts[,c('tumour_id', 'SBS4')], by = 'tumour_id')
colnames(total_df) <- c('tumour_id', 'weights', 'counts')

clonal_df      <- clonal_weights[,c('tumour_id', 'SBS4')] %>% left_join(clonal_counts[,c('tumour_id', 'SBS4')], by = 'tumour_id')
colnames(clonal_df) <- c('tumour_id', 'weights', 'counts')

subclonal_df      <- subclonal_weights[,c('tumour_id', 'SBS4')] %>% left_join(subclonal_counts[,c('tumour_id', 'SBS4')], by = 'tumour_id')
colnames(subclonal_df) <- c('tumour_id', 'weights', 'counts')

#split never smoker and smoker
neverSmokers <- all_combined_df$tumour_id[all_combined_df$smoking_status_merged == 'Never Smoked']
neverSmokers <- neverSmokers[!is.na(neverSmokers)]
neverSmokers <- neverSmokers[!neverSmokers %in% problematic_neverSmokers]
neverSmokers <- unique(neverSmokers)

neverSmk_clonal    <- clonal_df[clonal_df$tumour_id %in% neverSmokers,]
neverSmk_subclonal <- subclonal_df[subclonal_df$tumour_id %in% neverSmokers,]

smk_clonal    <- clonal_df[!clonal_df$tumour_id %in% neverSmokers,]
smk_subclonal <- subclonal_df[!subclonal_df$tumour_id %in% neverSmokers,]

#plot activities and thresholds --> Extended Figure 7d
plot_data <- rbind(data.frame(total_df, clonality = 'total'),
                   data.frame(clonal_df, clonality = 'clonal'),
                   data.frame(subclonal_df, clonality = 'subclonal'))
plot_data$smkStatus <- 'Smokers'
plot_data$smkStatus[plot_data$tumour_id %in% neverSmokers] <- 'neverSmokers'

plot_data$active <- 'MAYBE'
if(active == 'and'){
  plot_data$active[plot_data$weights > active_cutoff & plot_data$counts > minMut_active] <- 'YES'
} else if(active == 'or'){
  plot_data$active[plot_data$weights > active_cutoff | plot_data$counts > minMut_active] <- 'YES'
}

if(inactive == 'and'){
  plot_data$active[plot_data$weights < inactive_cutoff & plot_data$counts < maxMut_inactive] <- 'NO'
} else if(inactive == 'or'){
  plot_data$active[plot_data$weights < inactive_cutoff | plot_data$counts < maxMut_inactive] <- 'NO'
}

pdf(paste0(output_dir, 'scatter_deconstructSigs_SBS4activity.pdf'), width = 6, height = 4, useDingbats = F)
ggplot(plot_data, aes(x = weights, y = counts, colour = active)) + 
  geom_point() + 
  geom_vline(xintercept = active_cutoff) + 
  geom_vline(xintercept = inactive_cutoff) + 
  geom_hline(yintercept = maxMut_inactive, linetype = 'dotted', colour = 'black') + 
  geom_hline(yintercept = minMut_active, linetype = 'dotted', colour = 'black') + 
  facet_grid(smkStatus ~ clonality, scales = 'free_y') +
  scale_y_continuous(trans = 'log10') +
  scale_color_manual(name = 'SBS4 detected', values = c('YES' = '#b2182b', 'NO' = '#2166ac', 'MAYBE' = '#bababa')) +
  xlab('SBS4 weights') + ylab('SBS4 counts') + 
  ggtitle('DeconstructSigs for clonal and subclonal Mutations') +
  theme_bw()
dev.off()


#classify patients into SBS4 detected and not detected based on deconstructSigs outputs --> mainly use clonal classification
SBS4detection_class <- plot_data %>%
  dplyr::select(tumour_id, clonality, active, smkStatus) %>%
  spread(clonality, active, fill = 'NO')
colnames(SBS4detection_class) <- c('tumour', 'smkStatus', 'clonal_SBS4detected', 'subclonal_SBS4detected', 'total_SBS4detected')

overall_SBS4detetced <- sapply(1:nrow(SBS4detection_class), function(i){
  clonal    <- SBS4detection_class$clonal_SBS4detected[i]
  subclonal <- SBS4detection_class$subclonal_SBS4detected[i]
  total     <- SBS4detection_class$total_SBS4detected[i]
  if(clonal == 'YES'){
    return('YES')
  }
  if(any(total %in% c('MAYBE', 'NO')) & clonal == 'NO' & any(subclonal %in% c('NO', 'MAYBE'))){
    return('NO')
  }
  return('MAYBE')
})
SBS4detection_class$overall_SBS4detetced <- overall_SBS4detetced
write.table(SBS4detection_class, file = paste0(output_dir, 'SBS4detection_classification_perTumour.txt'), sep = '\t', row.names = F, quote = F)



#----------------------------#
# Analysis of Smoking Groups #
#----------------------------#
SBS4detection_df <- SBS4detection_class

#abbreviate histology
all_combined_df$histology <- 'Other'
all_combined_df$histology[all_combined_df$Histology_per_tumour_id_muttable %in% 'Invasive adenocarcinoma'] <- 'LUAD'
all_combined_df$histology[all_combined_df$Histology_per_tumour_id_muttable %in% 'Squamous cell carcinoma'] <- 'LUSC'

#plot fraction of detected and non detected SBS4 tumours of Smokers relative to smokeYrs
smoke_vars <- all_combined_df[, c('tumour_id', "years_smoking","cigs_perday","age", 'histology', 'smoking_status_merged')]
smoke_vars <- SBS4detection_df %>%
  left_join(smoke_vars, by = c('tumour' = 'tumour_id')) %>%
  filter(histology == 'LUAD' & smoking_status_merged != 'Never Smoked') %>%
  mutate(overall_SBS4detetced = factor(overall_SBS4detetced, levels = c('YES', 'MAYBE', 'NO')),
         clonal_SBS4detected = factor(clonal_SBS4detected, levels = c('YES', 'MAYBE', 'NO')),
         subclonal_SBS4detected = factor(subclonal_SBS4detected, levels = c('YES', 'MAYBE', 'NO')))
colnames(smoke_vars) <- c('tumour', "Patient", 'clonal_SBS4detected', 'subclonal_SBS4detected', 'total_SBS4detected', "overall_SBS4detetced", 'smokeYrs', 'cigsPerDay', 'age', 'histology','smoking_status_merged')

frac_smokeYrs_overall <- frac_SBS4detection_valueGroups(group_df = smoke_vars[smoke_vars$overall_SBS4detetced != 'MAYBE',c('overall_SBS4detetced', 'smokeYrs')], gap = 5)

add_data <- frac_smokeYrs_overall[1,,drop = F]
add_data$SBS4detected <- 'YES'
add_data$frac_group   <- 0

plot_data <- rbind(frac_smokeYrs_overall, add_data)

pdf(paste0(output_dir, 'bar_SBS4detection_smokeYrs_LUAD.pdf'), width = 4.5, height = 3.6,  useDingbats = F)
ggplot(plot_data, aes(x = factor(max_groupValue), y = frac_group, fill = SBS4detected)) +
  geom_bar(stat = 'identity', position = "dodge") +
  geom_smooth(method = 'loess', aes(group = SBS4detected, colour = SBS4detected), linetype = 'dashed', se = F, size = 0.8) +
  scale_fill_manual(name = 'SBS4', values = c('NO' = '#4d4d4d', 'MAYBE' = '#fdae61', 'YES' = '#b2182b'),
                    labels = c('NO' = 'not detected', 'MAYBE' = 'low confidence\ndetected', 'YES' = 'high confidence\ndetected')) +
  scale_colour_manual(name = 'SBS4', values = c('NO' = '#4d4d4d', 'MAYBE' = '#fdae61', 'YES' = '#b2182b'),
                      labels = c('NO' = 'not detected', 'MAYBE' = 'low confidence\ndetected', 'YES' = 'high confidence\ndetected')) +
  scale_y_continuous( breaks = seq(0, 1, 0.2), labels = seq(0,100, 20)) +
  xlab('upper limit smoking years') + ylab('% tumours') +
  ggtitle('Probability of detecting SBS4 by smoking duration') +
  theme_bw() + 
  theme(legend.position = 'none', 
        plot.title = element_text(hjust = 0.5, size = 11.5), 
        panel.grid = element_blank(),
        text = element_text( size = 16.5 ))
dev.off()

#comparison of  pack years between high, low and non detected SBS4 groups
smoke_vars <- all_combined_df[, c('tumour_id', "years_smoking","cigs_perday","age", 'histology', 'smoking_status_merged', 'packyears')]
smoke_vars <- SBS4detection_df %>%
  left_join(smoke_vars, by = c('tumour' = 'tumour_id')) %>%
  filter(histology == 'LUAD' & smoking_status_merged != 'Never Smoked') %>%
  mutate(overall_SBS4detetced = factor(overall_SBS4detetced, levels = c('YES', 'MAYBE', 'NO')),
         clonal_SBS4detected = factor(clonal_SBS4detected, levels = c('YES', 'MAYBE', 'NO')),
         subclonal_SBS4detected = factor(subclonal_SBS4detected, levels = c('YES', 'MAYBE', 'NO')))
smoke_vars <- smoke_vars %>% 
  left_join(clonal_df, by = c('tumour' = 'tumour_id')) %>%
  left_join(subclonal_df, by = c('tumour' = 'tumour_id')) %>%
  left_join(total_df, by = c('tumour' = 'tumour_id'))
colnames(smoke_vars) <- c('tumour', "smkStatus", 'clonal_SBS4detected', 'subclonal_SBS4detected', 'total_SBS4detected', "overall_SBS4detetced", 'smokeYrs', 'cigsPerDay', 'age', 'histology','smoking_status_merged', 'packyears',
                          'clonal_SBS4weights', 'clonal_SBS4counts', 'subclonal_SBS4weights', 'subclonal_SBS4counts', 'total_SBS4weights', 'total_SBS4counts')


nClonal         <- data.frame(tumour_id = clonal_counts$tumour_id, counts = rowSums(clonal_counts[,-1]))
plot_data       <- smoke_vars[,c('tumour', 'overall_SBS4detetced', 'clonal_SBS4weights','packyears')]
plot_data       <- plot_data[!is.na(plot_data$packyears),]
plot_data       <- plot_data %>% left_join(nClonal, by = c('tumour' = 'tumour_id'))

pdf(paste0(output_dir, 'box_packyears_SBS4status_LUADeverSmokers.pdf'), width = 5, height = 5, useDingbats = F)
ggplot(plot_data, aes(x = overall_SBS4detetced, y = packyears)) + 
  ggbeeswarm::geom_quasirandom(aes(size = counts, color = clonal_SBS4weights), alpha = 0.5) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  scale_colour_viridis_c(name = '% truncal \nSBS4-assigned\n mutations', option = "plasma") +
  scale_size_continuous(name = '# truncal mutations') +
  scale_x_discrete(labels = c('NO' = 'not detected', 'MAYBE' = 'low confidence\ndetected', 'YES' = 'high confidence\ndetected')) +
  ylab('pack years') + xlab('') +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#how many non detected SBS4 patients have smoked more than median pack years of high confidence SBS4 detected patients?
thresh <- smoke_vars %>% 
  filter(overall_SBS4detetced == 'YES') %>%
  filter(histology == 'LUAD') %>%
  summarise(median(packyears, na.rm = T))

smoke_vars %>% 
  filter(overall_SBS4detetced == 'NO') %>% 
  filter(packyears > as.numeric(thresh)) %>%
  filter(histology == 'LUAD')


# Association with Fusions & EGFR 
EGFRmut <- muttable_df[muttable_df$Hugo_Symbol == 'EGFR',]
EGFRmut <- EGFRmut[as.logical(EGFRmut$DriverMut),]
EGFRmut_tumours <- unique(EGFRmut$tumour_id)

fusion_data <- SBS4detection_df %>% left_join(all_combined_df[,c('tumour_id', 'patient_id', 'histology', 'smoking_status_merged')], by = c('tumour' = 'tumour_id'))
fusion_data <- fusion_data %>% left_join(fusions[!duplicated(fusions$patient_id), c('patient_id', 'Genes', 'Cat')], by = 'patient_id')
fusion_data$Cat[fusion_data$tumour %in% EGFRmut_tumours] <- 'EGFR'

count_neverSmoker <- fusion_data[fusion_data$smoking_status_merged == 'Never Smoked',] %>%
  filter(histology == 'LUAD') %>%
  group_by(Cat) %>%
  summarise(n = n())
count_neverSmoker$fraction <- count_neverSmoker$n / sum(count_neverSmoker$n)

count_smoker <- fusion_data[fusion_data$smoking_status_merged != 'Never Smoked',] %>%
  filter(histology == 'LUAD') %>%
  group_by(Cat, overall_SBS4detetced) %>%
  summarise(n = n())
count_smoker$fraction <- NA
count_smoker$fraction[count_smoker$overall_SBS4detetced == 'YES'] <- count_smoker$n[count_smoker$overall_SBS4detetced == 'YES'] / sum(count_smoker$n[count_smoker$overall_SBS4detetced == 'YES'])
count_smoker$fraction[count_smoker$overall_SBS4detetced == 'NO'] <- count_smoker$n[count_smoker$overall_SBS4detetced == 'NO'] / sum(count_smoker$n[count_smoker$overall_SBS4detetced == 'NO'])
count_smoker$fraction[count_smoker$overall_SBS4detetced == 'MAYBE'] <- count_smoker$n[count_smoker$overall_SBS4detetced == 'MAYBE'] / sum(count_smoker$n[count_smoker$overall_SBS4detetced == 'MAYBE'])

plot_data <- rbind(data.frame(Cat = count_neverSmoker$Cat, overall_SBS4detetced = 'Never Smoked', n = count_neverSmoker$n, fraction = count_neverSmoker$fraction, smkstatus = 'Never Smoked'),
                   data.frame(count_smoker, smkstatus = 'Smoked'))
plot_data$Cat[is.na(plot_data$Cat)] <- 'No Event'
plot_data <- plot_data %>%
  mutate(Cat = factor(Cat, levels = rev(c('No Event', 'Fusion', 'Oncogenic Isoform', 'EGFR'))),
         overall_SBS4detetced = ifelse(overall_SBS4detetced == 'NO', 'not detected', overall_SBS4detetced),
         overall_SBS4detetced = ifelse(overall_SBS4detetced == 'MAYBE', 'low confidence\ndetected', overall_SBS4detetced),
         overall_SBS4detetced = ifelse(overall_SBS4detetced == 'YES', 'high confidence\ndetected', overall_SBS4detetced),
         overall_SBS4detetced = factor(overall_SBS4detetced, levels = c('Never Smoked', 'not detected', 'low confidence\ndetected', 'high confidence\ndetected')))


pdf(paste0(output_dir, 'bar_fusion_EGFRmut_SBS4detection_LUAD.pdf'), width = 5, height = 5, useDingbats = F)
ggplot(plot_data, aes(x = overall_SBS4detetced, y = fraction, fill = Cat)) + 
  geom_bar(stat = 'identity', colour = 'black') +
  scale_fill_manual(name = '', values = c('EGFR' = '#8c6bb1',  'Fusion' = '#01665e', 'Oncogenic Isoform' = '#35978f', 'No Event' = '#e0e0e0'),
                    labels = c('EGFR' = 'EGFRmut', 'Fusions' = 'Fusions', 'Oncogenic Isoform' = 'Met exon 14 skipping', 'No Event' = 'No Event')) +
  geom_text(aes(label = n, colour = Cat), position = position_stack(vjust = 0.5), show.legend = F) +
  scale_colour_manual(values = c(c('EGFR' = 'white', 'Fusion' = 'white', 'Oncogenic Isoform' = 'white', 'No Event' = 'black'))) +
  facet_grid(.~smkstatus, scales = 'free_x', space = 'free_x') +
  scale_y_continuous(expand = c(0,0, 0, 0.02), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  xlab('') + ylab('% tumours') + 
  ggtitle('Fusions and Onocgenic Isoforms') +
  theme_bw() + theme(legend.position = 'top', plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 30, hjust = 1),
                     panel.grid = element_blank(), text = element_text( size = 15 ))
dev.off()

#fishers test
EGFR_test <- plot_data %>%
  filter(overall_SBS4detetced %in% c('not detected', 'high confidence\ndetected')) %>%
  dplyr::select(overall_SBS4detetced, Cat, n) %>%
  tidyr::spread(Cat, n) 
EGFR_test$nonEGFR <- rowSums(EGFR_test[,3:ncol(EGFR_test)], na.rm = T)
fisher.test(as.matrix(EGFR_test[,c('EGFR', 'nonEGFR')]))

fusion_test <- plot_data %>%
  filter(overall_SBS4detetced %in% c('not detected', 'high confidence\ndetected')) %>%
  dplyr::select(overall_SBS4detetced, Cat, n) %>%
  tidyr::spread(Cat, n, fill = 0) 
fusion_test$Fusion <- fusion_test$Fusion + fusion_test$`Oncogenic Isoform`
fusion_test$nonFusion <- rowSums(fusion_test[,c(2,5)], na.rm = T)
fisher.test(as.matrix(fusion_test[,c('Fusion', 'nonFusion')]))



#Association GD events 
gds <- evo_metrics[,c('tumour_id', 'First_GD', 'Second_GD', 'num_first_gd', 'num_second_gd', 'First_GD_homogen', 'Second_GD_homogen',
                      'GD_status_homogen', 'GD_statuses', 'frac_0_gd_regions', 'frac_1_gd_regions', 'frac_2_gd_regions', 'num_clonal_gds', 'num_subclonal_gds')]
gds <- data.table(gds)
gds[ num_clonal_gds == 0 & num_subclonal_gds == 0, gd_status := 'No GD']
gds[ num_clonal_gds == 1  & num_subclonal_gds == 0,  gd_status := '1 Truncal GD']
gds[ num_clonal_gds == 2 & num_subclonal_gds == 0,  gd_status := '2 Truncal GD']
gds[ num_clonal_gds == 0 & num_subclonal_gds == 1, gd_status := '1 Subclonal GD']
gds[ num_clonal_gds == 0 & num_subclonal_gds > 1,  gd_status := 'Multi-Subclonal GD']
gds[ num_clonal_gds == 1 & num_subclonal_gds == 1,  gd_status := '1 Truncal GD & Subclonal GD']
gds[ num_clonal_gds == 1 & num_subclonal_gds > 1,   gd_status := '1 Truncal GD & Multi-Subclonal GD']
colours <- c('No GD'='#e0e0e0', '1 Truncal GD'='#9999FF', '2 Truncal GD'='#0000CC', 
             '1 Truncal GD & Subclonal GD'='#CC66FF', '1 Truncal GD & Multi-Subclonal GD'='#330066', 
             '1 Subclonal GD'='#FF6666', 'Multi-Subclonal GD'='#990000')

gds[, gd_status := factor(gd_status, levels = rev(names(colours)))]
gds <- gds %>%
  left_join(all_combined_df[, c('tumour_id', 'histology', 'smoking_status_merged')]) %>%
  left_join(SBS4detection_df[,c('tumour', 'overall_SBS4detetced')], by = c('tumour_id' = 'tumour')) %>%
  filter(histology == 'LUAD') %>%
  filter(!is.na(gd_status))

plot_data <- gds
plot_data$class <- plot_data$overall_SBS4detetced
plot_data$class[plot_data$smoking_status_merged == 'Never Smoked'] <- 'Never Smoked'
plot_data$class <- sub('NO', 'not detected', plot_data$class)
plot_data$class <- sub('YES', 'high confidence\ndetected', plot_data$class)
plot_data$class <- sub('MAYBE', 'low confidence\ndetected', plot_data$class)
plot_data$class <- factor(plot_data$class, levels = c('Never Smoked', 'not detected', 'low confidence\ndetected', 'high confidence\ndetected'))

plot_data <- plot_data %>%
  group_by(class, gd_status) %>%
  summarise(count = n()) %>%
  group_by(class) %>%
  mutate(total = sum(count),
         fraction = count / total) %>%
  data.frame()

plot_data$group <- 'Smoked'
plot_data$group[plot_data$class == 'Never Smoked'] <- 'Never Smoked'

pdf(paste0(output_dir, 'bar_GD_SBS4detection_LUAD.pdf'), width = 8, height = 5, useDingbats = F)
ggplot(plot_data, aes(x = class, y = fraction, fill = gd_status)) + 
  geom_bar(stat = 'identity', colour = 'black') +
  scale_fill_manual(name = '', values = rev(colours)) +
  geom_text(aes(label = count, colour = gd_status), position = position_stack(vjust = 0.5), show.legend = F) +
  scale_colour_manual(values = c('No GD'='black', '1 Truncal GD'='white', '2 Truncal GD'='white', 
                                 '1 Truncal GD & Subclonal GD'='white', '1 Truncal GD & Multi-Subclonal GD'='white', 
                                 '1 Subclonal GD'='white', 'Multi-Subclonal GD'='white')) +
  facet_grid(.~group, scales = 'free_x', space = 'free_x') +
  scale_y_continuous(expand = c(0,0, 0, 0.02), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  xlab('') + ylab('% tumours') + 
  ggtitle('') +
  theme_bw() + theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 30, hjust = 1),
                     panel.grid = element_blank(), text = element_text( size = 15))
dev.off()

#fishers test
GD_test <- plot_data %>%
  mutate(cat = ifelse(class %in% c('Never Smoked', 'not detected'), 'noSmoker', 'Smoker')) %>%
  mutate(group = ifelse(gd_status == 'No GD', 'No GD', 'GD')) %>%
  group_by(cat, group) %>%
  summarise(n = sum(count)) %>%
  tidyr::spread(group, n) 
fisher.test(as.matrix(GD_test[,c('GD', 'No GD')]))




#-----------------------------------------#
# Analysis of lobe specific SBS4 exposure #
#-----------------------------------------#

#organise input data 
clonal_SBS4_df <- clonal_counts %>%
  dplyr::select(tumour_id, SBS4) %>%
  dplyr::rename(SBS4_muts_clonal = SBS4)


all_combined_df  <- all_combined_df %>%
  left_join(clonal_SBS4_df, by = "tumour_id") %>%
  mutate(histology = histology_3 ) %>%
  mutate(lobe = case_when(grepl("Upper", site_per_lesion) ~ "Upper",
                          grepl("Mid", site_per_lesion) ~ "Upper",
                          grepl("Lower", site_per_lesion) ~ "Lower"),
         side = case_when(grepl("Right", site_per_lesion) ~ "Right",
                          grepl("Left", site_per_lesion) ~ "Left")) %>%
  mutate(lobe = factor(lobe, levels = c("Lower", "Upper")),
         side = factor(side, levels = c("Left", "Right"))) %>%
  mutate(LungLobe = case_when(site_per_lesion == "Left Upper lobe" ~ "Left Upper Lobe", 
                              TRUE ~ site_per_lesion))  %>%
  mutate(LungLobe = factor(LungLobe, levels = c("Left Lower Lobe", "Left Upper Lobe","Right Lower Lobe","Right Middle Lobe", "Right Upper Lobe"  )))%>%
  dplyr::rename(Age = age,
                Lobe_UM = lobe, # ref = Lower
                Side_R = side, # ref = left
                CigsPerDay = cigs_perday,
                YearsSmoking = years_smoking,
                packyears = packyears,
                FamilyHistory = is.family.lung)


#lobe specific clonal SBS4 counts in LUAD & LUSC
in_df <- all_combined_df %>%
  filter(smoking_status_merged != "Never Smoked") %>%
  mutate(SmokingYears = YearsSmoking,
         ClonalSBS4counts = round(SBS4_muts_clonal) )

table_all <- in_df %>% 
  filter(histology %in% c("LUAD", "LUSC")) %>%
  group_by(LungLobe, histology) %>%
  summarise(MedianClonalSBS4 = median(ClonalSBS4counts, na.rm=T),
            TumourCount = n()) %>%
  mutate(mean.TumourCount.per.segment = case_when(LungLobe == "Left Upper Lobe" ~ TumourCount/4,
                                                  LungLobe == "Left Lower Lobe" ~ TumourCount/4,
                                                  LungLobe == "Right Upper Lobe" ~ TumourCount/3,
                                                  LungLobe == "Right Middle Lobe" ~ TumourCount/2,
                                                  LungLobe == "Right Lower Lobe" ~ TumourCount/5))

#plot lobe specific clonal SBS4 counts
pdf(paste0(output_dir,"/median.clonal.SBS4_by_lobe.pdf"), width = 12, height = 5)
ggplot(table_all, aes(x= LungLobe, y = 1)) +
  geom_point(aes(fill = MedianClonalSBS4, size = TumourCount), colour="black", pch=21, alpha = 0.99) +
  theme_cowplot() +
  ggtitle("(Ex)Smokers") +
  scale_fill_gradient(low = "yellow", high = "red")+
  labs(size = "Mean tumour count") +
  theme(axis.text.x = element_text(size=16, angle = 45,  hjust=1))+
  facet_grid(. ~ histology)
dev.off()
#--> Add  dots to lung schematic from Biorender manually


#--------------------------------------------#
# Regression analysis for clonal SBS4 counts #
#--------------------------------------------#
#--> Extended Figure 7 a and b

#organise input data 
clinicaldata_selected <- clinicalData %>%
  dplyr::select(cruk_id , age, sex, cigs_perday, years_smoking ,is.family.lung ,smoking_status_merged )

input_regression  <- all_combined_df %>%
  dplyr::select(tumour_id,cruk_id, site_per_lesion, histology_3 ) %>%
  left_join(clonal_SBS4_df, by = "tumour_id") %>%
  left_join(clinicaldata_selected , by = "cruk_id" ) %>%
  mutate(histology = histology_3 ) %>%
  mutate(lobe = case_when(grepl("Upper", site_per_lesion) ~ "Upper",
                          grepl("Mid", site_per_lesion) ~ "Upper",
                          grepl("Lower", site_per_lesion) ~ "Lower"),
         side = case_when(grepl("Right", site_per_lesion) ~ "Right",
                          grepl("Left", site_per_lesion) ~ "Left")) %>%
  mutate(lobe = factor(lobe, levels = c("Lower", "Upper")),
         side = factor(side, levels = c("Left", "Right"))) %>%
  mutate(LungLobe = case_when(site_per_lesion == "Left Upper lobe" ~ "Left Upper Lobe", 
                              TRUE ~ site_per_lesion))  %>%
  mutate(LungLobe = factor(LungLobe, levels = c("Left Lower Lobe", "Left Upper Lobe","Right Lower Lobe","Right Middle Lobe", "Right Upper Lobe"  )))%>%
  dplyr::rename(Age = age,
                Lobe_UM = lobe, # ref = Lower
                Side_R = side, # ref = left
                Gender_M = sex, # ref = female
                CigsPerDay = cigs_perday,
                YearsSmoking = years_smoking,
                FamilyHistory = is.family.lung) %>%
  filter(smoking_status_merged != "Never Smoked")


### regression analysis ###

# set dependant variable
dep_var <- "SBS4_muts_clonal"

# set explanatory variables
variables_MVA <- c("Age" , "Gender_M", "Side_R","Lobe_UM", "CigsPerDay", "YearsSmoking","FamilyHistory") 


## run for LUAD ##

# exclude NAs
variables_MVA.ori <- c(dep_var, variables_MVA) 
in_df_nona <- input_regression %>%
  filter(histology == "LUAD") 

in_df_nona$dep_var <- in_df_nona[[dep_var]]

in_df_nona <- as.data.table(in_df_nona) %>%
  .[,..variables_MVA.ori]  %>% 
  na.omit() 

## formula for model
f <- as.formula(paste(dep_var,
                      paste(variables_MVA, collapse = " + "), 
                      sep = " ~ "))

## beta binomial regression 
res_glm <- glm.nb(f , data = in_df_nona) 

summary(res_glm)

#rsq <- rsq(res_glm, adj=FALSE, type=c('n'))  # type=c('v','kl','sse','lr','n')  defalt : "v"  * 'n' will not accept adj = TRUE
#rsq   # [1]  0.3412369

# 95%CI 
Confint(res_glm)

# rate ratio 
rate_ratio.ori <- as.data.frame(t(exp(Confint(res_glm))))

rate_ratio <- rate_ratio.ori %>%
  mutate(Age = Age^10,
         CigsPerDay = CigsPerDay^10,
         YearsSmoking = YearsSmoking^10) %>%
  dplyr::rename(Gender_M = Gender_MMale,
                Side_R = Side_RRight,
                Lobe_UM = Lobe_UMUpper,
                FamilyHistory = FamilyHistoryTRUE) %>%
  dplyr::select(Age  ,Gender_M ,Side_R, Lobe_UM, CigsPerDay ,YearsSmoking, FamilyHistory ) %>%
  t()

# type 1 anova
anova_table <- anova(res_glm )  
anova_table
# coefficient in the reduced model
coef_vec <- coef(summary(res_glm))[,1][-1]
names(coef_vec) <- c(rownames(anova_table)[-1])
# coef pval 
coef_pval_vec <- coef(summary(res_glm))[,4][-1]
names(coef_pval_vec) <- c(rownames(anova_table)[-1])

# type 2 anova  
anova_table2 <- car::Anova(res_glm)  
anova_table2

anova_pval2 <- anova_table2$`Pr(>Chisq)`  
names(anova_pval2) <- c(rownames(anova_table2)) 


## plot data for LUAD
res_df_LUAD <- data.frame(variable = factor(names(anova_pval2) ,levels = rev(names(anova_pval2))),
                          coefficient = coef_vec, 
                          rate.ratio = rate_ratio[,1],
                          lower_rate.ratio = rate_ratio[,2],
                          upper_rate.ratio = rate_ratio[,3],
                          coef_pval = coef_pval_vec,
                          pval2 = anova_pval2,
                          stringsAsFactors = F
) 

res_df_LUAD <- res_df_LUAD %>%
  mutate(col_coef = case_when(coefficient > 0 ~ "pos",
                              coefficient < 0 ~ "neg") ,
         col_coef_sig2 = case_when(coefficient > 0 & pval2 < 0.05 ~ "pos",
                                   coefficient < 0 & pval2 < 0.05 ~ "neg"),
         sig2 = case_when(pval2 < 0.001 ~ "***",
                          pval2 < 0.01 ~ "**",
                          pval2 < 0.05 ~ "*"),
         sig3 = case_when(coef_pval < 0.001 ~ "***",
                          coef_pval < 0.01 ~ "**",
                          coef_pval < 0.05 ~ "*")
         
  )

write.csv(res_df_LUAD, paste0(output_dir,"/regression_clonal.SBS4_forest.plot_LUAD.csv"), row.names = F)

g1 <- ggplot(data = res_df_LUAD, aes(x = variable, y = rate.ratio, ymin=lower_rate.ratio, ymax=upper_rate.ratio)) +
  geom_pointrange(aes(colour= col_coef_sig2) ,size=1.5 ) + 
  geom_hline(yintercept=1, lty=2) +
  theme_cowplot() +
  labs( y = "Rate Ratio (95%CI)" ) + 
  scale_color_manual(values = c("pos" = "red", "neg" = "blue"), na.value = "#bababa",  drop = F ) +
  theme(axis.title.x = element_text(size=18),  
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),  
        axis.text.y = element_text(size=14),
        title = element_text(size=18)) +
  geom_text(aes(label=sig2), size = 10, y= log10(3.5), nudge_x = -0.2  ) +  #y= 3,
  scale_y_continuous(trans = "log10",  breaks = c(0.5,1,2,4),lim=c(0.4,4)) +
  ggtitle("Clonal SBS4 counts (LUAD)") + 
  coord_flip() +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.title.y=element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size=18),
    plot.title = element_text(size=22)
  ) 

# output in pdf
plot_out <- paste0(output_dir,"/regression_clonal.SBS4_forest.plot_LUAD.pdf")

pdf(plot_out, width = 6, height = 4.5)
print(g1)
dev.off()

# output in eps
out_eps <- paste0(output_dir,"/regression_clonal.SBS4_forest.plot_LUAD.eps")
g1
ggsave(out_eps,  width = 6, height = 4.5)
dev.off()


## run for LUSC ##

# exclude NAs
variables_MVA.ori <- c(dep_var, variables_MVA) 
in_df_nona <- input_regression %>%
  filter(histology == "LUSC")  

in_df_nona$dep_var <- in_df_nona[[dep_var]]

in_df_nona <- as.data.table(in_df_nona) %>%
  .[,..variables_MVA.ori]  %>% 
  na.omit() 

## formula for model
f <- as.formula(paste(dep_var,
                      paste(variables_MVA, collapse = " + "), 
                      sep = " ~ "))

## beta binomial regression 
res_glm <- glm.nb(f , data = in_df_nona) 

summary(res_glm)

#rsq <- rsq(res_glm, adj=FALSE, type=c('n'))  # type=c('v','kl','sse','lr','n')  defalt : "v"  * 'n' will not accept adj = TRUE
#rsq   # [1]  0.1089322

# 95%CI 
Confint(res_glm)

# rate ratio 
rate_ratio.ori <- as.data.frame(t(exp(Confint(res_glm))))

rate_ratio <- rate_ratio.ori %>%
  mutate(Age = Age^10,
         CigsPerDay = CigsPerDay^10,
         YearsSmoking = YearsSmoking^10) %>%
  dplyr::rename(Gender_M = Gender_MMale,
                Side_R = Side_RRight,
                Lobe_UM = Lobe_UMUpper,
                FamilyHistory = FamilyHistoryTRUE) %>%
  dplyr::select(Age  ,Gender_M ,Side_R, Lobe_UM, CigsPerDay ,YearsSmoking, FamilyHistory ) %>%
  t()

# type 1 anova
anova_table <- anova(res_glm )  
anova_table
# coefficient in the reduced model
coef_vec <- coef(summary(res_glm))[,1][-1]
names(coef_vec) <- c(rownames(anova_table)[-1])
# coef pval 
coef_pval_vec <- coef(summary(res_glm))[,4][-1]
names(coef_pval_vec) <- c(rownames(anova_table)[-1])

# type 2 anova  
anova_table2 <- car::Anova(res_glm) 
anova_table2

anova_pval2 <- anova_table2$`Pr(>Chisq)`  # no NULL in Anova()
names(anova_pval2) <- c(rownames(anova_table2)) # remove NULL

## plot data for LUSC
res_df_LUSC <- data.frame(variable = factor(names(anova_pval2) ,levels = rev(names(anova_pval2))),
                          coefficient = coef_vec, 
                          rate.ratio = rate_ratio[,1],
                          lower_rate.ratio = rate_ratio[,2],
                          upper_rate.ratio = rate_ratio[,3],
                          coef_pval = coef_pval_vec,
                          pval2 = anova_pval2,
                          stringsAsFactors = F
) 

res_df_LUSC <- res_df_LUSC %>%
  mutate(col_coef = case_when(coefficient > 0 ~ "pos",
                              coefficient < 0 ~ "neg") ,
         col_coef_sig2 = case_when(coefficient > 0 & pval2 < 0.05 ~ "pos",
                                   coefficient < 0 & pval2 < 0.05 ~ "neg"),
         sig2 = case_when(pval2 < 0.001 ~ "***",
                          pval2 < 0.01 ~ "**",
                          pval2 < 0.05 ~ "*"),
         sig3 = case_when(coef_pval < 0.001 ~ "***",
                          coef_pval < 0.01 ~ "**",
                          coef_pval < 0.05 ~ "*")
  )

write.csv(res_df_LUSC, paste0(output_dir,"/regression_clonal.SBS4_forest.plot_LUSC.csv"), row.names = F)

g1 <- ggplot(data = res_df_LUSC, aes(x = variable, y = rate.ratio, ymin=lower_rate.ratio, ymax=upper_rate.ratio)) +
  geom_pointrange(aes(colour= col_coef_sig2) ,size=1.5 ) + 
  geom_hline(yintercept=1, lty=2) +
  theme_cowplot() +
  labs( y = "Rate Ratio (95%CI)" ) + 
  scale_color_manual(values = c("pos" = "red", "neg" = "blue"), na.value = "#bababa",  drop = F ) +
  theme(axis.title.x = element_text(size=18),   
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),  
        axis.text.y = element_text(size=14),
        title = element_text(size=18)) +
  geom_text(aes(label=sig2), size = 10, y= log10(3.5), nudge_x = -0.2  ) +  #y= 3,
  scale_y_continuous(trans = "log10",  breaks = c(0.5,1,2,4),lim=c(0.4,4)) +
  ggtitle("Clonal SBS4 counts (LUSC)") + 
  coord_flip() +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.title.y=element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size=18),
    plot.title = element_text(size=22)
  ) 

# output in pdf
plot_out <- paste0(output_dir,"/regression_clonal.SBS4_forest.plot_LUSC.pdf")

pdf(plot_out, width = 6, height = 4.5)
print(g1)
dev.off()

# output in eps
out_eps <- paste0(output_dir,"/regression_clonal.SBS4_forest.plot_LUSC.eps")
g1
ggsave(out_eps,  width = 6, height = 4.5)
dev.off()






