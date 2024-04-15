######################################################################################################################
#############                Analysis of frequency of copy number driver events              ############# 
######################################################################################################################
# written by Alexander Frankell (alexander.frankell@crick.ac.uk)

# Description:
# Script to create Extended Figure 8bc of the manuscript "The natural history of NSCLC in TRACERx"
#options and libraries
options(stringsAsFactors = F)
suppressWarnings( suppressPackageStartupMessages( library(fst) ) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(amfFunctions) ) #www.github.com/amf71/amfFunctions
suppressPackageStartupMessages( library(biomaRt) )
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(ggrepel) )
suppressPackageStartupMessages( library(ggpubr) )

#setwd('/Volumes/proj-tracerx-lung/tctProjects/frankella/repos/shared_repos/TRACERx421_data_code/20220808_Tumour_evoutionary_histories//')

#parameters
outputs.folder  <- './'
date <- gsub("-","",Sys.Date())

##############################################
#### Get Inputs required for all analyses ####
##############################################

annot_genes_path <- '../20221109_Tumour_evo_histories_DATA/gene_anno_df.fst'
gistic_path <- '../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_GISTIC_peak_SCNAs.tsv'
clinical_data_path   <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_tumour_df.rds"
all_events_path   <- "../20221109_Tumour_evo_histories_DATA/20221130_TRACERx421_gd_driver_signature_events.tsv"
cn_path   <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_scna_table.fst"


### load cn drivers ###
annot_genes <- fst::read_fst( annot_genes_path, as.data.table=T )
cn_drivers <- annot_genes[ (is_cn_driver), gene_name ]

# load GISTIC data
gistic <- data.table::fread( gistic_path )

### load clinical data ###
clin <- as.data.table( readRDS( clinical_data_path) )

## load genomic events ##
all_events <- fread( all_events_path )
cn_events <- all_events[ consequence %in% c('Amp', 'Del') ]
cn <- read_fst( cn_path, as.data.table=T)

##########################
#### Sup fig 9 b & c ####
##########################

# Just need to label the data with Our CN Driver Genes

#get locations of cn driver genes
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
attributes <- c("hgnc_symbol",'chromosome_name', 
                'start_position', 'end_position')
gene_positions <- as.data.table( getBM(attributes=attributes,mart=ensembl, 
                                       useCache = FALSE) )

# filter for annotated nuclear genes 1-22
gene_positions <- gene_positions[ chromosome_name  %in% 1:22 ]

#annotate onto GISTIC output
chosen_gene_positions <- makeGRangesFromDataFrame( gene_positions[ hgnc_symbol %in% cn_drivers, .(chr = chromosome_name, 
                                                                                                  start = start_position, 
                                                                                                  end = end_position, 
                                                                                                  hgnc_symbol) ], keep.extra.columns = TRUE )

gistic[, chromosome := gsub( 'chr', '', tstrsplit(ID, split = ':')[[1]]) ]
gistic[, start := tstrsplit( tstrsplit(ID, split = ':')[[2]], split = '-' )[[1]] ]
gistic[, end := gsub( '\\(.*$', '', tstrsplit( tstrsplit(ID, split = ':')[[2]], split = '-' )[[2]]) ]

gistic_gr <- makeGRangesFromDataFrame( gistic, keep.extra.columns = TRUE )

overlaps <- findOverlaps(gistic_gr, chosen_gene_positions)

gistic_matched <- cbind( as.data.table(gistic_gr)[ overlaps@from ], 
                         as.data.table(chosen_gene_positions)[ overlaps@to ] )

gistic[, driver := gistic_matched[ match( gistic$ID, gistic_matched$ID ), hgnc_symbol ] ]


# get the frequency of each peak clonally and subclonally in different subsets
# Add Histology - the 'Histology' on there refer to what histology the peaks come from
gistic[, histology_sample := clin[ match( patient_id, cruk_id), Histology_per_tumour_id_muttable ] ]

# Ensure TSGs are only on deletion peaks and OGs on ampifications

gistic[ Type == 'DEL' & driver %in% annot_genes[ Bailey_TSG_Oncogene %in% c('oncogene', 'possible oncogene') & (is_cn_driver), gene_name ], driver := NA ]
gistic[ Type == 'AMP' & driver %in% annot_genes[ Bailey_TSG_Oncogene %in% c('tsg', 'possible tsg') & (is_cn_driver), gene_name], driver := NA ]

all_tumours <- gistic[, length(unique(patient_id))]
gistic_all <- gistic[, .( clonal_freq = sum(Alteration == 'clonal')/all_tumours,
                          subclonal_freq = sum(Alteration == 'subclonal')/all_tumours,
                          overall_freq = sum(!Alteration == 'no_alteration')/all_tumours),
                     by = .(Peak_ID, chromosome, start, end, Source, driver, Type)]

adeno_tumours <- gistic[ histology_sample == 'Invasive adenocarcinoma', length(unique(patient_id))]
gistic_adeno <- gistic[ histology_sample == 'Invasive adenocarcinoma' & Histology == 'LUAD' , 
                        .( clonal_freq = sum(Alteration == 'clonal')/adeno_tumours,
                           subclonal_freq = sum(Alteration == 'subclonal')/adeno_tumours,
                           overall_freq = sum(!Alteration == 'no_alteration')/adeno_tumours),
                        by = .(Peak_ID, chromosome, start, end, Source, driver, Type)]

squam_tumours <- gistic[histology_sample == 'Squamous cell carcinoma', length(unique(patient_id))]
gistic_lusc <- gistic[ histology_sample == 'Squamous cell carcinoma' & Histology == 'LUSC', 
                       .( clonal_freq = sum(Alteration == 'clonal')/squam_tumours,
                          subclonal_freq = sum(Alteration == 'subclonal')/squam_tumours,
                          overall_freq = sum(!Alteration == 'no_alteration')/squam_tumours ),
                       by = .(Peak_ID, chromosome, start, end, Source, driver, Type)]

# define overlapping peaks and remove duplicates if sam ehistology and peak type
gistic_all <- gistic_all[ is.na(driver) | !duplicated( paste(driver, Type) )]
gistic_adeno <- gistic_adeno[ is.na(driver) | !duplicated( paste(driver, Type) )]
gistic_lusc <- gistic_lusc[ is.na(driver) | !duplicated( paste(driver, Type) )]


# add classes
gistic_all[, clonal_subclonal_ratio := clonal_freq / subclonal_freq ]
gistic_adeno[, clonal_subclonal_ratio := clonal_freq / subclonal_freq ]
gistic_lusc[, clonal_subclonal_ratio := clonal_freq / subclonal_freq ]

gistic_all[ clonal_subclonal_ratio < 0.5, class := 'subclonal favoured']
gistic_all[ clonal_subclonal_ratio >= 0.5, class := 'subclonal and truncal']
gistic_all[ clonal_subclonal_ratio > 1.5, class := 'truncal favoured']

gistic_adeno[ clonal_subclonal_ratio < 0.5, class := 'subclonal favoured']
gistic_adeno[ clonal_subclonal_ratio >= 0.5, class := 'subclonal and truncal']
gistic_adeno[ clonal_subclonal_ratio > 1.5, class := 'truncal favoured']

gistic_lusc[ clonal_subclonal_ratio < 0.5, class := 'subclonal favoured']
gistic_lusc[ clonal_subclonal_ratio >= 0.5, class := 'subclonal and truncal']
gistic_lusc[ clonal_subclonal_ratio > 1.5, class := 'truncal favoured']

## plot

plot_function_classes <- function(data, title, add_legend = TRUE){
  
  colours_classes <- c('truncal favoured' = "#377EB8",
                       'subclonal and truncal' = "#984EA3",
                       'subclonal favoured' = "#E41A1C")
  
  clonal_nudge <- data[ class == 'truncal favoured', subclonal_freq ]
  clonal_nudge[ clonal_nudge < 1.5 ] <- 1.5
  
  plot <- ggplot(data, aes( x = clonal_freq, y = subclonal_freq )) +
    geom_point( aes(colour = class, size = overall_freq ) ) +
    scale_colour_manual( values = colours_classes ) +
    #geom_abline(slope = 1) +
    scale_size( limits = c(0.01, 1), 
                breaks = c(0.05, 0.1, 0.2, 0.5, 0.8),
                labels = c(5, 10, 20, 50, 80), range = c(1,5) ) +
    scale_y_continuous( breaks = seq(0,0.4,0.1), labels = seq(0, 40, 10), limits = c(0, 0.4)) +
    scale_x_continuous( breaks = seq(0,0.8,0.1), labels = seq(0, 80, 10), limits = c(0, 0.8)) +
    labs( x =  "Clonal CN event frequency (%)",
          y = "Subclonal CN event frequency (%)",
          fill = '',
          colour = '',
          size = '% Cases Amp/Del',
          title = title) +
    theme_classic() +
    theme( text = element_text( size = 20 ),
           plot.title = element_text( size = 35, hjust = 0.5 )) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    geom_text_repel( data = data[ class == 'truncal favoured' ], direction = "x", segment.alpha = 0.3, xlim = c(0.02, 0.8),
                     aes( label = driver ), box.padding = 1, force_pull = 1, angle = 90, segment.linetype = 2,
                     size = 5,  nudge_y  =   - clonal_nudge ) +
    geom_text_repel( data = data[ !class == 'truncal favoured' ], force_pull = 0, xlim = c(0, 0.8),
                     aes( label = driver, angle = 90 ), box.padding = 1, direction    = "x",  segment.alpha = 0.3, segment.linetype = 2,
                     size = 5, nudge_y  = 40*1.1 -  data[ !class == 'truncal favoured',  subclonal_freq ]  )
  
  if(!add_legend){
    plot <- plot + theme(legend.position = "none") 
  }
  
  print(plot)
}

pdf( paste0( outputs.folder, date, '_lusc_gistic_by_clonality.pdf'), width = 10.5, height = 8, useDingbats = FALSE  )
plot_function_classes( gistic_lusc, title = 'LUSC'  )
invisible( dev.off() )

pdf( paste0( outputs.folder, date, '_luad_gistic_by_clonality.pdf'), width = 8, height = 8, useDingbats = FALSE   )
plot_function_classes( gistic_adeno, title = 'LUAD', add_legend = FALSE  )
invisible( dev.off() )


#############
#### Extended Figure 8c ####
#############

cn[, histology := clin[ match( tumour_id, tumour_id_muttable_cruk), Histology_per_tumour_id_muttable ] ]
cn[ grepl('adenocarcinoma', histology), histology := 'LUAD' ]
cn[ grepl('Squamous', histology), histology := 'LUSC' ]
cn[ !histology %in% c('LUAD', 'LUSC'), histology := 'Other' ]

all_num <- cn[, length(unique(tumour_id)) ]
adeno_num <- cn[ histology == 'LUAD', length(unique(tumour_id)) ]
squam_num <- cn[ histology == 'LUSC', length(unique(tumour_id)) ]
other_num <- cn[ histology == 'Other', length(unique(tumour_id)) ]

cn_events[, histology := clin[ match( tumour_id, tumour_id_muttable_cruk), Histology_per_tumour_id_muttable ] ]
cn_events[ grepl('adenocarcinoma', histology), histology := 'LUAD' ]
cn_events[ grepl('Squamous', histology), histology := 'LUSC' ]
cn_events[ !histology %in% c('LUAD', 'LUSC'), histology := 'Other' ]

driver_summary <- unique(cn_events[, .(total_drivers = .N, histology), 
                                   by = .(tumour_id, clonality = ifelse(is_clonal, 'Truncal', 'Subclonal'))])
driver_summary[ histology == 'LUAD', drivers_per_case := total_drivers / adeno_num ]
driver_summary[ histology == 'LUSC', drivers_per_case := total_drivers / squam_num ]
driver_summary[ histology == 'Other', drivers_per_case := total_drivers / other_num ]

# add missing
missing <- setdiff( cn_events[, apply(expand.grid(unique(paste(tumour_id, histology, sep = '--' )), c('Truncal', 'Subclonal')), 1, paste, collapse="--") ],
                    driver_summary[, paste( tumour_id, histology, clonality, sep = '--' )]) 
missing <- data.table( tumour_id = tstrsplit(missing, split = '--')[[1]],
                       clonality = tstrsplit(missing, split = '--')[[3]],
                       total_drivers = 0,
                       histology = tstrsplit(missing, split = '--')[[2]],
                       drivers_per_case = 0)
driver_summary <- rbindlist( list(driver_summary, missing) )


colours <- c('Truncal' = "#377EB8",
             'Subclonal' = "#E41A1C")

my_comparisons <- list( c("Truncal", "Subclonal") )

pdf(paste0(outputs.folder, date, "_overall_cn_drivers.pdf"), width = 6, useDingbats = FALSE)

ggplot(driver_summary[ !histology == 'Other' ], 
       aes( x = histology, y = total_drivers, fill = clonality)) + 
  geom_boxplot( position = "dodge", outlier.shape = NA ) +
  geom_point( position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0.3) ) +
  stat_compare_means( size = 3.5, method = 'wilcox.test') +
  scale_fill_manual( values = colours ) +
  scale_y_continuous( limits = c(-0.3,8)) +
  scale_x_discrete( labels = c('LUAD', 'LUSC')) +
  
  theme_classic() +
  labs( x = '',
        y = 'SCNA driver events / tumour',
        fill = '') +
  theme( text = element_text( size = 20 ),
         axis.text.x = element_text( angle = 45, hjust = 1 ),
         plot.margin = margin(0,0,0,2.5, 'cm'))

dev.off()


#############
#### END ####
#############

