######################################################################################################################
#############                Analysis of timing of selection and evolutionary dependencies               ############# 
######################################################################################################################
# written by Alexander Frankell (alexander.frankell@crick.ac.uk)

# Description:
# Script to create Figure 3 of the manuscript "The natural history of NSCLC in TRACERx"

#options and libraries
options(stringsAsFactors = F)
suppressWarnings( suppressPackageStartupMessages( library(fst) ) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(amfFunctions) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(ggrepel) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(ggpubr) )
suppressPackageStartupMessages( library(dndscv) )
suppressPackageStartupMessages( library(metap) )
suppressPackageStartupMessages( library(discover) )
suppressPackageStartupMessages( library(GeneAccord) )
suppressPackageStartupMessages( library(GenomicRanges) )

#setwd('/Volumes/proj-tracerx-lung/tctProjects/frankella/repos/shared_repos/TRACERx421_data_code/20220808_Tumour_evoutionary_histories/')

#parameters 
outputs.folder  <- './'
#outputs.folder  <- '../../../../Tx_exome/Plots/Fig3_suppig9a/20221130/'
date <- gsub("-","",Sys.Date())

##############################################
#### Get Inputs required for all analyses ####
##############################################


# these input file paths will already exist if speciified in following args, if not use defaults specified here #
dnds_histology_path         <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_dndscv_output_by_histology.rds"
dnds_clonal_illusion_path   <- "../20221109_Tumour_evo_histories_DATA/20221122_TRACERx421_dndscv_output_clonal_illusions.tsv"
mutTable_path               <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst"
clinical_data_path          <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_all_tumour_df.rds"
mut_driver_lung_path        <- "../20221109_Tumour_evo_histories_DATA/gene_anno_df.fst"
events_path                 <- "../20221109_Tumour_evo_histories_DATA/20221130_TRACERx421_gd_driver_signature_events.tsv" 
events_path_background      <- "../20221109_Tumour_evo_histories_DATA/20221130_TRACERx421_gd_driver_signature_and_background_events.tsv"
tree_path     <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_phylogenetic_trees.rds"

dnds_hist <- readRDS(  dnds_histology_path )

dnds_clonal_illusion <- fread( dnds_clonal_illusion_path )

### load mut table ###
mutTable <- fst::read_fst( mutTable_path, as.data.table = TRUE )

### load clinical data ###
clin <- as.data.table( readRDS( clinical_data_path) )

#load lung driver genes
mut_drivers <- fst::read.fst( mut_driver_lung_path,as.data.table = T )
mut_drivers <- mut_drivers[ (is_lung_mut_driver), gene_name ]

# Load gd, driver and signatures for each tumour
events_all <- data.table::fread( events_path )
events_all_back <- data.table::fread( events_path_background )

# load trees
trees <- readRDS( tree_path )

###################
#### Figure 3A ####
###################

mround <- function( x, base ) base * round( x / base )
mfloor <- function( x, base ) base * floor( x / base )
mceiling <- function( x, base ) base * ceiling( x / base )

plot_dnds_susinct <- function( dnds_plot, type, title, limits = c(0.5, 500), driver_excess_lim = 2, total_cases,
                               add_legend = TRUE){
  
  cols <- paste( type, c( 'w', 'w_low', 'w_high'), sep = '_' )
  
  dnds_plot[, clonal_subclonal_diff_w_OR := get(cols[1])[ type == 'clonal'] / get(cols[1])[ type == 'subclonal'], by = gene_name ]
  dnds_plot[, gene_name := factor( gene_name, levels = dnds_plot[ order(clonal_subclonal_diff_w_OR), unique(gene_name) ]) ]
  dnds_plot[, cat := 'subclonal and clonal' ]
  dnds_plot[ clonal_subclonal_diff_w_OR < 0.5, cat := 'subclonal favoured' ]
  dnds_plot[ clonal_subclonal_diff_w_OR > 2, cat := 'clonal favoured' ]
  
  muts_obs_col <- paste('total', type, sep = "_")
  dnds_plot[, dnds_limited := get(cols[1]) ]
  dnds_plot[ get(cols[1]) < 1, dnds_limited := 1 ]
  dnds_plot[, excess_mutations := get(muts_obs_col) * ((dnds_limited - 1)/dnds_limited)]
  dnds_plot[, excess_mutations := sum(excess_mutations), by = gene_name]
  dnds_plot[, cohort_perc := excess_mutations / total_cases * 100 ]
  dnds_plot[ cohort_perc > 100, cohort_perc := 100 ]
  dnds_plot <- dnds_plot[ excess_mutations > driver_excess_lim ]
  siggenes  <- dnds_plot[ qallsubs_cv < 0.05 | get(cols[2]) > 1, gene_name ]
  dnds_plot <- dnds_plot[ gene_name %in% siggenes ]
  
  dnds_plot[, (cols) := lapply(.SD, function(col){ col[ col < limits[1] ] <- limits[1] ; col }), .SDcols = cols]
  dnds_plot[, (cols) := lapply(.SD, function(col){ col[ col > limits[2] ] <- limits[2] ; col }), .SDcols = cols]
  
  dnds_plot_cl <- dnds_plot[ type == 'clonal', .(gene_name, cat, clonal_subclonal_diff_w_OR, cohort_perc,  get(cols[1]), get(cols[2]), get(cols[3]) )]
  dnds_plot_subcl <- dnds_plot[ type == 'subclonal', .( get(cols[1]), get(cols[2]), get(cols[3]) )]
  
  setnames(dnds_plot_cl, 5:7, paste0( cols, "_clonal"))
  setnames(dnds_plot_subcl, paste0( cols, "_subclonal"))
  
  dnds_to_plot <- cbind(dnds_plot_cl, dnds_plot_subcl)
  
  cols <- c( paste0(cols, "_clonal"), paste0(cols, "_subclonal") )
  
  colours <- c('clonal favoured' = "#377EB8",
               'subclonal and clonal' = "#984EA3",
               'subclonal favoured' = "#E41A1C")
  
  
  w_max_x_lim <- mceiling(dnds_plot_cl[, max( get(paste0( type, "_w_high_clonal" )) )], 50)
  w_min_x_lim <- mfloor(dnds_plot_cl[, min( get(paste0( type, "_w_low_clonal" )) )], 0.1)
  w_max_y_lim <- mceiling(dnds_plot_subcl[, max( get(paste0( type, "_w_high_subclonal" )) )], 50)
  w_min_y_lim <- mfloor(dnds_plot_subcl[, min( get(paste0( type, "_w_low_subclonal" )) )], 0.1)
  w_max_y_lim <- max(w_max_x_lim, w_max_y_lim)*0.8
  w_max_x_lim <- max(w_max_x_lim, w_max_y_lim)
  
  both_to_sub_ratio <- dnds_to_plot[, table(cat)['subclonal favoured']/sum(table(cat)[c('subclonal favoured', 'subclonal and clonal')])  ]
  xlim_labels <- as.numeric( ((log10(w_max_x_lim)-log10(w_min_x_lim))*both_to_sub_ratio)+log10(w_min_x_lim) )
  xlims <- list( c( NA, xlim_labels), c( xlim_labels, NA) )
  xlims <- list( c( log10(w_min_x_lim), xlim_labels), c( xlim_labels, log10(w_max_x_lim)) )
  #if( dnds_to_plot[, table(cat)['subclonal favoured'] < 5 ] ) xlims[[1]][1] <- 0
  
  plot <- ggplot(dnds_to_plot, aes( x = get(cols[1]), y = get(cols[4]), colour = cat )) +
    geom_point( aes( size = cohort_perc) ) +
    geom_hline( yintercept = 1, colour = 'black', alpha = 0.2 ) +
    geom_vline( xintercept = 1, colour = 'black', alpha = 0.2 ) +
    geom_errorbar(  aes( ymin = get(cols[5]), ymax = get(cols[6]) ), alpha = 0.2 ) +
    geom_errorbarh( aes( xmin = get(cols[2]), xmax = get(cols[3]) ), alpha = 0.2 ) +
    geom_point(  aes( size = cohort_perc ) ) +
    scale_colour_manual( values = colours ) +
    scale_y_continuous( trans = 'log10', limits = c(w_min_y_lim*0.5, w_max_y_lim),
                        breaks = 10^(ceiling(log10(w_min_y_lim)):floor(log10(w_max_y_lim))),
                        expand = c(0.1, 0)) +
    scale_x_continuous( trans = 'log10', limits = c(w_min_x_lim, w_max_x_lim),
                        breaks = 10^(ceiling(log10(w_min_x_lim)):floor(log10(w_max_x_lim))) ) +
    scale_size( limits = c(0.5, 100), 
                breaks = c(1, 2, 5, 10, 20, 50), range = c(1,10) ) +
    guides(colour = guide_legend("",override.aes = aes(label = "", size = 4)),
           size = guide_legend('% Cases') ) +
    labs( x =  "dN/dS Clonal",
          y = "dN/dS Subclonal",
          title = title) +
    theme_classic() +
    theme( text = element_text( size = 20 ) ,
           plot.title = element_text( hjust = 0.5, size = 35))
  
  if( dnds_to_plot[, any(cat == 'clonal favoured') ]){
    clonal_nudge <- dnds_to_plot[ cat == 'clonal favoured',  get(cols[4]) ]
    clonal_nudge[ clonal_nudge < 1.5 ] <- 1.5
    plot <- plot + geom_text_repel( data = dnds_to_plot[ cat == 'clonal favoured' ], direction    = "x", segment.alpha = 0.17,
                                    aes( label = gene_name ), box.padding = 1, force_pull = 1, angle = 90, segment.linetype = 2,
                                    size = 6,  nudge_y  =   - clonal_nudge )
  }
  if( dnds_to_plot[, any(cat == 'subclonal and clonal') ]){
    plot <- plot + geom_text_repel( data = dnds_to_plot[ cat == 'subclonal and clonal' ], force_pull = 0, xlim = xlims[[2]],
                                    aes( label = gene_name, angle = 90 ), box.padding = 1, direction    = "x",  segment.alpha = 0.17, segment.linetype = 2,
                                    size = 6, nudge_y  = w_max_y_lim )
  }
  if( dnds_to_plot[, any(cat == 'subclonal favoured' ) ]){
    plot <- plot + geom_text_repel( data = dnds_to_plot[ cat == 'subclonal favoured' ], force_pull = 0,  xlim = xlims[[1]],
                                    aes( label = gene_name, angle = 90 ), box.padding = 1, direction    = "x",  segment.alpha = 0.17, segment.linetype = 2,
                                    size = 6, nudge_y  = w_max_y_lim )
  }
  
  if(!add_legend){
    plot <- plot + theme(legend.position = "none") 
  }
  
  print( plot )
}

pdf(useDingbats = FALSE, paste0( outputs.folder, date, "_Adenocarcinoma_dnds_drivers_point_mutations_scatter.pdf"), width = 7.7, height = 8 )

plot_dnds_susinct(dnds_hist[[1]][ histology == 'LUAD' &
                                    element_type == "Gene" ],
                  type = 'point', title = "LUAD", total_cases = 248, add_legend = FALSE )

dev.off()

data <- dnds_hist[[1]][ histology == 'LUAD' &
                          element_type == "Gene" ]
sig_genes <- data[ point_w_low > 1 | qallsubs_cv < 0.05, unique(gene_name) ]
length(sig_genes)
ors <- data[ gene_name %in% sig_genes, 
            .(or = point_w[ type == 'subclonal' ] / point_w[ type == 'clonal' ]), 
            by = gene_name ]
ors[, sum(or < 1)/.N ]
data[ type == 'subclonal' & point_w_low > 1, .N]
data[, intersect( gene_name[type == 'subclonal' & point_w_low > 1], 
                  gene_name[type == 'clonal' & point_w_low < 1])]
# [1] "CNTNAP2"  "CTNND2"   "HIST1H1C" "KMT2D"    "PTEN"     "RUNX1"    "SMAD4"   

pdf(useDingbats = FALSE, paste0( outputs.folder, date, "_squamous_dnds_drivers_point_mutations_scatter.pdf"), width = 10, height = 8 )

plot_dnds_susinct(dnds_hist[[1]][ histology == "LUSC" &
                                    element_type == "Gene" ],
                  type = 'point', title = "LUSC", total_cases = 138)

dev.off()


data <- dnds_hist[[1]][ histology == 'LUSC' &
                          element_type == "Gene" ]
sig_genes <- data[ point_w_low > 1 | qallsubs_cv < 0.05, unique(gene_name) ]
length(sig_genes)
ors <- data[ gene_name %in% sig_genes, 
             .(or = point_w[ type == 'subclonal' ] / point_w[ type == 'clonal' ]), 
             by = gene_name ]
ors[, sum(or < 1)/.N ]
data[ type == 'subclonal' & point_w_low > 1, .N]
data[ type == 'subclonal' & point_w_low > 1, gene_name ]

data[, intersect( gene_name[type == 'subclonal' & point_w_low > 1], 
                  gene_name[type == 'clonal' & point_w_low < 1])]
#[1] "ATM"    "B2M"    "BCLAF1" "DUSP22" "SETD2" 


##################################
#### Supplementary figure 8A ####
##################################

pdf(useDingbats = FALSE, paste0( outputs.folder, date, "_Adenocarcinoma_dnds_pathways_point_mutations_scatter.pdf"), width = 10, height = 7 )

plot_dnds_susinct(dnds_hist[[1]][ histology == 'LUAD' &
                                    element_type == "Pathway" ],
                  type = 'point', title = "", total_cases = 248 )

dev.off()

pdf(useDingbats = FALSE, paste0( outputs.folder, date, "_squamous_dnds_pathways_point_mutations_scatter.pdf"), width = 10, height = 7 )

plot_dnds_susinct(dnds_hist[[1]][ histology == "LUSC" &
                                    element_type == "Pathway" ],
                  type = 'point', title = "", total_cases = 138 )

dev.off()



###################
#### Figure 3B ####
###################

run_discover <- function( events, freq_theshold = NA, tumour_id = 'tumour_id', 
                          event_id = 'gene', genes = NA ){
  
  if( is.na(freq_theshold) & all(is.na(genes)) ) stop( 'provide either genes or a freq theshold')
  
  events <- as.data.table(events)
  
  setnames(events, event_id, 'event_id')
  setnames(events, tumour_id, 'tumour_id')
  
  if( !is.na(freq_theshold)){
    # get only relatively frequent events
    num_tumours  <- events[, length(unique(tumour_id))]
    freq_drivers <- events[ (is_driver), .(freq = .N/num_tumours), by = event_id ][ freq > freq_theshold , event_id ]
    freq_drivers <- as.character(freq_drivers)
    
  }
  
  if( !all(is.na(genes)) ) freq_drivers <- genes
  
  empty_output <- data.table( gene1 = character(),
                              gene2= character(),
                              p.value = numeric(),
                              q.value= numeric(),
                              relationship = character(),
                              neg_log_q = numeric())
  if( length(freq_drivers) == 0 ) return( empty_output )
  
  # make patient level interaction table of events
  events[, tumour_id := factor(tumour_id) ]
  events[, event_id := as.character(event_id) ]
  
  patient_driver_matrix <- suppressMessages( suppressWarnings( dcast(events[ event_id %in% freq_drivers & (is_driver) ], 
                                                                     event_id ~ tumour_id, fun.aggregate = length, drop = FALSE ) ) )
  patient_passenger_matrix <- suppressMessages( suppressWarnings( dcast(events[ !event_id %in% freq_drivers ], 
                                                                        event_id ~ tumour_id, fun.aggregate = length, drop = FALSE ) ) )
  patient_matrix <- rbind( patient_driver_matrix, patient_passenger_matrix )
  genes <- patient_matrix[, event_id ]
  patient_matrix <- as.matrix(patient_matrix[, -1])
  rownames(patient_matrix) <- genes
  patient_matrix[ patient_matrix > 1] <- 1
  
  events <- discover.matrix(patient_matrix)
  
  subset <- rownames(patient_matrix) %in% freq_drivers
  names(subset)<- rownames(patient_matrix)
  
  result.mutex <- pairwise.discover.test(events[ subset, ], alternative = 'less') 
  result.co <- pairwise.discover.test(events[ subset, ], alternative = 'greater') 
  
  result.mutex <- suppressWarnings( as.data.table( cbind(melt(result.mutex$p.values), melt(result.mutex$q.values)[,3]) ) )
  names(result.mutex) <- c("gene1", "gene2", "p.value", "q.value")
  result.co <- suppressWarnings( as.data.table( cbind(melt(result.co$p.values), melt(result.co$q.values)[,3] ) ) )
  setnames(result.co, c("gene1", "gene2", "p.value", "q.value"))
  
  result.mutex[, relationship := 'mutual exclusivity' ]
  result.co[, relationship := 'co-occurence' ]
  
  results <- rbindlist( list(result.mutex, result.co) )
  
  # remove the pair duplicates or where the same gene is being compared (these will be NA)
  results <- results[ !is.na(p.value) ]
  
  results[, neg_log_q := -log10(q.value) ]
  
  return(results)
  
}

## some functions

# add in missing pairs
add_missing <- function(data, all_pairs){
  all_pairs_concat <- apply(all_pairs, 2, function(x) paste(x, collapse = '_'))
  missing_pairs <- all_pairs[, !all_pairs_concat %in% data[, gene_pair]]
  data_missing <- data.table( gene1 = missing_pairs[1,], 
                              gene2 = missing_pairs[2,],
                              p.value = NA, q.value= NA, relationship = NA,
                              neg_log_q = NA, clonality_subset = NA, 
                              histology_subset = NA)
  data_missing[, gene_pair := paste(gene1, gene2, sep = '_') ]
  data <- rbind(data, data_missing)
  gene_order <- unique(c(all_pairs[1, ], all_pairs[2, ]))
  data <- data[ order( match(gene2, gene_order)), ]
  data <- data[ order( match(gene1, gene_order) ), ]
  data[, gene1 := factor(gene1, ordered = TRUE, levels = gene_order )]
  data[, gene2 := factor(gene2, ordered = TRUE, levels = gene_order )]
  return(data)
}

plot_me_co <- function( data, all_pairs, invert = FALSE, margin = NA, text_size = 25 ){
  
  data <- add_missing(data, all_pairs)
  
  if( invert ){
    gene_order <- rev(unique(c(all_pairs[1, ], all_pairs[2, ])))
    data[, gene1 := factor(gene1, ordered = TRUE, levels = gene_order )]
    data[, gene2 := factor(gene2, ordered = TRUE, levels = gene_order )]
  }
  
  data[ neg_log_q > 4, neg_log_q := 4 ]
  
  colours <- brewer.pal(10, 'Paired')[c(7,9)]
  names(colours) <- c('mutual exclusivity', 'co-occurence')
  
  plot <- ggplot(data, aes(x = gene1, y = gene2, fill = neg_log_q, colour = relationship)) + 
    geom_tile( data = data[ is.na(p.value) ], size = 0.1 ) +
    geom_tile( data = data[ !is.na(p.value) ], size = 2 ) +
    labs( fill = '-log(q)',
          x = "",
          y = "",
          colour = '') + 
    scale_fill_distiller( palette = 'Blues', limits = c(0,4), breaks = 1:4,
                          labels = c(seq(1, 3, 1), '>4'), direction = 1, na.value = 'white' ) +
    scale_colour_manual( values = colours, na.value = "grey"  ) +
    theme_classic() 
  
  if( !invert ){
    plot <- plot + scale_x_discrete(position = "top")+
      theme( text = element_text( size = text_size),
             axis.text.x = element_text( angle = 45, hjust = 0))
  }
  
  if( invert ){
    plot <- plot + scale_y_discrete(position = "right")+
      theme( text = element_text( size = text_size),
             axis.text.x = element_text( angle = 45, hjust = 1))
  }
  
  if( !all(is.na(margin)) ){
    plot <- plot +
      theme(plot.margin = unit(margin,"cm") )
  }
  
  print(plot)
  
}

# Function to run discover accross various mutations msubsets (histology and clonality)

run_all_meco <- function( events, clin, name, outdir, freq_multiplier = 1 ){
  
  
  clin <- as.data.frame( clin )
  
  
  is_adeno <- events$tumour_id %in% clin[clin$Histology_per_tumour_id_muttable == "Invasive adenocarcinoma", "tumour_id_muttable_cruk"]
  is_squam <- events$tumour_id %in% clin[clin$Histology_per_tumour_id_muttable == "Squamous cell carcinoma", "tumour_id_muttable_cruk"]
  is_other <- !is_adeno & !is_squam
  
  # freq <- 0.20
  # freq_adeno <- freq
  # freq_squam <- freq
  # freq_other <- freq
  # freq_all <- freq
  freq_adeno <- 0.05 * freq_multiplier
  freq_squam <- 0.05 * freq_multiplier
  freq_other <- 0.05 * freq_multiplier
  freq_all <- 0.05 * freq_multiplier
  
  run_discover_subsets <- function(events, subset, subsetname, freq_theshold ){
    me_co_sub <- run_discover( events[ subset ], freq_theshold = freq_theshold)
    me_co_sub[, clonality_subset := 'all']
    me_co_sub_clonal <- run_discover( events[ subset & is_clonal == TRUE ], freq_theshold = freq_theshold )
    me_co_sub_clonal[, clonality_subset := 'clonal']
    me_co_sub_subclonal <- run_discover( events[ subset & is_clonal == FALSE ], freq_theshold = freq_theshold )
    me_co_sub_subclonal[, clonality_subset := 'subclonal']
    me_co_sub <- rbindlist( list(me_co_sub, me_co_sub_clonal, me_co_sub_subclonal) )
    me_co_sub[, histology_subset := subsetname ]
    return(me_co_sub)
  }
  
  me_co_adeno <- run_discover_subsets( events, is_adeno, 'adeno', freq_theshold = freq_adeno )
  me_co_squam <- run_discover_subsets( events, is_squam, 'squam', freq_theshold = freq_squam )
  me_co_other <- run_discover_subsets( events, is_other, 'other', freq_theshold = freq_other )
  
  me_co_all   <- run_discover_subsets( events, rep(TRUE,events[,.N]), 'all', freq_theshold = freq_all )
  
  me_co_tumour_level <- rbindlist( list(me_co_adeno, me_co_squam, me_co_other, me_co_all) )
  me_co_tumour_level[, gene_pair := paste( gene1, gene2, sep = '__' ) ]
  
  return(me_co_tumour_level)
  
}

# Separate out different event types
events_all_specific <- copy( events_all_back )
events_all_specific[, gene := paste( gene, consequence, sep = '_' ) ]
events_all_specific[, gene := gsub('_', ' ', gene ) ]
events_all_specific[, gene := gsub('non-truncating|truncating', 'Mut', gene ) ]
events_all_specific[, gene := gsub('signature', '', gene ) ]
events_all_specific[, gene := gsub('Homdel', 'Del', gene ) ]
events_all_specific[, gene := gsub(' $', '', gene ) ]
events_all_specific[, gene := gsub(' GD', '', gene ) ]

# Run discover on subsets
me_co <- run_all_meco( events_all_specific, clin, outdir = outputs.folder, name = 'gene_sigs_specific' )
#fwrite(me_co, paste0(outputs.folder, date, '_all_me_co_tests_by_type_discover.tsv'))

## Remove any association that there is no evidence for this in either histology
all_tumours_me_co <- me_co[ histology_subset == 'all' & clonality_subset %in% c('clonal', 'subclonal', 'all') ]
all_tumours_me_co_all <- all_tumours_me_co[ clonality_subset == 'all' & q.value < 0.1 ]
all_tumours_me_co_cl <- all_tumours_me_co[ clonality_subset == 'clonal' & q.value < 0.1 ]
all_tumours_me_co_subcl <- all_tumours_me_co[ clonality_subset == 'subclonal'  & q.value < 0.1 ]

all_sig_genes <- all_tumours_me_co_all[, unique(c(as.character(gene1), as.character(gene2))) ]
all_sig_genes_cl <- all_tumours_me_co_cl[, unique(c(as.character(gene1), as.character(gene2))) ]
all_sig_genes_subcl <- all_tumours_me_co_subcl[, unique(c(as.character(gene1), as.character(gene2))) ]
all_tumours_me_co_cl[, gene_pair_type := paste(gene_pair, relationship, sep = ':')]
all_tumours_me_co_subcl[, gene_pair_type := paste(gene_pair, relationship, sep = ':')]
all_tumours_me_co_all[, gene_pair_type := paste(gene_pair, relationship, sep = ':')]

events_all_specific$SampleID <- sapply( strsplit(events_all_specific$tumour_id, split = '_'), function(x) x[1])
is_adeno <- events_all_specific$tumour_id %in% clin[ grepl("adenocarcinoma", clin$Histology_per_tumour_id_muttable), tumour_id_muttable_cruk ]
is_squam <- events_all_specific$tumour_id %in% clin[ grepl("Squamous", clin$Histology_per_tumour_id_muttable), tumour_id_muttable_cruk ]

adeno_clonal_sig_genes <- run_discover( events_all_specific[ is_adeno & is_clonal == TRUE ], genes = all_sig_genes_cl )
adeno_subclonal_sig_genes <- run_discover( events_all_specific[ is_adeno & is_clonal == FALSE ], genes = all_sig_genes_subcl )
squam_clonal_sig_genes <- run_discover( events_all_specific[ is_squam & is_clonal == TRUE ], genes = all_sig_genes_cl )
squam_subclonal_sig_genes <- run_discover( events_all_specific[ is_squam & is_clonal == FALSE ], genes = all_sig_genes_subcl )
adeno_all_sig_genes <- run_discover( events_all_specific[ is_adeno ], genes = all_sig_genes )
squam_all_sig_genes <- run_discover( events_all_specific[ is_squam ], genes = all_sig_genes )

adeno_all_sig_genes[, hist := 'adeno']
adeno_clonal_sig_genes[, hist := 'adeno']
adeno_subclonal_sig_genes[, hist := 'adeno']
squam_all_sig_genes[, hist := 'squam']
squam_clonal_sig_genes[, hist := 'squam']
squam_subclonal_sig_genes[, hist := 'squam']

clonal_hist_specific <- rbind(adeno_clonal_sig_genes, squam_clonal_sig_genes)
clonal_hist_specific[, gene_pair := paste( gene1, gene2, sep = '__')]
clonal_hist_specific[, gene_pair_type := paste(gene_pair, relationship, sep = ':')]
subclonal_hist_specific <- rbind(adeno_subclonal_sig_genes, squam_subclonal_sig_genes)
subclonal_hist_specific[, gene_pair := paste( gene1, gene2, sep = '__')]
subclonal_hist_specific[, gene_pair_type := paste(gene_pair, relationship, sep = ':')]
all_hist_specific <- rbind(adeno_all_sig_genes, squam_all_sig_genes)
all_hist_specific[, gene_pair := paste( gene1, gene2, sep = '__')]
all_hist_specific[, gene_pair_type := paste(gene_pair, relationship, sep = ':')]

clonal_hist_specific <- clonal_hist_specific[ gene_pair_type %in% all_tumours_me_co_cl$gene_pair_type ][, q_restricted := p.adjust(p.value, n = .N, method = 'fdr'), by = hist ]
clonal_hist_specific <- clonal_hist_specific[, .(min_q = min(q_restricted)), by = gene_pair ]
subclonal_hist_specific <- subclonal_hist_specific[ gene_pair_type %in% all_tumours_me_co_subcl$gene_pair_type][, q_restricted := p.adjust(p.value, n = .N, method = 'fdr'), by = hist  ]
subclonal_hist_specific <- subclonal_hist_specific[, .(min_q = min(q_restricted)), by = gene_pair ]
all_hist_specific <- all_hist_specific[ gene_pair_type %in% all_tumours_me_co_all$gene_pair_type][, q_restricted := p.adjust(p.value, n = .N, method = 'fdr'), by = hist  ]
all_hist_specific <- all_hist_specific[, .(min_q = min(q_restricted)), by = gene_pair ]

all_tumours_me_co_cl[, min_hist_specific_q := clonal_hist_specific[ match(all_tumours_me_co_cl$gene_pair, clonal_hist_specific$gene_pair), min_q ] ]
all_tumours_me_co_subcl[, min_hist_specific_q := subclonal_hist_specific[ match(all_tumours_me_co_subcl$gene_pair, subclonal_hist_specific$gene_pair), min_q ] ]
all_tumours_me_co_all[, min_hist_specific_q := all_hist_specific[ match(all_tumours_me_co_all$gene_pair, all_hist_specific$gene_pair), min_q ]]
all_tumours_me_co_cl <- all_tumours_me_co_cl[ min_hist_specific_q < 0.1  ]
all_tumours_me_co_subcl <- all_tumours_me_co_subcl[ min_hist_specific_q < 0.1 ]
all_tumours_me_co_all <- all_tumours_me_co_all[ min_hist_specific_q < 0.1 ]

# remove TERC & SOX2 and CCNE1 and ATK2 - these are just because they are in the same genomic location
all_tumours_me_co_cl <- all_tumours_me_co_cl[ !(gene1 %in% c('TERC Amp', 'SOX2 Amp') &  gene2 %in% c('TERC Amp', 'SOX2 Amp')) ]
all_tumours_me_co_cl <- all_tumours_me_co_cl[ !(gene1 %in% c('CCNE1 Amp', 'AKT2 Amp') &  gene2 %in% c('CCNE1 Amp', 'AKT2 Amp')) ]
all_tumours_me_co_subcl <- all_tumours_me_co_subcl[ !(gene1 %in% c('TERC Amp', 'SOX2 Amp') &  gene2 %in% c('TERC Amp', 'SOX2 Amp')) ]
all_tumours_me_co_subcl <- all_tumours_me_co_subcl[ !(gene1 %in% c('CCNE1 Amp', 'AKT2 Amp') &  gene2 %in% c('CCNE1 Amp', 'AKT2 Amp')) ]
all_tumours_me_co_all <- all_tumours_me_co_all[ !(gene1 %in% c('TERC Amp', 'SOX2 Amp') &  gene2 %in% c('TERC Amp', 'SOX2 Amp')) ]
all_tumours_me_co_all <- all_tumours_me_co_all[ !(gene1 %in% c('CCNE1 Amp', 'AKT2 Amp') &  gene2 %in% c('CCNE1 Amp', 'AKT2 Amp')) ]

all_genes_cl <- all_tumours_me_co_cl[, unique(c(as.character(gene1), as.character(gene2)))]
all_genes_cl <- all_genes_cl[ order( !grepl('WGD', all_genes_cl), !grepl('SBS', all_genes_cl), all_genes_cl ) ]
all_pairs_cl <- combn(all_genes_cl, 2)
all_pairs_concat <- apply(all_pairs_cl, 2, function(x) paste(x, collapse = '__'))

#make sure all pairs are the right way around
all_tumours_me_co_cl[, gene_pair_rev := paste(gene2, gene1, sep = '__') ]
all_tumours_me_co_cl[ gene_pair_rev %in% all_pairs_concat, 
                      `:=`(gene1 = tstrsplit(gene_pair_rev, split = '__')[[1]], 
                           gene2 = tstrsplit(gene_pair_rev, split = '__')[[2]])]
all_tumours_me_co_cl[, gene_pair := paste( gene1, gene2, sep = '__' ) ]
all_tumours_me_co_cl[, gene_pair_rev := NULL ]

all_tumours_me_co_cl[, gene1 := factor(gene1, levels = all_genes_cl )]
all_tumours_me_co_cl[, gene2 := factor(gene2, levels = all_genes_cl )]


all_genes_subcl <- all_tumours_me_co_subcl[, unique(c(as.character(gene1), as.character(gene2)))]
all_genes_subcl <- all_genes_subcl[ order( !grepl('WGD', all_genes_subcl), !grepl('SBS', all_genes_subcl), all_genes_subcl ) ]
all_pairs_subcl <- combn(all_genes_subcl, 2)
all_pairs_concat <- apply(all_pairs_subcl, 2, function(x) paste(x, collapse = '__'))

#make sure all pairs are the right way around
all_tumours_me_co_subcl[, gene_pair_rev := paste(gene2, gene1, sep = '__') ]
all_tumours_me_co_subcl[ gene_pair_rev %in% all_pairs_concat, 
                         `:=`(gene1 = tstrsplit(gene_pair_rev, split = '__')[[1]], 
                              gene2 = tstrsplit(gene_pair_rev, split = '__')[[2]])]
all_tumours_me_co_subcl[, gene_pair := paste( gene1, gene2, sep = '__' ) ]
all_tumours_me_co_subcl[, gene_pair_rev := NULL ]

all_tumours_me_co_subcl[, gene1 := factor(gene1, levels = all_genes_subcl )]
all_tumours_me_co_subcl[, gene2 := factor(gene2, levels = all_genes_subcl )]


all_genes_all <- all_tumours_me_co_all[, unique(c(as.character(gene1), as.character(gene2)))]
all_genes_all <- all_genes_all[ order( !grepl('WGD', all_genes_all), !grepl('SBS', all_genes_all), all_genes_all ) ]
all_pairs_all <- combn(all_genes_all, 2)
all_pairs_concat <- apply(all_pairs_all, 2, function(x) paste(x, collapse = '__'))

#make sure all pairs are the right way around
all_tumours_me_co_all[, gene_pair_rev := paste(gene2, gene1, sep = '__') ]
all_tumours_me_co_all[ gene_pair_rev %in% all_pairs_concat, 
                       `:=`(gene1 = tstrsplit(gene_pair_rev, split = '__')[[1]], 
                            gene2 = tstrsplit(gene_pair_rev, split = '__')[[2]])]
all_tumours_me_co_all[, gene_pair := paste( gene1, gene2, sep = '__' ) ]
all_tumours_me_co_all[, gene_pair_rev := NULL ]

all_tumours_me_co_all[, gene1 := factor(gene1, levels = all_genes_all )]
all_tumours_me_co_all[, gene2 := factor(gene2, levels = all_genes_all )]

# Plot relationships

pdf( paste0(outputs.folder, date, '_clonal_me_all_hist_hist_specific_p_filtered.pdf'), width = 10.7 )
plot_me_co( data = all_tumours_me_co_cl[, 1:9], all_pairs = all_pairs_cl, text_size = 30 )
dev.off()

pdf( paste0(outputs.folder, date, '_subclonal_me_all_hist_hist_specific_p_filtered.pdf'), width = 9.5, height = 4.6 )
plot_me_co( all_tumours_me_co_subcl[, 1:9], all_pairs_subcl, invert = T, margin = c(0,0,0,3), text_size = 30 )
dev.off()


###################
#### Figure 3C ####
###################

ordering_test <- function( events, clonal_freq_min = 0.1, subclonal_freq_min = 0.1, 
                           event_id = 'gene', tumour_id = 'tumour_id' ){
  
  events <- as.data.table(events)
  
  setnames(events, event_id, 'event_id')
  setnames(events, tumour_id, 'tumour_id')
  
  # to have power & limit hypothesis only test gene interactions A -> B
  # when A is at a reasonably high clonal frequecy and B is at a reasonably high
  # subclonal frequency
  total_tumours <-  events[, length(unique(tumour_id)) ]
  event_freqs <- events[, .(clonal_freq = .SD[ (is_clonal), length(unique(tumour_id)) / total_tumours],
                            subclonal_freq =  .SD[ (!is_clonal), length(unique(tumour_id)) / total_tumours]), 
                        by = event_id ]
  A_acceptable_genes <-  event_freqs[ clonal_freq > clonal_freq_min, event_id ]
  B_acceptable_genes <-  event_freqs[ clonal_freq > subclonal_freq_min, event_id ]
  
  out <- rbindlist( lapply( A_acceptable_genes, function(eventA){
    
    eventA_clonal_tumours <- events[ event_id == eventA & (is_clonal), unique(tumour_id)]
    
    out_geneA <- rbindlist( lapply( B_acceptable_genes, function(eventB){
      
      eventB_clonal_tumours <- events[ event_id == eventB & (is_clonal), unique(tumour_id) ]
      
      After_A_B_subcl_tumours <- events[ tumour_id %in% eventA_clonal_tumours & event_id == eventB & (!is_clonal), unique(tumour_id) ]
      Not_After_A_B_subcl_tumours <- events[ !tumour_id %in% eventA_clonal_tumours & event_id == eventB & (!is_clonal), unique(tumour_id) ]
      After_A_B_not_subcl_tumours <- setdiff( events[ tumour_id %in% eventA_clonal_tumours, unique(tumour_id) ], After_A_B_subcl_tumours )
      Not_After_A_B_not_subcl_tumours <- setdiff( events[ !tumour_id %in% eventA_clonal_tumours, unique(tumour_id) ], Not_After_A_B_subcl_tumours )
      
      EventB_after_table <- matrix( c( length(After_A_B_subcl_tumours), 
                                       length(Not_After_A_B_subcl_tumours),
                                       length(After_A_B_not_subcl_tumours),
                                       length(Not_After_A_B_not_subcl_tumours)), byrow = TRUE, ncol = 2 )
      
      After_B_A_subcl_tumours <- events[ tumour_id %in% eventB_clonal_tumours & event_id == eventA & (!is_clonal), unique(tumour_id) ]
      Not_After_B_A_subcl_tumours <- events[ !tumour_id %in% eventB_clonal_tumours & event_id == eventA & (!is_clonal), unique(tumour_id) ]
      After_B_A_not_subcl_tumours <- setdiff( events[ tumour_id %in% eventB_clonal_tumours, unique(tumour_id) ], After_B_A_subcl_tumours )
      Not_After_B_A_not_subcl_tumours <- setdiff( events[ !tumour_id %in% eventB_clonal_tumours, unique(tumour_id) ], Not_After_B_A_subcl_tumours )
      
      EventA_after_table <- matrix( c( length(After_B_A_subcl_tumours), 
                                       length(After_B_A_not_subcl_tumours),
                                       length(Not_After_B_A_subcl_tumours),
                                       length(Not_After_B_A_not_subcl_tumours)), byrow = TRUE, ncol = 2 )
      
      EventA_after_test <- fisher.test( EventA_after_table )
      EventB_after_test <- fisher.test( EventB_after_table )
      
      out <- data.table( Event_before = c(eventB, eventA), 
                         Event_after = c(eventA, eventB),
                         p = c(EventA_after_test$p.value, EventB_after_test$p.value),
                         OR =  c(EventA_after_test$estimate, EventB_after_test$estimate),
                         Event_before_event_after_subcl = c(length(After_B_A_subcl_tumours), length(After_A_B_subcl_tumours)),
                         Event_before_event_after_not_subcl =  c(length(After_B_A_not_subcl_tumours), length(After_A_B_not_subcl_tumours)),
                         Event_not_before_event_after_subcl =  c(length(Not_After_B_A_subcl_tumours), length(Not_After_A_B_subcl_tumours)),
                         Event_not_before_event_after_not_subcl =  c(length(Not_After_B_A_not_subcl_tumours), length(Not_After_A_B_not_subcl_tumours)) )
      
      return(out)
      
    } ))
    
    return(out_geneA)
    
  } ))
  
  out <- out[ !duplicated(paste( Event_before, Event_after )) ]
  
  out[, q := p.adjust( p, .N, method = 'fdr') ]
  
  return(out)
  
}

# Overlay histology
events_all[, histology := clin[ match( tumour_id, gsub('^.{1}_', '', gsub('-Cl', '+Cl', tumour_id_muttable_cruk))), Histology_per_tumour_id_muttable ] ]
events_all[ grepl('adenocarinoma', histology), histology := 'Invasive adenocarcinoma' ]
events_all[ grepl('Squamous', histology), histology := 'Squamous cell carcinoma' ]
events_all[ !histology %in% c('Invasive adenocarcinoma', 'Squamous cell carcinoma'), histology := 'Other']

# combine APOBEC
events_all[ gene %in% c('SBS2', 'SBS13'), gene := 'SBS2 & 13']
# combine adjacent CN events (on the same arm)
events_all[ gene %in% c('TERC', 'SOX2'), gene := 'TERC & SOX2']
events_all[ gene %in% c('CCNE1', 'AKT2'), gene := 'CCNE1 & AKT2']


events_all[, gene_type := paste( gene, consequence, sep = '_' ) ]
events_all[, gene_type := gsub('_', ' ', gene_type ) ]
events_all[, gene_type := gsub('non-truncating|truncating', 'Mut', gene_type ) ]
events_all[, gene_type := gsub('signature', '', gene_type ) ]
events_all[, gene_type := gsub('Homdel', 'Del', gene_type ) ]
events_all[, gene_type := gsub(' $', '', gene_type ) ]
events_all[, gene_type := gsub(' GD', '', gene_type ) ]

results <- ordering_test(events = events_all[ (is_driver) ], 
                         event_id = 'gene_type',
                         clonal_freq_min = 0.1, 
                         subclonal_freq_min = 0.1)
results[, type := 'All']
results_adeno <- ordering_test(events_all[ (is_driver) & histology == 'Invasive adenocarcinoma' ], 
                               clonal_freq_min = 0.1, subclonal_freq_min = 0.1 , event_id = 'gene_type')
results_adeno[, type := 'LUAD']
results_squam <- ordering_test(events_all[ (is_driver) &  histology == 'Squamous cell carcinoma' ], 
                               clonal_freq_min = 0.1, subclonal_freq_min = 0.1, event_id = 'gene_type')
results_squam[, type := 'LUSC']

all_results <- rbind( results, results_adeno, results_squam )

all_results_sig <- all_results[ q < 0.1 ]

all_results_sig[, label := paste(Event_before, '->', Event_after)]

# plot just significant hits with numbers

plot_hits <- function( data ){
  
  data <- as.data.table(data)
  
  #rename to capture correct cols in melt
  setnames(data, 1:2, c('event_before', 'event_after') )
  
  # if see an effect with the same gene amplificied this is an artefact 
  # (can't amplify more by our definitions)
  data <- data[ !(event_before == event_before & grepl('Amp', event_before) ) ]
  
  #supress warning on different column classes - defaults to character corrected below
  data <- melt(data[ q < 0.1 ], 
               measure = patterns('^Event_'), 
               id.vars = c('q', 'OR', 'event_before', 'event_after'),
               value.name = 'count',
               variable.name = 'subset' ) 
  
  
  data[, Event_before_label := ifelse( grepl( 'Event_before', subset), paste( 'Truncal',event_before), paste( 'No truncal',event_before)) ]
  data[, Event_before := ifelse( grepl( 'Event_before', subset), 'Truncal present', 'Truncal absent') ]
  data[, ylab := paste('% subclonal', event_after)]
  
  data <- data[, .(percentage = count[ !grepl('not_subcl', subset) ] / ( count[ grepl('not_subcl', subset) ] + count[ !grepl('not_subcl', subset) ] ) ),
               by = .(q, OR, Event_before, Event_before_label, ylab, event_before, event_after)]
  
  
  data[, gene_pair := paste(event_before, event_after) ]
  #data[, label :=  paste('Subclonal', event_after, sep = '\n')]
  data[, label := gsub('& SOX2 ', '', event_after) ]
  
  data[, test_label :=  paste0( 'q =\n', round(q, 3))]
  
  data[, gene_pair := factor(gene_pair, levels = data[ order(q), unique(gene_pair)])]
  data[, Event_before_label := factor(Event_before_label, levels = data[ rev(order(grepl('No Truncal', Event_before_label))), unique(Event_before_label)])]
  
  # Amplificaiton a bit too long..
  data[, label :=  gsub('Amplification', 'Amp', label)]
  data[, label :=  gsub(' & SOX', ' &\nSOX', label)]
  data[, label :=  gsub('mutation', 'Mut', label)]
  
  labels = data[Event_before=='Truncal present', label]
  names(labels) = data[Event_before=='Truncal present', gene_pair]
  
  colours <- c(adjustcolor( "#1550D6", alpha.f = 0.1), '#1550D6')
  
  print(
    ggplot(data, aes( x = Event_before_label, y = percentage  )) +
      geom_bar(stat = 'identity', aes(fill = Event_before), colour = 'black') + 
      stat_compare_means() +
      scale_size( limits = c(0.5, 100), 
                  breaks = c(1, 2, 5, 10, 20, 50), range = c(1,10) ) +
      scale_fill_manual( values = colours ) +
      geom_text( aes(label = test_label, x = 1.5, y = 1.25), size = 4, lineheight = 0.85 ) +
      geom_segment(aes(x = 1, y = 1.035, xend = 2, yend = 1.035)) +
      geom_segment(aes(x = 1, y = 1.035, xend = 1, yend = 1.02)) +
      geom_segment(aes(x = 2, y = 1.035, xend = 2, yend = 1.02)) +
      facet_grid(~ gene_pair, scales = "free",
                 labeller = labeller(gene_pair = labels) ) +
      scale_y_continuous( breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20), limits = c(0,1.4) ) +
      labs( x = '',
            y = '% tumours with subclonal event') +
      theme_classic() +
      theme( text = element_text( size = 17),
             axis.text.x = element_text( angle = 45, hjust = 1),
             strip.text.x.top = element_text(size = 13, angle = 90),
             strip.background = element_rect( color="#D63515", fill=adjustcolor( "#D63515", alpha.f = 0.35), 
                                              size=0.7, linetype="solid"),
             legend.position = "none",
             plot.margin=unit(c(0,0.2,0,2), "cm"))
  )
  
}

pdf(paste0(outputs.folder, date, "_figure3C.pdf"), width = 14.3, height = 6)
plot_hits( data = all_results_sig[ type == 'All' ] )
invisible( dev.off() )


###################
#### Figure 3D ####
###################


mutTable[, histology := clin[ match( tumour_id, tumour_id_muttable_cruk ), Histology_per_tumour_id_muttable ] ]

run_dndscv <- function(mutTable, remove.hypermutants.cutoff.frac = NA, gene_list = NA, max_muts = 10000){

  mutTable <- mutTable[,c("tumour_id","chr","start","ref","var")]
  
  mutTable$chr <- gsub("chr","",mutTable$chr)
  
  names(mutTable) <- c("sampleID","chr","start","ref","mut")
  
  if(!is.na(remove.hypermutants.cutoff.frac)){
    hypermuts <- sapply(unique(mutTable$sampleID), function(sample) sum(mutTable$sampleID == sample))
    hypermuts <- unique(mutTable$sampleID)[hypermuts>quantile(hypermuts,remove.hypermutants.cutoff.frac)]
    mutTable <- mutTable[!mutTable$sampleID %in% hypermuts,]
  }
  
  if(all(is.na(gene_list))) return(dndscv(mutTable,
                                          outmats = T, 
                                          max_coding_muts_per_sample = max_muts))
  
  if(!all(is.na(gene_list))) return(dndscv(mutTable,
                                           gene_list=gene_list,
                                           outmats = T, 
                                           max_coding_muts_per_sample = max_muts))
}



mut_drivers_cd <- mut_drivers
mut_drivers_cd <- c(mut_drivers_cd, "CDKN2A.p14arf",   "CDKN2A.p16INK4a")
mut_drivers_cd <- mut_drivers_cd[ !mut_drivers_cd == 'CDKN2A' ]

LUAD_sel_cl <- run_dndscv( mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'C' ])
LUAD_sel_cl <- as.data.table( genesetdnds( LUAD_sel_cl, gene_list = mut_drivers_cd  )[[1]], keep.rownames = 'name' )
LUAD_sel_cl[, `:=`(histology = 'LUAD',
                   group = 'All truncal') ]

LUAD_sel_subcl <- run_dndscv( mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' ] )
LUAD_sel_subcl <- as.data.table( genesetdnds( LUAD_sel_subcl, gene_list = mut_drivers_cd  )[[1]], keep.rownames = 'name'  )
LUAD_sel_subcl[, `:=`(histology = 'LUAD',
                      group = 'All subclonal') ]

LUSC_sel_cl <- run_dndscv( mutTable[ grepl('Squamous', histology) & PyCloneClonal_SC == 'C'] )
LUSC_sel_cl <- as.data.table( genesetdnds( LUSC_sel_cl, gene_list = mut_drivers_cd  )[[1]], keep.rownames = 'name' )
LUSC_sel_cl[, `:=`(histology = 'LUSC',
                   group = 'All truncal') ]

LUSC_sel_subcl <- run_dndscv( mutTable[ grepl('Squamous', histology) & PyCloneClonal_SC == 'S'] )
LUSC_sel_subcl <- as.data.table( genesetdnds( LUSC_sel_subcl, gene_list = mut_drivers_cd  )[[1]],  keep.rownames = 'name' )
LUSC_sel_subcl[, `:=`(histology = 'LUSC',
                      group = 'All subclonal') ]

dnds <- rbindlist( list(LUAD_sel_cl,  LUAD_sel_subcl,
                        LUSC_sel_cl, LUSC_sel_subcl) )
dnds[, name := gsub('_geneset$', '', name )]

dnds_plot <- dnds[ name == 'wall' ] 
dnds_clonal_illusion[, group := is_CI ]
dnds_clonal_illusion[, is_CI := NULL ]
setnames(dnds_clonal_illusion, 'Histology', 'histology')
dnds_clonal_illusion[, total_subcl_muts := sum(unique(no.mutations)), by = .(histology)]
dnds_clonal_illusion[, perc.muts := paste0(round(no.mutations *100 / total_subcl_muts,0), '%') ]
dnds_clonal_illusion[, group := gsub('Non-', 'No ', group) ]


dnds_plot <- merge(dnds_plot, dnds_clonal_illusion[ name == 'wall_geneset' ], all = TRUE) 

dnds_plot[ , group := factor(group, levels = c('All truncal', 'All subclonal', 'Clonal illusion','No clonal illusion'))]

dnds_plot[, clonality := 'Subclonal']
dnds_plot[ group == 'All truncal', clonality := 'Truncal']
dnds_plot[ is.na(perc.muts), perc.muts := '' ]

# value a for expansion lines
min_y = 3.7
max_y = 4
min_x = 2
x2 <- min_x + 1
max_x = 4
mid_x = ((max_x - 1) / 2) + min_x
mid_y = ((max_y - min_y) / 2 ) + min_y

## add in the dnds values from different subclones from Kristiana

colours <- c('Truncal' = "#377EB8",
             'Subclonal' = "#E41A1C")

pdf(useDingbats = FALSE, paste0( outputs.folder, date, "_overall_selection_adeno_vs_squam.pdf"), width = 5, height = 6 )

ggplot( dnds_plot ) +
  geom_pointrange( aes( x = group, y = mle, ymin = cilow, 
                        ymax = cihigh, colour = clonality ),
                   size = 1.2) +
  scale_color_manual( values = colours ) +
  geom_hline( yintercept = 1 ) +
  scale_y_continuous(  limits = c(0.1, 4) ) +
  labs( x = "",
        y = paste0("dN/dS Lung driver genes"),
        colour = '' ) +
  geom_text( aes( label = perc.muts, x = group, y = 3.555 ), size = 4.5) +
  geom_segment(aes(x = min_x, y = min_y, xend = min_x, yend = max_y)) +
  geom_segment(aes(x = min_x, y = max_y, xend = mid_x, yend = max_y)) +
  geom_segment(aes(x = mid_x, y = max_y, xend = mid_x, yend = mid_y)) +
  geom_segment(aes(x = max_x, y = mid_y, xend = x2, yend = mid_y)) +
  geom_segment(aes(x = max_x, y = mid_y, xend = max_x, yend = mid_y)) +
  geom_segment(aes(x = x2, y = mid_y, xend = x2, yend = min_y)) +
  geom_segment(aes(x = max_x, y = mid_y, xend = max_x, yend = min_y)) +
  theme_classic() +
  theme( text = element_text( size = 20 ),
         axis.text.x = element_text( angle = 45, hjust = 1),
         legend.title = element_blank()) +
  facet_wrap( ~ histology )


dev.off()

#####################################
#### Check for stat significance ####
#####################################

##### Function from sanger
### https://zenodo.org/record/3966023#.YanjS_HP2cZ
variable_dNdS_twodatasets_overall = function(dnds1, dnds2, genestotest) {
  
  # We can implement a simple LRT model based on the uniform dNdS model
  # This is different from the Fisher test in that it uses synonymous mutations (i.e. dN/dS ratios)
  # instead of comparing the contribution of nonsyn muts of a gene *relative* to other genes.
  # Being a uniform model it assumes no considerable changes in the mutation rate variation or coverage
  # across genes in both datasets. But takes into account signature and rate variation between two
  # datasets.
  # H0: wmis1==wmis2 & wtru1==wtru2
  # H1: wmis1!=wmis2 & wtru1!=wtru2
  # This is simply done using obs1, exp1, obs2, exp2 (y1 and y2 vectors below)
  
  y1 = as.numeric(c(NA, colSums( dnds1$genemuts[ dnds1$genemuts$gene_name %in% genestotest, 2:9 ] ) ))
  y2 = as.numeric(c(NA, colSums( dnds2$genemuts[ dnds2$genemuts$gene_name %in% genestotest, 2:9 ] ) ))
  
  # Global dN/dS ratios from all other genes (to normalise the differences for the gene being tested)        
  ind1 = !dnds1$genemuts$gene %in% genestotest
  ind2 = !dnds2$genemuts$gene %in% genestotest
  wmis1_global = sum(dnds1$genemuts$n_mis[ind1])/sum(dnds1$genemuts$exp_mis[ind1])
  wmis2_global = sum(dnds2$genemuts$n_mis[ind2])/sum(dnds2$genemuts$exp_mis[ind2])
  wtru1_global = sum(dnds1$genemuts$n_non[ind1]+dnds1$genemuts$n_spl[ind1])/sum(dnds1$genemuts$exp_non[ind1]+dnds1$genemuts$exp_spl[ind1])
  wtru2_global = sum(dnds2$genemuts$n_non[ind2]+dnds2$genemuts$n_spl[ind2])/sum(dnds2$genemuts$exp_non[ind2]+dnds2$genemuts$exp_spl[ind2])
  
  # MLE dN/dS ratios using the uniform model under H0 and H1
  wmis_mle0 = (y1[3]+y2[3])/(y1[7]*wmis1_global+y2[7]*wmis2_global)
  wtru_mle0 = sum(y1[4:5]+y2[4:5])/sum(y1[8:9]*wtru1_global+y2[8:9]*wtru2_global)
  wmis_mle1 = c(y1[3],y2[3])/c(y1[7]*wmis1_global,y2[7]*wmis2_global)
  wtru_mle1 = c(sum(y1[4:5]),sum(y2[4:5]))/c(sum(y1[8:9]*wtru1_global),sum(y2[8:9]*wtru2_global))
  
  # Observed and predicted counts under H0 and H1
  obs = as.numeric(c(y1[3], sum(y1[4:5]), y2[3], sum(y2[4:5])))
  exp0 = as.numeric(c(y1[7]*wmis1_global*wmis_mle0, sum(y1[8:9])*wtru1_global*wtru_mle0, y2[7]*wmis2_global*wmis_mle0, sum(y2[8:9])*wtru2_global*wtru_mle0))
  exp1 = as.numeric(c(y1[7]*wmis1_global*wmis_mle1[1], sum(y1[8:9])*wtru1_global*wtru_mle1[1], y2[7]*wmis2_global*wmis_mle1[2], sum(y2[8:9])*wtru2_global*wtru_mle1[2])) # Note that exp1 == obs (we only have this line here for confirmation purposes)
  ll0 = c(sum(dpois(x=obs[c(1,3)], lambda=exp0[c(1,3)], log=T)), sum(dpois(x=obs[c(2,4)], lambda=exp0[c(2,4)], log=T)))
  ll1 = c(sum(dpois(x=obs[c(1,3)], lambda=exp1[c(1,3)], log=T)), sum(dpois(x=obs[c(2,4)], lambda=exp1[c(2,4)], log=T)))
  
  # One-sided p-values
  pvals = (1-pchisq(2*(ll1-ll0), df=1))
  if (wmis_mle1[1]<wmis_mle1[2]) { pvals[1] = 1 } else { pvals[1] = pvals[1]/2 }
  if (wtru_mle1[1]<wtru_mle1[2]) { pvals[2] = 1 } else { pvals[2] = pvals[2]/2 }
  
  # Saving the results
  pvec = 1 - pchisq(-2 * sum(log(pvals)), df = 4) # Fisher combined p-value
  rmisvec = wmis_mle1[1]/wmis_mle1[2]
  rtruvec = wtru_mle1[1]/wtru_mle1[2]
  
  return( data.table(p = pvec, rmis = rmisvec, rtru =rtruvec) )

}

## Which mutations have clonal illusion
for( tumour in mutTable[, unique(tumour_id)[ unique(tumour_id) %in% names(trees) ]]){
  clone_table <- trees[[tumour]]$clonality_out$clonality_table
  is_LN <- grepl('LN', names(clone_table))
  clonal_illusion_prim_clones <- rownames(clone_table)[ apply(clone_table, 1, function(row) any(row == 'clonal' & !is_LN)) ]
  is_LN_unqiue <- rownames(clone_table)[ apply(clone_table, 1, function(row) !any(!row == 'absent' & !is_LN)) ]
  
  tree_clones <- rownames(clone_table)
  mutTable[ tumour_id == tumour, is_tree_clone := PyCloneCluster_SC %in% tree_clones ]
  mutTable[ tumour_id == tumour, clonal_illusion_prim := PyCloneCluster_SC %in% clonal_illusion_prim_clones ]
  mutTable[ tumour_id == tumour, is_LN_unqiue := PyCloneCluster_SC %in% is_LN_unqiue ]
  
}

dnds_CI <- run_dndscv( mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' & (is_tree_clone & clonal_illusion_prim & !is_LN_unqiue) ])
dnds_NCI <- run_dndscv( mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' & (is_tree_clone & !clonal_illusion_prim & !is_LN_unqiue) ])

dnds1 <- dnds_CI
dnds2 <- dnds_NCI


dnds_CIg <- run_dndscv( mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' & (is_tree_clone & clonal_illusion_prim & !is_LN_unqiue) ], gene_list = mut_drivers_cd)
dnds_NCIg <- run_dndscv( mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' & (is_tree_clone & !clonal_illusion_prim & !is_LN_unqiue) ], gene_list = mut_drivers_cd)
dnds_CIg[[1]]
dnds_NCIg[[1]]
mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' & (is_tree_clone & clonal_illusion_prim & !is_LN_unqiue), .N ]
mutTable[ grepl('adenocarcinoma', histology) & PyCloneClonal_SC == 'S' & (is_tree_clone & !clonal_illusion_prim & !is_LN_unqiue), .N ]
out <- variable_dNdS_twodatasets_overall(dnds_CI, dnds_NCI, mut_drivers_cd)



#############
#### END ####
#############

