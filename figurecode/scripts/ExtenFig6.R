######################################################################################################################
#############                Analysis of timing of selection and evolutionary dependencies               ############# 
######################################################################################################################
# written by Alexander Frankell (alexander.frankell@crick.ac.uk)

# Description:
# Script to create Extended Figure 6 of the manuscript "The natural history of NSCLC in TRACERx"
#options and libraries
options(stringsAsFactors = F)
suppressWarnings( suppressPackageStartupMessages( library(fst) ) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(amfFunctions) ) #www.github.com/amf71/amfFunctions
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(ParrallelGDDetect) )
# setwd('/Volumes/proj-tracerx-lung/tctProjects/frankella/repos/shared_repos/TRACERx421_data_code/20220407_Tumour_evoutionary_histories/')
source('../../../my_R_packages/ParrallelGDDetect/R/parrallel_gd_caller.R')

#parameters
outputs.folder  <- '.'
# outputs.folder  <- '../../../../Tx_exome/Plots/ExtendedFig6/'
date <- gsub("-","",Sys.Date())

##############################################
#### Get Inputs required for all analyses ####
##############################################

# these input file paths will already exist if speciified in following args, if not use defaults specified here #
mutTable_path         <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_mutation_table.fst"
mutTable_region_path  <- "../20221109_Tumour_evo_histories_DATA/20221123_TRACERx421_mutation_table_region.fst"
gd_data_path          <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_gd_table_raw.tsv"
purity_data_path      <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_manual_qc_sheet.tsv"
TRACERx_sims_path      <- "../20221109_Tumour_evo_histories_DATA/20221109_TRACERx421_500_simulated_TRACERx_tumours.fst"
TRACERx_sim_trees_path      <- "../20221109_Tumour_evo_histories_DATA/20220427_TRACERx421_simulated_trees.rds"
TRACERx_sim_trees_table_path      <- "../20221109_Tumour_evo_histories_DATA/20220427_TRACERx421_simulated_trees.table.tsv"
TRACERx_sim_trees_table_path      <- "../20221109_Tumour_evo_histories_DATA/20220427_TRACERx421_simulated_trees.table.tsv"
evo_metrics_path   <- "../20221109_Tumour_evo_histories_DATA/20221110_TRACERx421_evolutionary_metrics.tsv"

### load gd data raw ###
gd_data <- fread( gd_data_path )

### load purity qc table ###
purity_data <- fread( purity_data_path )
purity_data[, Sample := gsub('\\.', '-', Sample) ]

### load mut table ###
mutTable <- fst::read_fst( mutTable_path, as.data.table = TRUE )
mutTable_region <- fst::read_fst( mutTable_region_path, as.data.table = TRUE )

### load TRACERx simulations ###
tracerx_sims <- fst::read_fst( TRACERx_sims_path, as.data.table = TRUE )
tracerx_sim_trees <- readRDS( TRACERx_sim_trees_path )
tracerx_sim_trees_table <- fread( TRACERx_sim_trees_table_path )

## load evo metrics
evo_metrics <- fread( evo_metrics_path )

#######################
#### Extended Figure 6a & b ####
#######################

# adjust the nMaj3 cut off to 50%
gd_data[, second_gd := genome_frac_nMaj3 > 0.5 ]

# overlay purity info
gd_data[, purity := purity_data[ match( sample, Sample), ASCAT.purity ] ]
gd_data[, ploidy := purity_data[ match( sample, Sample), ASCAT_psi ] ]
reruns <- purity_data[ grepl("RERUN", ManualCall) | is.na(ASCAT.purity), Sample ]
gd_data[ sample %in% reruns, purity := purity_data[ match( sample, Sample ), Adjust_ASCAT.purity ] ]
gd_data[ sample %in% reruns, ploidy := purity_data[ match( sample, Sample ), Adjust_ASCAT.ploidy ] ]

#remove any CN_FAILS
cn_fails <- purity_data[ grepl("CN_FAIL|ALL_FAIL", ManualCall), Sample ]
gd_data <- gd_data[ !sample %in% cn_fails ]

# create col of gd status
gd_data[, gd_status := 0 ]
gd_data[ first_gd == TRUE, gd_status := 1 ]
gd_data[ second_gd == TRUE, gd_status := 2 ]

# recapitulate ASBOLUTE paper plot 

pdf( paste0( outputs.folder, date, "_ploidy_vs_mean_allelic_diff_0.5_nMaj3.pdf"), width = 8 )

ggplot( gd_data, aes( x = ploidy, y = mean_allelic_difference ) ) +
  geom_point( shape = 21, aes( colour = as.factor(gd_status) ), size = 2.5, stroke = 1.2 ) +
  scale_color_brewer( palette = "Set1" )+
  scale_y_continuous( limits = c( 0, 3.5 )) +
  scale_x_continuous( limits = c( 1, 7 )) +
  theme_classic() +
  theme( text = element_text( size = 15 ) )

invisible( dev.off() )


tum_ids <- mutTable[ ! is.na(MinorCPN_SC) ][ ! duplicated(tumour_id), .(tumour_id, MinorCPN_SC)]
tum_ids <- rbindlist(lapply(tum_ids[,unique(tumour_id)], function(tum){
  regions <- tum_ids[ tumour_id == tum, MinorCPN_SC ]
  regions <- sapply(strsplit(regions, split = ';')[[1]], function(x) gsub(':.*$', '', x))
  return( data.table( tumour_id = tum,
                      SampleID = gsub('_Cl.*$', '', tum),
                      regions = gsub('\\.','-', regions ) ) )
}))
tum_ids[, sample_full := paste(SampleID, regions, sep = '_')]
gd_data[, tumour := tum_ids[ match(sample, sample_full), tumour_id ]]

gd_data[, subclonal_gd := !length( unique( gd_status ) ) == 1, by = tumour ] 


# clonal WGD we'd be mroe confident in

# most useful for 1 vs 2 gd seems to be mean allelic diff, ploidy and nMaj3 look at them all

pdf( paste0( outputs.folder, date, "_nMaj3_vs_ploidy_vs_mean_allelic_diff.pdf") )

ggplot( gd_data, aes( x = genome_frac_nMaj3, y = mean_allelic_difference ) ) +
  geom_point( shape = 21,aes( colour = subclonal_gd, fill = ploidy, size = purity ), stroke = 1.5) +
  scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 4 ) +
  scale_colour_manual( values = c("black", "grey")) +
  scale_size( range = c(0.5, 6) ) +
  geom_vline( xintercept = 0.5,  linetype="dotted" ) +
  ylab("Mean allelic difference") +
  xlab("Genome fraction with major allele >= 3") +
  theme_classic() +
  theme( text = element_text( size = 20 ) )

invisible( dev.off() )

# most useful for 0 vs 1 gd seems to be mean allelic diff, ploidy and nMaj2 look at them all

pdf( paste0( outputs.folder, date, "_nMaj2_vs_ploidy_vs_mean_allelic_diff.pdf") )

ggplot( gd_data, aes( x = genome_frac_nMaj2, y = mean_allelic_difference ) ) +
  geom_point( shape = 21,aes( colour = subclonal_gd, fill = ploidy, size = purity ), stroke = 1.5) +
  scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 4 ) +
  scale_colour_manual( values = c("black", "grey")) +
  geom_vline( xintercept = 0.5,  linetype="dotted" ) +
  scale_size( range = c(0.5, 6) ) +
  ylab("Mean allelic difference") +
  xlab("Genome fraction with major allele >= 2") +
  theme_classic() +
  theme( text = element_text( size = 20 ) )

invisible( dev.off() )


#######################
#### Extended Figure 6d & e ####
#######################


# callect together data as ParallelGDDetect input
# merge data
gd_data[, num_gds := sum(first_gd, second_gd), 1:nrow(gd_data) ]
gd_data[, patient_id := gsub('^.{1}_', '', gsub('_SU.*$', '', sample)) ]
gd_data[, sample_short := gsub('^.{9}', '', sample) ]
mutTable_region[, patient_id := tstrsplit(tumour_id, split = '_Tum')[[1]] ]
mutTable_region[, tumour_id_short := paste0(patient_id, ifelse(grepl('Cluster', tumour_id), paste0('_', Tumour.cluster), '')) ]
mutTable_region[, patient_id := tstrsplit(tumour_id_short, split ='_Cl')[[1]] ]
mutTable_region[, region := gsub('\\.', '-', tstrsplit(RegionID, split = ':')[[2]]) ]
mutTable_region[, num_gds := gd_data[ match(paste(mutTable_region$patient_id, mutTable_region$region),
                                            paste(gd_data$patient_id, gd_data$sample_short)), num_gds] ]

#get long format Major/Minor CN
mutTable_region[, MajorCPN_SC_long := str_extract(MajorCPN_SC, paste0('(?<=', tstrsplit(RegionID, split = ':')[[2]], '.{1}?).{1}') ) ]
mutTable_region[, MinorCPN_SC_long := str_extract(MinorCPN_SC, paste0('(?<=', tstrsplit(RegionID, split = ':')[[2]], '.{1}?).{1}') ) ]
input_all <- copy( mutTable_region[ cleanCluster_SC == TRUE,
                                    .(tumour_id = tumour_id_short,
                                      cluster_id = PyCloneCluster_SC,
                                      is_clonal_cluster = PyCloneClonal_SC == 'C',
                                      sample_id = gsub('\\.', '-', region),
                                      CN_A = MajorCPN_SC_long, CN_B = MinorCPN_SC_long,
                                      mut_cpn = mut_cpn_SC, num_gds) ] )
input_all[, MajCN := max(CN_A, CN_B), 1:nrow(input_all) ]
data <- input_all[ !is.na(num_gds) ]




## benchmark with clonal mutations in gd/non-gd regions
data[ (is_clonal_cluster & num_gds == 0), benchmark := 'Clonal cluster 0 gd region' ]
data[ (is_clonal_cluster & num_gds == 1), benchmark := 'Clonal cluster 1 gd region' ]
data[ (is_clonal_cluster & num_gds == 2), benchmark := 'Clonal cluster 2 gd region' ]

cluster_metrics <- unique( data[, .(frac_cpn_2 = sum(mut_cpn > 1.5)/.N,
                                    ncpn_2 =  sum(mut_cpn > 1.5),
                                    frac_cpn_2_no_losses = sum(mut_cpn > 1.5 &  MajCN == 2^num_gds)/sum(MajCN == 2^num_gds),
                                    ncpn_2_no_losses = sum(mut_cpn > 1.5 &  MajCN == 2^num_gds),
                                    nmut_no_losses = sum(MajCN == 2^num_gds),
                                    nmut = .N,
                                    num_gds, 
                                    is_clonal_cluster,
                                    is_present = mean(mut_cpn) > 0.01,
                                    mean_mut_cpn = mean(mut_cpn),
                                    benchmark),
                                by = .(cluster_id, tumour_id, sample_id)])

pdf( paste0(outputs.folder, date, '_benchmark_all_mutations.pdf'), width = 9 )

ggplot( cluster_metrics[!is.na(benchmark)], aes(x = frac_cpn_2)) +
  geom_histogram( aes( fill = benchmark ), position = 'stack', colour = 'grey', bins = 50 ) +
  scale_fill_brewer( palette = 'Set1') +
  labs( x = 'fraction of mutations mut_cpn > 1.5') +
  theme_classic() +
  theme( text = element_text( size = 20 ),
         axis.text.x = element_text( size = 20 ))

dev.off()

pdf( paste0(outputs.folder, date, '_benchmark_majCN_filter.pdf'), width = 9 )

ggplot( cluster_metrics[!is.na(benchmark)], aes(x = frac_cpn_2_no_losses)) +
  geom_histogram( aes( fill = benchmark ), position = 'stack', colour = 'grey', bins = 50 ) +
  scale_fill_brewer( palette = 'Set1') +
  labs( x = 'fraction of mutations with MajCN = 2^number gds with mut_cpn > 1.5') +
  theme_classic() +
  theme( text = element_text( size = 20 ),
         axis.title.x = element_text( size = 15 ))

dev.off()


#######################
#### Exntended igure 6f ####
#######################

data <- copy( tracerx_sims )
truth_edges <- tracerx_sim_trees_table
trees <- tracerx_sim_trees
setnames(data,
         c('PATIENT', 'SAMPLE', 'MUT_CN', 'CHR', 'POS', "REF", 'ALT', 'CLUSTER'),
         c('tumour_id', 'sample_id', 'multiplicity', 'chromosome', 'position',
           'ref', 'alt', 'cluster_id'))

# Get average mutant copy number which would be observed from the vaf (ie with the same noise)
data[, vaf := VAR_COUNT / DEPTH ]
data[, mut_cpn := (vaf/PURITY)*(PURITY*(CN_A+CN_B)+2*(1-PURITY)) ]

# determine clonal cluster
data[, `:=`(mean_cluster_ccf_sample = mean(TRUE_CCF),
            nmuts = .N),
     by = .(cluster_id, sample_id)]
data[, `:=`(mean_cluster_ccf_tumour = mean(mean_cluster_ccf_sample)),
     by = .(cluster_id, tumour_id)]
data[, `:=`(max_mean_cluster_ccf_tumour = max(mean_cluster_ccf_sample)),
     by = .(cluster_id, tumour_id)]
data[, is_clonal_cluster := mean_cluster_ccf_tumour == max_mean_cluster_ccf_tumour,
     by = tumour_id ]
data[, mut_id := paste(tumour_id, chromosome, position, ref, alt) ]
data[, mean_CNA := mean(CN_A), by = mut_id ]
data[, mean_CNB := mean(CN_B), by = mut_id ]
data[, Maj_CN := ifelse(mean_CNA > mean_CNB, CN_A, CN_B)]
data[, frac_majCN_2 := sum(  round(Maj_CN) >= 2) / .N, by = sample_id]
data[, frac_majCN_3 := sum(  round(Maj_CN) >= 3) / .N, by = sample_id]
# Plot these across simulations to determine thesholds
data_plot <- unique( data[, .(frac_majCN_2, frac_majCN_3, PLOIDY, PURITY,
                              allelic_imbalence = mean(abs(CN_A - CN_B)[(is_clonal_cluster)]),
                              perc_clonal= sum(is_clonal_cluster)/.N,
                              num_clonal= sum(is_clonal_cluster) ),
                          by = .(sample_id, tumour_id) ] )

##########################################
#### Extract ground truth from trees  ####
##########################################
detectabliity_CCF_threshold <- 0.75 # Threshold for detecting a GD event in reality

## Overlay all the parents for each node
extract_all_parents <- function( child, tree ){
  parents <- unique( c(tree[ tree[,2] %in% child, 1]) )
  repeat{
    all.parents <- unique( c( parents, tree[ tree[,2] %in% parents, 1] ) )
    if( all(all.parents %in% parents ) ) break
    parents <- all.parents
  }
  return(parents)
}
truth_edges[, all_parents := paste( extract_all_parents( Child, tree = trees[[ which(names(trees) == tumour_id)]]), collapse = ','),
            1:nrow(truth_edges) ]

#iflter down just to the GD clones
truth_GD_events <- truth_edges[ (WGD),
                                .(tumour_id, all_parents, gd_clone = Child,
                                  is_clonal) ]

# Work out for each GD how many GD events have gone before in the lineage
# ie is this a 1st GD (from 2-4) or a second GD (from 4-8) etc
truth_GD_events[, num_gd := sum(paste(tumour_id, strsplit(all_parents, split = ',')[[1]]) %in% truth_GD_events[, paste(tumour_id, gd_clone)]) + 1, 1:nrow(truth_GD_events)]

# if its not the first GD note the parental cluster which were GD'd before
truth_GD_events[, parent_gds := paste(strsplit(all_parents, split = ',')[[1]][ paste(tumour_id, strsplit(all_parents, split = ',')[[1]]) %in% truth_GD_events[, paste(tumour_id, gd_clone)]], collapse = ','), 1:nrow(truth_GD_events)]

### Determine which regions and clones have GDs and which would be detectable (ie in at least 50% of cells)
# filter down to region - clone level (averaging over mutaions)
cluster_ccfs_sample <- data[, .(ccf = mean(TRUE_CCF)),
                            by = .(tumour_id, cluster_id, sample_id) ]
cluster_ccfs_sample[, is_gd_cluster := paste(tumour_id, cluster_id) %in% truth_GD_events[, paste(tumour_id, gd_clone)] ]
cluster_ccfs_sample[, num_gd := truth_GD_events[ match(paste(cluster_ccfs_sample$tumour_id, cluster_id),
                                                       paste(truth_GD_events$tumour_id, gd_clone)), num_gd ] ]
gd_cluster_ccfs_sample <- cluster_ccfs_sample[ (is_gd_cluster) ]

# if you add up all the GDs in a region that are in parallel would you still detect a GD (eve if no single GD has >50% CCF)
gd_cluster_ccfs_sample[ , is_region_gd_1_detectable := sum(ccf[num_gd == 1]) >= detectabliity_CCF_threshold,
                        by = .(tumour_id, sample_id) ]
gd_cluster_ccfs_sample[ , is_region_gd_2_detectable := sum(ccf[num_gd == 2]) >= detectabliity_CCF_threshold,
                        by = .(tumour_id, sample_id) ]
gd_cluster_ccfs_sample[, is_any_cluster_1_gd_detectable_region := any(ccf >= detectabliity_CCF_threshold & num_gd ==1), by = .(tumour_id, sample_id)]
gd_cluster_ccfs_sample[, is_clone_region_1_gd_detectable := ccf >= detectabliity_CCF_threshold & num_gd ==1 ]
gd_cluster_ccfs_sample[, is_gd_1_detectable := any(is_clone_region_1_gd_detectable), by = .(tumour_id, cluster_id) ]
gd_cluster_ccfs_sample[, is_any_cluster_2_gd_detectable_region := any(ccf >= detectabliity_CCF_threshold & num_gd ==2), by = .(tumour_id, sample_id)]
gd_cluster_ccfs_sample[, is_clone_region_2_gd_detectable := ccf >= detectabliity_CCF_threshold & num_gd ==2 ]
gd_cluster_ccfs_sample[, is_gd_2_detectable := any(is_clone_region_2_gd_detectable), by = .(tumour_id, cluster_id) ]

# Add that one of the undetectable GDs are detectable based on the addition of several low CCF GDs in a single region
gd_cluster_ccfs_sample[!is_any_cluster_1_gd_detectable_region & is_region_gd_1_detectable & !is_gd_1_detectable & num_gd == 1,
                       is_clone_region_1_gd_detectable := c(TRUE, rep(FALSE, .N-1)),
                       .(tumour_id, sample_id)]
gd_cluster_ccfs_sample[!is_any_cluster_2_gd_detectable_region & is_region_gd_2_detectable & !is_gd_2_detectable & num_gd == 2,
                       is_clone_region_2_gd_detectable := c(TRUE, rep(FALSE, .N-1)),
                       .(tumour_id, sample_id)]
gd_cluster_ccfs_sample[, is_clone_region_gd_detectable := is_clone_region_1_gd_detectable | is_clone_region_2_gd_detectable ]
gd_cluster_ccfs <- gd_cluster_ccfs_sample[, .(mean_ccf =  mean(ccf),
                                              is_clone_gd_detectable = any(is_clone_region_gd_detectable),
                                              is_gd_obs_clonal =  all(is_clone_region_gd_detectable)),
                                          by = .(tumour_id, cluster_id) ]

# add onto the tree info whether the GD events would be detectable based on this above
truth_GD_events[, is_gd_detectable := paste(tumour_id, gd_clone) %in% gd_cluster_ccfs[ (is_clone_gd_detectable), paste(tumour_id, cluster_id)] ]
truth_GD_events[, is_gd_obs_clonal := paste(tumour_id, gd_clone) %in% gd_cluster_ccfs[ (is_gd_obs_clonal), paste(tumour_id, cluster_id)] ]
truth_GD_events[, is_detectable_any_region := any(is_gd_detectable), .(tumour_id, gd_clone)]
truth_tumour <-  truth_GD_events[, .(num_clonal_gds = sum(is_clonal), num_subclonal_gds =  sum(!is_clonal),
                                     num_clonal_1_gds = sum(is_clonal & num_gd == 1), num_subclonal_1_gds =  sum(!is_clonal & num_gd == 1),
                                     num_clonal_2_gds = sum(is_clonal & num_gd == 2), num_subclonal_2_gds =  sum(!is_clonal & num_gd == 2),
                                     num_detectable_subclonal_gds = sum(!is_gd_obs_clonal & is_gd_detectable), num_detectable_clonal_gds =  sum(is_gd_obs_clonal & is_gd_detectable),
                                     num_detectable_subclonal_1_gds = sum(!is_gd_obs_clonal & is_gd_detectable & num_gd == 1), num_detectable_clonal_1_gds =  sum(is_gd_obs_clonal & is_gd_detectable & num_gd == 1),
                                     num_detectable_subclonal_2_gds = sum(!is_gd_obs_clonal & is_gd_detectable & num_gd == 2), num_detectable_clonal_2_gds =  sum(is_gd_obs_clonal & is_gd_detectable & num_gd == 2)),
                                 by = tumour_id ]

### Use ground truth number of detectable GDs per region as input for the tool ###
region_gds_truth <- gd_cluster_ccfs_sample[, .(num_detectable_gds = any(is_clone_region_1_gd_detectable) + any(is_clone_region_2_gd_detectable)),
                                           by = .(tumour_id, sample_id)]
data[, num_gds := region_gds_truth[ match( paste(data$tumour_id, data$sample_id),
                                           paste(region_gds_truth$tumour_id, region_gds_truth$sample_id) ), num_detectable_gds] ]
data[ is.na(num_gds), num_gds := 0 ] # for those tumours with no detectable gd
data[, sample_id := gsub('^.*SIM.{10}', '', sample_id) ]


#############
#### Run ####
#############

data_run <- data[, .(tumour_id, sample_id, chromosome, position, ref, alt, cluster_id,
                     is_clonal_cluster, mut_cpn, num_gds, MajCN = Maj_CN)]
output <- detect_par_gd( input = data_run, track = TRUE )
tumour <- copy( output[[1]] )
tumour[, truth_num_clonal_gds := truth_tumour[ match(tumour$tumour_id, truth_tumour$tumour_id), num_clonal_gds ]]
tumour[, truth_num_subclonal_gds := truth_tumour[ match(tumour$tumour_id, truth_tumour$tumour_id), num_subclonal_gds ]]
tumour[, truth_num_detectable_subclonal_gds := truth_tumour[ match(tumour$tumour_id, truth_tumour$tumour_id), num_detectable_subclonal_gds ]]
tumour[, truth_num_detectable_clonal_gds := truth_tumour[ match(tumour$tumour_id, truth_tumour$tumour_id), num_detectable_clonal_gds ]]
tumour[ is.na(truth_num_detectable_subclonal_gds), truth_num_detectable_subclonal_gds := 0 ]
tumour[ is.na(truth_num_detectable_clonal_gds), truth_num_detectable_clonal_gds := 0 ]

## Check determination of num gds against ground truth
data_plot[, num_detectable_gds := region_gds_truth[ match(data_plot[, paste(tumour_id, sample_id)],
                                                          region_gds_truth[, paste(tumour_id, sample_id)]), num_detectable_gds] ]
data_plot[ is.na(num_detectable_gds), num_detectable_gds := 0]
data_plot[, num_detectable_gds_lim := num_detectable_gds ]
data_plot[ num_detectable_gds > 2, num_detectable_gds_lim := 2 ]

####################
### Plot results ###
####################

# calculate what you would have got from older methods
NEJM_method <- function( gd_status, orig_method = TRUE ){
  gd_status <- as.numeric(gd_status)
  if(orig_method){
    gd_status[ gd_status > 1 ] <- 1
  }
  if( length(unique(gd_status)) == 1 ){
    subclonal_gd <- 0
    clonal_gd = unique(gd_status)
  } else {
    clonal_gd <- min(gd_status)
    subclonal_gd <- length(unique(gd_status)) - 1
  }
  return( c(clonal_gd, subclonal_gd))
}


tumour[, NEJM_method_num_cl_gds := NEJM_method(strsplit(GD_statuses, split=',')[[1]], orig_method = TRUE)[1] , 1:nrow(tumour)]
tumour[, NEJM_method_num_subcl_gds := NEJM_method(strsplit(GD_statuses, split=',')[[1]], orig_method = TRUE)[2] , 1:nrow(tumour)]
tumour[, NEJM_method_2gd_num_cl_gds := NEJM_method(strsplit(GD_statuses, split=',')[[1]], orig_method = FALSE)[1] , 1:nrow(tumour)]
tumour[, NEJM_method_2gd_num_subcl_gds := NEJM_method(strsplit(GD_statuses, split=',')[[1]], orig_method = FALSE)[2] , 1:nrow(tumour)]
data.plot.truth <- melt(tumour, id.vars = c('tumour_id'),
                        measure.vars = c('truth_num_detectable_subclonal_gds', 'truth_num_detectable_clonal_gds'),
                        variable.name = 'truth',
                        value = 'truth_GD_events')
data.plot.truth[, type := ifelse(grepl('_sub', truth), 'subclonal', 'clonal') ]
data.plot.measure <- melt(tumour, id.vars = c('tumour_id'),
                          measure.vars = c('num_clonal_gds', 'num_subclonal_gds',
                                           'NEJM_method_num_cl_gds', 'NEJM_method_num_subcl_gds',
                                           'NEJM_method_2gd_num_cl_gds', 'NEJM_method_2gd_num_subcl_gds'),
                          variable.name = 'method',
                          value = 'GD_events')


data.plot.measure[, type := ifelse(grepl('_sub', method), 'subclonal', 'clonal') ]
data.plot.measure[ grepl( 'al_gd', method ), method := 'ParallelGDDetect']
data.plot.measure[ grepl( 'NEJM_method_num', method ), method := 'NEJM method']
data.plot.measure[ grepl( 'NEJM_method_2gd', method ), method := 'NEJM 2nd method']
data.plot <- full_join(data.plot.measure, data.plot.truth[, .(tumour_id, type, truth_GD_events)] )
counts <- data.plot[, .(count = .N),
                    by = .(method, GD_events, truth_GD_events, type)]

pdf( paste0(outputs.folder, date, '_all_methods_vs_truth_facet.pdf'), width = 9 )

ggplot(data.plot[ type == 'subclonal' ], aes( x = method ) ) +
  geom_bar( aes(fill = method) ) +
  geom_text( data = counts[ type == 'subclonal'  ],
             aes(label = count, y = count - count *0.1, group = method))+
  scale_colour_brewer( palette = 'Paired' ) +
  labs( y = '# Tumours',
        x = '',
        title = 'Truth # subclonal GDs',
        fill = 'Method') +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Measured # subclonal GDs\n", labels = NULL, breaks = NULL)) +
  theme_classic() +
  theme(text = element_text( size = 20),
        axis.text.x = element_text( angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        plot.margin = margin(0.2,0.2,0.2,1.1, "cm")) +
  facet_grid( GD_events ~ truth_GD_events, scales = 'free')

dev.off()

#############################
## Extended figure 6 g, h and i ##
#############################

# Compare to SCNA ITH

evo_metrics[, wilcox.test( frac_abberant_genom_subcl[ num_subclonal_gds > 0],
                           frac_abberant_genom_subcl[ num_subclonal_gds == 0])]

evo_metrics[,num_subclonal_gds_plot :=  num_subclonal_gds]
evo_metrics[num_subclonal_gds > 2,num_subclonal_gds_plot :=  2]

pdf( paste(outputs.folder, '_subclonal_gd_vs_SCNA_ITH.pdf') )

ggplot(evo_metrics[ !is.na(num_subclonal_gds_plot)], 
       aes(x = as.factor(num_subclonal_gds_plot), y = frac_abberant_genom_subcl)) +
  geom_boxplot( aes(fill = as.factor(num_subclonal_gds_plot) ), outlier.shape = NA) +
  geom_jitter() +
  labs( x = '# Subclonal GDs',
        y = 'SCNA ITH', 
        fill = '# Subclonal GDs') +
  scale_fill_brewer( palette = 'Set1' )+
  theme_classic() +
  theme( text = element_text( size = 20))

dev.off()


# Compare to SNV ITH

evo_metrics[, wilcox.test( perc_subclonal[ num_subclonal_gds > 0],
                           perc_subclonal[ num_subclonal_gds == 0])]

pdf( paste(outputs.folder, date, '_subclonal_gd_vs_SNV_ITH.pdf') )

ggplot(evo_metrics[ !is.na(num_subclonal_gds_plot)], 
       aes(x = as.factor(num_subclonal_gds_plot), y = perc_subclonal )) +
  geom_boxplot( aes(fill = as.factor(num_subclonal_gds_plot) ), outlier.shape = NA) +
  geom_jitter() +
  labs( x = '# Subclonal GDs',
        y = 'SNV ITH', 
        fill = '# Subclonal GDs') +
  scale_fill_brewer( palette = 'Set1' )+
  theme_classic() +
  theme( text = element_text( size = 20))

dev.off()


# Compare to ABOPEC

evo_metrics[, wilcox.test( (SBS13_overall + SBS2_overall)[ num_subclonal_gds > 0],
                           (SBS13_overall + SBS2_overall)[ num_subclonal_gds == 0])]

pdf( paste(outputs.folder, date, '_subclonal_gd_vs_APOBECfrac.pdf') )

ggplot(evo_metrics[ !is.na(num_subclonal_gds_plot)], 
       aes(x = as.factor(num_subclonal_gds_plot), y = SBS13_overall + SBS2_overall )) +
  geom_boxplot( aes(fill = as.factor(num_subclonal_gds_plot) ), outlier.shape = NA) +
  geom_jitter() +
  labs( x = '# Subclonal GDs',
        y = 'Fraction ABOBEC mutations', 
        fill = '# Subclonal GDs') +
  scale_fill_brewer( palette = 'Set1' )+
  theme_classic() +
  theme( text = element_text( size = 20))

dev.off()

#############
#### END ####
#############

