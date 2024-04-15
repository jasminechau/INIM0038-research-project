##########################################################
##########################################################
# Script to plot extended figures 9a/b
# Aug 2022
# Kristiana Grigoriadis
##########################################################
##########################################################

####################################################################################
#################### Script to plot supplementary figures 10a/b ####################
####################################################################################

# Author: Kristiana Grigoriadis
# Date: Aug 2022 // checked 22/11/22

# script to be run from root directory of github repo

####################################################################################
####################################### Setup ######################################
####################################################################################

# options
options(stringsAsFactors = F)

# libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fst))
suppressPackageStartupMessages(library(dndscv))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(rjson))

# working directory
mount_point = '/camp/project'

# tx421 repo path
tx421_release = file.path(mount_point, 'proj-tracerx-lung/tctProjects/grigork/clonal_mixing/Data/primaryEvoPaper/TRACERx421_data_code')
setwd(tx421_release)

tx421_data = '20221109_Tumour_evo_histories_DATA/'
data_outdir = file.path(mount_point, 'proj-tracerx-lung/tctProjects/grigork/clonal_mixing/Data/primaryEvoPaper/20221109/data/')
figures_outdir = file.path(mount_point, 'proj-tracerx-lung/tctProjects/grigork/clonal_mixing/Data/primaryEvoPaper/20221109/plots')

### load data
# trees
trees_path = file.path(tx421_data, '20221109_TRACERx421_phylogenetic_trees.rds')
trees421 = readRDS(trees_path)

# clinical info
clin_path = file.path(tx421_data, '20221109_TRACERx421_all_tumour_df.rds')
clin_df = as.data.table(readRDS(clin_path))
clin_df[, Histology := 'Other']
clin_df[grep('Invasive adenocarcinoma', Histology_per_tumour_id_muttable), Histology := 'LUAD']
clin_df[grep('Squamous cell carcinoma', Histology_per_tumour_id_muttable), Histology := 'LUSC']
clin_small = unique(clin_df[, c('tumour_id_muttable_cruk', 'cruk_id', 'Histology')])
colnames(clin_small) = c('tumour_id', 'patient_id', 'Histology')

# palette
palette_path = 'Palette/TRACERx_palette.json'
palette = fromJSON(file = palette_path)

##########################################################
### Functions
##########################################################
# function to process all alternative trees in Tx:
process_alt_trees = function(alt_trees){
  if (length(alt_trees) == 1) names(alt_trees) = 'Corrected_tree'
  else{
    names(alt_trees) = c('Corrected_tree', paste0('alt_tree_', 1:(length(alt_trees) - 1)))
  }
  return(alt_trees)
}

# get tree level of a clone
get_tree_graph_level = function(tree_graph, clone){
  ### Function to get tree level of a clone ###
  # trunk has level 1
  if (length(unique(as.numeric(tree_graph))) == 1){
    return(1)
  } else{
    trunk = unique(tree_graph[, 1][!tree_graph[, 1] %in% tree_graph[, 2]])
    if (clone == trunk) return(1)
    else {
      clones_in_tree = unique(as.numeric(as.matrix(tree_graph)))
      colnames(tree_graph) = c('Parent', 'Child')
      tree_graph = apply(tree_graph, c(1, 2), as.numeric)
      tree_graph = as.data.frame(tree_graph)
      
      level = 1
      if (clone %in% clones_in_tree){
        current_clone = as.numeric(clone)
        while (current_clone != trunk) {
          parent = tree_graph[tree_graph$Child==current_clone, 'Parent']
          level = level + 1
          current_clone = parent
        }  
        return(level)
      } else return(NA)
    }
  }
}


# compute clone proportions
calculateCloneProportionsTopDown = function(ccf_cluster_table, tree_graph, as_fraction = T){
  ### Function to get clone proportions ###
  # this function takes in a CCF table, and a tree matrix and returns a table of clone proportions
  clones_in_tree = unique(as.numeric(as.matrix(tree_graph)))
  if (length(clones_in_tree) == 1) {
    trunk = clones_in_tree
  } else {
    trunk = unique(tree_graph[, 1][!tree_graph[, 1] %in% tree_graph[, 2]])
  }
  
  colnames(tree_graph) = c('Parent', 'Child')
  tree_graph = as.data.frame(tree_graph)
  ccf_cluster_table = ccf_cluster_table[rownames(ccf_cluster_table) %in% clones_in_tree, , drop = F] # only include
  
  if (all(ccf_cluster_table < 50)) ccf_cluster_table = ccf_cluster_table*100
  # set trunk to CCF == 100 
  ccf_cluster_table[rownames(ccf_cluster_table) == trunk, ] = 100
  
  ccf_cluster_df = as.data.frame(ccf_cluster_table)
  ccf_cluster_df$clones = rownames(ccf_cluster_df)
  
  # this function will return proportion_df:
  # we fill up the empty data frame, starting from the trunk then stepping down the tree
  proportion_df = as.data.frame(ccf_cluster_df)
  region_IDs = colnames(ccf_cluster_table)
  proportion_df[, region_IDs] = 0
  
  for (r in region_IDs){
    cols = c(r, 'clones')
    region.CCF = ccf_cluster_df[, cols]
    colnames(region.CCF)[1] = c('ccf')
    clones.present = region.CCF[region.CCF$ccf != 0, 'clones']
    parents.present = unique(tree_graph$Parent[tree_graph$Parent %in% clones.present])
    
    # order the parent subclones by tree level (trunk = level 1)
    parent.df = data.frame(parents = parents.present)
    parent.df$level = unlist(lapply(parent.df$parents, function(p){
      return(get_tree_graph_level(as.matrix(tree_graph), p))
    } ))
    setorder(parent.df, level) 
    
    for (p in parent.df$parents){
      child_clusters = tree_graph[tree_graph$Parent == p, 'Child']
      parent.ccf = as.numeric(region.CCF[region.CCF$clones == p, 'ccf'])
      sum.daughters.ccf = sum(region.CCF[region.CCF$clones %in% child_clusters, 'ccf'])
      if (sum.daughters.ccf > parent.ccf){
        # print('sum.daughters.ccf > parent.ccf')
        # print(p)
        parent_proportion = 0
        # force children of trunk to be proportional values of parent ccf so that they sum to parent
        region.CCF[region.CCF$clones %in% child_clusters, ]$ccf = parent.ccf*region.CCF[region.CCF$clones %in% child_clusters, 'ccf']/sum(region.CCF[region.CCF$clones %in% child_clusters, 'ccf'])
      } else {
        parent_proportion = parent.ccf - sum.daughters.ccf
      }
      proportion_df[proportion_df$clones == p, r] = parent_proportion
      # if the clones are terminal, add them to proportion_df as well
      for (d in child_clusters){
        if (d %in% tree_graph$Child & !(d %in% tree_graph$Parent)) {
          proportion_df[proportion_df$clones == d, r] = region.CCF[region.CCF$clones == d, 'ccf']
        }
      }
    }  
  }
  proportion_df = proportion_df[, region_IDs, drop = F]
  # fractions or percentages?
  if (as_fraction == T){
    if (max(as.matrix(proportion_df)) > 10) proportion_df = proportion_df/100
  }
  return(proportion_df)   
}


# function to get terminal clusters in the tree
get_terminal_clusters = function(tree_graph){
  term_clusters = tree_graph[, 2][!tree_graph[, 2] %in% tree_graph[, 1]]
  # if tree is a clonal tree (ie no subclones), return the clonal cluster
  if (length(unique(as.numeric(tree_graph))) == 1){
    term_clusters = as.character(unique(as.numeric(tree_graph)))
  }
  return(as.character(term_clusters))
}

# function to get maximum ccf of the leaf nodes in each region (there might not be any)
get_max_ccf_terminal_clusters_reg = function(tree_graph, ccf_cluster_table, as_fraction = T){
  regions = colnames(ccf_cluster_table)
  # get terminal clusters in tumour
  term_clusters = get_terminal_clusters(tree_graph)
  
  max_ccf_terminal_in_regs = sapply(regions, function(r) {
    return(max(ccf_cluster_table[rownames(ccf_cluster_table) %in% term_clusters, r]))
  })
  # max_ccf_terminal_tum = max(max_ccf_terminal_in_regs)
  if (as_fraction == T){
    if (any(as.matrix(max_ccf_terminal_in_regs) > 2)) max_ccf_terminal_in_regs = max_ccf_terminal_in_regs/100
  }
  
  # check whether tumour is clonal
  if (length(unique(as.numeric(tree_graph))) == 1){
    max_ccf_terminal_in_regs = 1* (max_ccf_terminal_in_regs == max_ccf_terminal_in_regs)
  }
  
  # cap at 1:
  max_ccf_terminal_in_regs[max_ccf_terminal_in_regs > 1] = 1
  
  return(max_ccf_terminal_in_regs)
}

# get which leaf node has maximum ccf in each region
get_leafnodes_with_max_ccf_reg = function(tree_graph, ccf_cluster_table){
  regions = colnames(ccf_cluster_table)
  # get terminal clusters in tumour
  # term_clusters = tree_graph[, 2][!tree_graph[, 2] %in% tree_graph[, 1]]
  term_clusters = get_terminal_clusters(tree_graph)
  
  which_max_ccf_terminal_in_regs = sapply(regions, function(r) {
    term_ccfs = ccf_cluster_table[rownames(ccf_cluster_table) %in% term_clusters, r, drop = F]
    if (all(as.numeric(term_ccfs) == 0)) return(NA)
    else {
      maxCCF_term_clusters = rownames(term_ccfs)[term_ccfs == max(term_ccfs)]
      return(paste(maxCCF_term_clusters, collapse = ';'))
    }
  })
  return(which_max_ccf_terminal_in_regs)
}

get_max_ccf_ancestral_clusters_reg = function(tree_graph, ccf_cluster_table, as_fraction = T){
  regions = colnames(ccf_cluster_table)
  
  clones_in_tree = unique(as.numeric(tree_graph))
  if (length(clones_in_tree) == 1) {
    trunk = clones_in_tree
  } else {
    trunk = unique(tree_graph[, 1][!tree_graph[, 1] %in% tree_graph[, 2]])
  }
  
  # get terminal clusters in tumour
  term_clusters = get_terminal_clusters(tree_graph)
  ancestral_subclones = clones_in_tree[!clones_in_tree %in% term_clusters & !clones_in_tree == trunk]
  
  if (length(ancestral_subclones) == 0){
    max_ccf_ancestral_in_regs = rep(NA, length(regions))
    names(max_ccf_ancestral_in_regs) = regions
  } else{
    max_ccf_ancestral_in_regs = sapply(regions, function(r) {
      return(max(ccf_cluster_table[rownames(ccf_cluster_table) %in% ancestral_subclones, r]))
    })
    
    if (as_fraction == T){
      if (any(as.matrix(max_ccf_ancestral_in_regs) > 2)) max_ccf_ancestral_in_regs = max_ccf_ancestral_in_regs/100
    }
  }
  
  # cap at 1:
  max_ccf_ancestral_in_regs[max_ccf_ancestral_in_regs > 1] = 1
  
  return(max_ccf_ancestral_in_regs)
}

process_mean_cluster_ccfs = function(ccf_cluster_table, phyloCCF_df){
  phyloCCF_columns = colnames(phyloCCF_df)[grep('PhyloCCF', colnames(phyloCCF_df))]
  regions = gsub('_PhyloCCF', '', phyloCCF_columns)
  
  meanCCF_dt = as.data.table(phyloCCF_df[, c(phyloCCF_columns, 'PycloneCluster')])
  meanCCF_dt[, (regions) := lapply(.SD, mean), by = PycloneCluster, .SDcols = phyloCCF_columns]
  meanCCF_cols = c(regions, 'PycloneCluster')
  meanCCF_dt = unique(meanCCF_dt[, ..meanCCF_cols])
  setorder(meanCCF_dt, 'PycloneCluster')
  
  mean_ccf_cluster_table = as.data.frame(meanCCF_dt)
  rownames(mean_ccf_cluster_table) = mean_ccf_cluster_table$PycloneCluster
  mean_ccf_cluster_table = as.matrix(mean_ccf_cluster_table[, regions, drop = F])
  
  if (any(abs(ccf_cluster_table - mean_ccf_cluster_table*100) > 1)) print('Difference in mean CCF and cluster CCF > 1 !!!')
  return(mean_ccf_cluster_table * 100)
}

##########################################################
### Extended figure 9a
##########################################################

### get clonal illusion clusters
cluster_df = rbindlist(lapply(names(trees421), function(n){
  t = trees421[[n]]
  corrected_tree = t$graph_pyclone$Corrected_tree
  clusters_in_tree = unique(as.numeric(corrected_tree))
  clonality_table = t$clonality_out$clonality_table_corrected
  ccf_cluster_table = t$nested_pyclone$ccf_cluster_table
  
  out_dt = data.table(tumour_id = n, 
                      cluster = sort(clusters_in_tree))
  
  # add column: whether cluster is clonal in any region (from clonality table)
  out_dt$is_CI = sapply(out_dt$cluster, function(c){
    any(clonality_table[rownames(clonality_table) == c,] == 'clonal')
  })
  
  # add column: whether cluster is trunk
  out_dt$is_trunk = out_dt$cluster == t$graph_pyclone$trunk
  
  return(out_dt)
}))

# extract only subclonal clusters
subclone_df = cluster_df[is_trunk == F]
subclone_df$is_trunk = NULL

# merge with histology data
subclone_df = merge.data.table(subclone_df, clin_small, by = 'tumour_id', all.x = T, all.y = F)

# get number of clonal illusion per tumour
subclone_df[, no_subclones := .N, by = 'tumour_id']
subclone_df[, no_CI_tumour := sum(is_CI), by = 'tumour_id']
subclone_df[, fraction_CI_tumour := no_CI_tumour / no_subclones]


# compare number of clonal illusion clusters in LUAD vs LUSC
# get tumour level table
tumour_columns = c('tumour_id', 'patient_id', 'Histology', 'no_subclones', 'no_CI_tumour', 'fraction_CI_tumour')
tumour_df = unique(subclone_df[, ..tumour_columns])

luad_lusc_palette = palette$histology_subtypes[c('LUAD', 'LUSC')]

g = ggplot(tumour_df[Histology %in% c('LUAD', 'LUSC')], aes(Histology, no_CI_tumour))+
  geom_boxplot(aes(fill = Histology))+
  stat_compare_means()+
  theme_classic()+
  labs(y = 'Number of clonal illusion subclones')+
  scale_fill_manual(values = luad_lusc_palette, breaks = names(luad_lusc_palette))
pdf(file.path(figures_outdir, 'suppfig10a.pdf'), width = 5, height = 4)
print(g)
dev.off()


# write table
# write.csv(subclone_df, file.path(tx421_data, 'subclone_df.20220816.csv'), row.names = F)
##########################################################
### Extended figure 9b
##########################################################

# Generate diversity metrics table
# Metrics to calculate:

# 1 Maximum CCF of only terminal (leaf) clones in each region. (If there are none, max ccf == 0)
# 1.1 Maximum of above across all regions of tumour

# 2 Maximum CCF of only ancestral (parental) clones in each region. (If there are none, max ccf == NAs)
# 2.1 Maximum of above across all regions of tumour

# 3. Simpson's diversity per region
# 3.1 Minimum regional Simpson's diversity across the tumour

# 4. Is region clonally pure based on clonality table?

### Note: exclude LN regions

cloneDivAlttreesDf = rbindlist(lapply(names(trees421), function(n){
  print(n)
  t = trees421[[n]]
  trunk = t$graph_pyclone$trunk
  clusters_in_tree = unique(as.numeric(t$graph_pyclone$Corrected_tree))
  
  # get data
  ccf_cluster_table = t$nested_pyclone$ccf_cluster_table
  clonality_table = t$clonality_out$clonality_table_corrected
  phyloCCF_df = as.data.frame(t$ccf_table_pyclone_clean)
  
  # restrict to only clusters in the tree
  ccf_cluster_table = ccf_cluster_table[rownames(ccf_cluster_table) %in% clusters_in_tree, , drop = F] # only use clusters in tree
  clonality_table = clonality_table[rownames(clonality_table) %in% clusters_in_tree, , drop = F] # only use clusters in tree
  phyloCCF_df = phyloCCF_df[phyloCCF_df$PycloneCluster %in% clusters_in_tree, , drop = F]
  
  # 0. remove LN regions:
  ccf_cluster_table = ccf_cluster_table[, !grepl('LN', colnames(ccf_cluster_table)), drop = F]
  clonality_table = clonality_table[, !grepl('LN', colnames(clonality_table)), drop = F]
  phyloCCF_df = phyloCCF_df[, !grepl('LN', colnames(phyloCCF_df)), drop = F]
  
  if (ncol(ccf_cluster_table) > 0){
    
    # process mean CCF of mutations within each cluster:
    mean_ccf_cluster_table = process_mean_cluster_ccfs(ccf_cluster_table, phyloCCF_df)
    
    ####################################################################################################
    # get all alternative trees:
    alt_trees = process_alt_trees(t$graph_pyclone$alt_trees)
    
    alt_trees_df = rbindlist(lapply(names(alt_trees), function(cur_tree){
      tree_graph = alt_trees[[cur_tree]]
      ##################################################
      ###### Steps 1-3 using ccf_cluster_table
      # 1. max CCF of only terminal clones
      max_ccf_terminalonly_reg = get_max_ccf_terminal_clusters_reg(tree_graph, ccf_cluster_table, as_fraction = T)
      # which clones are these
      max_ccf_leafnodes = get_leafnodes_with_max_ccf_reg(tree_graph, ccf_cluster_table)
      maxCCF_terminalonly_tum = max(max_ccf_terminalonly_reg)
      
      # 2. max CCF of any ancestral clone
      max_ccf_ancestralonly_reg = get_max_ccf_ancestral_clusters_reg(tree_graph, ccf_cluster_table, as_fraction = T)
      maxCCF_ancestralonly_tum = max(max_ccf_ancestralonly_reg, na.rm = T)
      
      # 3. simpson's diversity per region
      # get clone proportions:
      clone_props = calculateCloneProportionsTopDown(ccf_cluster_table, tree_graph, as_fraction = T)
      simpson = apply(clone_props, c(2), function(x) diversity(x, index = 'simpson'))
      minSimpson_tum = min(simpson)
      ##################################################
      ##### Repeat steps 1-3 but using mean CCF of mutations per cluster:
      # 1. max CCF of only terminal clones
      max_ccf_terminalonly_reg_meanccf = get_max_ccf_terminal_clusters_reg(tree_graph, mean_ccf_cluster_table, as_fraction = T)
      # which clones are these
      max_ccf_leafnodes_meanccf = get_leafnodes_with_max_ccf_reg(tree_graph, mean_ccf_cluster_table)
      maxCCF_terminalonly_tum_meanccf = max(max_ccf_terminalonly_reg_meanccf)
      
      # 2. max CCF of any ancestral clone
      max_ccf_ancestralonly_reg_meanccf = get_max_ccf_ancestral_clusters_reg(tree_graph, mean_ccf_cluster_table, as_fraction = T)
      maxCCF_ancestralonly_tum_meanccf = max(max_ccf_ancestralonly_reg_meanccf, na.rm = T)
      
      # 3. simpson's diversity per region
      # get clone proportions:
      clone_props_meanccf = calculateCloneProportionsTopDown(mean_ccf_cluster_table, tree_graph, as_fraction = T)
      simpson_meanccf = apply(clone_props_meanccf, c(2), function(x) diversity(x, index = 'simpson'))
      minSimpson_tum_meanccf = min(simpson_meanccf)
      
      
      ####################################################################################################
      # 4. is region clonally pure based on clonality table? are there subclones present?
      clonally_pure_reg = apply(clonality_table, 2, function(r) all(r %in% c("clonal", "absent")))
      subclones_present_reg = apply(X = clonality_table[rownames(clonality_table) != trunk, , drop = F], 
                                    MARGIN = 2, 
                                    FUN = function(r) any(r %in% c("subclonal", "clonal")))
      
      ####################################################################################################
      # 5. summarise in data frame
      clone_div_df = data.table(tumour_id = n,
                                tree_name = cur_tree,
                                sample = names(max_ccf_terminalonly_reg),
                                max_ccf_terminalonly_reg = max_ccf_terminalonly_reg,
                                max_ccf_leafnodes = max_ccf_leafnodes,
                                maxCCF_terminalonly_tum = maxCCF_terminalonly_tum,
                                max_ccf_ancestralonly_reg = max_ccf_ancestralonly_reg,
                                maxCCF_ancestralonly_tum = maxCCF_ancestralonly_tum,
                                simpson_div_reg = simpson,
                                min_simpson_tum = minSimpson_tum,
                                
                                max_ccf_terminalonly_reg_meanccf = max_ccf_terminalonly_reg_meanccf,
                                max_ccf_leafnodes_meanccf = max_ccf_leafnodes_meanccf,
                                maxCCF_terminalonly_tum_meanccf = maxCCF_terminalonly_tum_meanccf,
                                max_ccf_ancestralonly_reg_meanccf = max_ccf_ancestralonly_reg_meanccf,
                                maxCCF_ancestralonly_tum_meanccf = maxCCF_ancestralonly_tum_meanccf,
                                simpson_div_reg_meanccf = simpson_meanccf,
                                min_simpson_tum_meanccf = minSimpson_tum_meanccf,
                                
                                is_clonal_tree = (length(unique(as.numeric(tree_graph))) == 1),
                                clonally_pure_reg = clonally_pure_reg,
                                subclones_present_reg = subclones_present_reg)
      
      return(clone_div_df)
    }))
    
    return(alt_trees_df)
  }
}))

# fix trees with only clonal cluster: simpson_div_reg == min_simpson_tum == 0
cloneDivAlttreesDf[is_clonal_tree==T, simpson_div_reg := 0]
cloneDivAlttreesDf[is_clonal_tree==T, min_simpson_tum := 0]
cloneDivAlttreesDf[is_clonal_tree==T, simpson_div_reg := 0]
cloneDivAlttreesDf[is_clonal_tree==T, min_simpson_tum := 0]

write.csv(cloneDivAlttreesDf, file.path(data_outdir, 'cloneDivAlttreesDf.20221118.csv'), row.names = F)

write.csv(cloneDivCorrectedTree, file.path(data_outdir, 'cloneDivCorrectedTree.20221118.csv'), row.names = F)


##########################################################
### Save data frame generated above for survival analysis
##########################################################
# get tumour level metrics
tum_columns = c("tumour_id", 
                'tree_name', 
                'maxCCF_terminalonly_tum', 
                'maxCCF_ancestralonly_tum', 
                'maxCCF_terminalonly_tum_meanccf',
                'maxCCF_ancestralonly_tum_meanccf',
                'min_simpson_tum', 
                'min_simpson_tum_meanccf',
                'is_clonal_tree')
tumourDivAlttreesDf = unique(cloneDivAlttreesDf[, ..tum_columns])

# to get patient level metrics, take average across multiple tumour clusters:
tumourDivAlttreesDf[, patient_id := gsub('_Tumour.*', '', tumour_id)]

# summarise average score of maxCCF_terminalonly_tum across multiple trees:
#### Cluster CCF
tumourDivAlttreesDf[, meantree_maxCCF_terminalonly_tum := mean(maxCCF_terminalonly_tum), by = c('tumour_id')]
tumourDivAlttreesDf[, maxtree_maxCCF_terminalonly_tum := max(maxCCF_terminalonly_tum), by = c('tumour_id')]
tumourDivAlttreesDf[, mintree_maxCCF_terminalonly_tum := min(maxCCF_terminalonly_tum), by = c('tumour_id')]
tumourDivAlttreesDf[, mediantree_maxCCF_terminalonly_tum := median(maxCCF_terminalonly_tum), by = c('tumour_id')]
tumourDivAlttreesDf[, sdtree_maxCCF_terminalonly_tum := sd(maxCCF_terminalonly_tum), by = c('tumour_id')]
tumourDivAlttreesDf[is.na(sdtree_maxCCF_terminalonly_tum), sdtree_maxCCF_terminalonly_tum := 0]
tumourDivAlttreesDf[, no_uniq_scores := length(unique(maxCCF_terminalonly_tum)), by = c('tumour_id')] # check how many unique scores per tree we get

#### Mean CCF of mutations making up a cluster
tumourDivAlttreesDf[, meantree_maxCCF_terminalonly_tum_meanccf := mean(maxCCF_terminalonly_tum_meanccf), by = c('tumour_id')]
tumourDivAlttreesDf[, maxtree_maxCCF_terminalonly_tum_meanccf := max(maxCCF_terminalonly_tum_meanccf), by = c('tumour_id')]
tumourDivAlttreesDf[, mintree_maxCCF_terminalonly_tum_meanccf := min(maxCCF_terminalonly_tum_meanccf), by = c('tumour_id')]
tumourDivAlttreesDf[, mediantree_maxCCF_terminalonly_tum_meanccf := median(maxCCF_terminalonly_tum_meanccf), by = c('tumour_id')]
tumourDivAlttreesDf[, sdtree_maxCCF_terminalonly_tum_meanccf := sd(maxCCF_terminalonly_tum_meanccf), by = c('tumour_id')]
tumourDivAlttreesDf[is.na(sdtree_maxCCF_terminalonly_tum_meanccf), sdtree_maxCCF_terminalonly_tum_meanccf := 0]
tumourDivAlttreesDf[, no_uniq_scores_meanccf := length(unique(maxCCF_terminalonly_tum_meanccf)), by = c('tumour_id')] # check how many unique scores per tree we get

# summarise average score of diversity across multiple trees:
tumourDivAlttreesDf[, meantree_min_simpson_tum := mean(min_simpson_tum), by = c('tumour_id')]
tumourDivAlttreesDf[, maxtree_min_simpson_tum := max(min_simpson_tum), by = c('tumour_id')]
tumourDivAlttreesDf[, meantree_min_simpson_tum_meanccf := mean(min_simpson_tum_meanccf), by = c('tumour_id')]
tumourDivAlttreesDf[, maxtree_min_simpson_tum_meanccf := max(min_simpson_tum_meanccf), by = c('tumour_id')]
# write
write.csv(tumourDivAlttreesDf, file.path(data_outdir, '20221102_tumourDivAlttreesDf.csv'), row.names = F)

########################################################## Write score for one tree per tumours
### Corrected tree only - tumour level metrics
tumour_columns = c("patient_id", "tumour_id", 
                  "maxCCF_terminalonly_tum", "maxCCF_ancestralonly_tum", "min_simpson_tum",
                  "maxCCF_terminalonly_tum_meanccf", "maxCCF_ancestralonly_tum_meanccf", "min_simpson_tum_meanccf")

tumourDiv_corrected = tumourDivAlttreesDf[tree_name == 'Corrected_tree']
write.csv(tumourDiv_corrected[, ..tumour_columns], file.path(data_outdir, '20221102_tumourDivCorrectedtree.csv'), row.names = F)

### Min tree only - tumour level metrics
tumourDiv_mintree = tumourDivAlttreesDf[maxCCF_terminalonly_tum_meanccf == mintree_maxCCF_terminalonly_tum_meanccf]
# subclonal expansion score:
maxCCFterminal_cols = c("patient_id", "tumour_id", "mintree_maxCCF_terminalonly_tum", "mintree_maxCCF_terminalonly_tum_meanccf")
tumourDiv_mintree = unique(tumourDivAlttreesDf[, ..maxCCFterminal_cols])
write.csv(tumourDiv_mintree, file.path(data_outdir, '20221102_tumourDivMintree_subclonalExpScore.csv'), row.names = F)
# Simpson's diversity score:
simpson_cols = c("patient_id", "tumour_id", "maxtree_min_simpson_tum", "maxtree_min_simpson_tum_meanccf")
tumourDiv_simpson = unique(tumourDivAlttreesDf[, ..simpson_cols])
write.csv(tumourDiv_simpson, file.path(data_outdir, '20221102_tumourDivMintree_Simpson.csv'), row.names = F)

##########################################################
### Extended figure 9c
##########################################################

# corrected tree only
cloneDivCorrectedTree = cloneDivAlttreesDf[tree_name == 'Corrected_tree']
cloneDiv_corrected_sc = cloneDivCorrectedTree[subclones_present_reg == T]

# whether region is higher or lower than median for max ccf of leaf node
cloneDiv_corrected_sc[, large_recent_expansion_region := ifelse(max_ccf_terminalonly_reg_meanccf > median(max_ccf_terminalonly_reg_meanccf),
                                                      '> median',
                                                      '<= median')]

# recent subclonal expansion size vs simpson diversity
g = ggplot(cloneDiv_corrected_sc, aes(large_recent_expansion_region, simpson_div_reg))+
  stat_compare_means()+
  theme_classic()+
  labs(x = 'Size of largest recently expanded subclone in region',
       y = 'Clonal diveristy of region (Simpson)')+
  theme(axis.title = element_text(size = 12))+
  geom_jitter(alpha = .5, size = 2)+
  geom_boxplot(alpha = .6)
pdf(file.path(figures_outdir, 'suppfig10b.pdf'), width = 5, height = 4)
print(g)
dev.off()


##########################################################
### Extended figure 9f
##########################################################
corrected_long = melt.data.table(tumourDiv_corrected, 
                                 measure.vars = c('maxCCF_terminalonly_tum', 'maxCCF_ancestralonly_tum'))
corrected_long[variable == 'maxCCF_ancestralonly_tum', `Subclone type` := 'Ancestral']
corrected_long[variable == 'maxCCF_terminalonly_tum', `Subclone type` := 'Recent']

g = ggplot(corrected_long, aes(value))+
  geom_density(alpha = .6, aes(fill = `Subclone type`))+
  theme_classic()+
  labs(x = 'Largest subclonal expansion\nin any region')+
  scale_fill_manual(values = brewer.pal(3, 'Paired')[c(2,3)])+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        text = element_text(family = 'Helvetica'))
pdf(file.path(figures_outdir, 'extendFig9/extendfig9e.pdf'), width = 8, height = 6)
print(g)
dev.off()

















