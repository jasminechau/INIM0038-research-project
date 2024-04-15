#####################################################################################################
#############               Extended Figure Featuring Parallel Events              #############
#####################################################################################################
# written by Emilia Lim (emilia.lim@crick.ac.uk)


library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(cowplot)
library(reshape2)
library(tidyr)
library(DiagrammeR)
library(magrittr)
library(fst)
library(stringr)

addGeneAndCytoBandInfo = function(df,cytoband_df=NA){
  # CytoBand Information
  library(GenomicRanges)
  if(is.na(cytoband_df)){
    cytoband_df = getGeneAndCytoBandInfo()
  }
  cytoband_df = cytoband_df %>%
    dplyr::select(chromosome,start,end,gene_name,cytoband_name)
  gr_cytoband = with(cytoband_df, GRanges(chromosome, IRanges(start=start, end=end)))
  
  gr_cn = with(df, GRanges(chr, IRanges(start=startpos, end=endpos)))
  overlapping_regions = findOverlaps(gr_cn,gr_cytoband)
  overlapping_regions = data.frame(overlapping_regions)
  df_switch_cyto = cbind(data.frame(df)[overlapping_regions$queryHits,],
                         data.frame(cytoband_df)[overlapping_regions$subjectHits,])
  df_switch_cyto
  
}

areMutClustersParallel = function(tree_data,tree_name,mut_clusters){
  print(tree_name)
  clusters_that_muts_are_in = unique(unlist(strsplit(mut_clusters,"\\;")))
  clusters_that_muts_are_in = clusters_that_muts_are_in[clusters_that_muts_are_in!="NA"]
  if(length(clusters_that_muts_are_in)<2){
    return(NA)
  }
  
  c_auto_tree = tree_data[[tree_name]]$graph_pyclone$Corrected_tree
  c_auto_tree_trunk = tree_data[[tree_name]]$graph_pyclone$trunk
  
  merged_clusters_df = tree_data[[tree_name]]$merged_clusters
  if(!is.na(merged_clusters_df)&class(merged_clusters_df)=="character"){
    clusters_that_muts_are_in = ifelse(clusters_that_muts_are_in%in%merged_clusters_df[1],merged_clusters_df[3],clusters_that_muts_are_in)
    
  }
  
  if(!is.na(merged_clusters_df)&class(merged_clusters_df)=="matrix"){
    print("Merge Cluster - Matrix")
    clusters_that_muts_are_in = ifelse(clusters_that_muts_are_in%in%merged_clusters_df[1,1],merged_clusters_df[1,3],clusters_that_muts_are_in)
    
  }
  
  if(class(c_auto_tree)=="try-error"|class(c_auto_tree_trunk)=="try-error"|class(c_auto_tree)=="NULL"|class(c_auto_tree_trunk)=="NULL"){
    rr = NA
  }else{
    
    
    
    
    from_nodes <- as.numeric(as.character(c_auto_tree[,1]))
    to_nodes <- as.numeric(as.character(c_auto_tree[,2]))
    nodes <- sort(unique(c(from_nodes, to_nodes)))
    nodes_df <- tibble::tibble(id = nodes, type = NA, label = nodes)
    edges_df <- create_edge_df(from = from_nodes, to = to_nodes)
    
    all_paths_to_leafs = create_graph(nodes_df = nodes_df, edges_df = edges_df) %>%
      get_paths(from = c_auto_tree_trunk)
    
    all_paths_to_leafs_muts_in_path = lapply(all_paths_to_leafs,function(c_path){
      path_muts = clusters_that_muts_are_in[which(clusters_that_muts_are_in%in%c_path)]
    })
    all_paths_to_leafs_muts_in_path = unique(all_paths_to_leafs_muts_in_path)
    all_paths_to_leafs_num_muts_in_path = lapply(all_paths_to_leafs_muts_in_path,length)
    num_single_mut_paths = length(which(all_paths_to_leafs_num_muts_in_path==1))
    rr = num_single_mut_paths>1&num_single_mut_paths<=length(clusters_that_muts_are_in)
  }
  rr
}

findBestFitCluster = function(ccf_cluster_table,tree_removed_clusters,has_alteration_group){
  has_alteration_group = unlist(strsplit(has_alteration_group,"\\;"))
  
  if(length(setdiff(colnames(ccf_cluster_table),has_alteration_group))==0){
    best_fit_cluster = NA
  }else{
    ccf_cluster_table_summary = ccf_cluster_table %>%
      melt() %>%
      filter(!Var1%in%tree_removed_clusters) %>%
      mutate(has_alteration_group=ifelse(Var2%in%has_alteration_group,"has_alteration_group","no_alteration_group")) %>%
      group_by(Var1,has_alteration_group) %>%
      summarise(group_mean=mean(value)) %>%
      spread(has_alteration_group,group_mean,fill=0) %>%
      mutate(diff_between_groups=has_alteration_group-no_alteration_group)
    
    best_fit_cluster = ccf_cluster_table_summary$Var1[which(ccf_cluster_table_summary$diff_between_groups==max(ccf_cluster_table_summary$diff_between_groups))[1]]
  }
  best_fit_cluster
}

runfindBestFitCluster = function(tumour_oi,CN_Het_details_by_region_df_cyto,gene_anno_df,new_tree_input){
  print(tumour_oi)
  CN_Het_subclonal_by_patient_df_cyto2 = CN_Het_details_by_region_df_cyto %>%
    filter(tumour_id==tumour_oi) %>%
    distinct()
  
  CN_Het_subclonal_by_patient_gene_df = addGeneAndCytoBandInfo(CN_Het_subclonal_by_patient_df_cyto2,gene_anno_df) %>%
    select(-cytoband_name,-chromosome,-start,-end) %>%
    distinct()
  
  
  CN_Het_subclonal_by_patient_event_df = CN_Het_subclonal_by_patient_gene_df %>%
    mutate(region_id=str_sub(sample,-2,-1)) %>%
    mutate(CN_State=cpn_event_vs_ploidy) %>%
    mutate(seg_name=paste(chr,startpos,endpos,sep="_")) %>%
    filter(CN_State!="neutral") %>%
    mutate(event_name=paste(seg_name,gene_name,CN_State,sep="-"))
  
  ccf_cluster_table = new_tree_input$nested_pyclone$ccf_cluster_table
  colnames(ccf_cluster_table) = gsub("_LN",".LN",colnames(ccf_cluster_table))
  colnames(ccf_cluster_table) = do.call(rbind,strsplit(colnames(ccf_cluster_table),"\\."))[,2]
  tree_removed_clusters = c(new_tree_input$tree_removed_clusters,new_tree_input$cpn_removed_clusters)
  tree_removed_clusters = NA
  
  all_regions = paste(sort(unique(CN_Het_subclonal_by_patient_event_df$region_id)),collapse=";")
  region_list_df = CN_Het_subclonal_by_patient_event_df %>%
    group_by(event_name) %>%
    summarise(
      region_id_list = paste(unique(sort(region_id)),collapse=";")
    )
  
  region_list_best_fit_cluster_df = region_list_df %>%
    distinct(region_id_list) %>%
    filter(region_id_list!=all_regions) %>%
    rowwise() %>%
    mutate(best_fit_cluster = findBestFitCluster(ccf_cluster_table,tree_removed_clusters,region_id_list))
  
  CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster = region_list_df %>%
    mutate(tumour_id=tumour_oi) %>%
    left_join(.,region_list_best_fit_cluster_df,by=c("region_id_list"))
}

# setwd("/Volumes/lab-swantonc/working/lime/tct_lime/projects/tx_analysis/TRACERx421_data_code/20220808_Data/")

#Gene Annotation
gene_anno_df = fst::read.fst("gene_anno_df.fst")

bailey_gene_df = gene_anno_df %>%
  filter(!is.na(Bailey_TSG_Oncogene)) %>%
  distinct(gene_name,Bailey_TSG_Oncogene) %>%
  mutate(Bailey_TSG_Oncogene=gsub("possible ","",Bailey_TSG_Oncogene)) %>%
  arrange(Bailey_TSG_Oncogene,gene_name)

#Tree Data
tree_data_file = "20221109_TRACERx421_phylogenetic_trees.rds"
tree_data = readRDS(tree_data_file)
names(tree_data) = gsub("Cluster","Tumour",names(tree_data))

#SCNA Data
copy_number_file = "20221109_TRACERx421_scna_table.fst"
CN_Het_details_by_region_df_cyto = fst::read.fst(copy_number_file)
CN_Het_details_by_region_df_cyto$tumour_id = gsub("Cluster","Tumour",CN_Het_details_by_region_df_cyto$tumour_id)
tumour_id_list = unique(CN_Het_details_by_region_df_cyto$tumour_id)

patient_level_copy_number_file = "20221109_TRACERx421_scna_table_per_tumour.fst"
CN_Het_details_by_patient_df_cyto = read_fst(patient_level_copy_number_file)


#Finding best fit pyclone cluster for each SCNA
CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster_df = lapply((tumour_id_list),function(tumour_oi){
  print(tumour_oi)
  new_tree_input = try(tree_data[[tumour_oi]])
  if(class(new_tree_input)=="try-error"){
    scan_best_fit_cluster_df = NA
  }else{
    scan_best_fit_cluster_df = try(runfindBestFitCluster(tumour_oi,CN_Het_details_by_region_df_cyto,gene_anno_df,new_tree_input))
    if(class(scan_best_fit_cluster_df)=="try-error"){
      scan_best_fit_cluster_df = NA
    }
  }
  scan_best_fit_cluster_df
})

CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster_df2 = do.call(rbind,CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster_df) %>%
  filter(!is.na(best_fit_cluster))

#Summarising parallel SCNA events
MSAI_df = CN_Het_details_by_patient_df_cyto %>%
  filter(pt_seg_parallel_MSAI_gain_or_amp==TRUE|pt_seg_parallel_MSAI_loss_or_LOH==TRUE) %>%
  addGeneAndCytoBandInfo(.,cytoband_df=gene_anno_df)

MSAI_summary_df = MSAI_df %>%
  group_by(gene_name,tumour_id) %>%
  summarise(
    parallel_state=ifelse(pt_seg_parallel_MSAI_gain_or_amp==TRUE,"GainPara",""),
    parallel_state=ifelse(pt_seg_parallel_MSAI_loss_or_LOH==TRUE,"LossPara",parallel_state)
  ) %>%
  dplyr::rename(Hugo_Symbol=gene_name)

# Mut Table
muttable_file = "20221109_TRACERx421_mutation_table.fst"
muttable_df = fst::read.fst(muttable_file) 

driver_gene_df = muttable_df %>%
  filter(DriverMut==TRUE) %>%
  distinct(Hugo_Symbol)

mut_gene_tumour_count_df = muttable_df %>%
  filter(DriverMut==TRUE) %>%
  mutate(mutation_id_with_cluster=paste(PyCloneCluster_SC,mutation_id,sep="-")) %>%
  group_by(tumour_id,Hugo_Symbol) %>%
  summarise(mut_count=n(),
            mut_list=paste(((mutation_id_with_cluster)),collapse=";"),
            mut_cluster_count=length(unique(PyCloneCluster_SC[!is.na(PyCloneCluster_SC)])),
            mut_cluster_list=paste((unique(PyCloneCluster_SC)),collapse=";"))

CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster_df3 = CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster_df2 %>%
  separate(event_name,sep="-",into=c("coord","Hugo_Symbol","CN_State","CN_Extra","CN_Extra2")) %>%
  mutate(CN_State=ifelse(!is.na(CN_Extra),paste(CN_State,CN_Extra,sep="_"),CN_State)) %>%
  mutate(CN_State=ifelse(!is.na(CN_Extra2),paste(CN_State,CN_Extra2,sep="_"),CN_State)) %>%
  select(-CN_Extra,-CN_Extra2) %>%
  mutate(CN_State=ifelse(grepl("loss|LOH",CN_State),"loss",CN_State)) %>%
  mutate(CN_State=ifelse(grepl("amp|gain",CN_State),"gain",CN_State))


mut_gene_tumour_count_SCNA_df = mut_gene_tumour_count_df %>%
  left_join(.,CN_Het_subclonal_by_patient_df_cyto_best_fit_cluster_df3,by=c("Hugo_Symbol"="Hugo_Symbol","tumour_id"="tumour_id"))

mut_gene_tumour_count_SCNA_cluster_df = mut_gene_tumour_count_SCNA_df %>%
  mutate(scna_mut_cluster_list=ifelse(!is.na(mut_cluster_list)&!is.na(best_fit_cluster),paste(mut_cluster_list,best_fit_cluster,sep=";"),NA)) %>%
  rowwise() %>%
  mutate(isParallelwithSCNA=areMutClustersParallel(tree_data,tumour_id,scna_mut_cluster_list),
         isParallel=areMutClustersParallel(tree_data,tumour_id,mut_cluster_list)) %>%
  mutate(isParallel=ifelse(is.na(isParallel),FALSE,isParallel)) %>%
  mutate(isParallelwithSCNA=ifelse(is.na(isParallelwithSCNA),FALSE,isParallelwithSCNA)) %>%
  mutate(parallel_state=ifelse(isParallel==TRUE," SNV_Para","")) %>%
  mutate(parallel_state=ifelse(isParallelwithSCNA==TRUE,paste0(CN_State,"ParaWithSNV",parallel_state),parallel_state))

mut_gene_tumour_count_SCNA_cluster_SCNAonly_df = mut_gene_tumour_count_SCNA_cluster_df %>%
  mutate(parallel_state=ifelse(mut_count>1&parallel_state=="","multiSNV",parallel_state))



#Adding SCNA Parallel
mut_gene_tumour_count_SCNA_cluster_SCNAonly_df = bind_rows(mut_gene_tumour_count_SCNA_cluster_SCNAonly_df,MSAI_summary_df)

mut_gene_tumour_count_SCNA_cluster_SCNAonly_data = mut_gene_tumour_count_SCNA_cluster_SCNAonly_df %>%
  filter(parallel_state!="") %>%
  inner_join(.,bailey_gene_df,by=c("Hugo_Symbol"="gene_name")) %>%
  filter(Bailey_TSG_Oncogene%in%c("oncogene","tsg")) %>%
  group_by(Hugo_Symbol,Bailey_TSG_Oncogene,tumour_id) %>%
  summarise(parallel_state=paste(unique(parallel_state),collapse=";")) %>%
  mutate(parallel_state=gsub("gainParaWithSNV;GainPara","gainParaWithSNV",parallel_state)) %>%
  mutate(parallel_state=gsub("lossParaWithSNV;LossPara","lossParaWithSNV",parallel_state)) %>%
  mutate(parallel_state=gsub("gainParaWithSNV;lossParaWithSNV","lossParaWithSNV",parallel_state)) %>%
  mutate(parallel_state=gsub("lossParaWithSNV SNV_Para"," SNV_Para",parallel_state)) %>%
  mutate(parallel_state=gsub("gainParaWithSNV SNV_Para;lossParaWithSNV SNV_Para"," SNV_Para",parallel_state)) %>%
  mutate(parallel_state=gsub("gainParaWithSNV SNV_Para; SNV_Para"," SNV_Para",parallel_state)) %>%
  group_by(Hugo_Symbol,Bailey_TSG_Oncogene,parallel_state) %>%
  summarise(n=length(unique(tumour_id)))  %>%
  filter(Bailey_TSG_Oncogene=="oncogene"&parallel_state%in%c(" SNV_Para","GainPara","gainParaWithSNV","multiSNV")|
           Bailey_TSG_Oncogene=="tsg"&parallel_state%in%c(" SNV_Para","LossPara","lossParaWithSNV","multiSNV"))

gene_order = mut_gene_tumour_count_SCNA_cluster_SCNAonly_data%>%
  filter(parallel_state%in%c(" SNV_Para","GainPara","gainParaWithSNV","LossPara","lossParaWithSNV")) %>%
  group_by(Hugo_Symbol) %>%
  summarise(n=sum(n)) %>%
  arrange(-n) %>%
  pull(Hugo_Symbol)

color_bar = c("#ca0020","#f4a582","grey20","#92c5de","#0571b0","#fee08b")
names(color_bar) = c("gainParaWithSNV","GainPara"," SNV_Para","LossPara","lossParaWithSNV","multiSNV")

mut_gene_tumour_count_SCNA_cluster_SCNAonly_plot = mut_gene_tumour_count_SCNA_cluster_SCNAonly_data %>%
  filter(parallel_state%in%c(" SNV_Para","GainPara","gainParaWithSNV","LossPara","lossParaWithSNV")) %>%
  mutate(Hugo_Symbol=factor(Hugo_Symbol,levels=gene_order)) %>%
  mutate(parallel_state=factor(parallel_state,levels=c(" SNV_Para","GainPara","gainParaWithSNV","LossPara","lossParaWithSNV","multiSNV"))) %>%
  ggplot(aes(x=Hugo_Symbol,y=n,fill=parallel_state)) +
  geom_bar(stat="identity") +
  facet_wrap(~Bailey_TSG_Oncogene,scales="free",ncol=1) +
  ylab("Tumour Count") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values=color_bar)

mut_gene_tumour_count_SCNA_cluster_SCNAonly_plot
