##############################################################################################################
#############                  Extended Figure Featuring TX421 Cohort of Patients                #############
##############################################################################################################
# written by Emilia Lim (emilia.lim@crick.ac.uk)
# The male and female figures are part of the FontAwesome font package which must be installed. The font file
# "fontawesome-webfont.ttf" has been included with the upload of this script.

library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(fontawesome)
library(ggtext)

my_blank_gg = ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + theme_cowplot() + theme(
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.line = element_blank()
)

# load clinical data
path_to_clinical_data <-"20221109_TRACERx421_all_patient_df.rds"
clinical_df <- readRDS(path_to_clinical_data)

## tx421 patients
tx421_df <- clinical_df  %>%
  mutate(histology_group=ifelse(grepl("adenocarcinoma",histology_lesion1_merged),"Adenocarcinoma",as.character(histology_lesion1_merged))) %>%
  mutate(histology_group=ifelse(histology_group%in%c("Adenocarcinoma","Squamous cell carcinoma"),histology_group,"Other"))

hist_order = tx421_df %>%
  group_by(histology_group) %>%
  tally() %>%
  arrange(n) %>%
  pull(histology_group) %>%
  rev()

histology_colours = c("#c6dbf0","#e7b8bf","#8ed3c7","#fffdb3","#beb9da","#e7d3b7","grey")
names(histology_colours) = c("Adenocarcinoma","Squamous cell carcinoma","Adenosquamous carcinoma","Carcinosarcoma","Large cell carcinoma","Pleomorphic carcinoma","Other")

smoking_status_colours = c('#7993c9', "#a17db6", "#dea94a","#dea94a")
names(smoking_status_colours) = c("Never Smoked","Ex-Smoker","Recent Ex-Smoker","Smoker")

stage_colours = c("#eff3ff","#bdd7e7","#6baed6","#3182bd","#08519c","#022e56")
names(stage_colours) = c("IA","IB","IIA","IIB","IIIA","IIIB")


summary_anno_plot_df = tx421_df %>%
  mutate(smoking_status_merged=factor(smoking_status_merged,levels=rev(names(smoking_status_colours)))) %>%
  mutate(histology_group=factor(histology_group,levels=hist_order)) %>%
  arrange(histology_group,smoking_status_merged,pathologyTNM,sex)

tib <- tibble(
  family = c('firasans', 'lora', 'lobster', 'anton', 'syne'),
  x = 0,
  y = seq(0.0, 1, length.out = 5),
  label = "Let's talk cash <span style='font-family:fa-solid'></span>"
)
tib %>%
  ggplot(aes(x, y, label = label)) +
  geom_richtext(family = tib$family, size = 16, hjust = 0, col = 'dodgerblue4', label.colour = NA) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1.1)) +
  theme_void()

plotLittlePeople = function(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey90"){
  y_coord = unlist(lapply(1:nrows,function(x){rep(x,ncols)}))[1:nrow(summary_anno_plot_oi_df)]
  x_coord = rep(1:ncols,nrows)[1:nrow(summary_anno_plot_oi_df)]
  summary_anno_plot_oi_df$y_coord = y_coord
  summary_anno_plot_oi_df$x_coord = x_coord
  cohort_oi_plot = summary_anno_plot_oi_df %>%
    mutate(fa_sex=ifelse(sex=="Male",intToUtf8(61827),intToUtf8(61826))) %>%
    ggplot(aes(x=x_coord,y=y_coord)) +
    geom_tile(fill=bgcol,size=0.8) +
    scale_fill_manual(values=histology_colours) +
    geom_text(aes(label=fa_sex,color=smoking_status_merged),family="FontAwesome",size=8) +
    scale_color_manual(values=smoking_status_colours) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background = element_rect(fill = bgcol),
      plot.margin = margin(0,0,0,0, "cm")
    ) +
    xlim(0,ncols+1) +
    scale_y_reverse(limits=c(nrows+1,0))
}




stages_to_plot = summary_anno_plot_df$pathologyTNM %>% unique()
hists_to_plot = summary_anno_plot_df$histology_group %>% unique()

hist_oi = "Adenocarcinoma"
nrows = 7
ncols = 10

stage_oi = "IA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUADIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUADIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")

stage_oi = "IIA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUADIIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IIB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUADIIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")

stage_oi = "IIIA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUADIIIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IIIB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUADIIIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")




hist_oi = "Squamous cell carcinoma"
stage_oi = "IA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUSCIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUSCIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")

stage_oi = "IIA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUSCIIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IIB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUSCIIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")

stage_oi = "IIIA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUSCIIIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IIIB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_LUSCIIIB = try(plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85"))
if(class(cohort_oi_plot_LUSCIIIB)=="try-error"){cohort_oi_plot_LUSCIIIB = my_blank_gg}




hist_oi = "Other"
stage_oi = "IA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_OtherIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_OtherIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")

stage_oi = "IIA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_OtherIIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IIB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_OtherIIB = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85")

stage_oi = "IIIA"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_OtherIIIA = plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey95")

stage_oi = "IIIB"
summary_anno_plot_oi_df = summary_anno_plot_df %>%
  filter(pathologyTNM==stage_oi) %>%
  filter(histology_group==hist_oi)
cohort_oi_plot_OtherIIIB = try(plotLittlePeople(summary_anno_plot_oi_df,nrows,ncols,bgcol="grey85"))
if(class(cohort_oi_plot_OtherIIIB)=="try-error"){cohort_oi_plot_OtherIIIB = my_blank_gg}

patch_plot = (cohort_oi_plot_LUADIA+
    cohort_oi_plot_LUADIB+
    cohort_oi_plot_LUADIIA+
    cohort_oi_plot_LUADIIB+
    cohort_oi_plot_LUADIIIA+
    cohort_oi_plot_LUADIIIB+
    plot_layout(nrow = 1)) /
  (cohort_oi_plot_LUSCIA+
     cohort_oi_plot_LUSCIB+
     cohort_oi_plot_LUSCIIA+
     cohort_oi_plot_LUSCIIB+
     cohort_oi_plot_LUSCIIIA+
     cohort_oi_plot_LUSCIIIB+
     plot_layout(nrow = 1)) /
  (cohort_oi_plot_OtherIA+
     cohort_oi_plot_OtherIB+
     cohort_oi_plot_OtherIIA+
     cohort_oi_plot_OtherIIB+
     cohort_oi_plot_OtherIIIA+
     cohort_oi_plot_OtherIIIB+
     plot_layout(nrow = 1))



patch_plot

summary_anno_plot_df %>%
  group_by(histology_group) %>%
  tally()

