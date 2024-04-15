########################################################################
#############               Germline Drivers               ############# 
########################################################################
# written by James Black (james.black@ucl.ac.uk)

library(dplyr)
library(ggplot2)
library(readr)

# setwd('/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/TRACERx_421/TRACERx421_data_code/20220808_Tumour_evoutionary_histories')

germline_strict <- read.table('../20221109_Tumour_evo_histories_DATA/20220808_TRACERx421_GermlineDrivers.tsv', header = T, sep = '\t')

germline_plot <- germline_strict %>% 
  group_by(Gene, Second_hit) %>% 
  summarise(count = n()) 
for (i in 1:nrow(germline_plot)){
  germline_plot$total_count[i] <- sum(germline_plot$count[germline_plot$Gene==germline_plot$Gene[i]])
}
custom_palette <- c(`2nd hit detected` = "chartreuse4", "1" = "white")

germline_plot <- as.data.frame(germline_plot)
germline_plot <- germline_plot[order(germline_plot$Gene, decreasing = FALSE),]
germline_plot$order_1 <- 1:nrow(germline_plot)
germline_plot <- germline_plot[order(germline_plot$total_count,germline_plot$order_1, decreasing = c(TRUE, FALSE)),]
germline_plot$order_2 <- 1:nrow(germline_plot)

ggplot(germline_plot, aes(reorder(Gene, -order_2), count, fill = Second_hit)) + geom_col(alpha = 0.5,colour = "black")+ theme_bw()+ 
  theme(axis.text.x = element_text(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title = "Cancer predisposition germline mutations in TRACERx 421", x = "", y = "Number of patients")+ scale_fill_manual(name = "", labels= c("", "Second hit detected"),values = custom_palette,
                                                                                                                               guide = guide_legend(override.aes = list(fill = c("transparent","chartreuse4"), color = c("transparent", "black"))))+
  coord_flip()
