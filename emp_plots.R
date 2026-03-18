# library
library(ggplot2)
library(dplyr)
library(forcats)
theme_set(theme_light())
library(VennDiagram)
library(ComplexHeatmap)
library(ggtext)
library(UpSetR)
library(ComplexUpset)
library(tidyr)
library(stringr)
library(purrr)

######################################## STEADY STATES ##########################################
steady_states <- read.csv('steady_states.csv')
## all nodes
column_ha= HeatmapAnnotation(df=steady_states['Phenotype'],
                             annotation_legend_param = list(nrow=4), 
                             col = list(Phenotype = c ("M"="#9966ff", "E"="#99ccff", "pEM"="#ff9933", "Naive"="black")))
df=t(data.matrix(subset(steady_states, select = -Phenotype)))

colors = structure(c(2,3), names = c("0", "1")) 

png("steady states.png", width = 1000, height = 1500, res = 120)
h1=Heatmap(df, 
           top_annotation = column_ha,
           column_split = steady_states$Phenotype,
           show_column_names = TRUE, col=colors, 
           heatmap_legend_param = list(title='value',legend_direction = "vertical",   
                                       legend_width = unit(5, "cm")))

draw(h1, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# subset of nodes

subset_nodes=c("GRHL2", "P63", "miR200", "miR205", "OVOL1", "OVOL2", 
               "TightJunc","AdhJunc", "ECAD", "TWIST", "SNAI1","SNAI2", "ZEB1", "ZEB2", 
               "PRRX1", "BCAT_TCF", "NCAD","VIM", "MMP2", "Phenotype", "BMP", "TGFB_L", "GF")

short_ss=steady_states %>% subset(select = subset_nodes)
column_ha= HeatmapAnnotation(df=short_ss[c('Phenotype', 'BMP', 'TGFB_L', 'GF')],
                             annotation_legend_param = list(nrow=4), 
                             col = list(Phenotype = c ("M"="#9966ff", "E"="#99ccff", "pEM"="#ff9933", "Naive"="black"),
                                        BMP = c("0"="blue", "1"= "brown1"), 
                                        TGFB_L = c("0"="blue", "1"="brown1"), 
                                        GF= c("0"="blue", "1"="brown1")))
df=t(data.matrix(subset(short_ss, select = -c(Phenotype, GF, TGFB_L, BMP))))

colors = structure(c(2,3), names = c("0", "1")) 

row_split = rep("Epithelial factors", nrow(df))
row_split[10:nrow(df)]= "Mesenchymal factors"
h1=Heatmap(df, 
           top_annotation = column_ha,
           column_split = short_ss$Phenotype,
           row_split = row_split,
           cluster_rows = FALSE,
           show_column_names = FALSE, col=colors, 
           heatmap_legend_param = list(title='value',legend_direction = "vertical",   
                                       legend_width = unit(5, "cm")))

png("steady states_summary.png", width = 1000, height = 1000, res = 120)
draw(h1, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

######################################## FIXED POINTS EMT  ##########################################
fp <- read.csv("final_states_emt")
column_ha= HeatmapAnnotation(df=fp[c('Phenotype', 'BMP', 'TGFB_L', 'GF')],
                             annotation_legend_param = list(nrow=2), 
                             col = list(Phenotype = c ("M"="#9966ff", "E"="#99ccff", "pEM"="#ff9933"),
                                        BMP = c("0"="blue", "1"= "brown1"), 
                                        TGFB_L = c("0"="blue", "1"="brown1"), 
                                        GF= c("0"="blue", "1"="brown1")))
#df=t(unique(data.matrix(subset(fp, select = -c(Phenotype, X, input, Proba, BMP, GF, TGFB_L)))))
df=t(data.matrix(subset(fp, select = -c(Phenotype, X, input, Proba, BMP, GF, TGFB_L, E, M, pEM))))

colors = structure(c(2,3), names = c("0", "1")) 

row_split = rep("Mesenchymal factors", nrow(df))
row_split[10:nrow(df)]= "Epithelial factors"
#column_split[9:16] = "group2"
#column_split[17:24] = "group3"

png("final_states_from_EpiIC.png", width = 1000, height = 1000, res = 130)
Heatmap(df, 
           top_annotation = column_ha,
            column_split = fp$Phenotype,
            row_split = row_split,
            cluster_rows = FALSE,
           show_column_names = FALSE, col=colors, 
           heatmap_legend_param = list(title='value',legend_direction = "vertical",   
                                       legend_width = unit(5, "cm")))
dev.off()

#################################################################################################
#################################### MUTATIONS ON MET DRIVERS (MEsIC) ####################################

#################################### TGFBL dependent 
data <- read.csv('specific_mutations_mesenchymalcells_externalTGFBL.txt')
# Save original mutation order
mut_levels <- unique(data$mutation)

data <- data %>%
  mutate(
    mutation = as.character(mutation),  # Convert factor to character
    mutation = ifelse(mutation == "Wild type",
                      paste0("<b>", mutation, "</b>"),
                      mutation),
    mutation = factor(mutation, levels = unique(mutation))  # Use unique levels only
  )
# Plot
ggplot(data, aes(fill = phenotypes, y = mutation, x = proportion)) +
  geom_bar(position = "fill", stat = "identity", width = 0.4) +
  scale_fill_manual(values = c("#99ccff", "#9966ff", "#808080", "#ff9933")) +
  theme(
    text = element_text(size = 12, family = "Open Sans"),
    axis.text.y = element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  ggtitle('MesIC under TGFB_L=1')
ggsave("MesICmut_under_GFB_input_signal.png", width = 6, height = 5, device='png', dpi=700)

#################################### TGFBL independent 
data <- read.csv('specific_mutations_mesenchymalcells_no_external_TGFBL.txt')

# Save original mutation order
mut_levels <- unique(data$mutation)

data <- data %>%
  mutate(
    mutation = as.character(mutation),  # Convert factor to character
    mutation = ifelse(mutation == "Wild type",
                      paste0("<b>", mutation, "</b>"),
                      mutation),
    mutation = factor(mutation, levels = unique(mutation))  # Use unique levels only
  )
# Plot
ggplot(data, aes(fill = phenotypes, y = mutation, x = proportion)) +
  geom_bar(position = "fill", stat = "identity", width = 0.4) +
  scale_fill_manual(values = c("#99ccff", "#9966ff", "#808080", "#ff9933")) +
  theme(
    text = element_text(size = 12, family = "Open Sans"),
    axis.text.y = element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  ggtitle('MesIC under TGFB_L=0')
ggsave("MesICmut_no_input.png", width = 6, height = 5, device='png', dpi=700)

#################################### TGFBL and GF
data <- read.csv('specific_mutations_mesenchymalcells_TGFBL_GF.txt')

# Save original mutation order
mut_levels <- unique(data$mutation)

data <- data %>%
  mutate(
    mutation = as.character(mutation),  # Convert factor to character
    mutation = ifelse(mutation == "Wild type",
                      paste0("<b>", mutation, "</b>"),
                      mutation),
    mutation = factor(mutation, levels = unique(mutation))  # Use unique levels only
  )
# Plot
ggplot(data, aes(fill = phenotypes, y = mutation, x = proportion)) +
  geom_bar(position = "fill", stat = "identity", width = 0.4) +
  scale_fill_manual(values = c("#99ccff", "#9966ff", "#808080", "#ff9933")) +
  theme(
    text = element_text(size = 12, family = "Open Sans"),
    axis.text.y = element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  ggtitle('MesIC under TGFB_L=1 and GF=1')
ggsave("MesICmut_TGFBL_GF.png", width = 6, height = 5, device='png', dpi=700)

#################################### GF
data <- read.csv('specific_mutations_mesenchymalcells_GF.txt')

# Save original mutation order
mut_levels <- unique(data$mutation)

data <- data %>%
  mutate(
    mutation = as.character(mutation),  # Convert factor to character
    mutation = ifelse(mutation == "Wild type",
                      paste0("<b>", mutation, "</b>"),
                      mutation),
    mutation = factor(mutation, levels = unique(mutation))  # Use unique levels only
  )
# Plot
ggplot(data, aes(fill = phenotypes, y = mutation, x = proportion)) +
  geom_bar(position = "fill", stat = "identity", width = 0.4) +
  scale_fill_manual(values = c("#99ccff", "#9966ff", "#808080", "#ff9933")) +
  theme(
    text = element_text(size = 12, family = "Open Sans"),
    axis.text.y = element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  ggtitle('MesIC under GF=1')
ggsave("MesICmut_GF.png", width = 6, height = 5, device='png', dpi=700)

#################################################################################################
#################################################################################################

data <- read.csv('MesIC_MutationType_OFF.txt')

data_t <-data[which(data$phenotypes == "M"),]
data_t =data_t %>% arrange(proportion)
data_t[1]

data=data %>% arrange(factor(mutation, levels = data_t$mutation))
data$mutation <- factor(data$mutation, levels = unique(data$mutation))
data

ggplot(data, aes(fill=phenotypes, y=mutation, x=proportion)) + 
  geom_bar(position="fill", stat="identity", alpha=.6, width=.4) +
  scale_fill_manual(values = c("#99ccff", "#9966ff", "#808080", "#ff9933")) +
  theme(text=element_text(size=12,  family="Open Sans")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle('MesIC KO mutations under TGFB input')
ggsave("MesIC_KO.png", width = 6, height = 6, device='png', dpi=700)

#################################### PLOT IN MANUSCRIPT ####################################
data_t %>%
  mutate(mutation = fct_reorder(mutation, proportion)) %>%
  ggplot(aes(x=mutation, y=proportion, fill=tag)) +
  geom_bar(stat="identity", alpha=.6, width=.4) +
  scale_fill_manual(values = c("#f68060", 'red')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") +
  ggtitle("KO mutations on MesIC under TGFB_L input signal") +
  ylab("proportion of mesenchymal cells") + xlab("mutation (off)")
ggsave("MesIC_KO_TGFB_mesenchymal.png", width = 10, height = 6, device='png', dpi=700)


#################################################################################################
#################################################################################
data <- read.csv("modelaTGFB_mutations_intnodes.txt")
data_t <-data[which(data$phenotypes == "M"),]
aTGFB=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("modelaRTK_mutations_intnodes.txt")
data_t <-data[which(data$phenotypes == "M"),]
aRTK=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_intnodes.txt")
data_t <-data[which(data$phenotypes == "M"),]
original=data_t[which(data_t$proportion ==0),]$mutation

## common to all
RTK_TGFB_common=intersect(aRTK, aTGFB)
allmodels=intersect(RTK_TGFB_common, original)
write.csv(allmodels, "allmodels.csv", row.names=FALSE, quote = FALSE)

## in aRTK but not in aTGFB or original
TGFBandoriginal=c(aTGFB,original)
only_RTK=setdiff(aRTK, TGFBandoriginal)
write.csv(only_RTK, "only_aRTK.csv", row.names=FALSE, quote = FALSE)

## in aTGFB but not in RTK or original
RTKandoriginal=c(aRTK,original)
only_TGFB=setdiff(aTGFB, RTKandoriginal)
write.csv(only_TGFB, "only_aTGFB.csv", row.names=FALSE, quote = FALSE)


#################################### PLOT IN MANUSCRIPT (4) ####################################
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(aTGFB, aRTK, original),
  category.names = c("aTGFB only" , "aRTK only", "original"),
  filename = 'venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
#################################################################################################

### MUTATIONS ON ORIGINAL MODEL UNDER TGFBL, GF AND BOTH TGFBL AND GF ###

data <- read.csv("originalmodel_mutations_intnodes_TGFBL.txt")
data_t <-data[which(data$phenotypes == "M"),]
TGFB_L=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_intnodes_GF.txt")
data_t <-data[which(data$phenotypes == "M"),]
GF=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_intnodes_GF_TGFBL.txt")
data_t <-data[which(data$phenotypes == "M"),]
GF_and_TGFB_L=data_t[which(data_t$proportion ==0),]$mutation

input_list <- list(
  TGFB_L_input_original = TGFB_L,
  GF_input_original = GF,
  GF_and_TGFB_L_inputs_original = GF_and_TGFB_L,
  aTGFB = aTGFB,
  aRTK = aRTK,
  original = original
)

upset_data <- fromList(input_list)

png("intersection_inputs_no_inputs.png", width = 1000, height = 1000, res = 130)

ComplexUpset::upset(
  upset_data,
  intersect = names(input_list),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 5)
    )
  )
)
dev.off()

#################################################################################################

## COMBINATORIAL MUTATIONS 

data <- read.csv("modelaTGFB_mutations_mut_combination.txt")
data_t <-data[which(data$phenotypes == "M"),]
aTGFB=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("modelaRTK_mutations_mut_combination.txt")
data_t <-data[which(data$phenotypes == "M"),]
aRTK=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_mut_combination.txt")
data_t <-data[which(data$phenotypes == "M"),]
original=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_mut_combination_TGFBL.txt")
data_t <-data[which(data$phenotypes == "M"),]
TGFB_L=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_mut_combination_GF.txt")
data_t <-data[which(data$phenotypes == "M"),]
GF=data_t[which(data_t$proportion ==0),]$mutation

data <- read.csv("originalmodel_mutations_mut_combination_GF_TGFBL.txt")
data_t <-data[which(data$phenotypes == "M"),]
GF_and_TGFB_L=data_t[which(data_t$proportion ==0),]$mutation


input_list <- list(
  TGFB_L_input_original = TGFB_L,
  GF_input_original = GF,
  GF_and_TGFB_L_inputs_original = GF_and_TGFB_L,
  aTGFB = aTGFB,
  aRTK = aRTK,
  original = original
)

upset_data <- fromList(input_list)

# upset plot

png("combinatorial.png", width = 1000, height = 1000, res = 130)

ComplexUpset::upset(
  upset_data,
  intersect = names(input_list),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 5)
    )
  )
)

dev.off()

# heatmap

m=fromList(input_list)
m_tbl <- m %>%
  tibble::rownames_to_column("element")

df_long <- stack(input_list)

colnames(df_long) <- c("element", "set")

m_df <- df_long %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = set,
    values_from = value,
    values_fill = 0
    
  )


# heatmap 

m_df=as.data.frame(m_df)
rownames(m_df)<-m_df$element 
df=(data.matrix(subset(m_df, select = -element)))

rownames(df) <- rownames(df) %>%
  str_remove_all("^\\('") %>%     # remove starting ('
  str_remove_all("'\\)$") %>%     # remove ending ')
  str_replace("', '", " & ")  

png("combinatorial_heatmap.png", width = 1600, height = 1400, res = 130)

Heatmap(
  df,
  name = "val",
  col = c("0" = "white", "1" = "steelblue"),
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # add separation between squares
  rect_gp = gpar(col = "grey70", lwd = 1),
  
  # rotate column names
  column_names_rot = 45,
  
  # optional: improve readability
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 8),
  
)

dev.off()
