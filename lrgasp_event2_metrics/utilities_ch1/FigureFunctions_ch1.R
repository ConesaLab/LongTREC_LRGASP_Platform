#######################################
######### FIGURES FUNCTIONS CHALLENGE 1 #########
#######################################
## author: Francisco J. Pardo-Palacios, f.pardo.palacios@gmail.com
## author:Ana Conesa, ana.conesa@csic.es
## author: Wouter Maessen, woutermaessen2511@gmail.com
## Last modified: April 25th 2023
#######################################
list_of_packages <- c("ggplot2", "tidyverse", "ggpubr", "scales", "patchwork", "gridExtra", "grid", "RColorConesa", "fmsb", "MetBrewer","plotly",
                      "DT","knitr","optparse",'rmarkdown',"ggridges","reshape","ComplexHeatmap","jaccard","corrplot","huxtable","hrbrthemes","gt",
                      "ggpmisc","UpSetR","ComplexUpset","ggplot2movies","ggbreak","gg.gap","ggcorrplot")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

# LIBRARIES
#################
library(ggplot2)
library(ggrepel)
library(patchwork)
library(grid)
library(MetBrewer)
library(RColorConesa)
library(ggpubr)
library(scales)
library(huxtable)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(gridExtra)
library(reshape2)
library(gt)
library(ggpmisc)
library(UpSetR)
library(ComplexUpset)
library(ggplot2movies)
library(tidyr)
library(grid)
library(ggbreak) 
library(gg.gap)
library(ggcorrplot)
library(viridis)
library(fmsb)
library(plotly)
library(openxlsx)

outdir = "output/main"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

setwd(utilities.path)
setwd("/Users/woutermaessen/PycharmProjects/lrgasp_platform/lrgasp_event2_metrics/utilities_ch1")

# PATTERNS
#################

pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")

pub_theme_flip <- theme_pubclean(base_family = "Helvetica", flip=T) +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")


old.libplat.palette = c( "cDNA-PacBio"="#e66adb", "CapTrap-PacBio"="#ab0202", "cDNA-Illumina"="#FFCF71",  "Freestyle-Freestyle"="#75b562",
                         "cDNA-ONT"="#005C75", "CapTrap-ONT"="#7482F0", "R2C2-ONT"="#74CDF0", "dRNA-ONT"="#1b36d1"
)

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)

cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

# FUNCTIONS
#################
source("Functions_Supplementary_Figures_Challenge1_v4.R")

# FIGURES
#################

# Figure 2a. Detection
#########################
fig_2a_challenge1 <- function(new_data, 
                              new_data2 = NULL, 
                              data_sample = "WTC11", 
                              outdir = "./", 
                              ylims = c(160000, 250000), 
                              xlims = c(18000, 22000),
                              comparison,
                              bambu,
                              FLAIR,
                              Lyric,
                              IsoTools,
                              Mandalorion,
                              Iso_IB,
                              FLAMES,
                              IsoQuant,
                              Spectra,
                              TALON_LAPA,
                              StringTie2) {
  setwd(utilities.path)
  data_sample2 <- paste0("Challenge1_Figures_Data/", data_sample, "_results/", data_sample)
  
  genes_file <- paste0(data_sample2,".genes_SJ_table.csv")
  genes_SJ <- read.csv(genes_file, header = T, sep = ",")%>% t()
  genes_SJ <- rbind(genes_SJ, new_data[1:6])
  rownames(genes_SJ)[nrow(genes_SJ)] <- new_data[7]
  if(!is.null(new_data2)){
    genes_SJ <- rbind(genes_SJ, new_data2[1:6])
    rownames(genes_SJ)[nrow(genes_SJ)] <- new_data2[7]
  }
  
  num_trx_file <- paste0(data_sample2, ".summary_table_SC.csv") 
  num_trx <- read.csv(num_trx_file, header = T, sep = ",")
  num_trx <- rbind(num_trx, new_data[7:17])
  if(!is.null(new_data2)){
    num_trx <- rbind(num_trx, new_data2[7:17])
  }
  num_trx <- num_trx %>%
    mutate(across(2:11, as.integer))
  
  code_file <- paste0( data_sample2, ".code_updated.txt")
  code=read.csv(code_file, header = T, sep=",")
  
  code$Lib_Plat <- apply(code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  code$Lib_DC=apply(cbind(code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  code$Label <-apply(cbind(code[,c("Platform","Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  
  code <- rbind(code, meta_data_1)
  if(!is.null(new_data2)){
    code <- rbind(code, meta_data_2)
  }
  
  big_df <- merge(num_trx, genes_SJ, by.x="ID",by.y=0)
  
  big_df$FSM_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["FSM"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$ISM_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["ISM"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$NIC_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["NIC"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$NNC_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["NNC"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$knowTrx <- apply(big_df, 1, function(x){
    as.numeric(x["FSM"])+as.numeric(x["ISM"])+as.numeric(x["NIC"])+as.numeric(x["NNC"])
  })
  
  big_df$avgTrxGene <- apply(big_df, 1, function(x){
    as.numeric(x["knowTrx"])/as.numeric(x["Num.KnGenes"])
  })
  
  
  big_df <- merge(big_df, code, by.x="ID", by.y="pipelineCode")
  
  ## First we need to read all the data. 
  fsm_metrics <- read.csv(paste0(data_sample2,".FSM_metrics.csv"), sep=",", header = T) %>% t()
  ism_metrics <- read.csv(paste0(data_sample2,".ISM_metrics.csv"), sep=",", header = T) %>% t()
  nic_metrics <- read.csv(paste0(data_sample2,".NIC_metrics.csv"), sep=",", header=T) %>% t()
  nnc_metrics <- read.csv(paste0(data_sample2,".NNC_metrics.csv"), sep=",", header=T) %>% t()
  
  fsm_metrics <- rbind(fsm_metrics, temp_env$all.results$FSM$`Absolute value`)
  rownames(fsm_metrics)[nrow(fsm_metrics)] <- pipelineCode
  ism_metrics <- rbind(ism_metrics, temp_env$all.results$ISM$`Absolute value`)
  rownames(ism_metrics)[nrow(ism_metrics)] <- pipelineCode
  nic_metrics <- rbind(nic_metrics, temp_env$all.results$NIC$`Absolute value`)
  rownames(nic_metrics)[nrow(nic_metrics)] <- pipelineCode
  nnc_metrics <- rbind(nnc_metrics, temp_env$all.results$NNC$`Absolute value`)
  rownames(nnc_metrics)[nrow(nnc_metrics)] <- pipelineCode
  
  if(!is.null(new_data2)){
    fsm_metrics <- rbind(fsm_metrics, temp_env2$all.results$FSM$`Absolute value`)
    rownames(fsm_metrics)[nrow(fsm_metrics)] <- pipelineCode2
    ism_metrics <- rbind(ism_metrics, temp_env2$all.results$ISM$`Absolute value`)
    rownames(ism_metrics)[nrow(ism_metrics)] <- pipelineCode2
    nic_metrics <- rbind(nic_metrics, temp_env2$all.results$NIC$`Absolute value`)
    rownames(nic_metrics)[nrow(nic_metrics)] <- pipelineCode2
    nnc_metrics <- rbind(nnc_metrics, temp_env2$all.results$NNC$`Absolute value`)
    rownames(nnc_metrics)[nrow(nnc_metrics)] <- pipelineCode2
  }
  
  support_df <- data.frame(FSM_5ref=fsm_metrics[,"5' reference supported (gene)"], FSM_5illumina=fsm_metrics[,"5' CAGE supported"], 
                           FSM_3ref=fsm_metrics[,"3' reference supported (gene)"], FSM_3illumina=fsm_metrics[,"3' QuantSeq supported"],
                           ISM_5ref=ism_metrics[,"5' reference supported (gene)"], ISM_5illumina=ism_metrics[,"5' CAGE supported"], 
                           ISM_3ref=ism_metrics[,"3' reference supported (gene)"], ISM_3illumina=ism_metrics[,"3' QuantSeq supported"],
                           NIC_5ref=nic_metrics[,"5' reference supported (gene)"], NIC_5illumina=nic_metrics[,"5' CAGE supported"], 
                           NIC_3ref=nic_metrics[,"3' reference supported (gene)"], NIC_3illumina=nic_metrics[,"3' QuantSeq supported"],
                           NNC_5ref=nnc_metrics[,"5' reference supported (gene)"], NNC_5illumina=nnc_metrics[,"5' CAGE supported"], 
                           NNC_3ref=nnc_metrics[,"3' reference supported (gene)"], NNC_3illumina=nnc_metrics[,"3' QuantSeq supported"])
  support_df$total_5ref <- apply(support_df,1, function(x){
    as.numeric(x["FSM_5ref"])+as.numeric(x["ISM_5ref"])+as.numeric(x["NIC_5ref"])+as.numeric(x["NNC_5ref"])
  })
  
  support_df$total_3ref <- apply(support_df,1, function(x){
    as.numeric(x["FSM_3ref"])+as.numeric(x["ISM_3ref"])+as.numeric(x["NIC_3ref"])+as.numeric(x["NNC_3ref"])
  })
  
  support_df$total_5illumina <- apply(support_df,1, function(x){
    as.numeric(x["FSM_5illumina"])+as.numeric(x["ISM_5illumina"])+as.numeric(x["NIC_5illumina"])+as.numeric(x["NNC_5illumina"])
  })
  
  support_df$total_3illumina <- apply(support_df,1, function(x){
    as.numeric(x["FSM_3illumina"])+as.numeric(x["ISM_3illumina"])+as.numeric(x["NIC_3illumina"])+as.numeric(x["NNC_3illumina"])
  })
  
  support_df <- merge(support_df, big_df[,c("ID", "knowTrx")], by.x=0, by.y="ID")
  
  support_df$perc_5ref <- apply(support_df,1,function(x){
    as.numeric(x["total_5ref"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df$perc_5illumina <- apply(support_df,1,function(x){
    as.numeric(x["total_5illumina"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df$perc_3ref <- apply(support_df,1,function(x){
    as.numeric(x["total_3ref"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df$perc_3illumina <- apply(support_df,1,function(x){
    as.numeric(x["total_3illumina"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df <- merge(support_df,code, by.x="Row.names", by.y="pipelineCode")
  support_df$opacity <- ifelse(support_df$Row.names %in% c(pipelineCode, pipelineCode2), 1, 0.3)
  
  if(comparison == 'Custom'){
    if(bambu == 'NA'){
      support_df <- support_df[support_df$Tool != 'Bambu',]
    }
    if(FLAIR == 'NA'){
      support_df <- support_df[support_df$Tool != 'FLAIR',]
    }
    if(Lyric == 'NA'){
      support_df <- support_df[!support_df$Tool %in% c('Lyric*','Lyric'),]
    }
    if(IsoTools == 'NA'){
      support_df <- support_df[support_df$Tool != 'IsoTools',]
    }
    if(Mandalorion == 'NA'){
      support_df <- support_df[support_df$Tool != 'Mandalorion',]
    }
    if(Iso_IB == 'NA'){
      support_df <- support_df[support_df$Tool != 'Iso_IB',]
    }
    if(FLAMES == 'NA'){
      support_df <- support_df[support_df$Tool != 'FLAMES',]
    }
    if(IsoQuant == 'NA'){
      support_df <- support_df[support_df$Tool != 'IsoQuant',]
    }
    if(Spectra == 'NA'){
      support_df <- support_df[support_df$Tool != 'Spectra',]
    }
    if(TALON_LAPA == 'NA'){
      support_df <- support_df[support_df$Tool != 'TALON_LAPA',]
    }
    if(StringTie2 == 'NA'){
      support_df <- support_df[support_df$Tool != 'StringTie2',]
    }
  }
  
  p5 <- ggplot(support_df, aes(x=perc_5ref, y=perc_5illumina, color=Lib_Plat, alpha=opacity)) +
    geom_point(aes(size=knowTrx, shape=Data_Category)) +
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1))) +
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1))) +
    scale_color_manual(values = libplat.palette) +
    scale_alpha_identity() +
    pub_theme +
    xlab("5' end supported by reference") + 
    ylab("5' end supported by CAGE") +
    theme(axis.text.y = element_text(size=18),
          axis.text.x = element_text(size=18),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  p6 <- ggplot(support_df, aes(x=perc_3ref, y=perc_3illumina, color=Lib_Plat, alpha=opacity))+
    geom_point(aes(size=knowTrx, shape=Data_Category) )+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    scale_alpha_identity() +
    xlab("3' end supported by reference") + ylab("3' end supported by QuantSeq") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  sj_info <- read.csv(paste0(data_sample2,".SJ_info.csv"), sep=",", header=T)
  sj_info <- rbind(sj_info, c(new_data[7], 
                              new_data[18], 
                              new_data[19], 
                              new_data[20], 
                              new_data[21],
                              new_data[22],
                              new_data[23]))
  
  if(!is.null(new_data2)){
    sj_info <- rbind(sj_info, c(new_data2[7], 
                                new_data2[18], 
                                new_data2[19], 
                                new_data2[20], 
                                new_data2[21],
                                new_data2[22],
                                new_data2[23]))
  }
  
  sj_info <- merge(sj_info,code, by="pipelineCode")
  sj_info$perc_known <- apply(sj_info,1,function(x){
    as.numeric(x["known"])*100/as.numeric(x["total"])
  })
  sj_info$perc_canonical <- apply(sj_info,1,function(x){
    as.numeric(x["canonical"])*100/as.numeric(x["total"])
  })
  sj_info$perc_cov <- apply(sj_info,1,function(x){
    as.numeric(x["with_Coverage"])*100/as.numeric(x["total"])
  })
  
  sj_info$perc_RTS <- apply(sj_info,1,function(x){
    as.numeric(x["RTS"])*100/as.numeric(x["total"])
  })
  
  sj_info$opacity <- ifelse(sj_info$pipelineCode %in% c(pipelineCode,pipelineCode2), 1, 0.3)
  
  if(comparison == 'Custom'){
    if(bambu == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'Bambu',]
    }
    if(FLAIR == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'FLAIR',]
    }
    if(Lyric == 'NA'){
      sj_info <- sj_info[!sj_info$Tool %in% c('Lyric*','Lyric'),]
    }
    if(IsoTools == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'IsoTools',]
    }
    if(Mandalorion == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'Mandalorion',]
    }
    if(Iso_IB == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'Iso_IB',]
    }
    if(FLAMES == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'FLAMES',]
    }
    if(IsoQuant == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'IsoQuant',]
    }
    if(Spectra == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'Spectra',]
    }
    if(TALON_LAPA == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'TALON_LAPA',]
    }
    if(StringTie2 == 'NA'){
      sj_info <- sj_info[sj_info$Tool != 'StringTie2',]
    }
  }
  
  p8 <- ggplot(sj_info, aes(x=perc_canonical, y=perc_cov, color=Lib_Plat, alpha=opacity))+
    geom_point(aes(size=as.numeric(total), shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    scale_alpha_identity() +
    xlab("% canonical SJ detected") + ylab("% SJ with Short-Read coverage") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  #### zoom p8
  p8_zoom <- p8 + xlab("") + ylab("") +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)), limits = c(92,100))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)), limits = c(98.5,100))
  
  p2a <- grid.arrange(p5, p6, p8, p8_zoom, ncol = 2, nrow = 2)
  grid.rect(x = unit(0.45, "npc"),  # X position (adjust as needed)
            y = unit(0.415, "npc"),  # Y position (adjust as needed)
            width = unit(0.05, "npc"),  # Rectangle width
            height = unit(0.08, "npc"), # Rectangle height
            gp = gpar(fill = NA, col = "black", lwd = 2, lty = "solid"))  # Dashed border, transparent fill
  grid.lines(x = unit(c(0.475, 0.5675), "npc"), 
             y = unit(c(0.455, 0.48125), "npc"),  # Slight upward tilt
             gp = gpar(col = "black", lwd = 2, lty = "dashed"))
  grid.lines(x = unit(c(0.475, 0.568), "npc"), 
             y = unit(c(0.375, 0.075), "npc"),  # Slight upward tilt
             gp = gpar(col = "black", lwd = 2, lty = "dashed"))
  
  
  return(p2a)
  

}

# Figure 2c. SIRVs
#########################

fig_2c_challenge1 <- function(data_sample = "WTC11", outdir = outdir){
  setwd(utilities.path)
  data_sample2 <-paste0("Challenge1_Figures_Data/", data_sample, "_results/", data_sample)
  
  code_file <- paste0(data_sample2, ".code_updated.txt")
  code=read.csv(code_file, header = T, sep=",")
  
  code$Lib_Plat <- apply(code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  code$Lib_DC=apply(cbind(code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  code$Label <-apply(cbind(code[,c("Platform","Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  
  code <- rbind(code, meta_data_1)
  if(dataset2){
    code <- rbind(code, meta_data_2)
  }
  
  spliced_SIRV_metrics <- read.csv(paste0(data_sample2,".splicedSIRVS_metrics.csv"), sep=",", header=T) %>% t() %>% as.data.frame()
  spliced_SIRV_metrics <- rbind(spliced_SIRV_metrics, temp_env$all.results$spliced_SIRV$Value)
  rownames(spliced_SIRV_metrics)[nrow(spliced_SIRV_metrics)] <- pipelineCode
  if(dataset2){
    spliced_SIRV_metrics <- rbind(spliced_SIRV_metrics, temp_env2$all.results$spliced_SIRV$Value)
    rownames(spliced_SIRV_metrics)[nrow(spliced_SIRV_metrics)] <- pipelineCode2
  }
  spliced_SIRV_metrics[] <- lapply(spliced_SIRV_metrics, as.numeric)
  
  spliced_SIRV_metrics <- merge(spliced_SIRV_metrics, code, by.x=0, by.y="pipelineCode")
  spliced_SIRV_metrics[, "F1 score"] <- apply(spliced_SIRV_metrics, 1, function(x){
    s=as.numeric(x["Sensitivity"])
    p=as.numeric(x["Precision"])
    (2*s*p)/(s+p)
  })
  pivoted_spliced <- pivot_longer(spliced_SIRV_metrics,
                                  cols=c("Sensitivity","Precision", "F1 score"),
                                  names_to = "metrics",
                                  values_to = "value")
  
  pivoted_spliced$metrics  <- pivoted_spliced$metrics %>% factor(levels = c("Sensitivity","Precision", "F1 score"),
                                                                 labels = c("Sensitivity","Precision", "F1 score"))
  
  pivoted_spliced$full_Label <- paste0(pivoted_spliced$Label,"-", pivoted_spliced$Data_Category )
  
  if(comparison == 'Custom'){
    if(bambu == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'Bambu',]
    }
    if(FLAIR == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'FLAIR',]
    }
    if(Lyric == 'NA'){
      pivoted_spliced <- pivoted_spliced[!pivoted_spliced$Tool %in% c('Lyric*','Lyric'),]
    }
    if(IsoTools == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'IsoTools',]
    }
    if(Mandalorion == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'Mandalorion',]
    }
    if(Iso_IB == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'Iso_IB',]
    }
    if(FLAMES == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'FLAMES',]
    }
    if(IsoQuant == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'IsoQuant',]
    }
    if(Spectra == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'Spectra',]
    }
    if(TALON_LAPA == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'TALON_LAPA',]
    }
    if(StringTie2 == 'NA'){
      pivoted_spliced <- pivoted_spliced[pivoted_spliced$Tool != 'StringTie2',]
    }
  }
  
  pivoted_spliced$Alias <- factor(pivoted_spliced$Alias, levels = c(sort(setdiff(unique(pivoted_spliced$Alias), alias)), alias))
  
  
  pSIRVspliced <- ggplot(pivoted_spliced, aes(x=full_Label, y=value)) +
    geom_segment( aes(x=full_Label, xend=full_Label, y=0, yend=value, color=Lib_Plat), size=0.8) +
    geom_point( size=2, aes( shape=Data_Category, color=Lib_Plat ))  +
    facet_grid( metrics ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
    pub_theme +
    scale_color_manual(values = libplat.palette) +
    labs(x="", y="") +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, NA), position="right")
  
  #### same for unspliced
  
  unspliced_SIRV_metrics <- read.csv(paste0(data_sample2,".unsplicedSIRVS_metrics.csv"), sep=",", header=T) %>% t() %>% as.data.frame()
  unspliced_SIRV_metrics <- rbind(unspliced_SIRV_metrics, temp_env$all.results$unspliced_SIRV$Value)
  rownames(unspliced_SIRV_metrics)[nrow(unspliced_SIRV_metrics)] <- pipelineCode
  if(dataset2){
    unspliced_SIRV_metrics <- rbind(unspliced_SIRV_metrics, temp_env2$all.results$unspliced_SIRV$Value)
    rownames(unspliced_SIRV_metrics)[nrow(unspliced_SIRV_metrics)] <- pipelineCode2
  }
  unspliced_SIRV_metrics[] <- lapply(unspliced_SIRV_metrics, as.numeric)
  
  unspliced_SIRV_metrics <- merge(unspliced_SIRV_metrics, code, by.x=0, by.y="pipelineCode")
  
  unspliced_SIRV_metrics[, "F1 score"] <- apply(unspliced_SIRV_metrics, 1, function(x){
    s=as.numeric(x["Sensitivity"])
    p=as.numeric(x["Precision"])
    (2*s*p)/(s+p)
  })
  
  pivoted_unspliced <- pivot_longer(unspliced_SIRV_metrics,
                                    cols=c("Sensitivity","Precision", "F1 score"),
                                    names_to = "metrics",
                                    values_to = "value")
  pivoted_unspliced$metrics  <- pivoted_unspliced$metrics %>% factor(levels = c("Sensitivity","Precision", "F1 score"),
                                                                     labels = c("Sensitivity","Precision", "F1 score"))
  
  pivoted_unspliced$full_Label <- paste0(pivoted_unspliced$Label,"-", pivoted_unspliced$Data_Category )
  
  if(comparison == 'Custom'){
    if(bambu == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'Bambu',]
    }
    if(FLAIR == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'FLAIR',]
    }
    if(Lyric == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[!pivoted_unspliced$Tool %in% c('Lyric*','Lyric'),]
    }
    if(IsoTools == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'IsoTools',]
    }
    if(Mandalorion == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'Mandalorion',]
    }
    if(Iso_IB == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'Iso_IB',]
    }
    if(FLAMES == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'FLAMES',]
    }
    if(IsoQuant == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'IsoQuant',]
    }
    if(Spectra == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'Spectra',]
    }
    if(TALON_LAPA == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'TALON_LAPA',]
    }
    if(StringTie2 == 'NA'){
      pivoted_unspliced <- pivoted_unspliced[pivoted_unspliced$Tool != 'StringTie2',]
    }
  }
  
  pivoted_unspliced$Alias <- factor(pivoted_unspliced$Alias, levels = c(sort(setdiff(unique(pivoted_unspliced$Alias), alias)), alias))
  
  pSIRVunspliced <- ggplot(pivoted_unspliced, aes(x=full_Label, y=value)) +
    geom_segment( aes(x=full_Label, xend=full_Label, y=0, yend=value, color=Lib_Plat), size=0.8) +
    geom_point( size=2, aes( shape=Data_Category, color=Lib_Plat ))  +
    facet_grid( metrics ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
    pub_theme +
    scale_color_manual(values = libplat.palette) +
    labs(x="", y="") +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, NA), position="right")
  
  p2c <- grid.arrange(pSIRVspliced, pSIRVunspliced, ncol = 1, nrow = 2)
  
  return(p2c)
}



# Figure 2d. Radarplot
#########################
radar.simulation (species = "human", directory = "Challenge1_Figures_Data/Simulations/", pdf = paste0(outdir, "/Fig2d"))

radar.simulation (species = "mouse", directory = "Challenge1_Figures_Data/Simulations/", pdf = paste0(outdir, "/Fig2d_mouse"))


## Figure 2e. Evaluation against GENCODE
########################################
fig_2e_challenge1(new_data = c(temp_env$sqanti_data))

fig_2e_challenge1 <- function(new_data){
  setwd(utilities.path)
  pa_GENCODE <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/presence_absence.GENCODE_loci_2.csv", sep=",", header = T) [,1:51] # Presence absence analysis of all transcripits of the 50 loci in pipelines evaluated against manual annotation. 
  GENCODE_new <- pa_GENCODE$Row.names %in% new_data$LRGASP_id_2
  pa_GENCODE <- cbind(pa_GENCODE, as.numeric(GENCODE_new))
  colnames(pa_GENCODE)[ncol(pa_GENCODE)] <- pipelineCode
  
  pa.WTC11 <- read.csv("Challenge1_Figures_Data/WTC11_results/WTC11_comparison.pa.csv", as.is = TRUE)
  WTC11_new <- pa.WTC11$TAGS %in% new_data$LRGASP_id
  pa.WTC11 <- cbind(pa.WTC11, as.numeric(WTC11_new))
  colnames(pa.WTC11)[ncol(pa.WTC11)] <- pipelineCode
  
  gencode_eval_results <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/new_GENCODE_manualAnnot_evaluation.csv", header = T) # evaluation result
  gencode_genes <- read.xlsx("Challenge1_Figures_Data/GENCODE_manualAnnot/gencode_50_genes.xlsx", startRow = 2)
  colnames(gencode_genes) <- c('human_gene', 'human_gene_id', 'mouse_gene', 'mouse_gene_id', 'Ortholog')
  
  code_file <- paste0(data_sample2, ".code_updated.txt")
  code=read.csv(code_file, header = T, sep=",")
  
  code$Lib_Plat <- apply(code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  code$Lib_DC=apply(cbind(code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  code$Label <-apply(cbind(code[,c("Platform","Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  
  code <- rbind(code, meta_data_1)
  
  genocode_eval_WTC11 <- performance.genecode (gencode.pa = pa_GENCODE, ID_UIC = NULL,
                                               pa = pa.WTC11,  code = code, 
                                               selection = NULL, evaluation = gencode_eval_results,
                                               mypattern = "SQ3_human",
                                               directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/human/")
  
  pivoted_gencode_known <- pivot_longer(genocode_eval_WTC11 , cols = c("Sensitivity_known", "Precision_known", "F1_known"))
  pivoted_gencode_known$name <- pivoted_gencode_known$name %>% factor(levels = c("Sensitivity_known", "Precision_known", "F1_known"),
                                                                      labels = c("Sensitivity", "Precision", "F1-score"))
  
  pivoted_gencode_novel <- pivot_longer(genocode_eval_WTC11 , cols = c("Sensitivity_novel", "Precision_novel", "F1_novel"))
  pivoted_gencode_novel$name <- pivoted_gencode_novel$name %>% factor(levels = c("Sensitivity_novel", "Precision_novel", "F1_novel"),
                                                                      labels = c("Sensitivity", "Precision", "F1-score"))
  
  pC.known.mid <- ggplot(pivoted_gencode_known, aes(x=Label, y=value)) +
    geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Sample_code), size=0.8) +
    geom_point( size=2, aes( shape=Data_Category, color=Sample_code ))  +
    facet_grid( name ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
    pub_theme +
    scale_color_manual(values = libplat.palette) +
    labs(x="", y="",
         title="Known transcript level") +
    theme(legend.position="bottom") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y = element_text(size = 10))+
    scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, 1), position="right")
  
  pC.novel <- ggplot(pivoted_gencode_novel, aes(x=Label, y=value)) +
    geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Sample_code), size=0.8) +
    geom_point( size=2, aes( shape=Data_Category, color=Sample_code ))  +
    facet_grid( name ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
    pub_theme +
    scale_color_manual(values = libplat.palette) +
    labs(x="", y="",
         title="Novel transcript level") +
    theme(legend.position="bottom") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y = element_text(size = 10))+
    scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, 1), position="right")
  
  figure2d <- ggarrange(pC.known.mid,
                        pC.novel,
                        ncol = 2, nrow = 1, common.legend = TRUE)
}



###############
### Summary figure

realData_simulation_metrics <- get_realData_simulation_metrics(data_sample = "WTC11",
                                                               outdir = outdir, species = "human",
                                                               sim_directory = "Challenge1_Figures_Data/Simulations/")

summary_GENCODE_eval <- genocode_eval_WTC11 %>% select(Row.names, Sensitivity_known, Precision_known, 
                                                       Sensitivity_novel, Precision_novel, Total_detections)

df_summary_metrics <- merge(realData_simulation_metrics, summary_GENCODE_eval, by="Row.names")

summary_SIRVs <- spliced_SIRV_metrics %>% select(Row.names, Sensitivity, Precision)

df_summary_metrics <- merge(df_summary_metrics, summary_SIRVs, by="Row.names")

df_summary_metrics$perc_3illumina <- df_summary_metrics$perc_3illumina/100
df_summary_metrics$perc_5illumina <- df_summary_metrics$perc_5illumina/100
df_summary_metrics$perc_SRTM <- df_summary_metrics$perc_SRTM/100
df_summary_metrics$perc_SNTM <- df_summary_metrics$perc_SNTM/100
df_summary_metrics$perc_cov <- df_summary_metrics$perc_cov/100

df_pivoted_summary <- pivot_longer(df_summary_metrics, 
                                   cols=c("perc_SRTM", "perc_SNTM", "perc_5illumina", "perc_3illumina","perc_cov",
                                          "Sen_kn","Pre_kn","Sen_no", "Pre_no",
                                          "Sensitivity_known", "Precision_known",
                                          "Sensitivity_novel", "Precision_novel",
                                          "Sensitivity", "Precision"))

df_pivoted_summary <- df_pivoted_summary %>%
  mutate(value=ifelse(value=="NaN",
                      0,
                      value)) %>% 
  group_by(name, Lib_Plat) %>%
  mutate(Ranking = rank(value, na.last = "keep")) %>% 
  mutate(max_rank=max(Ranking, na.rm = T))

df_pivoted_summary$Ranking_adj <- df_pivoted_summary$Ranking*11/df_pivoted_summary$max_rank


df_pivoted_summary$name <-df_pivoted_summary$name %>% 
  factor(levels=rev(c("perc_SRTM","perc_SNTM","perc_5illumina","perc_3illumina","perc_cov",
                      "Sensitivity","Precision",
                      "Sensitivity_known","Precision_known",
                      "Sensitivity_novel","Precision_novel",
                      "Sen_kn", "Pre_kn","Sen_no","Pre_no")),
         labels=rev(c("% SRTM", "% SNTM", "% CAGE-Seq", "% Quant-Seq","% SJ cov.",
                      "SIRV Sensitivity", "SIRV Precision",
                      "GENCODE Sensit. (known)", "GENCODE Prec. (known)",
                      "GENCODE Sensit. (novel)", "GENCODE Prec. (novel)",
                      "SIMULATION Sensit. (known)", "SIMULATION Prec. (known)",
                      "SIMULATION Sensit. (novel)", "SIMULATION Prec. (novel)")))

df_pivoted_summary <- df_pivoted_summary %>% mutate(type_metric=ifelse(
  name %in% c("% SRTM", "% SNTM", "% CAGE-Seq", "% Quant-Seq", "% SJ cov."),
  "Real data",
  ifelse(
    name %in% c("SIRV Sensitivity", "SIRV Precision"),
    "SIRVs",
    ifelse(
      name %in% c("GENCODE Sensit. (known)", "GENCODE Prec. (known)",
                  "GENCODE Sensit. (novel)", "GENCODE Prec. (novel)"),
      "GENCODE manual annot",
      "Simulation"
    )
  )
)) 


top_labels <- c("Rest", "Bronze", "Silver", "Gold")
top_shapes <- c("Rest","Top3","Top2", "Top1")
df_pivoted_summary$Rank_top <- cut(df_pivoted_summary$Ranking_adj,
                                   breaks = c(0,8,9,10,11),
                                   labels = top_shapes,
                                   include.lowest = T)

df_pivoted_summary$Rank_top <- df_pivoted_summary$Rank_top %>% factor(levels=rev(top_shapes))

custom_palette <- colorRampPalette(c("#F5962A", "#ECF52A"))
top <- colorRampPalette(c("#B82D2D", "#F12626"))(3)
bottom <- colorRampPalette(c( "#8AE826", "#5CAE05"))(3)

paleta_final <- c(top, custom_palette(10), bottom)

palette_top <- c("Gold"="#FDDA3A",
                 "Silver"="#DAD8D2",
                 "Bronze"="#CA8B2B",
                 "Rest"="white")
shapes_top <- c("Top1"=24,
                "Top2"=23,
                "Top3"=22,
                "Rest"=21)

pdf(paste0(outdir,"/summary_challenge1_metrics.reorder_shapes.pdf"), width=9, height = 6)
for (i in c("cDNA-PacBio", "cDNA-ONT",
            "CapTrap-PacBio", "CapTrap-ONT",
            "R2C2-ONT", "dRNA-ONT")){
  
  if (startsWith(i, "cDNA")){
    sections <- c(4.5, 8.5, 10.5)
  }else{
    sections <- c(4.5, 6.5)
  }
  
  p.traffic_Rank <- ggplot(df_pivoted_summary  %>% filter(Lib_Plat==i) %>% na.omit(),
                           aes(y=name, x=Data_Category))+
    geom_point(size = 4, stroke=0.4, aes(shape=Rank_top, fill=Ranking_adj)) +  
    geom_hline(yintercept = sections,
               linetype="dashed", color="#0578C2") + 
    facet_grid(Lib_Plat~Alias, drop = T, scales="free")+
    scale_fill_viridis(breaks = c(2, 10), labels = c("Bottom", "Top"))+
    #scale_fill_gradientn(colors = paleta_final, values = rescale(c(1:12)) ,
    #                      na.value = "grey50",
    #                      breaks = c(2, 10), labels = c("Bottom", "Top"))  +
    scale_shape_manual(values = shapes_top)+
    pub_theme+
    theme(axis.text.x = element_text(angle=0, size=8),
          axis.text.y = element_text(size=8),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    labs(colour="Ranking\n")
  
  
  p.traffic_Value <- ggplot(df_pivoted_summary  %>% filter(Lib_Plat==i)%>% na.omit(),
                            aes(y=name, x=Data_Category))+
    geom_point(size = 4, stroke=0.4, aes(shape=Rank_top, fill=value)) +  
    geom_hline(yintercept = sections,
               linetype="dashed", color="#0578C2") +
    facet_grid(Lib_Plat~Alias, drop = T, scales="free")+
    scale_fill_viridis(breaks = c(0.1, 0.5, 0.9))  +
    #scale_fill_gradientn(colors = paleta_final, values = rescale(c(0:100)) ,
    #                      na.value = "grey50",
    #                      breaks = c(0.1, 0.5, 0.9))  +
    scale_shape_manual(values=shapes_top)+
    pub_theme+
    theme(axis.text.x = element_text(angle=0, size=8),
          axis.text.y = element_text(size=8),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    labs(colour="Value\n")
  
  p.bar_summary <- ggplot(df_summary_metrics  %>% filter(Lib_Plat==i),
                          aes(y=total, x=Data_Category, fill=Tool))+
    geom_bar(stat="identity")+
    geom_text(aes(label=paste(round(total*0.001,digits = 0), "K")),
              vjust=0, size=3) +
    facet_grid(Lib_Plat~Alias, drop = T, scales="free")+
    scale_color_gradientn(colors = paleta_final, values = rescale(c(0:100)) ,
                          na.value = "grey50",
                          breaks = c(0.1, 0.5, 0.9))  +
    scale_fill_met_d("Cross")+
    pub_theme+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          strip.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = -20) ) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1),
                       expand = expansion(mult=c(0,0.2)))+
    labs(y="Num. transcripts")
  
  sum1 <- p.bar_summary / p.traffic_Value + 
    plot_layout(heights = c(1, 3))
  sum2 <- p.bar_summary / p.traffic_Rank + 
    plot_layout(heights = c(1, 3))
  
  #ggsave(file=paste0(outdir, "/summary_figure.rank.",i,".svg"), plot=sum2, width=9, height=6)
  #ggsave(file=paste0(outdir, "/summary_figure.value.",i,".svg"), plot=sum1, width=9, height=6)
  
  
  print(sum1)
  print(sum2)
  
  fig2e <- ggarrange(sum1,
            sum2,
            ncol = 2, nrow = 1, common.legend = TRUE)
  
  return(fig2e)
  
}
dev.off()
