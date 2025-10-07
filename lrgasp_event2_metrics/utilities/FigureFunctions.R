################################################################################
### Script: Functions for figure consolidation #################################
### Author: Wouter Maessen, Ana Conesa Lab #####################################
### Date: 14-10-2024 ###########################################################
################################################################################

list_of_packages <- c("ggplot2", "tidyverse", "ggpubr", "scales", "patchwork", "gridExtra", "grid", "RColorConesa", "fmsb", "MetBrewer","plotly")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")


# Install and load packages
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)
library(RColorConesa)
library(fmsb)
library(MetBrewer)
library(plotly)


outdir = "output/main"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

#### set theme for plots
pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", linewidth = 0.4),
        axis.line.y = element_line(color="black", linewidth = 0.4)) +
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(vjust=0.5, size=14) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=14)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "none")

old.libplat.palette = c( "cDNA+ONT"="#74CDF0", "cDNA+PacBio"="#EE446F", "cDNA+Illumina"="#FFCF71", 
                         "CapTrap+ONT"="#7482F0", "R2C2+ONT"="#74F0D9", "dRNA+ONT"="#13BF5E", "CapTrap+PacBio"="#d14141")

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)
cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

# Load existing data
# Load and modify ES_code
ES_code <- read.csv("Challenge3_Figures_Data/ES/ES_code.txt", sep=",", header = T )
ES_summary_table <- read.csv("Challenge3_Figures_Data/ES_challenge1/ES_challenge1_metrics.summary_table_SC.csv", header = T)
manatee_code <- read.csv("Challenge3_Figures_Data/manatee/manatee_code.txt", sep=",", header = T )
manatee_metrics <- read.csv("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.challenge3_metrics.csv", sep=",",header = T) %>% t()
ES_metrics_perc <- read.csv("Challenge3_Figures_Data/ES/ES_challenge3_metrics.challenge3_metrics_perc.csv", sep=",",header = T) %>% t()
manatee_SIRV_metrics <- read.csv("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.SIRVS_metrics.csv", sep=",",header = T) %>% t()
df_SC_monoexons <- read.csv("Challenge3_Figures_Data/ES_challenge1/ES_challenge1_metrics.monoexons_SC.csv", sep=",", header = T) %>% t() %>% as.data.frame()

# Check correct structure of meta_data
if(!all(colnames(manatee_code) == colnames(meta_data))){
  print('ERROR: The column names of the meta_data is expected to match the column names of manatee_code: "pipelineCode" "Library_Preps" "Platform" "Data_Category" "Lab" "Tool" "Alias"')
  break
} 

manatee_code <- rbind(manatee_code, meta_data)
if (dataset2_provided != 'NA'){
  manatee_code <- rbind(manatee_code, meta_data2)
}


manatee_code$Lib_Plat <- apply(manatee_code, 1, function(x){
  paste(x["Library_Preps"], x["Platform"], sep = "-")
})
manatee_code$Lib_DC=apply(cbind(manatee_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
manatee_code$Label <-apply(cbind(manatee_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")

ES_code <- rbind(ES_code, meta_data)
if (dataset2_provided != 'NA'){
  ES_code <- rbind(ES_code, meta_data2)
}

ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
  paste(x["Library_Preps"], x["Platform"], sep = "-")
})
ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")

mouse_4a <- function(new_data, new_data_monoexons, new_data2 = NULL, new_data_monoexons2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  ### Panel 4a
  #############
  
  ES_summary_table <- merge(ES_summary_table, ES_code, by.x="ID", by.y="pipelineCode")
  ES_summary_table$Tool <- gsub("-", "\n", ES_summary_table$Tool)
  
  ES_summary_table <- rbind(ES_summary_table, new_data)
  if(!is.null(new_data2)){
    ES_summary_table <- rbind(ES_summary_table, new_data2)
  }
  
  melted_ES_summary <- ES_summary_table %>%
    pivot_longer(c("FSM","ISM","NIC","NNC","Antisense","Fusion","GenicGenomic","GenicIntron","Intergenic"))
  melted_ES_summary$Tool <- gsub("-", "\n", melted_ES_summary$Tool)
  
  melted_ES_summary$name <- melted_ES_summary$name %>%  factor(levels=c("FSM", "ISM", "NIC", "NNC",
                                                                        "Antisense", "Fusion", "GenicGenomic", "GenicIntron",
                                                                        "Intergenic"),
                                                               labels=c("FSM", "ISM", "NIC", "NNC",
                                                                        "Antisense", "Fusion", "GenicGenomic", "GenicIntron",
                                                                        "Intergenic"))
  
  p.A2.Labels= c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                 "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", unique(melted_ES_summary$Label)[length(unique(melted_ES_summary$Label))-1], if (!is.null(new_data2)) unique(melted_ES_summary$Label)[length(unique(melted_ES_summary$Label))] else NULL)
  melted_ES_summary$value <- as.integer(melted_ES_summary$value)
  
  # Convert the Tool column to a factor with the desired levels
  melted_ES_summary$Tool <- factor(melted_ES_summary$Tool, levels = c(sort(setdiff(unique(melted_ES_summary$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    melted_ES_summary <- melted_ES_summary[melted_ES_summary$Tool == tool_name,]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      melted_ES_summary <- melted_ES_summary[melted_ES_summary$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      melted_ES_summary <- melted_ES_summary[melted_ES_summary$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      melted_ES_summary <- melted_ES_summary[melted_ES_summary$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      melted_ES_summary <- melted_ES_summary[melted_ES_summary$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  p.A2 <- ggplot(melted_ES_summary, aes(x=Label, y=value, fill=name)) + 
    geom_bar(position = "fill", stat="identity", width = 0.7) +
    facet_grid( .~Tool , scales="free_x", space="free_x") +
    scale_fill_manual(values =  cat.palette, name="" ) +
    scale_x_discrete(breaks=p.A2.Labels,
                     labels=c("SO", rep("LO", 6), rep("LS", 3), data_category, if (!is.null(new_data2)) data_category2 else NULL))+
    xlab("")+ 
    ylab("Structural Category \nDistribution")+
    pub_theme +
    theme(text = element_text(size = 9)) +
    theme( strip.text.x = element_blank(),
           legend.text = element_text(size=14)) +
    scale_y_continuous(label = unit_format(unit = "%", scale = 100), expand = expansion(mult = c(0,0.1)))+
    #labs(tag = "LO:Long Only\nLS:Long and Short", hjust = 0) +
    theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1 )) +
    theme(plot.margin=unit(c(-1,0,0,0), "cm"), plot.tag.position = c(0.95,0.95), ) 
  
  df_SC_monoexons$total_mono <- apply(df_SC_monoexons,1,sum)
  df_SC_monoexons$ID <- row.names(df_SC_monoexons)
  df_SC_monoexons <- rbind(df_SC_monoexons, new_data_monoexons)
  if(!is.null(new_data_monoexons2)){
    df_SC_monoexons <- rbind(df_SC_monoexons, new_data_monoexons2)
  }
  
  total_mono <- df_SC_monoexons$total_mono
  
  ES_summary_table <- merge(ES_summary_table, df_SC_monoexons[, c("ID", "total_mono")], by="ID")
  ES_summary_table$total_multi <- apply(ES_summary_table,1, function(x){
    as.numeric(x["total"]) - as.numeric(x["total_mono"])
  })
  
  ES_summary_table$total_mono <- as.numeric(ES_summary_table$total_mono)
  melted_ES_summary_exons <- ES_summary_table %>%
    pivot_longer(c("total_multi", "total_mono"))
  melted_ES_summary_exons$Tool <- gsub("-", "\n", melted_ES_summary_exons$Tool)
  
  # Convert the Tool column to a factor with the desired levels
  melted_ES_summary_exons$Tool <- factor(melted_ES_summary_exons$Tool, levels = c(sort(setdiff(unique(melted_ES_summary_exons$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    melted_ES_summary_exons <- melted_ES_summary_exons[melted_ES_summary_exons$Tool == tool_name,]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      melted_ES_summary_exons <- melted_ES_summary_exons[melted_ES_summary_exons$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      melted_ES_summary_exons <- melted_ES_summary_exons[melted_ES_summary_exons$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      melted_ES_summary_exons <- melted_ES_summary_exons[melted_ES_summary_exons$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      melted_ES_summary_exons <- melted_ES_summary_exons[melted_ES_summary_exons$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  p.A1_exons <- ggplot(melted_ES_summary_exons, aes(x=Label, y=value, alpha=name, fill=Lib_Plat))+
    geom_bar(stat="identity", width=0.7)+
    facet_grid( .~Tool , scales="free_x", space="free_x") +
    scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
    scale_alpha_manual(values=c(0.6,1), breaks=c("total_mono", "total_multi"), labels=c("Mono-exons", "Multi-exons"))+
    pub_theme +
    ylab("Num. total isoforms") +
    xlab("")+
    theme( axis.text.x = element_blank(),
           strip.text.x = element_text(size=18)) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "none",
          axis.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0), "cm"))
  
  p.A1_exons <- p.A1_exons + labs(color = NULL)
  ply.A1_exons <- ggplotly(p.A1_exons)
  ply.A2 <- ggplotly(p.A2)
  
  for (i in seq_along(ply.A1_exons$x$data)) {
    ply.A1_exons$x$data[[i]]$showlegend <- FALSE
  }
  
  pA <- suppressWarnings(
    subplot(
      ply.A1_exons, ply.A2,
      nrows = 2,    # Arrange plots in 2 rows for stacking vertically
      titleX = TRUE, # Keep titles for x-axes
      titleY = TRUE  # Keep titles for y-axes
    ) %>% layout(
      width = 1000,  # Set plot width (in pixels)
      height = 500,
      legend = list(
        orientation = 'h', 
        x = 0.5, 
        xanchor = 'center', 
        y = -0.2, 
        title = list(text = NULL)
      ),
      yaxis = list(
        title = list(
          text = "No. total \nisoforms",
          font = list(size = 14),
          standoff = 20
        )
      ),
      yaxis2 = list(
        title = list(
          text = "Structural category \n distribution",
          font = list(size = 14),
          standoff = 20
        )
      )
    )
  )
  
  
  return(pA)
  
  ### End Panel 4a
  ################
}

### Script for figure 4b
append_new_data_4b <- function(new_data, new_data2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){ 

  manatee_metrics <- rbind(manatee_metrics, new_data)
  
  if(!is.null(new_data2)){
    manatee_metrics <- rbind(manatee_metrics, new_data2)
    rownames(manatee_metrics)[nrow(manatee_metrics)-1] <- as.character(manatee_code$pipelineCode[nrow(manatee_code)-1])
    rownames(manatee_metrics)[nrow(manatee_metrics)] <- as.character(manatee_code$pipelineCode[nrow(manatee_code)])
  } else {
    rownames(manatee_metrics)[nrow(manatee_metrics)] <- as.character(manatee_code$pipelineCode[nrow(manatee_code)])
  }
  
  # -- Create mono_exons_df base
  base_rows   <- c("ONT1","PB1","ONT2", "illumina1","ONT3","ONT4","PB2","PB3","ONT5")
  base_counts <- c(0,7,79174,635028,3703,470801,52,246565,85291)
  
  mono_exons_df <- data.frame(
    Row.names     = base_rows,
    Num.monoexons = base_counts
  )
  
  if(!is.null(new_data2)){
    new_data_row <- data.frame(
      Row.names     = manatee_code$pipelineCode[nrow(manatee_code)-1],
      Num.monoexons = sum(temp_env$sqanti_data$subcategory == 'mono-exon')
    )
    mono_exons_df <- rbind(mono_exons_df, new_data_row)
    new_data_row2 <- data.frame(
      Row.names     = manatee_code$pipelineCode[nrow(manatee_code)],
      Num.monoexons = sum(temp_env2$sqanti_data$subcategory == 'mono-exon')
    )
    mono_exons_df <- rbind(mono_exons_df, new_data_row2)
  } else {
    new_data_row <- data.frame(
      Row.names     = manatee_code$pipelineCode[nrow(manatee_code)],
      Num.monoexons = sum(temp_env$sqanti_data$subcategory == 'mono-exon')
    )
    mono_exons_df <- rbind(mono_exons_df, new_data_row)
  }
  
  manatee_metrics <- merge(manatee_metrics, mono_exons_df, by.x=0, by.y="Row.names")
  
  manatee_metrics$Num.multiexons <- apply(manatee_metrics, 1, function(x){
    as.numeric(x["Number of transcripts"]) - as.numeric(x["Num.monoexons"])
  })
  
  manatee_metrics <- merge(manatee_metrics, manatee_code, by.x="Row.names", by.y="pipelineCode")
  
  manatee_metrics$Tool <- gsub("-", "\n", manatee_metrics$Tool)
  colnames(manatee_metrics) <- make.names(colnames(manatee_metrics))
  
  melted_manatee <- manatee_metrics %>%
    pivot_longer(c("Num.monoexons", "Num.multiexons"))
  
  melted_manatee$Tool <- gsub("-", "\n", melted_manatee$Tool)
  melted_manatee <- melted_manatee %>% filter(Row.names!="ONT1")
  melted_manatee$Row.names <- as.factor(melted_manatee$Row.names)
  
  # -- Extract # transcripts from new_data
  number_of_tr_new_tool <- as.integer(
    manatee_metrics[manatee_metrics$Row.names == meta_data$pipelineCode, ]$Number.of.transcripts
  )
  example_label <- sprintf("%.1fK", number_of_tr_new_tool / 1000)
  
  example_label2 <- NULL
  if(!is.null(new_data2)){
    number_of_tr_new_tool2 <- as.integer(
      manatee_metrics[manatee_metrics$Row.names == meta_data2$pipelineCode, ]$Number.of.transcripts
    )
    example_label2 <- sprintf("%.1fK", number_of_tr_new_tool2 / 1000)
  }
  
  # -- Adjust factor levels so that tool_name is last
  all_tools <- sort(setdiff(unique(melted_manatee$Tool), tool_name))
  melted_manatee$Tool <- factor(melted_manatee$Tool, levels = c(all_tools, tool_name))
  
  # -- Build x-axis breaks and labels for pF4.1
  breaks_vec <- c("PB1","ONT3","PB2","ONT2","ONT5","illumina1","ONT4","PB3", meta_data$pipelineCode)
  labels_vec <- c("1.9K","63K","25K","179K","177K","916K","543K","294K", example_label)
  if(!is.null(new_data2)){
    breaks_vec <- c(breaks_vec, meta_data2$pipelineCode)
    labels_vec <- c(labels_vec, example_label2)
  }
  
  if(comparison == 'NA'){
    melted_manatee <- melted_manatee[melted_manatee$Tool == tool_name,]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      melted_manatee <- melted_manatee[melted_manatee$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      melted_manatee <- melted_manatee[melted_manatee$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      melted_manatee <- melted_manatee[melted_manatee$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      melted_manatee <- melted_manatee[melted_manatee$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  pF4.1 <- ggplot(melted_manatee, aes(x=Row.names, y=value, alpha=name, fill=Lib_Plat)) +
    geom_bar(stat="identity", position="stack", width=0.7) +
    facet_grid(. ~ Tool, scales='free_x', space="free_x") +
    scale_fill_manual(values = libplat.palette, name="Library-Platoform") +
    scale_alpha_manual(
      values=c(0.6,1),
      breaks=c("Num.monoexons", "Num.multiexons"),
      labels=c("Mono-exons", "Multi-exons")
    ) +
    scale_x_discrete(breaks=breaks_vec, labels=labels_vec) +
    pub_theme +
    ylab("Num. total isoforms") +
    xlab("") +
    theme(
      strip.background = element_rect(fill="gray90", color="gray50", linewidth=0.8),
      strip.text.x     = element_text(size=10, margin=margin(t=10, b=10)),
      axis.text.x      = element_text(size=8)
    ) +
    scale_y_continuous(
      label = unit_format(unit="K", scale=0.001),
      expand = expansion(mult = c(0, 0.1))
    ) +
    theme(
      legend.position = "none",
      plot.margin     = unit(c(0,0,-0.5,0), "cm")
    ) +
    labs(color = NULL)
  
  # ---------- Transcript categories for new_data only
  transcripts_per_gene <- temp_env$sqanti_data %>%
    group_by(associated_gene)  %>%
    summarize(num_transcripts = n_distinct(associated_gene))
  
  transcript_categories <- transcripts_per_gene %>%
    mutate(
      category = case_when(
        num_transcripts == 1 ~ "1",
        num_transcripts %in% c(2,3) ~ "2-3",
        num_transcripts %in% c(4,5) ~ "4-5",
        num_transcripts >= 6 ~ ">6"
      )
    ) %>%
    count(category)
  
  transcript_categories <- as.data.frame(transcript_categories)
  colnames(transcript_categories) <- c('Cat', 'Num_locus')
  
  if(!is.null(new_data2)){
    transcripts_per_gene2 <- temp_env2$sqanti_data %>%
      group_by(associated_gene) %>%
      summarize(num_transcripts = n_distinct(associated_transcript))
    
    transcript_categories2 <- transcripts_per_gene2 %>%
      mutate(
        category = case_when(
          num_transcripts == 1 ~ "1",
          num_transcripts %in% c(2,3) ~ "2-3",
          num_transcripts %in% c(4,5) ~ "4-5",
          num_transcripts >= 6 ~ ">6"
        )
      ) %>%
      count(category)
    
    transcript_categories2 <- as.data.frame(transcript_categories2)
    colnames(transcript_categories2) <- c('Cat', 'Num_locus')
  }
  
  # ---------- Process info number isoforms per locus
  trx_locus.list <- list()
  for(i in c("ONT1","ONT2","ONT3","ONT4","ONT5","PB1","PB2","PB3","illumina1")){
    trx_file <- file.path("Challenge3_Figures_Data", "manatee", 
                          "manatee_isoforms_per_gene", paste0(i, "_cpm_vs_trans.tsv"))
    df <- read.csv(trx_file, sep="\t", header=TRUE)
    num_trx_locus <- data.frame(
      Cat       = c("1","2-3","4-5",">6"),
      Num_locus = c(
        sum(df$n_isoforms == 1),
        sum(df$n_isoforms %in% c(2,3)),
        sum(df$n_isoforms %in% c(4,5)),
        sum(df$n_isoforms > 5)
      )
    )
    trx_locus.list[[i]] <- num_trx_locus
  }
  
  # -- Append transcript_categories (new_data)
  trx_locus.list[[length(trx_locus.list) + 1]] <- transcript_categories
  names(trx_locus.list)[[length(trx_locus.list)]] <- pipelineCode

  if(!is.null(new_data2)){
    trx_locus.list[[length(trx_locus.list) + 1]] <- transcript_categories2
    names(trx_locus.list)[[length(trx_locus.list)]] <- pipelineCode2
  }
  
  trx_locus.df <- bind_rows(trx_locus.list, .id = "pipeline")
  trx_locus.df <- merge(trx_locus.df, manatee_code, by.x="pipeline", by.y="pipelineCode")
  trx_locus.df$Tool <- gsub("-", "\n", trx_locus.df$Tool)
  
  # -- Adjust factor levels
  all_tools2 <- sort(setdiff(unique(trx_locus.df$Tool), tool_name))
  trx_locus.df$Tool <- factor(trx_locus.df$Tool, levels = c(all_tools2, tool_name))
  
  trx_locus.df <- trx_locus.df %>% filter(pipeline!="ONT1")
  
  # -- Build x-axis breaks and labels for pF4.2
  breaks_vec2 <- c("PB1","ONT3","PB2","ONT2","ONT5","illumina1","ONT4","PB3", pipelineCode)
  labels_vec2 <- c("LO","LO","LO","LO","LS","SO","LS","LS", data_category)
  if(!is.null(new_data2)){
    breaks_vec2 <- c(breaks_vec2, pipelineCode2)
    labels_vec2 <- c(labels_vec2, data_category2)
  }
  
  if(comparison == 'NA'){
    trx_locus.df <- trx_locus.df[trx_locus.df$Tool == tool_name,]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      trx_locus.df <- trx_locus.df[trx_locus.df$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      trx_locus.df <- trx_locus.df[trx_locus.df$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      trx_locus.df <- trx_locus.df[trx_locus.df$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      trx_locus.df <- trx_locus.df[trx_locus.df$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  trx_locus.df[24,] <- NA
  trx_locus.df[22:24, 'pipeline'] <- pipelineCode
  trx_locus.df$Cat[22:24] <- c('2-3', '4-5', '>6')
  trx_locus.df$Num_locus[21:24] <- trx_locus.df$Num_locus[13:16]
  trx_locus.df[22:24,4:12] <- trx_locus.df[21,4:12]
  
  pF4.2 <- ggplot(trx_locus.df, aes(x=pipeline, y=Num_locus, fill=Cat)) +
    geom_bar(position="fill", stat="identity", width=0.7) +
    facet_grid(. ~ Tool, scales="free_x", space="free_x", drop=TRUE) +
    scale_fill_met_d(name="Cassatt1", direction=-1) +
    scale_x_discrete(breaks=breaks_vec2, labels=labels_vec2) +
    ylab("Number of transcripts \n per gene") +
    pub_theme +
    theme(
      strip.text.x      = element_blank(),
      axis.text.x       = element_text(size=8),
      axis.title.x      = element_blank(),
      legend.text       = element_text(size=14),
      legend.title      = element_blank()
    ) +
    scale_y_continuous(
      label = unit_format(unit="%", scale=100),
      expand = expansion(mult = c(0,0.1))
    ) +
    theme(
      legend.position = "bottom",
      plot.margin     = unit(c(-1,0,0,0), "cm")
    )
  
  plyF4.1 <- ggplotly(pF4.1)
  plyF4.2 <- ggplotly(pF4.2)
  
  # -- Remove duplicate legends
  for(i in seq_along(plyF4.1$x$data)){
    plyF4.1$x$data[[i]]$showlegend <- FALSE
  }
  
  pF <- suppressWarnings(
    subplot(
      plyF4.1, plyF4.2,
      nrows  = 2,    # stacked vertically
      titleX = TRUE, # keep x-axis titles
      titleY = TRUE  # keep y-axis titles
    ) %>% 
      layout(
        width  = 1000,
        height = 500,
        legend = list(
          orientation = 'h',
          x = 0.5,
          xanchor = 'center',
          y = -0.2,
          title = list(text = NULL)
        ),
        yaxis = list(
          title = list(
            text  = "No. total \nisoforms",
            font  = list(size = 14),
            standoff = 20
          )
        ),
        yaxis2 = list(
          title = list(
            text  = "Number of transcripts \nper gene",
            font  = list(size = 14),
            standoff = 20
          )
        )
      )
  )
  
  return(pF)
}

### Script for figure 4c
append_new_data_4c <- function(new_data, new_data2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
  ######### length distribution mouse
  dist.list_ES <- list()
  for (i in c("ONT_1", "ONT_2","ONT_3","ONT_4","ONT_5",
              "ONT_6", "ONT_7","ONT_8","ONT_9","ONT_10",
              "ONT_11","PB_1","PB_2","PB_3",
              "PB_4", "PB_5", "illumina_1")){
    length_file=paste0("Challenge3_Figures_Data/ES/length_distributions/",i,"_length_dist.tsv")
    df=read.csv(length_file,sep="\t",header = T)
    dist.list_ES[[i]] <- as.numeric(df$length) %>% as.data.frame()
  }
  
  
  dist.list_ES[[ES_code$pipelineCode[nrow(ES_code)]]] <- as.numeric(new_data$length) %>% as.data.frame()
  if (!is.null(new_data2)) {
    dist.list_ES[[ES_code$pipelineCode[nrow(ES_code)-1]]] <- as.numeric(new_data$length) %>% as.data.frame()
    dist.list_ES[[ES_code$pipelineCode[nrow(ES_code)]]] <- as.numeric(new_data2$length) %>% as.data.frame()
  }
  
  dist_df_ES <- bind_rows(dist.list_ES, .id = "pipeline")
  colnames(dist_df_ES)<-c("pipeline","length")
  dist_df_ES <- merge(dist_df_ES, ES_code, by.x="pipeline", by.y="pipelineCode")
  dist_df_ES$Tool <- gsub("-", "\n", dist_df_ES$Tool)
  
  label_breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                   "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", ES_code$Label[nrow(ES_code)-1], 
                   if (!is.null(new_data2)) ES_code$Label[nrow(ES_code)] else NULL)
  labels = c("SO", rep("LO", 6), rep("LS", 3), data_category, if (!is.null(new_data2)) data_category2 else NULL)
  
  # print counts to use in caption, order to match plots
  count_df_ES <- dist_df_ES %>%
    group_by(pipeline, Label, Tool) %>%
    summarise(Count = n(), .groups = 'drop')
  count_df_ES <- count_df_ES[order(tolower(count_df_ES$Tool), tolower(count_df_ES$Label)),]
  
  # Convert the Tool column to a factor with the desired levels
  dist_df_ES$Tool <- factor(dist_df_ES$Tool, levels = c(sort(setdiff(unique(dist_df_ES$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    dist_df_ES <- dist_df_ES[dist_df_ES$Tool == tool_name,]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      dist_df_ES <- dist_df_ES[dist_df_ES$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      dist_df_ES <- dist_df_ES[dist_df_ES$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      dist_df_ES <- dist_df_ES[dist_df_ES$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      dist_df_ES <- dist_df_ES[dist_df_ES$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  pF3 <- ggplot(dist_df_ES, aes(x=Label, y=length, fill=Lib_Plat))+
    geom_violin(alpha = 0.6) +
    geom_boxplot(width=0.15, color="violet", alpha=0.8, outlier.shape = NA) +
    facet_grid( .~Tool , scales="free_x", space="free_x") +
    scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
    scale_x_discrete(breaks=label_breaks, labels=labels)+
    xlab("")+ 
    ylab("Length (bp), log10")+
    pub_theme +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size=10)) +
    scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)) )
  
  return(suppressWarnings(ggplotly(pF3, width = 1000, height = 500)))
  
}

### Script for figure 4d
append_new_data_4d <- function(new_data, new_data2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
  dist.list_manatee <- list()
  for (i in c("ONT1", "ONT2","ONT3","ONT4","ONT5","PB1","PB2","PB3","illumina1")){
    length_file=paste0("Challenge3_Figures_Data/manatee/length_distributions/",i,"_length_dist.tsv")
    df=read.csv(length_file,sep="\t",header = T)
    dist.list_manatee[[i]] <- as.numeric(df$length) %>% as.data.frame()
  }
  
  
  dist.list_manatee[[manatee_code$pipelineCode[nrow(manatee_code)]]] <- as.numeric(new_data$length) %>% as.data.frame()
  if (!is.null(new_data2)) {
    dist.list_manatee[[manatee_code$pipelineCode[nrow(manatee_code)-1]]] <- as.numeric(new_data$length) %>% as.data.frame()
    dist.list_manatee[[manatee_code$pipelineCode[nrow(manatee_code)]]] <- as.numeric(new_data2$length) %>% as.data.frame()
  }
  
  dist_df_manatee <- bind_rows(dist.list_manatee, .id = "pipeline")
  colnames(dist_df_manatee)<-c("pipeline","length")
  dist_df_manatee <- merge(dist_df_manatee, manatee_code, by.x="pipeline", by.y="pipelineCode")
  dist_df_manatee$Tool <- gsub("-", "\n", dist_df_manatee$Tool)
  
  sample_size = dist_df_manatee %>% group_by(pipeline) %>% summarize(num=n())
  label_breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                   "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", if (!is.null(new_data2)) manatee_code$Label[nrow(manatee_code)-1] else manatee_code$Label[nrow(manatee_code)], 
                   if (!is.null(new_data2)) manatee_code$Label[nrow(manatee_code)] else NULL)
  labels = c("SO", rep("LO", 6), rep("LS", 3), data_category, if (!is.null(new_data2)) data_category2 else NULL)
  
  # print counts to use in caption, order to match plots
  count_df_manatee <- dist_df_manatee %>%
    group_by(pipeline, Label, Tool) %>%
    summarise(Count = n(), .groups = 'drop')
  count_df_manatee <- count_df_manatee[order(tolower(count_df_manatee$Tool), tolower(count_df_manatee$Label)),]
  
  # Convert the Tool column to a factor with the desired levels
  dist_df_manatee$Tool <- factor(dist_df_manatee$Tool, levels = c(sort(setdiff(unique(dist_df_manatee$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    dist_df_manatee <- dist_df_manatee[dist_df_manatee$Tool == tool_name,]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      dist_df_manatee <- dist_df_manatee[dist_df_manatee$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      dist_df_manatee <- dist_df_manatee[dist_df_manatee$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      dist_df_manatee <- dist_df_manatee[dist_df_manatee$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      dist_df_manatee <- dist_df_manatee[dist_df_manatee$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  pF2 <- ggplot(dist_df_manatee, aes(x=Label, y=length, fill=Lib_Plat))+
    geom_violin(alpha = 0.6) +
    geom_boxplot(width=0.15, color="violet", alpha=0.8, outlier.shape = NA) +
    facet_grid( .~Tool , scales="free_x", space="free_x") +
    scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
    scale_x_discrete(breaks=label_breaks, labels=labels)+
    xlab("")+ 
    ylab("Length (bp), log10")+
    pub_theme +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size=10)) +
    scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)) )
  plyF2 <- suppressWarnings(ggplotly(pF2, width = 1000, height = 500))
  return(plyF2)
  
}

### Script for figure 4e1
append_new_data_4e1 <- function(new_data, new_data2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
  if(organism_like == 'mouse'){
    ES_metrics_perc <- read.csv(paste0("Challenge3_Figures_Data/ES/ES_challenge3_metrics.challenge3_metrics_perc.csv"), sep=",",header = T) %>% t()
  } else {
    ES_metrics_perc <- read.csv(paste0("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.challenge3_metrics_perc.csv"), sep=",",header = T) %>% t()
  }
  
  new_row <- rep(NA, ncol(ES_metrics_perc)) 
  names(new_row) <- colnames(ES_metrics_perc)  
  new_row[colnames(new_data)] <- new_data[2,]
  
  if (!is.null(new_data2)) {
    new_row2 <- rep(NA, ncol(ES_metrics_perc)) 
    names(new_row2) <- colnames(ES_metrics_perc)  
    new_row2[colnames(new_data2)] <- new_data2[2,]
    suppressWarnings(ES_metrics_perc <- rbind(ES_metrics_perc, new_row, new_row2))
  } else {
    suppressWarnings(ES_metrics_perc <- rbind(ES_metrics_perc, new_row))
  }
  
  rownames(ES_metrics_perc)[nrow(ES_metrics_perc)] <- pipelineCode
  if (!is.null(new_data2)) {
    rownames(ES_metrics_perc)[nrow(ES_metrics_perc)-1] <- pipelineCode
    rownames(ES_metrics_perc)[nrow(ES_metrics_perc)] <- pipelineCode2
  }
  
  if(organism_like == 'mouse'){
    ES_code <- read.csv(paste0("Challenge3_Figures_Data/ES/ES_code.txt"), sep=",", header = T )
  } else {
    ES_code <- read.csv(paste0("Challenge3_Figures_Data/manatee/manatee_code.txt"), sep=",", header = T )
  }
  ES_code <- rbind(ES_code, meta_data)
  if (!is.null(new_data2)) {
    ES_code <- rbind(ES_code, meta_data2)
  }
  
  ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")
  
  ES_metrics_perc <- merge(ES_metrics_perc, ES_code, by.x=0, by.y="pipelineCode")
  ES_metrics_perc$Tool <- gsub("-", "\n", ES_metrics_perc$Tool)
  colnames(ES_metrics_perc) <- make.names(colnames(ES_metrics_perc))
  ES_metrics_perc[nrow(ES_metrics_perc), "Mapping.transcripts"] <- 100
  
  # Convert the Tool column to a factor with the desired levels
  ES_metrics_perc$Tool <- factor(ES_metrics_perc$Tool, levels = c(sort(setdiff(unique(ES_metrics_perc$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    ES_metrics_perc <- ES_metrics_perc[-c(1:nrow(ES_metrics_perc)-1),]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  pC_SJ <- ggplot(ES_metrics_perc, aes(x=Label, y=as.numeric(Splice.Junctions.with.short.read.coverage), color=Lib_Plat, shape=Data_Category)) +  
    geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(Splice.Junctions.with.short.read.coverage), color=Lib_Plat), linewidth = 2) +
    geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
    pub_theme+
    theme(axis.title.y = element_text(size=12),
          axis.text.y  = element_text(size=12) ) +
    scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
    scale_color_manual(values =  libplat.palette, name="Library-Platform") +
    scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                              "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", if(!is.null(new_data2)) ES_metrics_perc$Label[2] else ES_metrics_perc$Label[1],
                              if (!is.null(new_data2)) ES_metrics_perc$Label[1] else NULL),
                     labels=c("SO", rep("LO", 6), rep("LS", 3), if(!is.null(new_data2)) ES_metrics_perc$Data_Category[2] else ES_metrics_perc$Data_Category[1],
                              if (!is.null(new_data2)) ES_metrics_perc$Data_Category[1] else NULL))+
    facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
    xlab("") + ylab("")+ 
    scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "none",
          strip.text.x = element_text(size = 10))
  
  return(ggplotly(pC_SJ, width = 1000, height = 500) %>% layout(yaxis=list(title=list(text="% Splice Junctions with short read coverage", standoff=20))))
}

### Script for figure 4e3
append_new_data_4e3 <- function(new_data, new_data2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
 if(organism_like == 'mouse'){
    ES_metrics_perc <- read.csv(paste0("Challenge3_Figures_Data/ES/ES_challenge3_metrics.challenge3_metrics_perc.csv"), sep=",",header = T) %>% t()
  } else {
    ES_metrics_perc <- read.csv(paste0("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.challenge3_metrics_perc.csv"), sep=",",header = T) %>% t()
  }
  
  new_row <- rep(NA, ncol(ES_metrics_perc)) 
  names(new_row) <- colnames(ES_metrics_perc)  
  new_row[colnames(new_data)] <- new_data[2,]
  
  if (!is.null(new_data2)) {
    new_row2 <- rep(NA, ncol(ES_metrics_perc)) 
    names(new_row2) <- colnames(ES_metrics_perc)  
    new_row2[colnames(new_data2)] <- new_data2[2,]
    suppressWarnings(ES_metrics_perc <- rbind(ES_metrics_perc, new_row, new_row2))
  } else {
    suppressWarnings(ES_metrics_perc <- rbind(ES_metrics_perc, new_row))
  }
  
  rownames(ES_metrics_perc)[nrow(ES_metrics_perc)] <- pipelineCode
  if (!is.null(new_data2)) {
    rownames(ES_metrics_perc)[nrow(ES_metrics_perc)-1] <- pipelineCode
    rownames(ES_metrics_perc)[nrow(ES_metrics_perc)] <- pipelineCode2
  }

  if(organism_like == 'mouse'){
    ES_code <- read.csv(paste0("Challenge3_Figures_Data/ES/ES_code.txt"), sep=",", header = T )
  } else {
    ES_code <- read.csv(paste0("Challenge3_Figures_Data/manatee/manatee_code.txt"), sep=",", header = T )
  }
  
  ES_code <- rbind(ES_code, meta_data)
  if (!is.null(new_data2)) {
    ES_code <- rbind(ES_code, meta_data2)
  }
  
  ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")
  
  ES_metrics_perc <- merge(ES_metrics_perc, ES_code, by.x=0, by.y="pipelineCode")
  ES_metrics_perc$Tool <- gsub("-", "\n", ES_metrics_perc$Tool)
  colnames(ES_metrics_perc) <- make.names(colnames(ES_metrics_perc))
  ES_metrics_perc[nrow(ES_metrics_perc), "Mapping.transcripts"] <- 100
  
  # Convert the Tool column to a factor with the desired levels
  ES_metrics_perc$Tool <- factor(ES_metrics_perc$Tool, levels = c(sort(setdiff(unique(ES_metrics_perc$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    ES_metrics_perc <- ES_metrics_perc[-c(1:nrow(ES_metrics_perc)-1),]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      ES_metrics_perc <- ES_metrics_perc[ES_metrics_perc$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  pC_nonCan <- ggplot(ES_metrics_perc, aes(x=Label, y=as.numeric(Non.canonical.Splice.Junctions), color=Lib_Plat, shape=Data_Category))  + 
    geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(Non.canonical.Splice.Junctions), color=Lib_Plat), linewidth=2) +
    geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
    pub_theme+
    theme(axis.title.y = element_text(size=12),
          axis.text.y  = element_text(size=12) ) +
    scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
    scale_color_manual(values =  libplat.palette, name="Library-Platform") +
    scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                              "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", if(!is.null(new_data2)) ES_metrics_perc$Label[2] else ES_metrics_perc$Label[1],
                              if (!is.null(new_data2)) ES_metrics_perc$Label[1] else NULL),
                     labels=c("SO", rep("LO", 6), rep("LS", 3), if(!is.null(new_data2)) ES_metrics_perc$Data_Category[2] else ES_metrics_perc$Data_Category[1],
                              if (!is.null(new_data2)) ES_metrics_perc$Data_Category[1] else NULL))+
    facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
    xlab("") + ylab("")+ 
    scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "none",
          strip.text.x = element_text(size = 10))
  
  return(ggplotly(pC_nonCan, width = 1000, height = 500) %>% layout(yaxis=list(title=list(text="% Non-canonical Splice Junctions", standoff=20))))
}

### Script for figure 4f
append_new_data_4f1_ES <- function(new_data_ES, new_data_ES2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
  ES_BUSCO <- read.csv(paste0("Challenge3_Figures_Data/ES/ES_challenge3_metrics.BUSCO_metrics.csv"), sep=",",header = T) %>% t()
  new_ES_BUSCO <- array(NA, dim = c(nrow(ES_BUSCO) + ifelse(is.null(new_data_ES2), 1, 2), ncol(ES_BUSCO)))
  new_ES_BUSCO[1:nrow(ES_BUSCO), ] <- ES_BUSCO
  new_ES_BUSCO[nrow(new_ES_BUSCO)-ifelse(is.null(new_data_ES2), 0, 1), ] <- new_data_ES$`Absolute value`
  
  if (!is.null(new_data_ES2)) {
    new_ES_BUSCO[nrow(new_ES_BUSCO), ] <- new_data_ES2$`Absolute value`
  }
  
  colnames(new_ES_BUSCO) <- colnames(ES_BUSCO)
  rownames(new_ES_BUSCO) <- c(rownames(ES_BUSCO), pipelineCode, if (!is.null(new_data_ES2)) pipelineCode2 else NULL)
  ES_BUSCO <- new_ES_BUSCO
  
  ES_BUSCO <- merge(ES_BUSCO, ES_code, by.x=0, by.y="pipelineCode")
  ES_BUSCO$Tool <- gsub("-", "\n", ES_BUSCO$Tool)
  colnames(ES_BUSCO) <- make.names(colnames(ES_BUSCO))
  
  total_busco <- 11366
  get_perc_busco_found <- function(x){
    r=(total_busco - as.numeric(x["Missing.BUSCOs"]))*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_complete <- function(x){
    c <- as.numeric(x["Complete.and.single.copy.BUSCOs"]) + as.numeric(x["Complete.and.duplicated.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_fragmented <- function(x){
    c <- as.numeric(x["Fragmented.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_missing <- function(x){
    c <- as.numeric(x["Missing.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  ES_BUSCO$BUSCO_found <- apply(ES_BUSCO,1,get_perc_busco_found)
  ES_BUSCO$BUSCO_complete <- apply(ES_BUSCO,1,get_perc_busco_complete)
  
  ES_BUSCO$BUSCO_fragmented <- apply(ES_BUSCO,1,get_perc_busco_fragmented)
  ES_BUSCO$BUSCO_missing <- apply(ES_BUSCO,1,get_perc_busco_missing)
  
  ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(sort(setdiff(unique(ES_BUSCO$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    ES_BUSCO <- ES_BUSCO[-c(1:nrow(ES_BUSCO)-1),]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  # Modify `pD7.2` plot
  pD7.1 <- ggplot(ES_BUSCO, aes(x=Label, y=as.numeric(BUSCO_found), color=Lib_Plat, shape=Data_Category)) +  
    geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_found), color=Lib_Plat), size=2) +
    geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
    pub_theme +
    theme_pubclean(flip=TRUE) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=14),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0,0,0,0)  # Remove margins around the plot
    ) +
    scale_fill_manual(values = libplat.palette, name="Library-Platform") +
    scale_color_manual(values = libplat.palette, name="Library-Platform") +
    scale_x_discrete(breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                                "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
                     labels = c("SO", rep("LO", 6), rep("LS", 3))) +
    facet_grid(Tool ~., drop = TRUE, scales = "free_y") +
    xlab("BUSCO genes found (%)") +
    scale_y_continuous("", sec.axis = sec_axis(~ ., breaks = NULL, name = "Mouse"),
                       label = unit_format(unit = "%"), limits = c(-1, 80), expand = expansion(mult = c(0, 0.1))) +
    theme(
      strip.text.y = element_blank(), 
      axis.line.y = element_blank(),
      axis.title = element_text(size=16)
    ) +
    theme(legend.position = "none") +
    coord_flip()
  
  # Calculate the appropriate height for the boxes in p.mid7 based on the number of levels
  num_tools <- length(unique(ES_BUSCO$Tool))  # Number of unique tools
  ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(tool_name ,sort(setdiff(unique(ES_BUSCO$Tool), tool_name), decreasing = TRUE)))
  
  # Create the left plot
  p.mid7 <- ggplot(ES_BUSCO, aes(y=Tool, x=1)) +
    geom_tile(fill="gray", color="white", width=0.5, height=0.95) +
    geom_text(aes(label=Tool), size=5) +
    theme_void() +
    theme(
      plot.margin=margin(10, 0, 25, -60),
      legend.position="none",
    ) +
    scale_x_continuous(expand=c(0, 0))
  
  # Combine the plots with adjusted widths
  pD7 <- suppressWarnings(subplot(
    ggplotly(p.mid7, width = 200, height = 500),  # Left plot: 200 units wide
    ggplotly(pD7.1, width = 800, height = 500), # Right plot: 800 units wide
    nrows = 1,
    widths = c(0.2, 0.8),  # Proportions for left and right plots
    titleX = TRUE,
    titleY = TRUE
  ) %>% layout(
    width = 1000,  # Total plot width
    height = 500,   # Total plot height
    xaxis2 = list(title = list(text = 'BUSCO Genes Found (%)', font = 14, standoff = 20))
  ))
  
  return(pD7)
  
}

append_new_data_4f2_ES <- function(new_data_ES, new_data_ES2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
  ES_BUSCO <- read.csv(paste0("Challenge3_Figures_Data/ES/ES_challenge3_metrics.BUSCO_metrics.csv"), sep=",",header = T) %>% t()
  new_ES_BUSCO <- array(NA, dim = c(nrow(ES_BUSCO) + ifelse(is.null(new_data_ES2), 1, 2), ncol(ES_BUSCO)))
  new_ES_BUSCO[1:nrow(ES_BUSCO), ] <- ES_BUSCO
  new_ES_BUSCO[nrow(new_ES_BUSCO)-ifelse(is.null(new_data_ES2), 0, 1), ] <- new_data_ES$`Absolute value`
  
  if (!is.null(new_data_ES2)) {
    new_ES_BUSCO[nrow(new_ES_BUSCO), ] <- new_data_ES2$`Absolute value`
  }
  
  colnames(new_ES_BUSCO) <- colnames(ES_BUSCO)
  rownames(new_ES_BUSCO) <- c(rownames(ES_BUSCO), pipelineCode, if (!is.null(new_data_ES2)) pipelineCode2 else NULL)
  ES_BUSCO <- new_ES_BUSCO
  
  ES_BUSCO <- merge(ES_BUSCO, ES_code, by.x=0, by.y="pipelineCode")
  ES_BUSCO$Tool <- gsub("-", "\n", ES_BUSCO$Tool)
  colnames(ES_BUSCO) <- make.names(colnames(ES_BUSCO))
  
  total_busco <- 11366
  get_perc_busco_found <- function(x){
    r=(total_busco - as.numeric(x["Missing.BUSCOs"]))*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_complete <- function(x){
    c <- as.numeric(x["Complete.and.single.copy.BUSCOs"]) + as.numeric(x["Complete.and.duplicated.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_fragmented <- function(x){
    c <- as.numeric(x["Fragmented.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_missing <- function(x){
    c <- as.numeric(x["Missing.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  ES_BUSCO$BUSCO_found <- apply(ES_BUSCO,1,get_perc_busco_found)
  ES_BUSCO$BUSCO_complete <- apply(ES_BUSCO,1,get_perc_busco_complete)
  
  ES_BUSCO$BUSCO_fragmented <- apply(ES_BUSCO,1,get_perc_busco_fragmented)
  ES_BUSCO$BUSCO_missing <- apply(ES_BUSCO,1,get_perc_busco_missing)
  
  ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(sort(setdiff(unique(ES_BUSCO$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    ES_BUSCO <- ES_BUSCO[-c(1:nrow(ES_BUSCO)-1),]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      ES_BUSCO <- ES_BUSCO[ES_BUSCO$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  pD7.1f <- ggplot(ES_BUSCO, aes(x = Label, y = BUSCO_fragmented, color = Lib_Plat, shape = Data_Category)) +
    geom_segment(aes(x = Label, xend = Label, y = 0, yend = BUSCO_fragmented, color = Lib_Plat), size = 2) +
    geom_point(position = position_dodge(width = 1), size = 5, aes(fill = Lib_Plat)) +
    scale_fill_manual(values = libplat.palette, name = "Library-Platform") +
    scale_color_manual(values = libplat.palette, name = "Library-Platform") +
    pub_theme +
    theme_pubclean(flip = TRUE) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, -20, 0, 0),
      legend.position = "none",
      strip.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title = element_text(size = 16)
    ) +
    scale_x_discrete(
      breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                 "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
      labels = c("SO", rep("LO", 6), rep("LS", 3))
    ) +
    scale_y_continuous(
      limits = c(-1, 20), label = unit_format(unit = "%"),
      expand = expansion(mult = c(0, 0.1))
    ) +
    facet_grid(Tool ~ ., drop = TRUE, scales = "free_y") +
    xlab("% BUSCO Fragmented") +
    ylab("") +
    coord_flip()
  
  # Calculate the appropriate height for the boxes in p.mid7 based on the number of levels
  num_tools <- length(unique(ES_BUSCO$Tool))  # Number of unique tools
  ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(tool_name ,sort(setdiff(unique(ES_BUSCO$Tool), tool_name), decreasing = TRUE)))
  
  # Create the left plot
  p.mid7 <- ggplot(ES_BUSCO, aes(y=Tool, x=1)) +
    geom_tile(fill="gray", color="white", width=0.5, height=0.95) +
    geom_text(aes(label=Tool), size=5) +
    theme_void() +
    theme(
      plot.margin=margin(10, 0, 25, -60),
      legend.position="none",
    ) +
    scale_x_continuous(expand=c(0, 0))
  
  # Combine the plots with adjusted widths
  pD7.2 <- suppressWarnings(subplot(
    ggplotly(p.mid7, width = 200, height = 500),  # Left plot: 200 units wide
    ggplotly(pD7.1f, width = 800, height = 500), # Right plot: 800 units wide
    nrows = 1,
    widths = c(0.2, 0.8),  # Proportions for left and right plots
    titleX = TRUE,
    titleY = TRUE
  ) %>% layout(
    width = 1000,  # Total plot width
    height = 500,   # Total plot height
    xaxis2 = list(title = list(text = 'BUSCO Genes Fragmented (%)', font = 14, standoff = 20))
  ))
  
  return(pD7.2)
  
}

# For using organism: Manatee
append_new_data_4f1_manatee <- function(new_data_manatee, new_data_manatee2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant){
  
  manatee_BUSCO <- read.csv(paste0("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.BUSCO_metrics.csv"), sep=",",header = T) %>% t()
  new_manatee_BUSCO <- array(NA, dim = c(nrow(manatee_BUSCO) + ifelse(is.null(new_data_manatee2), 1, 2), ncol(manatee_BUSCO)))
  new_manatee_BUSCO[1:nrow(manatee_BUSCO), ] <- manatee_BUSCO
  new_manatee_BUSCO[nrow(new_manatee_BUSCO)-ifelse(is.null(new_data_manatee2), 0, 1), ] <- new_data_manatee$`Absolute value`
  
  if (!is.null(new_data_manatee2)) {
    new_manatee_BUSCO[nrow(new_manatee_BUSCO), ] <- new_data_manatee2$`Absolute value`
  }
  
  colnames(new_manatee_BUSCO) <- colnames(manatee_BUSCO)
  rownames(new_manatee_BUSCO) <- c(rownames(manatee_BUSCO), pipelineCode, if (!is.null(new_data_manatee2)) pipelineCode2 else NULL)
  manatee_BUSCO <- new_manatee_BUSCO
  
  manatee_BUSCO <- merge(manatee_BUSCO, manatee_code, by.x=0, by.y="pipelineCode")
  manatee_BUSCO$Tool <- gsub("-", "\n", manatee_BUSCO$Tool)
  colnames(manatee_BUSCO) <- make.names(colnames(manatee_BUSCO))
  manatee_BUSCO <- manatee_BUSCO %>% filter(Row.names!="ONT1")
  
  total_busco <- 11366
  get_perc_busco_found <- function(x){
    r=(total_busco - as.numeric(x["Missing.BUSCOs"]))*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_complete <- function(x){
    c <- as.numeric(x["Complete.and.single.copy.BUSCOs"]) + as.numeric(x["Complete.and.duplicated.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_fragmented <- function(x){
    c <- as.numeric(x["Fragmented.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  get_perc_busco_missing <- function(x){
    c <- as.numeric(x["Missing.BUSCOs"]) 
    r=c*100/total_busco 
    round(r, digits = 2)
  }
  
  manatee_BUSCO$BUSCO_found <- apply(manatee_BUSCO,1,get_perc_busco_found)
  manatee_BUSCO$BUSCO_complete <- apply(manatee_BUSCO,1,get_perc_busco_complete)
  manatee_BUSCO$BUSCO_fragmented <- apply(manatee_BUSCO,1,get_perc_busco_fragmented)
  manatee_BUSCO$BUSCO_missing <- apply(manatee_BUSCO,1,get_perc_busco_missing)
  
  manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(sort(setdiff(unique(manatee_BUSCO$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    manatee_BUSCO <- manatee_BUSCO[-c(1:nrow(manatee_BUSCO)-1),]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  # Modify `pD7.2` plot
  pD7.2 <- ggplot(manatee_BUSCO, aes(x=Label, y=as.numeric(BUSCO_found), color=Lib_Plat, shape=Data_Category)) + 
    geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_found), color=Lib_Plat), size=2) +
    geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
    pub_theme +
    theme_pubclean(flip=TRUE) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=14),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0,0,0,0)  # Remove margins around the plot
    ) +
    scale_fill_manual(values = libplat.palette, name="Library-Platform") +
    scale_color_manual(values = libplat.palette, name="Library-Platform") +
    scale_x_discrete(breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                                "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
                     labels = c("SO", rep("LO", 6), rep("LS", 3))) +
    facet_grid(Tool ~., drop = TRUE, scales = "free_y") +
    xlab("BUSCO genes found (%)") +
    scale_y_continuous("", sec.axis = sec_axis(~ ., breaks = NULL, name = "Manatee"),
                       label = unit_format(unit = "%"), limits = c(-1, 80), expand = expansion(mult = c(0, 0.1))) +
    theme(
      strip.text.y = element_blank(), 
      axis.line.y = element_blank(),
      axis.title = element_text(size=16)
    ) +
    theme(legend.position = "none") +
    coord_flip()
  
  # Calculate the appropriate height for the boxes in p.mid7 based on the number of levels
  num_tools <- length(unique(manatee_BUSCO$Tool))  # Number of unique tools
  manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(tool_name ,sort(setdiff(unique(manatee_BUSCO$Tool), tool_name), decreasing = TRUE)))
  
  # Create the left plot
  p.mid <- ggplot(manatee_BUSCO, aes(y=Tool, x=1)) +
    geom_tile(fill="gray", color="white", width=0.5, height=0.95) +
    geom_text(aes(label=Tool), size=5) +
    theme_void() +
    theme(
      plot.margin=margin(10, 0, 25, -60),
      legend.position="none",
    ) +
    scale_x_continuous(expand=c(0, 0))
  
  # Combine the plots with adjusted widths
  combined_plot <- suppressWarnings(subplot(
    ggplotly(p.mid, width = 200, height = 500),  # Left plot: 200 units wide
    ggplotly(pD7.2, width = 800, height = 500), # Right plot: 800 units wide
    nrows = 1,
    widths = c(0.2, 0.8),  # Proportions for left and right plots
    titleX = TRUE,
    titleY = TRUE
  ) %>% layout(
    width = 1000,  # Total plot width
    height = 500,   # Total plot height
    xaxis2 = list(title = list(text = 'BUSCO Genes Found (%)', font = 14, standoff = 20))
  ))
  
  return(combined_plot)
}

append_new_data_4f2_manatee <- function(new_data_manatee, new_data_manatee2 = NULL, comparison, bambu, rna_bloom, rnaspades, stringtie2_isoquant) {
  
  # Load and process the data
  manatee_BUSCO <- read.csv(paste0("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.BUSCO_metrics.csv"), sep=",", header=TRUE) %>% t()
  new_manatee_BUSCO <- array(NA, dim=c(nrow(manatee_BUSCO) + ifelse(is.null(new_data_manatee2), 1, 2), ncol(manatee_BUSCO)))
  new_manatee_BUSCO[1:nrow(manatee_BUSCO), ] <- manatee_BUSCO
  new_manatee_BUSCO[nrow(new_manatee_BUSCO)-ifelse(is.null(new_data_manatee2), 0, 1), ] <- new_data_manatee$`Absolute value`
  
  if (!is.null(new_data_manatee2)) {
    new_manatee_BUSCO[nrow(new_manatee_BUSCO), ] <- new_data_manatee2$`Absolute value`
  }
  
  colnames(new_manatee_BUSCO) <- colnames(manatee_BUSCO)
  rownames(new_manatee_BUSCO) <- c(rownames(manatee_BUSCO), pipelineCode, if (!is.null(new_data_manatee2)) pipelineCode2 else NULL)
  manatee_BUSCO <- new_manatee_BUSCO
  
  manatee_BUSCO <- merge(manatee_BUSCO, manatee_code, by.x=0, by.y="pipelineCode")
  manatee_BUSCO$Tool <- gsub("-", "\n", manatee_BUSCO$Tool)
  colnames(manatee_BUSCO) <- make.names(colnames(manatee_BUSCO))
  manatee_BUSCO <- manatee_BUSCO %>% filter(Row.names != "ONT1")
  
  # Compute percentages
  total_busco <- 11366
  manatee_BUSCO$BUSCO_fragmented <- apply(manatee_BUSCO, 1, function(x) round(as.numeric(x["Fragmented.BUSCOs"]) * 100 / total_busco, 2))
  
  # Ensure the same scale for Tool
  manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(sort(setdiff(unique(manatee_BUSCO$Tool), tool_name)), tool_name))
  
  if(comparison == 'NA'){
    manatee_BUSCO <- manatee_BUSCO[-c(1:nrow(manatee_BUSCO)-1),]
  } else if(comparison == 'Custom'){
    if(bambu == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'Bambu',]
    }
    if(rna_bloom == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'RNA\nBloom',]
    }
    if(rnaspades == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'rnaSPAdes',]
    }
    if(stringtie2_isoquant == 'NA'){
      manatee_BUSCO <- manatee_BUSCO[manatee_BUSCO$Tool != 'StringTie2\nIsoQuant',]
    }
  }
  
  p.right <- ggplot(manatee_BUSCO, aes(x = Label, y = BUSCO_fragmented, color = Lib_Plat, shape = Data_Category)) +
    geom_segment(aes(x = Label, xend = Label, y = 0, yend = BUSCO_fragmented, color = Lib_Plat), size = 2) +
    geom_point(position = position_dodge(width = 1), size = 5, aes(fill = Lib_Plat)) +
    scale_fill_manual(values = libplat.palette, name = "Library-Platform") +
    scale_color_manual(values = libplat.palette, name = "Library-Platform") +
    pub_theme +
    theme_pubclean(flip = TRUE) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, -20, 0, 0),
      legend.position = "none",
      strip.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title = element_text(size = 16)
    ) +
    scale_x_discrete(
      breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                 "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
      labels = c("SO", rep("LO", 6), rep("LS", 3))
    ) +
    scale_y_continuous(
      limits = c(-1, 20), label = unit_format(unit = "%"),
      expand = expansion(mult = c(0, 0.1))
    ) +
    facet_grid(Tool ~ ., drop = TRUE, scales = "free_y") +
    xlab("% BUSCO Fragmented") +
    ylab("") +
    coord_flip()
  
  # Calculate the appropriate height for the boxes in p.mid7 based on the number of levels
  num_tools <- length(unique(manatee_BUSCO$Tool))  # Number of unique tools
  manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(tool_name ,sort(setdiff(unique(manatee_BUSCO$Tool), tool_name), decreasing = TRUE)))
  
  # Create the left plot
  p.left <- ggplot(manatee_BUSCO, aes(y=Tool, x=1)) +
    geom_tile(fill="gray", color="white", width=0.5, height=0.95) +
    geom_text(aes(label=Tool), size=5) +
    theme_void() +
    theme(
      plot.margin=margin(10, 0, 25, -60),
      legend.position="none",
    ) +
    scale_x_continuous(expand=c(0, 0))
  
  # Combine the plots with adjusted widths
  combined_plot <- suppressWarnings(subplot(
    ggplotly(p.left, width = 200, height = 500),  # Left plot: 200 units wide
    ggplotly(p.right, width = 800, height = 500), # Right plot: 800 units wide
    nrows = 1,
    widths = c(0.2, 0.8),  # Proportions for left and right plots
    titleX = TRUE,
    titleY = TRUE
  ) %>% layout(
    width = 1000,  # Total plot width
    height = 500,   # Total plot height
    xaxis2 = list(title = list(text = 'BUSCO Genes Fragmented (%)', font = 14, standoff = 20))
  ))
  
  return(combined_plot)
}


print('Scripts loaded')