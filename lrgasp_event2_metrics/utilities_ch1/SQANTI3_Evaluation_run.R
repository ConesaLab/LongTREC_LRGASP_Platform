#####################################
##### SQANTI3 report generation ######
#####################################



### Francisco J Pardo-Palacios
### Last Modified: 05/03/2020 by francisco.pardo.palacios@gmail.com

#********************** Taking arguments from python script
args <- c("/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file1/lrgasp_platform_challenge_1_classification.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file1/lrgasp_platform_challenge_1_junctions.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/lrgasp_event2_metrics/utilities_ch1",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file1/lrgasp_platform_challenge_1_Rdata",
          "mouse",
          "Bambu_testyour_tool",
          "PacBio",
          "cDNA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "False",
          "Custom",
          "NA",
          "Brooks",
          "CRG",
          "MaxPlanck",
          "NA",
          "IB",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA")

args <- c("/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file1/lrgasp_platform_challenge_1_classification.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file1/lrgasp_platform_challenge_1_junctions.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/lrgasp_event2_metrics/utilities_ch1",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file2/lrgasp_platform_challenge_1_Rdata",
          "human",
          "Bambu_testyour_tool",
          "PacBio",
          "cDNA",
          "NA",
          "NA",
          "NA",
          "NA",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file2/lrgasp_platform_challenge_1_classification.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file2/lrgasp_platform_challenge_1_junctions.txt",
          "ONT",
          "R2C2",
          "LO",
          "NA",
          "NA",
          "NA",
          "True",
          "Custom",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NA",
          "NAd",
          "NA",
          "NA",
          "NA",
          "NA")

args <- commandArgs(trailingOnly = TRUE)

args

class.file <- args[1]
junc.file <- args[2]
utilities.path <- args[3]
rdata <- args[4]
organism <- args[5]
if (args[6] == 'your_tool'){
  tool_name <- "your_tool" 
} else{
  tool_name <- gsub("your_tool", "", args[6])
}
name <- tool_name
alias <- tool_name
platform <- args[7]
library_prep <- args[8]
data_category <- args[9]
sirv_file <- args[10]
ercc_file <- args[11]
sequin_file <- args[12]
class.file2 <- args[13]
junc.file2 <- args[14]
platform2 <- args[15]
library_prep2 <- args[16]
data_category2 <- args[17]
sirv_file2 <- args[18]
ercc_file2 <- args[19]
sequin_file2 <- args[20]
dataset2 <- ifelse(args[21] == "True", TRUE, FALSE)
comparison <- args[22]
bambu <- args[23]
FLAIR <- args[24]
Lyric <- args[25]
IsoTools <- args[26]
Mandalorion <- args[27]
Iso_IB <- args[28]
FLAMES <- args[29]
IsoQuant <- args[30]
Spectra <- args[31]
TALON_LAPA <- args[32]
StringTie2 <- args[33]


# Generate Alias
chars <- unlist(strsplit(tool_name, ""))
caps <- chars[grepl("[A-Z]", chars)]

if (length(caps) >= 2) {
  alias <- paste0(caps[1:2], collapse = "")  # Take first 3 capital letters
} else if (length(caps) == 1) {
  # Find the position of the first capital letter
  cap_pos <- regexpr("[A-Z]", name)
  if (cap_pos > 0 && cap_pos < nchar(name)) {
    alias <- substr(name, cap_pos, cap_pos + 1)  # Take the capital letter and the next letter
  } else {
    alias <- caps[1]  # If it's the last letter, just take it
  }
} else {
  # No capital letters, take first 3 letters of 'name' as they are
  alias <- substr(name, 1, 2)
}

list_alias <- c("Ba",  "FL",  "Ly*", "IT",  "Ma",  "IB",  "FM",  "IQ",  "Sp",  "TL",  "ST")

if (alias %in% list_alias){
  alias <- paste0(alias, "2")
}

# Generate pipelineCode
pipelineCode <- paste0(alias,'_',platform,'_',library_prep)
pipelineCode2 <- paste0(alias,'_',platform2,'_',library_prep2,"_2")

meta_data_1 <- c("pipelineCode" = pipelineCode, "Library_Preps" = library_prep, "Platform" = platform, "Data_Category" = data_category, "Lab" = paste(name, "'s lab"), "Tool" = name, "Alias" = alias, 'Lib_Plat' = paste0(library_prep, '-', platform), 'Lib_DC' = paste0(library_prep, '-', data_category), 'Label' = paste0(platform, '-', library_prep, '-', data_category))                                       
meta_data_2 <- c("pipelineCode" = pipelineCode2, "Library_Preps" = library_prep2, "Platform" = platform2, "Data_Category" = data_category2, "Lab" = paste(name, "'s lab"), "Tool" = name, "Alias" = alias, 'Lib_Plat' = paste0(library_prep2, '-', platform2), 'Lib_DC' = paste0(library_prep2, '-', data_category2), 'Label' = paste0(platform2, '-', library_prep2, '-', data_category2))                                       

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
report.file <- paste(report.prefix, "Evaluation_report.html", sep="_");


#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "knitr","rmarkdown", "data.table")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(knitr)
library(rmarkdown)
library(data.table)
library(tidyverse)
#********************* Run Calculation scripts

setwd(utilities.path)

source("LRGASP_calculations.R")

LRGASP_calculations(NAME = name, out.dir = rdata,
                    class.file=class.file, junc.file=junc.file,
                    functions.dir = utilities.path)

LRGASP_calculations(NAME = name, out.dir = rdata,
                    class.file=class.file, junc.file=junc.file,
                    functions.dir = utilities.path, sirv_file = sirv_file,
                    ercc_file = ercc_file, sequin_file = sequin_file)

#file 1
temp_env <- new.env()

load(paste0(rdata, '/', name, '_classification.RData'), envir = temp_env) #sqanti_data
load(paste0(rdata, '/', name, '_junctions.RData'), envir = temp_env)      #sqanti_data.junc
load(paste0(rdata, '/', name, '_results.RData'), envir = temp_env)        #all.results
load(paste0(rdata, '/', name, '_SIRVs_class.RData'), envir = temp_env)    #sirv_data
load(paste0(rdata, '/', name, '_SIRVs_junc.RData'), envir = temp_env)     #sirv_data.junc

if(dataset2 == TRUE){
  LRGASP_calculations(NAME = name, out.dir = rdata,
                      class.file=class.file2, junc.file=junc.file2,
                      functions.dir = utilities.path)
}


if(dataset2 == TRUE){
  LRGASP_calculations(NAME = name, out.dir = rdata,
                      class.file=class.file2, junc.file=junc.file2,
                      functions.dir = utilities.path, sirv_file = sirv_file2,
                      ercc_file = ercc_file2, sequin_file = sequin_file2)
}

#file 2
temp_env2 <- new.env()

load(paste0(rdata, '/', name, '_classification.RData'), envir = temp_env2) #sqanti_data
load(paste0(rdata, '/', name, '_junctions.RData'), envir = temp_env2)      #sqanti_data.junc
load(paste0(rdata, '/', name, '_results.RData'), envir = temp_env2)        #transcriptome.results
load(paste0(rdata, '/', name, '_SIRVs_class.RData'), envir = temp_env2)    #sirv_data
load(paste0(rdata, '/', name, '_SIRVs_junc.RData'), envir = temp_env2)     #sirv_data.junc

setwd(utilities.path)
source('FigureFunctions_ch1.R')


RMD = paste(utilities.path, "Evaluation_metrics.Rmd", sep = "/")


rmarkdown::render(RMD, params = list(
  output.directory = rdata,
  Name = name,
  Platform = platform  ), output_file = report.file
)



