#####################################
##### SQANTI3 report generation ######
#####################################



### Francisco J Pardo-Palacios
### Last Modified: 05/03/2020 by francisco.pardo.palacios@gmail.com

#********************** Taking arguments from python script

results_file1 <- "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file1/"
results_file1 <- '/Users/woutermaessen/Desktop/ConesaInternship/Platform_runs/Challenge3/manatee/StringTie2_cDNA_PB_LO.fasta/results_file1/'



args <- c(paste0(results_file1, "lrgasp_platform_challenge_3_classification.txt"),
          paste0(results_file1, "lrgasp_platform_challenge_3_junctions.txt"),
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/lrgasp_event2_metrics/utilities",
          paste0(results_file1, "lrgasp_platform_challenge_3_Rdata"),
          paste0(results_file1, "BUSCO_results.tsv"),
          "manatee",
          "Stringtie_testyour_tool",
          "PacBio",
          "cDNA",
          "LO",
          "Custom",
          "NA",
          "NA",
          "NA",
          "Stringtie2",
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
          "NA")

args <- c(paste0(results_file1, "lrgasp_platform_challenge_3_classification.txt"),
          paste0(results_file1, "lrgasp_platform_challenge_3_junctions.txt"),
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/lrgasp_event2_metrics/utilities",
          paste0(results_file1, "lrgasp_platform_challenge_3_Rdata"),
          paste0(results_file1, "BUSCO_results.tsv"),
          "manatee",
          "Stringtie_testyour_tool",
          "PacBio",
          "cDNA",
          "LO",
          "Custom",
          "NA",
          "NA",
          "NA",
          "Stringtie2",
          "NA",
          "NA",
          "NA",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file2/lrgasp_platform_challenge_3_classification.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file2/lrgasp_platform_challenge_3_junctions.txt",
          "/Users/woutermaessen/PycharmProjects/lrgasp_platform/sqanti_results/results_file2/BUSCO_results.tsv",
          "PacBio",
          "cDNA",
          "LO",
          "NA",
          "NA",
          "NA",
          "Provided")

args <- commandArgs(trailingOnly = TRUE)
args

# for file 1
class.file <- args[1]
junc.file <- args[2]
utilities.path <- args[3]
rdata <- args[4]
busco <- args[5]
organism <- args[6]
if (args[7] == 'your_tool'){
  tool_name <- "your_tool" 
} else{
  tool_name <- gsub("your_tool", "", args[7])
}
name <- tool_name
alias <- tool_name
platform <- args[8]
library_prep <- args[9]
data_category <- args[10]
comparison <- args[11]
bambu <- args[12]
rna_bloom <- args[13]
rnaspades <- args[14]
stringtie2_isoquant <- args[15]
sirv_file <- args[16]
ercc_file <- args[17]
sequin_file <- args[18]
class.file2 <- args[19]
junc.file2 <- args[20]
busco2 <- args[21]
platform2 <- args[22]
library_prep2 <- args[23]
data_category2 <- args[24]
sirv_file2 <- args[25]
ercc_file2 <- args[26]
sequin_file2 <- args[27]
dataset2_provided <- args[28]

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

list_alias <- c("Ba", "Bl", "rS", "IQ")

if (alias %in% list_alias){
  alias <- paste0(alias, "2")
}

# Generate pipelineCode
pipelineCode <- paste0(alias,'_',platform,'_',library_prep)
pipelineCode2 <- paste0(alias,'_',platform2,'_',library_prep2,"_2")

lab <- 'testing_lab'

meta_data <- data.frame(pipelineCode = pipelineCode, Library_Preps = library_prep, Platform = platform, Data_Category = data_category, Lab = lab, Tool = tool_name, Alias = alias)
meta_data2 <- data.frame(pipelineCode = pipelineCode2, Library_Preps = library_prep2, Platform = platform2, Data_Category = data_category2, Lab = lab, Tool = tool_name, Alias = alias)


report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
report.file <- paste(report.prefix, "Evaluation_report.html", sep="_");
bam_file <- paste(report.prefix, "corrected.bam", sep="_")


#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "knitr","rmarkdown", "Rsamtools")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(knitr)
library(rmarkdown)
library(Rsamtools)

#********************* Run Calculation scripts

setwd(utilities.path)
source("Platform_calculations.R")

#file 1
LRGASP_calculations_challenge3(NAME = name, out.dir = rdata,
                               class.file=class.file, junc.file=junc.file,
                               platform = platform, 
                               functions.dir = utilities.path,
                               bam = bam_file,
                               organism = organism,
                               sirv_list = sirv_file,
                               ercc_list = ercc_file,
                               sequin_list = sequin_file)

#file 2
if (dataset2_provided != 'NA'){
  bam_file2 <- paste(strsplit(class.file2, "_classification.txt")[[1]][1], "corrected.bam", sep="_")
  LRGASP_calculations_challenge3(NAME = paste0(name, '2'), out.dir = rdata,
                                 class.file=class.file2, junc.file=junc.file2,
                                 platform = platform2, 
                                 functions.dir = utilities.path,
                                 bam = bam_file2,
                                 organism = organism,
                                 sirv_list = sirv_file2,
                                 ercc_list = ercc_file2,
                                 sequin_list = sequin_file2)
}



#file 1
busco_table = read.table(busco, sep="\t", header=F)
busco_table = as.data.frame(t(busco_table))
busco_results = data.frame(row.names = busco_table[,1])
busco_results[,"Absolute value"]=apply(busco_table, 1, function(Y) as.integer(Y[2]))
total_BUSCO = sum(busco_results[,"Absolute value"])
busco_results[,"Relative value (%)"] = apply(busco_results,1, function(Z){
  round( ((Z[1]/total_BUSCO)*100), digits = 2)
})
save(busco_results , file = paste(name, "_BUSCO.RData", sep = ''))

#file 2
if (dataset2_provided != 'NA'){
  busco_table = read.table(busco2, sep="\t", header=F)
  busco_table = as.data.frame(t(busco_table))
  busco_results = data.frame(row.names = busco_table[,1])
  busco_results[,"Absolute value"]=apply(busco_table, 1, function(Y) as.integer(Y[2]))
  total_BUSCO = sum(busco_results[,"Absolute value"])
  busco_results[,"Relative value (%)"] = apply(busco_results,1, function(Z){
    round( ((Z[1]/total_BUSCO)*100), digits = 2)
  })
  save(busco_results , file = paste(paste0(name,'2'), "_BUSCO.RData", sep = ''))
}

temp_env <- new.env()

#file 1
load(paste0(rdata, '/', name, '_classification.RData'), envir = temp_env) #sqanti_data
load(paste0(rdata, '/', name, '_junctions.RData'), envir = temp_env)      #sqanti_data.junc
load(paste0(rdata, '/', name, '_Transcriptome_results.RData'), envir = temp_env)        #transcriptome.results
load(paste0(rdata, '/', name, '_SIRVs_class.RData'), envir = temp_env)    #sirv_data
load(paste0(rdata, '/', name, '_SIRVs_junc.RData'), envir = temp_env)     #sirv_data.junc
load(paste0(rdata, '/', name, '_BUSCO.RData'), envir = temp_env)          #busco_results
if(tail(strsplit(sirv_file, split='/')[[1]],1) != 'NA'){
  load(paste0(rdata, '/', name, '_SIRV_results.RData'), envir = temp_env)  #sirv type: SIRV          
} 
if(tail(strsplit(ercc_file, split='/')[[1]],1) != 'NA'){
  load(paste0(rdata, '/', name, '_ERCC_results.RData'), envir = temp_env)  #sirv type: ERCC  
}
if(tail(strsplit(sequin_file, split='/')[[1]],1) != 'NA'){
  load(paste0(rdata, '/', name, '_SEQUIN_results.RData'), envir = temp_env)  #sirv type: SEQUIN   
}
if(length(unique(temp_env$sqanti_data$structural_category)) == 1){
  temp_env$sqanti_data$structural_category <- NA
}


#file 2
if (dataset2_provided != 'NA'){
  temp_env2 <- new.env()
  load(paste0(rdata, '/', paste0(name,'2'), '_classification.RData'), envir = temp_env2) #sqanti_data
  load(paste0(rdata, '/', paste0(name,'2'), '_junctions.RData'), envir = temp_env2)      #sqanti_data.junc
  load(paste0(rdata, '/', paste0(name,'2'), '_Transcriptome_results.RData'), envir = temp_env2)        #transcriptome.results
  load(paste0(rdata, '/', paste0(name,'2'), '_SIRVs_class.RData'), envir = temp_env2)    #sirv_data
  load(paste0(rdata, '/', paste0(name,'2'), '_SIRVs_junc.RData'), envir = temp_env2)     #sirv_data.junc
  load(paste0(rdata, '/', paste0(name,'2'), '_BUSCO.RData'), envir = temp_env2)          #busco_results
  if(sirv_file2 != 'NA'){
    load(paste0(rdata, '/', paste0(name,'2'), '_SIRV_results.RData'), envir = temp_env2)  #sirv type: SIRV          
  } 
  if(ercc_file2 != 'NA'){
    load(paste0(rdata, '/', paste0(name,'2'), '_ERCC_results.RData'), envir = temp_env2)  #sirv type: ERCC  
  }
  if(sequin_file2 != 'NA'){
    load(paste0(rdata, '/', paste0(name,'2'), '_SEQUIN_results.RData'), envir = temp_env2)  #sirv type: SEQUIN   
  }
  if(length(unique(temp_env2$sqanti_data$structural_category)) == 1){
    temp_env2$sqanti_data$structural_category <- NA
  }
}

if(organism == 'mouse' | !unique(is.na(temp_env$sqanti_data$structural_category))){
  organism_like <- 'mouse'
} else {
  organism_like <- 'manatee'
}


setwd(utilities.path)
source('FigureFunctions.R')

if(dataset2_provided == 'NA'){
  RMD = paste(utilities.path, "Evaluation_metrics.Rmd", sep = "/")
  
  rmarkdown::render(RMD, params = list(
    output.directory = rdata,
    Name = name,
    Platform = platform, 
    library_prep = library_prep,
    data_category = data_category,
    pipelineCode = pipelineCode), output_file = report.file
  )
  
} else{
  print('file 2 scripts')
  
  
  RMD = paste(utilities.path, "Evaluation_metrics2.Rmd", sep = "/")
  
  rmarkdown::render(RMD, params = list(
    output.directory = rdata,
    Name = name,
    Platform = platform,
    platform2 = platform2,
    library_prep = library_prep,
    library_prep2 = library_prep2,
    data_category = data_category,
    data_category2 = data_category2,
    pipelineCode = pipelineCode,
    pipelineCode2 = pipelineCode2), output_file = report.file
  )
}




