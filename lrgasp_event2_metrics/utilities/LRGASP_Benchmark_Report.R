
LRGASP_Benchmark_Report <- function(NAME, class.file, junc.file, out.dir, platform, functions.dir, bam, organism, tool_name, alias, data_category, library_prep, lab){
  
  setwd(rdata)
  #source(paste0(utilities.path,"/FigureFunctions.R"))
  source('/Users/woutermaessen/Desktop/ConesaInternship/ScriptsToRecreateFiguresLRGASP_Manatee/FigureFunctions.R')
  
  temp_env <- new.env()
  
  load(paste0(rdata, '/', name, '_classification.RData'), envir = temp_env) #sqanti_data
  load(paste0(rdata, '/', name, '_junctions.RData'), envir = temp_env)      #sqanti_data.junc
  load(paste0(rdata, '/', name, '_results.RData'), envir = temp_env)        #all.results
  load(paste0(rdata, '/', name, '_SIRVs_class.RData'), envir = temp_env)    #sirv_data
  load(paste0(rdata, '/', name, '_SIRVs_junc.RData'), envir = temp_env)     #sirv_data.junc
  
  meta_data <- data.frame(pipelineCode = pipelineCode, Library_Preps = library_prep, Platform = platform, Data_Category = data_category, Lab = lab, Tool = tool_name, Alias = alias)
  
  # Generate fig 4a (in Progress...)
  
  # Generate fig 4b
  append_new_data_4b(new_data = t(as.data.frame(as.numeric(temp_env$all.results$`Transcriptome without reference`[1:10,1]))))
  
  # Generate fig 4c or 4d
  if(organism == 'mouse'){
    append_new_data_4c(new_data = temp_env$sqanti_data[,c(1,4)])
  } else{
    append_new_data_4d(new_data = temp_env$sqanti_data[,c(1,4)])
  }
  
  # Generate 4e
  ES_metrics_perc <- read.csv("Challenge3_Figures_Data/ES/ES_challenge3_metrics.challenge3_metrics_perc.csv", sep=",",header = T) %>% t()
  new_data <- t(as.data.frame(temp_env$all.results$`Transcriptome without reference`))
  new_row <- rep(NA, ncol(ES_metrics_perc)) 
  names(new_row) <- colnames(ES_metrics_perc)  
  new_row[colnames(new_data)] <- new_data[2,] 
  append_new_data_4e1(new_data = new_data)
  append_new_data_4e3(new_data = new_data)
  
  # Generate 4f
  if(organism == 'mouse'){
    append_new_data_4f1_ES(busco_results)
    append_new_data_4f2_ES(busco_results)
  } else {
    append_new_data_4f1_manatee(busco_results)
    append_new_data_4f2_manatee(busco_results)
  }
  

  # Generate 4g
  append_new_data_4g2(new_data = as.data.frame(t(temp_env$all.results$SIRV)))
}
