##############################
### R SCRIPT: Create function for SIRV Evaluation
### Author: Wouter Maessen
### Initial Day: Monday January 20th 2025
### Last Edit: -
##############################

SIRV_Evaluation <- function(sirv_data, sirv_list, sirv_name, out.dir, NAME){
  sirv_list <- trimws(readLines(sirv_list))
  sirv_data$TP=apply(sirv_data,1,TP_function)
  SIRVs_transcripts=as.integer(length(sirv_data$isoform))
  SIRVs_called=intersect(sirv_data[which(sirv_data$structural_category=="FSM" & 
                                           sirv_data$TP==TRUE),"associated_transcript"],
                         sirv_list)
  TP=length(SIRVs_called)
  RM_isoforms=sirv_data[which(sirv_data$structural_category=="FSM" & 
                                sirv_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  SIRVs_transcripts_incomplete=sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM") & 
                                                 sirv_data$TP==FALSE),"isoform"]
  
  SIRVs_called_wrong_ends=setdiff(intersect(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"],
                                            sirv_list),SIRVs_called)
  PTP=length(SIRVs_called_wrong_ends)
  SIRVs_not_detected=setdiff(sirv_list,sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"])
  FN=length(SIRVs_not_detected)
  FP_sirvs_detected=sirv_data[-which(sirv_data$structural_category %in% c("FSM","ISM")),"isoform"]
  FP=length(FP_sirvs_detected)
  
  SIRVs_redundancy=length(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"isoform"])/length(unique(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"]))
  
  # Write out results
  b.SIRVs_results=data.frame(row.names = c("SIRV transcripts", "True Positive detections (TP)", "SIRV transcripts associated to TP (Reference Match)",
                                           "Partial True Positive detections (PTP)", "SIRV transcripts associated to PTP",
                                           "False Negative (FN)", "False Positive (FP)", 
                                           "Sensitivity", "Precision",
                                           "Non Redundant Precision","Positive Detection Rate",
                                           "False Discovery Rate","Redundancy"))
  b.SIRVs_results[,"Value"]="-"
  b.SIRVs_results["SIRV transcripts","Value"]=SIRVs_transcripts
  b.SIRVs_results["True Positive detections (TP)","Value"]=as.integer(TP)
  b.SIRVs_results["SIRV transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  b.SIRVs_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  b.SIRVs_results["SIRV transcripts associated to PTP","Value"]=as.integer(length(SIRVs_transcripts_incomplete))
  b.SIRVs_results["False Negative (FN)","Value"]=as.integer(FN)
  b.SIRVs_results["False Positive (FP)","Value"]=as.integer(FP)
  b.SIRVs_results["Sensitivity","Value"]=round(TP/length(sirv_list), digits = 2)
  b.SIRVs_results["Precision","Value"]=round(RM/SIRVs_transcripts, digits = 2)
  b.SIRVs_results["Non Redundant Precision","Value"]=round(TP/SIRVs_transcripts, digits = 2)
  b.SIRVs_results["Positive Detection Rate", "Value"]=round(length(unique(c(SIRVs_called,SIRVs_called_wrong_ends)))/length(sirv_list), digits = 2)
  b.SIRVs_results["False Discovery Rate","Value"]=round((SIRVs_transcripts - RM)/SIRVs_transcripts, digits = 2)
  b.SIRVs_results["Redundancy","Value"]=round(SIRVs_redundancy, digits = 2)
  
  
  
  ####Create a list with all results and save all
  ###############################################
  
  files <- ls(pattern = "b.SIRVs_results")
  assign(paste0("sirv.results_", sirv_name), list())
  
  for ( i in 1: length(files) ) {
    assign(paste0("sirv.results_", sirv_name)[[i]], eval(parse(text = files[i])))
  }
  print(paste("out.dir:", out.dir))
  setwd(out.dir)
  
  save(list = paste0("sirv.results_", sirv_name), file = paste(NAME,sirv_name, "results.RData", sep = '_'))
}
