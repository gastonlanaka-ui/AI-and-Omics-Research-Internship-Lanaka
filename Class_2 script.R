# A function to classify genes using log2FC and padj
classify_gene <- function(logFC, padj){
  #operation
  classified_gene <- "Not_significant"
  
  classified_gene <- ifelse(logFC > 1 & padj < 0.05, "Upregulated", classified_gene)
  classified_gene <- ifelse(logFC < -1 & padj < 0.05, "Downregulated", classified_gene)
   
  return(classified_gene)
}
#calling the function

classified_gene1 <- classify_gene(logFC = 1.1, padj = 0.04)

#Automating my workflow with a For Loop, steps
#Import datasets 
#check for missing values
#clean the data
#calculate BMI and save results

#Using loops to automate the process
#Loop Workflow

#1 Define Input and output folders
getwd()
input_dir2 <- "data"
output_dir2 <- "results"

#2.List the files you want to process
files_to_process2 <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

#3. create empty results list
Results_list2 <- list()

#4.for each file,
#Import data for each file
#handle NA values, calulate BMI and save results, both csv and R list

for (file_names2 in files_to_process2) {
  cat("\nProcessing:", file_names2, "\n")
  
  input_file_path2 <- file.path(input_dir2, file_names2)
  #import dataset
  data2 <- read.csv(input_file_path2, header = TRUE)
  #checking for missing values
  
  if("logFC" %in% names(data2)){
    missing_count <- sum(is.na(data2$logFC))
    cat("missing values in 'logFC':", missing_count,"\n")
    data2$logFC[is.na(data2$logFC)] <- 1
  }
  if("padj" %in% names(data2)){
    missing_count <- sum(is.na(data2$padj))
    cat("missing values in 'padj':", missing_count,"\n")
    data2$padj[is.na(data2$padj)] <- 1
  }
  #processing status
  data2$status <- classify_gene(data2$logFC, data2$padj)
  cat("status has been processed successfully. \n")
  
  #Save the Results in R
  Results_list2[[file_names2]] <- data2
  
  #Save in results folder
  output_file_path2 <- file.path(output_dir2, paste0("Gene_status", file_names2))
  write.csv(data2, output_file_path2, row.names = FALSE)
  cat("Results saved to:", output_file_path2, "\n")
}

results2_1 <- Results_list2[[1]]
results2_2 <- Results_list2[[2]]  
 














