#' @title calculator of Pvalue
#' @importFrom survminer surv_cutpoint surv_categorize
#' @importFrom survival Surv survdiff
#' @description
#'     By using this function, you can calculate the Pvalue of all genes you have provided.
#' @param survival the dataframe that contains survival data
#' @param RNA the dataframe that contains the expression data of genes
#' @param result the dataframe which will contains the outcome
#'
#' @return this function will return a dataframe that contains either the gene's ensemble IDs and it's Pvalue.
#' @export
#'
#' @examples
#' library(Oncofilterfast)
#' result <- data.frame(gene = c("A"),Pvalue = c(1))
#' RNA_all_path=system.file("extdata", "TCGA-LGG.htseq_fpkm.tsv", package = "Oncofilterfast")
#' RNA_all=read.csv(RNA_all_path,header=TRUE,sep="\t")
#' rows_to_keep <- apply(RNA_all[, -1], 1, function(row) {
#'   non_zero_count <- sum(row != 0)
#'   total_elements <- length(row)
#'   (non_zero_count / total_elements) >= 0.5
#'   })
#' RNA <- RNA_all[rows_to_keep, ]
#' survival_path=system.file("extdata", "TCGA-LGG.survival.tsv", package = "Oncofilterfast")
#' survival=read.csv(survival_path,header=TRUE,sep="\t")
#' final=calculator(survival=survival,RNA=RNA,result=result)
#' print(nrow(final))
#' filtered_result <- final[final$Pvalue < 0.01, ]
#' print(nrow(filtered_result))
#' print(filtered_result)
#'
calculator=function(survival,RNA,result){
  numofgene=nrow(RNA)
  result <- data.frame(
    gene = c("A"),
    Pvalue = c(1))
  #prepare data(working_table)
  working_table=survival
  working_table=t(working_table)
  working_table.colnames=working_table[1,]
  working_table <- as.data.frame(working_table)
  colnames(working_table) <- working_table[1, ]
  working_table <- working_table[-1, ]
  colnames(working_table) <- gsub("-", ".", colnames(working_table))
  for (i in 1:numofgene){
    #prepare process_df(should be a loop)
    process_df=RNA[i,]
    result[i,1]=process_df[1,1]
    #Aligning tables(should be a loop)
    common_ids <- intersect(colnames(working_table), colnames(process_df))
    process_df_aligned <- process_df[, match(common_ids, colnames(process_df))]
    working_table_aligned <- working_table[, match(common_ids, colnames(working_table))]
    combined_table <- rbind(working_table_aligned, process_df_aligned)
    combined_table=t(combined_table)
    combined_table=as.data.frame(combined_table)
    colnames(combined_table)[4] <- "gene"
    combined_table$gene <- as.numeric(as.character(combined_table$gene))
    combined_table$OS.time <- as.numeric(as.character(combined_table$OS.time))
    combined_table$OS <- as.numeric(as.character(combined_table$OS))
    #calculate P_value
    res.cut <- surv_cutpoint(minprop = 0.05,combined_table,
                             time = "OS.time",
                             event = "OS",
                             variables = "gene"
    )
    res.cat <- surv_categorize(res.cut)
    surv_obj <- Surv(res.cat$OS.time, res.cat$OS)
    surv_diff_res <- survdiff(surv_obj ~ gene, data = res.cat)
    result[i,2]=surv_diff_res$pvalue
  }
  return(result)
}
