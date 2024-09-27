
# ------------------------------ diann -------------------------------
library(diann)

#' @param data DIANN report.tsv file
#' @param dcm data aquisition method, character, dia or dda
diannLFQ <- function(data, dam) {
  if (dam == "DIA") {
    lfq <- diann_maxlfq(data[data$Q.Value <= 0.01 & data$PG.Q.Value <= 0.01,], 
                                   group.header = "Protein.Group", 
                                   id.header = "Precursor.Id", 
                                   quantity.header = "Precursor.Normalised", 
                                   sample.header = "Run")
  } else if (dam == "DDA") {
    lfq <- diann_maxlfq(data[data$Q.Value <= 0.01 & data$PG.Q.Value <= 0.01,], 
                                   group.header = "Protein.Group", 
                                   id.header = "Precursor.Id", 
                                   quantity.header = "Ms1.Area", 
                                   sample.header = "Run")
  }
  lfq <- as.data.frame(lfq) %>% 
    mutate(Accession = row.names(lfq)) %>% 
    select(Accession, everything())
  return(lfq)
}


#' @param data.evidence evidence file
#' @param data.protein.group protein group file
#' @param normalization bool. Whether do median normalization.
makeProteinListDIANN <- function(data.report, normalization = T) {
  data.report <- data.report %>% 
    filter(!is.na(Q.Value) & Q.Value < 0.01 & 
           !is.na(PG.Q.Value) & PG.Q.Value < 0.01)
  norm.data <- preprocess(data.report, 
                          sample_id  = "Run", 
                          primary_id = "Protein.Group", 
                          secondary_id = "Precursor.Id", 
                          intensity_col = "Precursor.Normalised", 
                          pdf_out = NULL, 
                          median_normalization = normalization)
  protein.list <- create_protein_list(norm.data)
  return(protein.list)
}


library(openxlsx2)
library(dplyr)
library(stringr)
library(hash)
load("./modification_unimod.RData")

makeResultListDIANN <- function(protein.table, data.report) {
  # 定义表头
  protein.table <- as.data.table(protein.table)
  data.report <- as.data.table(data.report)
  protein.header1 <- c("First Protein ID", "IDs of protein", "Protein Descriptions", "Gene names", "Score")
  protein.header2 <- c("Unique peptides")
  protein.table.header <- colnames(protein.table)[-1]
  peptide.header <- c("Modified sequence", "Peptide sequence", "Charges", "RT (min)", "Intensity", 
                      "Best MS/MS scan", "Missed cleavages", "Modification")
  
  # >>> 初始化返回值 >>>
  result.table <- list(data.frame(t(c(protein.header1, protein.table.header, protein.header2))))
  rows.grouping <- list()
  
  protein.groups <- unique(data.report$Protein.Group)
  
  # >> 检查用谁做padding
  n.col <- max(9, 6 + length(protein.table.header))
  if (length(protein.table.header) < 3) {
    padding <- "protein"
    pad.columns <- setNames(replicate(3 - length(protein.table.header), "", simplify = FALSE),
                            paste0("X", 7 + length(protein.table.header): 9))
    result.table <- result.table %>% mutate(!!!pad.columns)
    peptide.header <- c("", peptide.header)
  } else if (length(protein.table.header) > 3) {
    padding <- "peptide"
    pad.columns <- setNames(replicate(length(protein.table.header) - 3, "", simplify = FALSE),
                            paste0("X", 10: 6 + length(protein.table.header)))
    peptide.header <- c("", peptide.header, rep("", length(protein.table.header) - 3))
  } else {
    padding <- "none"
    peptide.header <- c("", peptide.header, rep("", length(protein.table.header) - 3))
  }
  
  row.n <- 2
  round.n <- length(protein.groups)
  suppressWarnings({
    withProgress(message = "Extract info", {
      for (i in 1: round.n) {
        # >> 计算pep信息
        accession <- protein.groups[i]
        protein.group.info.tmp <- data.report[Protein.Group == accession]
        # 获取 intensity
        intensity <- protein.table[Accession == accession, -1, with = FALSE]
        if (nrow(intensity) == 0) {
          intensity <- as.data.table(matrix(NA, nrow = 1, ncol = length(protein.table.header)))
          colnames(intensity) <- protein.table.header
        }
        
        # 处理 protein.data
        protein.data <- protein.group.info.tmp[1, .(Protein.id = str_split_i(Protein.Group, ";", 1), 
                                                    protein.ids = Protein.Ids, 
                                                    description = First.Protein.Description, 
                                                    genenames = str_split_i(Genes, ";", 1), 
                                                    score = CScore)]
        
        protein.data <- cbind(protein.data, intensity)
        protein.data[, unique.peptide := NA]
        # 处理 peptide.data
        peptide.data <- protein.group.info.tmp[, {
          bestIndex <- which.max(Quantity.Quality)
          .(
            peptide.sequence = unique(Stripped.Sequence), 
            charges = paste(unique(Precursor.Charge), collapse = "; "), 
            RT = RT[bestIndex][1], 
            intensity = Precursor.Normalised[bestIndex][1], 
            bestmsms = paste(Run[bestIndex][1], 
                             MS2.Scan[bestIndex][1]), 
            missed.cleavage = as.double(sum(ifelse(gregexpr("[KR](?!P|$)", Stripped.Sequence[1], perl = TRUE)[[1]] > 0,
                                                   attr(gregexpr("[KR](?!P|$)", Stripped.Sequence[1], perl = TRUE)[[1]], "match.length"), 0))),
            modification = paste(sapply(gsub("\\(UniMod:(\\d+)\\)", "\\1", 
                                             str_extract_all(Modified.Sequence, "\\(.*?\\)", simplify = FALSE)[[1]]), 
                                        function(x) if (length(x)) modification.unimod[[x]] else "", 
                                        simplify = TRUE, USE.NAMES = FALSE), collapse = "; ")
          )}, by = Modified.Sequence]
        peptide.data <- peptide.data
        # 更新 unique.peptide
        protein.data$unique.peptide <- length(which(str_split(peptide.data$charges, ";", simplify = T) != ""))
      
        if (padding == "protein") {
          protein.data <- protein.data %>% 
            mutate(!!!pad.columns)
          
          peptide.data <- peptide.data %>% 
            mutate(X1 = "") %>% 
            select(X1, everything())
          
        } else if (padding == "peptide") {
          peptide.data <- peptide.data %>% 
            mutate(X1 = "", !!!pad.columns) %>% 
            select(X1, everything())
          
        } else {
          peptide.data <- peptide.data %>% 
            mutate(X1 = "") %>% 
            select(X1, everything())
        }
        
        # 更新列名
        setnames(protein.data, paste0("X", 1:n.col))
        setnames(peptide.data, paste0("X", 1:n.col))
        
        # 合并结果
        result.table <- c(result.table, list(protein.data, data.table(t(peptide.header)), peptide.data))
        row.n <- row.n + 1L
        rows.grouping[[i]] <- c(row.n, 0)
        
        # 添加 peptide.header
        row.n <- row.n + 1L
        
        # 合并 peptide.data
        row.n <- row.n + nrow(peptide.data)
        rows.grouping[[i]][2] <- row.n
        incProgress(1 / round.n)
      }
      result.table <- rbindlist(result.table, use.names = F)
    })
  })
  return(list(result.table = result.table, rows.grouping = rows.grouping))
}


# data.report <- readData("../LibrarySearchingCompare/data/new/dia/diann/report.tsv")
# protein.table1 <- diannLFQ(data.report, "DIA")
# result1 <- system.time({
#   res <- createResultTableDIANN(protein.table1, data.report)
# })
# result1
# result2 <- system.time({
#   groupRows(res$result.table, res$rows.grouping, "./diann_test.xlsx")
# })
# 
# result2
