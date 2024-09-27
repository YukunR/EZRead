library(iq)
library(data.table)
library(openxlsx2)
source("./spectronuat.R", local = T)
source("./diann.R", local = T)
source("./pd.R", local = T)
source("./maxquant.R", local = T)


#' @param data 各个来源的数据
#' @param dam.qh 数据采集方式或定量头, data aquisition method, character, one of DDA and DIA
#' @param format 数据格式, character, one of pd, mq, diann, sp, autodetect
checkFormat <- function(data, dam.qh, format = "autodetect") {
  if (all(colnames(data)[1: 2] == c("PSMs.Workflow.ID", "PSMs.Peptide.ID"))) {
    res <- list(data = data, dam = dam.qh, format = "Proteome Discoverer", comment = "")
  } else if (all(colnames(data)[1: 2] == c("Sequence", "Length"))) {
    res <- list(data = data, dam = dam.qh, format = "MaxQuant", comment = "")
  } else if (all(colnames(data)[1: 2] == c("File.Name", "Run"))) {
    res <- list(data = data, dam = dam.qh, format = "DIA-NN", comment = "")
  } else if (all(colnames(data)[1: 2] == c("R.Data.Points.per.Peak..MS2..EXT.", "R.MS1.Average.Tolerance..ppm."))) {   # TODO: 标准化所有SP结果
    res <- list(data = data, dam = dam.qh, format = "Spectronuat", comment = "")
  } else {
    print(colnames(data))
  }
  
  if (format == "antodetect") {
    res$comment <- paste0("Your data are in ", res$format, " format.")
  } else if (format != res$format) {
    res$comment <- paste0("Your data are in ", res$format, " format! Please check your data!")
  } else {
    res$comment <- "wow! You select correct format!"
  }
  return(res)
}


createProteinTable <- function(protein.list, method) {
  checkNAinList <- function(data) {
    i <- length(data)
    while (i > 0) {
      if (is.logical(data[[i]])) {
        data[[i]] <- NULL
      }
      i <- i - 1
    }
    return(data)
  }
  if (method == "Mean") {
    protein.list <- checkNAinList(protein.list)
    res <- create_protein_table(protein.list, method = "meanInt")
  } else if (method == "TopN") {
    protein.list <- checkNAinList(protein.list)
    res <- create_protein_table(protein.list, method = "topN", N = 3)
  } else if (method == "Median Polish") {
    res <- create_protein_table(protein.list, method = "median_polish")
  } else if (method == "MaxLFQ") {
    res <- create_protein_table(protein.list, method = "maxLFQ")
  }
  
  res <- cbind(Accession = rownames(res$estimate),
               res$estimate)
  res <- as.data.frame(res)
  res <- na.omit(res)
  res[2: 4] <- 2^apply(res[2: 4], MARGIN = 2, FUN = as.numeric)
  return(res)
}


readData <- function(data.path) {
  return(as.data.frame(fread(data.path, check.names = T)))
}


# ====================================== deprecated method ===========================================
#' #' @param data.long long format data
#' #' @param file.format file format
#' #' @param data.wide 'maxquant': wide format data
#' #' @param intensity.col 'spectronaut': PEP.MS1Quantity or PEP.MS2Quantity
#' makeProteinList <- function(data.long, file.format, normalization = T, data.wide = NULL, intensity.col = NULL) {
#'   if (file.format == "Proteome Discoverer") {
#'     makeProteinListPD(data.long, normalization = normalization)
#'   } else if (file.format == "MaxQuant") {
#'     protein.list <- makeProteinListMaxquant(data.protein.group = data.wide, 
#'                                             data.evidence = data.long, 
#'                                             normalization = normalization)
#'   } else if (file.format == "DIA-NN") {
#'     protein.list <- makeProteinListDIANN(data.long, 
#'                                          normalization = normalization)
#'   } else if (file.format == "Spectronuat") {
#'     protein.list <- makeProteinListSP(data.long, intensity.col = intensity.col)
#'   }
#' }

groupRows <- function(result.list, rows.grouping) {
  wb <- wb_workbook() %>% wb_add_worksheet("resultTable")
  wb$page_setup(summary_row = "Above")
  write_data(wb, "resultTable", result.list, col_names = F)
  withProgress(message = "Group rows", {
    i <- 1 / length(rows.grouping)
    for (rows in rows.grouping) {
      if (length(rows) == 2) {
        wb$group_rows("resultTable", rows = rows[1]: rows[2] - 1, collapsed = T)
        incProgress(i)
      }
    }
  })
  return(wb)
}
